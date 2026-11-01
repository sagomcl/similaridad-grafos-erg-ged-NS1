# ============================================================================
# calibracion_cluster.py
#
# Explicacion detallada de este script (que hace cada paso, como se corre y
# que hay que cambiar para extrapolar el metodo a otra proteina):
#     docs/02_explicacion_del_codigo.md
# Explicacion de los conceptos (nodo, grupo farmacoforico, colapso de atomos,
# papel de la geometria 3D, GED, algoritmo hungaro):
#     docs/01_explicacion_del_metodo.md
# ============================================================================
import sys, json, pickle
sys.path.insert(0, '.')
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from scipy.optimize import differential_evolution
from pipeline_final_completo import smiles_to_3d_mol, build_erg, ged_similarity, DEFAULT_THETA

CLUSTER_IDX = int(sys.argv[1])  # 0 = mas grande, 1 = segundo, etc.
np.random.seed(42)

with open('activos_clusters.pkl','rb') as fh:
    data = pickle.load(fh)
smiles_all = data['smiles']
clusters = data['clusters']
cluster = clusters[CLUSTER_IDX]
cluster_smiles = [smiles_all[i] for i in cluster]
print(f"Cluster #{CLUSTER_IDX}: {len(cluster_smiles)} moleculas")

inactivos = pd.read_csv("AID2326_inactivos_muestra5000.csv").sample(n=60, random_state=100+CLUSTER_IDX).reset_index(drop=True)

n_act = len(cluster_smiles)
n_train_act = max(3, int(n_act*0.6))
train_act_smiles = cluster_smiles[:n_train_act]
test_act_smiles = cluster_smiles[n_train_act:]

train_inact = inactivos.iloc[:30]
test_inact = inactivos.iloc[30:60]

def build_graphs(smiles_list, prefix):
    graphs = {}
    for i, s in enumerate(smiles_list):
        try:
            mol3d, _ = smiles_to_3d_mol(s)
            graphs[f"{prefix}_{i}"] = build_erg(mol3d)
        except Exception:
            continue
    return graphs

g_act_train = build_graphs(train_act_smiles, "act")
g_act_test = build_graphs(test_act_smiles, "actT")
g_inact_train = build_graphs(train_inact['Smiles'].tolist(), "inact")
g_inact_test = build_graphs(test_inact['Smiles'].tolist(), "inactT")
print(f"train: {len(g_act_train)} act + {len(g_inact_train)} inact | test: {len(g_act_test)} act + {len(g_inact_test)} inact")

keys = list(DEFAULT_THETA.keys())
train_names = list(g_act_train.keys()) + list(g_inact_train.keys())
train_labels = [1]*len(g_act_train) + [0]*len(g_inact_train)
all_train = {**g_act_train, **g_inact_train}

def score(name, theta_vec, act_graphs):
    theta = dict(zip(keys, theta_vec))
    sims = [ged_similarity(all_train[name], act_graphs[a], theta) for a in act_graphs if a != name]
    return np.mean(sims) if sims else 0.0

def objetivo(theta_vec):
    scores = [score(n, theta_vec, g_act_train) for n in train_names]
    try:
        return -roc_auc_score(train_labels, scores)
    except ValueError:
        return 0.0

bounds = [(0.05,2.5),(0.05,2.5),(0.3,3.0),(0.2,3.0),(0.02,1.5),(0.1,2.0)]
res = differential_evolution(objetivo, bounds, seed=42, maxiter=15, popsize=8, tol=1e-3)
theta_opt = dict(zip(keys, res.x))
auc_train = -res.fun
print(f"AUC train: {auc_train:.3f}")

test_names = list(g_act_test.keys()) + list(g_inact_test.keys())
test_labels = [1]*len(g_act_test) + [0]*len(g_inact_test)
all_test = {**g_act_test, **g_inact_test}

def score_test(name, theta):
    sims = [ged_similarity(all_test[name], g_act_train[a], theta) for a in g_act_train]
    return np.mean(sims) if sims else 0.0

test_scores = [score_test(n, theta_opt) for n in test_names]
try:
    auc_test = roc_auc_score(test_labels, test_scores)
except ValueError:
    auc_test = float('nan')
print(f"AUC test: {auc_test:.3f}  (n_act_test={len(g_act_test)}, n_inact_test={len(g_inact_test)})")

with open(f"cluster{CLUSTER_IDX}_resultado.json","w") as fh:
    json.dump({"cluster_idx": CLUSTER_IDX, "n_activos_cluster": n_act,
               "auc_train": auc_train, "auc_test": auc_test,
               "n_train_act": len(g_act_train), "n_test_act": len(g_act_test),
               "theta": theta_opt}, fh, indent=2)

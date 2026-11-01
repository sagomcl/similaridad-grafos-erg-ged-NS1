# ============================================================================
# calibracion_pubchem_real.py
#
# Explicacion detallada de este script (que hace cada paso, como se corre y
# que hay que cambiar para extrapolar el metodo a otra proteina):
#     docs/02_explicacion_del_codigo.md
# Explicacion de los conceptos (nodo, grupo farmacoforico, colapso de atomos,
# papel de la geometria 3D, GED, algoritmo hungaro):
#     docs/01_explicacion_del_metodo.md
# ============================================================================
import sys, json
sys.path.insert(0, '.')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve

from pipeline_final_completo import (smiles_to_3d_mol, build_erg, ged_similarity,
                                       DEFAULT_THETA)
from scipy.optimize import differential_evolution

np.random.seed(42)

activos = pd.read_csv("AID2326_activos_reales.csv")
inactivos = pd.read_csv("AID2326_inactivos_muestra5000.csv")

# muestras manejables (calcular ErG de miles de moleculas con embedding 3D es costoso)
act_sample = activos.sample(n=35, random_state=42).reset_index(drop=True)
inact_sample = inactivos.sample(n=70, random_state=42).reset_index(drop=True)

# split train/test para que el AUC final sea de GENERALIZACION, no memorizado
act_train, act_test = act_sample.iloc[:15], act_sample.iloc[15:35]
inact_train, inact_test = inact_sample.iloc[:30], inact_sample.iloc[30:70]

def build_graphs(df, label):
    graphs = {}
    for _, r in df.iterrows():
        try:
            mol3d, _ = smiles_to_3d_mol(r["Smiles"])
            graphs[f"{label}_{int(r['CID'])}"] = build_erg(mol3d)
        except Exception:
            continue
    return graphs

print("Construyendo grafos ErG (train)...")
g_act_train = build_graphs(act_train, "act")
g_inact_train = build_graphs(inact_train, "inact")
print(f"  activos train: {len(g_act_train)}  inactivos train: {len(g_inact_train)}")

print("Construyendo grafos ErG (test)...")
g_act_test = build_graphs(act_test, "act")
g_inact_test = build_graphs(inact_test, "inact")
print(f"  activos test: {len(g_act_test)}  inactivos test: {len(g_inact_test)}")

keys = list(DEFAULT_THETA.keys())
train_names = list(g_act_train.keys()) + list(g_inact_train.keys())
train_labels = [1]*len(g_act_train) + [0]*len(g_inact_train)
all_train_graphs = {**g_act_train, **g_inact_train}

def score(name, theta_vec, act_graphs):
    theta = dict(zip(keys, theta_vec))
    sims = [ged_similarity(all_train_graphs[name], act_graphs[a], theta)
            for a in act_graphs if a != name]
    return np.mean(sims) if sims else 0.0

def objetivo(theta_vec):
    scores = [score(n, theta_vec, g_act_train) for n in train_names]
    try:
        return -roc_auc_score(train_labels, scores)
    except ValueError:
        return 0.0

print("Calibrando theta (differential_evolution, esto tarda unos minutos)...")
bounds = [(0.05, 2.5), (0.05, 2.5), (0.3, 3.0), (0.2, 3.0), (0.02, 1.5), (0.1, 2.0)]
res = differential_evolution(objetivo, bounds, seed=42, maxiter=15, popsize=8, tol=1e-3)
theta_opt = dict(zip(keys, res.x))
auc_train = -res.fun
print(f"AUC en TRAIN (40 activos + 100 inactivos): {auc_train:.3f}")

# Evaluacion en TEST (nunca visto durante la calibracion) -- esto es lo que
# realmente dice si el metodo generaliza
test_names = list(g_act_test.keys()) + list(g_inact_test.keys())
test_labels = [1]*len(g_act_test) + [0]*len(g_inact_test)
all_test_graphs = {**g_act_test, **g_inact_test}

def score_test(name, theta):
    sims = [ged_similarity(all_test_graphs[name], g_act_train[a], theta)
            for a in g_act_train]
    return np.mean(sims)

test_scores = [score_test(n, theta_opt) for n in test_names]
auc_test = roc_auc_score(test_labels, test_scores)
print(f"AUC en TEST (60 activos + 200 inactivos, no vistos en calibracion): {auc_test:.3f}")

fpr, tpr, _ = roc_curve(test_labels, test_scores)
plt.figure(figsize=(5,5))
plt.plot(fpr, tpr, color="#378ADD", lw=2, label=f"ErG+GED (AUC test={auc_test:.3f})")
plt.plot([0,1],[0,1], color="#b4b2a9", lw=1, linestyle="--", label="azar (AUC=0.5)")
plt.xlabel("Tasa de falsos positivos"); plt.ylabel("Tasa de verdaderos positivos")
plt.title("Curva ROC en TEST -- PubChem AID 2326 (NS1 influenza)")
plt.legend(loc="lower right"); plt.tight_layout()
plt.savefig("pubchem_roc_test.png", dpi=150)
plt.close()

with open("pubchem_theta_calibrado.json", "w") as fh:
    json.dump({"theta": theta_opt, "auc_train": auc_train, "auc_test": auc_test,
               "n_train_activos": len(g_act_train), "n_train_inactivos": len(g_inact_train),
               "n_test_activos": len(g_act_test), "n_test_inactivos": len(g_inact_test)}, fh, indent=2)

pd.DataFrame({"Molecula": test_names, "Etiqueta": test_labels, "Score": test_scores}).to_csv(
    "pubchem_test_scores.csv", index=False)

print("Guardado: pubchem_theta_calibrado.json, pubchem_roc_test.png, pubchem_test_scores.csv")

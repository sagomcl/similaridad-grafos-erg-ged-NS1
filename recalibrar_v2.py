# ============================================================================
# recalibrar_v2.py
#
# Explicacion detallada de este script (que hace cada paso, como se corre y
# que hay que cambiar para extrapolar el metodo a otra proteina):
#     docs/02_explicacion_del_codigo.md
# Explicacion de los conceptos (nodo, grupo farmacoforico, colapso de atomos,
# papel de la geometria 3D, GED, algoritmo hungaro):
#     docs/01_explicacion_del_metodo.md
# ============================================================================
"""
Recalibracion y revalidacion de theta con el GED corregido (ged_mejorado.py).

Que hace (todo con train/test: el AUC que se reporta es siempre sobre moleculas
que NO se usaron para calibrar):
  1. Carga PubChem AID 2326 desde BasePubChem.xlsx (activos e inactivos reales).
  2. Agrupa los activos por similitud (Butina) para ver las familias quimicas.
  3. Calibra y evalua en: (a) activos mezclados al azar, (b) cada uno de los 3 clusters
     mas grandes, (c) los 3 clusters juntos  -> (c) produce theta_final_v2.json.
  4. Imprime una tabla comparando contra los valores anteriores (calculados con el
     emparejamiento defectuoso).

Uso (con el .venv activo, en la carpeta Universidad):
    python recalibrar_v2.py --rapido     # prueba rapida (pocas iteraciones, unos minutos)
    python recalibrar_v2.py              # corrida completa (mismos parametros que antes)

La primera vez lee BasePubChem.xlsx (270k filas, tarda unos minutos) y guarda
AID2326_cache.csv para las siguientes.
Salida: recalibracion_v2_resultados.json, theta_final_v2.json
"""
import os
import sys
import json
import time
import random
import pandas as pd

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

import ged_mejorado
ged_mejorado.parchear()               # ANTES de usar las funciones del pipeline

import pipeline_final_completo as pf
from pipeline_final_completo import (_limpiar_smiles_sal, clusterizar_activos,
                                     calibrar_y_evaluar_train_test)

OUT = os.path.dirname(os.path.abspath(__file__))
RAPIDO = "--rapido" in sys.argv
MAXITER, POPSIZE = (5, 5) if RAPIDO else (15, 8)
N_INACT = 60
SEED = 42

# ---------------------------------------------------------------- datos
cache = os.path.join(OUT, "AID2326_cache.csv")
if not os.path.exists(cache):
    print("Leyendo BasePubChem.xlsx (solo la primera vez, puede tardar varios minutos)...")
    df = pd.read_excel(os.path.join(OUT, "BasePubChem.xlsx"),
                       usecols=["PUBCHEM_CID", "PUBCHEM_EXT_DATASOURCE_SMILES", "PUBCHEM_ACTIVITY_OUTCOME"])
    df = df.dropna()
    df = df[df["PUBCHEM_ACTIVITY_OUTCOME"].isin(["Active", "Inactive"])]
    df.to_csv(cache, index=False)
df = pd.read_csv(cache)
df["cid"] = df["PUBCHEM_CID"].astype(float).astype(int).astype(str)

act_df = df[df["PUBCHEM_ACTIVITY_OUTCOME"] == "Active"]
inact_df = df[df["PUBCHEM_ACTIVITY_OUTCOME"] == "Inactive"].sample(n=5000, random_state=SEED)
print(f"AID 2326: {len(act_df)} activos, {len(df) - len(act_df)} inactivos "
      f"(se usa una muestra de {len(inact_df)} inactivos)")


def limpiar(sub):
    d = {}
    for cid, smi in zip(sub["cid"], sub["PUBCHEM_EXT_DATASOURCE_SMILES"]):
        s = _limpiar_smiles_sal(smi)
        if s:
            d[cid] = s
    return d


activos = limpiar(act_df)
inactivos = limpiar(inact_df)
clusters = clusterizar_activos(activos, cutoff_distancia=0.6)
print(f"(referencia anterior: 737 clusters y 588 singletons)")

rng = random.Random(SEED)
inact_ids = list(inactivos)


# ---------------------------------------------------------------- calibracion con progreso
import numpy as np
from scipy.optimize import differential_evolution
from sklearn.metrics import roc_auc_score


def calibrar_con_progreso(activos_smiles, decoys_smiles, frac_train=0.6,
                          maxiter=15, popsize=8, seed=42):
    """Igual que pipeline_final_completo.calibrar_y_evaluar_train_test (mismo split
    y misma funcion objetivo), pero imprime el avance de la optimizacion.
    Diferencia: polish=False (el AUC es una funcion escalonada, el pulido final
    no aporta y gasta evaluaciones extra)."""
    rng = np.random.RandomState(seed)
    act_keys = list(activos_smiles.keys()); rng.shuffle(act_keys)
    dec_keys = list(decoys_smiles.keys()); rng.shuffle(dec_keys)
    n_act_train = max(2, int(len(act_keys) * frac_train))
    n_dec_train = max(2, int(len(dec_keys) * frac_train))
    act_train, act_test = act_keys[:n_act_train], act_keys[n_act_train:]
    dec_train, dec_test = dec_keys[:n_dec_train], dec_keys[n_dec_train:]

    def build(keys, smiles_dict):
        graphs = {}
        for k in keys:
            try:
                mol3d, _ = pf.smiles_to_3d_mol(smiles_dict[k])
                graphs[k] = pf.build_erg(mol3d)
            except Exception:
                continue
        return graphs

    print("    construyendo grafos ErG...")
    g_act_train, g_act_test = build(act_train, activos_smiles), build(act_test, activos_smiles)
    g_dec_train, g_dec_test = build(dec_train, decoys_smiles), build(dec_test, decoys_smiles)
    print(f"    train: {len(g_act_train)} activos + {len(g_dec_train)} decoys | "
          f"test: {len(g_act_test)} activos + {len(g_dec_test)} decoys")

    keys_theta = list(pf.DEFAULT_THETA.keys())
    train_names = list(g_act_train) + list(g_dec_train)
    train_labels = [1] * len(g_act_train) + [0] * len(g_dec_train)
    all_train = {**g_act_train, **g_dec_train}
    bounds = [(0.05, 2.5), (0.05, 2.5), (0.3, 3.0), (0.2, 3.0), (0.02, 1.5), (0.1, 2.0)]
    total_ev = popsize * len(bounds) * (maxiter + 1)       # maximo de evaluaciones
    paso = max(1, total_ev // 20)                          # un aviso cada ~5%
    estado = {"n": 0, "mejor": 0.0, "t0": time.time()}

    def score(name, theta_vec):
        theta = dict(zip(keys_theta, theta_vec))
        sims = [pf.ged_similarity(all_train[name], g_act_train[a], theta)
                for a in g_act_train if a != name]
        return np.mean(sims) if sims else 0.0

    def objetivo(theta_vec):
        scores = [score(n, theta_vec) for n in train_names]
        try:
            auc = roc_auc_score(train_labels, scores)
        except ValueError:
            auc = 0.0
        estado["n"] += 1
        estado["mejor"] = max(estado["mejor"], auc)
        n, dt = estado["n"], time.time() - estado["t0"]
        if n == 1:
            print(f"    1 evaluacion tarda {dt:.1f} s -> como maximo ~{total_ev * dt / 60:.0f} min "
                  f"({total_ev} evaluaciones; puede terminar antes si converge)")
        elif n % paso == 0:
            falta = (total_ev - n) * dt / n / 60
            print(f"    ... {n}/{total_ev} evaluaciones ({100 * n / total_ev:.0f}%) | "
                  f"mejor AUC train hasta ahora: {estado['mejor']:.3f} | "
                  f"{dt / 60:.1f} min transcurridos, ~{falta:.0f} min restantes (maximo)")
        return -auc

    print(f"    optimizando los 6 pesos (evolucion diferencial)...")
    res = differential_evolution(objetivo, bounds, seed=seed, maxiter=maxiter,
                                 popsize=popsize, tol=1e-3, polish=False)
    theta_final = dict(zip(keys_theta, res.x))
    auc_train = -res.fun
    print(f"    optimizacion terminada tras {estado['n']} evaluaciones. AUC train = {auc_train:.3f}")

    test_names = list(g_act_test) + list(g_dec_test)
    test_labels = [1] * len(g_act_test) + [0] * len(g_dec_test)
    all_test = {**g_act_test, **g_dec_test}
    if len(set(test_labels)) < 2:
        auc_test = float("nan")
    else:
        print("    evaluando en moleculas de prueba (nunca vistas)...")
        test_scores = [np.mean([pf.ged_similarity(all_test[n], g_act_train[a], theta_final)
                                for a in g_act_train]) for n in test_names]
        auc_test = roc_auc_score(test_labels, test_scores)
    print(f"    AUC train = {auc_train:.3f} | AUC test = {auc_test:.3f}")

    # Linea base: Morgan+Tanimoto (similitud media a los activos de entrenamiento),
    # sobre exactamente las mismas moleculas de prueba.
    auc_tan = float("nan")
    if len(set(test_labels)) == 2:
        from rdkit import Chem, DataStructs
        from rdkit.Chem import AllChem
        todos = {**activos_smiles, **decoys_smiles}
        fp = lambda s: AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), radius=2, nBits=1024)
        fps_train = [fp(todos[k]) for k in g_act_train]
        tan_scores = [float(np.mean(DataStructs.BulkTanimotoSimilarity(fp(todos[n]), fps_train)))
                      for n in test_names]
        auc_tan = roc_auc_score(test_labels, tan_scores)
    print(f"    linea base Tanimoto: AUC test = {auc_tan:.3f}  (ErG+GED: {auc_test:.3f})")

    detalle = {"n_train_act": len(g_act_train), "n_train_dec": len(g_dec_train),
               "n_test_act": len(g_act_test), "n_test_dec": len(g_dec_test),
               "auc_test_tanimoto": auc_tan}
    return theta_final, auc_train, auc_test, detalle


def muestra_inactivos(n, semilla):
    r = random.Random(semilla)
    return {k: inactivos[k] for k in r.sample(inact_ids, n)}


# ---------------------------------------------------------------- experimentos
resultados = {}


def correr(nombre, act, inact, frac=0.6):
    print(f"\n>>> {nombre}: {len(act)} activos, {len(inact)} inactivos "
          f"(maxiter={MAXITER}, popsize={POPSIZE})")
    t0 = time.time()
    theta, auc_tr, auc_te, det = calibrar_con_progreso(
        act, inact, frac_train=frac, maxiter=MAXITER, popsize=POPSIZE, seed=SEED)
    resultados[nombre] = {"theta": theta, "auc_train": auc_tr, "auc_test": auc_te, **det,
                          "minutos": round((time.time() - t0) / 60, 1)}
    print(f"    listo en {resultados[nombre]['minutos']} min")
    return theta, auc_tr, auc_te, det


# (a) activos mezclados al azar
ids_mezcla = random.Random(SEED).sample(list(activos), 35)
correr("mezclado_azar", {k: activos[k] for k in ids_mezcla}, muestra_inactivos(N_INACT, 100))

# (b) cada uno de los 3 clusters mas grandes
for k in range(3):
    cl = clusters[k]
    correr(f"cluster_{k + 1}", {c: activos[c] for c in cl}, muestra_inactivos(N_INACT, 100 + k))

# (c) los 3 clusters juntos -> theta final
union = {c: activos[c] for cl in clusters[:3] for c in cl}
theta_f, auc_tr, auc_te, det = correr("3_clusters_juntos", union, muestra_inactivos(N_INACT, 200))

with open(os.path.join(OUT, "recalibracion_v2_resultados.json"), "w") as fh:
    json.dump(resultados, fh, indent=2)
with open(os.path.join(OUT, "theta_final_v2.json"), "w") as fh:
    json.dump({"theta": theta_f, "auc_train": auc_tr, "auc_test": auc_te,
               "n_train_act": det["n_train_act"], "n_test_act": det["n_test_act"],
               "nota": "calibrado con ged_mejorado (emparejamiento que considera aristas)"}, fh, indent=2)

# ---------------------------------------------------------------- resumen
anterior = {"mezclado_azar": 0.49, "cluster_1": 0.962, "cluster_2": 0.671,
            "cluster_3": 0.733, "3_clusters_juntos": 0.65}
print("\n" + "=" * 70)
print("RESUMEN (AUC en moleculas de PRUEBA, nunca vistas en la calibracion)")
print("=" * 70)
print(f"{'experimento':<22}{'ErG+GED nuevo':>15}{'Tanimoto':>10}{'ErG anterior':>14}{'n test act':>12}")
for k, v in resultados.items():
    print(f"{k:<22}{v['auc_test']:>15.3f}{v.get('auc_test_tanimoto', float('nan')):>10.3f}"
          f"{anterior.get(k, float('nan')):>14.3f}{v['n_test_act']:>12}")
print("\n(anterior = calculado con el emparejamiento defectuoso)")
print("\ntheta anterior vs nuevo (3 clusters juntos):")
viejo = json.load(open(os.path.join(OUT, "theta_final.json")))["theta"]
for k in theta_f:
    print(f"  {k:<24}{viejo[k]:>8.3f} -> {theta_f[k]:>8.3f}")
print("\nGuardado: recalibracion_v2_resultados.json, theta_final_v2.json")

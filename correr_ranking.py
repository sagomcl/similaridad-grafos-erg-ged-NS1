# ============================================================================
# correr_ranking.py
#
# Explicacion detallada de este script (que hace cada paso, como se corre y
# que hay que cambiar para extrapolar el metodo a otra proteina):
#     docs/02_explicacion_del_codigo.md
# Explicacion de los conceptos (nodo, grupo farmacoforico, colapso de atomos,
# papel de la geometria 3D, GED, algoritmo hungaro):
#     docs/01_explicacion_del_metodo.md
# ============================================================================
"""
Corre el ranking completo de Moleculas.xlsx (~33k) contra A9 y A22 usando el
theta validado (theta_final.json) y reutilizando las funciones de
pipeline_final_completo.py. Al terminar imprime y guarda el top-20 de cada uno.

Uso (desde VS Code o terminal, dentro de la carpeta Universidad):
    python correr_ranking.py

Requisitos: rdkit, scipy, scikit-learn, networkx, pandas, numpy, matplotlib, openpyxl
Tarda ~30-45 min. Si se interrumpe, vuelve a correrlo: retoma desde el checkpoint.
Para empezar de cero, borra ranking_33k_checkpoint.json y ranking_33k_completo.csv.
"""
import os
import json
import pandas as pd

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")  # silencia los DEPRECATION WARNING de RDKit

import sys
USAR_THETA_NUEVO = "--theta-nuevo" in sys.argv   # usa theta_final_v2.json (de recalibrar_v2.py)
USAR_V2 = "--v2" in sys.argv or USAR_THETA_NUEVO   # python correr_ranking.py --v2
if USAR_V2:
    import ged_mejorado
    ged_mejorado.parchear()           # debe ir ANTES de importar las funciones del pipeline
SUF = "_v2b" if USAR_THETA_NUEVO else ("_v2" if USAR_V2 else "")
THETA_FILE = "theta_final_v2.json" if USAR_THETA_NUEVO else "theta_final.json"

import time
from collections import defaultdict
import pipeline_final_completo as pf
from pipeline_final_completo import rankear_excel_en_bloques, INFLUENZA_ACTIVOS

# --- cronometro por etapa (envuelve las funciones del pipeline, sin modificarlo) ---
T, N = defaultdict(float), defaultdict(int)

def cronometrar(nombre, f):
    def envuelta(*a, **k):
        t0 = time.perf_counter()
        try:
            return f(*a, **k)
        finally:
            T[nombre] += time.perf_counter() - t0
            N[nombre] += 1
    return envuelta

pf.smiles_to_3d_mol = cronometrar("embedding_3D", pf.smiles_to_3d_mol)
pf.build_erg = cronometrar("construir_ErG", pf.build_erg)
pf.ged_similarity = cronometrar("GED", pf.ged_similarity)

BATCH = 1500   # moleculas por bloque (cada bloque vuelve a leer el Excel)

OUT = os.path.dirname(os.path.abspath(__file__))
EXCEL = os.path.join(OUT, "Moleculas.xlsx")
OUT_CSV = os.path.join(OUT, f"ranking_33k_completo{SUF}.csv")
CKPT = os.path.join(OUT, f"ranking_33k_checkpoint{SUF}.json")
TOP_N = 20

with open(os.path.join(OUT, THETA_FILE)) as fh:
    theta = json.load(fh)["theta"]
print(f"theta usado ({THETA_FILE}):", theta)

# 1) Procesar todas las moleculas, por bloques, hasta llegar al 100%
TIEMPOS_CSV = os.path.join(OUT, f"tiempos_por_bloque{SUF}.csv")
idx_prev = json.load(open(CKPT))["next_idx"] if os.path.exists(CKPT) else 0
n_refs = len(INFLUENZA_ACTIVOS)
t_total = time.perf_counter()
while True:
    T.clear(); N.clear()
    t0 = time.perf_counter()
    idx, total = rankear_excel_en_bloques(
        EXCEL, INFLUENZA_ACTIVOS, theta,
        out_csv=OUT_CSV, checkpoint_json=CKPT, batch_size=BATCH)
    dt = time.perf_counter() - t0
    n_mol = idx - idx_prev
    if n_mol > 0:
        ms = lambda k: 1000 * T[k] / n_mol
        resto_min = (total - idx) * dt / n_mol / 60
        print(f"  [tiempos] {n_mol} moleculas en {dt:.1f} s -> {1000*dt/n_mol:.1f} ms/molecula (total, incluye lectura del Excel y Tanimoto)")
        print(f"            embedding 3D: {ms('embedding_3D'):.1f} ms/mol | construir ErG: {ms('construir_ErG'):.1f} ms/mol | "
              f"GED ({n_refs} referencias): {ms('GED'):.1f} ms/mol = {1000*T['GED']/max(N['GED'],1):.2f} ms por comparacion")
        print(f"            tiempo restante estimado: {resto_min:.1f} min")
        nuevo = not os.path.exists(TIEMPOS_CSV)
        with open(TIEMPOS_CSV, "a") as fh:
            if nuevo:
                fh.write("fila_fin,moleculas,segundos,ms_por_molecula,ms_embedding,ms_ErG,ms_GED,ms_por_comparacion_GED\n")
            fh.write(f"{idx},{n_mol},{dt:.2f},{1000*dt/n_mol:.2f},{ms('embedding_3D'):.2f},{ms('construir_ErG'):.2f},"
                     f"{ms('GED'):.2f},{1000*T['GED']/max(N['GED'],1):.3f}\n")
    idx_prev = idx
    if idx >= total:
        break
print(f"\nTiempo total de esta sesion: {(time.perf_counter() - t_total)/60:.1f} min")

# 2) Top-20 por activo (por separado, nunca promediado)
df = pd.read_csv(OUT_CSV)
meta = pd.read_excel(EXCEL, usecols=["ChEMBL ID", "Name", "Smiles"]).rename(columns={"ChEMBL ID": "ID"})
df = df.merge(meta, on="ID", how="left")

tops = {}
for ref in INFLUENZA_ACTIVOS:
    col = f"ErG+GED__{ref}"
    tcol = f"Morgan+Tanimoto__{ref}"
    top = df.sort_values(col, ascending=False).head(TOP_N)[["ID", "Name", col, tcol, "Smiles"]]
    top.to_csv(os.path.join(OUT, f"top{TOP_N}_{ref}{SUF}.csv"), index=False)
    tops[ref] = top
    print(f"\n=== Top-{TOP_N} ErG+GED vs {ref} ===")
    print(top[["ID", "Name", col, tcol]].to_string(index=False))

refs = list(INFLUENZA_ACTIVOS)
if len(refs) >= 2:
    comunes = set(tops[refs[0]]["ID"]) & set(tops[refs[1]]["ID"])
    print(f"\nMoleculas en comun entre top-{TOP_N} de {refs[0]} y {refs[1]}: {len(comunes)}")
    if comunes:
        print(sorted(comunes))

print(f"\nGuardado: ranking_33k_completo{SUF}.csv, top{TOP_N}_A9{SUF}.csv, top{TOP_N}_A22{SUF}.csv")

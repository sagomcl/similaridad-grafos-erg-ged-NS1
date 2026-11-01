# ============================================================================
# comparar_A9_A22.py
#
# Explicacion detallada de este script (que hace cada paso, como se corre y
# que hay que cambiar para extrapolar el metodo a otra proteina):
#     docs/02_explicacion_del_codigo.md
# Explicacion de los conceptos (nodo, grupo farmacoforico, colapso de atomos,
# papel de la geometria 3D, GED, algoritmo hungaro):
#     docs/01_explicacion_del_metodo.md
# ============================================================================
"""
Comparacion A9 vs A22 y de sus top-20 (ErG+GED vs Tanimoto).

PARTE 1: A9 vs A22 en detalle -> similitud ErG+GED (theta calibrado y por defecto),
         Tanimoto, emparejamiento de nodos farmacoforicos, diferencias en distancias
         entre grupos, subestructura comun maxima (MCS) y propiedades.
PARTE 2: top-20 de A9 vs top-20 de A22 -> matrices cruzadas 20x20 con ErG+GED y
         Tanimoto, similitud dentro de cada lista, similitud de cada lista contra
         ambos activos, y heatmaps.

Requisitos: haber corrido correr_ranking.py (genera top20_A9.csv y top20_A22.csv)
Uso:  python comparar_A9_A22.py
Salida: comparacion_A9_A22_resumen.txt, comparacion_top20_ErG.csv,
        comparacion_top20_Tanimoto.csv, comparacion_heatmaps.png
"""
import os
import json
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment
from scipy.stats import spearmanr

from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, rdFMCS
RDLogger.DisableLog("rdApp.*")

import sys
USAR_V2 = "--v2" in sys.argv          # python comparar_A9_A22.py --v2
if USAR_V2:
    import ged_mejorado
    ged_mejorado.parchear()           # debe ir ANTES de importar ged_similarity
SUF = "_v2" if USAR_V2 else ""

from pipeline_final_completo import (smiles_to_3d_mol, build_erg, ged_similarity,
                                     cost_matrix_from_theta, INFLUENZA_ACTIVOS,
                                     DEFAULT_THETA, IDX)

OUT = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(OUT, "theta_final.json")) as fh:
    THETA = json.load(fh)["theta"]

lines = []
def log(s=""):
    print(s)
    lines.append(str(s))


def prep(smi):
    mol3d, mol2d = smiles_to_3d_mol(smi)
    return {"G": build_erg(mol3d), "mol": mol2d,
            "fp": AllChem.GetMorganFingerprintAsBitVect(mol2d, radius=2, nBits=1024)}


def tanimoto(a, b):
    return DataStructs.TanimotoSimilarity(a["fp"], b["fp"])


def erg_sim(a, b, theta=THETA):
    return ged_similarity(a["G"], b["G"], theta)


def nodo_str(P, i):
    d = P["G"].nodes[i]
    simbolos = [P["mol"].GetAtomWithIdx(k).GetSymbol() for k in d["atom_ids"]]
    return f"{d['label']}[{''.join(simbolos)}]"


def detalle_ged(P1, P2, n1, n2, theta):
    """Replica ged_similarity pero imprimiendo el emparejamiento y las diferencias."""
    G1, G2 = P1["G"], P2["G"]
    C_SUB = cost_matrix_from_theta(theta)
    a, b = G1.number_of_nodes(), G2.number_of_nodes()
    N = max(a, b)
    Cn = np.zeros((N, N))
    for i in range(a):
        for j in range(b):
            Cn[i, j] = C_SUB[IDX[G1.nodes[i]["label"]], IDX[G2.nodes[j]["label"]]]
    for i in range(a, N):
        Cn[i, :b] = theta["node_indel"]
    for i in range(a):
        Cn[i, b:] = theta["node_indel"]
    if USAR_V2:
        mapping = ged_mejorado.ged_detallado(G1, G2, theta)[1]
    else:
        r, c = linear_sum_assignment(Cn)
        mapping = {int(i): int(j) for i, j in zip(r, c) if i < a and j < b}

    log(f"  Nodos {n1}: {a}   |   Nodos {n2}: {b}   |   Aristas {n1}: {G1.number_of_edges()}   Aristas {n2}: {G2.number_of_edges()}")
    log("  Emparejamiento de nodos (grupos farmacoforicos):")
    for i, j in mapping.items():
        costo = Cn[i, j]
        marca = "igual tipo" if costo == 0 else f"tipo distinto, costo {costo:.2f}"
        log(f"    {n1}:{nodo_str(P1, i):<22} <-> {n2}:{nodo_str(P2, j):<22} ({marca})")
    for i in range(a):
        if i not in mapping:
            log(f"    {n1}:{nodo_str(P1, i)} sin pareja en {n2}")
    usados = set(mapping.values())
    for j in range(b):
        if j not in usados:
            log(f"    {n2}:{nodo_str(P2, j)} sin pareja en {n1}")

    log("  Diferencias de distancia (en enlaces) entre pares de grupos emparejados:")
    usadas, dif_total, n_dif = set(), 0, 0
    for (u, v, d) in G1.edges(data=True):
        if u in mapping and v in mapping:
            u2, v2 = mapping[u], mapping[v]
            if G2.has_edge(u2, v2):
                d2 = int(G2.edges[u2, v2]["dist"])
                usadas.add(tuple(sorted((u2, v2))))
                dif = abs(int(d["dist"]) - d2)
                if dif > 0:
                    n_dif += 1
                    dif_total += dif
                    log(f"    {nodo_str(P1, u)}-{nodo_str(P1, v)}: {n1}={int(d['dist'])} enlaces vs {n2}={d2}  (dif {dif})")
            else:
                log(f"    {nodo_str(P1, u)}-{nodo_str(P1, v)}: existe en {n1} ({int(d['dist'])}) pero no en {n2}")
        else:
            log(f"    arista {u}-{v} de {n1} sin equivalente")
    for (x, y, d) in G2.edges(data=True):
        if tuple(sorted((x, y))) not in usadas:
            log(f"    arista {nodo_str(P2, x)}-{nodo_str(P2, y)} solo existe en {n2}")
    log(f"  -> {n_dif} pares con distancia distinta (suma de diferencias = {dif_total} enlaces)")


# ======================================================================
# PARTE 1 — A9 vs A22 en detalle
# ======================================================================
log("=" * 78)
log("PARTE 1 — A9 vs A22")
log("=" * 78)
P9 = prep(INFLUENZA_ACTIVOS["A9"])
P22 = prep(INFLUENZA_ACTIVOS["A22"])

log(f"\nSimilitud ErG+GED (theta calibrado): {erg_sim(P9, P22):.3f}")
log(f"Similitud ErG+GED (theta por defecto, sin calibrar): {erg_sim(P9, P22, DEFAULT_THETA):.3f}")
log(f"Similitud Morgan+Tanimoto: {tanimoto(P9, P22):.3f}")

log("\n--- Como se emparejan los grafos ErG (theta calibrado) ---")
detalle_ged(P9, P22, "A9", "A22", THETA)

log("\n--- Subestructura comun maxima (MCS) ---")
res = rdFMCS.FindMCS([P9["mol"], P22["mol"]], timeout=10)
patt = Chem.MolFromSmarts(res.smartsString)
n9, n22 = P9["mol"].GetNumHeavyAtoms(), P22["mol"].GetNumHeavyAtoms()
log(f"  Atomos en comun: {res.numAtoms} de {n9} (A9) y {n22} (A22)  |  SMARTS: {res.smartsString}")
for nombre, P in (("A9", P9), ("A22", P22)):
    match = set(P["mol"].GetSubstructMatch(patt))
    fuera = [f"{P['mol'].GetAtomWithIdx(k).GetSymbol()}{k}" for k in range(P["mol"].GetNumHeavyAtoms()) if k not in match]
    log(f"  Atomos de {nombre} fuera de la subestructura comun: {', '.join(fuera) if fuera else 'ninguno'}")

log("\n--- Propiedades fisicoquimicas ---")
props = {}
for nombre, P in (("A9", P9), ("A22", P22)):
    m = P["mol"]
    props[nombre] = {"MW": Descriptors.MolWt(m), "logP": Descriptors.MolLogP(m),
                     "HBD": rdMolDescriptors.CalcNumHBD(m), "HBA": rdMolDescriptors.CalcNumHBA(m),
                     "TPSA": rdMolDescriptors.CalcTPSA(m), "Enlaces rotables": rdMolDescriptors.CalcNumRotatableBonds(m)}
log(pd.DataFrame(props).round(2).to_string())

# ======================================================================
# PARTE 2 — top-20 de A9 vs top-20 de A22
# ======================================================================
log("\n" + "=" * 78)
log("PARTE 2 — top-20 de A9 vs top-20 de A22")
log("=" * 78)
f9, f22 = os.path.join(OUT, "top20_A9.csv"), os.path.join(OUT, "top20_A22.csv")
if not (os.path.exists(f9) and os.path.exists(f22)):
    log("No encontre top20_A9.csv / top20_A22.csv. Corre primero correr_ranking.py.")
else:
    t9, t22 = pd.read_csv(f9), pd.read_csv(f22)

    def build_list(df):
        ids, preps = [], []
        for _, r in df.iterrows():
            try:
                preps.append(prep(r["Smiles"]))
                ids.append(r["ID"])
            except Exception:
                continue
        return ids, preps

    ids9, L9 = build_list(t9)
    ids22, L22 = build_list(t22)
    log(f"Moleculas procesadas: top A9 = {len(L9)}, top A22 = {len(L22)}")
    comunes = set(ids9) & set(ids22)
    log(f"Moleculas en comun entre ambos top-20: {len(comunes)} {sorted(comunes) if comunes else ''}")

    M_erg = np.array([[erg_sim(a, b) for b in L22] for a in L9])
    M_tan = np.array([[tanimoto(a, b) for b in L22] for a in L9])

    def interna(L, f):
        v = [f(L[i], L[j]) for i in range(len(L)) for j in range(i + 1, len(L))]
        return float(np.mean(v)) if v else float("nan")

    resumen = pd.DataFrame({
        "ErG+GED": {
            "cruzada A9-list x A22-list (media)": M_erg.mean(),
            "cruzada (maximo)": M_erg.max(),
            "interna lista A9 (media)": interna(L9, erg_sim),
            "interna lista A22 (media)": interna(L22, erg_sim),
        },
        "Tanimoto": {
            "cruzada A9-list x A22-list (media)": M_tan.mean(),
            "cruzada (maximo)": M_tan.max(),
            "interna lista A9 (media)": interna(L9, tanimoto),
            "interna lista A22 (media)": interna(L22, tanimoto),
        },
    })
    log("\nSimilitud entre las listas (0 a 1):")
    log(resumen.round(3).to_string())

    ref = pd.DataFrame({
        "ErG+GED": {
            "top A9 vs A9": np.mean([erg_sim(x, P9) for x in L9]),
            "top A9 vs A22": np.mean([erg_sim(x, P22) for x in L9]),
            "top A22 vs A22": np.mean([erg_sim(x, P22) for x in L22]),
            "top A22 vs A9": np.mean([erg_sim(x, P9) for x in L22]),
        },
        "Tanimoto": {
            "top A9 vs A9": np.mean([tanimoto(x, P9) for x in L9]),
            "top A9 vs A22": np.mean([tanimoto(x, P22) for x in L9]),
            "top A22 vs A22": np.mean([tanimoto(x, P22) for x in L22]),
            "top A22 vs A9": np.mean([tanimoto(x, P9) for x in L22]),
        },
    })
    log("\nSimilitud media de cada lista contra A9 y A22:")
    log(ref.round(3).to_string())

    rho, _ = spearmanr(M_erg.flatten(), M_tan.flatten())
    log(f"\nCorrelacion (Spearman) entre ErG+GED y Tanimoto en los 400 pares cruzados: {rho:.3f}")
    i, j = np.unravel_index(np.argmax(M_erg), M_erg.shape)
    log(f"Par cruzado mas parecido segun ErG+GED: {ids9[i]} (top A9) vs {ids22[j]} (top A22) = {M_erg[i, j]:.3f}")
    i, j = np.unravel_index(np.argmax(M_tan), M_tan.shape)
    log(f"Par cruzado mas parecido segun Tanimoto: {ids9[i]} (top A9) vs {ids22[j]} (top A22) = {M_tan[i, j]:.3f}")
    log("\nLectura: si la similitud cruzada es mucho menor que la interna de cada lista, los dos top-20")
    log("son vecindarios quimicos distintos. Si la cruzada es parecida a la interna, las listas no se")
    log("diferencian tanto como sugiere el 0/20 de moleculas en comun.")

    pd.DataFrame(M_erg, index=ids9, columns=ids22).round(3).to_csv(os.path.join(OUT, f"comparacion_top20_ErG{SUF}.csv"))
    pd.DataFrame(M_tan, index=ids9, columns=ids22).round(3).to_csv(os.path.join(OUT, f"comparacion_top20_Tanimoto{SUF}.csv"))

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    for ax, M, titulo in ((axes[0], M_erg, "ErG+GED"), (axes[1], M_tan, "Morgan+Tanimoto")):
        im = ax.imshow(M, vmin=0, vmax=1, cmap="viridis")
        ax.set_title(f"{titulo}: top-20 A9 (filas) vs top-20 A22 (columnas)")
        ax.set_xlabel("top A22"); ax.set_ylabel("top A9")
        fig.colorbar(im, ax=ax, fraction=0.046)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, f"comparacion_heatmaps{SUF}.png"), dpi=150)
    plt.close()

with open(os.path.join(OUT, f"comparacion_A9_A22_resumen{SUF}.txt"), "w") as fh:
    fh.write("\n".join(lines))
log(f"\nGuardado: comparacion_A9_A22_resumen{SUF}.txt, comparacion_top20_ErG{SUF}.csv, comparacion_top20_Tanimoto{SUF}.csv, comparacion_heatmaps{SUF}.png")

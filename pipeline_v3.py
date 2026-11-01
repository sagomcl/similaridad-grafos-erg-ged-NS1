# ============================================================================
# pipeline_v3.py
#
# Explicacion detallada de este script (que hace cada paso, como se corre y
# que hay que cambiar para extrapolar el metodo a otra proteina):
#     docs/02_explicacion_del_codigo.md
# Explicacion de los conceptos (nodo, grupo farmacoforico, colapso de atomos,
# papel de la geometria 3D, GED, algoritmo hungaro):
#     docs/01_explicacion_del_metodo.md
# ============================================================================

import os
import math
import json
from collections import deque

import numpy as np
import networkx as nx
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment, differential_evolution
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score, roc_curve

import rdkit.Chem as Chem
from rdkit.Chem import AllChem, ChemicalFeatures, rdReducedGraphs, DataStructs, Descriptors
from rdkit import RDConfig




fdef_path = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
feat_factory = ChemicalFeatures.BuildFeatureFactory(fdef_path)

LABEL_MAP = {"Donor": "HBD", "Acceptor": "HBA", "Aromatic": "Aromatico", "Hydrophobe": "Hidrofobo"}
PRIORITY = ["HBD", "HBA", "Aromatico", "Hidrofobo"]
TYPES = ["Aromatico", "Hidrofobo", "HBA", "HBD"]
IDX = {t: i for i, t in enumerate(TYPES)}


def smiles_to_3d_mol(smiles, seed=0xF00D):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES invalido: {smiles}")
    mol2d = mol
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError(f"Fallo el embedding 3D para: {smiles}")
    AllChem.UFFOptimizeMolecule(mol)
    return mol, mol2d


def build_erg(mol):
    feats = feat_factory.GetFeaturesForMol(mol)
    atom_type = {}
    for f in feats:
        fam = f.GetFamily()
        if fam not in LABEL_MAP:
            continue
        label = LABEL_MAP[fam]
        for a in f.GetAtomIds():
            cur = atom_type.get(a)
            if cur is None or PRIORITY.index(label) < PRIORITY.index(cur):
                atom_type[a] = label

    bond_graph = nx.Graph()
    for b in mol.GetBonds():
        bond_graph.add_edge(b.GetBeginAtomIdx(), b.GetEndAtomIdx())

    conf = mol.GetConformer()
    visited, nodes = set(), []
    for a, label in atom_type.items():
        if a in visited:
            continue
        comp, dq = {a}, deque([a])
        visited.add(a)
        while dq:
            cur = dq.popleft()
            for nb in bond_graph.neighbors(cur):
                if nb in atom_type and atom_type[nb] == label and nb not in visited:
                    visited.add(nb); comp.add(nb); dq.append(nb)
        center = np.array([list(conf.GetAtomPosition(i)) for i in comp]).mean(axis=0)
        nodes.append({"label": label, "atom_ids": sorted(comp), "center": center})

    def shortest_bond_distance(atomsA, atomsB):
        dq = deque([(i, 0) for i in atomsA])
        seen, targets = set(atomsA), set(atomsB)
        while dq:
            v, d = dq.popleft()
            if v in targets:
                return d
            for nb in bond_graph.neighbors(v):
                if nb not in seen:
                    seen.add(nb); dq.append((nb, d + 1))
        return math.inf

    G = nx.Graph()
    for i, n in enumerate(nodes):
        G.add_node(i, **n)
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            d = shortest_bond_distance(nodes[i]["atom_ids"], nodes[j]["atom_ids"])
            if 1 <= d <= 6:
                G.add_edge(i, j, dist=int(d))
    return G


def erg_visualization_data(smiles):
    """Devuelve la estructura de datos paso a paso (atomos, bonds, tipos,
    nodos ErG, aristas con distancia) lista para alimentar una visualizacion
    tipo la que ya tienes en molecula_a_erg_aspirina.html, para CUALQUIER
    molecula (no solo aspirina). Util para la presentacion."""
    mol3d, mol2d = smiles_to_3d_mol(smiles)
    G = build_erg(mol3d)
    AllChem.Compute2DCoords(mol2d)
    conf = mol2d.GetConformer()
    atoms = [{"idx": a.GetIdx(), "symbol": a.GetSymbol(),
              "x": round(conf.GetAtomPosition(a.GetIdx()).x, 3),
              "y": round(conf.GetAtomPosition(a.GetIdx()).y, 3)} for a in mol2d.GetAtoms()]
    bonds = [{"a": b.GetBeginAtomIdx(), "b": b.GetEndAtomIdx()} for b in mol2d.GetBonds()]
    feats = feat_factory.GetFeaturesForMol(mol2d)
    atom_type = {}
    for f in feats:
        fam = f.GetFamily()
        if fam not in LABEL_MAP:
            continue
        label = LABEL_MAP[fam]
        for a in f.GetAtomIds():
            cur = atom_type.get(a)
            if cur is None or PRIORITY.index(label) < PRIORITY.index(cur):
                atom_type[a] = label
    nodes_out = []
    for i, n in enumerate(G.nodes(data=True)):
        _, data = n
        xs = [conf.GetAtomPosition(a).x for a in data["atom_ids"]]
        ys = [conf.GetAtomPosition(a).y for a in data["atom_ids"]]
        nodes_out.append({"label": data["label"], "atoms": data["atom_ids"],
                           "x": round(np.mean(xs), 3), "y": round(np.mean(ys), 3)})
    edges_out = [{"i": i, "j": j, "dist": d["dist"]} for i, j, d in G.edges(data=True)]
    return {"atoms": atoms, "bonds": bonds,
            "atom_type": {str(k): v for k, v in atom_type.items()},
            "nodes": nodes_out, "edges": edges_out}


# =========================================================================
# SECCION 2 — Distancia de edicion de grafos (GED) con algoritmo Hungaro
# =========================================================================
DEFAULT_THETA = {
    "w_aromatico_hidrofobo": 0.50,
    "w_hba_hbd": 0.75,
    "w_polar_apolar": 1.50,
    "node_indel": 1.00,
    "edge_beta": 0.25,
    "edge_indel": 0.50,
}


def cost_matrix_from_theta(theta):
    w_ar_hy = theta["w_aromatico_hidrofobo"]
    w_hba_hbd = theta["w_hba_hbd"]
    w_pa = theta["w_polar_apolar"]
    return np.array([
        [0.0,      w_ar_hy,  w_pa,      w_pa],
        [w_ar_hy,  0.0,      w_pa,      w_pa],
        [w_pa,     w_pa,     0.0,       w_hba_hbd],
        [w_pa,     w_pa,     w_hba_hbd, 0.0],
    ])


def ged_similarity_v1(G1, G2, theta=DEFAULT_THETA, verbose=False):
    """VERSION ANTERIOR (pre-correccion). El húngaro solo mira nodos, por eso
    A9 y A22 (grafos identicos) dan 0,43. Se conserva solo para trazabilidad
    de resultados viejos; usar ged_similarity (ged_mejorado) en todo lo nuevo."""
    C_SUB = cost_matrix_from_theta(theta)
    node_indel, edge_beta, edge_indel = theta["node_indel"], theta["edge_beta"], theta["edge_indel"]

    n1, n2 = G1.number_of_nodes(), G2.number_of_nodes()
    N = max(n1, n2)
    Cn = np.zeros((N, N))
    for i in range(n1):
        li = G1.nodes[i]["label"]
        for j in range(n2):
            Cn[i, j] = C_SUB[IDX[li], IDX[G2.nodes[j]["label"]]]
    for i in range(n1, N):
        Cn[i, :n2] = node_indel
    for i in range(n1):
        Cn[i, n2:] = node_indel

    r, c = linear_sum_assignment(Cn)
    cost_nodes = Cn[r, c].sum()
    mapping = {int(i): int(j) for i, j in zip(r, c) if i < n1 and j < n2}

    used, cost_edges = set(), 0.0
    for (u, v, d) in G1.edges(data=True):
        if u in mapping and v in mapping:
            u2, v2 = mapping[u], mapping[v]
            if G2.has_edge(u2, v2):
                cost_edges += edge_beta * abs(int(d["dist"]) - int(G2.edges[u2, v2]["dist"]))
                used.add(tuple(sorted((u2, v2))))
            else:
                cost_edges += edge_indel
        else:
            cost_edges += edge_indel
    for (x, y, _) in G2.edges(data=True):
        if tuple(sorted((x, y))) not in used:
            cost_edges += edge_indel

    ged = cost_nodes + cost_edges
    denom = n1 + n2 + G1.number_of_edges() + G2.number_of_edges()
    n_ged = ged / denom if denom > 0 else 0.0
    similarity = math.exp(-n_ged)

    if verbose:
        print(f"  emparejamiento G1->G2: {mapping}")
        print(f"  costo_nodos={cost_nodes:.3f}  costo_aristas={cost_edges:.3f}  "
              f"GED={ged:.3f}  GED_norm={n_ged:.3f}  similitud={similarity:.3f}")
    return similarity


# Version corregida del GED (emparejamiento que considera aristas). Reemplaza a
# ged_similarity_v1 sin necesidad de llamar parchear(): hay una sola fuente.
from ged_mejorado import ged_similarity_v2 as ged_similarity  # noqa: E402


# =========================================================================
# SECCION 3 — Moleculas de referencia (verificadas, con fuente citada)
# =========================================================================

# --- 3a. FASE A (ACTUAL): activos de NS1 INFLUENZA = A9 y A22 ---
# Fuente: PDB reales de Salome, parseados con RDKit -> 2-hidroxi-N-fenetil-
# benzamida (A9, C15H15NO2, MW 241.3) y 2-amino-N-fenetilbenzamida
# (A22, C15H16N2O, MW 240.3). Confirmados por Kleinpeter et al. 2018 (JBC)
# como inhibidores del bolsillo CPSF30 de NS1 de INFLUENZA (no de dengue).
# ROL: sistema de control positivo para validar si ErG+GED funciona, antes
# de aplicarlo a dengue (FASE B).
INFLUENZA_ACTIVOS = {
    "A9":  "O=C(NCCc1ccccc1)c1ccccc1O",
    "A22": "Nc1ccccc1C(=O)NCCc1ccccc1",
}

# --- 3b. Benchmark adicional de influenza CON IC50 real (opcional, mas rico) ---
# Fuente: Basu et al. 2009, J. Virol., DOI 10.1128/JVI.01805-08 — Tabla 1.
# SMILES verificados via el resolvedor NCI Cactus (nombre NSC -> estructura).
INFLUENZA_BENCHMARK_NSC = {
    "NSC125044": "Oc1ccccc1C(=O)NN=Cc2ccc3CCCc3c2",
    "NSC128164": "CN(C)Cc1cc(ccc1c2ccc(cc2CN(C)C)C(C)(C)C)C(C)(C)C",
    "NSC109834": "CC1(C)N=C(N)N=C(N)N1c2cccc(OCc3cccc(N)c3)c2",
    "NSC95676":  "CC1CCCC2(CC[N+](C)(CCC[N+](C)(C)C)CC2)C1",
}
INFLUENZA_IC50_UM = {
    "NSC128164": np.mean([14, 5, 13]),
    "NSC109834": np.mean([10, 2, 12]),
    "NSC95676":  np.mean([12, 20, 19]),
    "NSC125044": np.mean([11, 7, 12]),
}

# --- 3c. FASE B (futura, solo si FASE A da resultados razonables): dengue ---
# Fuente: Ansari, Alisha, Baig et al. 2025/2026, Molecular Diversity,
# DOI 10.1007/s11030-025-11447-5 — validados in vitro contra NS1 dengue.
DENGUE_ACTIVOS = {
    "Berberina":        "O1c2c(OC1)cc5c(c2)c4cc3ccc(OC)c(OC)c3c[n+]4CC5",
    "Acido_betulinico": "CC(=C)C1CCC2(C1C3CCC4C5(CCC(C(C5CCC4(C3(CC2)C)C)(C)C)O)C)C(=O)O",
    "Acido_ursolico":   "CC1CCC2(CCC3(C(=CCC4C3(CCC5C4(CCC(C5(C)C)O)C)C)C2C1C)C)C(=O)O",
}


# =========================================================================
# SECCION 4 — Metodos de comparacion (para saber si ErG+GED sirve o no)
# =========================================================================
def tanimoto_bits(fp1, fp2):
    return DataStructs.TanimotoSimilarity(fp1, fp2)


def hamming_similarity_bits(fp1, fp2):
    a1 = np.zeros((len(fp1),), dtype=int); DataStructs.ConvertToNumpyArray(fp1, a1)
    a2 = np.zeros((len(fp2),), dtype=int); DataStructs.ConvertToNumpyArray(fp2, a2)
    return 1 - np.sum(a1 != a2) / len(a1)


def erg_fp_cosine(fp1, fp2):
    a, b = np.array(fp1), np.array(fp2)
    num = np.dot(a, b)
    den = np.dot(a, a) + np.dot(b, b) - num
    return num / den if den > 0 else 0.0


def _prep(smi):
    mol3d, mol2d = smiles_to_3d_mol(smi)
    return {
        "erg_graph": build_erg(mol3d),
        "erg_fp": rdReducedGraphs.GetErGFingerprint(mol2d),
        "morgan_fp": AllChem.GetMorganFingerprintAsBitVect(mol2d, radius=2, nBits=1024),
        "mw": Descriptors.MolWt(mol2d),
    }


def compute_all_similarities(smiles_dict, referencia_smiles_dict, theta=DEFAULT_THETA,
                              por_activo=False):
    """Para cada molecula en smiles_dict, calcula su similitud contra el set
    de referencia con los 4 metodos.
    Si por_activo=True, devuelve una columna de similitud POR CADA referencia
    por separado (recomendado: promediar sobre activos quimicamente distintos
    diluye la senal real, como se documento en el hallazgo del 2026-08-16)."""
    ref_data = {name: _prep(smi) for name, smi in referencia_smiles_dict.items()}
    rows = []
    for name, smi in smiles_dict.items():
        d = _prep(smi)
        row = {"Molecula": name}
        s_erg_ged, s_ergfp, s_morgan_tan, s_morgan_ham = [], [], [], []
        for ref_name, r in ref_data.items():
            if ref_name == name:
                continue
            v_ged = ged_similarity(d["erg_graph"], r["erg_graph"], theta)
            v_ergfp = erg_fp_cosine(d["erg_fp"], r["erg_fp"])
            v_tan = tanimoto_bits(d["morgan_fp"], r["morgan_fp"])
            v_ham = hamming_similarity_bits(d["morgan_fp"], r["morgan_fp"])
            s_erg_ged.append(v_ged); s_ergfp.append(v_ergfp)
            s_morgan_tan.append(v_tan); s_morgan_ham.append(v_ham)
            if por_activo:
                row[f"ErG+GED__{ref_name}"] = v_ged
                row[f"Morgan+Tanimoto__{ref_name}"] = v_tan
        row["ErG+GED"] = np.mean(s_erg_ged)
        row["ErG-fp+coseno"] = np.mean(s_ergfp)
        row["Morgan+Tanimoto"] = np.mean(s_morgan_tan)
        row["Morgan+Hamming"] = np.mean(s_morgan_ham)
        rows.append(row)
    return pd.DataFrame(rows)


# =========================================================================
# SECCION 5 — Decoys reales + calibracion de pesos CON curva ROC guardada
# =========================================================================
# Que hace: en vez de usar farmacos comunes (aspirina, cafeina) como decoys
# "faciles", construye decoys REALES sacados de tu propio dataset de
# candidatos (Moleculas.xlsx / la muestra de 500), filtrados por:
#   (a) peso molecular parecido al de los activos (+-100 Da aprox), y
#   (b) similitud de Tanimoto BAJA (<0.08) frente a TODOS los activos
# -- este es el criterio tipo DUD-E: decoys que se parecen en propiedades
# fisicoquimicas generales pero no en estructura real, para que la tarea
# de discriminacion sea justa (no trivial) y no un maquillaje facil.
#
# Luego calibra theta (differential_evolution maximizando AUC) y esta vez
# SI guarda: (1) el JSON de pesos, (2) el CSV con las etiquetas y scores
# usados, (3) una figura con la curva ROC y el AUC, para que puedas citarla
# directamente en la tesis con evidencia reproducible.
# =========================================================================

def construir_decoys_reales(df_candidatos, activos_smiles, smiles_col="Smiles",
                             mw_tolerancia=100, tanimoto_max=0.08, n_decoys=20, seed=42):
    """df_candidatos: DataFrame con al menos la columna smiles_col (y MW si
    la tienes ya calculada en 'Molecular Weight', si no la calcula con RDKit).
    Devuelve dict {id: smiles} de decoys reales, no farmacos de ejemplo."""
    activos_prep = {n: _prep(s) for n, s in activos_smiles.items()}
    mw_activos = [Descriptors.MolWt(Chem.MolFromSmiles(s)) for s in activos_smiles.values()]
    mw_lo, mw_hi = min(mw_activos) - mw_tolerancia, max(mw_activos) + mw_tolerancia

    candidatos = []
    for i, row in df_candidatos.iterrows():
        smi = row[smiles_col]
        if pd.isna(smi):
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        mw = Descriptors.MolWt(mol)
        if not (mw_lo <= mw <= mw_hi):
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
        max_tan = max(tanimoto_bits(fp, activos_prep[n]["morgan_fp"]) for n in activos_prep)
        if max_tan < tanimoto_max:
            candidatos.append((row.get("ChEMBL_ID", f"decoy_{i}"), smi, max_tan))

    candidatos.sort(key=lambda x: x[2])  # los MAS distintos primero
    seleccion = candidatos[:n_decoys]
    print(f"  [decoys reales] candidatos que cumplen filtro: {len(candidatos)}  "
          f"-> seleccionados: {len(seleccion)}  (MW activos: {mw_lo:.0f}-{mw_hi:.0f}, "
          f"Tanimoto<{tanimoto_max})")
    return {cid: smi for cid, smi, _ in seleccion}


def _verificar_no_existe(*rutas):
    """Frena antes de sobrescribir un resultado guardado. Si quieres rehacer un
    resultado, usa un nombre nuevo (por ejemplo con la fecha o la version)."""
    existentes = [r for r in rutas if os.path.exists(r)]
    if existentes:
        raise FileExistsError(
            "No se sobrescribe un resultado existente. Usa otro nombre para: "
            + ", ".join(os.path.basename(r) for r in existentes))


def calibrar_theta_con_roc(activos_smiles, decoys_smiles, out_prefix,
                            bounds=None, maxiter=60, seed=42):
    """Igual que calibrar_theta, pero ademas guarda:
      - {out_prefix}_theta.json   (pesos calibrados)
      - {out_prefix}_roc.png      (curva ROC con el AUC en el titulo)
      - {out_prefix}_scores.csv   (molecula, etiqueta, score -- trazabilidad)
    """
    _verificar_no_existe(f"{out_prefix}_theta.json", f"{out_prefix}_roc.png",
                         f"{out_prefix}_scores.csv")
    if bounds is None:
        bounds = [(0.05, 2.5), (0.05, 2.5), (0.3, 3.0), (0.2, 3.0), (0.02, 1.5), (0.1, 2.0)]
    keys = list(DEFAULT_THETA.keys())

    graphs = {}
    for name, smi in {**activos_smiles, **decoys_smiles}.items():
        mol3d, _ = smiles_to_3d_mol(smi)
        graphs[name] = build_erg(mol3d)

    cal_set = list(activos_smiles.keys()) + list(decoys_smiles.keys())
    cal_labels = [1] * len(activos_smiles) + [0] * len(decoys_smiles)

    def score(name, theta_vec):
        theta = dict(zip(keys, theta_vec))
        sims = [ged_similarity(graphs[name], graphs[a], theta)
                for a in activos_smiles if a != name]
        return np.mean(sims)

    def objetivo(theta_vec):
        scores = [score(n, theta_vec) for n in cal_set]
        try:
            auc = roc_auc_score(cal_labels, scores)
        except ValueError:
            auc = 0.5
        return -auc

    res = differential_evolution(objetivo, bounds, seed=seed, maxiter=maxiter, popsize=15, tol=1e-4)
    theta_opt = dict(zip(keys, res.x))
    auc_opt = -res.fun

    scores_finales = [score(n, res.x) for n in cal_set]
    fpr, tpr, _ = roc_curve(cal_labels, scores_finales)

    plt.figure(figsize=(5, 5))
    plt.plot(fpr, tpr, color="#378ADD", lw=2, label=f"ErG+GED calibrado (AUC={auc_opt:.3f})")
    plt.plot([0, 1], [0, 1], color="#b4b2a9", lw=1, linestyle="--", label="azar (AUC=0.5)")
    plt.xlabel("Tasa de falsos positivos"); plt.ylabel("Tasa de verdaderos positivos")
    plt.title("Curva ROC — calibracion de pesos (activos vs decoys reales)")
    plt.legend(loc="lower right"); plt.tight_layout()
    plt.savefig(f"{out_prefix}_roc.png", dpi=150)
    plt.close()

    with open(f"{out_prefix}_theta.json", "w") as fh:
        json.dump(theta_opt, fh, indent=2)

    pd.DataFrame({"Molecula": cal_set, "Etiqueta": cal_labels, "Score": scores_finales}) \
      .to_csv(f"{out_prefix}_scores.csv", index=False)

    return theta_opt, auc_opt


# =========================================================================
# SECCION 6 — Carga de datos reales (33k) y de una futura base con relacion
#              real a NS1 influenza (PubChem AID 2326 u homologa)
# =========================================================================

def cargar_candidatos_excel(path, smiles_col="Smiles", id_col=None, sheet_name=0):
    """Lee un Excel/CSV con una columna de SMILES y devuelve un dict
    {id: smiles} listo para usar en compute_all_similarities."""
    if str(path).lower().endswith(".csv"):
        df = pd.read_csv(path)
    else:
        df = pd.read_excel(path, sheet_name=sheet_name)
    if smiles_col not in df.columns:
        raise KeyError(f"No se encontro la columna '{smiles_col}'. Columnas disponibles: {list(df.columns)}")
    ids = df[id_col] if id_col and id_col in df.columns else [f"mol_{i}" for i in range(len(df))]
    return {str(i): s for i, s in zip(ids, df[smiles_col]) if pd.notna(s)}


def _limpiar_smiles_sal(smi):
    """Los compuestos de PubChem suelen venir como sal (ej. 'C...N.Cl' con
    un contraion separado por '.'). Nos quedamos con el fragmento mas grande
    (la molecula real), no con el contraion."""
    frags = [f for f in str(smi).split(".") if Chem.MolFromSmiles(f) is not None]
    if not frags:
        return None
    return max(frags, key=len)


def cargar_pubchem_ns1_assay(csv_path):
    """Carga el CSV real descargado de PubChem BioAssay AID 2326 (qHTS
    primario contra la funcion de NS1 de influenza, NIH Chemical Genomics
    Center / Engel lab -- el mismo screening que origino A9/A22/ML303).
    Como conseguirlo:
      1) https://pubchem.ncbi.nlm.nih.gov/bioassay/2326 -> boton "Download"
      2) En el modal, seccion "Deposited Data Table Only" -> CSV -> Save
         (las otras opciones son metadatos del ensayo, no resultados).
      3) Ese CSV YA trae SMILES (columna PUBCHEM_EXT_DATASOURCE_SMILES),
         no hace falta pedirlos aparte via PUG REST.
    Composicion real (confirmada 2026-08-16): 270,315 filas -- 1,073
    Active, 260,684 Inactive, 8,410 Inconclusive, 130 Unspecified.
    Devuelve (activos, inactivos): dos dict {CID: smiles_limpio (sin sales)}."""
    # Antes: pd.read_csv(..., skiprows=[1,2,3,4]) asumia 4 filas de metadatos.
    # AID2326_cache.csv no las trae y perderia 4 moleculas reales. Ahora se
    # filtran por contenido: solo se conservan filas con un desenlace valido.
    df = pd.read_csv(csv_path, low_memory=False)
    df = df[df["PUBCHEM_ACTIVITY_OUTCOME"].isin(["Active", "Inactive", "Inconclusive", "Unspecified"])]
    df = df.dropna(subset=["PUBCHEM_CID", "PUBCHEM_EXT_DATASOURCE_SMILES", "PUBCHEM_ACTIVITY_OUTCOME"])

    activos_df = df[df["PUBCHEM_ACTIVITY_OUTCOME"] == "Active"].copy()
    inactivos_df = df[df["PUBCHEM_ACTIVITY_OUTCOME"] == "Inactive"].copy()

    activos_df["smiles_limpio"] = activos_df["PUBCHEM_EXT_DATASOURCE_SMILES"].apply(_limpiar_smiles_sal)
    inactivos_df["smiles_limpio"] = inactivos_df["PUBCHEM_EXT_DATASOURCE_SMILES"].apply(_limpiar_smiles_sal)
    activos_df = activos_df.dropna(subset=["smiles_limpio"])
    inactivos_df = inactivos_df.dropna(subset=["smiles_limpio"])

    activos = {str(int(r.PUBCHEM_CID)): r.smiles_limpio for r in activos_df.itertuples()}
    inactivos = {str(int(r.PUBCHEM_CID)): r.smiles_limpio for r in inactivos_df.itertuples()}
    print(f"  [PubChem AID 2326] activos reales cargados: {len(activos)}  "
          f"inactivos reales cargados: {len(inactivos)}")
    return activos, inactivos


# =========================================================================
# SECCION 6b — Clustering de activos por similitud estructural (Butina)
# =========================================================================
# Por que hace falta: un cribado fenotipico primario como AID 2326 mezcla
# compuestos que "funcionan" por mecanismos completamente distintos entre
# si (no comparten farmacoforo). Promediar/calibrar contra TODOS esos
# activos mezclados diluye cualquier señal real (mismo problema ya visto
# con berberina vs. terpenoides en dengue). Clusterizar los activos permite
# aislar sub-grupos quimicamente coherentes (mismo chemotipo) para medir,
# de forma justa, si ErG+GED discrimina bien DENTRO de un chemotipo.
# =========================================================================

def clusterizar_activos(activos_smiles_dict, cutoff_distancia=0.6):
    """activos_smiles_dict: {id: smiles}. cutoff_distancia=0.6 equivale
    aproximadamente a Tanimoto>0.4 dentro del mismo cluster (Butina, RDKit).
    Devuelve: lista de clusters (cada uno es una lista de ids), ordenada
    de mayor a menor tamaño."""
    from rdkit.ML.Cluster import Butina
    ids = list(activos_smiles_dict.keys())
    fps = []
    ids_validos = []
    for cid in ids:
        mol = Chem.MolFromSmiles(activos_smiles_dict[cid])
        if mol is None:
            continue
        fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024))
        ids_validos.append(cid)

    n = len(fps)
    dists = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - s for s in sims])

    clusters_idx = Butina.ClusterData(dists, n, cutoff_distancia, isDistData=True)
    clusters_idx = sorted(clusters_idx, key=len, reverse=True)
    clusters_ids = [[ids_validos[i] for i in c] for c in clusters_idx]
    n_singletons = sum(1 for c in clusters_ids if len(c) == 1)
    print(f"  [clustering] {n} activos -> {len(clusters_ids)} clusters "
          f"({n_singletons} sin ningun parecido a otro activo). "
          f"Tamaños top 5: {[len(c) for c in clusters_ids[:5]]}")
    return clusters_ids


# =========================================================================
# SECCION 6c — Calibracion train/test genérica (evita el sobreajuste que
#              vimos con solo A9/A22 -- AUC=1.000 no confiable con n=2)
# =========================================================================
# Diferencia clave con calibrar_theta_con_roc (Seccion 5): esa calibra Y
# evalua sobre el MISMO set (facil de sobreajustar con pocos activos).
# Esta separa un grupo de calibracion (train) de un grupo NUNCA VISTO
# (test) para reportar un AUC de GENERALIZACION real -- la unica forma
# honesta de decir "el metodo funciona" cuando tienes pocos activos.
# =========================================================================

def calibrar_y_evaluar_train_test(activos_smiles, decoys_smiles, frac_train=0.65,
                                   maxiter=15, popsize=8, seed=42,
                                   act_train_ids=None, act_test_ids=None):
    """activos_smiles / decoys_smiles: dict {id: smiles}. Divide cada grupo
    en train/test, calibra theta SOLO con train, evalua AUC en train y en
    test (nunca visto). Devuelve (theta, auc_train, auc_test, detalle).
    Si se pasan act_train_ids / act_test_ids, la particion de activos ya viene
    decidida (por ejemplo por elegir_chemotipos_solo_train) y no se reparte al azar."""
    rng = np.random.RandomState(seed)

    # El orden de los shuffle (activos y luego decoys) se conserva del original
    # para que la particion por defecto siga siendo la misma con la misma semilla.
    if act_train_ids is None:
        act_keys = list(activos_smiles.keys()); rng.shuffle(act_keys)
        n_act_train = max(2, int(len(act_keys) * frac_train))
        act_train, act_test = act_keys[:n_act_train], act_keys[n_act_train:]
    else:
        act_train, act_test = list(act_train_ids), list(act_test_ids)

    dec_keys = list(decoys_smiles.keys()); rng.shuffle(dec_keys)
    n_dec_train = max(2, int(len(dec_keys) * frac_train))
    dec_train, dec_test = dec_keys[:n_dec_train], dec_keys[n_dec_train:]

    def build(keys, smiles_dict):
        graphs = {}
        for k in keys:
            try:
                mol3d, _ = smiles_to_3d_mol(smiles_dict[k])
                graphs[k] = build_erg(mol3d)
            except Exception:
                continue
        return graphs

    g_act_train = build(act_train, activos_smiles)
    g_act_test = build(act_test, activos_smiles)
    g_dec_train = build(dec_train, decoys_smiles)
    g_dec_test = build(dec_test, decoys_smiles)
    print(f"  train: {len(g_act_train)} activos + {len(g_dec_train)} decoys | "
          f"test: {len(g_act_test)} activos + {len(g_dec_test)} decoys")

    keys_theta = list(DEFAULT_THETA.keys())
    train_names = list(g_act_train.keys()) + list(g_dec_train.keys())
    train_labels = [1] * len(g_act_train) + [0] * len(g_dec_train)
    all_train = {**g_act_train, **g_dec_train}

    def score(name, theta_vec, act_graphs):
        theta = dict(zip(keys_theta, theta_vec))
        sims = [ged_similarity(all_train[name], act_graphs[a], theta)
                for a in act_graphs if a != name]
        return np.mean(sims) if sims else 0.0

    def objetivo(theta_vec):
        scores = [score(n, theta_vec, g_act_train) for n in train_names]
        try:
            return -roc_auc_score(train_labels, scores)
        except ValueError:
            return 0.0

    bounds = [(0.05, 2.5), (0.05, 2.5), (0.3, 3.0), (0.2, 3.0), (0.02, 1.5), (0.1, 2.0)]
    res = differential_evolution(objetivo, bounds, seed=seed, maxiter=maxiter, popsize=popsize, tol=1e-3)
    theta_final = dict(zip(keys_theta, res.x))
    auc_train = -res.fun

    test_names = list(g_act_test.keys()) + list(g_dec_test.keys())
    test_labels = [1] * len(g_act_test) + [0] * len(g_dec_test)
    all_test = {**g_act_test, **g_dec_test}

    def score_test(name, theta):
        sims = [ged_similarity(all_test[name], g_act_train[a], theta) for a in g_act_train]
        return np.mean(sims) if sims else 0.0

    if len(set(test_labels)) < 2:
        auc_test = float("nan")
    else:
        test_scores = [score_test(n, theta_final) for n in test_names]
        auc_test = roc_auc_score(test_labels, test_scores)

    print(f"  AUC train={auc_train:.3f}  AUC test (nunca visto)={auc_test:.3f}")
    detalle = {"n_train_act": len(g_act_train), "n_train_dec": len(g_dec_train),
               "n_test_act": len(g_act_test), "n_test_dec": len(g_dec_test)}
    return theta_final, auc_train, auc_test, detalle


def elegir_chemotipos_solo_train(activos_smiles, frac_train=0.65, n_clusters=3,
                                 cutoff_distancia=0.6, seed=42):
    """Corrige la fuga del paso 6 original: antes los 3 clusters mas grandes se
    elegian mirando TODOS los activos, incluidos los de prueba.
    Aqui: 1) se reparten activos en train/test; 2) se clusteriza SOLO train y se
    eligen los n_clusters mas grandes; 3) un activo de prueba entra solo si su
    Tanimoto maximo contra un miembro de esos clusters supera 1 - cutoff (el mismo
    umbral de Butina). Devuelve (ids_train, ids_test)."""
    rng = np.random.RandomState(seed)
    keys = list(activos_smiles.keys()); rng.shuffle(keys)
    n_train = max(2, int(len(keys) * frac_train))
    train_ids, test_ids = keys[:n_train], keys[n_train:]

    clusters = clusterizar_activos({k: activos_smiles[k] for k in train_ids}, cutoff_distancia)
    ids_train = [c for cl in clusters[:n_clusters] for c in cl]

    fp = lambda smi: AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), radius=2, nBits=1024)
    fps_chemo = [fp(activos_smiles[c]) for c in ids_train]
    ids_test = []
    for k in test_ids:
        sims = DataStructs.BulkTanimotoSimilarity(fp(activos_smiles[k]), fps_chemo)
        if max(sims) >= 1 - cutoff_distancia:
            ids_test.append(k)
    print(f"  [chemotipos solo-train] train: {len(ids_train)} activos en {n_clusters} clusters; "
          f"test asignados a esos chemotipos: {len(ids_test)} de {len(test_ids)} de prueba")
    return ids_train, ids_test


# =========================================================================
# SECCION 7 — Ranking sobre las 33k (o una muestra) contra A9/A22
# =========================================================================

def rankear_candidatos(candidatos_smiles, activos_smiles, theta, top_n=20):
    """Corre compute_all_similarities con por_activo=True y devuelve, ademas
    de la tabla completa, el top-N de cada metodo para inspeccion rapida."""
    df = compute_all_similarities(candidatos_smiles, activos_smiles, theta=theta, por_activo=True)
    top_erg = df.sort_values("ErG+GED", ascending=False).head(top_n)
    top_tan = df.sort_values("Morgan+Tanimoto", ascending=False).head(top_n)
    overlap = len(set(top_erg["Molecula"]) & set(top_tan["Molecula"]))
    print(f"  top-{top_n} ErG+GED vs top-{top_n} Morgan+Tanimoto: {overlap}/{top_n} en comun")
    return df, top_erg, top_tan


# =========================================================================
# SECCION 7b — Correr TODA una base grande (33k+) en bloques, con
#              checkpoint en disco para poder pausar y retomar
# =========================================================================
# Por que hace falta: procesar cada molecula implica generar su conformero
# 3D (ETKDG+UFF) y correr el Hungaro -- a ~70-90 ms/molecula, 33 mil
# moleculas toman 30-45 minutos. Si tu entorno corta ejecuciones largas
# (como el sandbox usado para la tesis), esta funcion permite invocarla
# muchas veces seguidas: cada llamada retoma justo donde quedo la anterior
# (guarda el indice en un .json y va agregando filas al CSV de salida).
# =========================================================================

def rankear_excel_en_bloques(path_excel, activos_smiles, theta, out_csv,
                              checkpoint_json, batch_size=1500, smiles_col="Smiles",
                              id_col="ChEMBL ID", sheet_name=0):
    """Procesa un bloque de `batch_size` moleculas del Excel/CSV y AGREGA los
    resultados a out_csv. Llamala repetidamente (por ejemplo en un loop de
    shell, o una vez por cada turno de conversacion) hasta que imprima
    Progreso: 100.0% -- ya para entonces el CSV tiene el ranking completo."""
    ref_data = {}
    for name, smi in activos_smiles.items():
        mol3d, mol2d = smiles_to_3d_mol(smi)
        ref_data[name] = {
            "erg_graph": build_erg(mol3d),
            "morgan_fp": AllChem.GetMorganFingerprintAsBitVect(mol2d, radius=2, nBits=1024),
        }

    err_csv = out_csv.replace(".csv", "_errores.csv")
    if os.path.exists(checkpoint_json):
        with open(checkpoint_json) as fh:
            start_idx = json.load(fh)["next_idx"]
    else:
        start_idx = 0
        # Antes el codigo borraba out_csv aqui. Ahora se niega a empezar encima
        # de un resultado existente sin checkpoint.
        _verificar_no_existe(out_csv, err_csv)

    if str(path_excel).lower().endswith(".csv"):
        df_full = pd.read_csv(path_excel)
    else:
        df_full = pd.read_excel(path_excel, sheet_name=sheet_name)
    df_full = df_full.dropna(subset=[smiles_col]).reset_index(drop=True)
    total = len(df_full)
    end_idx = min(start_idx + batch_size, total)
    print(f"  Procesando filas {start_idx} a {end_idx} de {total}")

    rows, errores = [], []
    for i in range(start_idx, end_idx):
        r = df_full.iloc[i]
        smi = r[smiles_col]
        cid = r.get(id_col, f"mol_{i}") if id_col in df_full.columns else f"mol_{i}"
        try:
            mol3d, mol2d = smiles_to_3d_mol(smi)
            g = build_erg(mol3d)
            morgan_fp = AllChem.GetMorganFingerprintAsBitVect(mol2d, radius=2, nBits=1024)
            row = {"ID": cid}
            for ref_name, rd in ref_data.items():
                row[f"ErG+GED__{ref_name}"] = ged_similarity(g, rd["erg_graph"], theta)
                row[f"Morgan+Tanimoto__{ref_name}"] = tanimoto_bits(morgan_fp, rd["morgan_fp"])
            rows.append(row)
        except Exception as exc:
            # Antes: 'except Exception: continue' descartaba la molecula sin rastro.
            errores.append({"indice": i, "ID": cid, "SMILES": smi, "error": str(exc)})

    # Las moleculas descartadas quedan registradas en un CSV aparte.
    if errores:
        pd.DataFrame(errores).to_csv(err_csv, mode="a", header=not os.path.exists(err_csv), index=False)
    pd.DataFrame(rows).to_csv(out_csv, mode="a", header=(start_idx == 0), index=False)
    with open(checkpoint_json, "w") as fh:
        json.dump({"next_idx": end_idx, "total": total}, fh)
    print(f"  Progreso: {end_idx}/{total} ({100*end_idx/total:.1f}%)  "
          f"| calculadas en el lote: {len(rows)}  | descartadas: {len(errores)} (ver {os.path.basename(err_csv)})")
    return end_idx, total


# =========================================================================
# SECCION 8 — Ejecucion completa (lo que puedes correr hoy mismo)
# =========================================================================
if __name__ == "__main__":
    OUT = os.path.dirname(os.path.abspath(__file__))

    print("=" * 78)
    print("PASO 1 — Datos de visualizacion ErG (para tu presentacion)")
    print("=" * 78)
    viz = erg_visualization_data(INFLUENZA_ACTIVOS["A9"])
    _verificar_no_existe(os.path.join(OUT, "viz_A9.json"))
    with open(os.path.join(OUT, "viz_A9.json"), "w") as fh:
        json.dump(viz, fh, indent=2, default=str)
    print(f"  Grafo ErG de A9: {len(viz['nodes'])} nodos, {len(viz['edges'])} aristas "
          f"(guardado en viz_A9.json, reutiliza la logica de molecula_a_erg_aspirina.html)")

    print("\n" + "=" * 78)
    print("PASO 2 — Decoys reales sacados de tu dataset (500-muestra de Moleculas.xlsx)")
    print("=" * 78)
    ruta_500 = os.path.join(OUT, "merged_500_con_props.csv")
    if not os.path.exists(ruta_500):
        # Este archivo no esta en la carpeta de trabajo. Sin el no hay decoys, y
        # sin decoys no hay calibracion: se detiene con un mensaje claro.
        raise FileNotFoundError(f"Falta {ruta_500}. Genera el archivo de 500 "
                                "moleculas con propiedades antes de correr el paso 2.")
    df_500 = pd.read_csv(ruta_500)
    decoys_influenza = construir_decoys_reales(df_500, INFLUENZA_ACTIVOS, mw_tolerancia=120,
                                                tanimoto_max=0.08, n_decoys=20)

    print("\n" + "=" * 78)
    print("PASO 3 — Calibracion de pesos (A9/A22 vs decoys reales) + curva ROC")
    print("=" * 78)
    theta_flu, auc_flu = calibrar_theta_con_roc(
        INFLUENZA_ACTIVOS, decoys_influenza,
        out_prefix=os.path.join(OUT, "influenza"), maxiter=50)
    print(f"  AUC (A9/A22 vs {len(decoys_influenza)} decoys reales): {auc_flu:.3f}")
    for k, v in theta_flu.items():
        print(f"    {k}: {DEFAULT_THETA[k]:.3f} -> {v:.3f}")
    print("  Guardado: influenza_theta.json, influenza_roc.png, influenza_scores.csv")

    print("\n" + "=" * 78)
    print("PASO 4 — Comparacion de metodos: candidatos (500) vs A9/A22, theta calibrado")
    print("=" * 78)
    candidatos_500 = {r.ChEMBL_ID: r.Smiles for r in df_500.itertuples()}
    df_rank, top_erg, top_tan = rankear_candidatos(candidatos_500, INFLUENZA_ACTIVOS,
                                                    theta=theta_flu, top_n=20)
    _verificar_no_existe(os.path.join(OUT, "ranking_influenza_500.csv"))
    df_rank.to_csv(os.path.join(OUT, "ranking_influenza_500.csv"), index=False)
    print("  Guardado: ranking_influenza_500.csv (incluye columnas por-activo: "
          "ErG+GED__A9, ErG+GED__A22, etc.)")

    print("\n" + "=" * 78)
    print("PASO 5 (referencia opcional) — benchmark IC50 real de Basu et al. 2009")
    print("=" * 78)
    graphs_flu = {n: build_erg(smiles_to_3d_mol(s)[0]) for n, s in INFLUENZA_BENCHMARK_NSC.items()}
    ref = "NSC109834"
    for n in INFLUENZA_BENCHMARK_NSC:
        if n == ref:
            continue
        s = ged_similarity(graphs_flu[n], graphs_flu[ref], theta=theta_flu)
        print(f"  {n} vs {ref}: similitud={s:.3f}  IC50_prom={INFLUENZA_IC50_UM[n]:.1f} uM")

    print("\n" + "=" * 78)
    print("PASO 6 — VALIDACION REAL con PubChem AID 2326 (opcional, requiere el CSV)")
    print("=" * 78)
    # Antes buscaba AID_2326_datatable.csv, que no existe; el archivo de la
    # carpeta es AID2326_cache.csv (mismo formato de PubChem).
    csv_pubchem = os.path.join(OUT, "AID2326_cache.csv")
    if os.path.exists(csv_pubchem):
        activos_reales, inactivos_reales = cargar_pubchem_ns1_assay(csv_pubchem)

        # 6a) chemotipos elegidos SOLO con los activos de entrenamiento (la
        # version anterior los elegia con todos los activos: fuga de datos).
        # Un cribado fenotipico mezcla mecanismos distintos (ver Seccion 6b).
        ids_train, ids_test = elegir_chemotipos_solo_train(activos_reales, frac_train=0.65,
                                                           n_clusters=3, seed=42)

        # 6b) muestra manejable de inactivos reales (miles son demasiados
        # para el embedding 3D uno por uno en tiempo razonable)
        import random
        random.seed(42)
        inact_sample_ids = random.sample(list(inactivos_reales.keys()), min(80, len(inactivos_reales)))
        inactivos_sample = {k: inactivos_reales[k] for k in inact_sample_ids}

        # 6c) calibracion train/test HONESTA con la particion ya decidida
        theta_final, auc_train, auc_test, detalle = calibrar_y_evaluar_train_test(
            activos_reales, inactivos_sample, maxiter=8, popsize=5,
            act_train_ids=ids_train, act_test_ids=ids_test)
        print(f"  RESULTADO VALIDACION: AUC train={auc_train:.3f}  AUC test={auc_test:.3f}  {detalle}")
        ruta_theta = os.path.join(OUT, "theta_validacion_v3.json")
        _verificar_no_existe(ruta_theta)  # no pisar theta_final.json ni otros resultados
        with open(ruta_theta, "w") as fh:
            json.dump({"theta": theta_final, "auc_train": auc_train, "auc_test": auc_test,
                       "detalle": detalle}, fh, indent=2)
    else:
        print(f"  (Salta este paso: no se encontro {csv_pubchem}.)")
        theta_final = theta_flu  # usa el theta de PASO 3 como respaldo

    print("\n" + "=" * 78)
    print("PASO 7 — Aplicacion real: ranking de TODA la base de candidatos (33k+)")
    print("=" * 78)
    # Antes la ruta era data/Moleculas.xlsx; el archivo esta en la raiz de la carpeta.
    excel_candidatos = os.path.join(OUT, "Moleculas.xlsx")
    if os.path.exists(excel_candidatos):
        # Esto puede tardar 30-45 min para 33k moleculas reales. Si tu
        # entorno corta ejecuciones largas, llama esta misma linea varias
        # veces seguidas -- retoma solo donde quedo gracias al checkpoint.
        # Nombre de salida nuevo: no se pisa ranking_33k_completo*.csv.
        idx, total = rankear_excel_en_bloques(
            excel_candidatos, INFLUENZA_ACTIVOS, theta_final,
            out_csv=os.path.join(OUT, "ranking_33k_v3.csv"),
            checkpoint_json=os.path.join(OUT, "ranking_33k_v3_checkpoint.json"),
            batch_size=1500, id_col="ChEMBL ID")
        if idx < total:
            print(f"  Progreso parcial: {idx}/{total}. Vuelve a correr este script para continuar.")
        else:
            print("  Ranking completo. Reportar SIEMPRE por activo por separado "
                  "(ErG+GED__A9, ErG+GED__A22), nunca promediado. Revisa el CSV de "
                  "errores para ver cuantas moleculas no se pudieron calcular.")
    else:
        print(f"  (Salta este paso: no se encontro {excel_candidatos}.)")

    print("\n" + "=" * 78)
    print("PASO 8 — FASE B (dengue, el target real) -- correr solo si FASE A convence")
    print("=" * 78)
    print("  df_dengue = compute_all_similarities(candidatos_500, DENGUE_ACTIVOS, "
          "theta=theta_final, por_activo=True)")
    print("  # Idealmente: repetir Secciones 6b/6c con activos reales de dengue una vez")
    print("  # que existan mas que 3 compuestos confirmados (ver matriz de literatura).")

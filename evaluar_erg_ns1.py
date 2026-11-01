# ============================================================================
# evaluar_erg_ns1.py
#
# Explicacion detallada de este script (que hace cada paso, como se corre y
# que hay que cambiar para extrapolar el metodo a otra proteina):
#     docs/02_explicacion_del_codigo.md
# Explicacion de los conceptos (nodo, grupo farmacoforico, colapso de atomos,
# papel de la geometria 3D, GED, algoritmo hungaro):
#     docs/01_explicacion_del_metodo.md
# ============================================================================
"""
Evaluacion justa de ErG+GED frente a Tanimoto, con los criterios farmacoforicos del
sitio de union de NS1.  Script aparte: no modifica pipeline_final_completo.py.

QUE HACE (todo sobre moleculas de PRUEBA, nunca usadas para entrenar):
  0. Chequea que A9 y A22 cumplan los criterios NS1 (ver ns1_utils.py) y lo reporta.
  1. Evalua con PubChem AID 2326 en tres escenarios:
       - mezclado:      activos al azar (sin familias)
       - intra-cluster: entrenar y probar dentro de una misma familia
       - entre-clusters: entrenar con una familia y probar con OTRAS (aqui Tanimoto
                         no tiene ventaja por construccion: es la prueba justa)
     contra dos tipos de inactivos: faciles (al azar) y dificiles (los de propiedades
     fisicoquimicas mas parecidas a las de los activos: peso, logP, donadores, aceptores,
     enlaces rotables, al estilo DUD-E; criterio independiente de Tanimoto y de ErG).
  2. Compara 8 metodos: Tanimoto (max/media), ErG-fp, ErG+GED (max/media),
     fusion ErG+GED con Tanimoto, ErG+GED con filtro NS1 y solo el filtro NS1.
     'max' = similitud al activo mas parecido (vecino mas cercano); 'media' = promedio.
  3. Metricas: AUC con intervalo 95% (bootstrap) y enriquecimiento EF10/EF20.
  4. Opcional: --calibrar ajusta theta SOLO con el set de entrenamiento de cada prueba
     (sin fuga) y con la restriccion  w_polar_apolar >= w_aromatico_hidrofobo.
  5. Opcional: --cribado aplica los criterios NS1 y el ErG+GED a las ~33k moleculas.

USO (con el .venv activo, en la carpeta Universidad):
    python evaluar_erg_ns1.py --solo-referencias     # solo el chequeo de A9/A22
    python evaluar_erg_ns1.py                         # evaluacion (unos minutos)
    python evaluar_erg_ns1.py --calibrar              # con theta calibrado por prueba (lento)
    python evaluar_erg_ns1.py --cribado               # cribado de las 33k (~30 min)
Opciones de criterios: --tam-min 15  --ratio-extremos 0.7  --A-por-enlace 1.5 2.0  --exigir-tamano
Requiere haber corrido antes recalibrar_v2.py (crea AID2326_cache.csv) o tener BasePubChem.xlsx.
"""
import os
import sys
import time
import random
import argparse
import numpy as np
import pandas as pd

from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import AllChem, rdReducedGraphs, Descriptors, rdMolDescriptors
RDLogger.DisableLog("rdApp.*")

import ged_mejorado
import pipeline_final_completo as pf
import ns1_utils as U

ap = argparse.ArgumentParser(description="Evaluacion ErG+GED vs Tanimoto con criterios NS1")
ap.add_argument("--solo-referencias", action="store_true")
ap.add_argument("--cribado", action="store_true")
ap.add_argument("--calibrar", action="store_true")
ap.add_argument("--theta-archivo", default=None, help="json con {'theta': {...}} (ojo: puede tener fuga)")
ap.add_argument("--semillas", type=int, default=3)
ap.add_argument("--n-act", type=int, default=20)
ap.add_argument("--n-dec", type=int, default=40)
ap.add_argument("--iter", type=int, default=8, help="iteraciones de la calibracion")
ap.add_argument("--tam-min", type=float, default=15.0)
ap.add_argument("--ratio-extremos", type=float, default=0.70)
ap.add_argument("--A-por-enlace", type=float, nargs=2, default=[1.5, 2.0])
ap.add_argument("--exigir-tamano", action="store_true",
                help="incluir el criterio de tamano en el filtro NS1 (A9/A22 pueden no cumplirlo)")
args = ap.parse_args()

OUT = os.path.dirname(os.path.abspath(__file__))
P = {"tam_min": args.tam_min, "ratio_extremos": args.ratio_extremos,
     "A_por_enlace": tuple(args.A_por_enlace)}

if args.theta_archivo:
    import json
    THETA = json.load(open(os.path.join(OUT, args.theta_archivo)))["theta"]
    print(f"theta: {args.theta_archivo} (puede tener fuga de informacion entre entrenamiento y prueba)")
else:
    THETA = dict(pf.DEFAULT_THETA)
    print("theta: valores por defecto (sin calibrar, sin fuga)" if not args.calibrar
          else "theta: se calibra por prueba, solo con el set de entrenamiento")

METODOS = ["Tanimoto (max)", "Tanimoto (media)", "ErG-fp (max)", "ErG+GED (max)",
           "ErG+GED (media)", "Fusion ErG+GED + Tanimoto",
           "ErG+GED (max) + filtro NS1", "Solo filtro NS1"]


# ------------------------------------------------------------------ analisis de moleculas
def analizar(smi):
    mol3d, mol2d = pf.smiles_to_3d_mol(smi)
    G = pf.build_erg(mol3d)
    conf = mol3d.GetConformer()
    heavy = np.array([list(conf.GetAtomPosition(a.GetIdx())) for a in mol3d.GetAtoms()
                      if a.GetAtomicNum() > 1])
    nodos = [(G.nodes[i]["label"], np.asarray(G.nodes[i]["center"], float))
             for i in range(G.number_of_nodes())]
    aristas = [(i, j, int(d["dist"])) for i, j, d in G.edges(data=True)]
    return {"G": G,
            "fp": AllChem.GetMorganFingerprintAsBitVect(mol2d, radius=2, nBits=1024),
            "ergfp": np.asarray(rdReducedGraphs.GetErGFingerprint(mol2d), float),
            "crit": U.criterios_ns1(nodos, aristas, heavy, P)}


CACHE = {}


def item(cid, smi):
    if cid not in CACHE:
        try:
            CACHE[cid] = analizar(smi)
        except Exception:
            CACHE[cid] = None
    return CACHE[cid]


def largo_multiconf(smi, n=30):
    """Largo maximo (A) entre atomos pesados sobre n conformeros: da una idea de la
    extension real, porque un solo conformero puede estar plegado."""
    molH = Chem.AddHs(Chem.MolFromSmiles(smi))
    prm = AllChem.ETKDGv3(); prm.randomSeed = 0xF00D
    cids = list(AllChem.EmbedMultipleConfs(molH, numConfs=n, params=prm))
    idx = [a.GetIdx() for a in molH.GetAtoms() if a.GetAtomicNum() > 1]
    return max((U.largo_max(molH.GetConformer(c).GetPositions()[idx]) for c in cids), default=0.0)


def erg_sim(a, b, theta):
    if a["G"].number_of_nodes() == 0 or b["G"].number_of_nodes() == 0:
        return 0.0
    return ged_mejorado.ged_similarity_v2(a["G"], b["G"], theta)


def ergfp_sim(a, b):
    num = float(np.dot(a, b)); den = float(np.dot(a, a) + np.dot(b, b) - num)
    return num / den if den > 0 else 0.0


def propiedades(smi):
    m = Chem.MolFromSmiles(smi)
    return np.array([Descriptors.MolWt(m), Descriptors.MolLogP(m), rdMolDescriptors.CalcNumHBD(m),
                     rdMolDescriptors.CalcNumHBA(m), rdMolDescriptors.CalcNumRotatableBonds(m)], float)


def filtro_ok(c):
    return bool(c["estructura"] and (c["tam_ok"] or not args.exigir_tamano))


# ------------------------------------------------------------------ paso 0: referencias
def chequeo_referencias():
    print("\n" + "=" * 78)
    print("PASO 0 - ¿A9 y A22 cumplen los criterios del sitio de union de NS1?")
    print("=" * 78)
    print(f"Parametros: tamano > {P['tam_min']} A | aromaticos en extremos (ratio >= {P['ratio_extremos']}) | "
          f"{P['A_por_enlace'][0]}-{P['A_por_enlace'][1]} A por enlace")
    todo_ok = True
    for nombre, smi in pf.INFLUENZA_ACTIVOS.items():
        it = item(nombre, smi)
        c = it["crit"]
        largo_mc = largo_multiconf(smi)
        print(f"\n{nombre}:")
        print(f"  composicion (>=2 aromaticos, >=1 hidrofobico, >=1 donador, >=1 aceptor): {c['composicion']} "
              f"[arom={c['n_arom']} hidrof={c['n_hidrof']} HBD={c['n_hbd']} HBA={c['n_hba']}]")
        print(f"  grafo conectado (grupos a 1-6 enlaces): {c['conectado']}")
        print(f"  aromaticos en los extremos: {c['extremos']} (ratio {c['ratio_extremos']:.2f})")
        print(f"  largo de la molecula: {c['largo_A']:.1f} A (1 conformero) | {largo_mc:.1f} A (max de 30 conformeros) "
              f"-> {'cumple' if largo_mc > P['tam_min'] else 'NO supera'} {P['tam_min']} A")
        print(f"  A por enlace entre grupos: media {c['A_por_enlace_medio']:.2f} | "
              f"{100 * c['ajuste_dist']:.0f}% de las aristas dentro de {P['A_por_enlace']}")
        todo_ok &= c["estructura"] and largo_mc > P["tam_min"]
    if not todo_ok:
        print("\n AVISO: A9/A22 no cumplen TODOS los criterios (ver arriba). Si lo que buscamos son")
        print(" moleculas parecidas a ellos, conviene revisar que mide exactamente cada criterio")
        print(" (por ejemplo si el 15 son Angstrom, atomos pesados u otra medida). El tamano solo")
        print(" entra al filtro con --exigir-tamano.")


# ------------------------------------------------------------------ datos AID 2326
def cargar_aid():
    cache = os.path.join(OUT, "AID2326_cache.csv")
    if not os.path.exists(cache):
        print("Leyendo BasePubChem.xlsx (solo la primera vez, tarda unos minutos)...")
        df = pd.read_excel(os.path.join(OUT, "BasePubChem.xlsx"),
                           usecols=["PUBCHEM_CID", "PUBCHEM_EXT_DATASOURCE_SMILES", "PUBCHEM_ACTIVITY_OUTCOME"])
        df = df.dropna()
        df[df["PUBCHEM_ACTIVITY_OUTCOME"].isin(["Active", "Inactive"])].to_csv(cache, index=False)
    df = pd.read_csv(cache)
    df["cid"] = df["PUBCHEM_CID"].astype(float).astype(int).astype(str)
    act = df[df["PUBCHEM_ACTIVITY_OUTCOME"] == "Active"]
    ina = df[df["PUBCHEM_ACTIVITY_OUTCOME"] == "Inactive"].sample(n=5000, random_state=42)

    def limpiar(sub):
        d = {}
        for cid, smi in zip(sub["cid"], sub["PUBCHEM_EXT_DATASOURCE_SMILES"]):
            s = pf._limpiar_smiles_sal(smi)
            if s:
                d[cid] = s
        return d
    return limpiar(act), limpiar(ina)


# ------------------------------------------------------------------ calibracion (sin fuga)
def calibrar(train, decoys, seed):
    from scipy.optimize import differential_evolution
    y = [1] * len(train) + [0] * len(decoys)
    pool = train + decoys
    popsize = 6
    total = popsize * 6 * (args.iter + 1)
    st = {"n": 0, "t0": time.time(), "mejor": 0.0}

    def a_theta(v):
        return {"w_aromatico_hidrofobo": v[0], "w_hba_hbd": v[1], "w_polar_apolar": v[0] + v[2],
                "node_indel": v[3], "edge_beta": v[4], "edge_indel": v[5]}

    def objetivo(v):
        th = a_theta(v)
        sc = []
        for i, x in enumerate(pool):
            sims = [erg_sim(x, a, th) for j, a in enumerate(train) if not (i < len(train) and i == j)]
            sc.append(max(sims))
        a_ = U.auc(y, sc)
        st["n"] += 1; st["mejor"] = max(st["mejor"], a_)
        if st["n"] == 1 or st["n"] % max(1, total // 10) == 0:
            dt = time.time() - st["t0"]
            print(f"      calibracion {st['n']}/{total} | mejor AUC train {st['mejor']:.3f} | "
                  f"{dt / 60:.1f} min (~{(total - st['n']) * dt / st['n'] / 60:.0f} min restantes, maximo)")
        return -a_

    bounds = [(0.05, 2.5), (0.05, 2.5), (0.0, 2.0), (0.2, 3.0), (0.02, 1.5), (0.1, 2.0)]
    res = differential_evolution(objetivo, bounds, seed=seed, maxiter=args.iter, popsize=popsize,
                                 tol=1e-3, polish=False)
    return a_theta(res.x)


# ------------------------------------------------------------------ evaluacion
def puntajes(test, train, theta):
    tan = [[DataStructs.TanimotoSimilarity(t["fp"], a["fp"]) for a in train] for t in test]
    egp = [[ergfp_sim(t["ergfp"], a["ergfp"]) for a in train] for t in test]
    erg = [[erg_sim(t, a, theta) for a in train] for t in test]
    mx = lambda M: np.array([max(r) for r in M])
    md = lambda M: np.array([np.mean(r) for r in M])
    ok = np.array([1.0 if filtro_ok(t["crit"]) else 0.0 for t in test])
    erg_max, tan_max = mx(erg), mx(tan)
    return {"Tanimoto (max)": tan_max, "Tanimoto (media)": md(tan),
            "ErG-fp (max)": mx(egp), "ErG+GED (max)": erg_max, "ErG+GED (media)": md(erg),
            "Fusion ErG+GED + Tanimoto": (U.rango_norm(erg_max) + U.rango_norm(tan_max)) / 2,
            "ErG+GED (max) + filtro NS1": erg_max + ok, "Solo filtro NS1": ok}, ok


def evaluar(grupo, ids_train, ids_test, seed):
    t0 = time.time()
    train = [x for x in (item(c, activos[c]) for c in ids_train) if x]
    test_act = [x for x in (item(c, activos[c]) for c in ids_test) if x]
    if len(train) < 3 or len(test_act) < 2:
        print(f"\n>>> {grupo}: pocos datos ({len(train)} train, {len(test_act)} test), se omite")
        return
    print(f"\n>>> {grupo}: {len(train)} activos de entrenamiento, {len(test_act)} de prueba")

    rng = random.Random(seed)
    # inactivos dificiles: propiedades fisicoquimicas (z-score) mas cercanas a algun activo de entrenamiento
    z_tr = np.array([(propiedades(activos[c]) - PROP_MU) / PROP_SD for c in ids_train])
    dist_pool = {c: float(np.min(np.linalg.norm(z_tr - Z_INACT[c], axis=1))) for c in inact_ids}
    dificiles = sorted(inact_ids, key=lambda c: dist_pool[c])[:args.n_dec]
    resto = [c for c in inact_ids if c not in set(dificiles)]
    rng.shuffle(resto)
    faciles, calib_ids = resto[:args.n_dec], resto[args.n_dec:args.n_dec + 30]
    fps_tr = [a["fp"] for a in train]
    tan_medio = lambda ids: np.mean([max(DataStructs.BulkTanimotoSimilarity(fp_inact[c], fps_tr)) for c in ids])
    print(f"    inactivos dificiles (propiedades parecidas): distancia media {np.mean([dist_pool[c] for c in dificiles]):.2f} "
          f"vs faciles {np.mean([dist_pool[c] for c in faciles]):.2f} | "
          f"Tanimoto max medio al entrenamiento: dificiles {tan_medio(dificiles):.2f}, faciles {tan_medio(faciles):.2f}")

    theta = THETA
    if args.calibrar:
        dec_cal = [x for x in (item(c, inactivos[c]) for c in calib_ids) if x]
        print("    calibrando theta con el set de entrenamiento...")
        theta = calibrar(train, dec_cal, seed)
        print("    theta:", {k: round(v, 3) for k, v in theta.items()})

    for tipo, ids in (("faciles", faciles), ("dificiles", dificiles)):
        decs = [x for x in (item(c, inactivos[c]) for c in ids) if x]
        test = test_act + decs
        y = np.r_[np.ones(len(test_act)), np.zeros(len(decs))]
        sc, ok = puntajes(test, train, theta)
        pa = 100 * ok[:len(test_act)].mean(); pd_ = 100 * ok[len(test_act):].mean()
        for m in METODOS:
            lo, hi = U.auc_ci(y, sc[m], seed=seed)
            FILAS.append({"grupo": grupo, "decoys": tipo, "metodo": m, "AUC": U.auc(y, sc[m]),
                          "IC95_lo": lo, "IC95_hi": hi, "EF10": U.ef(y, sc[m], 0.10),
                          "EF20": U.ef(y, sc[m], 0.20), "n_train": len(train),
                          "n_act_test": len(test_act), "n_decoys": len(decs),
                          "pct_filtro_act": pa, "pct_filtro_dec": pd_})
        print(f"    [{tipo}] AUC  Tanimoto(max)={U.auc(y, sc['Tanimoto (max)']):.3f} | "
              f"ErG+GED(max)={U.auc(y, sc['ErG+GED (max)']):.3f} | "
              f"Fusion={U.auc(y, sc['Fusion ErG+GED + Tanimoto']):.3f} | "
              f"pasan filtro NS1: {pa:.0f}% de activos vs {pd_:.0f}% de inactivos")
    print(f"    ({time.time() - t0:.0f} s)")


def evaluacion():
    global activos, inactivos, fp_inact, inact_ids, FILAS, PROP_MU, PROP_SD, Z_INACT
    FILAS = []
    activos, inactivos = cargar_aid()
    print(f"\nAID 2326: {len(activos)} activos, muestra de {len(inactivos)} inactivos")
    clusters = pf.clusterizar_activos(activos, cutoff_distancia=0.6)
    inact_ids = list(inactivos)
    fp_inact = {c: AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), radius=2, nBits=1024)
                for c, s in inactivos.items()}
    props_inact = {c: propiedades(s) for c, s in inactivos.items()}
    M = np.array([props_inact[c] for c in inact_ids])
    PROP_MU, PROP_SD = M.mean(axis=0), M.std(axis=0) + 1e-9
    Z_INACT = {c: (props_inact[c] - PROP_MU) / PROP_SD for c in inact_ids}

    n = args.n_act
    for s in range(args.semillas):
        ids = random.Random(1000 + s).sample(list(activos), 2 * n)
        evaluar("mezclado", ids[:n], ids[n:], seed=s)
    for k in range(3):
        cl = list(clusters[k]); random.Random(k).shuffle(cl)
        corte = max(3, int(0.6 * len(cl)))
        evaluar(f"intra-cluster {k + 1}", cl[:corte], cl[corte:], seed=k)
        otros = [c for j in range(3) if j != k for c in clusters[j]]
        evaluar(f"entre-clusters: train {k + 1} -> test otros", list(clusters[k]), otros, seed=10 + k)

    df = pd.DataFrame(FILAS)
    df.to_csv(os.path.join(OUT, "evaluacion_erg_ns1_resultados.csv"), index=False)
    for tipo in ("faciles", "dificiles"):
        sub = df[df["decoys"] == tipo]
        for met, titulo in (("AUC", "AUC"), ("EF20", "Enriquecimiento EF20 (top 20% de la lista)")):
            print(f"\n{'=' * 100}\n{titulo} en moleculas de PRUEBA - inactivos {tipo.upper()}\n{'=' * 100}")
            print(sub.pivot_table(index="grupo", columns="metodo", values=met, aggfunc="mean")[METODOS]
                  .round(2).to_string())
    sub = df[(df["decoys"] == "dificiles") & (df["metodo"] == "Solo filtro NS1")]
    print("\n% de moleculas que pasan el filtro NS1 (activos de prueba vs inactivos dificiles):")
    print(sub.groupby("grupo")[["pct_filtro_act", "pct_filtro_dec"]].mean().round(0).to_string())
    print("\nGuardado: evaluacion_erg_ns1_resultados.csv (incluye intervalos de confianza 95%)")


# ------------------------------------------------------------------ cribado de las 33k
def cribado():
    print("\n" + "=" * 78 + "\nCRIBADO de Moleculas.xlsx con criterios NS1 + ErG+GED\n" + "=" * 78)
    refs = {k: item(k, s) for k, s in pf.INFLUENZA_ACTIVOS.items()}
    m = pd.read_excel(os.path.join(OUT, "Moleculas.xlsx"), usecols=["ChEMBL ID", "Smiles"]).dropna()
    total, filas, t0 = len(m), [], time.time()
    for n, (cid, smi) in enumerate(zip(m["ChEMBL ID"], m["Smiles"]), 1):
        try:
            it = analizar(smi)
        except Exception:
            continue
        c = it["crit"]
        fila = {"ID": cid, "Smiles": smi, **{k: c[k] for k in (
            "composicion", "conectado", "extremos", "estructura", "largo_A", "tam_ok",
            "A_por_enlace_medio", "ajuste_dist", "ratio_extremos")}}
        for k, r in refs.items():
            fila[f"ErG+GED__{k}"] = erg_sim(it, r, THETA)
            fila[f"Tanimoto__{k}"] = DataStructs.TanimotoSimilarity(it["fp"], r["fp"])
        fila["ErG_max"] = max(fila[f"ErG+GED__{k}"] for k in refs)
        fila["Tanimoto_max"] = max(fila[f"Tanimoto__{k}"] for k in refs)
        fila["cumple_NS1"] = filtro_ok(c)
        filas.append(fila)
        if n % 1500 == 0 or n == total:
            dt = time.time() - t0
            print(f"  {n}/{total} ({100 * n / total:.0f}%) | {dt / 60:.1f} min | "
                  f"~{(total - n) * dt / n / 60:.0f} min restantes")
    df = pd.DataFrame(filas)
    df.to_csv(os.path.join(OUT, "cribado_ns1.csv"), index=False)
    print("\nEmbudo de criterios (moleculas que los cumplen):")
    for nombre, mask in (("total analizadas", df["ID"].notna()),
                         ("composicion (2 arom + hidrof + donador + aceptor)", df["composicion"]),
                         ("+ grafo conectado", df["composicion"] & df["conectado"]),
                         ("+ aromaticos en extremos (= estructura NS1)", df["estructura"]),
                         (f"+ largo > {P['tam_min']} A", df["estructura"] & df["tam_ok"])):
        print(f"  {nombre:<55}{int(mask.sum()):>7}")
    top = df[df["cumple_NS1"]].sort_values("ErG_max", ascending=False).head(20)
    top.to_csv(os.path.join(OUT, "top20_ns1.csv"), index=False)
    print("\nTop-20 por ErG+GED entre las que cumplen el filtro NS1:")
    print(top[["ID", "ErG_max", "Tanimoto_max", "largo_A", "ratio_extremos"]].round(3).to_string(index=False))
    top_t = df.sort_values("Tanimoto_max", ascending=False).head(20)
    print(f"\nEn comun con el top-20 de Tanimoto (sin filtro): {len(set(top['ID']) & set(top_t['ID']))}/20")
    print("Guardado: cribado_ns1.csv, top20_ns1.csv")


if __name__ == "__main__":
    chequeo_referencias()
    if args.solo_referencias:
        sys.exit(0)
    if args.cribado:
        cribado()
    else:
        evaluacion()

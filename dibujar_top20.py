# ============================================================================
# dibujar_top20.py
#
# Explicacion detallada de este script (que hace cada paso, como se corre y
# que hay que cambiar para extrapolar el metodo a otra proteina):
#     docs/02_explicacion_del_codigo.md
# Explicacion de los conceptos (nodo, grupo farmacoforico, colapso de atomos,
# papel de la geometria 3D, GED, algoritmo hungaro):
#     docs/01_explicacion_del_metodo.md
# ============================================================================
"""
Dibuja las estructuras del top-20 (con A9 y A22 de referencia) para revisarlas con la
quimica: puntaje ErG+GED y Tanimoto contra A9 y A22 en cada molecula.
Tambien cuenta cuantas moleculas empatan con el puntaje del puesto 20.

Uso (con el .venv activo):
    python dibujar_top20.py            # usa los archivos *_v2  (ranking con GED corregido)
    python dibujar_top20.py _v2b       # usa los archivos *_v2b (GED corregido + theta recalibrado)
    python dibujar_top20.py ""         # usa los archivos originales (GED anterior)
Salida: top20_estructuras<sufijo>.png
"""
import os
import sys
import pandas as pd
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import AllChem, Draw
RDLogger.DisableLog("rdApp.*")

from pipeline_final_completo import INFLUENZA_ACTIVOS

SUF = sys.argv[1] if len(sys.argv) > 1 else "_v2"
OUT = os.path.dirname(os.path.abspath(__file__))

t9 = pd.read_csv(os.path.join(OUT, f"top20_A9{SUF}.csv"))
t22 = pd.read_csv(os.path.join(OUT, f"top20_A22{SUF}.csv"))
rank = pd.read_csv(os.path.join(OUT, f"ranking_33k_completo{SUF}.csv")).set_index("ID")

smiles = pd.concat([t9[["ID", "Smiles"]], t22[["ID", "Smiles"]]]).drop_duplicates("ID").set_index("ID")["Smiles"]
ids = list(t9["ID"]) + [i for i in t22["ID"] if i not in set(t9["ID"])]


def fp(smi):
    return AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), radius=2, nBits=1024)


ref_fp = {k: fp(v) for k, v in INFLUENZA_ACTIVOS.items()}

mols, legends = [], []
for k, smi in INFLUENZA_ACTIVOS.items():
    mols.append(Chem.MolFromSmiles(smi))
    legends.append(f"REFERENCIA {k}")

filas = []
for pos, cid in enumerate(ids, 1):
    smi = smiles[cid]
    f = fp(smi)
    e9, e22 = rank.loc[cid, "ErG+GED__A9"], rank.loc[cid, "ErG+GED__A22"]
    t_9 = DataStructs.TanimotoSimilarity(f, ref_fp["A9"])
    t_22 = DataStructs.TanimotoSimilarity(f, ref_fp["A22"])
    mols.append(Chem.MolFromSmiles(smi))
    legends.append(f"{pos}. {cid}\nErG A9 {e9:.3f} | A22 {e22:.3f}\nTanimoto A9 {t_9:.2f} | A22 {t_22:.2f}")
    filas.append({"puesto": pos, "ID": cid, "ErG_A9": e9, "ErG_A22": e22,
                  "Tanimoto_A9": t_9, "Tanimoto_A22": t_22, "Smiles": smi})

img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(380, 300), legends=legends)
ruta = os.path.join(OUT, f"top20_estructuras{SUF}.png")
img.save(ruta)
pd.DataFrame(filas).round(3).to_csv(os.path.join(OUT, f"top20_tabla{SUF}.csv"), index=False)
print(f"Guardado: {ruta} y top20_tabla{SUF}.csv")

# Empates en el borde del top-20
for ref in INFLUENZA_ACTIVOS:
    col = f"ErG+GED__{ref}"
    corte = t9[col].min() if ref == "A9" else t22[col].min()
    n = int((rank[col] >= corte - 1e-9).sum())
    print(f"{ref}: puntaje del puesto 20 = {corte:.4f} -> {n} moleculas tienen ese puntaje o mas"
          + ("  (hay empates en el borde: el top-20 es parcialmente arbitrario)" if n > 20 else ""))
print("Puntajes ErG distintos entre los 33k:", rank["ErG+GED__A9"].round(6).nunique())

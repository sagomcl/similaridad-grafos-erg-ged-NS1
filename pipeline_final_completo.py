"""
=============================================================================
pipeline_final_completo.py  --  NUCLEO del metodo ErG + GED + algoritmo hungaro
=============================================================================

Este archivo contiene todo el metodo, en el orden en que se usa:

  SECCION 1  smiles_to_3d_mol / build_erg   -> de SMILES a grafo ErG
  SECCION 2  ged_similarity                  -> similitud entre dos grafos ErG
             (version 1, conservada por trazabilidad; en produccion se usa la
              de ged_mejorado.py -- ver nota al final de la Seccion 2)
  SECCION 3  INFLUENZA_ACTIVOS / DENGUE_ACTIVOS -> moleculas de referencia
  SECCION 4  Tanimoto, Hamming, ErG-fp       -> metodos de comparacion
  SECCION 5  decoys reales + calibracion con curva ROC
  SECCION 6  carga de datos reales (PubChem AID 2326) y clustering de activos
  SECCION 7  ranking de bases grandes, por bloques y con checkpoint
  SECCION 8  ejemplo de ejecucion completa (ver aviso mas abajo)

PUNTOS DONDE SE CAMBIA ALGO PARA EXTRAPOLAR A OTRA DIANA
--------------------------------------------------------
  (a) LABEL_MAP / PRIORITY / TYPES  (Seccion 1): que familias farmacoforicas
      se consideran y con que prioridad se resuelven los solapamientos.
  (b) El rango 1 <= d <= 6 en build_erg: hasta que distancia (en enlaces) se
      conectan dos grupos farmacoforicos.
  (c) DEFAULT_THETA y los `bounds` de la calibracion (Secciones 2 y 5): los
      seis pesos del costo de edicion. No se editan a mano para una diana
      nueva: se recalibran y se guardan en un JSON.
  (d) INFLUENZA_ACTIVOS / DENGUE_ACTIVOS (Seccion 3): los inhibidores de
      referencia. Para otra proteina, se reemplaza este diccionario.
  (e) cutoff_distancia en clusterizar_activos (Seccion 6b): que tan parecidas
      deben ser dos moleculas para considerarse del mismo quimiotipo.
  (f) ns1_utils.PARAMS: los criterios geometricos del bolsillo (archivo aparte).

AVISO SOBRE LA SECCION 8 (bloque `if __name__ == "__main__"`)
-------------------------------------------------------------
Ese bloque es la demostracion original del pipeline y apunta a rutas que ya no
existen tal cual (merged_500_con_props.csv, data/Moleculas.xlsx,
AID_2326_datatable.csv). Se conserva sin cambios porque documenta el orden
logico de los pasos, pero NO es la forma de reproducir los resultados de la
tesis. Para eso estan los scripts:
    recalibrar_v2.py     -> calibrar theta (produce theta_final_v2.json)
    correr_ranking.py    -> aplicar el metodo a las 33 000 candidatas
    evaluar_erg_ns1.py   -> evaluacion justa frente a Tanimoto
    comparar_A9_A22.py   -> analisis detallado de un par de moleculas
    dibujar_top20.py     -> figura con las estructuras del top-20
"""
import os
import math
import json
from collections import deque

import numpy as np
import networkx as nx
import pandas as pd
import matplotlib
matplotlib.use("Agg")          # backend sin ventana: necesario para guardar PNG en servidor
import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment, differential_evolution
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score, roc_curve

import rdkit.Chem as Chem
from rdkit.Chem import AllChem, ChemicalFeatures, rdReducedGraphs, DataStructs, Descriptors
from rdkit import RDConfig


# =========================================================================
# SECCION 1 — De SMILES a grafo ErG (Extended Reduced Graph)
# =========================================================================
# BaseFeatures.fdef es el catalogo de definiciones farmacoforicas que trae
# RDKit: un archivo de reglas SMARTS que dice, por ejemplo, "un nitrogeno con
# hidrogeno es un donador", "un anillo aromatico es Aromatic", etc. Al
# construir la "feature factory" una sola vez al importar el modulo, se evita
# releer y recompilar esas reglas en cada una de las 33 000 moleculas.
fdef_path = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
feat_factory = ChemicalFeatures.BuildFeatureFactory(fdef_path)

# Traduccion de las familias de RDKit a las cuatro etiquetas del ErG.
# RDKit define mas familias (PosIonizable, NegIonizable, ZnBinder...) que aqui
# simplemente se ignoran: cualquier familia que no este en LABEL_MAP no genera
# nodo.
# >>> PARA EXTRAPOLAR (a): si la diana nueva depende de cargas (por ejemplo un
#     bolsillo con un residuo de aspartato), agregar aqui
#     "PosIonizable": "Cation", "NegIonizable": "Anion", y ampliar PRIORITY,
#     TYPES y la matriz de costos.
LABEL_MAP = {"Donor": "HBD", "Acceptor": "HBA", "Aromatic": "Aromatico", "Hydrophobe": "Hidrofobo"}

# PRIORITY resuelve el problema central de la construccion del ErG: un mismo
# atomo puede pertenecer a varias familias a la vez. El oxigeno de un fenol es
# donador Y aceptor; los carbonos de un anillo de benceno son aromaticos Y
# tambien hidrofobos. Si se creara un nodo por cada familia, el grafo CRECERIA
# en vez de reducirse, que es lo contrario del proposito de un ErG.
# La regla adoptada: cada atomo recibe UN solo tipo, el primero de esta lista
# que le aplique. El orden va de la propiedad mas especifica e informativa
# (donar un puente de hidrogeno) a la mas generica (ser apolar).
PRIORITY = ["HBD", "HBA", "Aromatico", "Hidrofobo"]

# Orden canonico de los tipos: define filas y columnas de la matriz de costos.
TYPES = ["Aromatico", "Hidrofobo", "HBA", "HBD"]
IDX = {t: i for i, t in enumerate(TYPES)}


def smiles_to_3d_mol(smiles, seed=0xF00D):
    """SMILES -> (molecula con conformero 3D y con hidrogenos, molecula 2D).

    Pasos:
      1. MolFromSmiles: interpreta la cadena SMILES y construye el grafo
         molecular (atomos y enlaces). Si el SMILES esta mal escrito devuelve
         None y aqui se lanza ValueError.
      2. AddHs: agrega los hidrogenos explicitos. Son indispensables para que
         la geometria 3D tenga sentido (los angulos de enlace dependen de
         ellos), aunque despues no se usen como nodos del ErG.
      3. ETKDGv3 + EmbedMolecule: genera coordenadas 3D. ETKDG ("Experimental
         Torsion Knowledge Distance Geometry") usa geometria de distancias con
         preferencias de torsion tomadas de estructuras cristalograficas
         reales, por lo que produce conformeros plausibles y no solo
         geometricamente validos. randomSeed fija la semilla: con la misma
         semilla, la misma molecula da siempre el mismo conformero, lo cual es
         necesario para que los resultados sean reproducibles.
      4. UFFOptimizeMolecule: relaja la geometria con el campo de fuerzas UFF
         ("Universal Force Field") para quitar tensiones residuales.

    Devuelve DOS objetos porque cumplen papeles distintos:
      - mol (3D, con H): de aqui salen los centros de los nodos y las
        coordenadas para los criterios geometricos.
      - mol2d (sin H): de aqui salen los fingerprints (Morgan, ErG-fp) y las
        propiedades fisicoquimicas, que se calculan sobre la topologia.

    Es el paso mas lento del pipeline (del orden de decenas de milisegundos
    por molecula) y la unica fuente de fallos habituales: moleculas muy
    grandes o con estereoquimica imposible hacen que EmbedMolecule devuelva
    -1 y se lance RuntimeError. Los scripts de ranking capturan esa excepcion
    y registran la molecula en un CSV de errores (17 de 33 047 en esta tesis).
    """
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
    """Construye el grafo ErG de una molecula con conformero 3D.

    ESTA ES LA FUNCION CENTRAL DEL METODO. Hace cuatro cosas:

    PASO A -- Asignar un tipo farmacoforico a cada atomo, sin solapamientos.
      Se piden a RDKit todas las "features" de la molecula. Cada feature trae
      una familia (Donor, Acceptor, Aromatic, Hydrophobe, ...) y la lista de
      atomos que la componen. Como un atomo puede aparecer en varias features,
      se conserva la de mayor prioridad segun PRIORITY. Resultado: el
      diccionario atom_type = {indice de atomo: etiqueta}. Los atomos que no
      pertenecen a ninguna de las cuatro familias (por ejemplo un carbono de
      una cadena entre dos grupos) simplemente no quedan en el diccionario y
      no generan nodo: son "pegamento" estructural.

    PASO B -- Colapsar atomos contiguos del mismo tipo en un solo nodo.
      Se recorre la molecula con una busqueda en anchura (BFS) sobre el grafo
      de enlaces, agrupando atomos que (i) tienen el mismo tipo y (ii) estan
      unidos por enlaces. Un anillo de benceno, que son seis carbonos
      aromaticos contiguos, se vuelve UN nodo "Aromatico". Esto es la
      "reduccion" del Extended Reduced Graph: se pasa de ~20 atomos a ~6
      nodos.
      Como no se pierde la informacion atomica: cada nodo guarda
      `atom_ids`, la lista exacta de atomos que lo componen. Siempre se puede
      volver de un nodo a sus atomos originales (comparar_A9_A22.py lo usa
      para imprimir los simbolos quimicos de cada nodo, y evaluar_erg_ns1.py
      para recuperar las coordenadas de los atomos pesados). El ErG es una
      vista agregada de la molecula, no una destruccion de la molecula.

    PASO C -- Guardar la geometria 3D de cada nodo.
      `center` es el centroide de las coordenadas 3D de los atomos del nodo.
      Es la posicion fisica del grupo farmacoforico en el espacio. No se usa
      en el calculo del GED (ver PASO D) pero si en los criterios geometricos
      de ns1_utils.py: distancias reales en angstroms entre grupos, largo de
      la molecula, aromaticos en los extremos.

    PASO D -- Conectar los nodos con aristas etiquetadas por distancia.
      Para cada par de nodos se calcula `shortest_bond_distance`: la longitud,
      EN NUMERO DE ENLACES, del camino mas corto entre cualquier atomo del
      primer nodo y cualquier atomo del segundo. Se crea arista solo si esa
      distancia esta entre 1 y 6.
        - Por que en enlaces y no en angstroms: la distancia topologica no
          depende del conformero. Una misma molecula embebida con otra semilla
          cambia sus angstroms pero no sus enlaces; usar enlaces hace que el
          grafo (y por tanto el ranking) sea estable y reproducible.
        - Por que el tope de 6: dos grupos separados por mas de 6 enlaces
          tienen tantos grados de libertad entre ellos que su distancia real
          en el espacio deja de estar determinada; la arista no informaria
          nada. Ademas acota el numero de aristas y por tanto el costo del
          GED.
      >>> PARA EXTRAPOLAR (b): este rango es un parametro del metodo. Para
          dianas con bolsillos mas grandes puede ampliarse (por ejemplo a 8);
          hay que recalibrar theta despues de cambiarlo, porque cambia el
          numero de aristas y con ello la escala del costo.

    Devuelve un `networkx.Graph` donde cada nodo tiene los atributos
    label / atom_ids / center, y cada arista el atributo dist.
    """
    # ---------- PASO A: un tipo por atomo, resolviendo solapamientos
    feats = feat_factory.GetFeaturesForMol(mol)
    atom_type = {}
    for f in feats:
        fam = f.GetFamily()
        if fam not in LABEL_MAP:
            continue                       # familia que no nos interesa
        label = LABEL_MAP[fam]
        for a in f.GetAtomIds():
            cur = atom_type.get(a)
            # se queda la etiqueta de mayor prioridad (indice menor en PRIORITY)
            if cur is None or PRIORITY.index(label) < PRIORITY.index(cur):
                atom_type[a] = label

    # grafo de enlaces de la molecula completa: se usa para el BFS del paso B
    # y para las distancias topologicas del paso D
    bond_graph = nx.Graph()
    for b in mol.GetBonds():
        bond_graph.add_edge(b.GetBeginAtomIdx(), b.GetEndAtomIdx())

    # ---------- PASO B y C: colapsar atomos contiguos del mismo tipo
    conf = mol.GetConformer()              # coordenadas 3D generadas por ETKDG
    visited, nodes = set(), []
    for a, label in atom_type.items():
        if a in visited:
            continue
        comp, dq = {a}, deque([a])         # componente conexa del mismo tipo
        visited.add(a)
        while dq:
            cur = dq.popleft()
            for nb in bond_graph.neighbors(cur):
                if nb in atom_type and atom_type[nb] == label and nb not in visited:
                    visited.add(nb); comp.add(nb); dq.append(nb)
        center = np.array([list(conf.GetAtomPosition(i)) for i in comp]).mean(axis=0)
        nodes.append({"label": label, "atom_ids": sorted(comp), "center": center})

    # ---------- PASO D: distancia topologica entre dos grupos de atomos
    def shortest_bond_distance(atomsA, atomsB):
        """BFS multi-origen: menor numero de enlaces entre el conjunto A y B."""
        dq = deque([(i, 0) for i in atomsA])
        seen, targets = set(atomsA), set(atomsB)
        while dq:
            v, d = dq.popleft()
            if v in targets:
                return d
            for nb in bond_graph.neighbors(v):
                if nb not in seen:
                    seen.add(nb); dq.append((nb, d + 1))
        return math.inf                    # moleculas con fragmentos separados

    G = nx.Graph()
    for i, n in enumerate(nodes):
        G.add_node(i, **n)
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            d = shortest_bond_distance(nodes[i]["atom_ids"], nodes[j]["atom_ids"])
            if 1 <= d <= 6:                # <<< rango de conexion (extrapolable)
                G.add_edge(i, j, dist=int(d))
    return G


def erg_visualization_data(smiles):
    """Exporta el grafo ErG de CUALQUIER molecula en un dict listo para la web.

    Devuelve atomos con coordenadas 2D, enlaces, el tipo asignado a cada atomo,
    los nodos ErG con su posicion y las aristas con su distancia. Es lo que
    alimenta explorador_erg_3d.html y las figuras paso a paso de la
    presentacion (molecula -> atomos tipificados -> nodos colapsados -> grafo).

    Ojo: aqui los tipos se recalculan sobre mol2d (sin hidrogenos) porque las
    coordenadas que se exportan son 2D, de dibujo. Para el metodo se usa
    siempre build_erg(mol3d).
    """
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
# Los seis pesos que definen el costo de edicion. theta es el UNICO lugar
# donde entra criterio quimico ajustable en la medida de similitud:
#
#   w_aromatico_hidrofobo : costo de confundir un aromatico con un hidrofobo
#   w_hba_hbd             : costo de confundir un aceptor con un donador
#   w_polar_apolar        : costo de confundir un grupo polar con uno apolar
#   node_indel            : costo de insertar o borrar un nodo entero
#                           (una molecula tiene un grupo que la otra no tiene)
#   edge_beta             : costo por cada enlace de diferencia entre dos
#                           aristas emparejadas (estirar o encoger la
#                           separacion entre dos grupos)
#   edge_indel            : costo de insertar o borrar una arista
#
# Estos valores por defecto son una conjetura inicial razonable, NO los
# usados en los resultados: los finales estan en theta_final_v2.json y se
# obtuvieron con recalibrar_v2.py.
# >>> PARA EXTRAPOLAR (c): no editar estos numeros a mano para una diana
#     nueva. Recalibrarlos con activos de esa diana.
DEFAULT_THETA = {
    "w_aromatico_hidrofobo": 0.50,
    "w_hba_hbd": 0.75,
    "w_polar_apolar": 1.50,
    "node_indel": 1.00,
    "edge_beta": 0.25,
    "edge_indel": 0.50,
}


def cost_matrix_from_theta(theta):
    """Matriz 4x4 de costos de sustitucion de tipos de nodo.

    El orden de filas y columnas es TYPES = [Aromatico, Hidrofobo, HBA, HBD].
    Es simetrica y con diagonal cero: sustituir un tipo por si mismo es gratis.
    """
    w_ar_hy = theta["w_aromatico_hidrofobo"]
    w_hba_hbd = theta["w_hba_hbd"]
    w_pa = theta["w_polar_apolar"]
    return np.array([
        [0.0,      w_ar_hy,  w_pa,      w_pa],
        [w_ar_hy,  0.0,      w_pa,      w_pa],
        [w_pa,     w_pa,     0.0,       w_hba_hbd],
        [w_pa,     w_pa,     w_hba_hbd, 0.0],
    ])


def ged_similarity(G1, G2, theta=DEFAULT_THETA, verbose=False):
    """GED aproximado entre dos grafos ErG, resuelto con el algoritmo hungaro.

    VERSION 1 -- CONSERVADA POR TRAZABILIDAD. Tiene un defecto conocido: el
    hungaro elige el emparejamiento mirando SOLO los tipos de nodo, de modo
    que con nodos del mismo tipo hay empates y puede quedar un emparejamiento
    cruzado; el costo de aristas se calcula despues sobre ese emparejamiento
    ya fijado. Por eso A9 y A22, con grafos ErG identicos, daban 0,43 en vez
    de 1,00. La version corregida esta en ged_mejorado.py y es la que produce
    todos los resultados reportados. Esta se mantiene para poder reproducir
    los numeros antiguos y documentar la correccion.

    COMO FUNCIONA (es el esquema de Riesen & Bunke 2009):

    1. Se arma una matriz cuadrada Cn de tamano N x N con N = max(n1, n2).
       El bloque real (n1 x n2) tiene los costos de sustituir el tipo del
       nodo i de G1 por el tipo del nodo j de G2. Las filas y columnas de
       relleno representan nodos ficticios y valen `node_indel`: emparejar un
       nodo real con un ficticio equivale a borrarlo o insertarlo.

    2. `linear_sum_assignment` resuelve el problema de asignacion: elegir una
       pareja para cada fila, sin repetir columnas, minimizando la suma de
       costos. Es exactamente el ALGORITMO HUNGARO (Kuhn-Munkres), que lo
       resuelve de forma optima en tiempo polinomico O(N^3) -- frente al GED
       exacto, que es NP-dificil. Aqui esta el aporte del hungaro: convierte
       un problema intratable en uno que corre en milisegundos.

    3. Con el emparejamiento ya decidido se recorren las aristas de G1. Si la
       arista correspondiente existe en G2, se cobra edge_beta por cada enlace
       de diferencia en la distancia; si no existe, se cobra edge_indel. Al
       final se cobra edge_indel por cada arista de G2 que nadie uso.

    4. Normalizacion y conversion a similitud:
           GED_norm = GED / (n1 + n2 + m1 + m2)
           similitud = exp(-GED_norm)
       El denominador es el numero maximo de operaciones de edicion posibles,
       asi que divide el costo por "lo mucho que habia para editar" y quita el
       sesgo por tamano de molecula. exp(-x) lleva el costo a una similitud en
       (0, 1] de forma suave y monotona.
    """
    C_SUB = cost_matrix_from_theta(theta)
    node_indel, edge_beta, edge_indel = theta["node_indel"], theta["edge_beta"], theta["edge_indel"]

    # ---- 1. matriz de costos de nodos, cuadrada con nodos ficticios
    n1, n2 = G1.number_of_nodes(), G2.number_of_nodes()
    N = max(n1, n2)
    Cn = np.zeros((N, N))
    for i in range(n1):
        li = G1.nodes[i]["label"]
        for j in range(n2):
            Cn[i, j] = C_SUB[IDX[li], IDX[G2.nodes[j]["label"]]]
    for i in range(n1, N):
        Cn[i, :n2] = node_indel        # insertar un nodo que G1 no tiene
    for i in range(n1):
        Cn[i, n2:] = node_indel        # borrar un nodo que G2 no tiene

    # ---- 2. ALGORITMO HUNGARO: emparejamiento optimo de nodos
    r, c = linear_sum_assignment(Cn)
    cost_nodes = Cn[r, c].sum()
    mapping = {int(i): int(j) for i, j in zip(r, c) if i < n1 and j < n2}

    # ---- 3. costo de aristas sobre ese emparejamiento
    used, cost_edges = set(), 0.0
    for (u, v, d) in G1.edges(data=True):
        if u in mapping and v in mapping:
            u2, v2 = mapping[u], mapping[v]
            if G2.has_edge(u2, v2):
                # misma relacion entre los mismos grupos: se paga la
                # diferencia de distancia, en enlaces
                cost_edges += edge_beta * abs(int(d["dist"]) - int(G2.edges[u2, v2]["dist"]))
                used.add(tuple(sorted((u2, v2))))
            else:
                cost_edges += edge_indel      # la relacion no existe en G2
        else:
            cost_edges += edge_indel          # alguno de los nodos fue borrado
    for (x, y, _) in G2.edges(data=True):
        if tuple(sorted((x, y))) not in used:
            cost_edges += edge_indel          # aristas de G2 que hay que insertar

    # ---- 4. normalizacion y similitud
    ged = cost_nodes + cost_edges
    denom = n1 + n2 + G1.number_of_edges() + G2.number_of_edges()
    n_ged = ged / denom if denom > 0 else 0.0
    similarity = math.exp(-n_ged)

    if verbose:
        print(f"  emparejamiento G1->G2: {mapping}")
        print(f"  costo_nodos={cost_nodes:.3f}  costo_aristas={cost_edges:.3f}  "
              f"GED={ged:.3f}  GED_norm={n_ged:.3f}  similitud={similarity:.3f}")
    return similarity


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
#
# DATO IMPORTANTE PARA LA TESIS: A9 y A22 difieren en un solo atomo (OH vs
# NH2 en la posicion orto del anillo bencenico). Con la reduccion ErG ese
# cambio desaparece, porque tanto -OH como -NH2 son donadores de puente de
# hidrogeno: los dos grafos ErG resultan identicos (6 nodos, 14 aristas) y su
# similitud es 1,000. Es decir, en Fase A las dos referencias aportan la
# misma informacion al metodo.
#
# >>> PARA EXTRAPOLAR (d): aqui se cambian las referencias. Un diccionario
#     {nombre: SMILES} con los inhibidores conocidos de la diana nueva es
#     todo lo que necesita el resto del pipeline.
INFLUENZA_ACTIVOS = {
    "A9":  "O=C(NCCc1ccccc1)c1ccccc1O",
    "A22": "Nc1ccccc1C(=O)NCCc1ccccc1",
}

# --- 3b. Benchmark adicional de influenza CON IC50 real (opcional, mas rico) ---
# Fuente: Basu et al. 2009, J. Virol., DOI 10.1128/JVI.01805-08 — Tabla 1.
# SMILES verificados via el resolvedor NCI Cactus (nombre NSC -> estructura).
# Sirve para una comprobacion blanda: ¿moleculas con IC50 parecido salen con
# similitud parecida? Con n=4 no da para una metrica, solo para inspeccion.
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
# Son tres productos naturales quimicamente muy distintos entre si (un
# alcaloide de isoquinolina y dos triterpenos pentaciclicos). Por el hallazgo
# central de la Fase A, NO deben promediarse: cada uno se usa como referencia
# por separado y produce su propio ranking.
DENGUE_ACTIVOS = {
    "Berberina":        "O1c2c(OC1)cc5c(c2)c4cc3ccc(OC)c(OC)c3c[n+]4CC5",
    "Acido_betulinico": "CC(=C)C1CCC2(C1C3CCC4C5(CCC(C(C5CCC4(C3(CC2)C)C)(C)C)O)C)C(=O)O",
    "Acido_ursolico":   "CC1CCC2(CCC3(C(=CCC4C3(CCC5C4(CCC(C5(C)C)O)C)C)C2C1C)C)C(=O)O",
}


# =========================================================================
# SECCION 4 — Metodos de comparacion (para saber si ErG+GED sirve o no)
# =========================================================================
# Un metodo nuevo solo se puede defender contra una linea base. La linea base
# estandar en cribado por similitud es el fingerprint de Morgan (ECFP4) con
# coeficiente de Tanimoto. Se incluyen tambien Hamming y el fingerprint ErG
# que trae RDKit, para separar dos efectos distintos: el de usar una
# representacion reducida (ErG-fp) y el de usar el grafo explicito con GED.

def tanimoto_bits(fp1, fp2):
    """Coeficiente de Tanimoto entre dos fingerprints de bits.

    T = c / (a + b - c), con a y b los bits encendidos de cada molecula y c
    los comunes. Es la medida mas usada en quimioinformatica. Mide
    coincidencia de SUBESTRUCTURAS exactas: dos moleculas con el mismo
    farmacoforo pero anillos distintos pueden tener Tanimoto bajo.
    """
    return DataStructs.TanimotoSimilarity(fp1, fp2)


def hamming_similarity_bits(fp1, fp2):
    """1 - (fraccion de bits en que difieren).

    A diferencia de Tanimoto, premia tambien los ceros coincidentes. Con
    fingerprints dispersos (la mayoria de los 1024 bits en cero) esto infla la
    similitud de cualquier par de moleculas, por lo que suele ser una medida
    pobre; se incluye porque el anteproyecto la contemplaba como comparacion.
    """
    a1 = np.zeros((len(fp1),), dtype=int); DataStructs.ConvertToNumpyArray(fp1, a1)
    a2 = np.zeros((len(fp2),), dtype=int); DataStructs.ConvertToNumpyArray(fp2, a2)
    return 1 - np.sum(a1 != a2) / len(a1)


def erg_fp_cosine(fp1, fp2):
    """Similitud tipo Tanimoto continua entre dos fingerprints ErG de RDKit.

    El ErG-fp de RDKit es un vector de flotantes (no de bits): cuenta pares de
    tipos farmacoforicos a cada distancia topologica. Esta funcion aplica la
    forma generalizada de Tanimoto para vectores continuos,
    <a,b> / (<a,a> + <b,b> - <a,b>).
    Es la comparacion mas informativa de todas, porque usa la MISMA idea de
    reduccion farmacoforica que el metodo propuesto pero sin grafo explicito
    ni GED: la diferencia entre ErG-fp y ErG+GED aisla exactamente cuanto
    aporta el grafo.
    """
    a, b = np.array(fp1), np.array(fp2)
    num = np.dot(a, b)
    den = np.dot(a, a) + np.dot(b, b) - num
    return num / den if den > 0 else 0.0


def _prep(smi):
    """Prepara una molecula una sola vez: grafo ErG + los dos fingerprints + MW.

    Se usa para no repetir el embedding 3D (lo mas costoso) cada vez que una
    molecula se compara con otra referencia.
    """
    mol3d, mol2d = smiles_to_3d_mol(smi)
    return {
        "erg_graph": build_erg(mol3d),
        "erg_fp": rdReducedGraphs.GetErGFingerprint(mol2d),
        "morgan_fp": AllChem.GetMorganFingerprintAsBitVect(mol2d, radius=2, nBits=1024),
        "mw": Descriptors.MolWt(mol2d),
    }


def compute_all_similarities(smiles_dict, referencia_smiles_dict, theta=DEFAULT_THETA,
                              por_activo=False):
    """Para cada molecula de smiles_dict, su similitud a las referencias con
    los cuatro metodos.

    Si por_activo=True, agrega una columna por CADA referencia por separado
    (ErG+GED__A9, ErG+GED__A22, ...). Es el modo recomendado: promediar la
    similitud sobre activos quimicamente distintos diluye la senal real. Ese
    fue el hallazgo central de la Fase A y por eso todos los resultados de la
    tesis se reportan por activo.

    Las columnas sin sufijo (ErG+GED, Morgan+Tanimoto, ...) son el promedio
    sobre referencias; se conservan por compatibilidad pero no deben usarse
    para concluir nada cuando las referencias no comparten quimiotipo.
    """
    ref_data = {name: _prep(smi) for name, smi in referencia_smiles_dict.items()}
    rows = []
    for name, smi in smiles_dict.items():
        d = _prep(smi)
        row = {"Molecula": name}
        s_erg_ged, s_ergfp, s_morgan_tan, s_morgan_ham = [], [], [], []
        for ref_name, r in ref_data.items():
            if ref_name == name:
                continue                 # no compararse consigo misma
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
#
# ADVERTENCIA METODOLOGICA: calibrar_theta_con_roc calibra y evalua sobre el
# MISMO conjunto, asi que su AUC mide memorizacion, no generalizacion. Para
# resultados publicables usar calibrar_y_evaluar_train_test (Seccion 6c) o
# recalibrar_v2.py.
# =========================================================================

def construir_decoys_reales(df_candidatos, activos_smiles, smiles_col="Smiles",
                             mw_tolerancia=100, tanimoto_max=0.08, n_decoys=20, seed=42):
    """Selecciona decoys (supuestos inactivos) con el criterio tipo DUD-E.

    Se queda con moleculas de la base propia que tengan peso molecular en el
    rango de los activos (+-mw_tolerancia) pero Tanimoto maximo contra todos
    los activos por debajo de tanimoto_max. Las ordena de menor a mayor
    parecido y toma las n_decoys mas distintas.

    Por que importa: si los decoys son demasiado faciles (moleculas enormes o
    diminutas), cualquier metodo da AUC alto y la validacion no vale nada. Si
    son demasiado parecidos, se estan tratando como inactivos moleculas que
    probablemente serian activas. El filtro por propiedades + disimilitud
    estructural es el compromiso estandar.

    LIMITACION: son decoys PRESUNTOS (nunca se midieron en el laboratorio).
    En la Fase A de esta tesis se reemplazaron por inactivos CONFIRMADOS de
    PubChem AID 2326, que es mucho mas solido. Esta funcion se conserva para
    el caso de una diana nueva sin cribado disponible -- situacion exacta del
    dengue en la Fase B.
    """
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


def calibrar_theta_con_roc(activos_smiles, decoys_smiles, out_prefix,
                            bounds=None, maxiter=60, seed=42):
    """Calibra theta maximizando el AUC y guarda pesos, curva ROC y scores.

    COMO SE CALIBRA (este es el nucleo del objetivo especifico 3):
      - Se construyen los grafos ErG de activos y decoys una sola vez.
      - El "puntaje" de una molecula es su similitud media a los activos
        (excluyendose a si misma, para no inflar el resultado de los activos).
      - La funcion objetivo es -AUC: separar activos de decoys lo mejor
        posible.
      - `differential_evolution` (Storn & Price 1997) es un optimizador global
        por poblaciones. Se usa en vez de un metodo de gradiente porque el AUC
        es una funcion ESCALONADA de theta: cambia solo cuando dos moleculas
        intercambian su posicion en el orden, asi que su gradiente es cero
        casi en todas partes y los metodos de descenso no sirven.

    `bounds` fija el rango admisible de cada peso, en el mismo orden que las
    claves de DEFAULT_THETA. Los limites no son arbitrarios: mantienen los
    pesos en una escala comparable y evitan soluciones degeneradas (por
    ejemplo todos los costos en cero, que haria que todo se parezca a todo).
    >>> PARA EXTRAPOLAR (c): si en la diana nueva se tiene la intuicion de que,
        por ejemplo, confundir polar con apolar debe ser siempre mas grave que
        confundir aromatico con hidrofobo, se puede imponer como restriccion.
        evaluar_erg_ns1.py lo hace: parametriza w_polar_apolar = v[0] + v[2]
        con v[2] >= 0, lo que garantiza w_polar_apolar >= w_aromatico_hidrofobo.
    """
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
    """Lee un Excel/CSV con una columna de SMILES -> {id: smiles}.

    >>> PARA EXTRAPOLAR: para cribar otra biblioteca de compuestos basta
        apuntar `path` al archivo nuevo y ajustar smiles_col / id_col. No hay
        nada especifico de ChEMBL en el resto del pipeline.
    """
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
    (la molecula real), no con el contraion.

    Sin esta limpieza, el embedding 3D y el ErG se calcularian sobre un
    sistema de dos fragmentos desconectados: el grafo saldria desconectado y
    la similitud no tendria sentido quimico. Se usa el fragmento con el SMILES
    mas largo como aproximacion al mas grande.
    """
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
    Devuelve (activos, inactivos): dos dict {CID: smiles_limpio (sin sales)}.

    Este ensayo es el `ground truth` de toda la Fase A: sin el, no habria con
    que medir si el metodo discrimina o no. Solo se usan Active e Inactive;
    Inconclusive y Unspecified se descartan por no aportar etiqueta fiable.

    OJO con `skiprows=[1, 2, 3, 4]`: el CSV que PubChem entrega trae cuatro
    filas de metadatos debajo del encabezado. Si el archivo no las tiene (por
    ejemplo AID2326_cache.csv, generado por recalibrar_v2.py a partir del
    Excel), este salto borraria cuatro moleculas reales. pipeline_v3.py filtra
    por contenido en vez de por posicion, que es mas seguro.
    """
    df = pd.read_csv(csv_path, skiprows=[1, 2, 3, 4], low_memory=False)
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
    """Agrupa activos en quimiotipos con el algoritmo de Butina (1999).

    Como funciona Butina: calcula todas las distancias (1 - Tanimoto) entre
    pares, elige como centro de cluster la molecula con mas vecinos dentro del
    umbral, le asigna esos vecinos, los retira y repite. Es el metodo estandar
    en quimioinformatica porque no exige fijar el numero de clusters de
    antemano y respeta la idea de "vecindario quimico".

    cutoff_distancia=0.6 significa que dos moleculas quedan en el mismo
    cluster si su Tanimoto es mayor que 0,4.
    >>> PARA EXTRAPOLAR (e): bajar el cutoff (por ejemplo 0,4, o Tanimoto>0,6)
        produce quimiotipos mas estrictos y mas pequenos; subirlo produce
        grupos mas laxos. Es un parametro con efecto directo sobre el AUC que
        se va a reportar, asi que conviene dejarlo fijo y declararlo.

    Devuelve la lista de clusters (listas de ids) ordenada de mayor a menor.
    En esta tesis: 1 073 activos -> 737 clusters, 588 de ellos de una sola
    molecula. Ese numero es el resultado que explica por que el AUC agregado
    daba 0,49.
    """
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

    # matriz de distancias en formato triangular condensado, como la pide Butina
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
                                   maxiter=15, popsize=8, seed=42):
    """Calibra theta solo con el conjunto de entrenamiento y evalua en prueba.

    Detalle que importa: el puntaje de una molecula de PRUEBA se calcula como
    su similitud media a los activos de ENTRENAMIENTO (nunca a los de prueba).
    Asi se reproduce la situacion real de un cribado: se conocen unos pocos
    inhibidores y se quiere ordenar una biblioteca desconocida.

    Devuelve (theta, auc_train, auc_test, detalle). La diferencia entre
    auc_train y auc_test es la medida directa del sobreajuste: en esta tesis,
    con n=2 activos, auc_train era 1,00 y el auc_test sobre activos reales
    cayo a 0,49.
    """
    rng = np.random.RandomState(seed)

    def split(d):
        keys = list(d.keys())
        rng.shuffle(keys)
        k = max(1, int(len(keys) * frac_train))
        return keys[:k], keys[k:]

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
                mol3d, _ = smiles_to_3d_mol(smiles_dict[k])
                graphs[k] = build_erg(mol3d)
            except Exception:
                continue           # moleculas que no se pudieron embeber en 3D
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
        # similitud a los activos de ENTRENAMIENTO: nada de informacion de prueba
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


# =========================================================================
# SECCION 7 — Ranking sobre las 33k (o una muestra) contra A9/A22
# =========================================================================

def rankear_candidatos(candidatos_smiles, activos_smiles, theta, top_n=20):
    """Ranking de un diccionario de candidatos, con los cuatro metodos.

    Imprime cuantas moleculas comparten el top-N de ErG+GED y el de
    Morgan+Tanimoto: es la forma mas directa de ver si el metodo propuesto
    aporta algo distinto o solo reordena lo mismo que un fingerprint.
    """
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
    """Procesa un bloque de `batch_size` moleculas y AGREGA filas a out_csv.

    Patron de checkpoint: el .json guarda {"next_idx": ..., "total": ...}. Al
    arrancar, si el .json existe se retoma desde next_idx; si no existe, se
    empieza de cero y se borra el CSV anterior. Hay que llamarla repetidamente
    (correr_ranking.py lo hace en un while) hasta que next_idx == total.

    Las referencias (activos) se preparan UNA vez al inicio: su grafo ErG y su
    fingerprint se reutilizan para las 33 000 comparaciones.

    `except Exception: continue` descarta en silencio las moleculas cuyo
    embedding 3D falla. pipeline_v3.py cambia esto por un registro en un CSV
    de errores, que es lo que permitio contar exactamente 17 moleculas
    descartadas de 33 047 y reportarlo en la tesis.
    """
    ref_data = {}
    for name, smi in activos_smiles.items():
        mol3d, mol2d = smiles_to_3d_mol(smi)
        ref_data[name] = {
            "erg_graph": build_erg(mol3d),
            "morgan_fp": AllChem.GetMorganFingerprintAsBitVect(mol2d, radius=2, nBits=1024),
        }

    if os.path.exists(checkpoint_json):
        with open(checkpoint_json) as fh:
            start_idx = json.load(fh)["next_idx"]
    else:
        start_idx = 0
        if os.path.exists(out_csv):
            os.remove(out_csv)

    if str(path_excel).lower().endswith(".csv"):
        df_full = pd.read_csv(path_excel)
    else:
        df_full = pd.read_excel(path_excel, sheet_name=sheet_name)
    df_full = df_full.dropna(subset=[smiles_col]).reset_index(drop=True)
    total = len(df_full)
    end_idx = min(start_idx + batch_size, total)
    print(f"  Procesando filas {start_idx} a {end_idx} de {total}")

    rows = []
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
                # una columna por referencia: NUNCA promediadas
                row[f"ErG+GED__{ref_name}"] = ged_similarity(g, rd["erg_graph"], theta)
                row[f"Morgan+Tanimoto__{ref_name}"] = tanimoto_bits(morgan_fp, rd["morgan_fp"])
            rows.append(row)
        except Exception:
            continue

    pd.DataFrame(rows).to_csv(out_csv, mode="a", header=(start_idx == 0), index=False)
    with open(checkpoint_json, "w") as fh:
        json.dump({"next_idx": end_idx, "total": total}, fh)
    print(f"  Progreso: {end_idx}/{total} ({100*end_idx/total:.1f}%)")
    return end_idx, total


# =========================================================================
# SECCION 8 — Ejecucion completa (DEMOSTRACION HISTORICA, ver aviso arriba)
# =========================================================================
# Este bloque documenta el orden logico de los ocho pasos del proyecto, pero
# apunta a rutas que ya no corresponden a la carpeta actual. Para reproducir
# los resultados de la tesis, usar los scripts (recalibrar_v2.py,
# correr_ranking.py, evaluar_erg_ns1.py). Se conserva sin modificaciones para
# no alterar el historial del codigo.
if __name__ == "__main__":
    OUT = os.path.dirname(os.path.abspath(__file__))

    print("=" * 78)
    print("PASO 1 — Datos de visualizacion ErG (para tu presentacion)")
    print("=" * 78)
    viz = erg_visualization_data(INFLUENZA_ACTIVOS["A9"])
    with open(os.path.join(OUT, "viz_A9.json"), "w") as fh:
        json.dump(viz, fh, indent=2, default=str)
    print(f"  Grafo ErG de A9: {len(viz['nodes'])} nodos, {len(viz['edges'])} aristas "
          f"(guardado en viz_A9.json, reutiliza la logica de molecula_a_erg_aspirina.html)")

    print("\n" + "=" * 78)
    print("PASO 2 — Decoys reales sacados de tu dataset (500-muestra de Moleculas.xlsx)")
    print("=" * 78)
    df_500 = pd.read_csv(os.path.join(OUT, "merged_500_con_props.csv"))
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
    csv_pubchem = os.path.join(OUT, "AID_2326_datatable.csv")
    if os.path.exists(csv_pubchem):
        activos_reales, inactivos_reales = cargar_pubchem_ns1_assay(csv_pubchem)

        # 6a) clusterizar los activos reales: aisla chemotipos coherentes
        # (un cribado fenotipico mezcla mecanismos MUY distintos -- ver
        # Seccion 6b). Sin esto, calibrar contra los ~1073 mezclados da
        # AUC~0.5 (azar): no es que el metodo falle, es que promedias
        # sobre activos sin relacion estructural entre si.
        clusters = clusterizar_activos(activos_reales, cutoff_distancia=0.6)
        activos_coherentes = {}
        for cluster in clusters[:3]:  # los 3 clusters mas grandes/coherentes
            for cid in cluster:
                activos_coherentes[cid] = activos_reales[cid]
        print(f"  activos coherentes combinados (3 clusters mas grandes): {len(activos_coherentes)}")

        # 6b) muestra manejable de inactivos reales (miles son demasiados
        # para el embedding 3D uno por uno en tiempo razonable)
        import random
        random.seed(42)
        inact_sample_ids = random.sample(list(inactivos_reales.keys()), min(80, len(inactivos_reales)))
        inactivos_sample = {k: inactivos_reales[k] for k in inact_sample_ids}

        # 6c) calibracion train/test HONESTA (AUC de generalizacion, no de
        # memorizacion -- ver advertencia en Seccion 6c)
        theta_final, auc_train, auc_test, detalle = calibrar_y_evaluar_train_test(
            activos_coherentes, inactivos_sample, maxiter=8, popsize=5)
        print(f"  RESULTADO VALIDACION: AUC train={auc_train:.3f}  AUC test={auc_test:.3f}  {detalle}")
        with open(os.path.join(OUT, "theta_final.json"), "w") as fh:
            json.dump({"theta": theta_final, "auc_train": auc_train, "auc_test": auc_test}, fh, indent=2)
    else:
        print(f"  (Salta este paso: no se encontro {csv_pubchem}.")
        print("   Descargalo desde https://pubchem.ncbi.nlm.nih.gov/bioassay/2326 ->")
        print("   Download -> 'Deposited Data Table Only' -> CSV, y ponlo junto a este script.)")
        theta_final = theta_flu  # usa el theta de PASO 3 como respaldo

    print("\n" + "=" * 78)
    print("PASO 7 — Aplicacion real: ranking de TODA la base de candidatos (33k+)")
    print("=" * 78)
    excel_candidatos = os.path.join(OUT, "data", "Moleculas.xlsx")
    if os.path.exists(excel_candidatos):
        # Esto puede tardar 30-45 min para 33k moleculas reales. Si tu
        # entorno corta ejecuciones largas, llama esta misma linea varias
        # veces seguidas -- retoma solo donde quedo gracias al checkpoint.
        idx, total = rankear_excel_en_bloques(
            excel_candidatos, INFLUENZA_ACTIVOS, theta_final,
            out_csv=os.path.join(OUT, "ranking_33k_completo.csv"),
            checkpoint_json=os.path.join(OUT, "ranking_33k_checkpoint.json"),
            batch_size=1500)
        if idx < total:
            print(f"  Progreso parcial: {idx}/{total}. Vuelve a correr este script para continuar.")
        else:
            print("  Ranking completo. Recuerda: reportar SIEMPRE por activo por separado "
                  "(ErG+GED__A9, ErG+GED__A22), nunca promediado -- ver hallazgo de A9 vs A22 "
                  "con 0/20 de overlap en el top-20 pese a ser quimicamente muy parecidos.")
    else:
        print(f"  (Salta este paso: no se encontro {excel_candidatos}. Ajusta la ruta a tu Excel real.)")

    print("\n" + "=" * 78)
    print("PASO 8 — FASE B (dengue, el target real) -- correr solo si FASE A convence")
    print("=" * 78)
    print("  df_dengue = compute_all_similarities(candidatos_500, DENGUE_ACTIVOS, "
          "theta=theta_final, por_activo=True)")
    print("  # Idealmente: repetir Secciones 6b/6c con activos reales de dengue una vez")
    print("  # que existan mas que 3 compuestos confirmados (ver matriz de literatura).")

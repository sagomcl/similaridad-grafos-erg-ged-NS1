"""
=============================================================================
ged_mejorado.py  --  Distancia de edicion de grafos (GED) para grafos ErG,
                     con emparejamiento de nodos que SI tiene en cuenta las
                     aristas (distancias entre grupos farmacoforicos).
=============================================================================

QUE RESUELVE ESTE ARCHIVO
-------------------------
La version original del GED (`pipeline_final_completo.ged_similarity`) decide
el emparejamiento de nodos con el algoritmo hungaro mirando UNICAMENTE el
costo de sustitucion de tipos (HBD, HBA, Aromatico, Hidrofobo). Cuando hay
nodos del mismo tipo (por ejemplo dos anillos aromaticos, o dos donadores)
la matriz de costos tiene empates, y el hungaro devuelve cualquiera de los
emparejamientos optimos: a veces cruza los nodos (empareja el aromatico de
la izquierda de una molecula con el de la derecha de la otra). Solo DESPUES
se calcula el costo de aristas sobre ese emparejamiento ya fijado, que puede
ser muy malo.

Consecuencia concreta medida en este proyecto: A9 y A22, cuyos grafos ErG son
IDENTICOS (6 nodos y 14 aristas iguales), daban similitud 0,43 en vez de 1,00.

COMO LO ARREGLA (estandar en GED bipartito, Riesen & Bunke 2009)
---------------------------------------------------------------
1) Matriz de costos "aumentada": al costo de emparejar el nodo i con el nodo j
   se le suma la mitad del costo de alinear sus estructuras locales de aristas
   (las distancias de las aristas incidentes a i contra las incidentes a j).
   Asi el hungaro ya no decide a ciegas entre nodos del mismo tipo: prefiere
   el que ademas esta a distancias parecidas de sus vecinos.
2) Refinamiento local: se intercambian pares de nodos del emparejamiento
   mientras el costo TOTAL (nodos + aristas) siga bajando.
3) Se prueban los dos puntos de partida (matriz de solo nodos y matriz
   aumentada) y se conserva el mejor, de modo que el costo nunca sale peor
   que el de la version anterior.

DOS FORMAS DE USARLO
--------------------
a) Reemplazo explicito (recomendado en codigo nuevo):
       from ged_mejorado import ged_similarity_v2
       sim = ged_similarity_v2(G1, G2, theta)

b) Parche en caliente (asi lo hacen los scripts de este repositorio, para no
   duplicar el resto del pipeline):
       import ged_mejorado
       ged_mejorado.parchear()                      # ANTES del import siguiente
       from pipeline_final_completo import ged_similarity   # ya es la v2

LIMITACION IMPORTANTE (documentarla en la tesis)
------------------------------------------------
El paso 2 es una busqueda local (greedy): su resultado puede depender del
orden en que vienen numerados los nodos. Por eso dos grafos identicos pero
con nodos numerados distinto (como A9 y A22) pueden terminar con similitudes
que difieren en el tercer decimal frente a una misma molecula candidata. Es
un optimo local, no el GED exacto -- que es NP-dificil.
"""
import math
import numpy as np
from scipy.optimize import linear_sum_assignment

# -----------------------------------------------------------------------------
# Los cuatro tipos de nodo de un grafo ErG. El orden define las filas/columnas
# de la matriz de costos de sustitucion.
#
# >>> PARA EXTRAPOLAR: si se quieren mas tipos de nodo (por ejemplo separar
#     "cargado positivo" y "cargado negativo", que es parte de la propuesta
#     original de Stiefl et al. 2006), hay que: (1) agregarlos a esta lista,
#     (2) ampliar la matriz de _cost_matrix con los costos nuevos, y
#     (3) ampliar LABEL_MAP y PRIORITY en pipeline_final_completo.py.
# -----------------------------------------------------------------------------
TYPES = ["Aromatico", "Hidrofobo", "HBA", "HBD"]
IDX = {t: i for i, t in enumerate(TYPES)}


def _cost_matrix(theta):
    """Matriz 4x4 de costo de SUSTITUIR un tipo de nodo por otro.

    Leer asi: C[i, j] = cuanto cuesta decir que un nodo de tipo TYPES[i] de la
    molecula 1 "es" un nodo de tipo TYPES[j] de la molecula 2.

      - Diagonal = 0: cambiar un aromatico por un aromatico no cuesta nada.
      - a = w_aromatico_hidrofobo: confundir aromatico con hidrofobo. Barato,
        porque ambos son grupos apolares que ocupan volumen.
      - h = w_hba_hbd: confundir aceptor con donador de puente de hidrogeno.
        Son el mismo "idioma" (polar) pero con el papel invertido.
      - p = w_polar_apolar: confundir un grupo polar (HBA/HBD) con uno apolar
        (aromatico/hidrofobo). Es el error quimicamente mas grave, por eso en
        los pesos calibrados suele ser el mas alto.

    >>> PARA EXTRAPOLAR a otra proteina: estos tres numeros son el "criterio
        quimico" del metodo. No hace falta reescribir codigo para cambiarlos:
        se calibran con recalibrar_v2.py sobre los activos de la nueva diana y
        quedan guardados en un JSON (ver theta_final_v2.json).
    """
    a = theta["w_aromatico_hidrofobo"]
    h = theta["w_hba_hbd"]
    p = theta["w_polar_apolar"]
    return np.array([
        #        Aromatico  Hidrofobo  HBA   HBD
        [0.0, a,   p,   p],     # Aromatico
        [a,   0.0, p,   p],     # Hidrofobo
        [p,   p,   0.0, h],     # HBA
        [p,   p,   h,   0.0],   # HBD
    ])


def _edges(G):
    """Aristas del grafo como diccionario {(u, v) ordenado: distancia}.

    La "distancia" de una arista ErG es el numero de enlaces del camino mas
    corto entre los dos grupos farmacoforicos (1 a 6). Se normaliza la clave
    con sorted() para que (2, 5) y (5, 2) sean la misma arista.
    """
    return {tuple(sorted((u, v))): int(d["dist"]) for u, v, d in G.edges(data=True)}


def _incident(n, E):
    """Para cada nodo, la lista de distancias de las aristas que salen de el.

    Es la "huella local" del nodo: un nodo que esta a 2, 3 y 5 enlaces de sus
    vecinos tiene huella [2, 3, 5]. Comparar huellas es lo que permite
    distinguir dos nodos del mismo tipo (paso 1 del arreglo).
    """
    inc = [[] for _ in range(n)]
    for (u, v), d in E.items():
        inc[u].append(d)
        inc[v].append(d)
    return inc


def _local_edge_cost(a, b, beta, indel):
    """Costo minimo de alinear dos huellas locales de aristas, a y b.

    Es un problema de asignacion en miniatura: cada distancia de `a` se empareja
    con una de `b` pagando beta * |diferencia|; si a una no le alcanza pareja
    (porque un nodo tiene mas vecinos que el otro) se paga `indel` (insertar o
    borrar una arista). Se resuelve otra vez con el hungaro.
    """
    if not a and not b:
        return 0.0
    k = max(len(a), len(b))
    M = np.full((k, k), indel, dtype=float)   # relleno = costo de insertar/borrar
    for x, da in enumerate(a):
        for y, db in enumerate(b):
            M[x, y] = beta * abs(da - db)
    r, c = linear_sum_assignment(M)
    return float(M[r, c].sum())


def _total_cost(p, Cn, n1, n2, E1, E2, beta, indel):
    """Costo TOTAL de edicion para un emparejamiento concreto `p`.

    `p` es una lista: p[i] = indice del nodo de G2 con el que se empareja el
    nodo i de G1 (los indices >= n2 son nodos "ficticios", es decir borrados).

    El costo tiene dos partes:
      1. Nodos: suma de Cn[i, p[i]] (sustituciones + inserciones/borrados).
      2. Aristas: para cada arista (u, v) de G1 se mira si la arista
         (p[u], p[v]) existe en G2.
           - si existe: se paga beta * |dist1 - dist2|  (estirar/encoger la
             distancia entre dos grupos farmacoforicos);
           - si no existe: se paga `indel` (borrar la arista).
         Al final se cobra `indel` por cada arista de G2 que quedo sin usar
         (hay que insertarla).
    """
    cost_nodes = sum(Cn[i, p[i]] for i in range(len(p)))
    used, cost_edges = set(), 0.0
    for (u, v), d1 in E1.items():
        u2, v2 = p[u], p[v]
        if u2 < n2 and v2 < n2:                      # ambos nodos tienen pareja real
            key = (u2, v2) if u2 < v2 else (v2, u2)
            if key in E2:
                cost_edges += beta * abs(d1 - E2[key])
                used.add(key)
                continue
        cost_edges += indel
    cost_edges += indel * sum(1 for k in E2 if k not in used)
    return cost_nodes + cost_edges


def _refine(p, Cn, n1, n2, E1, E2, beta, indel, max_iter=50):
    """Busqueda local: intercambia pares de asignaciones mientras mejore.

    Recorre todos los pares (i, k), prueba el intercambio p[i] <-> p[k] y lo
    conserva solo si el costo total baja. Repite hasta que una pasada completa
    no encuentre ninguna mejora (o hasta max_iter pasadas).

    Es el paso que corrige los emparejamientos "cruzados" que el hungaro de
    solo-nodos dejaba pasar por empate.
    """
    p = list(p)
    best = _total_cost(p, Cn, n1, n2, E1, E2, beta, indel)
    N = len(p)
    for _ in range(max_iter):
        mejora = False
        for i in range(N):
            for k in range(i + 1, N):
                p[i], p[k] = p[k], p[i]                      # probar el intercambio
                c = _total_cost(p, Cn, n1, n2, E1, E2, beta, indel)
                if c < best - 1e-12:
                    best, mejora = c, True                   # se queda
                else:
                    p[i], p[k] = p[k], p[i]                  # se deshace
        if not mejora:
            break
    return p, best


def ged_detallado(G1, G2, theta, refine=True):
    """Calcula el GED entre dos grafos ErG y devuelve todo el detalle.

    Devuelve (similitud, mapping, ged, ged_normalizado):
      - similitud: exp(-GED_normalizado), en (0, 1]. 1 = grafos identicos.
      - mapping: {nodo de G1: nodo de G2} del emparejamiento elegido. Sirve
        para explicar el resultado a un quimico (ver comparar_A9_A22.py).
      - ged: costo de edicion crudo, en las unidades de theta.
      - ged_normalizado: ged / (n1 + n2 + m1 + m2).

    POR QUE SE NORMALIZA: sin normalizar, dos moleculas grandes acumulan mas
    costo que dos pequenas solo por tener mas nodos y aristas, y el ranking
    quedaria sesgado por tamano. Se divide por el total de elementos de los dos
    grafos (nodos + aristas), que es el maximo numero de operaciones de edicion
    posibles, y asi el costo queda en una escala comparable entre pares.

    POR QUE exp(-x) Y NO 1-x: exp(-x) convierte un costo en [0, inf) en una
    similitud en (0, 1] de forma suave y monotona, sin necesidad de recortar
    valores, y penaliza mucho las primeras diferencias y cada vez menos las
    siguientes (dos moleculas ya muy distintas apenas se distinguen entre si,
    que es justo lo que interesa en un cribado por similitud).
    """
    C_SUB = _cost_matrix(theta)
    node_indel, beta, indel = theta["node_indel"], theta["edge_beta"], theta["edge_indel"]
    n1, n2 = G1.number_of_nodes(), G2.number_of_nodes()
    N = max(n1, n2)                      # la matriz se cuadra con nodos ficticios
    E1, E2 = _edges(G1), _edges(G2)
    inc1, inc2 = _incident(n1, E1), _incident(n2, E2)

    # ---- (1) Costo de nodos "puro": solo tipos. Es el que se usa para
    #      evaluar el costo final, igual que en la version original.
    Cn = np.zeros((N, N))
    for i in range(n1):
        li = IDX[G1.nodes[i]["label"]]
        for j in range(n2):
            Cn[i, j] = C_SUB[li, IDX[G2.nodes[j]["label"]]]
    for i in range(n1, N):               # filas ficticias = insertar nodo en G1
        Cn[i, :n2] = node_indel
    for i in range(n1):                  # columnas ficticias = borrar nodo de G1
        Cn[i, n2:] = node_indel

    # ---- (2) Costo "aumentado": tipos + huella local de aristas. Solo se usa
    #      para ELEGIR el emparejamiento, no para puntuarlo. El factor 0,5
    #      evita contar dos veces cada arista (una por cada uno de sus extremos).
    Ca = Cn.copy()
    for i in range(n1):
        for j in range(n2):
            Ca[i, j] += 0.5 * _local_edge_cost(inc1[i], inc2[j], beta, indel)
        for j in range(n2, N):
            Ca[i, j] += 0.5 * indel * len(inc1[i])
    for i in range(n1, N):
        for j in range(n2):
            Ca[i, j] += 0.5 * indel * len(inc2[j])

    # ---- (3) Dos candidatos (hungaro sobre Cn y sobre Ca), cada uno refinado;
    #      se conserva el de menor costo total real.
    candidatos = []
    for M in (Cn, Ca):
        r, c = linear_sum_assignment(M)          # <-- ALGORITMO HUNGARO
        p = [0] * N
        for i, j in zip(r, c):
            p[int(i)] = int(j)
        if refine:
            p, costo = _refine(p, Cn, n1, n2, E1, E2, beta, indel)
        else:
            costo = _total_cost(p, Cn, n1, n2, E1, E2, beta, indel)
        candidatos.append((costo, p))
    costo, p = min(candidatos, key=lambda x: x[0])

    denom = n1 + n2 + len(E1) + len(E2)
    n_ged = costo / denom if denom > 0 else 0.0
    mapping = {i: p[i] for i in range(n1) if p[i] < n2}
    return math.exp(-n_ged), mapping, costo, n_ged


def ged_similarity_v2(G1, G2, theta, verbose=False, refine=True):
    """Reemplazo directo de pipeline_final_completo.ged_similarity.

    Devuelve solo la similitud (un numero en (0, 1]), que es lo que necesitan
    el ranking y la calibracion. Con verbose=True imprime el emparejamiento,
    util para revisar un par concreto de moleculas con la quimica.
    """
    sim, mapping, ged, n_ged = ged_detallado(G1, G2, theta, refine=refine)
    if verbose:
        print(f"  emparejamiento G1->G2: {mapping}")
        print(f"  GED={ged:.3f}  GED_norm={n_ged:.3f}  similitud={sim:.3f}")
    return sim


def parchear():
    """Sustituye `ged_similarity` dentro de pipeline_final_completo por la v2.

    Hay que llamarla ANTES de hacer `from pipeline_final_completo import ...`,
    porque ese import copia la referencia a la funcion en el momento en que se
    ejecuta: si se parchea despues, el script seguiria usando la version vieja.

    Nota: pipeline_v3.py ya importa directamente la v2 y no necesita parche.
    """
    import pipeline_final_completo as pf
    default = pf.DEFAULT_THETA

    def ged_similarity(G1, G2, theta=default, verbose=False):
        return ged_similarity_v2(G1, G2, theta, verbose=verbose)

    pf.ged_similarity = ged_similarity
    return ged_similarity

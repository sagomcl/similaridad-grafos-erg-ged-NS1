"""
=============================================================================
ns1_utils.py  --  (1) criterios farmacoforicos del sitio de union de NS1 y
                  (2) metricas de evaluacion (AUC, enriquecimiento, bootstrap)
=============================================================================

Este archivo no depende de RDKit a proposito: recibe el grafo ErG ya construido
(nodos, aristas y coordenadas) como listas y arreglos de numpy. Asi se puede
probar y reutilizar sin instalar nada de quimioinformatica.

PARTE 1 -- QUE SE LE PIDE A UNA MOLECULA PARA PODER INHIBIR NS1
---------------------------------------------------------------
El bolsillo de NS1 donde se unen A9 y A22 (el sitio CPSF30, descrito por
Kleinpeter et al. 2018) es un canal alargado y mayoritariamente hidrofobo, con
dos cavidades en los extremos donde encajan anillos aromaticos, y con puntos
de anclaje polar en el medio. De ahi salen los cinco criterios que implementa
`criterios_ns1`:

  1. COMPOSICION. Al menos 2 nodos aromaticos, >= 1 hidrofobico, >= 1 donador
     (HBD) y >= 1 aceptor (HBA). Es decir: la molecula debe tener con que
     llenar el canal apolar Y con que formar puentes de hidrogeno.
  2. CONECTIVIDAD. El grafo ErG debe quedar conectado. Como solo se dibujan
     aristas entre grupos separados por 1 a 6 enlaces, un grafo desconectado
     significa que hay grupos farmacoforicos a mas de 6 enlaces de todos los
     demas: la molecula es demasiado laxa y sus grupos no pueden ocupar el
     bolsillo a la vez.
  3. AROMATICOS EN LOS EXTREMOS. La separacion 3D entre los dos aromaticos mas
     alejados debe ser al menos `ratio_extremos` (70 % por defecto) de la
     separacion maxima entre cualquier par de grupos. Traduce la idea de que
     los anillos deben estar en las puntas de la molecula, no apilados en el
     centro.
  4. TAMANO. El largo de la molecula (distancia maxima entre dos atomos
     pesados) debe superar `tam_min` angstroms, para que alcance a cubrir el
     canal.
  5. AJUSTE DE DISTANCIAS. Por cada arista se mide (distancia 3D entre los
     centros de los dos grupos) / (numero de enlaces que los separa). Si ese
     cociente cae entre 1,5 y 2,0 A por enlace, la cadena esta razonablemente
     extendida; si es mucho menor, esta plegada sobre si misma. Se reporta la
     fraccion de aristas que cumplen.

>>> PARA EXTRAPOLAR A OTRA PROTEINA: este es el bloque que hay que cambiar.
    Los criterios 1 a 5 son una descripcion geometrica del bolsillo de NS1.
    Para otra diana se ajustan los umbrales en PARAMS (o por linea de comandos
    en evaluar_erg_ns1.py: --tam-min, --ratio-extremos, --A-por-enlace) y, si
    el bolsillo tiene otra forma (por ejemplo una cavidad compacta en vez de un
    canal), se reescribe `criterios_ns1` con las condiciones que correspondan.
    El resto del pipeline (ErG, GED, hungaro, calibracion) no cambia.

PARTE 2 -- METRICAS
-------------------
  - auc:    area bajo la curva ROC, calculada a mano (equivale al estadistico
            de Mann-Whitney): probabilidad de que un activo tomado al azar
            reciba un puntaje mayor que un inactivo tomado al azar.
            0,5 = azar; 1,0 = separacion perfecta.
  - ef:     factor de enriquecimiento. Cuantas veces mas activos hay en el
            top `frac` de la lista que los que habria ordenando al azar. Es la
            metrica que de verdad le importa a un cribado: si EF20 = 4, en el
            primer 20 % de la lista hay 4 veces mas activos que por azar.
  - auc_ci: intervalo de confianza del 95 % para el AUC por bootstrap. Con
            muestras de 5 a 30 activos, reportar un AUC sin intervalo es
            enganoso; esta funcion es la que permite decir si la diferencia
            entre dos metodos es real o es ruido.
  - rango_norm: convierte puntajes a rangos en [0, 1]. Se usa para fusionar
            dos metodos de escalas distintas (ErG+GED y Tanimoto) sumando
            rangos en vez de valores crudos.
"""
import itertools
import numpy as np

# -----------------------------------------------------------------------------
# Umbrales de los criterios geometricos. Valores por defecto deducidos del
# bolsillo CPSF30 de NS1 de influenza.
#
# >>> PARA EXTRAPOLAR: cambiar estos tres numeros (o pasarlos en el argumento
#     `p` de criterios_ns1) es la forma mas rapida de adaptar el filtro a otra
#     diana sin tocar codigo.
# -----------------------------------------------------------------------------
PARAMS = {
    "tam_min": 15.0,            # largo minimo de la molecula, en angstroms
    "ratio_extremos": 0.70,     # que tan "en las puntas" deben estar los aromaticos
    "A_por_enlace": (1.5, 2.0), # rango aceptable de angstroms por enlace
}


def es_conexo(n, aristas):
    """¿El grafo de n nodos con esas aristas es conexo? (union-find simple).

    Se usa para el criterio 2. Implementado a mano para no depender de networkx
    en este archivo.
    """
    if n == 0:
        return False
    padre = list(range(n))

    def raiz(x):
        while padre[x] != x:
            padre[x] = padre[padre[x]]   # compresion de camino
            x = padre[x]
        return x

    for i, j, _ in aristas:
        padre[raiz(i)] = raiz(j)
    return len({raiz(i) for i in range(n)}) == 1


def largo_max(coords):
    """Distancia euclidiana maxima entre dos puntos de `coords` (Nx3).

    Con las coordenadas de los atomos pesados del conformero 3D, esto es el
    "largo" de la molecula en angstroms (criterio 4).
    """
    coords = np.asarray(coords, float)
    if len(coords) < 2:
        return 0.0
    return float(np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2).max())


def criterios_ns1(nodos, aristas, coords_pesados, p=None):
    """Evalua los cinco criterios sobre una molecula ya convertida a ErG.

    Parametros
    ----------
    nodos : lista de (etiqueta, centro_xyz)
        Etiqueta es "Aromatico" / "Hidrofobo" / "HBA" / "HBD"; centro_xyz es el
        centroide 3D de los atomos que se colapsaron en ese nodo.
    aristas : lista de (i, j, n_enlaces)
        n_enlaces = longitud en enlaces del camino mas corto entre los grupos.
    coords_pesados : arreglo Nx3
        Coordenadas 3D de todos los atomos pesados (sin hidrogenos).
    p : dict opcional
        Umbrales que sobreescriben PARAMS.

    Devuelve un diccionario con el detalle de cada criterio. Las claves
    importantes son:
      - "estructura": True si cumple composicion + conectado + extremos.
        Es el filtro "duro" de forma farmacoforica.
      - "tam_ok": True si supera el largo minimo. Se deja aparte porque
        depende mucho del conformero generado (ver nota abajo).
      - "ajuste_dist": fraccion de aristas con angstroms/enlace en rango. Es
        un numero continuo, util como variable descriptiva mas que como filtro.

    NOTA SOBRE EL TAMANO: un unico conformero puede salir plegado y subestimar
    el largo real. Por eso en evaluar_erg_ns1.py existe `largo_multiconf`, que
    genera 30 conformeros y toma el maximo, y por eso el criterio de tamano
    solo entra al filtro con la bandera --exigir-tamano.
    """
    p = {**PARAMS, **(p or {})}
    labels = [l for l, _ in nodos]

    # --- criterio 1: composicion farmacoforica minima
    comp = (labels.count("Aromatico") >= 2 and labels.count("Hidrofobo") >= 1
            and labels.count("HBD") >= 1 and labels.count("HBA") >= 1)

    # --- criterio 2: el grafo no se parte en pedazos
    conectado = es_conexo(len(nodos), aristas)

    # --- criterio 3: los dos aromaticos mas separados deben marcar los extremos
    centros = [np.asarray(c, float) for _, c in nodos]
    d_max_all, d_max_ar = 0.0, 0.0
    for a, b in itertools.combinations(range(len(nodos)), 2):
        d = float(np.linalg.norm(centros[a] - centros[b]))   # distancia 3D real
        d_max_all = max(d_max_all, d)
        if labels[a] == "Aromatico" and labels[b] == "Aromatico":
            d_max_ar = max(d_max_ar, d)
    ratio_ext = d_max_ar / d_max_all if d_max_all > 0 else 0.0
    extremos = labels.count("Aromatico") >= 2 and ratio_ext >= p["ratio_extremos"]

    # --- criterio 4: largo de la molecula
    largo = largo_max(coords_pesados)

    # --- criterio 5: angstroms por enlace (cadena extendida vs. plegada)
    ratios = [float(np.linalg.norm(centros[i] - centros[j])) / d for i, j, d in aristas if d > 0]
    lo, hi = p["A_por_enlace"]
    ajuste = float(np.mean([lo <= r <= hi for r in ratios])) if ratios else 0.0

    estructura = bool(comp and conectado and extremos)
    return {
        "composicion": bool(comp), "conectado": bool(conectado),
        "ratio_extremos": ratio_ext, "extremos": bool(extremos),
        "estructura": estructura,
        "largo_A": largo, "tam_ok": bool(largo > p["tam_min"]),
        "A_por_enlace_medio": float(np.mean(ratios)) if ratios else 0.0,
        "ajuste_dist": ajuste,
        "n_arom": labels.count("Aromatico"), "n_hidrof": labels.count("Hidrofobo"),
        "n_hbd": labels.count("HBD"), "n_hba": labels.count("HBA"),
    }


# ============================================================ metricas
def auc(y, s):
    """AUC-ROC calculada como estadistico de Mann-Whitney.

    y: etiquetas 0/1 (1 = activo). s: puntajes (mayor = mas parecido a activo).
    Cuenta, sobre todos los pares (activo, inactivo), la fraccion en que el
    activo tiene mayor puntaje; los empates cuentan medio. Equivale al area
    bajo la curva ROC y no necesita ordenar umbrales.
    """
    y = np.asarray(y); s = np.asarray(s, float)
    pos, neg = s[y == 1], s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    mayor = (pos[:, None] > neg[None, :]).sum() + 0.5 * (pos[:, None] == neg[None, :]).sum()
    return float(mayor / (len(pos) * len(neg)))


def ef(y, s, frac):
    """Factor de enriquecimiento en el primer `frac` de la lista ordenada.

    EF = (proporcion de activos en el top k) / (proporcion de activos global).
    EF = 1 es lo mismo que ordenar al azar; EF = 5 es cinco veces mejor.
    """
    y = np.asarray(y); s = np.asarray(s, float)
    if y.sum() == 0:
        return float("nan")
    k = max(1, int(round(frac * len(y))))
    idx = np.argsort(-s, kind="stable")[:k]
    return float((y[idx].sum() / k) / y.mean())


def auc_ci(y, s, n=300, seed=0):
    """Intervalo de confianza del 95 % del AUC por bootstrap estratificado.

    Remuestrea con reemplazo los activos y los inactivos por separado `n`
    veces, recalcula el AUC y devuelve los percentiles 2,5 y 97,5. Con pocos
    activos el intervalo sale ancho, y eso es precisamente la informacion que
    hay que reportar.
    """
    y = np.asarray(y); s = np.asarray(s, float)
    pos, neg = s[y == 1], s[y == 0]
    if len(pos) < 2 or len(neg) < 2:
        return (float("nan"), float("nan"))
    rng = np.random.RandomState(seed)
    vals = []
    for _ in range(n):
        p = rng.choice(pos, len(pos)); q = rng.choice(neg, len(neg))
        vals.append(auc(np.r_[np.ones(len(p)), np.zeros(len(q))], np.r_[p, q]))
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))


def rango_norm(s):
    """Convierte puntajes en rangos normalizados a [0, 1].

    Sirve para fusionar metodos con escalas distintas: en vez de sumar una
    similitud ErG+GED (que vive cerca de 0,8-1,0) con un Tanimoto (que vive
    cerca de 0,1-0,5), se suman sus rangos, que si son comparables.
    """
    s = np.asarray(s, float)
    if len(s) < 2:
        return np.zeros(len(s))
    return np.argsort(np.argsort(s, kind="stable"), kind="stable") / (len(s) - 1)

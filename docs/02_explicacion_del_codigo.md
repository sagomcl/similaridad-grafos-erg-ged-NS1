# 02 — Explicación del código, archivo por archivo

Antes de leer esto conviene leer [`01_explicacion_del_metodo.md`](01_explicacion_del_metodo.md),
que explica los conceptos. Aquí se explica **qué hace cada archivo, cada función y qué
hay que cambiar** para adaptar el método a otra proteína. Los lugares extrapolables
están marcados con ▶.

## Mapa del repositorio

| Archivo | Qué es | ¿Se toca para extrapolar? |
|---|---|---|
| `pipeline_final_completo.py` | **núcleo**: ErG, GED, referencias, calibración, ranking | sí, varias constantes |
| `ged_mejorado.py` | GED corregido (emparejamiento que mira las aristas) | solo si se agregan tipos de nodo |
| `ns1_utils.py` | criterios geométricos del bolsillo + métricas | **sí, es el bloque específico de la diana** |
| `pipeline_v3.py` | misma funcionalidad, endurecida y sin fuga de datos | igual que el núcleo |
| `recalibrar_v2.py` | calibra y valida θ → `theta_final_v2.json` | sí, el bloque de datos |
| `correr_ranking.py` | aplica el método a las 33 047 candidatas | sí, rutas y referencias |
| `evaluar_erg_ns1.py` | evaluación justa frente a Tanimoto | sí, criterios y referencias |
| `comparar_A9_A22.py` | análisis explicativo de un par de moléculas | sí, el par a comparar |
| `dibujar_top20.py` | figura con las estructuras del top‑20 | no |
| `calibracion_cluster.py`, `calibracion_pubchem_real.py` | scripts históricos | no (no reproducibles) |

Dependencias: `rdkit`, `numpy`, `scipy`, `scikit-learn`, `pandas`, `networkx`,
`matplotlib`, `openpyxl`. Python 3.10–3.12.

---

# 1. `pipeline_final_completo.py` — el núcleo

Es el archivo que importan todos los demás. Está organizado en ocho secciones.

## Sección 1 — De SMILES a grafo ErG

### Constantes del módulo

```python
fdef_path   = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
feat_factory = ChemicalFeatures.BuildFeatureFactory(fdef_path)
```
`BaseFeatures.fdef` es el catálogo de definiciones farmacofóricas que trae RDKit: un
archivo de reglas SMARTS que dice "un nitrógeno con hidrógeno es un donador", "un
anillo aromático es `Aromatic`", etc. Se construye la *feature factory* **una sola vez
al importar el módulo**, para no recompilar esas reglas en cada una de las 33 047
moléculas. Es una optimización importante.

```python
LABEL_MAP = {"Donor": "HBD", "Acceptor": "HBA",
             "Aromatic": "Aromatico", "Hydrophobe": "Hidrofobo"}
```
▶ **Extrapolable.** Traduce las familias de RDKit a las cuatro etiquetas del ErG.
Cualquier familia que no esté en este diccionario (RDKit también define `PosIonizable`,
`NegIonizable`, `ZnBinder`...) simplemente **no genera nodo**. Si la diana nueva depende
de cargas —por ejemplo un bolsillo con un residuo de aspartato que exige un catión—, hay
que añadir `"PosIonizable": "Cation"` aquí y además ampliar `PRIORITY`, `TYPES` y la
matriz de costos.

```python
PRIORITY = ["HBD", "HBA", "Aromatico", "Hidrofobo"]
```
Resuelve el problema de los solapamientos (un átomo en varias familias). Cada átomo
recibe **un solo** tipo: el primero de esta lista que le aplique. El orden va de la
propiedad más específica a la más genérica. ▶ Cambiar el orden cambia la topología de
todos los grafos; si se cambia, hay que recalibrar θ.

```python
TYPES = ["Aromatico", "Hidrofobo", "HBA", "HBD"]
IDX   = {t: i for i, t in enumerate(TYPES)}
```
Orden canónico: define las filas y columnas de la matriz de costos. Tiene que ser el
mismo en `pipeline_final_completo.py` y en `ged_mejorado.py`.

### `smiles_to_3d_mol(smiles, seed=0xF00D)`

Convierte texto SMILES en molécula con geometría 3D. Cuatro pasos:

1. `Chem.MolFromSmiles` — interpreta el SMILES y construye el grafo molecular. Devuelve
   `None` si el SMILES está mal escrito → se lanza `ValueError`.
2. `Chem.AddHs` — agrega los hidrógenos explícitos. Son imprescindibles para que la
   geometría 3D tenga sentido (los ángulos de enlace dependen de ellos), aunque después
   no sean nodos del ErG.
3. `AllChem.ETKDGv3()` + `EmbedMolecule` — genera las coordenadas 3D. **`randomSeed` fija
   la semilla**: la misma molécula produce siempre el mismo conformero. Sin esto nada
   sería reproducible.
4. `UFFOptimizeMolecule` — relaja la geometría con el campo de fuerzas UFF.

Devuelve **dos** objetos con papeles distintos:
- `mol` (3D, con hidrógenos): de aquí salen los centros de los nodos y las coordenadas
  para los criterios geométricos.
- `mol2d` (sin hidrógenos): de aquí salen los *fingerprints* y las propiedades
  fisicoquímicas, que se calculan sobre la topología.

Es **el paso más lento** del pipeline (decenas de ms por molécula) y la única fuente
habitual de fallos: moléculas muy grandes o con estereoquímica imposible hacen que
`EmbedMolecule` devuelva −1 y se lance `RuntimeError`. En el ranking final fallaron
**17 de 33 047** (registradas en `resultados/ranking_33k_thetav2_errores.csv`).

### `build_erg(mol)` — la función central

Devuelve un `networkx.Graph` donde cada nodo tiene `label`, `atom_ids` y `center`, y
cada arista tiene `dist`. Cuatro pasos:

**Paso A — un tipo por átomo.**
```python
for f in feats:
    if fam not in LABEL_MAP: continue
    for a in f.GetAtomIds():
        if cur is None or PRIORITY.index(label) < PRIORITY.index(cur):
            atom_type[a] = label
```
Se queda la etiqueta de mayor prioridad. Los átomos que no caen en ninguna de las
cuatro familias no entran en `atom_type` y no generan nodo: son "pegamento" estructural.

**Paso B — colapsar contiguos del mismo tipo.** Búsqueda en anchura (`deque`) sobre
`bond_graph`, agrupando átomos que (i) tienen el mismo tipo y (ii) están unidos por
enlaces. Un benceno (6 carbonos aromáticos contiguos) se vuelve **un** nodo.
`comp` es el conjunto de átomos del nodo y queda guardado en `atom_ids`: **así no se
pierde la información atómica**.

**Paso C — geometría del nodo.**
```python
center = np.array([list(conf.GetAtomPosition(i)) for i in comp]).mean(axis=0)
```
Centroide 3D de los átomos del nodo. No se usa en el GED; sí en los criterios
geométricos de `ns1_utils.py`.

**Paso D — aristas por distancia topológica.** La función interna
`shortest_bond_distance(atomsA, atomsB)` hace un BFS multiorigen y devuelve el **número
de enlaces** del camino más corto entre los dos grupos. Se crea arista solo si:
```python
if 1 <= d <= 6:
    G.add_edge(i, j, dist=int(d))
```
▶ **Extrapolable (el rango 1–6).** Por qué en enlaces y no en ángstroms: la distancia
topológica no depende del conformero, así que el grafo y el ranking son estables. Por
qué el tope de 6: más allá hay tantos grados de libertad que la distancia real en el
espacio deja de estar determinada, y además acota el número de aristas. Para bolsillos
más grandes puede ampliarse (por ejemplo a 8), pero **hay que recalibrar θ** porque
cambia el número de aristas y con ello la escala del costo.

### `erg_visualization_data(smiles)`

Exporta todo el paso a paso (átomos con coordenadas 2D, enlaces, tipo de cada átomo,
nodos ErG con posición, aristas con distancia) en un diccionario listo para la web.
Alimenta `figuras/explorador_erg_3d.html` y las figuras didácticas de la presentación.
Ojo: aquí los tipos se recalculan sobre `mol2d` porque las coordenadas exportadas son de
dibujo, 2D. Para el método se usa siempre `build_erg(mol3d)`.

## Sección 2 — GED con algoritmo húngaro

### `DEFAULT_THETA`

Los seis pesos. ▶ **Extrapolable**, pero **no editándolos a mano**: estos valores son una
conjetura inicial y **no** son los usados en los resultados. Los finales están en
`resultados/theta_final_v2.json` y salen de `recalibrar_v2.py`.

| Peso | Qué penaliza |
|---|---|
| `w_aromatico_hidrofobo` | confundir aromático con hidrofóbico |
| `w_hba_hbd` | confundir aceptor con donador |
| `w_polar_apolar` | confundir polar con apolar (el error más grave) |
| `node_indel` | insertar o borrar un nodo entero |
| `edge_beta` | por cada enlace de diferencia entre aristas emparejadas |
| `edge_indel` | insertar o borrar una arista |

### `cost_matrix_from_theta(theta)`

Matriz 4×4 simétrica con diagonal cero (sustituir un tipo por sí mismo es gratis). El
orden de filas/columnas es `TYPES`.

### `ged_similarity(G1, G2, theta)` — **versión 1, conservada por trazabilidad**

Esta es la versión con el defecto. Se mantiene para poder reproducir los números
antiguos y documentar la corrección. **Los resultados reportados usan la de
`ged_mejorado.py`.**

Pasos:
1. Matriz `Cn` de tamaño N×N con N = max(n₁, n₂). Bloque real = costos de sustitución
   de tipos; filas/columnas de relleno = `node_indel` (nodos ficticios).
2. `linear_sum_assignment(Cn)` ← **el algoritmo húngaro**. Óptimo en O(N³).
3. Con el emparejamiento fijado, se recorren las aristas de G1: si la correspondiente
   existe en G2 se cobra `edge_beta × |Δdist|`; si no, `edge_indel`. Al final,
   `edge_indel` por cada arista de G2 que nadie usó.
4. `GED_norm = GED / (n₁+n₂+m₁+m₂)`; `similitud = exp(−GED_norm)`.

**El defecto:** el paso 1 mira solo tipos, así que con nodos del mismo tipo la matriz
empata y el húngaro puede devolver un emparejamiento cruzado. Síntoma: A9 vs A22 = 0,43
en lugar de 1,00.

## Sección 3 — Moléculas de referencia

```python
INFLUENZA_ACTIVOS = {"A9":  "O=C(NCCc1ccccc1)c1ccccc1O",
                     "A22": "Nc1ccccc1C(=O)NCCc1ccccc1"}
```
▶ **El punto de extrapolación más importante.** Un diccionario `{nombre: SMILES}` con los
inhibidores conocidos de la diana nueva es todo lo que necesita el resto del pipeline.

A9 = 2‑hidroxi‑N‑fenetilbenzamida; A22 = 2‑amino‑N‑fenetilbenzamida. Difieren en **un
átomo** (–OH frente a –NH₂) y los dos se unen al bolsillo CPSF30 de NS1 de influenza
(Kleinpeter et al. 2018). Con la reducción ErG sus grafos resultan **idénticos**, porque
–OH y –NH₂ son los dos donadores: sobre influenza, por tanto, **las dos referencias aportan la
misma información al método**. Es un dato que hay que declarar en la tesis.

`INFLUENZA_BENCHMARK_NSC` / `INFLUENZA_IC50_UM`: cuatro compuestos de Basu et al. (2009)
con IC50 medido. Con n = 4 no da para una métrica, solo para inspección.

`DENGUE_ACTIVOS`: berberina, ácido betulínico y ácido ursólico (Ansari et al. 2025).
Son un alcaloide y dos triterpenos: químicamente muy distintos entre sí, así que **no
deben promediarse** — cada uno se usa como referencia por separado.

## Sección 4 — Métodos de comparación

- `tanimoto_bits(fp1, fp2)` — T = c/(a+b−c) sobre fingerprints de bits. La línea base
  estándar. Mide coincidencia de **subestructuras exactas**.
- `hamming_similarity_bits(fp1, fp2)` — 1 − fracción de bits distintos. Premia también
  los ceros coincidentes, lo que con vectores dispersos infla la similitud de cualquier
  par; medida pobre, incluida porque el anteproyecto la contemplaba.
- `erg_fp_cosine(fp1, fp2)` — Tanimoto generalizado para vectores continuos sobre el
  ErG‑fingerprint de RDKit. **La comparación más informativa**: usa la misma idea de
  reducción farmacofórica pero sin grafo ni GED, así que la diferencia entre ErG‑fp y
  ErG+GED aísla exactamente cuánto aporta el grafo.
- `_prep(smi)` — prepara una molécula **una sola vez** (grafo ErG + los dos fingerprints
  + peso molecular) para no repetir el embedding 3D en cada comparación.
- `compute_all_similarities(..., por_activo=False)` — tabla de similitudes. Con
  `por_activo=True` agrega una columna **por cada referencia por separado**
  (`ErG+GED__A9`, `ErG+GED__A22`, ...). **Es el modo que hay que usar**: promediar sobre
  activos de quimiotipos distintos diluye la señal. Las columnas sin sufijo son el
  promedio y se conservan solo por compatibilidad.

## Sección 5 — Decoys y calibración con curva ROC

### `construir_decoys_reales(...)`

Selecciona supuestos inactivos con el criterio tipo **DUD‑E**: moléculas de la biblioteca
propia con peso molecular en el rango de los activos (±`mw_tolerancia`) pero Tanimoto
máximo contra todos los activos por debajo de `tanimoto_max` (0,08). Se ordenan de menor
a mayor parecido y se toman las `n_decoys` más distintas.

Por qué importa: con decoys demasiado fáciles (moléculas enormes o diminutas) cualquier
método da AUC alto y la validación no vale nada; con decoys demasiado parecidos se están
etiquetando como inactivas moléculas que probablemente serían activas.

**Limitación:** son inactivos **presuntos**, nunca medidos. En este trabajo se reemplazaron
por los inactivos **confirmados** de AID 2326, mucho más sólido. Esta función queda para
el caso de una diana sin cribado disponible — que es exactamente la situación del dengue.

### `calibrar_theta_con_roc(...)`

Calibra θ maximizando el AUC y guarda tres cosas: el JSON de pesos, la curva ROC en PNG
y el CSV de puntajes con etiquetas (trazabilidad).

- El puntaje de una molécula es su similitud **media a los activos**, excluyéndose a sí
  misma (para no inflar el resultado de los propios activos).
- La función objetivo es **−AUC**.
- `differential_evolution` (Storn & Price 1997): optimizador global por poblaciones. Se
  usa **en vez de un método de gradiente** porque el AUC es una función **escalonada** de
  θ —solo cambia cuando dos moléculas se intercambian en el orden—, de modo que su
  gradiente es cero casi en todas partes.
- `bounds` fija el rango admisible de cada peso, en el mismo orden que las claves de
  `DEFAULT_THETA`. No son arbitrarios: evitan soluciones degeneradas (por ejemplo todos
  los costos en cero, que haría que todo se parezca a todo).

⚠️ **Advertencia:** esta función calibra **y evalúa sobre el mismo conjunto**. Su AUC
mide memorización, no generalización. Para resultados publicables usar
`calibrar_y_evaluar_train_test` o `recalibrar_v2.py`.

## Sección 6 — Datos reales

- `cargar_candidatos_excel(path, smiles_col, id_col)` — lee Excel/CSV → `{id: smiles}`.
  ▶ Para cribar otra biblioteca basta apuntar `path` al archivo nuevo.
- `_limpiar_smiles_sal(smi)` — los compuestos de PubChem suelen venir como sal
  (`"C...N.Cl"`, con el contraión separado por un punto). Se queda con el fragmento más
  grande. **Sin esta limpieza** el ErG se calcularía sobre dos fragmentos desconectados
  y la similitud no tendría sentido químico.
- `cargar_pubchem_ns1_assay(csv_path)` — carga AID 2326. El docstring incluye las
  instrucciones exactas para descargarlo (Download → *Deposited Data Table Only* → CSV).
  Solo se usan `Active` e `Inactive`; `Inconclusive` y `Unspecified` se descartan.
  ⚠️ `skiprows=[1,2,3,4]` asume las cuatro filas de metadatos que trae el CSV de
  PubChem; si el archivo no las tiene (como `AID2326_cache.csv`, generado desde el
  Excel), ese salto **borraría cuatro moléculas reales**. `pipeline_v3.py` filtra por
  contenido, que es lo correcto.

### `clusterizar_activos(activos, cutoff_distancia=0.6)`

Agrupa activos en quimiotipos con el algoritmo de **Butina (1999)**: calcula todas las
distancias (1 − Tanimoto), elige como centro la molécula con más vecinos dentro del
umbral, le asigna esos vecinos, los retira y repite. Es el estándar en
quimioinformática porque no exige fijar el número de grupos de antemano.

`cutoff_distancia = 0.6` ⇔ mismo cluster si Tanimoto > 0,4. ▶ **Extrapolable**: bajarlo
da quimiotipos más estrictos; subirlo, más laxos. Tiene efecto directo sobre el AUC que
se va a reportar, así que conviene dejarlo fijo y **declararlo**.

Resultado en esta tesis: **1 073 activos → 737 clusters, 588 de ellos singletons**. Ese
número es la explicación de por qué el AUC agregado daba 0,49.

### `calibrar_y_evaluar_train_test(...)`

La versión honesta de la calibración. Parte activos y decoys en entrenamiento/prueba,
calibra **solo** con entrenamiento y reporta AUC en ambos.

Detalle que importa: el puntaje de una molécula de **prueba** se calcula como su
similitud media a los activos de **entrenamiento**, nunca a los de prueba. Así se
reproduce la situación real de un cribado.

La diferencia entre `auc_train` y `auc_test` **es** la medida del sobreajuste: con n = 2
activos, 1,00 frente a 0,49.

## Sección 7 — Ranking

- `rankear_candidatos(...)` — ranking de un diccionario de candidatos con los cuatro
  métodos. Imprime cuántas moléculas comparten el top‑N de ErG+GED y el de
  Morgan+Tanimoto: la forma más directa de ver si el método aporta algo distinto.
- `rankear_excel_en_bloques(...)` — procesa un bloque de `batch_size` moléculas y
  **agrega** filas al CSV. Patrón de *checkpoint*: el JSON guarda
  `{"next_idx":…, "total":…}`; al arrancar, si existe se retoma desde ahí, y si no se
  empieza de cero borrando el CSV anterior. Hay que llamarla repetidamente (lo hace
  `correr_ranking.py` en un `while`) hasta que `next_idx == total`.
  Las referencias se preparan **una** vez y se reutilizan en las 33 047 comparaciones.
  ⚠️ `except Exception: continue` descarta en silencio las moléculas que fallan;
  `pipeline_v3.py` las registra en un CSV de errores — de ahí sale el conteo exacto de 17.

## Sección 8 — `if __name__ == "__main__"`

Demostración histórica de los ocho pasos. **No es la forma de reproducir los resultados**:
apunta a rutas que ya no existen (`merged_500_con_props.csv`, `data/Moleculas.xlsx`,
`AID_2326_datatable.csv`). Se conserva sin cambios porque documenta el orden lógico.

---

# 2. `ged_mejorado.py` — el GED corregido

Es la versión que produce **todos los resultados reportados**.

| Función | Qué hace |
|---|---|
| `_cost_matrix(theta)` | matriz 4×4 de sustitución de tipos (igual que en el núcleo) |
| `_edges(G)` | aristas como `{(u,v) ordenado: dist}` |
| `_incident(n, E)` | para cada nodo, la lista de distancias de sus aristas = su **huella local** |
| `_local_edge_cost(a, b, beta, indel)` | costo mínimo de alinear dos huellas locales; es un problema de asignación en miniatura, resuelto otra vez con el húngaro |
| `_total_cost(p, ...)` | costo total (nodos + aristas) de un emparejamiento `p` concreto |
| `_refine(p, ...)` | búsqueda local: intercambia pares de asignaciones mientras el costo baje |
| `ged_detallado(G1, G2, theta)` | devuelve `(similitud, mapping, ged, ged_normalizado)` |
| `ged_similarity_v2(...)` | reemplazo directo de `ged_similarity`; devuelve solo la similitud |
| `parchear()` | sustituye `ged_similarity` dentro del núcleo por la v2 |

### Cómo corrige el defecto, en tres pasos

1. **Matriz aumentada `Ca`.** Al costo de emparejar i con j se le suma
   `0.5 × _local_edge_cost(huella_i, huella_j)`. Así, entre dos nodos aromáticos
   empatados por tipo, el húngaro prefiere el que además está a distancias parecidas de
   los demás grupos. El factor ½ evita contar dos veces cada arista (una por cada
   extremo).
2. **Refinamiento** con `_refine`.
3. **Dos puntos de partida** (`Cn` y `Ca`), cada uno refinado, y se conserva el de menor
   costo real → **nunca sale peor que la versión anterior**.

Resultado: A9 vs A22 pasa de 0,43 a **1,000**.

### Cómo usarlo

```python
# opción A (recomendada en código nuevo)
from ged_mejorado import ged_similarity_v2
sim = ged_similarity_v2(G1, G2, theta)

# opción B (lo que hacen los scripts de este repo)
import ged_mejorado
ged_mejorado.parchear()                               # ANTES del import siguiente
from pipeline_final_completo import ged_similarity    # ya es la v2
```
El orden de la opción B **no es opcional**: `from … import ged_similarity` copia la
referencia a la función en ese instante, así que parchear después no tendría efecto.

### ⚠️ Limitación que hay que declarar en la tesis

`_refine` es una búsqueda local *greedy*: su resultado puede depender del **orden en que
estén numerados los nodos**. Dos grafos idénticos con numeración distinta —como A9 y
A22, cuyos átomos vienen en otro orden— pueden dar similitudes que difieren en el tercer
decimal frente a una misma candidata. Es un **óptimo local**, no el GED exacto (que es
NP‑difícil). Consecuencia medible: los top‑20 de A9 y de A22 coinciden en **17 de 20** y
no en 20 de 20, a pesar de que sus grafos ErG son el mismo grafo.

▶ Extrapolable: si se añaden tipos de nodo hay que ampliar `TYPES` y `_cost_matrix` aquí
**y** en el núcleo, con el mismo orden.

---

# 3. `ns1_utils.py` — criterios del bolsillo y métricas

No depende de RDKit a propósito: recibe el grafo ya construido como listas y arrays.

## Criterios geométricos

```python
PARAMS = {"tam_min": 15.0, "ratio_extremos": 0.70, "A_por_enlace": (1.5, 2.0)}
```
▶ **El bloque más específico de la diana.** Cambiar estos tres números es la forma más
rápida de adaptar el filtro a otra proteína; también se pueden pasar por línea de
comandos en `evaluar_erg_ns1.py` (`--tam-min`, `--ratio-extremos`, `--A-por-enlace`).

| Función | Qué hace |
|---|---|
| `es_conexo(n, aristas)` | union‑find con compresión de camino; criterio 2 |
| `largo_max(coords)` | distancia euclidiana máxima entre dos puntos = largo de la molécula en Å; criterio 4 |
| `criterios_ns1(nodos, aristas, coords_pesados, p)` | evalúa los cinco criterios y devuelve el detalle |

Claves importantes del diccionario que devuelve:
- `estructura` — `True` si cumple composición **y** conexo **y** extremos (el filtro duro).
- `tam_ok` — `True` si supera el largo mínimo. Va aparte porque un único conformero
  puede salir plegado y subestimar el largo real.
- `ajuste_dist` — fracción de aristas con Å/enlace en rango. Variable continua, más útil
  como descriptor que como filtro.

▶ Si el bolsillo de la diana nueva tiene otra forma (una cavidad compacta en vez de un
canal alargado), hay que **reescribir `criterios_ns1`** con las condiciones que
correspondan. El resto del pipeline no cambia.

## Métricas

| Función | Qué calcula |
|---|---|
| `auc(y, s)` | AUC‑ROC como estadístico de Mann‑Whitney: fracción de pares (activo, inactivo) en que el activo puntúa más alto; los empates cuentan medio. No necesita ordenar umbrales. |
| `ef(y, s, frac)` | factor de enriquecimiento = (proporción de activos en el top k) / (proporción global). EF = 1 es azar; EF = 5 es cinco veces mejor. **Es la métrica que de verdad le importa a un cribado.** |
| `auc_ci(y, s, n=300)` | intervalo de confianza del 95 % por bootstrap estratificado. Con 5–30 activos, reportar un AUC sin intervalo es engañoso. |
| `rango_norm(s)` | convierte puntajes a rangos en [0,1], para fusionar métodos de escalas distintas (ErG+GED vive en 0,8–1,0; Tanimoto en 0,1–0,5: sumar los valores crudos no tendría sentido, sumar los rangos sí). |

---

# 4. `pipeline_v3.py` — variante endurecida

Mismo método, cuatro correcciones de robustez. **Es la versión recomendada para una
corrida desde cero.**

1. Importa directamente `ged_similarity_v2`: no hay que acordarse de llamar a
   `parchear()`. `ged_similarity_v1` queda disponible para reproducir lo antiguo.
2. `_verificar_no_existe()`: se niega a sobrescribir un resultado ya guardado. Evita
   perder una corrida de 45 minutos por un descuido.
3. Registra en un CSV aparte las moléculas cuyo embedding 3D falla, en vez de
   descartarlas en silencio (de ahí salen las 17 de 33 047).
4. `elegir_chemotipos_solo_train()`: corrige una **fuga de datos**. Antes, los tres
   quimiotipos más grandes se elegían mirando **todos** los activos, incluidos los que
   después se usaban como prueba. Ahora se clusteriza **solo** con el conjunto de
   entrenamiento y una molécula de prueba entra únicamente si su Tanimoto contra algún
   miembro de esos quimiotipos supera el mismo umbral de Butina. También arregla la
   lectura de `AID2326_cache.csv` (filtrado por contenido en vez de saltar cuatro filas
   por posición) y la ruta de `Moleculas.xlsx`.

Los demás scripts importan `pipeline_final_completo` porque así se corrieron los
resultados reportados. Para trabajo nuevo: cambiar ese import por `pipeline_v3` y quitar
la llamada a `parchear()`.

---

# 5. `recalibrar_v2.py` — calibración y validación de θ

Produce `theta_final_v2.json` (los pesos finales) y
`recalibracion_v2_resultados.json` (los cinco experimentos).

### Qué hace, en orden

1. Carga AID 2326 desde `BasePubChem.xlsx` y deja `AID2326_cache.csv` para las
   siguientes corridas (la primera lectura de 270 000 filas tarda varios minutos).
2. Limpia sales con `_limpiar_smiles_sal`.
3. Agrupa los activos con Butina.
4. Corre **cinco experimentos**, todos con partición entrenamiento/prueba:

| Experimento | Qué mide |
|---|---|
| `mezclado_azar` | 35 activos al azar de los 1 073: qué pasa al mezclar mecanismos de acción |
| `cluster_1`, `cluster_2`, `cluster_3` | cada quimiotipo grande por separado: si discrimina **dentro** de una familia |
| `3_clusters_juntos` | los tres juntos → **produce el θ final** |

5. Calcula, sobre **exactamente las mismas moléculas de prueba**, la línea base
   Morgan+Tanimoto. Sin esa comparación no se puede afirmar que el método aporte algo.
6. Imprime una tabla comparando contra los AUC obtenidos con el GED defectuoso (el
   diccionario `anterior` está ahí solo para eso).

### `calibrar_con_progreso(...)`

Es `calibrar_y_evaluar_train_test` con dos diferencias: imprime el avance de la
optimización con estimación de tiempo restante, y usa `polish=False` (el pulido final de
`differential_evolution` no aporta nada cuando la función objetivo es escalonada y gasta
evaluaciones).

### Cómo se corre y qué cambiar

```bash
python recalibrar_v2.py --rapido    # prueba rápida, pocas iteraciones
python recalibrar_v2.py             # corrida completa (~1–2 h en total)
```

▶ **Extrapolable:**
- el bloque de datos: reemplazar `BasePubChem.xlsx` por el cribado etiquetado de la diana
  nueva (necesita columnas CID / SMILES / Active‑Inactive);
- `N_INACT = 60` — inactivos por experimento. Subirlo estabiliza las estimaciones y
  multiplica el tiempo de cómputo;
- `MAXITER`, `POPSIZE` — esfuerzo de la evolución diferencial;
- `cutoff_distancia = 0.6` — definición de quimiotipo;
- `SEED = 42` — fijo, para que la partición sea reproducible.

---

# 6. `correr_ranking.py` — aplicar el método a la biblioteca

Recorre `Moleculas.xlsx` (33 047 compuestos de ChEMBL) en bloques de 1 500, calcula para
cada molécula su ErG+GED y su Tanimoto frente a **cada** referencia por separado, y al
final imprime y guarda el top‑20 de cada una.

### Las tres configuraciones

```bash
python correr_ranking.py                # GED v1 + theta_final.json      (resultados viejos)
python correr_ranking.py --v2           # GED corregido + theta_final.json
python correr_ranking.py --theta-nuevo  # GED corregido + theta_final_v2.json  ← la final
```
El sufijo de los archivos de salida cambia con la opción (`""`, `"_v2"`, `"_v2b"`), así
que las tres corridas conviven sin pisarse.

> Nota de trazabilidad: el ranking final del repositorio se llama
> `ranking_33k_thetav2.csv` y trae un CSV de errores, porque se generó con la función de
> ranking de `pipeline_v3.py` usando `theta_final_v2.json`. Es la misma configuración
> que `--theta-nuevo`, solo con otro nombre de archivo.

### Puede interrumpirse

Si se corta, se vuelve a lanzar y retoma desde `ranking_33k_checkpoint*.json`. Para
empezar de cero hay que borrar el checkpoint **y** el CSV de salida.

### El cronómetro

```python
pf.smiles_to_3d_mol = cronometrar("embedding_3D", pf.smiles_to_3d_mol)
pf.build_erg        = cronometrar("construir_ErG", pf.build_erg)
pf.ged_similarity   = cronometrar("GED", pf.ged_similarity)
```
Envuelve las funciones del pipeline **sin modificarlo** y escribe
`tiempos_por_bloque*.csv`. De ahí salen los números de eficiencia computacional que pide
el objetivo específico 4.

▶ **Extrapolable:** `EXCEL` (biblioteca a cribar), las referencias importadas
(`INFLUENZA_ACTIVOS` → `DENGUE_ACTIVOS`), `THETA_FILE` (pesos de esa diana) y `BATCH`
(solo afecta memoria y frecuencia de guardado, no los resultados).

---

# 7. `evaluar_erg_ns1.py` — la evaluación justa

Es el experimento más exigente del proyecto. Produce
`evaluacion_erg_ns1_resultados.csv`.

### Tres escenarios

| Escenario | Qué mide |
|---|---|
| `mezclado` | entrenar y probar con activos al azar (mecanismos distintos) |
| `intra-cluster k` | entrenar y probar dentro del mismo quimiotipo |
| `entre-clusters: train k -> test otros` | entrenar con una familia, probar con **otras**. **La prueba realmente justa**: por construcción, aquí Tanimoto no tiene la ventaja de estar buscando análogos cercanos |

### Dos tipos de inactivos

- **fáciles**: al azar.
- **difíciles**: los de propiedades fisicoquímicas más parecidas a los activos (peso,
  logP, donadores, aceptores, enlaces rotables, en espacio z‑score), al estilo DUD‑E. El
  criterio es **independiente** tanto de Tanimoto como de ErG, para no favorecer a
  ninguno de los dos.

### Ocho métodos

`Tanimoto (max)`, `Tanimoto (media)`, `ErG-fp (max)`, `ErG+GED (max)`,
`ErG+GED (media)`, `Fusion ErG+GED + Tanimoto`, `ErG+GED (max) + filtro NS1`,
`Solo filtro NS1`.

- **"max"** = similitud al activo de entrenamiento **más** parecido (vecino más
  cercano). Es lo correcto cuando los activos no comparten quimiotipo.
- **"media"** = promedio sobre todos los activos de entrenamiento.

### Métricas

AUC con intervalo de confianza del 95 % por bootstrap, y EF10 / EF20.

### Funciones propias

| Función | Qué hace |
|---|---|
| `analizar(smi)` | grafo ErG + Morgan‑fp + ErG‑fp + criterios NS1, todo de una vez |
| `item(cid, smi)` | caché en memoria: cada molécula se analiza una sola vez por corrida |
| `largo_multiconf(smi, n=30)` | genera 30 conformeros y toma el largo máximo: evita subestimar el tamaño por un conformero plegado |
| `propiedades(smi)` | el vector de 5 descriptores para elegir los inactivos difíciles |
| `calibrar(train, decoys, seed)` | calibra θ **solo** con el entrenamiento de cada prueba (sin fuga) |
| `puntajes(test, train, theta)` | calcula los ocho métodos de una vez |
| `cribado()` | aplica criterios NS1 + ErG+GED a las 33 047 y produce el embudo de criterios |

Detalle importante de `calibrar`:
```python
def a_theta(v):
    return {..., "w_polar_apolar": v[0] + v[2], ...}   # v[2] >= 0
```
Esta parametrización **impone** `w_polar_apolar ≥ w_aromatico_hidrofobo` por
construcción: confundir polar con apolar nunca puede costar menos que confundir
aromático con hidrofóbico. Es una restricción química razonable y reduce el espacio de
búsqueda. ▶ Si en otra diana no tiene sentido, se quita de `a_theta`.

```bash
python evaluar_erg_ns1.py --solo-referencias   # ¿cumplen A9/A22 los criterios?
python evaluar_erg_ns1.py                      # evaluación (minutos)
python evaluar_erg_ns1.py --calibrar           # con θ calibrado por prueba (lento)
python evaluar_erg_ns1.py --cribado            # cribado NS1 de las 33k (~30 min)
```

---

# 8. `comparar_A9_A22.py` — abrir la caja negra

En vez de un número de similitud, imprime **nodo por nodo** cómo se emparejaron los dos
grafos, qué grupos quedaron sin pareja y en cuántos enlaces difiere cada distancia. Es
el script para discutir un resultado con un químico.

- **Parte 1:** A9 vs A22 en detalle — similitud con θ calibrado y por defecto, Tanimoto,
  emparejamiento de nodos (`detalle_ged`, que replica el GED imprimiendo el paso a
  paso), diferencias de distancia, subestructura común máxima (`rdFMCS`) y propiedades
  fisicoquímicas.
- **Parte 2:** top‑20 de A9 frente a top‑20 de A22 — matrices cruzadas 20×20 con ErG+GED
  y con Tanimoto, similitud interna de cada lista, similitud de cada lista contra ambos
  activos, correlación de Spearman entre las dos medidas y mapas de calor.

```bash
python comparar_A9_A22.py        # GED v1
python comparar_A9_A22.py --v2   # GED corregido  ← usar este
```

`nodo_str(P, i)` es la función que materializa el "no se pierde la información atómica":
usa `atom_ids` para imprimir los símbolos químicos de cada nodo, de modo que se lee
`Aromatico[CCCCCC]` o `HBD[O]`.

---

# 9. `dibujar_top20.py` — la figura

Rejilla con las referencias y las 20 moléculas, cada una rotulada con su ErG+GED y su
Tanimoto frente a A9 y A22 (`Draw.MolsToGridImage`). Además imprime dos diagnósticos que
importan para la tesis:

1. cuántas moléculas de las 33 047 **empatan** con el puntaje del puesto 20 (si hay
   empates, el corte exacto del top‑20 es en parte arbitrario);
2. cuántos valores de similitud **distintos** existen entre las 33 047, que es una
   medida de la granularidad del método.

```bash
python dibujar_top20.py        # archivos _v2
python dibujar_top20.py _v2b   # GED corregido + θ recalibrado
python dibujar_top20.py ""     # archivos originales
```

---

# 10. Scripts históricos

`calibracion_cluster.py` y `calibracion_pubchem_real.py` son las primeras versiones de
la calibración, con el GED sin corregir. **No son reproducibles tal cual**: dependen de
archivos intermedios que ya no están (`activos_clusters.pkl`,
`AID2326_activos_reales.csv`, `AID2326_inactivos_muestra5000.csv`).

Se conservan por trazabilidad: ahí aparecieron por primera vez los AUC por cluster
(0,962 / 0,671 / 0,733) que cita el anteproyecto, y
`calibracion_pubchem_real.py` es el que produjo `figuras/pubchem_roc_test.png`, la curva
ROC con AUC_test = 0,49 que motivó todo el análisis de clustering.

Para calibrar hoy: `recalibrar_v2.py`.

---

# 11. Orden de ejecución recomendado

```bash
# 0) entorno
pip install -r requirements.txt

# 1) calibrar y validar los pesos  (necesita BasePubChem.xlsx; ~1-2 h)
python recalibrar_v2.py
#    -> theta_final_v2.json, recalibracion_v2_resultados.json, AID2326_cache.csv

# 2) evaluación justa frente a Tanimoto  (minutos)
python evaluar_erg_ns1.py
#    -> evaluacion_erg_ns1_resultados.csv

# 3) aplicar a la biblioteca completa  (necesita Moleculas.xlsx; ~45 min)
python correr_ranking.py --theta-nuevo
#    -> ranking_33k_*.csv, top20_A9_*.csv, top20_A22_*.csv, tiempos_por_bloque_*.csv

# 4) análisis e interpretación
python comparar_A9_A22.py --v2
python dibujar_top20.py _v2b
```

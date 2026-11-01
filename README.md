# ErG + GED + algoritmo húngaro para priorizar inhibidores de NS1

Prefiltro computacional que ordena una biblioteca de **33 047 moléculas** según su
parecido farmacofórico con inhibidores conocidos de la proteína **NS1**, usando
**Extended Reduced Graphs** (ErG), **distancia de edición de grafos** (GED) y el
**algoritmo húngaro** de asignación óptima.

Código del trabajo de grado de **Salomé Castro Gómez** — Programa de Ciencia de Datos,
Departamento de Matemáticas, Universidad Externado de Colombia. Tutor: PhD. Germán
Andrés Combariza.

---

## El problema

El dengue sigue siendo endémico en decenas de países y **no existe un tratamiento
antiviral específico aprobado**. La proteína no estructural 1 (NS1) es un blanco
terapéutico atractivo por su papel en la replicación viral y la evasión inmunitaria,
pero el número de inhibidores de NS1 de dengue confirmados experimentalmente es
mínimo — un puñado de compuestos.

Esa escasez de datos es, a la vez, el problema y el obstáculo para resolverlo: sin
inhibidores conocidos no hay con qué entrenar ni validar un método computacional.

**Estrategia de dos fases:**

- **Lo que hace este repositorio.** Desarrollar, poner a prueba y aplicar el método sobre la **NS1 de influenza**, una
  proteína homóloga para la que sí existen inhibidores caracterizados (A9 y A22,
  Kleinpeter et al. 2018) y un cribado experimental público de 270 315 compuestos con
  1 073 activos confirmados (PubChem AID 2326). Ahí se puede medir honestamente si el
  método discrimina.
- **Extrapolación pendiente.** Aplicar el método ya calibrado a NS1 de dengue, con los pocos
  inhibidores confirmados que existen.

El método es **extrapolable**: cambiar de diana requiere cambiar las moléculas de
referencia, recalibrar seis pesos y ajustar los criterios geométricos del bolsillo. Todo
lo demás es genérico. Ver [`docs/03_como_extrapolar.md`](docs/03_como_extrapolar.md).

---

## Cómo funciona, en cinco líneas

1. Cada molécula (SMILES) se convierte en un conformero 3D (ETKDG + UFF).
2. Sus átomos se tipifican en cuatro familias farmacofóricas (donador, aceptor,
   aromático, hidrofóbico), **una sola por átomo**, y los contiguos del mismo tipo se
   **colapsan en un nodo**: un benceno pasa de 6 átomos a 1 nodo.
3. Los nodos se conectan con aristas etiquetadas por la **distancia topológica** (número
   de enlaces, de 1 a 6). El resultado es un grafo de ~6 nodos en vez de ~20 átomos.
4. Dos grafos se comparan con la **distancia de edición**: el emparejamiento óptimo de
   nodos se resuelve con el **algoritmo húngaro** en O(N³), convirtiendo un problema
   NP‑difícil en uno de microsegundos. El costo se normaliza y se convierte en similitud
   con `exp(−GED/(n₁+n₂+m₁+m₂))`.
5. Los seis pesos del costo de edición (**θ**) se **calibran con datos** por evolución
   diferencial, maximizando el AUC sobre un conjunto etiquetado, con partición
   entrenamiento/prueba.

Explicación completa, sin suponer conocimientos de quimioinformática:
[`docs/01_explicacion_del_metodo.md`](docs/01_explicacion_del_metodo.md).

---

## Resultados principales

Todos los números están verificables en [`resultados/`](resultados) y detallados en
[`docs/04_resultados.md`](docs/04_resultados.md).

| | |
|---|---|
| Moléculas cribadas | **33 013** de 33 047 (17 descartadas por fallo de embedding 3D) |
| Tiempo total del cribado | **≈ 24 min** en un portátil, sin paralelizar |
| Costo por comparación GED | **0,86 ms** |
| AUC de prueba, dentro de quimiotipos coherentes | **0,88 – 0,98** |
| AUC de prueba, con activos de mecanismos mezclados | **0,48** (azar) |
| Similitud A9 vs A22 tras corregir el GED | **1,000** (antes: 0,43) |

### Los tres hallazgos del trabajo

**1. Un defecto de implementación que falseaba todo.** La primera versión dejaba que el
algoritmo húngaro emparejara nodos mirando solo su tipo; con nodos del mismo tipo la
matriz de costos empata y el emparejamiento salía cruzado. Síntoma: A9 y A22, cuyos
grafos ErG son **idénticos**, daban similitud 0,43. Corregido con el esquema estándar de
GED bipartito (costo aumentado con la estructura local de aristas + refinamiento local),
la similitud pasa a 1,000 y el AUC de prueba sube de 0,65 a **0,865**.

**2. El método sirve para buscar análogos, no para predecir actividad.** Calibrado solo
contra A9 y A22 (n = 2), el AUC parecía perfecto: 1,000. Al recalibrar contra los 1 073
activos reales de PubChem con partición entrenamiento/prueba, el AUC de prueba cayó a
**0,479**: azar. La explicación apareció al agrupar esos activos con el algoritmo de
Butina: **737 quimiotipos, 588 de ellos de una sola molécula**. Un cribado fenotípico
mezcla mecanismos de acción sin farmacóforo común. Al validar **dentro** de los tres
quimiotipos más grandes, el AUC subió a **0,98 / 0,91 / 0,88**.
→ Regla operativa: **reportar siempre por activo de referencia, nunca promediado.**

**3. El método recupera el quimiotipo correcto, y con moléculas que un fingerprint
dejaría abajo.** De las 20 moléculas mejor rankeadas frente a A9 (una arilamida de
fenetilamina), **17 son arilamidas de fenetilamina o variantes cercanas**, con el anillo
sustituido o reemplazado por furano, benzofurano, tiofeno, indol o piridina. Sus
coeficientes de Tanimoto van de **0,125 a 0,525**: en términos de fingerprint no son
"muy parecidas", pero sí comparten la organización farmacofórica. La correlación de
Spearman entre ErG+GED y Tanimoto es **0,166** — las dos medidas miden cosas distintas.

### Y un resultado que conviene leer con cuidado

En la evaluación más estricta (θ **sin** calibrar, inactivos difíciles por propiedades
fisicoquímicas al estilo DUD‑E, siete escenarios), **Morgan+Tanimoto es igual o mejor que
ErG+GED en 5 de los 7 escenarios**, y ningún método generaliza entre quimiotipos
distintos. Lo que sí resulta consistentemente bueno es la **fusión por rangos de ambos**.
El aporte defendible del método es la complementariedad, no la superioridad. Está
reportado tal cual en [`docs/04_resultados.md`](docs/04_resultados.md), con intervalos de
confianza.

---

## Estructura del repositorio

```
├── pipeline_final_completo.py   NÚCLEO: ErG, GED, referencias, calibración, ranking
├── ged_mejorado.py              GED corregido (el que produce los resultados)
├── ns1_utils.py                 criterios geométricos del bolsillo + métricas
├── pipeline_v3.py               misma funcionalidad, endurecida y sin fuga de datos
│
├── recalibrar_v2.py             calibra y valida θ        → theta_final_v2.json
├── evaluar_erg_ns1.py           evaluación justa vs Tanimoto
├── correr_ranking.py            aplica el método a las 33 047 candidatas
├── comparar_A9_A22.py           análisis explicativo de un par de moléculas
├── dibujar_top20.py             figura con las estructuras del top-20
├── calibracion_*.py             scripts históricos (trazabilidad, no reproducibles)
│
├── docs/
│   ├── 01_explicacion_del_metodo.md    conceptos, desde cero
│   ├── 02_explicacion_del_codigo.md    archivo por archivo, función por función
│   ├── 03_como_extrapolar.md           receta para otra proteína
│   └── 04_resultados.md                todos los números, con su archivo de origen
│
├── resultados/                  pesos calibrados, rankings, métricas (archivos chicos)
├── figuras/                     curva ROC, mapas de calor, explorador ErG interactivo
└── data/                         (vacía) cómo obtener los datos — ver data/README.md
```

Los datos de entrada **no están en el repositorio** por tamaño (más de 60 MB).
[`data/README.md`](data/README.md) explica cómo obtener cada archivo.

---

## Instalación

```bash
git clone <url-del-repositorio>
cd erg-ns1
python -m venv .venv

# Windows
.venv\Scripts\activate
# Linux / macOS
source .venv/bin/activate

pip install -r requirements.txt
```

Python 3.10–3.13. La dependencia principal es **RDKit** (código abierto), que trae el
catálogo de definiciones farmacofóricas (`BaseFeatures.fdef`) sobre el que se construyen
los grafos.

### Con conda

```bash
conda create -n erg-ns1 -c conda-forge python=3.13 rdkit
conda activate erg-ns1
pip install scipy scikit-learn pandas networkx matplotlib openpyxl
```

### Entorno verificado

Los resultados de [`docs/04_resultados.md`](docs/04_resultados.md) se reprodujeron con
Python 3.13.2, rdkit 2025.03.6, numpy 2.2.2, scipy 1.18.1, scikit-learn 1.9.1,
pandas 2.2.3, networkx 3.5, matplotlib 3.10.0 y openpyxl 3.1.5.

Comprobación rápida de que el entorno funciona: los dos inhibidores de referencia deben
dar similaridad exactamente 1,000, porque sus grafos ErG son idénticos.

```python
import ged_mejorado, json
ged_mejorado.parchear()
from pipeline_final_completo import INFLUENZA_ACTIVOS, smiles_to_3d_mol, build_erg

theta = json.load(open("resultados/theta_final_v2.json"))["theta"]
G = {k: build_erg(smiles_to_3d_mol(s)[0]) for k, s in INFLUENZA_ACTIVOS.items()}
print(ged_mejorado.ged_similarity_v2(G["A9"], G["A22"], theta))   # -> 1.0
```

---

## Uso

```bash
# 1) Calibrar y validar los pesos θ   (necesita BasePubChem.xlsx; ~1-2 h)
python recalibrar_v2.py --rapido     # prueba rápida primero
python recalibrar_v2.py              # corrida completa
#    → theta_final_v2.json, recalibracion_v2_resultados.json, AID2326_cache.csv

# 2) Evaluación justa frente a Tanimoto   (minutos)
python evaluar_erg_ns1.py
#    → evaluacion_erg_ns1_resultados.csv

# 3) Aplicar a la biblioteca completa   (necesita Moleculas.xlsx; ~25-45 min)
python correr_ranking.py --theta-nuevo
#    → ranking_33k_*.csv, top20_A9_*.csv, top20_A22_*.csv, tiempos_por_bloque_*.csv

# 4) Análisis e interpretación
python comparar_A9_A22.py --v2      # emparejamiento nodo por nodo, MCS, propiedades
python dibujar_top20.py _v2b        # figura con las estructuras
```

El paso 3 **puede interrumpirse**: guarda un *checkpoint* en disco y al volver a
lanzarlo retoma donde quedó.

---

## Qué cambiar para aplicarlo a otra proteína

Resumen; la receta completa, con ejemplos, está en
[`docs/03_como_extrapolar.md`](docs/03_como_extrapolar.md).

| # | Qué | Dónde | ¿Obligatorio? |
|---|---|---|---|
| 1 | Inhibidores de referencia | `pipeline_final_completo.py` → `INFLUENZA_ACTIVOS` | **sí** |
| 2 | Pesos θ (recalibrar, no editar) | `recalibrar_v2.py` → nuevo `theta_*.json` | **sí**, si hay datos etiquetados |
| 3 | Criterios geométricos del bolsillo | `ns1_utils.py` → `PARAMS`, `criterios_ns1` | recomendable |
| 4 | Tipos de nodo (p. ej. cargas) | `LABEL_MAP`, `PRIORITY`, `TYPES`, `_cost_matrix` | solo si hace falta |
| 5 | Rango de conexión / umbral de quimiotipo | `build_erg` (`1<=d<=6`), `clusterizar_activos` | casos especiales |

**Dos reglas que el trabajo sobre influenza dejó claras y que hay que respetar al extrapolar:**

- Un ranking **por cada** inhibidor de referencia, nunca promediado entre varios.
- Declarar siempre de dónde salió θ: calibrado con datos de esa diana, transferido de
  otra o sin calibrar.

---

## Fuentes de los datos y de las referencias

| Recurso | Fuente |
|---|---|
| A9, A22 (inhibidores de NS1 de influenza) | Kleinpeter et al. 2018, *J. Biol. Chem.* — bolsillo CPSF30 |
| Compuestos con IC50 de influenza | Basu et al. 2009, *J. Virol.*, DOI 10.1128/JVI.01805-08 |
| Conjunto etiquetado de validación | [PubChem BioAssay AID 2326](https://pubchem.ncbi.nlm.nih.gov/bioassay/2326) — NIH Chemical Genomics Center |
| Biblioteca de candidatos | ChEMBL (descarga de 33 047 compuestos) |
| Inhibidores de NS1 de dengue (extrapolación prevista) | Ansari et al. 2025, *Molecular Diversity*, DOI 10.1007/s11030-025-11447-5 |

## Métodos y algoritmos

| Algoritmo | Referencia |
|---|---|
| Extended Reduced Graphs | Stiefl, Watson, Baumann & Zaliani (2006), *J. Chem. Inf. Model.* |
| GED aproximado por asignación bipartita | Riesen & Bunke (2009), *Image and Vision Computing* |
| Algoritmo húngaro de asignación óptima | Kuhn (1955); Munkres — vía `scipy.optimize.linear_sum_assignment` |
| ErG + GED para cribado virtual por ligando | García‑Hernández, Fernández & Serratosa (2019), *J. Chem. Inf. Model.* |
| Agrupamiento de quimiotipos | Butina (1999), *J. Chem. Inf. Comput. Sci.* |
| Evolución diferencial | Storn & Price (1997), *Journal of Global Optimization* |
| Quimioinformática | RDKit — Landrum et al. |

---

## Limitaciones

Verificadas con los datos de este repositorio; el detalle está en
[`docs/04_resultados.md`](docs/04_resultados.md) §7.

- Es una **búsqueda de análogos**, no un predictor de actividad biológica.
- **No generaliza** a familias químicas no vistas (AUC ≤ 0,78 entre quimiotipos).
- El GED es **aproximado** y su refinamiento es *greedy*: depende del orden de numeración
  de los nodos.
- La reducción ErG **pierde resolución por diseño**: 22 064 de 33 013 moléculas reciben el
  mismo puntaje frente a las dos referencias.
- **Un solo conformero** por molécula.
- El filtro geométrico de NS1 **no mejoró** el ranking (AUC 0,30–0,63).
- **Ninguna molécula priorizada se ha ensayado experimentalmente.** El resultado es una
  hipótesis priorizada, pensada como prefiltro previo al docking o al ensayo, no un
  hallazgo biológico.

## Licencia

MIT — ver [`LICENSE`](LICENSE).

# 04 — Resultados

Todos los números de este documento salen de los archivos de la carpeta
[`../resultados/`](../resultados). Cada tabla indica su archivo de origen, de modo que
cualquier cifra se puede verificar.

Salvo cuando se diga lo contrario, los resultados corresponden a la configuración final:
**GED corregido** (`ged_mejorado.py`) con **θ recalibrado** (`theta_final_v2.json`).

---

## 1. Los datos

### Biblioteca de candidatos — `Moleculas.xlsx` (descarga de ChEMBL)

| | |
|---|---|
| Registros | 33 047 |
| Con SMILES válido | 33 030 |
| Tipo "Small molecule" | 33 034 |
| Con fase clínica declarada (`Max Phase`) | 180 (127 de ellas aprobadas, fase 4) |
| Descartadas por fallo del embedding 3D | **17** |
| **Procesadas en el ranking final** | **33 013** |

Descriptores (del propio archivo de ChEMBL):

| Descriptor | n | media | mín | mediana | máx |
|---|---|---|---|---|---|
| Peso molecular (Da) | 33 037 | 358,0 | 124,2 | 350,4 | 1 485,7 |
| AlogP | 32 998 | 3,41 | −3,48 | 3,44 | 13,59 |
| Átomos pesados | 32 998 | 24,9 | 6 | 24 | 73 |
| Anillos aromáticos | 32 998 | 2,51 | 0 | 2 | 8 |
| Enlaces rotables | 32 998 | 4,59 | 0 | 4 | 30 |
| Aceptores (HBA) | 32 998 | 4,84 | 0 | 5 | 22 |
| Donadores (HBD) | 32 998 | 1,04 | 0 | 1 | 13 |

### Conjunto etiquetado de validación — PubChem AID 2326

Cribado de alto rendimiento (qHTS) del NIH Chemical Genomics Center contra la función
de NS1 de influenza; es el mismo cribado del que salieron A9 y A22.

| Desenlace | Compuestos |
|---|---|
| Active | **1 073** |
| Inactive | 260 684 |
| Inconclusive | 8 410 |
| Unspecified | 130 |
| **Total** | **270 315** |

Solo se usan `Active` e `Inactive`.

---

## 2. Corrección del GED: el resultado que lo hizo todo medible

`resultados/comparacion_A9_A22_resumen_v2.txt`

A9 (2‑hidroxi‑N‑fenetilbenzamida) y A22 (2‑amino‑N‑fenetilbenzamida) difieren en **un
átomo**: –OH frente a –NH₂. Los dos se unen al bolsillo CPSF30 de NS1 de influenza.

| | Antes (GED v1) | Después (GED corregido) |
|---|---|---|
| Similitud ErG+GED A9–A22 | **0,43** | **1,000** |
| Nodos / aristas de cada grafo | 6 / 14 | 6 / 14 |
| Nodos emparejados a coste 0 | — | 6 de 6 |
| Pares con distancia distinta | — | **0** |
| Morgan + Tanimoto | 0,706 | 0,706 |
| Subestructura común máxima | 17 átomos de 18 | 17 de 18 |

Lectura: los dos grafos ErG son **idénticos**, porque la reducción trata –OH y –NH₂ como
el mismo tipo de nodo (donador de puente de hidrógeno). La representación acertó: las dos
moléculas se unen al mismo sitio. Un fingerprint de Morgan sí los distingue (0,706),
porque distingue algo que a la proteína le da igual.

**Consecuencia metodológica:** en Fase A las dos referencias aportan al método la misma
información. No son dos validaciones independientes.

Propiedades fisicoquímicas de las dos referencias:

| | A9 | A22 |
|---|---|---|
| Peso molecular | 241,29 | 240,31 |
| logP | 2,36 | 2,24 |
| HBD / HBA | 2 / 2 | 2 / 2 |
| TPSA | 49,33 | 55,12 |
| Enlaces rotables | 4 | 4 |

---

## 3. Calibración y validación de θ

`resultados/recalibracion_v2_resultados.json`, `resultados/theta_final_v2.json`

Cinco experimentos, todos con partición entrenamiento/prueba; el AUC reportado es
**siempre** sobre moléculas que no se usaron para calibrar. La última columna es la línea
base Morgan+Tanimoto **sobre exactamente las mismas moléculas de prueba**.

| Experimento | Activos entren. | Activos prueba | AUC entren. | **AUC prueba** | AUC prueba Tanimoto | AUC prueba antes de corregir el GED | min |
|---|---|---|---|---|---|---|---|
| Activos mezclados al azar | 21 | 14 | 0,663 | **0,479** | 0,563 | 0,49 | 19,4 |
| Quimiotipo 1 | 12 | 8 | 1,000 | **0,984** | 1,000 | 0,962 | 11,1 |
| Quimiotipo 2 | 9 | 7 | 1,000 | **0,911** | 1,000 | 0,671 | 3,4 |
| Quimiotipo 3 | 7 | 5 | 0,972 | **0,875** | 1,000 | 0,733 | 43,3 |
| Los 3 quimiotipos juntos | 28 | 20 | 0,902 | **0,865** | 1,000 | 0,65 | 35,5 |

Cada experimento usa 36 inactivos de entrenamiento y 24 de prueba, tomados de AID 2326.

### Los pesos finales

`theta_final_v2.json`, obtenidos del experimento "3 quimiotipos juntos":

| Peso | Antes de corregir el GED | **Final** |
|---|---|---|
| `w_aromatico_hidrofobo` | 0,244 | **1,893** |
| `w_hba_hbd` | 1,668 | **0,801** |
| `w_polar_apolar` | 1,754 | **0,421** |
| `node_indel` | 2,318 | **0,670** |
| `edge_beta` | 1,166 | **1,469** |
| `edge_indel` | 0,714 | **0,465** |
| AUC entrenamiento | 0,950 | 0,902 |
| **AUC prueba** | 0,650 | **0,865** |

Dos observaciones que conviene declarar:

1. La corrección del GED **sube el AUC de prueba de 0,65 a 0,865** y, a la vez, baja el de
   entrenamiento de 0,95 a 0,90: menos memorización y más generalización. Es el efecto
   que se espera de corregir un defecto real.
2. Los pesos calibrados **no coinciden con la intuición química**. Se esperaría que
   `w_polar_apolar` (confundir polar con apolar) fuera el costo más alto, y resulta el
   más bajo (0,421), mientras que `w_aromatico_hidrofobo` es el más alto (1,893). Esto
   quiere decir que, para discriminar en este conjunto, lo que más informa es **distinguir
   anillos aromáticos de cadenas apolares**, no distinguir polar de apolar. Con 28
   activos de entrenamiento es un resultado que hay que tomar con cautela: puede ser una
   propiedad del quimiotipo dominante y no un principio general. Por eso
   `evaluar_erg_ns1.py` ofrece una calibración con la restricción
   `w_polar_apolar ≥ w_aromatico_hidrofobo` impuesta por construcción.

### El hallazgo metodológico central

`resultados/recalibracion_v2_resultados.json` + `clusterizar_activos`

| Paso | AUC de prueba | Qué significa |
|---|---|---|
| Calibrar solo con A9 y A22 (n = 2) | 1,000 (de calibración) | **espejismo**: con dos activos casi idénticos, cualquier método memoriza |
| Recalibrar con los 1 073 activos mezclados | **0,479** | azar |
| Validar **dentro** de quimiotipos coherentes | **0,875 – 0,984** | el método sí funciona, pero solo dentro de una familia química |

La explicación está en el agrupamiento de Butina de los 1 073 activos:

| | |
|---|---|
| Clusters (quimiotipos) | **737** |
| De una sola molécula (*singletons*) | **588** |

Es decir: el cribado fenotípico primario mezcla compuestos que inhiben NS1 por
mecanismos completamente distintos, sin un farmacóforo común que el método pueda
aprender a reconocer de forma agregada.

**Conclusión:** ErG+GED+húngaro sirve para priorizar **análogos estructurales de un
inhibidor conocido** (búsqueda de "primos químicos"), no para detectar actividad
biológica en general. De ahí la regla operativa de todo el proyecto: **reportar siempre
por activo de referencia, nunca promediado**.

---

## 4. Evaluación justa frente a Tanimoto

`resultados/evaluacion_erg_ns1_resultados.csv` — 144 filas: 7 escenarios × 2 tipos de
inactivos × 8 métodos, con intervalos de confianza del 95 % por bootstrap.

**Importante:** esta evaluación se corrió con **θ sin calibrar** (`DEFAULT_THETA`), que es
la condición sin fuga de información más estricta posible. No usa los pesos
optimizados. Es el escenario menos favorable para el método propuesto, y así hay que
leerlo.

### AUC con inactivos difíciles (los de propiedades fisicoquímicas más parecidas)

| Escenario | Tanimoto (max) | ErG‑fp (max) | ErG+GED (max) | ErG+GED (media) | Fusión | Solo filtro NS1 |
|---|---|---|---|---|---|---|
| Mezclado | 0,477 | 0,544 | 0,487 | 0,487 | 0,471 | 0,521 |
| Intra‑quimiotipo 1 | **0,975** | 0,834 | 0,928 | 0,931 | 0,956 | 0,563 |
| Intra‑quimiotipo 2 | **1,000** | 0,989 | 0,875 | 0,879 | 0,950 | 0,300 |
| Intra‑quimiotipo 3 | **0,975** | 0,725 | 0,760 | 0,915 | 0,933 | 0,513 |
| Entre quimiotipos (1→otros) | 0,504 | 0,240 | 0,388 | 0,355 | 0,443 | 0,364 |
| Entre quimiotipos (2→otros) | 0,570 | **0,727** | 0,643 | 0,473 | 0,663 | 0,500 |
| Entre quimiotipos (3→otros) | **0,608** | 0,504 | 0,482 | 0,553 | 0,567 | 0,460 |

### AUC con inactivos fáciles (al azar)

| Escenario | Tanimoto (max) | ErG‑fp (max) | ErG+GED (max) | ErG+GED (media) | Fusión | Solo filtro NS1 |
|---|---|---|---|---|---|---|
| Mezclado | 0,527 | 0,608 | 0,486 | 0,548 | 0,490 | 0,558 |
| Intra‑quimiotipo 1 | 0,997 | 0,953 | 0,991 | 0,994 | **1,000** | 0,625 |
| Intra‑quimiotipo 2 | 0,996 | 0,996 | 0,879 | 0,875 | 0,932 | 0,400 |
| Intra‑quimiotipo 3 | **1,000** | 0,825 | 0,740 | 0,880 | 0,928 | 0,563 |
| Entre quimiotipos (1→otros) | **0,718** | 0,503 | 0,467 | 0,530 | 0,608 | 0,452 |
| Entre quimiotipos (2→otros) | 0,675 | 0,734 | 0,733 | 0,548 | **0,784** | 0,550 |
| Entre quimiotipos (3→otros) | 0,750 | 0,665 | 0,622 | 0,707 | 0,718 | 0,447 |

### Enriquecimiento EF20 (inactivos difíciles)

El EF20 dice cuántas veces más activos hay en el primer 20 % de la lista que ordenando al
azar. EF = 1 es azar.

| Escenario | Tanimoto (max) | ErG+GED (max) | ErG+GED (media) | Fusión |
|---|---|---|---|---|
| Mezclado | 1,08 | 1,25 | 1,00 | 1,17 |
| Intra‑quimiotipo 1 | 4,80 | 3,60 | 3,60 | 4,20 |
| Intra‑quimiotipo 2 | 5,22 | 4,48 | 4,48 | 4,48 |
| Intra‑quimiotipo 3 | 5,00 | 2,00 | 4,00 | 4,00 |
| Entre quimiotipos (1→otros) | 0,87 | 0,69 | 0,69 | 0,87 |
| Entre quimiotipos (2→otros) | 1,45 | 1,29 | 0,80 | 1,29 |
| Entre quimiotipos (3→otros) | 1,41 | 0,84 | 1,27 | 1,27 |

### Qué se concluye de esta tabla (sin adornos)

1. **Dentro de un quimiotipo, todos los métodos funcionan bien** (AUC 0,74–1,00; EF20 de
   2 a 5,2). La tarea es fácil: buscar análogos de moléculas muy parecidas.
2. **Con activos mezclados, ningún método funciona** (AUC 0,47–0,61). Confirma el
   hallazgo de la sección 3 con una metodología independiente.
3. **Entre quimiotipos —la prueba realmente justa— ningún método pasa de 0,78** y varios
   caen por debajo de 0,5. Ninguna de las dos aproximaciones generaliza a familias
   químicas no vistas.
4. **Con θ sin calibrar, Morgan+Tanimoto es igual o mejor que ErG+GED en 5 de los 7
   escenarios.** Este es el resultado más incómodo del trabajo y hay que reportarlo tal
   cual. Dos matices legítimos, no excusas: (a) esta evaluación no usa los pesos
   calibrados, y con θ calibrado el AUC de prueba sube a 0,865 en el experimento
   equivalente de la sección 3; (b) los intervalos de confianza del 95 % de los dos
   métodos se solapan en casi todos los escenarios (con 5–36 activos de prueba no hay
   potencia estadística para separarlos).
5. **La fusión por rangos de ErG+GED con Tanimoto es la mejor opción o está muy cerca de
   la mejor en varios escenarios** (0,956 / 0,950 / 0,933 intra‑quimiotipo con inactivos
   difíciles; 0,784 en entre‑quimiotipos 2 con inactivos fáciles). Es el hallazgo
   positivo más defendible: **las dos medidas son complementarias**, no sustitutas.
6. **El filtro geométrico de NS1 no discrimina** (AUC 0,30–0,63) y combinarlo con ErG+GED
   empeora el ranking. Las tasas de aprobación del filtro lo explican: pasa el filtro
   entre el 0 % y el 50 % de los activos y entre el 20 % y el 45 % de los inactivos — no
   separa. La razón es que los criterios describen el bolsillo CPSF30 concreto, mientras
   que los activos de AID 2326 actúan sobre NS1 por mecanismos variados que no pasan
   necesariamente por ese bolsillo.

La complementariedad del punto 5 se confirma de forma independiente: la correlación de
Spearman entre ErG+GED y Tanimoto en los 400 pares cruzados de los dos top‑20 es
**0,166** (`comparacion_A9_A22_resumen_v2.txt`). Las dos medidas están midiendo cosas
distintas.

---

## 5. Aplicación a las 33 013 moléculas

`resultados/top20_A9_thetav2.csv`, `resultados/top20_A22_thetav2.csv`

### Top‑20 frente a A9

| # | ChEMBL ID | ErG+GED | Tanimoto | MW | SMILES |
|---|---|---|---|---|---|
| 1 | CHEMBL1505940 | 0,9895 | 0,477 | 279,3 | `Cc1c(C(=O)NCCc2ccccc2)oc2ccccc12` |
| 2 | CHEMBL1428886 | 0,9895 | 0,479 | 374,3 | `Cc1oc(-c2cc(Cl)cc(Cl)c2)cc1C(=O)NCCc1ccccc1` |
| 3 | CHEMBL1417824 | 0,9895 | 0,422 | 263,7 | `Cc1occc1C(=O)NCCc1ccc(Cl)cc1` |
| 4 | CHEMBL1408523 | 0,9895 | 0,417 | 263,7 | `Cc1occc1C(=O)NCCc1cccc(Cl)c1` |
| 5 | CHEMBL1349270 | 0,9802 | 0,314 | 309,4 | `COc1ccccc1C(=O)NCC1(c2ccccc2)CCCC1` |
| 6 | CHEMBL1446779 | 0,9699 | 0,125 | 357,3 | `CN1/C(=C\C(=O)c2cncc(Br)c2)C(C)(C)c2ccccc21` |
| 7 | CHEMBL1453837 | 0,9541 | 0,328 | 404,3 | `COc1ccc(CCNC(=O)c2cc(-c3ccc(Cl)c(Cl)c3)oc2C)cc1` |
| 8 | CHEMBL1570735 | 0,9457 | 0,396 | 323,4 | `CC1CCN(c2ccc(C(=O)NCCc3ccccc3)cn2)CC1` |
| 9 | CHEMBL1498912 | 0,9457 | 0,412 | 323,4 | `COc1ccc(CCNC(=O)c2cc3ccccc3oc2=O)cc1` |
| 10 | CHEMBL1423811 | 0,9457 | 0,368 | 323,4 | `CC1CCCN(c2ccc(C(=O)NCCc3ccccc3)cn2)C1` |
| 11 | CHEMBL1565294 | 0,9451 | 0,360 | 259,3 | `COc1ccccc1CCNC(=O)c1ccoc1C` |
| 12 | CHEMBL1502464 | 0,9451 | 0,315 | 284,4 | `Cc1ccsc1C(=O)NCCc1c[nH]c2ccccc12` |
| 13 | CHEMBL1331497 | 0,9438 | 0,432 | 285,8 | `O=C(/C=C/c1ccccc1Cl)NCCc1ccccc1` |
| 14 | CHEMBL62044 | 0,9438 | 0,432 | 285,8 | `O=C(/C=C/c1ccc(Cl)cc1)NCCc1ccccc1` |
| 15 | CHEMBL1430612 | 0,9424 | 0,467 | 309,3 | `O=C(NCCc1ccccc1)c1ccc(-c2ccc(F)cc2)o1` |
| 16 | CHEMBL1522952 | 0,9366 | 0,444 | 277,8 | `Cc1cc(C(=O)NCCc2ccc(Cl)cc2)c(C)o1` |
| 17 | CHEMBL1545664 | 0,9339 | 0,136 | 317,3 | `O=C(/C=C/c1ccco1)NC(=S)Nc1cccc([N+](=O)[O-])c1` |
| 18 | CHEMBL1446411 | 0,9303 | 0,346 | 339,2 | `O=C(NCCc1cccc(Cl)c1)c1ccc(Cl)c([N+](=O)[O-])c1` |
| 19 | CHEMBL1399595 | 0,9303 | 0,381 | 271,4 | `O=C(NC12CC3CC(CC(C3)C1)C2)c1ccccc1O` |
| 20 | CHEMBL1565959 | 0,9303 | 0,525 | 304,2 | `O=C(NCCc1ccccc1)c1cccc(Br)c1` |

### Lo que hay que ver en esta tabla

**El método recupera el quimiotipo correcto.** A9 es
`O=C(NCCc1ccccc1)c1ccccc1O`: una **arilamida de fenetilamina**. De las 20 moléculas
seleccionadas, **17 son arilamidas de fenetilamina o variantes muy cercanas** — el mismo
esqueleto, con el anillo aromático sustituido o reemplazado por un heterociclo (furano,
benzofurano, tiofeno, indol, piridina). Sin ninguna información sobre la proteína más
allá de dos moléculas de referencia, el ranking reconstruyó la familia química.

**Y lo hace con moléculas que un fingerprint habría dejado más abajo.** Los Tanimoto del
top‑20 van de **0,125 a 0,525**, con mediana cerca de 0,40. En términos de fingerprint
esas moléculas no son "muy parecidas" a A9; en términos de farmacóforo sí lo son. Los
casos más llamativos son el puesto 6 (Tanimoto 0,125) y el 17 (0,136): comparten la
organización farmacofórica pero muy pocas subestructuras exactas. Eso es exactamente lo
que el método debía aportar sobre Tanimoto, y lo aporta.

**Advertencia sobre el corte.** Hay empates: cuatro moléculas comparten el puntaje
máximo (0,9895) y tres el del puesto 18–20 (0,9303). El corte exacto en el puesto 20 es
en parte arbitrario; lo que el método entrega es un **vecindario químico**, no un orden
estricto entre las 20.

### Corrección de un resultado del anteproyecto

El anteproyecto reporta que "el top‑20 de A9 y el top‑20 de A22 **no comparten ninguna
molécula**". **Ese resultado era un artefacto del GED defectuoso** y queda corregido:

| | Moléculas en común entre el top‑20 de A9 y el de A22 |
|---|---|
| GED v1 (defectuoso) | **0 / 20** |
| GED corregido + θ recalibrado | **17 / 20** |

Y así tenía que ser: los grafos ErG de A9 y A22 son el mismo grafo, así que las dos
listas deberían ser idénticas. Las 3 moléculas de diferencia se explican por los empates
del borde y por el carácter *greedy* del refinamiento del emparejamiento, que puede
llegar a óptimos locales distintos según el orden de numeración de los nodos (los átomos
de A9 y A22 vienen en orden distinto).

Esta corrección **refuerza** la conclusión del anteproyecto en lugar de debilitarla: la
recomendación de reportar siempre por activo de referencia sigue siendo válida —y es
imprescindible cuando los activos son de quimiotipos distintos, como los tres de
dengue—, pero el "0/20" no era evidencia de ello, era un error de implementación.

### Comparación de los dos top‑20 (listas del GED corregido con θ anterior)

`resultados/comparacion_A9_A22_resumen_v2.txt`

| | ErG+GED | Tanimoto |
|---|---|---|
| Similitud media cruzada (lista A9 × lista A22) | 0,727 | 0,162 |
| Similitud máxima cruzada | 0,891 | 0,455 |
| Similitud media interna de la lista de A9 | 0,808 | 0,246 |
| Similitud media interna de la lista de A22 | 0,734 | 0,187 |
| Similitud media de la lista de A9 contra A9 | 0,843 | 0,314 |
| Similitud media de la lista de A22 contra A22 | 0,798 | 0,207 |
| Correlación de Spearman entre las dos medidas (400 pares) | **0,166** | |

La similitud cruzada (0,727) es comparable a la interna de cada lista (0,808 y 0,734):
los dos top‑20 son **el mismo vecindario químico**, no dos vecindarios distintos.

---

## 6. Eficiencia computacional

`resultados/tiempos_por_bloque_v2.csv` — 16 bloques de 1 500 moléculas.

| Etapa | ms por molécula |
|---|---|
| Embedding 3D (ETKDG + UFF) | **19,6** |
| Construcción del grafo ErG | **20,5** |
| GED frente a las 2 referencias | **1,7** |
| **Total (incluye lectura del Excel y Tanimoto)** | **44,3** |
| Por comparación GED individual | **0,86** |

Tiempo total del ranking completo: **≈ 24 minutos** para 33 013 moléculas en un portátil,
sin paralelizar.

Dos lecturas que importan:

1. **El GED no es el cuello de botella.** Representa el 4 % del tiempo. El 90 % se lo
   llevan el embedding 3D y la construcción del ErG, que son *preprocesamiento* y se
   hacen **una sola vez por molécula** independientemente de cuántas referencias haya.
   Añadir una tercera referencia costaría 0,86 ms más por molécula (+2 %).
2. **Es el argumento de eficiencia del método.** Un prefiltro de 0,86 ms por comparación
   frente a minutos u horas por molécula del docking molecular: tres a cinco órdenes de
   magnitud de diferencia. Para eso sirve un prefiltro.

Si el preprocesamiento se guardara en caché (los grafos ErG en disco), un cribado
posterior contra nuevas referencias costaría solo la parte del GED: **menos de 30
segundos para las 33 000**.

---

## 7. Resumen de limitaciones

Todas verificadas con los datos de este repositorio:

1. **El método es una búsqueda de análogos, no un predictor de actividad.** AUC 0,48 con
   activos mezclados frente a 0,88–0,98 dentro de un quimiotipo.
2. **No generaliza a familias químicas no vistas.** Entre quimiotipos, ningún método
   evaluado pasa de 0,78 y varios caen por debajo del azar.
3. **Con θ sin calibrar, Morgan+Tanimoto es competitivo o mejor.** El aporte defendible
   del método es la **complementariedad** (fusión por rangos; Spearman 0,166), no la
   superioridad.
4. **El GED es aproximado y su refinamiento es *greedy*.** Depende del orden de
   numeración de los nodos: 17/20 de coincidencia entre los top‑20 de dos grafos
   idénticos, en lugar de 20/20.
5. **La reducción ErG pierde resolución por diseño.** 22 064 de las 33 013 moléculas
   reciben exactamente el mismo puntaje frente a A9 y frente a A22, y hay empates en el
   borde del top‑20.
6. **Un solo conformero por molécula.** Una molécula flexible puede quedar plegada y
   subestimar su tamaño real; de ahí que el criterio de tamaño se calcule aparte con 30
   conformeros y no entre al filtro por defecto.
7. **El filtro geométrico de NS1 no aportó.** AUC 0,30–0,63; pasa el filtro un 0–50 % de
   los activos y un 20–45 % de los inactivos.
8. **Las dos referencias de Fase A son una sola para el método**, porque sus grafos ErG
   son idénticos. La validación con `n = 2` referencias es, en realidad, con una.
9. **Pocos activos en cada experimento** (5 a 36 en prueba): los intervalos de confianza
   del 95 % son anchos y se solapan entre métodos.
10. **No hay validación experimental.** Ninguna de las 20 moléculas priorizadas se ha
    ensayado contra NS1. El resultado es una hipótesis priorizada, no un hallazgo
    biológico.

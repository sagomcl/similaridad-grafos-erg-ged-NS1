# 01 — Explicación del método, de cero

Este documento explica **qué hace el método y por qué**, sin suponer que quien lo lee
sabe quimioinformática. Está escrito en el orden en que ocurren las cosas dentro del
código. La explicación línea por línea del código está en
[`02_explicacion_del_codigo.md`](02_explicacion_del_codigo.md).

---

## 1. El problema, en una frase

Hay una proteína (NS1) que interesa bloquear. Se conocen **muy pocas** moléculas que
la bloqueen. Existe una biblioteca de **33 047 moléculas candidatas** y la pregunta es:
¿cuáles de esas 33 047 vale la pena llevar al laboratorio primero?

Probar una molécula en el laboratorio cuesta dinero y semanas. Hacer *docking*
molecular (simular el acoplamiento físico con la proteína) cuesta horas de cómputo por
molécula. Lo que hace este proyecto es un **prefiltro**: ordenar las 33 047 en
milisegundos por molécula, según cuánto se parecen —**en lo que importa para unirse a
la proteína**— a los inhibidores que ya se sabe que funcionan.

La frase clave es "en lo que importa". Todo el método consiste en definir
operacionalmente ese "lo que importa".

---

## 2. Qué es un grupo farmacofórico

Cuando una molécula pequeña (un fármaco) se une a una proteína, no se une "toda". Se
une por unos pocos puntos de contacto. Esos puntos son de cuatro clases, y son las
cuatro únicas cosas que una molécula puede ofrecerle a una proteína:

| Clase | Qué hace físicamente | Ejemplo en una molécula |
|---|---|---|
| **HBD** (donador de puente de hidrógeno) | tiene un hidrógeno unido a N o a O y lo "presta" | el –OH de un fenol, el –NH– de una amida |
| **HBA** (aceptor de puente de hidrógeno) | tiene un par de electrones libres y recibe ese hidrógeno | el O del grupo C=O, el N de una piridina |
| **Aromático** | anillo plano con electrones deslocalizados; se apila cara a cara con los anillos de la proteína (interacción π–π) | un anillo de benceno |
| **Hidrofóbico** | zona apolar que "huye" del agua y se mete en las cavidades grasas de la proteína | una cadena –CH₂–CH₃, un metilo |

El conjunto de estos puntos, **con su disposición en el espacio**, se llama
**farmacóforo**. La idea central de la farmacología molecular es que *dos moléculas con
el mismo farmacóforo se unen al mismo sitio, aunque sus átomos sean distintos*.

Esto es importante y vale la pena subrayarlo: la proteína **no ve átomos**, ve
funciones. Le da igual si el donador de hidrógeno es un –OH o un –NH₂; lo que le
importa es que haya un donador en ese lugar. El método se construye exactamente sobre
esa idea.

### Qué es un anillo aromático

Un anillo aromático es un anillo (normalmente de seis carbonos, como el benceno) cuyos
electrones no están fijos en enlaces individuales sino repartidos por todo el anillo,
por encima y por debajo del plano. Eso tiene dos consecuencias que importan aquí:

1. **Es plano y rígido.** No se dobla. Su geometría es predecible.
2. **Esa nube de electrones puede interactuar con la nube de otro anillo aromático de
   la proteína**, apilándose como dos cartas. Es una de las interacciones más fuertes y
   más frecuentes en el reconocimiento proteína–ligando.

Por eso los anillos aromáticos son el ancla principal en la mayoría de los bolsillos de
unión, y por eso en los criterios de NS1 (sección 7) se exige que haya **al menos dos**
y que estén **en los extremos** de la molécula.

---

## 3. Qué es un nodo del ErG, y por qué se colapsan los átomos

### El grafo molecular completo (y por qué no sirve)

La representación natural de una molécula como grafo es: **un nodo por átomo, una
arista por enlace**. Para el inhibidor A9 eso son 18 nodos y 18 aristas (sin contar
hidrógenos).

El problema es doble:

- **Demasiado detalle para la pregunta.** El grafo completo distingue un cloro de un
  bromo, o un anillo de seis carbonos de uno de cinco con un azufre. Pero la proteína
  no distingue eso: ve "zona apolar voluminosa" en los dos casos. El grafo completo
  codifica información que no es relevante para predecir unión, y esa información
  irrelevante actúa como ruido: hace que dos moléculas que se unirían igual parezcan
  diferentes.
- **Demasiado costoso.** Comparar dos grafos con la distancia de edición es un problema
  NP‑difícil: el costo crece factorialmente con el número de nodos. Con 18 nodos por
  molécula y 33 047 moléculas, no se termina nunca.

### La reducción: de átomos a grupos funcionales

El **Extended Reduced Graph** (ErG, Stiefl et al. 2006) resuelve las dos cosas con una
sola idea: **un nodo por grupo funcional, no por átomo**.

El procedimiento, tal como lo implementa `build_erg()`:

1. **Tipificar cada átomo.** Se le pregunta a RDKit, usando su catálogo de reglas
   SMARTS (`BaseFeatures.fdef`), a qué familias farmacofóricas pertenece cada átomo.

2. **Resolver los solapamientos.** Aquí aparece el primer problema real del proyecto:
   un mismo átomo pertenece a varias familias a la vez. El oxígeno de un fenol es
   donador **y** aceptor. Los carbonos de un benceno son aromáticos **y** también
   hidrofóbicos. Si se creara un nodo por cada familia, el grafo **crecería** en lugar
   de reducirse — exactamente lo contrario del propósito. Fue el primer error que hubo
   que corregir en este proyecto: un solo anillo de benceno terminaba representado por
   un nodo aromático *más* varios nodos hidrofóbicos.

   La regla que se adoptó: **un solo tipo por átomo**, con la prioridad

   ```
   HBD  >  HBA  >  Aromático  >  Hidrofóbico
   ```

   El orden va de lo más específico e informativo a lo más genérico. Ser donador de
   puente de hidrógeno es una propiedad muy particular y muy determinante; ser apolar
   es lo que le pasa a casi cualquier carbono.

3. **Colapsar los contiguos del mismo tipo.** Se recorre la molécula por sus enlaces
   (búsqueda en anchura) y todos los átomos vecinos que comparten tipo se fusionan en
   **un solo nodo**. Los seis carbonos aromáticos de un benceno se vuelven un nodo
   "Aromático". Una cadena de tres carbonos apolares se vuelve un nodo "Hidrofóbico".

   Los átomos que no pertenecen a ninguna de las cuatro familias (por ejemplo un
   carbono de la cadena que une dos grupos) **no generan nodo**: son pegamento
   estructural. Están en la molécula, cuentan para las distancias, pero no son puntos
   de contacto.

4. **Resultado.** A9 pasa de 18 átomos a **6 nodos**: dos aromáticos, dos donadores, un
   aceptor y un hidrofóbico.

### Por qué esto no es "perder información"

Es la pregunta más razonable que se le puede hacer al método, y tiene tres respuestas.

**Primera: la información atómica no se borra, se guarda.** Cada nodo del grafo
almacena el atributo `atom_ids`, la lista exacta de los átomos que se colapsaron en él.
Siempre se puede ir del nodo de vuelta a sus átomos. El script
`comparar_A9_A22.py` lo usa para imprimir, por cada nodo emparejado, los símbolos
químicos de los átomos que lo componen —así se lee `Aromatico[CCCCCC]` o `HBD[O]`—; y
`evaluar_erg_ns1.py` lo usa para recuperar las coordenadas de los átomos pesados y
medir el tamaño real de la molécula. El ErG es una **vista agregada** de la molécula,
no una destrucción de la molécula: la molécula completa sigue ahí, con su SMILES, su
conformero 3D y sus átomos.

**Segunda: la información que se descarta es, por diseño, la que no importa.** El
método descarta deliberadamente la identidad del átomo cuando no cambia la función. El
ejemplo perfecto está en el propio proyecto: **A9 y A22 difieren en un solo átomo**
(A9 tiene –OH donde A22 tiene –NH₂). Con el ErG los dos grafos resultan **idénticos**
(6 nodos, 14 aristas, similitud 1,000), porque –OH y –NH₂ son los dos donadores de
puente de hidrógeno. Y, en efecto, **las dos moléculas se unen al mismo bolsillo de
NS1 con potencia parecida** (Kleinpeter et al. 2018): la representación acertó. Un
fingerprint de Morgan, en cambio, les da similitud 0,706, porque sí distingue el
oxígeno del nitrógeno — distingue algo que a la proteína le da igual.

**Tercera: la reducción es lo que hace el cálculo posible.** Con 6 nodos en vez de 18,
la distancia de edición de grafos pasa de intratable a ~0,2 ms por comparación. Sin la
reducción no habría método.

El precio de esta decisión hay que decirlo con claridad, porque es una limitación real:
el ErG **no puede distinguir dos moléculas que compartan farmacóforo pero difieran en
algo que sí afecte la actividad** (un cloro que mejora la afinidad, un
estereocentro invertido, un anillo de cinco frente a uno de seis). En el ranking de las
33 047 moléculas, 22 064 reciben exactamente el mismo puntaje frente a A9 y frente a
A22, y hay empates en el borde del top‑20. El método propone *vecindarios químicos*,
no moléculas individuales.

---

## 4. Cómo entra la geometría tridimensional

Esto merece una sección propia porque es donde más fácil es confundirse, y porque
conviene ser preciso en la tesis.

Para cada molécula se genera un **conformero 3D**: una geometría concreta en el
espacio, con coordenadas (x, y, z) para cada átomo. Se hace con **ETKDG v3** (geometría
de distancias con preferencias de torsión tomadas de estructuras cristalográficas
reales, de modo que el conformero no solo es geométricamente válido sino químicamente
plausible) y después se relaja con el campo de fuerzas **UFF**. La semilla aleatoria
está fija (`seed=0xF00D`), así que la misma molécula produce siempre el mismo
conformero: el método es reproducible.

Esas coordenadas 3D se usan **en dos lugares y no en un tercero**:

**(a) Sí: el centro de cada nodo.** Cada nodo guarda `center`, el centroide de las
coordenadas 3D de sus átomos. Es la posición física del grupo farmacofórico en el
espacio, en ángstroms.

**(b) Sí: los criterios geométricos del bolsillo.** Los `center` y las coordenadas de
los átomos pesados alimentan los cinco criterios de `ns1_utils.criterios_ns1` (sección
7): la distancia real entre los dos anillos aromáticos más alejados, el largo total de
la molécula en ángstroms, y los ángstroms por enlace de cada arista. Ahí la geometría
3D es el dato principal.

**(c) No: las aristas del grafo y el cálculo del GED.** Las aristas del ErG **no** están
etiquetadas con la distancia en ángstroms, sino con la **distancia topológica**: el
número de enlaces del camino más corto entre los dos grupos (de 1 a 6). Y el costo de
editar una arista es `edge_beta × |diferencia en número de enlaces|`.

**Por qué se tomó esa decisión.** Una molécula flexible tiene miles de conformeros
posibles. Si las aristas llevaran ángstroms, el grafo —y por tanto el ranking— cambiaría
con cada conformero que se generara: la misma molécula embebida con otra semilla daría
otro puntaje. La distancia en **enlaces**, en cambio, es una propiedad de la molécula,
no de la pose: es invariante. Sacrificar precisión geométrica a cambio de
reproducibilidad y estabilidad es lo correcto para un prefiltro que tiene que ordenar
decenas de miles de moléculas de forma consistente. El refinamiento geométrico es lo
que hace el *docking* después, sobre las 20 moléculas que el prefiltro seleccionó.

**Cómo se recupera la dimensión 3D entonces.** De dos maneras complementarias:

1. El número de enlaces es un *proxy* acotado de la distancia física: entre 1,5 y 2,0 Å
   por enlace si la cadena está extendida. El criterio 5 de NS1 mide exactamente ese
   cociente (distancia 3D real entre centros ÷ número de enlaces) y permite detectar
   cuándo la molécula está plegada y el proxy deja de valer.
2. Los criterios geométricos (largo mínimo, aromáticos en los extremos) sí usan
   ángstroms reales y actúan como **filtro** sobre el ranking: una molécula puede tener
   el farmacóforo correcto en términos topológicos y quedar descartada porque en el
   espacio no alcanza a cubrir el canal del bolsillo.

Si en un trabajo futuro se quisiera un GED propiamente tridimensional, el cambio está
localizado: usar `|dist3D(i,j) − dist3D(i',j')|` en lugar de la diferencia en enlaces
dentro de `_total_cost` y `_local_edge_cost` de `ged_mejorado.py`, generando varios
conformeros por molécula y quedándose con la mejor coincidencia. Está señalado en
[`03_como_extrapolar.md`](03_como_extrapolar.md).

---

## 5. Cómo se comparan dos grafos: la distancia de edición

Ya tenemos dos moléculas convertidas en grafos pequeños. ¿Cómo se mide si se parecen?

La **distancia de edición de grafos** (Graph Edit Distance, GED) responde con una
pregunta muy intuitiva: *¿cuál es el conjunto más barato de modificaciones que convierte
el grafo 1 en el grafo 2?* Las modificaciones permitidas son seis, y cada una tiene un
precio:

| Operación | Precio | Significado químico |
|---|---|---|
| sustituir el tipo de un nodo: aromático ↔ hidrofóbico | `w_aromatico_hidrofobo` | confundir un anillo con una cadena grasa; los dos son apolares, error leve |
| sustituir el tipo de un nodo: HBA ↔ HBD | `w_hba_hbd` | confundir aceptor con donador; mismo "idioma" polar, papel invertido |
| sustituir el tipo de un nodo: polar ↔ apolar | `w_polar_apolar` | confundir un grupo polar con uno apolar; el error químicamente **más grave** |
| insertar o borrar un nodo | `node_indel` | una molécula tiene un grupo funcional que la otra no tiene |
| estirar o encoger una arista | `edge_beta` × (enlaces de diferencia) | los dos grupos están, pero a distinta separación |
| insertar o borrar una arista | `edge_indel` | la relación entre dos grupos existe en una molécula y no en la otra |

Esos seis números son el vector **θ** (*theta*). Son el único lugar del método donde
entra criterio químico ajustable, y **no se eligen a mano**: se calibran con datos
(sección 8).

La suma de los precios de todas las operaciones necesarias es el GED.

### Normalización

El GED crudo tiene un sesgo: dos moléculas grandes acumulan más costo que dos pequeñas
simplemente porque hay más nodos y aristas que editar. Sin corregirlo, el ranking
favorecería a las moléculas pequeñas. Por eso se divide por el número total de
elementos de los dos grafos:

```
GED_normalizado = GED / (n₁ + n₂ + m₁ + m₂)
```

donde n son nodos y m aristas. El denominador es el número máximo de operaciones de
edición posibles, así que el cociente responde a "qué fracción de todo lo que había para
editar hubo que editar".

### De distancia a similitud

```
similitud = exp( − GED_normalizado )
```

Esta transformación lleva un costo en [0, ∞) a una similitud en (0, 1], con 1 para
grafos idénticos. Se usa la exponencial y no `1 − x` por dos razones: no hace falta
recortar valores (nunca sale negativa), y penaliza mucho las primeras diferencias y
cada vez menos las siguientes — que es lo adecuado, porque entre dos moléculas ya muy
distintas da lo mismo cuánto más distintas sean.

---

## 6. Dónde entra el algoritmo húngaro

Hay un detalle que se pasó por alto arriba: para saber qué nodo hay que sustituir por
cuál, **primero hay que decidir qué nodo del grafo 1 se corresponde con qué nodo del
grafo 2**. Y ese emparejamiento es justamente lo difícil.

Si los dos grafos tienen *n* nodos, hay *n!* emparejamientos posibles. Probarlos todos
y quedarse con el más barato es el GED exacto, y es **NP‑difícil**: con 6 nodos son 720
combinaciones (aún viable), pero el costo crece factorialmente y para moléculas más
grandes se vuelve imposible — y hay que hacerlo 33 047 veces.

**Riesen y Bunke (2009)** propusieron la aproximación que usa este trabajo: plantear el
emparejamiento como un **problema de asignación**.

1. Se construye una **matriz de costos** cuadrada. La entrada (i, j) dice cuánto costaría
   decir que el nodo i del grafo 1 "es" el nodo j del grafo 2. Como los grafos pueden
   tener distinto número de nodos, la matriz se completa con filas y columnas de nodos
   ficticios cuyo costo es `node_indel`: emparejar un nodo real con uno ficticio
   equivale a borrarlo o insertarlo.

2. El problema es entonces: *elegir exactamente una pareja para cada fila, sin repetir
   columnas, minimizando la suma total*. Ese es el **problema de asignación óptima**, y
   el **algoritmo húngaro** (Kuhn 1955; Munkres) lo resuelve de forma **exacta** en
   tiempo **O(N³)**.

   En el código es una sola línea:
   ```python
   r, c = linear_sum_assignment(Cn)   # scipy.optimize
   ```

**Ese es el aporte del húngaro al proyecto**: convierte un problema combinatorio
intratable en uno que corre en microsegundos, y lo hace sin heurísticas — el
emparejamiento que devuelve es el óptimo *de la matriz que se le dio*. El GED resultante
es una cota superior del GED exacto.

### El defecto que hubo que corregir

La primera implementación armaba la matriz de costos **mirando solo los tipos de nodo**.
Eso tiene un problema: cuando hay dos nodos del mismo tipo (dos anillos aromáticos, dos
donadores), las entradas de la matriz **empatan**, y el húngaro devuelve cualquiera de
los emparejamientos óptimos — a veces uno "cruzado", que empareja el aromático de la
izquierda de una molécula con el de la derecha de la otra. El costo de aristas se
calculaba *después*, sobre ese emparejamiento ya fijado, y salía malísimo.

El síntoma fue inconfundible: **A9 y A22, con grafos ErG idénticos, daban similitud
0,43 en vez de 1,00**.

La corrección (`ged_mejorado.py`) es la estándar en la literatura de GED bipartito:

1. **Matriz de costos aumentada.** Al costo de emparejar i con j se le suma la mitad del
   costo de alinear sus "huellas locales de aristas": la lista de distancias a sus
   vecinos. Así, entre dos nodos aromáticos empatados por tipo, el húngaro prefiere el
   que además está a distancias parecidas de los demás grupos. (El factor ½ evita contar
   dos veces cada arista, una por cada extremo.)
2. **Refinamiento local.** Se intercambian pares de asignaciones mientras el costo total
   (nodos + aristas) siga bajando.
3. **Se prueban los dos puntos de partida** (matriz de solo tipos y matriz aumentada) y
   se conserva el mejor, de modo que el resultado nunca es peor que el de la versión
   anterior.

Después de la corrección: A9 vs A22 = **1,000**, con los 6 nodos emparejados a coste
cero y 0 diferencias de distancia.

**Limitación que hay que declarar:** el paso 2 es una búsqueda local (*greedy*), así que
su resultado puede depender del orden en que estén numerados los nodos. Dos grafos
idénticos pero con numeración distinta —como A9 y A22, cuyos átomos vienen en otro
orden— pueden terminar con similitudes que difieren en el tercer decimal frente a una
misma molécula candidata. Es un óptimo local, no el GED exacto. En la práctica esto
explica por qué los top‑20 de A9 y de A22 coinciden en 17 de 20 y no en 20 de 20, a
pesar de que sus grafos son el mismo grafo.

---

## 7. Qué se le pide a una molécula para poder inhibir NS1

Hasta aquí el método mide *parecido a un inhibidor conocido*. Pero se puede además
imponer directamente lo que el bolsillo necesita, sin pasar por la similitud.

El sitio de NS1 donde se unen A9 y A22 es el **bolsillo CPSF30** (Kleinpeter et al.
2018): un canal alargado, mayoritariamente hidrofóbico, con dos cavidades en los
extremos donde encajan anillos aromáticos y puntos de anclaje polar en el medio. De esa
descripción salen cinco criterios, implementados en `ns1_utils.criterios_ns1`:

1. **Composición.** Al menos 2 nodos aromáticos, ≥1 hidrofóbico, ≥1 donador y ≥1
   aceptor. Es decir: la molécula tiene con qué llenar el canal apolar **y** con qué
   formar puentes de hidrógeno.
2. **Conectividad.** El grafo ErG debe quedar **conexo**. Como solo se dibujan aristas
   entre grupos separados por 1–6 enlaces, un grafo desconectado significa que hay
   grupos a más de 6 enlaces de todos los demás: la molécula es demasiado laxa y sus
   grupos no pueden ocupar el bolsillo al mismo tiempo.
3. **Aromáticos en los extremos.** La separación 3D entre los dos aromáticos más
   alejados debe ser al menos el 70 % de la separación máxima entre cualquier par de
   grupos. Traduce la idea de que los anillos van en las puntas, no apilados en el
   centro.
4. **Tamaño.** El largo de la molécula (distancia máxima entre dos átomos pesados) debe
   superar 15 Å, para que alcance a cubrir el canal.
5. **Ajuste de distancias.** Por cada arista: (distancia 3D entre centros) ÷ (número de
   enlaces) debe caer entre 1,5 y 2,0 Å por enlace. Si es mucho menor, la cadena está
   plegada sobre sí misma y los grupos no están donde el grafo sugiere. Se reporta la
   fracción de aristas que cumplen.

Los criterios 1–3 forman el filtro "duro" (`estructura`). El 4 se deja aparte porque un
único conformero puede salir plegado y subestimar el largo real: por eso existe
`largo_multiconf`, que genera 30 conformeros y toma el máximo, y por eso el criterio de
tamaño solo entra al filtro con la bandera `--exigir-tamano`.

**Resultado honesto de este filtro:** aplicado sobre activos reales de PubChem, el
filtro NS1 por sí solo tiene un AUC de 0,30–0,63 según el escenario — es decir, apenas
discrimina, y combinarlo con ErG+GED empeora el ranking en vez de mejorarlo. La
explicación es que los criterios describen el bolsillo CPSF30 concreto, mientras que los
activos del cribado AID 2326 actúan sobre NS1 por mecanismos variados que no pasan
necesariamente por ese bolsillo. Es un resultado negativo informativo, y está reportado
como tal.

---

## 8. Cómo se calibran los pesos θ (y por qué eso es la parte delicada)

Los seis pesos de θ no se pueden elegir a ojo: son el criterio químico del método y
cambian el ranking por completo. Se **calibran con datos**.

### El procedimiento

1. Se necesita un conjunto **etiquetado**: moléculas que se sabe que son activas y
   moléculas que se sabe que no lo son. En este proyecto es **PubChem AID 2326**, el
   cribado de alto rendimiento del NIH Chemical Genomics Center contra la función de
   NS1 de influenza — el mismo cribado que dio origen a A9 y A22. Tiene **270 315
   compuestos: 1 073 activos confirmados y 260 684 inactivos confirmados**. Es el
   *ground truth* de todo el trabajo.
2. Se parte en **entrenamiento** y **prueba**.
3. Se busca el θ que maximice el **AUC** (área bajo la curva ROC) en el conjunto de
   entrenamiento. El AUC es la probabilidad de que un activo tomado al azar reciba un
   puntaje mayor que un inactivo tomado al azar: 0,5 es azar, 1,0 es separación
   perfecta.
4. Se reporta el AUC en el conjunto de **prueba**, que nunca se vio durante la
   calibración. El puntaje de una molécula de prueba se calcula como su similitud a los
   activos de **entrenamiento**: así se reproduce la situación real de un cribado, donde
   se conocen unos pocos inhibidores y hay que ordenar una biblioteca desconocida.

La búsqueda se hace con **evolución diferencial** (Storn & Price 1997), un optimizador
global por poblaciones, y no con un método de gradiente. La razón es que el AUC es una
función **escalonada** de θ: solo cambia cuando dos moléculas intercambian su posición
en el orden, así que su gradiente es cero casi en todas partes y los métodos de descenso
no sirven.

### El hallazgo metodológico central del proyecto

Esta parte es el resultado más importante de la Fase A, y conviene contarla como
sucedió:

- Calibrando **solo contra A9 y A22** (n = 2 activos), el AUC de calibración fue
  **1,000**. Parecía un éxito. No lo era: con dos activos químicamente casi idénticos,
  cualquier método razonable puede "memorizar" la diferencia frente a los decoys sin
  tener ninguna capacidad real de generalizar.
- Al recalibrar con los **1 073 activos reales** de AID 2326, con partición
  entrenamiento/prueba, el AUC de prueba fue **0,49**: azar.
- La explicación apareció al agrupar esos 1 073 activos por similitud estructural con el
  algoritmo de **Butina** (1999): salieron **737 grupos, de los cuales 588 son de una
  sola molécula**. Es decir, el cribado fenotípico mezcla compuestos que inhiben NS1 por
  mecanismos completamente distintos, sin un farmacóforo común que el método pueda
  aprender a reconocer de forma agregada.
- Al repetir la validación **dentro** de los tres quimiotipos más grandes, el AUC de
  prueba subió a **0,96 / 0,91 / 0,88** (con el GED corregido).

**Conclusión:** ErG+GED+húngaro sirve para priorizar **análogos estructurales de un
inhibidor conocido** — una búsqueda de "primos químicos" —, no para detectar actividad
biológica en general cuando los activos de referencia no comparten farmacóforo. De ahí
la regla que se aplica en todo el proyecto: **reportar siempre la similitud frente a
cada activo de referencia por separado, nunca promediada**.

---

## 9. Contra qué se compara el método

Un método nuevo solo se puede defender contra una línea base. Las líneas base usadas
son:

- **Morgan (ECFP4) + Tanimoto**: el estándar de la industria. Codifica la molécula como
  un vector de 1024 bits donde cada bit indica la presencia de una subestructura
  circular. Mide coincidencia de **subestructuras exactas**.
- **Morgan + Hamming**: incluida porque el anteproyecto la contemplaba. Premia también
  los ceros coincidentes, lo que con vectores dispersos infla la similitud de cualquier
  par; es una medida pobre y los resultados lo confirman.
- **ErG‑fingerprint + Tanimoto generalizado**: el fingerprint ErG que trae RDKit es un
  vector de flotantes que cuenta pares de tipos farmacofóricos a cada distancia
  topológica. Es la comparación más informativa de todas, porque usa la **misma** idea
  de reducción farmacofórica que el método propuesto **pero sin grafo explícito ni
  GED**. La diferencia entre ErG‑fp y ErG+GED aísla exactamente cuánto aporta el grafo.
- **Fusión por rangos** de ErG+GED con Tanimoto: se convierten los dos puntajes a rangos
  en [0, 1] y se promedian. Mide si los dos métodos son complementarios.

Los resultados de esta comparación están en
[`04_resultados.md`](04_resultados.md), incluido el resultado incómodo pero honesto:
con θ sin calibrar, **Tanimoto es igual o mejor que ErG+GED en casi todos los
escenarios**, y la fusión de ambos es la mejor opción en varios de ellos.

---

## 10. Resumen del flujo completo

```
SMILES  (texto: "O=C(NCCc1ccccc1)c1ccccc1O")
  │
  ├─ MolFromSmiles ─────────> grafo molecular (18 átomos, 18 enlaces)
  │
  ├─ AddHs + ETKDGv3 + UFF ─> conformero 3D con coordenadas (x,y,z)
  │
  ├─ feat_factory ──────────> cada átomo recibe familias farmacofóricas
  │                           (con solapamientos)
  ├─ PRIORITY ──────────────> un solo tipo por átomo
  │
  ├─ BFS sobre enlaces ─────> átomos contiguos del mismo tipo se COLAPSAN
  │                           → 6 nodos, cada uno con label, atom_ids, center 3D
  │
  └─ distancia topológica ──> aristas entre grupos a 1–6 enlaces
                              → GRAFO ErG (6 nodos, 14 aristas)

  dos grafos ErG
  │
  ├─ matriz de costos desde θ (tipos + huella local de aristas)
  ├─ ALGORITMO HÚNGARO ─────> emparejamiento óptimo de nodos
  ├─ refinamiento local ────> se corrigen emparejamientos cruzados
  ├─ costo de aristas ──────> edge_beta × |Δenlaces|, edge_indel
  ├─ normalizar ────────────> GED / (n₁+n₂+m₁+m₂)
  └─ exp(−·) ───────────────> SIMILITUD en (0,1]

  θ se calibra con evolución diferencial maximizando AUC sobre PubChem AID 2326,
  con partición entrenamiento/prueba y dentro de quimiotipos coherentes (Butina).

  33 047 candidatas × 2 referencias → ranking → top‑20 por referencia
  (opcionalmente filtrado por los criterios geométricos del bolsillo)
```

# 03 — Cómo extrapolar el método a otra proteína

El método no está atado a NS1. Lo que lo ata a una diana concreta son **cinco bloques
bien localizados**; todo lo demás (construcción del ErG, GED, algoritmo húngaro,
calibración, ranking, métricas) es genérico.

Este documento es la receta operativa. La justificación conceptual de cada decisión está
en [`01_explicacion_del_metodo.md`](01_explicacion_del_metodo.md).

---

## Resumen: los cinco puntos de ajuste

| # | Qué se cambia | Dónde | ¿Obligatorio? |
|---|---|---|---|
| 1 | **Inhibidores de referencia** | `pipeline_final_completo.py` → `INFLUENZA_ACTIVOS` / `DENGUE_ACTIVOS` | **sí** |
| 2 | **Pesos θ** (recalibrar, no editar a mano) | `recalibrar_v2.py` → nuevo `theta_*.json` | **sí**, si existe un conjunto etiquetado |
| 3 | **Criterios geométricos del bolsillo** | `ns1_utils.py` → `PARAMS` y `criterios_ns1` | opcional pero recomendable |
| 4 | **Tipos de nodo** del ErG | `LABEL_MAP`, `PRIORITY`, `TYPES` + `_cost_matrix` | solo si la diana exige tipos nuevos (p. ej. cargas) |
| 5 | **Rango de conexión** de aristas y **umbral de quimiotipo** | `build_erg` (`1 <= d <= 6`), `clusterizar_activos` (`cutoff_distancia`) | solo en casos especiales |

---

## Paso 1 — Cambiar los inhibidores de referencia

Es el único cambio estrictamente imprescindible.

```python
# pipeline_final_completo.py, Sección 3
MI_DIANA_ACTIVOS = {
    "inhibidor_1": "SMILES...",
    "inhibidor_2": "SMILES...",
}
```

**Reglas que hay que respetar**, y que son el resultado principal de la Fase A:

- **Un ranking por referencia, nunca promediado.** Si los inhibidores conocidos no
  comparten farmacóforo, promediar la similitud sobre todos ellos destruye la señal (es
  literalmente lo que hizo caer el AUC a 0,49). Usar siempre `por_activo=True` y, en la
  evaluación, la variante **"max"** (similitud al activo más parecido) en lugar de la
  media.
- **Verificar la procedencia de cada SMILES.** Conviene anotar en el propio código la
  fuente (DOI, PDB, nombre exacto del compuesto) como se hizo con A9/A22 y con los tres
  activos de dengue.
- **Limpiar sales.** Si los SMILES vienen de PubChem o de un proveedor, pasarlos por
  `_limpiar_smiles_sal`: un contraión separado por `.` produce un grafo desconectado y
  una similitud sin sentido.
- **Revisar si el ErG los distingue.** Antes de calibrar, correr
  `comparar_A9_A22.py --v2` adaptado al par de referencias nuevas. Si dos de ellas dan
  similitud 1,000, significa que el ErG las ve como la misma molécula y aportan la misma
  información: hay que saberlo antes de reportar resultados "con dos referencias
  independientes". Fue exactamente el caso de A9 y A22.

---

## Paso 2 — Recalibrar los pesos θ

**No editar `DEFAULT_THETA` a mano.** Los seis pesos son el criterio químico del método y
se obtienen de datos.

### Si existe un cribado etiquetado para la diana nueva (el caso bueno)

Es la situación de influenza, gracias a PubChem AID 2326. Se necesita un archivo con
tres columnas: identificador, SMILES y desenlace (`Active` / `Inactive`).

1. Reemplazar en `recalibrar_v2.py` el bloque de carga de datos:
   ```python
   df = pd.read_excel("MI_CRIBADO.xlsx",
                      usecols=["CID", "SMILES", "OUTCOME"])
   df = df[df["OUTCOME"].isin(["Active", "Inactive"])]
   ```
2. Correr `python recalibrar_v2.py --rapido` primero, para verificar que todo funciona y
   estimar el tiempo.
3. Correr `python recalibrar_v2.py` completo. Produce el nuevo `theta_*.json` y el JSON
   con los cinco experimentos.
4. **Leer los AUC de prueba antes de seguir.** Si el AUC de prueba con activos mezclados
   sale cerca de 0,5 pero dentro de quimiotipos sale alto, el diagnóstico es el mismo que
   en esta tesis: la diana tiene varios mecanismos de acción y el método solo sirve para
   buscar análogos de cada inhibidor por separado.

### Si **no** existe cribado etiquetado (el caso del dengue)

Es la situación de la Fase B y hay que decirlo con claridad en cualquier informe: **sin
conjunto etiquetado no se puede medir generalización**. Las opciones, de mejor a peor:

1. **Transferir θ de una diana homóloga** calibrada con datos reales. Es lo que propone
   este proyecto: usar `theta_final_v2.json` (calibrado sobre NS1 de influenza) como
   punto de partida para NS1 de dengue. Supone que las dos proteínas premian el mismo
   tipo de parecido farmacofórico — una hipótesis razonable entre homólogas, pero una
   hipótesis.
2. **Construir decoys presuntos** con `construir_decoys_reales` (criterio tipo DUD‑E:
   peso molecular parecido, Tanimoto bajo) y calibrar contra ellos. El AUC resultante es
   optimista, porque los "inactivos" nunca se midieron.
3. **Usar `DEFAULT_THETA` sin calibrar** y reportarlo así. Es la opción más honesta
   cuando no hay nada: el ranking sigue siendo útil como ordenamiento por parecido
   farmacofórico, pero no hay forma de afirmar nada sobre su poder predictivo.

En los tres casos: **reportar explícitamente de dónde salió θ** y no presentar el AUC de
una diana como si fuera el de la otra.

---

## Paso 3 — Adaptar los criterios geométricos del bolsillo

Los cinco criterios de `ns1_utils.criterios_ns1` son una descripción del **bolsillo
CPSF30 de NS1**: canal alargado, hidrofóbico, con anillos aromáticos en los extremos.
Para otra diana hay que rehacer esa descripción.

### Ajuste rápido: solo los umbrales

```python
# ns1_utils.py
PARAMS = {
    "tam_min": 15.0,            # largo mínimo en Å → depende del tamaño del bolsillo
    "ratio_extremos": 0.70,     # qué tan "en las puntas" van los aromáticos
    "A_por_enlace": (1.5, 2.0), # rango de Å por enlace (cadena extendida)
}
```
O sin tocar código, desde la línea de comandos:
```bash
python evaluar_erg_ns1.py --tam-min 12 --ratio-extremos 0.5 --A-por-enlace 1.4 2.2
```

### Ajuste de fondo: otras condiciones

Si el bolsillo tiene otra forma, hay que **reescribir** `criterios_ns1`. Ejemplos de
cómo se traduce una descripción estructural a código:

| Descripción del bolsillo | Criterio equivalente |
|---|---|
| cavidad pequeña y compacta | `largo_max(coords) < tam_max` (invertir el criterio 4) |
| exige un catión que interactúe con un aspartato | contar nodos de un tipo nuevo `Cation` (ver paso 4) |
| dos subsitios polares separados ~8 Å | distancia 3D entre el HBD y el HBA más alejados dentro de un rango |
| bolsillo plano (hendidura superficial) | exigir pocos enlaces rotables, o varios aromáticos coplanares |
| sitio alargado con tres puntos de anclaje | exigir ≥3 nodos y que el grafo sea un camino, no una estrella |

### De dónde se saca esa descripción

- La estructura cristalográfica del complejo proteína–inhibidor (PDB), midiendo
  directamente el bolsillo.
- El artículo que caracterizó el sitio de unión (en este proyecto, Kleinpeter et al.
  2018 para CPSF30).
- Herramientas de detección de cavidades (fpocket, SiteMap) si solo se tiene la
  estructura de la proteína sola.

### Y una advertencia del propio proyecto

En esta tesis el filtro geométrico **no mejoró el ranking**: por sí solo dio AUC
0,30–0,63 y combinarlo con ErG+GED lo empeoró. La razón es que los criterios describen un
bolsillo concreto mientras que los activos del cribado actúan por mecanismos variados. La
lección extrapolable es: **medir siempre el aporte del filtro de forma aislada antes de
incorporarlo**, con `--solo-referencias` (¿lo cumplen los propios inhibidores conocidos?)
y con la columna `Solo filtro NS1` de la evaluación.

---

## Paso 4 — Añadir tipos de nodo (solo si hace falta)

Si la diana exige una propiedad que las cuatro etiquetas actuales no capturan —lo más
común, cargas—, hay que tocar cuatro sitios **de forma coherente**:

```python
# 1) pipeline_final_completo.py
LABEL_MAP = {"Donor": "HBD", "Acceptor": "HBA",
             "Aromatic": "Aromatico", "Hydrophobe": "Hidrofobo",
             "PosIonizable": "Cation", "NegIonizable": "Anion"}   # nuevo

PRIORITY = ["Cation", "Anion", "HBD", "HBA", "Aromatico", "Hidrofobo"]
TYPES    = ["Aromatico", "Hidrofobo", "HBA", "HBD", "Cation", "Anion"]

# 2) cost_matrix_from_theta: ampliar a 6x6 y añadir los pesos nuevos a theta
# 3) ged_mejorado.py: TYPES y _cost_matrix, con el MISMO orden
# 4) recalibrar_v2.py: ampliar `bounds` con un rango por cada peso nuevo
```

Consecuencias que hay que asumir: θ gana dimensiones, así que la calibración necesita
más iteraciones y más datos; y **todos los resultados anteriores dejan de ser
comparables**, porque los grafos cambian.

Pensar bien la prioridad: poner `Cation` antes que `HBD` significa que un amonio
protonado dejará de contarse como donador. Esa decisión cambia la topología de muchos
grafos.

---

## Paso 5 — Rango de conexión y umbral de quimiotipo

```python
# pipeline_final_completo.py, build_erg
if 1 <= d <= 6:          # ← rango de conexión, en enlaces
```
Ampliarlo (por ejemplo a 8) para bolsillos grandes donde grupos muy separados sí
interactúan a la vez con la proteína; reducirlo (a 4) para sitios compactos. **Cambia el
número de aristas y por tanto la escala del costo: obliga a recalibrar θ.**

```python
clusterizar_activos(activos, cutoff_distancia=0.6)   # 0,6 ⇔ Tanimoto > 0,4
```
Define qué cuenta como "mismo quimiotipo". Bajarlo da grupos más estrictos y más
pequeños; subirlo, grupos más laxos. Tiene efecto directo sobre el AUC que se reporta, así
que conviene **fijarlo, declararlo y no ajustarlo mirando el resultado**.

---

## Caso concreto: pasar de NS1 de influenza a NS1 de dengue (Fase B)

Lo que ya está hecho y lo que falta:

```python
# 1) Las referencias ya están en el código (Ansari et al. 2025):
DENGUE_ACTIVOS = {
    "Berberina":        "...",   # alcaloide de isoquinolina
    "Acido_betulinico": "...",   # triterpeno pentacíclico
    "Acido_ursolico":   "...",   # triterpeno pentacíclico
}
```

2. **Un ranking por cada uno, por separado.** Son tres quimiotipos distintos (un
   alcaloide y dos triterpenos): promediarlos reproduciría exactamente el error que la
   Fase A documentó. Esperar tres listas de candidatos muy diferentes entre sí, y
   presentarlas como tres listas.

3. **θ transferido de influenza**, declarándolo como supuesto. No hay cribado masivo
   etiquetado para NS1 de dengue, así que no se puede recalibrar ni medir
   generalización.

4. **Criterios geométricos**: los de CPSF30 **no aplican**. El bolsillo de dengue donde
   actúan berberina y los triterpenos no está caracterizado con el mismo detalle. Lo
   razonable es correr la Fase B **sin filtro geométrico** y decirlo, o derivar criterios
   nuevos de la estructura de NS1 de dengue si se dispone de ella.

5. **Comando:**
   ```bash
   # en correr_ranking.py: importar DENGUE_ACTIVOS en lugar de INFLUENZA_ACTIVOS
   python correr_ranking.py --theta-nuevo
   ```

6. **Cómo hay que presentar el resultado.** Un ranking priorizado de "primos químicos"
   de cada inhibidor conocido de dengue, útil como **prefiltro previo al docking o al
   ensayo experimental**, no como predicción de actividad. Esa es la afirmación que la
   evidencia de la Fase A permite sostener, y no más.

---

## Lista de verificación antes de reportar resultados en una diana nueva

- [ ] Los SMILES de referencia están verificados y con fuente citada en el código.
- [ ] Se comprobó si el ErG distingue entre sí a las referencias (`comparar_*.py`).
- [ ] Se declaró de dónde salió θ: calibrado con datos de esta diana, transferido de
      otra, o sin calibrar.
- [ ] Todo AUC reportado es de un conjunto de **prueba** que no se usó para calibrar.
- [ ] Hay una **línea base** (Morgan+Tanimoto) sobre exactamente las mismas moléculas de
      prueba.
- [ ] Los AUC van con **intervalo de confianza** (`auc_ci`), sobre todo con pocos activos.
- [ ] Los rankings se reportan **por referencia**, nunca promediados.
- [ ] Se contaron y registraron las moléculas descartadas por fallo de embedding 3D.
- [ ] Se revisó si hay **empates** en el borde del top‑N (`dibujar_top20.py` lo imprime).
- [ ] Se declaró que el resultado es un prefiltro, no una predicción de actividad.

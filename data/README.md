# Datos de entrada

Los archivos de datos **no se versionan** en este repositorio: suman más de 60 MB y
todos son de acceso público, así que se pueden reconstruir. Esta carpeta explica cómo
obtener cada uno y dónde debe quedar.

Por defecto los scripts buscan los archivos **en la raíz del repositorio** (junto a
`pipeline_final_completo.py`), no dentro de `data/`. Si se prefiere tenerlos aquí, hay
que ajustar las rutas `EXCEL` en `correr_ranking.py` y la ruta de `BasePubChem.xlsx` en
`recalibrar_v2.py`.

---

## 1. `Moleculas.xlsx` — biblioteca de candidatos (≈ 10 MB)

La biblioteca que se criba: **33 047 compuestos** descargados de ChEMBL.

Columnas que usa el código: `ChEMBL ID`, `Name`, `Smiles`. El archivo trae además 40
columnas de descriptores calculados por ChEMBL (peso molecular, AlogP, HBA, HBD, enlaces
rotables, anillos aromáticos, fase clínica, etc.) que se usan para caracterizar el
conjunto y para elegir los inactivos "difíciles".

**Cómo obtener un archivo equivalente:**

1. Entrar a <https://www.ebi.ac.uk/chembl/>.
2. Buscar o filtrar el conjunto de compuestos de interés.
3. Botón de descarga → formato **XLSX** (o CSV) con todas las columnas.

Cualquier archivo con una columna de SMILES y una de identificador funciona: basta
ajustar `smiles_col` e `id_col` en las llamadas a `rankear_excel_en_bloques` o
`cargar_candidatos_excel`.

---

## 2. `BasePubChem.xlsx` / `AID2326_cache.csv` — conjunto etiquetado (≈ 32 MB / 18 MB)

El *ground truth* de toda la Fase A: el cribado de alto rendimiento **PubChem BioAssay
AID 2326** del NIH Chemical Genomics Center contra la función de NS1 de influenza — el
mismo cribado del que salieron A9 y A22.

| Desenlace | Compuestos |
|---|---|
| Active | 1 073 |
| Inactive | 260 684 |
| Inconclusive | 8 410 |
| Unspecified | 130 |
| **Total** | **270 315** |

**Cómo obtenerlo:**

1. Entrar a <https://pubchem.ncbi.nlm.nih.gov/bioassay/2326>.
2. Botón **Download**.
3. En el modal, sección **"Deposited Data Table Only"** → **CSV** → guardar.
   (Las otras opciones descargan metadatos del ensayo, no los resultados.)
4. Ese CSV ya trae los SMILES en la columna `PUBCHEM_EXT_DATASOURCE_SMILES`; no hace
   falta pedirlos aparte vía PUG REST.

Columnas que usa el código: `PUBCHEM_CID`, `PUBCHEM_EXT_DATASOURCE_SMILES`,
`PUBCHEM_ACTIVITY_OUTCOME`.

**Sobre el caché:** la primera corrida de `recalibrar_v2.py` lee el Excel completo
(270 000 filas, varios minutos), se queda solo con las filas `Active`/`Inactive` y
guarda `AID2326_cache.csv` para las siguientes. Si solo se tiene el CSV original de
PubChem, se puede usar directamente con `cargar_pubchem_ns1_assay`.

⚠️ **Cuidado con el encabezado.** El CSV que entrega PubChem trae cuatro filas de
metadatos debajo del encabezado, y `cargar_pubchem_ns1_assay` las salta con
`skiprows=[1,2,3,4]`. El `AID2326_cache.csv` generado desde el Excel **no** las tiene, y
ese salto borraría cuatro moléculas reales. `pipeline_v3.py` filtra por contenido en vez
de por posición, que es lo correcto.

---

## 3. Archivos intermedios de los scripts históricos (no necesarios)

`calibracion_cluster.py` y `calibracion_pubchem_real.py` dependen de tres archivos
intermedios que ya no se generan:

- `activos_clusters.pkl`
- `AID2326_activos_reales.csv`
- `AID2326_inactivos_muestra5000.csv`

Esos scripts se conservan solo por trazabilidad y **no son reproducibles tal cual**. Para
calibrar hoy se usa `recalibrar_v2.py`, que parte directamente de `BasePubChem.xlsx`.

---

## 4. `merged_500_con_props.csv` (opcional)

Muestra de 500 moléculas con propiedades, usada por el bloque de demostración
(`__main__`) de `pipeline_final_completo.py` para construir decoys presuntos. No es
necesaria para reproducir los resultados: se puede generar tomando 500 filas de
`Moleculas.xlsx` con las columnas `ChEMBL_ID`, `Smiles` y `Molecular Weight`.

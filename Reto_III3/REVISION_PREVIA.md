# Revisión previa de la entrega — Reto III.3

Esta revisión es editorial y de consistencia interna. No sustituye una auditoría numérica
completa de los scripts ni certifica todos los resultados científicos.

## Cambios organizativos realizados

- Se renombraron las carpetas principales con numeración consecutiva y nombres descriptivos.
- Se corrigió la secuencia que saltaba del número 2 al 4.
- Se expandió el ZIP del repositorio de código en `02_Codigo_y_Datos/` para que los archivos
  fuente, pruebas, datos, figuras y documentación puedan navegarse en GitHub.
- Se conservó el material de informe, bitácora, presentación, videos y apoyo.
- En la diapositiva de portada se retiró la identificación de la universidad y la escuela;
  se conservó la fecha.

## Observaciones que deben verificarse antes de la entrega

### 1. Diferencia en los valores de Voyager 2

La diapositiva 2 informa `ΔE = 164.2 km²/s²`, `n medio = 173.2 km²/s²`
(5.5 %) y con velocidad angular instantánea `164.8 km²/s²` (0.4 %).
El README interno del repositorio de código presenta esos mismos valores.

En cambio, el informe PDF, en la sección de formulación/verificación, presenta
`ΔE = 172.42 km²/s²` y `n Δℓz = 173.38 km²/s²` (0.6 %).
Son cifras incompatibles entre sí. Hay que identificar qué definición, intervalo temporal,
estado inicial/final y archivo de resultados corresponde a la cifra correcta y unificar
informe, diapositivas, bitácora y README. No se modificó ninguna cifra para evitar
elegir arbitrariamente.

### 2. Límite de C3

El informe utiliza `C3 ≤ 100 km²/s²` para describir el óptimo, mientras que el pie de la
figura 1 menciona `C3 ≤ 105 km²/s²`. La diapositiva 4 usa `C3 ≤ 100 km²/s²`.
Conviene confirmar el umbral usado en la búsqueda y corregir el texto o el pie de figura.

### 3. Bitácora de IA

La bitácora declara que los estudiantes revisan cada entrada, pero en la tabla las celdas
«Revisado por» aparecen sin nombres. Además, algunos prompts se presentan como resúmenes.
Antes de entregar, deben completar los revisores reales y distinguir transcripciones literales
de resúmenes; no se deben inventar revisiones ni prompts.

### 4. Presentación

La presentación tiene seis diapositivas y tres imágenes insertadas. No se detectaron
imágenes idénticas por hash entre esos tres recursos. En esta reorganización no se rediseñó
la presentación ni se sustituyeron esas imágenes; únicamente se retiró la identificación
institucional de la portada y se normalizó el nombre del archivo.

## Recomendación

Antes de publicar, resolver como mínimo los puntos 1 y 2, porque son discrepancias numéricas
visibles entre materiales que deberían contar la misma historia. Después, abrir el ZIP
finalmente preparado y revisar que todos los archivos se visualicen correctamente en GitHub.

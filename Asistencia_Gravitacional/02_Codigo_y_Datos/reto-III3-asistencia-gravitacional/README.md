# Reto III.3: Diseñar una asistencia gravitacional

Mecánica Clásica 2026-B, UIS. Encuentro con Júpiter que lleva una sonda hacia Saturno,
resuelto con **cónicas empalmadas** (Modelo 1) e **integración de N cuerpos** (Modelo 2),
con **efemérides reales** (JPL DE421) y validación con la **trayectoria real de Voyager 2**.

**Autores:** Frankly Giovanni Poveda Pinzón (código 2221414) y Camila Castillo López (código 2221748). Universidad Industrial de Santander, Escuela de Física.

## Animaciones

`docs/animaciones/index.html` (se abre con doble clic; necesita `datos.js` en la misma carpeta) reúne cuatro escenas interactivas: el encuentro visto desde Júpiter y desde el Sol, la conservación de h, el error de las cónicas y Voyager 2 real. Los datos de las escenas 2 a 4 se regeneran con `python scripts/07_exportar_animaciones.py`.

### Videos de referencia

`docs/animaciones/video/` contiene la escena (`escena.html`, dibujada en canvas con los datos de `datos.js` y `mision.js`), los tiempos de cada clip (`cfg.py`, `guion.json`) y el script que la graba (`render.py`, requiere Playwright y ffmpeg). Son cinco clips de 25 a 29 s, sin voz ni subtítulos (misión, encuentro, conservación de h, error de las cónicas, Voyager 2) y una versión completa de unos 2 min. El enlace está en la introducción del informe. Los datos de `mision.js` salen de `scripts/08_datos_mision_video.py`.

## Reproducir todo

Requiere Python ≥ 3.12 para las versiones exactas de `requirements.txt` (con las que se
generaron los resultados). Con Python 3.10 u 3.11 use `pip install -r requirements-min.txt`.

```bash
pip install -r requirements.txt
python reproducir.py             # 40 pruebas + resultados + figuras (≈ 2 min)
python reproducir.py --reporte   # además, figuras para el informe en figuras/reporte/
```

No hay números aleatorios: los resultados son deterministas. Las tolerancias numéricas están
declaradas en `reproducir.py`.

## Estructura

| Carpeta | Contenido |
|---------|-----------|
| `src/asistencia/` | `kepler.py` (dos cuerpos), `lambert.py`, `sobrevuelo.py` (hipérbola 3D), `efemerides.py` (DE421 y circular), `voyager.py` (kernel SPICE), `modelo1.py`, `modelo2.py` (REBOUND), `jupiter_real.py` (Modelo 2 ampliado: satélites galileanos y J2) |
| `tests/` | `test_modelo1.py`, `test_modelo2.py` (Jacobi, convergencia, límites, maniobra de Voyager 2), `test_datos.py` (datos reales frente a valores publicados) |
| `scripts/` | `01_buscar_ventanas.py`, `02_figuras_modelo1.py`, `03_voyager2.py`, `04_verificar_lagrangiano.py`, `05_modelo2.py`, `06_voyager2_satelites.py`, `07_exportar_animaciones.py`, `descargar_horizons.py` (opcional) |
| `datos/` | Datos crudos y su procedencia (`LEEME.md`) |
| `resultados/` | Salidas numéricas en JSON y NPZ |
| `figuras/` | Figuras (y `figuras/reporte/` para LaTeX) |
| `docs/` | Informe, deducción lagrangiana, hipótesis, bitácora de IA (4D) y entregas previas |

## Resultados

**Modelo 1 con efemérides reales.** Sin usar las fechas históricas, la búsqueda encuentra la
ventana de 1977. Óptimo con C3 ≤ 100 km²/s²: salida 3 sep 1977, Júpiter 12 abr 1979,
Saturno 1 feb 1981 (3.42 años).

**Validación con datos no ajustados.** Evaluado en las fechas reales, el modelo predice:

| | r_p predicho (circular) | r_p predicho (DE421) | r_p real |
|---|---|---|---|
| Voyager 2 | 12.6 R_J | 9.72 R_J | 10.09 R_J |
| Voyager 1 | 5.75 R_J | 4.74 R_J | 4.88 R_J |

La hipótesis de órbitas circulares (H1) introducía un error del 25 % en r_p; con
efemérides reales baja a 3–4 %.

**Qué se conserva** (`docs/formulacion_lagrangiana_reto_III3.pdf`,
`scripts/04_verificar_lagrangiano.py`): en el sistema rotante Sol–Júpiter se conserva
h = E − n ℓ_z (constante de Jacobi), exacta si la órbita de Júpiter es circular. Con los datos
reales de Voyager 2 (`resultados/voyager2.json`), ΔE = 164.2 km²/s²; con el n medio de Júpiter
nΔℓ_z = 173.2 (falla 5.5 %) y con su velocidad angular instantánea ∫ω dℓ_z = 164.8 (0.4 %),
porque la órbita real de Júpiter es excéntrica. Por lo mismo, Tisserand cambia 3.6 %.

**Modelo 2** (`scripts/05_modelo2.py`, REBOUND/IAS15):
- Jacobi constante a 2.5×10⁻¹⁴ mientras la energía heliocéntrica cambia; convergencia a ≤ 2.1 m
  en 250 días; DOP853 (independiente) difiere 3 m.
- **Error de las cónicas** en la ganancia de energía (R_e = r_SOI, r_p = 10 R_J): 17.5 % con
  v∞ = 2.5 km/s, 5.7 % con 5, 2.0 % con 8, 0.13 % con 15. Laplace y Hill son el mejor
  compromiso; radios menores o mayores empeoran el resultado.

**Voyager 2 real** (`scripts/03_voyager2.py`, `scripts/06_voyager2_satelites.py`):
- El kernel reproduce lo publicado: máximo acercamiento al centro de Júpiter a las 22:29:01 UTC
  del 9 jul 1979, a 721 673 km (publicado: 22:29 UTC, 721 670 km); sobrevuelos de los cuatro
  satélites galileanos al 0.2 % o mejor.
- Desde el estado real del 4 jun 1979, el Modelo 2 con Júpiter puntual sigue la trayectoria con
  1.3×10³ km de error en el periapsis (la hipérbola: 1.0×10⁵ km), pero se separa 6.8×10⁴ km
  después y no reproduce la caída de la energía jovicéntrica (v∞ equivalente 7.650 → 7.578 km/s).
- La causa: los satélites galileanos y J2 (−0.029 km/s) y una **maniobra** de 76.6 minutos que
  empieza 130 minutos después del máximo acercamiento (2.83 mm/s², Δv = 13.0 m/s, −0.033 km/s),
  documentada en NASA SP-439. Con todo incluido: v∞ final 7.5780 frente a 7.5777 km/s, y 6 km
  de error en los 11 días posteriores a la maniobra.
- Lejos de Júpiter la energía heliocéntrica pasa de −121.8 a +26.3 km²/s²: la sonda queda en
  escape (e = 1.30).

## Documentos

- `docs/informe/informe_reto_III3.pdf`: informe técnico (8 páginas, modelito UIS). Sus figuras
  son las de `figuras/reporte/` (`python reproducir.py --reporte`).
- `docs/formulacion_lagrangiana_reto_III3.pdf`: deducción completa de la formulación.
- `docs/hipotesis.md` y `docs/bitacora_ia.md`.
- `docs/banco_preguntas_defensa.pdf`: preguntas probables de la sustentación.
- `docs/entregas_previas/`: avance de identificación del problema y marco teórico (octubre).

## Pendiente

1. (Opcional) comparar con los vectores de JPL Horizons (`scripts/descargar_horizons.py`).
2. Repetir el barrido del error de las cónicas con la órbita excéntrica de Júpiter y en 3D.

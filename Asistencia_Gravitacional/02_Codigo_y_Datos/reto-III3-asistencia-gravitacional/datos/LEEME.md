# Datos crudos

## `vgr2_jup230.bsp` (trayectoria real de Voyager 2)

- **Qué es:** kernel SPICE (formato SPK, tipo 1) con la trayectoria reconstruida de
  Voyager 2 respecto al baricentro de Júpiter, del 4 de junio al 20 de julio de 1979.
- **Origen:** NAIF/JPL, carpeta `pub/naif/VOYAGER/kernels/spk/` del servidor de NAIF.
- **Cómo se obtuvo:** el servidor de NAIF no era accesible desde el entorno de trabajo,
  así que se descargó la copia idéntica que distribuye el paquete `asteRisk`
  (github.com/Rafael-Ayala/asteRisk, carpeta `inst/`). Tamaño: 330 752 bytes.
- **Cómo se comprobó que es el archivo correcto:** con valores publicados por fuentes
  independientes (`tests/test_datos.py`, `resultados/voyager2_satelites.json`):

  | Máximo acercamiento a | Este trabajo (kernel) | Publicado |
  |---|---|---|
  | Júpiter (centro del planeta) | 9 jul 1979, 22:29:01 UTC, 721 673 km | 22:29 UTC, 721 670 km |
  | Calisto | 8 jul, 12:21:33 UTC, 215 281 km | 12:21 UTC, 214 930 km |
  | Ganímedes | 9 jul, 07:14:44 UTC, 62 078 km | 07:14 UTC, 62 130 km |
  | Europa | 9 jul, 17:51:25 UTC, 205 641 km | 17:53 UTC, 205 720 km |
  | Ío | 9 jul, 23:17:45 UTC, 1 129 728 km | 23:17 UTC, 1 129 900 km |

  El kernel da la posición respecto al baricentro del sistema joviano (721 552 km a las
  22:28:52 UTC); el centro del planeta está a menos de 200 km de ese punto y se obtiene con
  las posiciones de los satélites. Valores publicados: tabla de la misión (NSSDC/JPL),
  reproducida por el PDS Ring-Moon Systems Node.
- **El kernel contiene una maniobra.** Entre las 00:38:44 y las 01:55:23 UTC del 10 de julio
  (de 130 a 206 minutos después del máximo acercamiento) la sonda tiene una aceleración
  constante de 2.83 mm/s² que la gravedad no explica: Δv = 13.0 m/s en 76.6 minutos. Es la
  maniobra de corrección de trayectoria que describe NASA SP-439 (Morrison y Samz, *Voyage to
  Jupiter*, 1980, cap. 7): un encendido de 76 minutos unas dos horas después del máximo
  acercamiento. Ver `scripts/06_voyager2_satelites.py`.
- **Pendiente:** comparar con JPL Horizons (`scripts/descargar_horizons.py`, ID −32).
- **Uso:** solo validación. No se ajusta ningún parámetro con estos datos.

## Satélites galileanos

- **Posiciones:** función `JupiterMoons` del paquete de PyPI `astronomy-engine==2.1.19`
  (D. Cross, licencia MIT), que da los vectores de Ío, Europa, Ganímedes y Calisto respecto al
  centro de Júpiter en el ecuador J2000. No requiere internet.
- **Constantes:** GM de los satélites (5959.916, 3202.739, 9887.834 y 7179.289 km³/s²),
  J2 = 14696.43×10⁻⁶ con radio de referencia 71 492 km (solución JUP230 de JPL, R. A. Jacobson)
  y polo de Júpiter de la IAU (α = 268.057°, δ = 64.495°).
- **Comprobación:** la tabla de arriba. Las distancias de sobrevuelo coinciden con las
  publicadas al 0.2 % o mejor, lo que valida a la vez el kernel, el modelo de los satélites y
  los marcos de referencia.

## Efeméride planetaria JPL DE421

- **Qué es:** efeméride numérica del JPL (Folkner, Williams y Boggs 2009), 1900–2050.
- **Cómo se obtiene:** paquete de PyPI `de421==2008.1` (está en `requirements.txt`), leído con
  `jplephem`. No se copia aquí porque ocupa 27 MB y se instala con pip.
- **Marco:** DE421 da posiciones baricéntricas ecuatoriales (ICRF). El código las pasa a
  heliocéntricas y a la eclíptica J2000 (oblicuidad 84381.448″).
- **Escala de tiempo:** DE421 y el kernel usan TDB. Las fechas UTC se convierten con
  TDB = UTC + 32.184 s + (TAI − UTC), con la tabla de segundos intercalares del IERS
  (50.184 s en julio de 1979).
- **Comprobación:** `tests/test_datos.py` (Tierra en J2000, distancias e inclinaciones de
  Júpiter y Saturno dentro de sus valores publicados).

## Fechas históricas de las Voyager

PDS Ring-Moon Systems Node, https://pds-rings.seti.org/voyager/mission/ (lanzamientos y máximos
acercamientos a Júpiter y Saturno). Maniobra de Voyager 2 en Júpiter: D. Morrison y J. Samz,
*Voyage to Jupiter*, NASA SP-439 (1980), capítulo 7.

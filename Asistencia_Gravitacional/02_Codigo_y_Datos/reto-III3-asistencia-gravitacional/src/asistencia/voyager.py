"""Trayectoria real de Voyager 2 durante el encuentro con Júpiter.

Fuente: kernel SPICE `vgr2_jup230.bsp` de NAIF/JPL (trayectoria reconstruida de
Voyager 2 respecto al baricentro de Júpiter, 4 jun – 21 jul 1979). Ver datos/LEEME.md.
Estos datos NO se usan para ajustar nada: solo para validar.
"""
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from .constantes import DIA
from .efemerides import estado_planeta

KERNEL = Path(__file__).resolve().parents[2] / "datos" / "vgr2_jup230.bsp"
_CARGADO = False
# El segmento es SPK tipo 1 (diferencias modificadas), que jplephem no lee: se usa
# CSPICE a través de spiceypy. El marco ECLIPJ2000 viene incluido en CSPICE.


def _cargar():
    global _CARGADO
    if not _CARGADO:
        import spiceypy as sp
        sp.furnsh(str(KERNEL))
        _CARGADO = True


def _et(jd):
    return (jd - 2451545.0) * DIA          # segundos TDB desde J2000


def cobertura():
    """(jd_inicio, jd_fin) en TDB del segmento de Voyager 2."""
    import spiceypy as sp
    _cargar()
    ventana = sp.spkcov(str(KERNEL), -32)
    a, b = sp.wnfetd(ventana, 0)
    return 2451545.0 + a / DIA, 2451545.0 + b / DIA


def estado_relativo_jupiter(jd):
    """Posición (km) y velocidad (km/s) de Voyager 2 respecto al baricentro de Júpiter,
    eclíptica J2000. jd en TDB."""
    import spiceypy as sp
    _cargar()
    x, _ = sp.spkez(-32, _et(jd), "ECLIPJ2000", "NONE", 5)
    return np.array(x[:3]), np.array(x[3:])


def estado_heliocentrico(jd):
    """Estado heliocéntrico = Júpiter (DE421) + Voyager 2 relativo a Júpiter (SPICE)."""
    rJ, vJ = estado_planeta("Jupiter", jd, "de421")
    r, v = estado_relativo_jupiter(jd)
    return rJ + r, vJ + v


def instante_de_minimo(f, jd_a, jd_b, n=4000):
    """Instante (JD) y valor del mínimo de f(jd) en [jd_a, jd_b]: malla de n puntos y
    refinado con Brent.

    El refinado se hace en segundos respecto al mejor punto de la malla y no en días
    julianos: la tolerancia del método acotado de scipy es relativa al valor de la variable
    (1.5e-8 |x|) y con x ~ 2.4e6 días equivale a 50 minutos. La primera versión refinaba
    en días julianos y daba el máximo acercamiento a Júpiter con 3 minutos de error
    (bitácora, entrada 25)."""
    ts = np.linspace(jd_a, jd_b, n)
    t0 = ts[int(np.argmin([f(t) for t in ts]))]
    paso = (jd_b - jd_a) / (n - 1) * DIA
    r = minimize_scalar(lambda s: f(t0 + s / DIA), bounds=(-2 * paso, 2 * paso), method="bounded",
                        options={"xatol": 1e-3})
    return t0 + r.x / DIA, float(r.fun)


def maximo_acercamiento():
    """(jd, distancia en km) del máximo acercamiento de Voyager 2 al baricentro del sistema
    de Júpiter."""
    a, b = cobertura()
    return instante_de_minimo(lambda t: np.linalg.norm(estado_relativo_jupiter(t)[0]), a + 0.01, b - 0.01, 8000)

"""Constantes físicas y orbitales usadas en todo el proyecto.

Unidades internas: km, s, km/s. Las fechas se manejan como días julianos en TDB.

Fuentes:
- Parámetros gravitacionales: constantes de la efeméride JPL DE421 (Folkner et al. 2009),
  que es la misma efeméride que usamos para las posiciones de los planetas.
- Radio de Júpiter: IAU (radio ecuatorial a 1 bar).
- Elementos medios de J2000 para el modelo circular: Standish, "Keplerian Elements for
  Approximate Positions of the Major Planets", JPL (tabla 1, 1800–2050).
"""
import numpy as np

AU = 1.495978707e8            # km (valor IAU 2012; DE421 usa 149597870.691 km, diferencia 1e-9)
DIA = 86400.0                 # s
JD_J2000 = 2451545.0

# DE421: GMS = 2.959122082855911e-4 UA^3/día^2 ; GM5 (sistema joviano) = 2.825345909524226e-7
_UA3_DIA2 = 149597870.691**3 / DIA**2
GM_SOL = 2.959122082855911e-4 * _UA3_DIA2        # ≈ 1.32712440e11 km^3/s^2
GM_JUPITER = 2.825345909524226e-7 * _UA3_DIA2    # ≈ 1.26712765e8 km^3/s^2 (Júpiter + satélites)
R_JUPITER = 71492.0           # km
OBLICUIDAD_J2000 = np.radians(84381.448 / 3600)  # rad, para pasar de ecuatorial a eclíptica

# Semiejes mayores (UA) y longitudes medias en J2000 (grados), Standish tabla 1.
PLANETAS = {
    "Tierra":  {"a_ua": 1.00000261, "L0_deg": 100.46457166},
    "Jupiter": {"a_ua": 5.20288700, "L0_deg": 34.39644051},
    "Saturno": {"a_ua": 9.53667594, "L0_deg": 49.95424423},
}


def radio_esfera_influencia(a_km, gm_planeta, gm_sol=GM_SOL):
    """Esfera de influencia de Laplace: r = a (m_p / M_sol)^(2/5)."""
    return a_km * (gm_planeta / gm_sol) ** 0.4


def radio_hill(a_km, gm_planeta, gm_sol=GM_SOL):
    """Radio de Hill: r = a (m_p / 3 M_sol)^(1/3)."""
    return a_km * (gm_planeta / (3 * gm_sol)) ** (1 / 3)


A_JUPITER_KM = PLANETAS["Jupiter"]["a_ua"] * AU
R_SOI_JUPITER = radio_esfera_influencia(A_JUPITER_KM, GM_JUPITER)
R_HILL_JUPITER = radio_hill(A_JUPITER_KM, GM_JUPITER)
V_CIRC_JUPITER = np.sqrt(GM_SOL / A_JUPITER_KM)

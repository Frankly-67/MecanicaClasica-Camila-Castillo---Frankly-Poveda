"""Posiciones y velocidades heliocéntricas de los planetas, en el marco eclíptico J2000.

Dos fuentes, para poder medir el efecto de la hipótesis H1 (ver docs/hipotesis.md):

- "de421"   : efeméride real JPL DE421 (Folkner, Williams y Boggs 2009), distribuida en
              el paquete de Python `de421` y leída con `jplephem`. Es la fuente por omisión.
- "circular": órbitas circulares y coplanares con la longitud media de J2000 (Standish).
              Es el modelo simplificado de la primera versión.

Las fechas son días julianos TDB. Para Júpiter y Saturno DE421 da el baricentro del
sistema (planeta + satélites); la diferencia con el centro del planeta es de cientos de km.
"""
from datetime import datetime, timezone

import numpy as np

from .constantes import AU, DIA, GM_SOL, JD_J2000, OBLICUIDAD_J2000, PLANETAS

FUENTE = "de421"
_EFEM = None

# TDB - UTC = 32.184 s + (TAI - UTC). TAI - UTC cambia con cada segundo intercalar (tabla del
# IERS, Bulletin C). Se desprecia TDB - TT, que es periódica y menor que 2 ms.
_TAI_MENOS_UTC = [("1972-01-01", 10), ("1972-07-01", 11), ("1973-01-01", 12), ("1974-01-01", 13),
                  ("1975-01-01", 14), ("1976-01-01", 15), ("1977-01-01", 16), ("1978-01-01", 17),
                  ("1979-01-01", 18), ("1980-01-01", 19), ("1981-07-01", 20), ("1982-07-01", 21),
                  ("1983-07-01", 22), ("1985-07-01", 23), ("1988-01-01", 24), ("1990-01-01", 25),
                  ("1991-01-01", 26), ("1992-07-01", 27), ("1993-07-01", 28), ("1994-07-01", 29),
                  ("1996-01-01", 30), ("1997-07-01", 31), ("1999-01-01", 32), ("2006-01-01", 33),
                  ("2009-01-01", 34), ("2012-07-01", 35), ("2015-07-01", 36), ("2017-01-01", 37)]
_SALTOS = [(datetime.fromisoformat(f).replace(tzinfo=timezone.utc).timestamp(), n)
           for f, n in _TAI_MENOS_UTC]


def tdb_menos_utc(unix_utc):
    """Segundos que hay que sumar a UTC para tener TDB, en el instante dado (segundos Unix)."""
    n = 10
    for t, k in _SALTOS:
        if unix_utc >= t:
            n = k
    return 32.184 + n


def fecha_a_jd(fecha):
    """'YYYY-MM-DD[THH:MM[:SS]]' (UTC) -> día juliano TDB."""
    if isinstance(fecha, str):
        fecha = datetime.fromisoformat(fecha)
    if fecha.tzinfo is None:
        fecha = fecha.replace(tzinfo=timezone.utc)
    return (fecha.timestamp() + tdb_menos_utc(fecha.timestamp())) / DIA + 2440587.5


def jd_a_fecha(jd, hora=False, segundos=False):
    """Día juliano TDB -> fecha UTC 'YYYY-MM-DD', con hora y minutos (hora=True) o hasta
    segundos (segundos=True)."""
    tdb = (jd - 2440587.5) * DIA
    utc = tdb - tdb_menos_utc(tdb - 50.0)
    utc = tdb - tdb_menos_utc(utc)
    f = datetime.fromtimestamp(round(utc, 3), tz=timezone.utc)   # evita 14:28:59.99996 -> '14:28'
    if segundos:
        return f.strftime("%Y-%m-%d %H:%M:%S")
    return f.strftime("%Y-%m-%d %H:%M" if hora else "%Y-%m-%d")


_ROT = np.array([[1, 0, 0],
                 [0, np.cos(OBLICUIDAD_J2000), np.sin(OBLICUIDAD_J2000)],
                 [0, -np.sin(OBLICUIDAD_J2000), np.cos(OBLICUIDAD_J2000)]])


def ecuatorial_a_ecliptica(x):
    return _ROT @ x


def _efemeride():
    global _EFEM
    if _EFEM is None:
        import de421
        from jplephem import Ephemeris
        _EFEM = Ephemeris(de421)
    return _EFEM


_NOMBRE_JPL = {"Jupiter": "jupiter", "Saturno": "saturn", "Tierra": "earthmoon"}


def _baricentrico(nombre, jd):
    e = _efemeride()
    r, v = e.position_and_velocity(_NOMBRE_JPL[nombre], jd)
    r, v = r[:, 0], v[:, 0] / DIA
    if nombre == "Tierra":       # Tierra = baricentro Tierra-Luna - Luna/(1 + M_T/M_L)
        rm, vm = e.position_and_velocity("moon", jd)
        f = 1 / (1 + e.EMRAT)
        r, v = r - f * rm[:, 0], v - f * vm[:, 0] / DIA
    return r, v


def estado_planeta(nombre, jd, fuente=None):
    """Posición (km) y velocidad (km/s) heliocéntricas, eclíptica J2000."""
    fuente = fuente or FUENTE
    if fuente == "circular":
        p = PLANETAS[nombre]
        a = p["a_ua"] * AU
        n = np.sqrt(GM_SOL / a**3)
        L = np.radians(p["L0_deg"]) + n * (jd - JD_J2000) * DIA
        return (a * np.array([np.cos(L), np.sin(L), 0.0]),
                a * n * np.array([-np.sin(L), np.cos(L), 0.0]))
    if fuente == "de421":
        rp, vp = _baricentrico(nombre, jd)
        rs, vs = _efemeride().position_and_velocity("sun", jd)
        return (ecuatorial_a_ecliptica(rp - rs[:, 0]),
                ecuatorial_a_ecliptica(vp - vs[:, 0] / DIA))
    raise ValueError(f"fuente desconocida: {fuente}")

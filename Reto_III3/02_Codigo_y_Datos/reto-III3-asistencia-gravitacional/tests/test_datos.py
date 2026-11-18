"""Verificación de los datos reales: que se leen bien y en el marco correcto.

Cada prueba compara con un valor publicado por una fuente distinta de la que generó
los datos (no con otro resultado de nuestro propio código).
"""
import numpy as np
import pytest

from asistencia.constantes import AU, R_JUPITER
from asistencia.efemerides import estado_planeta, fecha_a_jd, jd_a_fecha
from asistencia.voyager import estado_relativo_jupiter, instante_de_minimo, maximo_acercamiento


def test_tierra_en_j2000():
    """El 1 ene 2000 12h TDB la Tierra está a 0.9833 UA del Sol (perihelio fue el 3 ene)
    y su longitud eclíptica heliocéntrica es la del Sol geocéntrico + 180° ≈ 100.4°."""
    r, _ = estado_planeta("Tierra", 2451545.0, "de421")
    assert abs(np.linalg.norm(r) / AU - 0.9833) < 0.0005
    assert abs(np.degrees(np.arctan2(r[1], r[0])) - 100.4) < 0.2


@pytest.mark.parametrize("nombre,q,Q,imax", [("Jupiter", 4.95, 5.46, 1.31), ("Saturno", 9.0, 10.1, 2.49)])
def test_planetas_dentro_de_su_perihelio_afelio_e_inclinacion(nombre, q, Q, imax):
    """Durante 1977–1982 la distancia debe quedar entre perihelio y afelio publicados,
    y la latitud eclíptica no puede superar la inclinación orbital."""
    for jd in np.linspace(fecha_a_jd("1977-01-01"), fecha_a_jd("1982-01-01"), 60):
        r, _ = estado_planeta(nombre, jd, "de421")
        d = np.linalg.norm(r) / AU
        assert q - 0.02 < d < Q + 0.02
        assert abs(np.degrees(np.arcsin(r[2] / np.linalg.norm(r)))) < imax + 0.01


def test_circular_y_de421_se_parecen():
    """Las dos fuentes de efemérides deben diferir solo por la excentricidad e inclinación
    (unos pocos por ciento), no por un error de marco o de fechas."""
    jd = fecha_a_jd("1979-07-09")
    r1, v1 = estado_planeta("Jupiter", jd, "circular")
    r2, v2 = estado_planeta("Jupiter", jd, "de421")
    assert np.degrees(np.arccos(r1 @ r2 / np.linalg.norm(r1) / np.linalg.norm(r2))) < 6
    assert abs(np.linalg.norm(v1) - np.linalg.norm(v2)) < 0.7


def test_escala_de_tiempo_utc_a_tdb():
    """La época J2000 es, por definición, el 1 ene 2000 a las 12:00:00 TT = 11:58:55.816 UTC
    (TT - UTC = 64.184 s ese día). Comprueba la tabla de segundos intercalares, y que en
    1979 la diferencia era de 50.184 s (TAI - UTC = 18 s)."""
    assert abs(fecha_a_jd("2000-01-01T11:58:55.816") - 2451545.0) < 1e-9
    assert abs((fecha_a_jd("1979-07-09T00:00") - 2444063.5) * 86400 - 50.184) < 1e-4
    assert jd_a_fecha(fecha_a_jd("1977-08-20T14:29"), hora=True) == "1977-08-20 14:29"


def test_voyager2_maximo_acercamiento_coincide_con_lo_publicado():
    """PDS Ring-Moon Systems Node: máximo acercamiento de Voyager 2 a Júpiter el
    9 jul 1979 a las 22:29 UTC, a 721 670 km del centro. El kernel da la posición respecto
    al baricentro del sistema joviano, que está a menos de 200 km del centro del planeta:
    se tolera 0.05 % en la distancia y un minuto en la hora."""
    t, d = maximo_acercamiento()
    assert abs(d - 721670) / 721670 < 5e-4
    assert abs((t - fecha_a_jd("1979-07-09T22:29")) * 1440) < 1.0
    assert 10.0 < d / R_JUPITER < 10.2


def test_voyager2_respecto_al_centro_de_jupiter():
    """Corrigiendo del baricentro al centro del planeta con las posiciones de los satélites,
    la distancia mínima debe acercarse aún más al valor publicado (721 670 km)."""
    from asistencia.jupiter_real import centro_jupiter
    t0, _ = maximo_acercamiento()
    _, d = instante_de_minimo(lambda t: np.linalg.norm(estado_relativo_jupiter(t)[0] - centro_jupiter(t)[0]),
                              t0 - 0.1, t0 + 0.1, 400)
    assert abs(d - 721670) < 30


@pytest.mark.parametrize("luna,fecha,d_pub", [("callisto", "1979-07-08T12:21", 214930),
                                              ("ganymede", "1979-07-09T07:14", 62130),
                                              ("europa", "1979-07-09T17:53", 205720),
                                              ("io", "1979-07-09T23:17", 1129900)])
def test_sobrevuelos_de_los_satelites_galileanos(luna, fecha, d_pub):
    """Valida a la vez el kernel de Voyager 2, el modelo de los satélites y los marcos de
    referencia: la distancia mínima a cada satélite debe coincidir con la tabla publicada
    de la misión (NSSDC/JPL) al 0.3 %, y la hora con 3 minutos."""
    from asistencia.jupiter_real import centro_jupiter, lunas

    def d(t):
        sat = lunas(t)
        return np.linalg.norm(estado_relativo_jupiter(t)[0] - (centro_jupiter(t, sat)[0] + sat[luna][0]))
    t, dmin = instante_de_minimo(d, fecha_a_jd(fecha) - 0.1, fecha_a_jd(fecha) + 0.1, 600)
    assert abs(dmin / d_pub - 1) < 3e-3
    assert abs((t - fecha_a_jd(fecha)) * 1440) < 3.0

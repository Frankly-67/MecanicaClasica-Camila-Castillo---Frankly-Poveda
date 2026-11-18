"""Verificaciones del Modelo 2 (N cuerpos) y de la comparación con las cónicas."""
import numpy as np

from asistencia.constantes import DIA, GM_JUPITER, R_JUPITER, R_SOI_JUPITER
from asistencia.kepler import propagar
from asistencia.modelo2 import (conicas_empalmadas, estado_periapsis_plano, integrar,
                                tiempo_hasta_radio)


def test_jacobi_se_conserva_en_rebound():
    """Invariante exacto del problema restringido circular: h = E - n l_z."""
    r_p, v_p = estado_periapsis_plano(7.6, 10 * R_JUPITER)
    _, _, s = integrar(r_p, v_p, 200 * DIA, muestras=100)
    assert np.ptp(s[:, 1]) / abs(s[0, 1]) < 1e-12


def test_convergencia_ias15():
    """Bajar la tolerancia 1000 veces no cambia la posición final en más de 1 m."""
    r_p, v_p = estado_periapsis_plano(7.6, 10 * R_JUPITER)
    r1, _, _ = integrar(r_p, v_p, 250 * DIA, 1e-8)
    r2, _, _ = integrar(r_p, v_p, 250 * DIA, 1e-11)
    assert np.linalg.norm(r1 - r2) < 1e-3


def test_tiempo_hasta_radio_coincide_con_propagacion():
    """La fórmula de la anomalía hiperbólica y la propagación universal (otro método)
    deben dar la misma distancia."""
    r_p, v_p = estado_periapsis_plano(5.0, 10 * R_JUPITER)
    t = tiempo_hasta_radio(r_p, v_p, R_SOI_JUPITER)
    r, _ = propagar(r_p, v_p, t, GM_JUPITER)
    assert abs(np.linalg.norm(r) / R_SOI_JUPITER - 1) < 1e-9


def test_conicas_y_nbody_coinciden_cerca_del_periapsis():
    """Límite: muy cerca del periapsis domina Júpiter y los dos modelos deben coincidir."""
    r_p, v_p = estado_periapsis_plano(10.0, 10 * R_JUPITER)
    t = 0.5 * DIA
    rN, _, _ = integrar(r_p, v_p, t)
    rP, _, _ = conicas_empalmadas(r_p, v_p, R_SOI_JUPITER, t)
    assert np.linalg.norm(rN - rP) < 50.0          # km, frente a ~4e5 km recorridos


def test_masa_de_jupiter_tendiendo_a_cero_no_hay_asistencia():
    """Sin la masa de Júpiter la energía heliocéntrica no cambia (Kepler puro)."""
    import asistencia.modelo2 as m2
    from asistencia.kepler import energia_especifica
    mu_original = m2.MU2
    try:
        m2.MU2 = 1e-6
        r_p, v_p = np.array([0, -10 * R_JUPITER, 0]), np.array([7.6, 0, 0])
        r1, v1, _ = m2.integrar(r_p, v_p, 100 * DIA)
        r0, v0, _ = m2.integrar(r_p, v_p, -100 * DIA)
        E1, E0 = energia_especifica(r1, v1, m2.MU1), energia_especifica(r0, v0, m2.MU1)
        assert abs(E1 - E0) / abs(E0) < 1e-8
    finally:
        m2.MU2 = mu_original


# ---------- Modelo 2 ampliado (satélites, J2) frente a la trayectoria real ----------

def test_aceleracion_residual_revela_la_maniobra_de_voyager2():
    """Aceleración real (kernel) menos la gravitacional del modelo con satélites y J2: debe
    ser nula salvo durante el encendido de 76 minutos que empezó unas dos horas después del
    máximo acercamiento (Morrison y Samz, NASA SP-439). Dato histórico independiente."""
    from asistencia.jupiter_real import aceleracion_residual
    from asistencia.voyager import maximo_acercamiento
    t_ca, _ = maximo_acercamiento()
    fuera = [np.linalg.norm(aceleracion_residual(t_ca + m / 1440)) for m in (-300, 60, 120, 215, 400)]
    dentro = [np.linalg.norm(aceleracion_residual(t_ca + m / 1440)) for m in (135, 170, 200)]
    assert max(fuera) < 5e-9                       # km/s^2
    assert all(2.5e-6 < a < 3.1e-6 for a in dentro)


def test_modelo_con_satelites_sigue_a_voyager2_tras_la_maniobra():
    """Tramo balístico de 11 días después de la maniobra, partiendo del estado real: con
    satélites y J2 el error final debe ser de pocos km; con Júpiter puntual, cientos."""
    from asistencia.jupiter_real import estado_voyager, integrar as integrar_real
    from asistencia.voyager import cobertura, maximo_acercamiento
    t_ca, _ = maximo_acercamiento()
    t_a, t_b = t_ca + 0.16, cobertura()[1] - 0.01
    err = {}
    for modelo in ("puntual", "satelites+J2"):
        _, y = integrar_real(estado_voyager(t_a), t_a, t_b, modelo)
        err[modelo] = np.linalg.norm(y[:3, -1] - estado_voyager(t_b)[:3])
    assert err["satelites+J2"] < 50.0
    assert err["puntual"] > 10 * err["satelites+J2"]

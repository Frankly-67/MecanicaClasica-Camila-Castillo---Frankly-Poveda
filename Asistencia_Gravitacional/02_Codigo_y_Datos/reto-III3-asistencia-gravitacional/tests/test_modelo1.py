"""Verificaciones independientes del Modelo 1.

Cada prueba compara el código con algo que NO sale del mismo cálculo: una solución
analítica, un invariante, un caso límite o un segundo método numérico.
Ejecutar con:  python -m pytest -v
"""
import numpy as np
import pytest
from scipy.integrate import solve_ivp

from asistencia.constantes import AU, DIA, GM_JUPITER, GM_SOL, R_JUPITER
from asistencia.kepler import energia_especifica, parametro_tisserand, propagar
from asistencia.lambert import lambert
from asistencia.sobrevuelo import angulo_entre, angulo_giro, estado_periapsis

Z = np.array([0.0, 0.0, 1.0])


def integrar_dos_cuerpos(r0, v0, t, mu):
    """Segundo método, independiente de las funciones de Stumpff: integración numérica."""
    def f(_, y):
        return np.concatenate([y[3:], -mu * y[:3] / np.linalg.norm(y[:3]) ** 3])
    sol = solve_ivp(f, (0, t), np.concatenate([r0, v0]), method="DOP853", rtol=1e-12, atol=1e-6)
    return sol.y[:3, -1], sol.y[3:, -1]


# ---------- Propagación kepleriana ----------

@pytest.mark.parametrize("v_factor", [0.8, 1.0, 1.3, 1.6])   # elipse, círculo, elipse, hipérbola
def test_propagacion_coincide_con_integracion_numerica(v_factor):
    r0 = np.array([AU, 0, 0])
    v0 = np.array([0, v_factor * np.sqrt(GM_SOL / AU), 0.3])
    t = 400 * DIA
    r_k, v_k = propagar(r0, v0, t, GM_SOL)
    r_n, v_n = integrar_dos_cuerpos(r0, v0, t, GM_SOL)
    assert np.linalg.norm(r_k - r_n) / AU < 1e-8
    assert np.linalg.norm(v_k - v_n) < 1e-6


def test_propagacion_conserva_energia_y_momento():
    r0 = np.array([AU, 0.2 * AU, 0])
    v0 = np.array([-5.0, 33.0, 0])
    E0, h0 = energia_especifica(r0, v0, GM_SOL), np.cross(r0, v0)
    for dias in [10, 200, 1000, 5000]:
        r, v = propagar(r0, v0, dias * DIA, GM_SOL)
        assert abs(energia_especifica(r, v, GM_SOL) / E0 - 1) < 1e-10
        assert np.linalg.norm(np.cross(r, v) - h0) / np.linalg.norm(h0) < 1e-10


def test_orbita_circular_periodo_completo():
    r0 = np.array([AU, 0, 0])
    v0 = np.array([0, np.sqrt(GM_SOL / AU), 0])
    r, _ = propagar(r0, v0, 2 * np.pi * np.sqrt(AU**3 / GM_SOL), GM_SOL)
    assert np.linalg.norm(r - r0) / AU < 1e-9


# ---------- Lambert ----------

@pytest.mark.parametrize("dias,ang", [(300, 100), (700, 150), (1200, 170), (900, 230)])
def test_lambert_llega_al_punto_pedido(dias, ang):
    """Se resuelve Lambert y se propaga v1 con otro método: debe llegar a r2 a tiempo.
    El caso de 230° es el que hizo fallar la primera versión (bitácora, entrada 4)."""
    r1 = np.array([AU, 0, 0])
    a = np.radians(ang)
    r2 = 5.2 * AU * np.array([np.cos(a), np.sin(a), 0.02])
    v1, v2 = lambert(r1, r2, dias * DIA, GM_SOL)
    r_fin, v_fin = integrar_dos_cuerpos(r1, v1, dias * DIA, GM_SOL)
    assert np.linalg.norm(r_fin - r2) / AU < 1e-7
    assert np.linalg.norm(v_fin - v2) < 1e-5


@pytest.mark.parametrize("dias,ang", [(60, 90), (30, 20)])
def test_lambert_transferencia_muy_rapida(dias, ang):
    """Transferencias hiperbólicas muy rápidas: la raíz cae junto a la zona donde la
    variable universal no tiene solución (y < 0). La versión anterior respondía, sin
    razón, que no había solución (bitácora, entrada 29)."""
    r1 = np.array([AU, 0, 0])
    a = np.radians(ang)
    r2 = 5.2 * AU * np.array([np.cos(a), np.sin(a), 0])
    v1, _ = lambert(r1, r2, dias * DIA, GM_SOL)
    r_fin, _ = integrar_dos_cuerpos(r1, v1, dias * DIA, GM_SOL)
    assert np.linalg.norm(r_fin - r2) / AU < 1e-7


def test_lambert_limite_hohmann():
    """Transferencia casi de 180° con el tiempo de Hohmann: v1 debe ser la de Hohmann."""
    r1n, r2n = AU, 5.2 * AU
    a_t = 0.5 * (r1n + r2n)
    t_h = np.pi * np.sqrt(a_t**3 / GM_SOL)
    v_hohmann = np.sqrt(GM_SOL * (2 / r1n - 1 / a_t))
    a = np.radians(179.99)
    v1, _ = lambert([r1n, 0, 0], r2n * np.array([np.cos(a), np.sin(a), 0]), t_h, GM_SOL)
    assert abs(np.linalg.norm(v1) - v_hohmann) < 0.01       # km/s


# ---------- Sobrevuelo ----------

@pytest.mark.parametrize("vinf,rp_RJ", [(5.6, 10), (8.0, 3), (12.0, 40)])
def test_hiperbola_integrada_confirma_giro_y_simetria(vinf, rp_RJ):
    """Se integra numéricamente la hipérbola desde el periapsis hacia atrás y adelante.
    Lejos de Júpiter |v_inf| debe ser igual en ambos extremos y el ángulo entre las
    velocidades debe ser el delta analítico."""
    mu, rp = GM_JUPITER, rp_RJ * R_JUPITER
    r_p, v_p = estado_periapsis(np.array([0.0, -vinf, 0.0]), rp, mu, Z)
    t = 400 * DIA
    r_a, v_a = integrar_dos_cuerpos(r_p, v_p, t, mu)
    r_b, v_b = integrar_dos_cuerpos(r_p, v_p, -t, mu)
    vinf_a = np.sqrt(v_a @ v_a - 2 * mu / np.linalg.norm(r_a))
    vinf_b = np.sqrt(v_b @ v_b - 2 * mu / np.linalg.norm(r_b))
    assert abs(vinf_a - vinf) < 1e-6 and abs(vinf_b - vinf) < 1e-6
    # a distancia finita la dirección aún no es la asintótica: tolerancia 0.5°
    assert abs(np.degrees(angulo_entre(v_b, v_a) - angulo_giro(vinf, rp, mu))) < 0.5


def test_periapsis_en_3d_gira_alrededor_de_la_normal():
    """Con un plano inclinado, la velocidad de salida debe quedar en el plano de
    v_in y la normal dada, y la distancia mínima debe ser r_p."""
    vin = np.array([1.0, -6.0, 2.0])
    n = np.cross(vin, np.array([3.0, 1.0, -1.0]))
    rp = 8 * R_JUPITER
    r_p, v_p = estado_periapsis(vin, rp, GM_JUPITER, n)
    assert abs(np.linalg.norm(r_p) - rp) < 1e-6
    assert abs(r_p @ v_p) < 1e-6 * rp                     # en el periapsis r ⟂ v
    assert abs(np.dot(v_p, n)) < 1e-9 * np.linalg.norm(n) * np.linalg.norm(v_p)


def test_limite_masa_cero_no_desvia():
    assert angulo_giro(6.0, 10 * R_JUPITER, 1e-12) < 1e-15


def test_limite_periapsis_lejano_no_desvia():
    """delta debe decrecer con r_p y tender a cero (la primera versión de esta prueba,
    escrita por la IA, estaba mal planteada: bitácora, entrada 5)."""
    deltas = [angulo_giro(6.0, rp, GM_JUPITER) for rp in np.logspace(0, 6, 30) * R_JUPITER]
    assert np.all(np.diff(deltas) < 0)
    assert np.degrees(deltas[-1]) < 0.01


def test_analisis_dimensional_excentricidad():
    """e - 1 = r_p v^2 / mu es adimensional: cambiar km -> m no cambia el resultado."""
    e_km = 1 + 10 * R_JUPITER * 5.6**2 / GM_JUPITER
    e_m = 1 + (10 * R_JUPITER * 1e3) * (5.6e3) ** 2 / (GM_JUPITER * 1e9)
    assert abs(e_km - e_m) < 1e-12


def test_tisserand_de_un_planeta_vale_tres():
    """Límite conocido: un cuerpo en la misma órbita circular del planeta tiene T = 3."""
    a = 5.2 * AU
    T = parametro_tisserand([a, 0, 0], [0, np.sqrt(GM_SOL / a), 0], GM_SOL, a)
    assert abs(T - 3) < 1e-12


def test_jacobi_se_conserva_en_el_sistema_rotante():
    """Integración del problema restringido en el marco rotante (otro método): durante un
    sobrevuelo E y l_z cambian, pero h = E - n l_z (Jacobi) se conserva."""
    from scipy.integrate import solve_ivp
    mu1, mu2, a = GM_SOL, GM_JUPITER, 5.2026 * AU
    n = np.sqrt((mu1 + mu2) / a**3)
    al = mu2 / (mu1 + mu2)
    P1, P2 = np.array([-al * a, 0, 0]), np.array([(1 - al) * a, 0, 0])

    def f(_, s):
        r, v = s[:3], s[3:]
        d1, d2 = r - P1, r - P2
        acc = -mu1 * d1 / np.linalg.norm(d1) ** 3 - mu2 * d2 / np.linalg.norm(d2) ** 3
        acc += [2 * n * v[1] + n**2 * r[0], -2 * n * v[0] + n**2 * r[1], 0]
        return np.concatenate([v, acc])

    def h(s):
        r, v = s[:3], s[3:]
        vin = v + np.cross([0, 0, n], r)
        E = vin @ vin / 2 - mu1 / np.linalg.norm(r - P1) - mu2 / np.linalg.norm(r - P2)
        return E, np.cross(r, vin)[2], E - n * np.cross(r, vin)[2]

    rp = 10 * R_JUPITER
    s0 = np.concatenate([P2 + [0, -rp, 0], [np.sqrt(7.6**2 + 2 * mu2 / rp), 0, 0]])
    s1 = solve_ivp(f, (0, 60 * DIA), s0, method="DOP853", rtol=1e-12, atol=1e-8).y[:, -1]
    s2 = solve_ivp(f, (0, -60 * DIA), s0, method="DOP853", rtol=1e-12, atol=1e-8).y[:, -1]
    E1, l1, h1 = h(s1)
    E2, l2, h2 = h(s2)
    assert abs(E1 - E2) > 50                     # la energía heliocéntrica sí cambia
    assert abs(h1 - h2) / abs(h1) < 1e-9          # la de Jacobi no

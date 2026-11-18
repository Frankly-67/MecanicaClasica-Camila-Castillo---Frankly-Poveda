"""Modelo 2: integración gravitacional de N cuerpos con REBOUND (integrador IAS15).

Unidades: km, s y G = 1 (las masas se dan como GM en km^3/s^2).

Para comparar con el Modelo 1 sin mezclar otras diferencias, el caso de estudio usa el
mismo escenario idealizado en ambos modelos: Júpiter en órbita circular alrededor del Sol
(problema restringido circular de tres cuerpos), plano. El encuentro parte del MISMO
estado en el periapsis y se propaga hacia adelante y hacia atrás:

- Modelo 1 (cónicas empalmadas): hipérbola kepleriana alrededor de Júpiter hasta el radio de
  empalme R_e; desde ahí, cónica heliocéntrica (solo el Sol).
- Modelo 2: Sol y Júpiter actúan siempre (REBOUND/IAS15).
"""
import numpy as np
import rebound

from .constantes import A_JUPITER_KM, GM_JUPITER, GM_SOL
from .kepler import propagar

MU1, MU2, A = GM_SOL, GM_JUPITER, A_JUPITER_KM
N_J = np.sqrt((MU1 + MU2) / A**3)          # movimiento medio de la órbita relativa Sol-Júpiter


def jupiter_heliocentrico(t):
    """Júpiter respecto al Sol en órbita circular (relativa, mu = G(M1+M2)); fase 0 en t = 0."""
    c, s = np.cos(N_J * t), np.sin(N_J * t)
    return A * np.array([c, s, 0.0]), A * N_J * np.array([-s, c, 0.0])


def estado_periapsis_plano(vinf, rp, beta_deg=0.0, sentido=+1):
    """Estado relativo a Júpiter en el periapsis. v_inf de llegada apunta en -y (opuesto al
    movimiento de Júpiter) girado beta grados; sentido=+1 gira v_inf hacia la dirección de
    movimiento de Júpiter (la sonda gana energía)."""
    from .sobrevuelo import estado_periapsis
    b = np.radians(beta_deg)
    vin = vinf * np.array([np.sin(b), -np.cos(b), 0.0])
    return estado_periapsis(vin, rp, MU2, np.array([0.0, 0.0, sentido]))


# ----------------------------- Modelo 1 ------------------------------------

def tiempo_hasta_radio(r_p_vec, v_p_vec, R, mu=MU2):
    """Tiempo desde el periapsis hasta |rho| = R en una hipérbola (anomalía hiperbólica F):
    cosh F = (1 + R/|a|)/e,   t = sqrt(|a|^3/mu) (e sinh F - F)."""
    rp = np.linalg.norm(r_p_vec)
    vp = np.linalg.norm(v_p_vec)
    a_abs = mu / (vp**2 - 2 * mu / rp)
    e = 1 + rp / a_abs
    F = np.arccosh((1 + R / a_abs) / e)
    return np.sqrt(a_abs**3 / mu) * (e * np.sinh(F) - F)


def conicas_empalmadas(r_p, v_p, R_e, t_final):
    """Propaga con cónicas empalmadas desde el periapsis (t = 0) hasta t_final (+ o -).
    Devuelve el estado heliocéntrico en t_final y el instante del empalme."""
    signo = np.sign(t_final)
    t_e = signo * tiempo_hasta_radio(r_p, v_p, R_e)
    if abs(t_final) <= abs(t_e):          # todavía dentro de la esfera: solo Júpiter
        rho, vrho = propagar(r_p, v_p, t_final, MU2)
        rJ, vJ = jupiter_heliocentrico(t_final)
        return rJ + rho, vJ + vrho, t_e
    rho, vrho = propagar(r_p, v_p, t_e, MU2)
    rJ, vJ = jupiter_heliocentrico(t_e)
    r, v = propagar(rJ + rho, vJ + vrho, t_final - t_e, MU1)
    return r, v, t_e


# ----------------------------- Modelo 2 ------------------------------------

def simulacion(r_p, v_p, epsilon=1e-9, saturno=False):
    """Sol + Júpiter (órbita circular) + sonda (partícula de prueba) en t = 0."""
    sim = rebound.Simulation()
    sim.G = 1.0
    sim.integrator = "ias15"
    try:                                   # REBOUND 5
        sim.integrator.epsilon = epsilon
    except AttributeError:                 # REBOUND 4
        sim.ri_ias15.epsilon = epsilon
    rJ, vJ = jupiter_heliocentrico(0.0)
    sim.add(m=MU1)
    sim.add(m=MU2, x=rJ[0], y=rJ[1], z=rJ[2], vx=vJ[0], vy=vJ[1], vz=vJ[2])
    if saturno:
        from .constantes import AU
        aS = 9.537 * AU
        GM_S = 8.459715185680659e-08 * 149597870.691**3 / 86400.0**2
        ang = np.radians(28)         # Saturno 28° por delante de Júpiter, como en julio de 1979 (DE421)
        vS = np.sqrt((MU1 + GM_S) / aS)
        sim.add(m=GM_S, x=aS * np.cos(ang), y=aS * np.sin(ang), vx=-vS * np.sin(ang), vy=vS * np.cos(ang))
    sim.N_active = sim.N
    r0, v0 = rJ + r_p, vJ + v_p
    sim.add(m=0.0, x=r0[0], y=r0[1], z=r0[2], vx=v0[0], vy=v0[1], vz=v0[2])
    sim.move_to_com()
    return sim


def estado(sim, i):
    p = sim.particles[i]
    return np.array([p.x, p.y, p.z]), np.array([p.vx, p.vy, p.vz])


def heliocentrico(sim):
    rs, vs = estado(sim, 0)
    r, v = estado(sim, sim.N - 1)
    return r - rs, v - vs


def jacobi(sim):
    """h = E - n l_z de la sonda en el sistema inercial baricéntrico (sección 5 de la
    formulación lagrangiana). Exacta solo para Sol + Júpiter en órbita circular."""
    r, v = estado(sim, sim.N - 1)
    r1, _ = estado(sim, 0)
    r2, _ = estado(sim, 1)
    E = 0.5 * v @ v - MU1 / np.linalg.norm(r - r1) - MU2 / np.linalg.norm(r - r2)
    return E - N_J * np.cross(r, v)[2]


def integrar(r_p, v_p, t_final, epsilon=1e-9, muestras=0, saturno=False):
    """Integra desde el periapsis hasta t_final. Devuelve el estado heliocéntrico final y,
    si muestras > 0, la serie temporal de (t, h, distancia a Júpiter)."""
    sim = simulacion(r_p, v_p, epsilon, saturno)
    if t_final < 0:
        sim.dt = -abs(sim.dt) if sim.dt else -1.0
    serie = []
    tiempos = np.linspace(0, t_final, muestras + 1)[1:] if muestras else [t_final]
    if muestras:
        serie.append((0.0, jacobi(sim), np.linalg.norm(r_p)))
    for t in tiempos:
        sim.integrate(t, exact_finish_time=1)
        if muestras:
            r, _ = estado(sim, sim.N - 1)
            r2, _ = estado(sim, 1)
            serie.append((t, jacobi(sim), np.linalg.norm(r - r2)))
    r, v = heliocentrico(sim)
    return r, v, np.array(serie)

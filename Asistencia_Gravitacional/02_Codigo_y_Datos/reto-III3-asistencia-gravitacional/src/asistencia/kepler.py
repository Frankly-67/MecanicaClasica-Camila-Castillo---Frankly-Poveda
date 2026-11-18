"""Problema de dos cuerpos: funciones de Stumpff, propagación con variable universal
y elementos orbitales.

La propagación usa la formulación de variable universal (Curtis, "Orbital Mechanics
for Engineering Students", alg. 3.3 y 3.4). El mismo código sirve para órbitas
elípticas, parabólicas e hiperbólicas.
"""
import numpy as np


def stumpff_C(z):
    if z > 1e-8:
        return (1 - np.cos(np.sqrt(z))) / z
    if z < -1e-8:
        return (np.cosh(np.sqrt(-z)) - 1) / (-z)
    return 1 / 2 - z / 24 + z * z / 720


def stumpff_S(z):
    if z > 1e-8:
        s = np.sqrt(z)
        return (s - np.sin(s)) / s**3
    if z < -1e-8:
        s = np.sqrt(-z)
        return (np.sinh(s) - s) / s**3
    return 1 / 6 - z / 120 + z * z / 5040


def propagar(r0, v0, dt, mu, tol=1e-10, max_iter=200):
    """Propaga el estado (r0, v0) un tiempo dt bajo un solo cuerpo central de parámetro mu.

    r0, v0: vectores 3D (km, km/s). dt en segundos (puede ser negativo).
    Devuelve (r, v).
    """
    r0 = np.asarray(r0, float)
    v0 = np.asarray(v0, float)
    if dt == 0:
        return r0.copy(), v0.copy()
    r0n = np.linalg.norm(r0)
    vr0 = np.dot(r0, v0) / r0n
    alpha = 2 / r0n - np.dot(v0, v0) / mu      # = 1/a
    sqmu = np.sqrt(mu)

    # Estimación inicial de la anomalía universal. Elipse: Curtis, alg. 3.3.
    # Hipérbola: Vallado (alg. 8), porque la de Curtis hace divergir a Newton cuando
    # dt es grande (lo detectó el barrido del Modelo 2; ver bitácora).
    if alpha < -1e-12:
        a = 1 / alpha
        s = np.sign(dt)
        arg = (-2 * mu * alpha * dt) / (np.dot(r0, v0) + s * np.sqrt(-mu * a) * (1 - r0n * alpha))
        chi = s * np.sqrt(-a) * np.log(arg) if arg > 0 else sqmu * abs(alpha) * dt
    else:
        chi = sqmu * abs(alpha) * dt if abs(alpha) > 1e-12 else sqmu * dt / r0n
    for _ in range(max_iter):
        z = alpha * chi * chi
        C, S = stumpff_C(z), stumpff_S(z)
        F = (r0n * vr0 / sqmu * chi * chi * C + (1 - alpha * r0n) * chi**3 * S
             + r0n * chi - sqmu * dt)
        dF = (r0n * vr0 / sqmu * chi * (1 - alpha * chi * chi * S)
              + (1 - alpha * r0n) * chi * chi * C + r0n)
        paso = F / dF
        chi -= paso
        if abs(paso) < tol * max(1.0, abs(chi)):
            break
    else:
        raise RuntimeError("propagar: Newton no convergió")

    z = alpha * chi * chi
    C, S = stumpff_C(z), stumpff_S(z)
    f = 1 - chi * chi / r0n * C
    g = dt - chi**3 / sqmu * S
    r = f * r0 + g * v0
    rn = np.linalg.norm(r)
    fdot = sqmu / (rn * r0n) * (alpha * chi**3 * S - chi)
    gdot = 1 - chi * chi / rn * C
    v = fdot * r0 + gdot * v0
    return r, v


def energia_especifica(r, v, mu):
    return 0.5 * np.dot(v, v) - mu / np.linalg.norm(r)


def elementos(r, v, mu):
    """Semieje a, excentricidad e, perihelio q, afelio Q (inf si no ligada), energía E,
    momento angular h (vector) e inclinación i (rad, respecto al plano xy)."""
    r = np.asarray(r, float)
    v = np.asarray(v, float)
    E = energia_especifica(r, v, mu)
    h = np.cross(r, v)
    hn = np.linalg.norm(h)
    evec = np.cross(v, h) / mu - r / np.linalg.norm(r)
    e = np.linalg.norm(evec)
    a = -mu / (2 * E)
    q = hn**2 / mu / (1 + e)
    Q = a * (1 + e) if E < 0 else np.inf
    return {"a": a, "e": e, "q": q, "Q": Q, "E": E, "h": h, "i": np.arccos(h[2] / hn)}


def parametro_tisserand(r, v, mu, a_planeta):
    """T = a_p/a + 2 cos(i) sqrt(p / a_p), con p = a (1 - e^2) el semilatus rectum
    (así la fórmula vale para elipses e hipérbolas)."""
    el = elementos(r, v, mu)
    p = np.linalg.norm(el["h"]) ** 2 / mu
    return a_planeta / el["a"] + 2 * np.cos(el["i"]) * np.sqrt(p / a_planeta)

"""Problema de Lambert: dados r1, r2 y el tiempo de vuelo, hallar v1 y v2.

Formulación de variable universal (Curtis, alg. 5.2), una sola revolución,
movimiento prógrado (antihorario visto desde el polo norte de la eclíptica).
Se busca primero un intervalo con cambio de signo sobre una malla de z y luego se
refina con el método de Brent, que no depende de una buena estimación inicial.
(La primera versión buscaba el intervalo de otra forma y fallaba para transferencias
de más de 180°; ver docs/bitacora_ia.md, entrada 4.)
"""
import numpy as np
from scipy.optimize import brentq

from .kepler import stumpff_C, stumpff_S

Z_MAX = (2 * np.pi) ** 2 - 1e-6   # límite superior para una revolución

_ZS = np.concatenate([-np.logspace(4, -3, 120), np.linspace(1e-6, Z_MAX, 300)])
_ZP = np.sqrt(np.clip(_ZS, 0, None))
_ZN = np.sqrt(np.clip(-_ZS, 0, None))
with np.errstate(all="ignore"):
    _C = np.where(_ZS > 0, (1 - np.cos(_ZP)) / _ZS, (np.cosh(_ZN) - 1) / -_ZS)
    _S = np.where(_ZS > 0, (_ZP - np.sin(_ZP)) / _ZP**3, (np.sinh(_ZN) - _ZN) / _ZN**3)


def _F_malla(r1n, r2n, A, sqmu_tof):
    y = r1n + r2n + A * (_ZS * _S - 1) / np.sqrt(_C)
    with np.errstate(invalid="ignore"):
        F = (y / _C) ** 1.5 * _S + A * np.sqrt(y) - sqmu_tof
    return np.where(y < 0, -1e30, F)


def lambert(r1, r2, tof, mu, progrado=True):
    r1 = np.asarray(r1, float)
    r2 = np.asarray(r2, float)
    r1n, r2n = np.linalg.norm(r1), np.linalg.norm(r2)
    dth = np.arccos(np.clip(np.dot(r1, r2) / (r1n * r2n), -1, 1))
    cz = np.cross(r1, r2)[2]
    if (progrado and cz < 0) or (not progrado and cz >= 0):
        dth = 2 * np.pi - dth
    if abs(np.sin(dth)) < 1e-10:
        raise ValueError("Transferencia de 0° o 180°: el plano de la órbita no está definido")
    A = np.sin(dth) * np.sqrt(r1n * r2n / (1 - np.cos(dth)))
    sqmu = np.sqrt(mu)

    def y(z):
        return r1n + r2n + A * (z * stumpff_S(z) - 1) / np.sqrt(stumpff_C(z))

    def F(z):
        yz = y(z)
        if yz < 0:
            return -1e30
        return (yz / stumpff_C(z)) ** 1.5 * stumpff_S(z) + A * np.sqrt(yz) - sqmu * tof

    Fs = _F_malla(r1n, r2n, A, sqmu * tof)
    ok = Fs > -1e29
    # Cambio de signo entre puntos vecinos. El de la izquierda puede caer en la zona sin
    # solución (y < 0, donde F vale -1e30): como F -> -sqrt(mu) tof < 0 cuando y -> 0+, también
    # ahí hay una raíz. (Antes se exigían los dos puntos válidos y se perdía la solución de
    # las transferencias muy rápidas; bitácora, entrada 29.)
    idx = np.where(ok[1:] & (np.sign(Fs[:-1]) != np.sign(Fs[1:])))[0]
    if len(idx) == 0:
        raise ValueError("Lambert: no hay solución de una revolución para ese tiempo de vuelo")
    i = idx[0]
    z = brentq(F, _ZS[i], _ZS[i + 1], xtol=1e-12, rtol=1e-12, maxiter=200)

    yz = y(z)
    f = 1 - yz / r1n
    g = A * np.sqrt(yz / mu)
    gdot = 1 - yz / r2n
    return (r2 - f * r1) / g, (gdot * r2 - r1) / g

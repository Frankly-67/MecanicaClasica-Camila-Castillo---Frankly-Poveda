"""Verificación independiente de la formulación lagrangiana (problema restringido de 3 cuerpos).

Parte A (simbólica, sympy):
  1. Ecuaciones de Euler–Lagrange del lagrangiano en el sistema rotante.
  2. dh/dt = 0 sobre las soluciones (h = función de energía de Jacobi).
  3. h = E - n*l_z (energía inercial menos n por momento angular inercial).
  4. Límite lejano: -2 h a / GM1 = parámetro de Tisserand.
Parte B (numérica): un sobrevuelo de Júpiter integrado en el sistema rotante con DOP853:
  E y l_z cambian mucho, E - n l_z se conserva a precisión numérica.

Salida: resultados/verificacion_lagrangiano.json
"""
import json

import numpy as np
import sympy as sp
from estilo import RES
from scipy.integrate import solve_ivp

from asistencia.constantes import AU, DIA, GM_JUPITER, GM_SOL, R_JUPITER

# ------------------------- Parte A: simbólica -------------------------
t = sp.symbols("t")
n, G1, G2, a, al = sp.symbols("n G_1 G_2 a alpha", positive=True)
x, y, z = (sp.Function(s)(t) for s in "xyz")
r1 = sp.sqrt((x + al * a) ** 2 + y**2 + z**2)          # Sol en (-alpha a, 0, 0)
r2 = sp.sqrt((x - (1 - al) * a) ** 2 + y**2 + z**2)    # Júpiter en ((1-alpha) a, 0, 0)
V = -G1 / r1 - G2 / r2
xd, yd, zd = x.diff(t), y.diff(t), z.diff(t)
T = sp.Rational(1, 2) * ((xd - n * y) ** 2 + (yd + n * x) ** 2 + zd**2)
L = T - V

el = [sp.simplify(sp.diff(L.diff(q.diff(t)), t) - L.diff(q)) for q in (x, y, z)]
Om = n**2 * (x**2 + y**2) / 2 - V
esperado = [x.diff(t, 2) - 2 * n * yd - Om.diff(x),
            y.diff(t, 2) + 2 * n * xd - Om.diff(y),
            z.diff(t, 2) - Om.diff(z)]
ok_el = all(sp.simplify(e - f) == 0 for e, f in zip(el, esperado))

h = sum(q.diff(t) * L.diff(q.diff(t)) for q in (x, y, z)) - L
h_esperado = (xd**2 + yd**2 + zd**2) / 2 - Om
ok_h = sp.simplify(h - h_esperado) == 0

# dh/dt sobre las soluciones: sustituir las aceleraciones de Euler–Lagrange
sol = sp.solve(esperado, [x.diff(t, 2), y.diff(t, 2), z.diff(t, 2)], dict=True)[0]
ok_conserva = sp.simplify(h.diff(t).subs(sol)) == 0

# h = E - n l_z con la velocidad inercial v = (xd - n y, yd + n x, zd)
vx, vy = xd - n * y, yd + n * x
E = (vx**2 + vy**2 + zd**2) / 2 + V
lz = x * vy - y * vx
ok_E_nl = sp.simplify(h - (E - n * lz)) == 0

# Límite lejano de Júpiter (G2 -> 0, alpha -> 0): elementos keplerianos del Sol
As, es, inc = sp.symbols("a_s e i", positive=True)
E_kep = -G1 / (2 * As)
lz_kep = sp.sqrt(G1 * As * (1 - es**2)) * sp.cos(inc)
n_kep = sp.sqrt(G1 / a**3)
TJ = sp.simplify(-2 * (E_kep - n_kep * lz_kep) * a / G1)
TJ_esperado = a / As + 2 * sp.cos(inc) * sp.sqrt(As * (1 - es**2) / a)
ok_tisserand = sp.simplify(TJ - TJ_esperado) == 0

# ------------------------- Parte B: numérica --------------------------
mu1, mu2 = GM_SOL, GM_JUPITER
aJ = 5.2026 * AU
nJ = np.sqrt((mu1 + mu2) / aJ**3)
alfa = mu2 / (mu1 + mu2)
P1, P2 = np.array([-alfa * aJ, 0, 0]), np.array([(1 - alfa) * aJ, 0, 0])


def f(_, s):
    r, v = s[:3], s[3:]
    d1, d2 = r - P1, r - P2
    acc = -mu1 * d1 / np.linalg.norm(d1) ** 3 - mu2 * d2 / np.linalg.norm(d2) ** 3
    acc += np.array([2 * nJ * v[1] + nJ**2 * r[0], -2 * nJ * v[0] + nJ**2 * r[1], 0.0])
    return np.concatenate([v, acc])


def cantidades(s):
    r, v = s[:3], s[3:]
    vin = v + np.cross([0, 0, nJ], r)
    Vp = -mu1 / np.linalg.norm(r - P1) - mu2 / np.linalg.norm(r - P2)
    E = vin @ vin / 2 + Vp
    lz = np.cross(r, vin)[2]
    return E, lz, E - nJ * lz


# Condición inicial: sonda que pasa a 10 R_J de Júpiter, empezando a 0.8 UA de él
rp = 10 * R_JUPITER
vinf = 7.6
vp = np.sqrt(vinf**2 + 2 * mu2 / rp)
s_p = np.concatenate([P2 + [0, -rp, 0], [vp, 0, 0]])   # periapsis (velocidad relativa al marco rotante ~ inercial local)
dur = 120 * DIA
ida = solve_ivp(f, (0, dur), s_p, method="DOP853", rtol=1e-13, atol=1e-9, dense_output=True)
vuelta = solve_ivp(f, (0, -dur), s_p, method="DOP853", rtol=1e-13, atol=1e-9, dense_output=True)
a_ini, a_fin = cantidades(vuelta.y[:, -1]), cantidades(ida.y[:, -1])
ts = np.linspace(-dur, dur, 2001)
hs = [cantidades(vuelta.sol(t) if t < 0 else ida.sol(t))[2] for t in ts]

salida = {
    "simbolica": {"euler_lagrange_coincide": ok_el, "h_igual_v2/2_menos_Omega": ok_h,
                  "dh_dt_cero_sobre_soluciones": ok_conserva, "h_igual_E_menos_n_lz": ok_E_nl,
                  "limite_lejano_es_tisserand": ok_tisserand},
    "numerica": {
        "duracion_dias": 240,
        "E_inicio_km2_s2": round(a_ini[0], 4), "E_fin_km2_s2": round(a_fin[0], 4),
        "lz_inicio_km2_s": float(f"{a_ini[1]:.6e}"), "lz_fin_km2_s": float(f"{a_fin[1]:.6e}"),
        "h_inicio": round(a_ini[2], 8), "h_fin": round(a_fin[2], 8),
        "variacion_relativa_maxima_h": float(f"{(max(hs) - min(hs)) / abs(hs[0]):.2e}"),
        "cambio_E_igual_n_por_cambio_lz": {
            "Delta_E": round(a_fin[0] - a_ini[0], 4),
            "n_Delta_lz": round(nJ * (a_fin[1] - a_ini[1]), 4)},
    },
}
(RES / "verificacion_lagrangiano.json").write_text(json.dumps(salida, indent=2, ensure_ascii=False))
print(json.dumps(salida, indent=2, ensure_ascii=False))

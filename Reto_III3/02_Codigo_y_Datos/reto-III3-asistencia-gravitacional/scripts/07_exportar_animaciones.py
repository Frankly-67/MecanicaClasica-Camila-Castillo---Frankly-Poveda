"""Exporta los datos que usa la página de animaciones (docs/animaciones/datos.js).

  1. Sobrevuelo simulado con REBOUND (v_inf = 7.6 km/s, r_p = 10 R_J, 250 días a cada lado):
     posición heliocéntrica de la sonda y de Júpiter, E, n*l_z y h = E - n*l_z.
  2. Voyager 2 real (kernel SPICE) y el Modelo 2 con y sin satélites, J2 y maniobra.
  3. Error de las cónicas (resultados/modelo2.json).
No se calcula nada nuevo: se reutilizan las mismas funciones que los scripts 05 y 06.
Salida: docs/animaciones/datos.js
"""
import json
from pathlib import Path

import numpy as np
from estilo import RES

from asistencia.constantes import AU, DIA, GM_JUPITER, R_JUPITER, R_SOI_JUPITER
from asistencia.efemerides import jd_a_fecha
from asistencia.jupiter_real import centro_jupiter, cuerpos, estado_voyager, integrar, lunas
from asistencia.kepler import energia_especifica
from asistencia.modelo2 import MU1, N_J, estado, estado_periapsis_plano, simulacion
from asistencia.voyager import cobertura, estado_relativo_jupiter, instante_de_minimo

SALIDA = Path(__file__).resolve().parent.parent / "docs" / "animaciones" / "datos.js"


def r(x, n=3):
    return [round(float(v), n) for v in np.ravel(x)]


def sobrevuelo_nbody():
    r_p, v_p = estado_periapsis_plano(7.6, 10 * R_JUPITER)
    T = 250 * DIA
    tiempos = np.linspace(-T, T, 501)
    filas = {k: [] for k in ("t", "px", "py", "jx", "jy", "E", "nLz", "h")}
    for sentido, ts in ((-1, tiempos[tiempos < 0][::-1]), (+1, tiempos[tiempos >= 0])):
        sim = simulacion(r_p, v_p, 1e-9, False)
        if sentido < 0:
            sim.dt = -abs(sim.dt) if sim.dt else -1.0
        for t in ts:
            sim.integrate(t, exact_finish_time=1)
            rp_, vp_ = estado(sim, sim.N - 1)
            rs, _ = estado(sim, 0)
            rj, _ = estado(sim, 1)
            E = 0.5 * vp_ @ vp_ - MU1 / np.linalg.norm(rp_ - rs) - GM_JUPITER / np.linalg.norm(rp_ - rj)
            nlz = N_J * np.cross(rp_, vp_)[2]
            for k, v in (("t", t / DIA), ("px", (rp_ - rs)[0] / AU), ("py", (rp_ - rs)[1] / AU),
                         ("jx", (rj - rs)[0] / AU), ("jy", (rj - rs)[1] / AU), ("E", E), ("nLz", nlz), ("h", E - nlz)):
                filas[k].append(float(v))
    orden = np.argsort(filas["t"])
    return {k: r(np.array(v)[orden], 5 if k in ("px", "py", "jx", "jy") else 4) for k, v in filas.items()}


def voyager():
    a, b = cobertura()
    t0, t1 = a + 0.01, b - 0.01
    t_ca, _ = instante_de_minimo(lambda t: np.linalg.norm(estado_relativo_jupiter(t)[0]), t0, t1, 8000)
    # muestreo más fino cerca del periapsis
    jd = np.unique(np.concatenate([np.linspace(t0, t1, 400), np.linspace(t_ca - 1.5, t_ca + 1.0, 500)]))
    rho0, v0 = estado_relativo_jupiter(t_ca)
    # plano de la órbita jovicéntrica (para proyectar)
    nrm = np.cross(rho0, v0); nrm /= np.linalg.norm(nrm)
    e1 = rho0 / np.linalg.norm(rho0); e2 = np.cross(nrm, e1)

    def proy(x):
        return [round(float(x @ e1 / R_JUPITER), 3), round(float(x @ e2 / R_JUPITER), 3)]
    y0 = estado_voyager(t0)
    ts = (jd - t0) * DIA
    ini, fin = t_ca + 129.9 / 1440, t_ca + 206.5 / 1440
    from asistencia.jupiter_real import aceleracion_residual
    tm = np.linspace(ini, fin, 153)[1:-1]
    a_med = np.mean([aceleracion_residual(t) for t in tm], axis=0)
    modelos = {}
    for nombre, modelo, man in (("puntual", "puntual", None), ("completo", "satelites+J2", (ini, fin, a_med))):
        _, y = integrar(y0, t0, t1, modelo, man, t_eval=ts)
        modelos[nombre] = y

    def vinf_eq(y, t):
        (_, _), (rj, vj), _ = cuerpos(t)
        return float(np.sqrt(2 * energia_especifica(y[:3] - rj, y[3:] - vj, GM_JUPITER)))
    real = np.array([estado_voyager(t) for t in jd]).T
    out = {"dias": r(jd - t_ca, 4), "real": [], "puntual": [], "completo": [], "io": [], "europa": [], "ganymede": [], "callisto": [],
           "vinf_real": [], "vinf_puntual": [], "vinf_completo": []}
    for k, t in enumerate(jd):
        c = cuerpos(t)[1][0]               # baricentro de Júpiter (posición)
        for nom, y in (("real", real), ("puntual", modelos["puntual"]), ("completo", modelos["completo"])):
            out[nom].append(proy(y[:3, k] - c))
        sat = lunas(t)
        cen = centro_jupiter(t, sat)[0]
        for n in ("io", "europa", "ganymede", "callisto"):
            out[n].append(proy(cen + sat[n][0]))
        out["vinf_real"].append(round(vinf_eq(real[:, k], t), 4))
        out["vinf_puntual"].append(round(vinf_eq(modelos["puntual"][:, k], t), 4))
        out["vinf_completo"].append(round(vinf_eq(modelos["completo"][:, k], t), 4))
    out["maniobra_dias"] = [round(float(ini - t_ca), 5), round(float(fin - t_ca), 5)]
    out["periapsis_utc"] = jd_a_fecha(t_ca, segundos=True)
    out["inicio_utc"] = jd_a_fecha(t0, hora=True)
    out["fin_utc"] = jd_a_fecha(t1, hora=True)
    out["r_soi_RJ"] = round(R_SOI_JUPITER / R_JUPITER, 1)
    return out


def main():
    m2 = json.loads((RES / "modelo2.json").read_text())
    err = m2["C_error_conicas_vs_rp_(Re=rSOI)_en_%"]
    datos = {"nbody": sobrevuelo_nbody(), "voyager": voyager(), "error_conicas": err,
             "constantes": {"GM_J": GM_JUPITER, "RJ": R_JUPITER, "rSOI_RJ": R_SOI_JUPITER / R_JUPITER, "n": N_J}}
    SALIDA.parent.mkdir(parents=True, exist_ok=True)
    SALIDA.write_text("window.DATOS = " + json.dumps(datos, separators=(",", ":")) + ";\n")
    print("escrito", SALIDA, SALIDA.stat().st_size // 1024, "KB")


if __name__ == "__main__":
    main()

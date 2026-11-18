"""Modelo 1: cónicas empalmadas Tierra -> Júpiter -> Saturno.

  1. Tramo heliocéntrico Tierra -> Júpiter (problema de Lambert, solo el Sol).
  2. Sobrevuelo en el sistema de Júpiter (hipérbola, solo Júpiter; instantáneo visto
     desde el Sol, porque la esfera de influencia es pequeña frente a la órbita).
  3. Tramo heliocéntrico Júpiter -> Saturno (Lambert, solo el Sol).

Condiciones para un sobrevuelo sin motor:
  |v_inf_in| = |v_inf_out|   (energía conservada en el sistema de Júpiter)
  r_p >= r_min               (el periapsis no puede estar dentro de Júpiter)
"""
import numpy as np
from scipy.optimize import minimize

from .constantes import DIA, GM_JUPITER, GM_SOL, R_JUPITER
from .efemerides import estado_planeta
from .lambert import lambert
from .sobrevuelo import analizar_sobrevuelo


def evaluar_mision(jd_salida, tof1_dias, tof2_dias, fuente=None):
    jd_J = jd_salida + tof1_dias
    jd_S = jd_J + tof2_dias
    rE, vE = estado_planeta("Tierra", jd_salida, fuente)
    rJ, vJ = estado_planeta("Jupiter", jd_J, fuente)
    rS, vS = estado_planeta("Saturno", jd_S, fuente)

    v1, v2 = lambert(rE, rJ, tof1_dias * DIA, GM_SOL)
    v3, v4 = lambert(rJ, rS, tof2_dias * DIA, GM_SOL)

    vinf_in, vinf_out = v2 - vJ, v3 - vJ
    fb = analizar_sobrevuelo(vinf_in, vinf_out, GM_JUPITER)
    return {
        "jd_salida": jd_salida, "jd_jupiter": jd_J, "jd_saturno": jd_S,
        "tof1": tof1_dias, "tof2": tof2_dias,
        "C3": float(np.dot(v1 - vE, v1 - vE)),
        "rE": rE, "rJ": rJ, "rS": rS, "vE": vE, "vJ": vJ, "vS": vS,
        "v_salida": v1, "v_llega_J": v2, "v_sale_J": v3, "v_llega_S": v4,
        "vinf_in_vec": vinf_in, "vinf_out_vec": vinf_out,
        "vinf_saturno": float(np.linalg.norm(v4 - vS)),
        **fb,
        "rp_RJ": fb["rp"] / R_JUPITER,
    }


def refinar(jd0, tof1, tof2, rp_min_RJ=2.0, c3_max=None, fuente=None):
    """Ajusta (fecha, tof1, tof2) para que el desajuste de |v_inf| sea exactamente cero.

    Sin c3_max minimiza el C3 de salida; con c3_max minimiza el tiempo total hasta
    Saturno sujeto a C3 <= c3_max. Siempre exige r_p >= rp_min_RJ.
    """
    x0 = np.array([jd0, tof1, tof2], float)
    escala = np.array([10.0, 10.0, 20.0])
    cache = {}

    def mision(x):
        k = tuple(np.round(x, 12))
        if k not in cache:
            jd, t1, t2 = x0 + escala * x
            try:
                cache[k] = evaluar_mision(jd, t1, t2, fuente=fuente)
            except ValueError:
                # Lambert sin solución en ese punto: valores que el optimizador rechaza
                cache[k] = {"C3": 1e4, "desajuste": 50.0, "rp_RJ": -100.0, "tof1": t1, "tof2": t2}
        return cache[k]

    if c3_max is None:
        objetivo = lambda x: mision(x)["C3"] / 100.0  # noqa: E731
    else:
        objetivo = lambda x: (mision(x)["tof1"] + mision(x)["tof2"]) / 1000.0  # noqa: E731
    restricciones = [
        {"type": "eq", "fun": lambda x: mision(x)["desajuste"]},
        {"type": "ineq", "fun": lambda x: mision(x)["rp_RJ"] - rp_min_RJ},
    ]
    if c3_max is not None:
        restricciones.append({"type": "ineq", "fun": lambda x: (c3_max - mision(x)["C3"]) / 10.0})
    res = minimize(objetivo, np.zeros(3), method="SLSQP", constraints=restricciones,
                   options={"ftol": 1e-10, "maxiter": 300})
    return mision(res.x), res

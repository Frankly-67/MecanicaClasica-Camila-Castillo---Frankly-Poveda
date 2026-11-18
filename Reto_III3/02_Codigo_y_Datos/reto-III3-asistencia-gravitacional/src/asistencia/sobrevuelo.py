"""Geometría del sobrevuelo en el sistema de referencia de Júpiter.

En ese sistema, durante el encuentro, la sonda sigue una hipérbola. Se conservan la
energía y el momento angular respecto a Júpiter, así que |v_inf| no cambia: el
encuentro solo gira el vector v_inf un ángulo delta.

    e = 1 + r_p v_inf^2 / mu,      sin(delta / 2) = 1 / e
"""
import numpy as np


def angulo_giro(v_inf, r_p, mu):
    e = 1 + r_p * v_inf**2 / mu
    return 2 * np.arcsin(1 / e)


def periapsis_requerido(delta, v_inf, mu):
    """Distancia de máximo acercamiento necesaria para girar v_inf un ángulo delta."""
    return mu * (1 / np.sin(delta / 2) - 1) / v_inf**2


def angulo_entre(u, w):
    c = np.dot(u, w) / (np.linalg.norm(u) * np.linalg.norm(w))
    return np.arccos(np.clip(c, -1, 1))


def rotar(v, eje, ang):
    """Rota v un ángulo ang alrededor del eje unitario (fórmula de Rodrigues)."""
    k = eje / np.linalg.norm(eje)
    return v * np.cos(ang) + np.cross(k, v) * np.sin(ang) + k * np.dot(k, v) * (1 - np.cos(ang))


def analizar_sobrevuelo(v_inf_in, v_inf_out, mu):
    """Dado lo que piden los tramos heliocéntricos (v_inf de llegada y de salida),
    devuelve lo que tendría que hacer el sobrevuelo sin propulsión."""
    vin, vout = np.linalg.norm(v_inf_in), np.linalg.norm(v_inf_out)
    delta = angulo_entre(v_inf_in, v_inf_out)
    v_med = 0.5 * (vin + vout)
    rp = periapsis_requerido(delta, v_med, mu) if delta > 1e-9 else np.inf
    return {"vinf_in": vin, "vinf_out": vout, "desajuste": vout - vin,
            "delta": delta, "rp": rp, "normal": np.cross(v_inf_in, v_inf_out),
            "e": 1 / np.sin(delta / 2) if delta > 1e-9 else np.inf}


def estado_periapsis(v_inf_in, r_p, mu, normal):
    """Estado (r, v) relativo a Júpiter en el periapsis de la hipérbola que llega con
    v_inf_in, pasa a distancia r_p y gira v_inf alrededor del eje `normal`
    (normal = v_inf_in x v_inf_out da el giro pedido)."""
    v_inf = np.linalg.norm(v_inf_in)
    delta = angulo_giro(v_inf, r_p, mu)
    u = rotar(v_inf_in / v_inf, normal, delta / 2)        # dirección de v en el periapsis
    r_dir = rotar(u, normal, -np.pi / 2)
    return r_p * r_dir, np.sqrt(v_inf**2 + 2 * mu / r_p) * u

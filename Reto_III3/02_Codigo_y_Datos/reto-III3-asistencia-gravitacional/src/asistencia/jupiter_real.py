"""Modelo 2 ampliado para comparar con la trayectoria real de Voyager 2.

El Modelo 2 básico trata a Júpiter como una masa puntual (GM del sistema joviano en su
baricentro). Aquí se levanta esa hipótesis (H5): Júpiter se separa en el planeta, con su
achatamiento (J2), y los cuatro satélites galileanos. El Sol, el baricentro de Júpiter y
Saturno se toman de la efeméride DE421 (posiciones prescritas), así que solo se integra
la sonda, en el sistema inercial del baricentro del Sistema Solar (eclíptica J2000).

Fuentes:
- Posiciones de los satélites: paquete `astronomy-engine` (D. Cross), función
  `JupiterMoons`, que da vectores respecto al centro de Júpiter. Se validan en
  `tests/test_datos.py` contra las distancias publicadas de los sobrevuelos de Voyager 2.
- GM de los satélites, J2 y radio de referencia: solución JUP230 de JPL (R. A. Jacobson).
- Polo de Júpiter: IAU (ascensión recta 268.057°, declinación 64.495°, J2000).
"""
import numpy as np
from scipy.integrate import solve_ivp

from .constantes import AU, DIA, GM_JUPITER, GM_SOL, R_JUPITER
from .efemerides import _efemeride, ecuatorial_a_ecliptica
from .voyager import estado_relativo_jupiter

GM_LUNAS = {"io": 5959.916, "europa": 3202.739, "ganymede": 9887.834, "callisto": 7179.289}  # km^3/s^2
NOMBRES = {"io": "Ío", "europa": "Europa", "ganymede": "Ganímedes", "callisto": "Calisto"}
GM_PLANETA = GM_JUPITER - sum(GM_LUNAS.values())      # Júpiter solo, sin satélites
GM_SATURNO = 8.459715185680659e-08 * 149597870.691**3 / DIA**2       # DE421
J2 = 14696.43e-6                                      # con radio de referencia 71 492 km
_AR, _DEC = np.radians(268.056595), np.radians(64.495303)
POLO = ecuatorial_a_ecliptica(np.array([np.cos(_DEC) * np.cos(_AR), np.cos(_DEC) * np.sin(_AR),
                                        np.sin(_DEC)]))
MODELOS = ("puntual", "satelites", "satelites+J2")


def _de421(nombre, jd):
    r, v = _efemeride().position_and_velocity(nombre, jd)
    return ecuatorial_a_ecliptica(r[:, 0]), ecuatorial_a_ecliptica(v[:, 0] / DIA)


def cuerpos(jd):
    """Estados (r, v) del Sol, del baricentro de Júpiter y de Saturno respecto al baricentro
    del Sistema Solar, eclíptica J2000 (km, km/s)."""
    return _de421("sun", jd), _de421("jupiter", jd), _de421("saturn", jd)


def lunas(jd):
    """Posición y velocidad de cada satélite respecto al CENTRO de Júpiter (km, km/s)."""
    import astronomy
    jm = astronomy.JupiterMoons(astronomy.Time.FromTerrestrialTime(jd - 2451545.0))
    out = {}
    for n in GM_LUNAS:
        s = getattr(jm, n)
        out[n] = (ecuatorial_a_ecliptica(np.array([s.x, s.y, s.z]) * AU),
                  ecuatorial_a_ecliptica(np.array([s.vx, s.vy, s.vz]) * AU / DIA))
    return out


def centro_jupiter(jd, sat=None):
    """Centro del planeta respecto al baricentro del sistema joviano (r, v). El baricentro
    es el promedio pesado por masa del planeta y sus satélites, de modo que
    r_centro = - sum(GM_i r_i) / GM_sistema, con r_i respecto al centro."""
    sat = sat or lunas(jd)
    r = -sum(GM_LUNAS[n] * sat[n][0] for n in GM_LUNAS) / GM_JUPITER
    v = -sum(GM_LUNAS[n] * sat[n][1] for n in GM_LUNAS) / GM_JUPITER
    return r, v


def aceleracion(jd, r, modelo="satelites+J2"):
    """Aceleración gravitacional de la sonda en la posición r (respecto al baricentro del
    Sistema Solar). modelo: 'puntual', 'satelites' o 'satelites+J2'."""
    (rs, _), (rj, _), (rsa, _) = cuerpos(jd)
    a = -GM_SOL * (r - rs) / np.linalg.norm(r - rs) ** 3
    a -= GM_SATURNO * (r - rsa) / np.linalg.norm(r - rsa) ** 3
    if modelo == "puntual":
        d = r - rj
        return a - GM_JUPITER * d / np.linalg.norm(d) ** 3
    sat = lunas(jd)
    rc = rj + centro_jupiter(jd, sat)[0]
    d = r - rc
    dn = np.linalg.norm(d)
    a -= GM_PLANETA * d / dn**3
    if modelo == "satelites+J2":
        z = d @ POLO
        a -= 1.5 * J2 * GM_PLANETA * R_JUPITER**2 / dn**5 * ((1 - 5 * z * z / dn**2) * d + 2 * z * POLO)
    for n, gm in GM_LUNAS.items():
        dl = r - (rc + sat[n][0])
        a -= gm * dl / np.linalg.norm(dl) ** 3
    return a


def estado_voyager(jd):
    """Estado real de Voyager 2 (kernel SPICE) respecto al baricentro del Sistema Solar."""
    rj, vj = _de421("jupiter", jd)
    rho, vrho = estado_relativo_jupiter(jd)
    return np.concatenate([rj + rho, vj + vrho])


def integrar(y0, jd0, jd1, modelo="satelites+J2", maniobra=None, t_eval=None):
    """Integra la sonda desde jd0 hasta jd1 con DOP853 (rtol = 1e-12).

    maniobra = (jd_inicio, jd_fin, a_vec) añade una aceleración constante a_vec (km/s^2)
    entre esas fechas; la integración se parte en tramos para que el integrador no cruce
    la discontinuidad del empuje. t_eval: instantes (segundos desde jd0) donde se quiere el
    estado. Devuelve (t, y): y tiene una columna por instante; sin t_eval, solo el final."""
    T = (jd1 - jd0) * DIA
    cortes = [0.0, T]
    if maniobra is not None and jd0 < maniobra[0] < maniobra[1] < jd1:
        cortes = [0.0, (maniobra[0] - jd0) * DIA, (maniobra[1] - jd0) * DIA, T]
    y = np.asarray(y0, float)
    ts, ys = [], []
    for k in range(len(cortes) - 1):
        a, b = cortes[k], cortes[k + 1]
        empuje = maniobra[2] if (len(cortes) == 4 and k == 1) else 0.0

        def f(t, y, empuje=empuje):
            return np.concatenate([y[3:], aceleracion(jd0 + t / DIA, y[:3], modelo) + empuje])
        s = solve_ivp(f, (a, b), y, method="DOP853", rtol=1e-12, atol=1e-6,
                      dense_output=t_eval is not None)
        if t_eval is not None:
            ultimo = k == len(cortes) - 2
            m = (t_eval >= min(a, b)) & ((t_eval <= max(a, b)) if ultimo else (t_eval < max(a, b)))
            ts.append(t_eval[m])
            ys.append(s.sol(t_eval[m]))
        y = s.y[:, -1]
    if t_eval is None:
        return np.array([T]), y[:, None]
    return np.concatenate(ts), np.concatenate(ys, axis=1)


def aceleracion_residual(jd, modelo="satelites+J2", h=5.0):
    """Aceleración real de la sonda (derivada numérica de la velocidad del kernel, diferencia
    centrada de paso h segundos) menos la gravitacional del modelo. Si el modelo de fuerzas
    es completo debe ser cero; lo que quede es una fuerza no gravitacional (o no modelada)."""
    a_real = (estado_voyager(jd + h / DIA)[3:] - estado_voyager(jd - h / DIA)[3:]) / (2 * h)
    return a_real - aceleracion(jd, estado_voyager(jd)[:3], modelo)

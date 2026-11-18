"""Voyager 2 con el Modelo 2 ampliado: satélites galileanos, J2 y la maniobra de corrección.

El Modelo 2 básico (05_modelo2.py, Júpiter puntual) no reproduce la caída de la energía
jovicéntrica real cerca del periapsis y se separa 68 000 km en 11 días. Aquí se busca la causa:

  A. Validación de las efemérides de los satélites: distancias de máximo acercamiento de
     Voyager 2 a cada satélite frente a los valores publicados.
  B. Aceleración residual (real menos gravitacional): si el modelo de fuerzas es completo debe
     ser cero. Aparece una aceleración constante durante 76 minutos: una maniobra.
  C. Integración desde el estado real del 4 jun 1979 con cada modelo de Júpiter, con y sin
     la maniobra, y por tramos balísticos (antes y después de la maniobra).
  D. De dónde sale la caída de la energía jovicéntrica.
  E. Prolongación de la trayectoria fuera de la esfera de influencia: energía heliocéntrica,
     excentricidad y parámetro de Tisserand asintóticos.

Salidas: resultados/voyager2_satelites.json, figuras/fig09_voyager2_satelites.png
"""
import json

import numpy as np
from estilo import COL, FIG, RES, figura, titulo
from scipy.optimize import brentq

from asistencia.constantes import A_JUPITER_KM, AU, DIA, GM_JUPITER, GM_SOL, R_JUPITER, R_SOI_JUPITER
from asistencia.efemerides import jd_a_fecha
from asistencia.jupiter_real import (NOMBRES, aceleracion_residual, centro_jupiter, cuerpos,
                                     estado_voyager, integrar, lunas)
from asistencia.kepler import elementos, energia_especifica, parametro_tisserand
from asistencia.voyager import cobertura, estado_relativo_jupiter, instante_de_minimo

# Tabla de la misión (NSSDC/JPL; la reproduce el PDS Ring-Moon Systems Node): hora UTC y km
PUBLICADO = {"callisto": ("1979-07-08 12:21", 214930), "ganymede": ("1979-07-09 07:14", 62130),
             "europa": ("1979-07-09 17:53", 205720), "io": ("1979-07-09 23:17", 1129900)}
UMBRAL = 1e-7          # km/s^2 (0.1 mm/s^2): por encima se considera empuje


def vinf_eq(y, jd):
    """sqrt(2 eps_J): energía respecto al baricentro de Júpiter como v_inf equivalente."""
    (_, _), (rj, vj), _ = cuerpos(jd)
    return np.sqrt(2 * energia_especifica(y[:3] - rj, y[3:] - vj, GM_JUPITER))


def error(y, jd):
    return float(np.linalg.norm(y[:3] - estado_voyager(jd)[:3]))


def main():
    a, b = cobertura()
    t0, t1 = a + 0.01, b - 0.01
    t_ca, d_ca = instante_de_minimo(lambda t: np.linalg.norm(estado_relativo_jupiter(t)[0]), t0, t1, 8000)

    # ---------------- A. Sobrevuelos de los satélites ----------------
    def dist_luna(n):
        def d(t):
            sat = lunas(t)
            return np.linalg.norm(estado_relativo_jupiter(t)[0] - (centro_jupiter(t, sat)[0] + sat[n][0]))
        return d
    sobrevuelos = {}
    for n, (fecha_pub, d_pub) in PUBLICADO.items():
        t, d = instante_de_minimo(dist_luna(n), t_ca - 2.5, t_ca + 1.0)
        sobrevuelos[NOMBRES[n]] = {"fecha_utc": jd_a_fecha(t, segundos=True), "distancia_km": round(float(d)),
                                   "publicado_utc": fecha_pub, "publicado_km": d_pub,
                                   "diferencia_%": round(100 * (d / d_pub - 1), 3)}
    t_centro, d_centro = instante_de_minimo(
        lambda t: np.linalg.norm(estado_relativo_jupiter(t)[0] - centro_jupiter(t)[0]), t_ca - 0.2, t_ca + 0.2, 2000)

    # ---------------- B. Aceleración residual y maniobra ----------------
    def exceso(minutos):
        return np.linalg.norm(aceleracion_residual(t_ca + minutos / 1440)) - UMBRAL
    mins = np.arange(-1500.0, 900.0, 1.0)      # incluye los sobrevuelos de Ganímedes y Europa
    res_fino = np.array([np.linalg.norm(aceleracion_residual(t_ca + m / 1440)) for m in mins])
    activo = mins[res_fino > UMBRAL]
    m_ini = brentq(exceso, activo[0] - 1, activo[0], xtol=1e-3)
    m_fin = brentq(exceso, activo[-1], activo[-1] + 1, xtol=1e-3)
    jd_ini, jd_fin = t_ca + m_ini / 1440, t_ca + m_fin / 1440
    tm = np.linspace(jd_ini, jd_fin, 153)[1:-1]
    acc = np.array([aceleracion_residual(t) for t in tm])
    a_media = acc.mean(axis=0)
    dv = a_media * (jd_fin - jd_ini) * DIA
    maniobra = (jd_ini, jd_fin, a_media)
    rho_m, vrho_m = estado_relativo_jupiter(0.5 * (jd_ini + jd_fin))
    fuera = res_fino[(mins < m_ini - 2) | (mins > m_fin + 2)]
    m_pico = mins[(mins < m_ini - 2) | (mins > m_fin + 2)][int(np.argmax(fuera))]

    # ---------------- C. Integración con cada modelo ----------------
    y0 = estado_voyager(t0)
    ts = np.unique(np.concatenate([np.linspace(0, (t1 - t0) * DIA, 800),
                                   (np.linspace(t_ca - 2.5, t_ca + 1.5, 1600) - t0) * DIA]))
    jds = t0 + ts / DIA
    i_ca = int(np.argmin(np.abs(jds - t_ca)))
    real = np.array([estado_voyager(t) for t in jds]).T
    series = {"real": np.array([vinf_eq(real[:, k], jds[k]) for k in range(len(jds))])}
    modelos = {}
    for nombre, modelo, man in [("puntual", "puntual", None), ("satelites", "satelites", None),
                                ("satelites+J2", "satelites+J2", None),
                                ("satelites+J2+maniobra", "satelites+J2", maniobra)]:
        _, y = integrar(y0, t0, t1, modelo, man, t_eval=ts)
        err = np.linalg.norm(y[:3] - real[:3], axis=0)
        series[nombre] = np.array([vinf_eq(y[:, k], jds[k]) for k in range(len(jds))])
        modelos[nombre] = {"error_en_periapsis_km": round(float(err[i_ca])),
                           "error_final_km": round(float(err[-1])),
                           "vinf_eq_final_km_s": round(float(series[nombre][-1]), 4)}
    # tramos balísticos: hasta el inicio de la maniobra, y desde su fin partiendo del estado real
    tramos = {}
    for modelo in ("puntual", "satelites+J2"):
        _, y = integrar(y0, t0, jd_ini, modelo)
        _, y2 = integrar(estado_voyager(jd_fin), jd_fin, t1, modelo)
        tramos[modelo] = {"tramo1_4jun_a_maniobra_error_km": round(error(y[:, -1], jd_ini)),
                          "tramo2_maniobra_a_20jul_error_km": round(error(y2[:, -1], t1))}

    # ---------------- D. Caída de la energía jovicéntrica ----------------
    v_ini, v_fin = series["real"][0], series["real"][-1]
    salto = {"vinf_eq_real_inicio_km_s": round(float(v_ini), 4), "vinf_eq_real_fin_km_s": round(float(v_fin), 4),
             "caida_total_km_s": round(float(v_fin - v_ini), 4),
             "marea_solar_km_s": round(float(series["puntual"][-1] - v_ini), 4),
             "satelites_y_J2_km_s": round(float(series["satelites+J2"][-1] - series["puntual"][-1]), 4),
             "maniobra_km_s": round(float(series["satelites+J2+maniobra"][-1] - series["satelites+J2"][-1]), 4),
             "sin_explicar_km_s": round(float(v_fin - series["satelites+J2+maniobra"][-1]), 4)}

    # ---------------- E. Prolongación fuera de la esfera de influencia ----------------
    def resumen(y, jd):
        (rs, vs), (rj, _), _ = cuerpos(jd)
        r, v = y[:3] - rs, y[3:] - vs
        el = elementos(r, v, GM_SOL)
        return {"fecha": jd_a_fecha(jd), "distancia_a_jupiter_rSOI": round(float(np.linalg.norm(y[:3] - rj) / R_SOI_JUPITER), 2),
                "energia_heliocentrica_km2_s2": round(float(el["E"]), 2), "excentricidad": round(float(el["e"]), 3),
                "semieje_UA": round(float(el["a"] / AU), 2),
                "tisserand": round(float(parametro_tisserand(r, v, GM_SOL, A_JUPITER_KM)), 4)}
    dias = 200
    _, ya = integrar(y0, t0, t0 - dias, "satelites+J2")
    _, yb = integrar(estado_voyager(t1), t1, t1 + dias, "satelites+J2")
    asint = {"antes": resumen(ya[:, -1], t0 - dias), "inicio_registro": resumen(y0, t0),
             "fin_registro": resumen(estado_voyager(t1), t1), "despues": resumen(yb[:, -1], t1 + dias)}
    asint["ganancia_de_energia_asintotica_km2_s2"] = round(
        asint["despues"]["energia_heliocentrica_km2_s2"] - asint["antes"]["energia_heliocentrica_km2_s2"], 1)
    asint["cambio_relativo_tisserand"] = round(asint["despues"]["tisserand"] / asint["antes"]["tisserand"] - 1, 4)

    out = {
        "A_sobrevuelos_satelites": sobrevuelos,
        "A_maximo_acercamiento_a_jupiter": {
            "baricentro_del_sistema": {"fecha_utc": jd_a_fecha(t_ca, segundos=True), "distancia_km": round(float(d_ca))},
            "centro_del_planeta": {"fecha_utc": jd_a_fecha(t_centro, segundos=True), "distancia_km": round(float(d_centro))},
            "publicado": {"fecha_utc": "1979-07-09 22:29", "distancia_km": 721670}},
        "B_maniobra": {
            "inicio_min_tras_periapsis": round(float(m_ini), 1), "fin_min_tras_periapsis": round(float(m_fin), 1),
            "inicio_utc": jd_a_fecha(jd_ini, segundos=True), "fin_utc": jd_a_fecha(jd_fin, segundos=True),
            "duracion_min": round(float(m_fin - m_ini), 1),
            "aceleracion_mm_s2": round(float(np.linalg.norm(a_media) * 1e6), 3),
            "delta_v_m_s": round(float(np.linalg.norm(dv) * 1e3), 2),
            "angulo_con_velocidad_jovicentrica_grados": round(float(np.degrees(np.arccos(
                dv @ vrho_m / np.linalg.norm(dv) / np.linalg.norm(vrho_m)))), 1),
            "distancia_a_jupiter_RJ": round(float(np.linalg.norm(rho_m) / R_JUPITER), 2),
            "rapidez_jovicentrica_km_s": round(float(np.linalg.norm(vrho_m)), 2),
            "cambio_energia_jovicentrica_km2_s2": round(float(vrho_m @ dv), 3),
            "residual_fuera_de_la_maniobra_mm_s2": {"mediana": float(f"{np.median(fuera) * 1e6:.2g}"),
                                                      "maximo": float(f"{fuera.max() * 1e6:.2g}"),
                                                      "maximo_en_min_tras_periapsis": float(m_pico)},
            "publicado": "encendido de 76 minutos, unas dos horas después del máximo acercamiento "
                         "(Morrison y Samz, Voyage to Jupiter, NASA SP-439, 1980)"},
        "C_integracion_desde_4_jun_1979": modelos,
        "C_tramos_balisticos": tramos,
        "D_caida_energia_jovicentrica": salto,
        "E_prolongacion_200_dias": asint,
    }
    (RES / "voyager2_satelites.json").write_text(json.dumps(out, indent=2, ensure_ascii=False))
    print(json.dumps(out, indent=2, ensure_ascii=False))

    # ---------------- Figura ----------------
    horas = np.arange(-24.0, 12.01, 2 / 60)
    res = {m: np.array([np.linalg.norm(aceleracion_residual(t_ca + h / 24, m)) for h in horas]) * 1e6
           for m in ("puntual", "satelites+J2")}
    fig, ax = figura(1, 2, figsize=(12, 4.6))
    for mod, col, lab, xy, ha in [("puntual", COL["alterno"], "Júpiter puntual", (-22.0, 0.018), "left"),
                                  ("satelites+J2", COL["sonda"], "con satélites y $J_2$", (1.2, 1.4e-3), "right")]:
        ax[0].semilogy(horas, res[mod], color=col, lw=2)
        ax[0].annotate(lab, xy, ha=ha, va="bottom", fontsize=9.5, color=col)
    ax[0].set_ylim(1e-4, 30)
    ax[0].set_xlabel("horas desde el máximo acercamiento")
    ax[0].set_ylabel("aceleración no explicada (mm/s²)")
    ax[0].annotate("maniobra\n76 min", (0.5 * (m_ini + m_fin) / 60, 4.2), ha="center", va="bottom", fontsize=9.5)
    ax[0].annotate("Ganímedes", (-15.3, 3.2), ha="center", va="bottom", fontsize=9.5, color=COL["gris"])
    ax[0].annotate("Europa", (-4.6, 0.17), ha="center", va="bottom", fontsize=9.5, color=COL["gris"])
    ax[0].set_title("(a) Aceleración real menos la del modelo", loc="left", fontsize=12)

    dias_ca = jds - t_ca
    m = (dias_ca > -2.2) & (dias_ca < 1.2)
    for nombre, col, lab in [("puntual", COL["alterno"], "Júpiter puntual"),
                             ("satelites+J2", COL["saturno"], "+ satélites y $J_2$"),
                             ("satelites+J2+maniobra", COL["sonda"], "+ maniobra")]:
        ax[1].plot(dias_ca[m], series[nombre][m], color=col, lw=2)
        ax[1].annotate(lab, (1.2, series[nombre][m][-1]), xytext=(0, 4), textcoords="offset points",
                       ha="right", va="bottom", fontsize=9.5, color=col)
    ax[1].plot(dias_ca[m], series["real"][m], color=COL["tinta"], lw=1.6, ls=(0, (4, 3)), label="Voyager 2 real")
    ax[1].set_xlabel("días desde el máximo acercamiento")
    ax[1].set_ylabel(r"$\sqrt{2\varepsilon_J}$ (km/s)")
    ax[1].legend(fontsize=9.5, loc="lower left")
    ax[1].set_title("(b) Energía respecto a Júpiter", loc="left", fontsize=12)
    titulo(fig, "Voyager 2: qué le faltaba al Modelo 2", x=0.01, ha="left", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "fig09_voyager2_satelites.png")


if __name__ == "__main__":
    main()

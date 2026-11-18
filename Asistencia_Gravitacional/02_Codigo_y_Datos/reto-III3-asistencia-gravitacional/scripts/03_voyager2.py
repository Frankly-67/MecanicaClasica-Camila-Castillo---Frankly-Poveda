"""Validación con datos reales: el sobrevuelo de Voyager 2 por Júpiter (julio de 1979).

Datos: kernel SPICE vgr2_jup230.bsp (trayectoria reconstruida) + efeméride DE421.
Nada de esto se usó para ajustar el Modelo 1.

Calcula:
  1. Máximo acercamiento y elementos jovicéntricos osculantes en el periapsis.
  2. Qué tan constante es la energía jovicéntrica real (¿|v_inf| entra = sale?).
  3. Error de la hipérbola kepleriana (cónica) frente a la trayectoria real.
  4. Cambio de energía y de momento angular respecto al Sol: relación de Jacobi
     Delta E = n Delta l_z con el n medio de Júpiter y con su velocidad angular instantánea.
  5. Predicción del Modelo 1 en las fechas reales de Voyager 2 (circular y DE421).

(El parámetro de Tisserand y los valores asintóticos están en 06_voyager2_satelites.py,
porque hay que prolongar la trayectoria fuera de la esfera de influencia.)

Salidas: resultados/voyager2.json, figuras/fig05_voyager2_real.png
"""
import json

import numpy as np
from estilo import COL, FIG, RES, figura, titulo

from asistencia.constantes import A_JUPITER_KM, DIA, GM_JUPITER, GM_SOL, R_JUPITER, R_SOI_JUPITER
from asistencia.efemerides import estado_planeta, fecha_a_jd, jd_a_fecha
from asistencia.kepler import elementos, energia_especifica, propagar
from asistencia.modelo1 import evaluar_mision
from asistencia.sobrevuelo import angulo_entre, angulo_giro
from asistencia.voyager import (cobertura, estado_heliocentrico, estado_relativo_jupiter,
                                instante_de_minimo)

V2 = ("1977-08-20T14:29", "1979-07-09T22:29", "1981-08-26T03:24")


def main():
    a, b = cobertura()
    t_a, t_b = a + 0.01, b - 0.01
    ts = np.linspace(t_a, t_b, 4000)
    d = np.array([np.linalg.norm(estado_relativo_jupiter(t)[0]) for t in ts])
    t_ca, _ = instante_de_minimo(lambda t: np.linalg.norm(estado_relativo_jupiter(t)[0]), t_a, t_b, 8000)
    r_ca, v_ca = estado_relativo_jupiter(t_ca)

    # 1. Hipérbola osculante en el periapsis (lo que usaría el modelo de cónicas)
    el = elementos(r_ca, v_ca, GM_JUPITER)
    vinf_ca = np.sqrt(2 * el["E"])
    delta_formula = angulo_giro(vinf_ca, np.linalg.norm(r_ca), GM_JUPITER)

    # 2. Energía jovicéntrica a lo largo del registro, expresada como v_inf equivalente
    eps = np.array([energia_especifica(*estado_relativo_jupiter(t), GM_JUPITER) for t in ts])
    vinf_eq = np.sqrt(2 * eps)

    # 3. Hipérbola kepleriana propagada desde el periapsis frente a los datos reales
    err = np.array([np.linalg.norm(propagar(r_ca, v_ca, (t - t_ca) * DIA, GM_JUPITER)[0]
                                   - estado_relativo_jupiter(t)[0]) for t in ts])
    v_a = estado_relativo_jupiter(t_a)[1]
    v_b = estado_relativo_jupiter(t_b)[1]
    va_k = propagar(r_ca, v_ca, (t_a - t_ca) * DIA, GM_JUPITER)[1]
    vb_k = propagar(r_ca, v_ca, (t_b - t_ca) * DIA, GM_JUPITER)[1]

    # 4. Heliocéntrico
    rh_a, vh_a = estado_heliocentrico(t_a)
    rh_b, vh_b = estado_heliocentrico(t_b)
    vh = np.array([np.linalg.norm(estado_heliocentrico(t)[1]) for t in ts])

    # 4b. Relación de Jacobi Delta E = n Delta l_z con datos reales.
    # E es la energía TOTAL de la sonda (cinética + potencial del Sol + potencial de Júpiter)
    # y l_z su momento angular a lo largo del polo de la órbita de Júpiter, ambos respecto al
    # baricentro Sol-Júpiter, como en la formulación. La órbita real de Júpiter no es circular
    # (e = 0.048): se compara con el n medio y con la velocidad angular instantánea
    # w = |r_J x v_J| / r_J^2.
    al = GM_JUPITER / (GM_SOL + GM_JUPITER)
    hJ = np.cross(*estado_planeta("Jupiter", t_ca, "de421"))
    k = hJ / np.linalg.norm(hJ)
    n_J = np.sqrt((GM_SOL + GM_JUPITER) / A_JUPITER_KM**3)

    def jacobi(t):
        rJ, vJ = estado_planeta("Jupiter", t, "de421")
        r, v = estado_heliocentrico(t)
        rb, vb = r - al * rJ, v - al * vJ               # respecto al baricentro Sol-Júpiter
        E = 0.5 * vb @ vb - GM_SOL / np.linalg.norm(r) - GM_JUPITER / np.linalg.norm(r - rJ)
        # potencia del término radial: (dR_J/dt) (a_J . R_J^), con a_J la aceleración que
        # Júpiter da a la sonda y R_J la posición de Júpiter respecto al baricentro
        Rh = rJ / np.linalg.norm(rJ)
        a_J = -GM_JUPITER * (r - rJ) / np.linalg.norm(r - rJ) ** 3
        return (E, np.cross(rb, vb) @ k, np.linalg.norm(np.cross(rJ, vJ)) / (rJ @ rJ),
                (1 - al) * (vJ @ Rh) * (a_J @ Rh))
    tj = np.unique(np.concatenate([ts, np.linspace(t_ca - 1.5, t_ca + 1.5, 6000)]))
    Ej, lzj, wj, prad = np.array([jacobi(t) for t in tj]).T
    radial = float(np.sum(0.5 * (prad[1:] + prad[:-1]) * np.diff(tj) * DIA))
    dE_total = Ej[-1] - Ej[0]
    n_dlz = n_J * (lzj[-1] - lzj[0])
    w_dlz = float(np.sum(0.5 * (wj[1:] + wj[:-1]) * np.diff(lzj)))      # integral de w dl_z
    dE_kepler = energia_especifica(rh_b, vh_b, GM_SOL) - energia_especifica(rh_a, vh_a, GM_SOL)

    # 5. Modelo 1 en las fechas reales
    j0, jJ, jS = (fecha_a_jd(x) for x in V2)
    pred = {f: evaluar_mision(j0, jJ - j0, jS - jJ, f) for f in ("circular", "de421")}

    def km_s(x):
        return round(float(x), 3)
    out = {
        "fuente": "kernel SPICE vgr2_jup230.bsp (NAIF) + DE421",
        "cobertura": [jd_a_fecha(t_a, True), jd_a_fecha(t_b, True)],
        "maximo_acercamiento_al_baricentro_joviano": {
            "fecha_utc": jd_a_fecha(t_ca, segundos=True),
            "distancia_km": round(float(np.linalg.norm(r_ca))),
            "distancia_RJ": round(float(np.linalg.norm(r_ca) / R_JUPITER), 3),
            "publicado_PDS_respecto_al_centro_del_planeta": "1979-07-09 22:29 UTC, 721670 km",
            "nota": "respecto al centro del planeta: resultados/voyager2_satelites.json"},
        "hiperbola_osculante_en_periapsis": {
            "e": round(float(el["e"]), 4), "vinf_km_s": km_s(vinf_ca),
            "delta_formula_grados": round(float(np.degrees(delta_formula)), 2),
            "inclinacion_respecto_ecliptica_grados": round(float(np.degrees(el["i"])), 2)},
        "energia_jovicentrica_real": {
            "vinf_equivalente_inicio_km_s": km_s(vinf_eq[0]),
            "vinf_equivalente_periapsis_km_s": km_s(vinf_ca),
            "vinf_equivalente_fin_km_s": km_s(vinf_eq[-1]),
            "distancia_inicio_en_rSOI": round(float(d[0] / R_SOI_JUPITER), 3),
            "distancia_fin_en_rSOI": round(float(d[-1] / R_SOI_JUPITER), 3)},
        "error_conica_vs_real": {
            "error_posicion_inicio_km": round(float(err[0])),
            "error_posicion_fin_km": round(float(err[-1])),
            "error_posicion_inicio_RJ": round(float(err[0] / R_JUPITER), 1),
            "error_posicion_fin_RJ": round(float(err[-1] / R_JUPITER), 1),
            "giro_real_entre_extremos_grados": round(float(np.degrees(angulo_entre(v_a, v_b))), 3),
            "giro_conica_entre_extremos_grados": round(float(np.degrees(angulo_entre(va_k, vb_k))), 3)},
        "heliocentrico_en_los_extremos_del_registro": {
            "nota": "ambos extremos están dentro de la esfera de influencia; los valores "
                    "asintóticos están en voyager2_satelites.json",
            "rapidez_inicio_km_s": km_s(np.linalg.norm(vh_a)),
            "rapidez_fin_km_s": km_s(np.linalg.norm(vh_b)),
            "energia_solo_Sol_inicio_km2_s2": round(float(energia_especifica(rh_a, vh_a, GM_SOL)), 2),
            "energia_solo_Sol_fin_km2_s2": round(float(energia_especifica(rh_b, vh_b, GM_SOL)), 2),
            "energia_total_inicio_km2_s2": round(float(Ej[0]), 2),
            "energia_total_fin_km2_s2": round(float(Ej[-1]), 2)},
        "relacion_de_jacobi": {
            "Delta_E_total_km2_s2": round(float(dE_total), 2),
            "n_medio_por_Delta_lz_km2_s2": round(float(n_dlz), 2),
            "integral_w_dlz_km2_s2": round(w_dlz, 2),
            "w_sobre_n_en_periapsis": round(float(wj[int(np.argmin(np.abs(tj - t_ca)))] / n_J), 4),
            "distancia_Sol_Jupiter_UA": round(float(np.linalg.norm(estado_planeta("Jupiter", t_ca, "de421")[0]) / 1.495978707e8), 3),
            "diferencia_relativa_con_n_medio": round(float(n_dlz / dE_total - 1), 4),
            "diferencia_relativa_con_w_instantanea": round(float(w_dlz / dE_total - 1), 4),
            "termino_radial_km2_s2": round(radial, 2),
            "resto_fuerzas_no_incluidas_km2_s2": round(float(dE_total - w_dlz - radial), 2),
            "nota_resto": "trabajo de lo que no está en el potencial Sol + Júpiter puntual: "
                          "satélites, J2 y maniobra (voyager2_satelites.json)",
            "Delta_E_sin_potencial_de_Jupiter_km2_s2": round(float(dE_kepler), 2)},
        "modelo1_en_fechas_reales": {
            f: {"vinf_entrada_km_s": km_s(m["vinf_in"]), "vinf_salida_km_s": km_s(m["vinf_out"]),
                "desajuste_km_s": km_s(m["desajuste"]), "delta_grados": round(float(np.degrees(m["delta"])), 2),
                "rp_RJ": round(float(m["rp_RJ"]), 2), "C3_km2_s2": round(float(m["C3"]), 1)}
            for f, m in pred.items()},
    }
    (RES / "voyager2.json").write_text(json.dumps(out, indent=2, ensure_ascii=False))
    print(json.dumps(out, indent=2, ensure_ascii=False))

    # ---------------- Figura ----------------
    dias = (ts - t_ca)
    fig, ax = figura(1, 3, figsize=(15, 4.8))
    # (a) trayectoria real vs cónica, sistema de Júpiter (plano eclíptico)
    rr = np.array([estado_relativo_jupiter(t)[0] for t in ts]) / R_JUPITER
    rk = np.array([propagar(r_ca, v_ca, (t - t_ca) * DIA, GM_JUPITER)[0] for t in ts]) / R_JUPITER
    ax[0].plot(rk[:, 0], rk[:, 1], color=COL["alterno"], lw=1.8, ls="--", label="Hipérbola kepleriana")
    ax[0].plot(rr[:, 0], rr[:, 1], color=COL["sonda"], lw=2, label="Voyager 2 real (SPICE)")
    ax[0].annotate("entrada (4 jun)", rr[0, :2], xytext=(-95, 0), textcoords="offset points", fontsize=9)
    ax[0].annotate("salida (20 jul)", rr[-1, :2], xytext=(-20, 10), textcoords="offset points", fontsize=9)
    ax[0].plot(0, 0, "o", color=COL["jupiter"], ms=6)
    ax[0].set_aspect("equal")
    ax[0].set_xlim(-200, 120)
    ax[0].set_ylim(-410, 60)
    ax[0].set_xlabel("x eclíptico ($R_J$)")
    ax[0].set_ylabel("y eclíptico ($R_J$)")
    ax[0].set_title("(a) Sistema de Júpiter", loc="left", fontsize=12)
    ax[0].legend(fontsize=9, loc="lower left", title="(coinciden a esta escala)", title_fontsize=8)
    # (b) energía jovicéntrica como v_inf equivalente
    ax[1].plot(dias, vinf_eq, color=COL["sonda"], lw=2)
    ax[1].axhline(vinf_ca, color=COL["gris"], ls="--", lw=1)
    ax[1].set_xlabel("días desde el máximo acercamiento")
    ax[1].set_ylabel(r"$\sqrt{2\varepsilon_J}$ (km/s)")
    ax[1].set_title("(b) Energía respecto a Júpiter", loc="left", fontsize=12)
    # (c) rapidez heliocéntrica
    ax[2].plot(dias, vh, color=COL["sonda"], lw=2)
    ax[2].axhline(np.sqrt(GM_SOL / np.linalg.norm(rh_b) * 2), color=COL["gris"], ls=":", lw=1)
    ax[2].annotate("escape del Sol", (dias[-1], np.sqrt(2 * GM_SOL / np.linalg.norm(rh_b))),
                   xytext=(-80, 4), textcoords="offset points", fontsize=9, color=COL["gris"])
    ax[2].set_xlabel("días desde el máximo acercamiento")
    ax[2].set_ylabel("rapidez respecto al Sol (km/s)")
    ax[2].set_title("(c) Rapidez heliocéntrica", loc="left", fontsize=12)
    titulo(fig, "Voyager 2 en Júpiter: datos reales", x=0.01, ha="left", fontsize=15, fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "fig05_voyager2_real.png")


if __name__ == "__main__":
    main()

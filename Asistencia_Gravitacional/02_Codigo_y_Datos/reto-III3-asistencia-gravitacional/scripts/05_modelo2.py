"""Modelo 2 (N cuerpos, REBOUND/IAS15) y comparación con el Modelo 1 (cónicas empalmadas).

Partes:
  A. Invariante: la constante de Jacobi h = E - n l_z se conserva en la simulación
     mientras la energía heliocéntrica cambia.
  B. Convergencia numérica: tolerancia de IAS15 y un integrador independiente (DOP853).
  C. Error de las cónicas: barrido de v_inf, radio de empalme R_e y periapsis r_p,
     partiendo del mismo periapsis en ambos modelos.
  D. Efecto de añadir Saturno.
  E. Voyager 2: integración desde el estado real (SPICE) con Sol, Júpiter y Saturno de DE421.

Salidas: resultados/modelo2.json, figuras/fig06_jacobi_convergencia.png,
         figuras/fig07_error_conicas.png, figuras/fig08_voyager2_modelo2.png
"""
import json

import numpy as np
import rebound
from estilo import COL, FIG, RES, figura, nota, titulo
from scipy.integrate import solve_ivp

from asistencia.constantes import DIA, GM_JUPITER, GM_SOL, R_HILL_JUPITER, R_JUPITER, R_SOI_JUPITER
from asistencia.efemerides import estado_planeta, jd_a_fecha
from asistencia.kepler import energia_especifica, propagar
from asistencia.modelo2 import (MU1, MU2, N_J, conicas_empalmadas, estado_periapsis_plano,
                                integrar, jupiter_heliocentrico)
from asistencia.voyager import cobertura, estado_heliocentrico, estado_relativo_jupiter

GM_SATURNO = 8.459715185680659e-08 * 149597870.691**3 / 86400.0**2   # DE421


def h_exacto(r_h, v_h, t):
    """Jacobi h = E - n l_z exacta del problema restringido, a partir del estado
    heliocéntrico en el instante t. Se pasa al sistema baricéntrico con la órbita circular
    de Júpiter: el Sol está en -alpha r_J y se mueve con -alpha v_J."""
    al = MU2 / (MU1 + MU2)
    rJ, vJ = jupiter_heliocentrico(t)
    r_b, v_b = r_h - al * rJ, v_h - al * vJ
    E = 0.5 * v_b @ v_b - MU1 / np.linalg.norm(r_h) - MU2 / np.linalg.norm(r_h - rJ)
    return E - N_J * np.cross(r_b, v_b)[2]


def comparar(vinf, rp_RJ, R_e, eps=1e-9, saturno=False):
    r_p, v_p = estado_periapsis_plano(vinf, rp_RJ * R_JUPITER)
    T = 3 * R_SOI_JUPITER / vinf + 30 * DIA            # lejos de Júpiter (> 3 r_SOI)
    T = max(T, 1.3 * abs(conicas_empalmadas(r_p, v_p, R_e, 1.0)[2]))
    rN, vN, _ = integrar(r_p, v_p, T, eps, saturno=saturno)
    rNb, vNb, _ = integrar(r_p, v_p, -T, eps, saturno=saturno)
    rP, vP, _ = conicas_empalmadas(r_p, v_p, R_e, T)
    rPb, vPb, _ = conicas_empalmadas(r_p, v_p, R_e, -T)
    dEN = energia_especifica(rN, vN, MU1) - energia_especifica(rNb, vNb, MU1)
    dEP = energia_especifica(rP, vP, MU1) - energia_especifica(rPb, vPb, MU1)
    return {"vinf": vinf, "rp_RJ": rp_RJ, "R_e_rSOI": R_e / R_SOI_JUPITER, "T_dias": T / DIA,
            "dE_nbody": dEN, "dE_conicas": dEP, "error_dE_rel": (dEP - dEN) / dEN,
            "error_v_salida_km_s": float(np.linalg.norm(vP - vN)),
            "error_r_salida_RJ": float(np.linalg.norm(rP - rN) / R_JUPITER),
            "jacobi_violada_conicas_rel": float((h_exacto(rP, vP, T) - h_exacto(rPb, vPb, -T)) / abs(h_exacto(rPb, vPb, -T))),
            "jacobi_nbody_rel": float((h_exacto(rN, vN, T) - h_exacto(rNb, vNb, -T)) / abs(h_exacto(rNb, vNb, -T)))}


def parte_a():
    r_p, v_p = estado_periapsis_plano(7.6, 10 * R_JUPITER)
    T = 250 * DIA
    _, _, s_f = integrar(r_p, v_p, T, muestras=400)
    _, _, s_b = integrar(r_p, v_p, -T, muestras=400)
    s = np.vstack([s_b[::-1], s_f[1:]])
    # energía heliocéntrica a lo largo de la misma integración
    E = []
    for t in s[:, 0]:
        r, v, _ = integrar(r_p, v_p, t) if t != 0 else (None, None, None)
        E.append(energia_especifica(r, v, MU1) if t != 0 else np.nan)
    return s, np.array(E)


def parte_b():
    r_p, v_p = estado_periapsis_plano(7.6, 10 * R_JUPITER)
    T = 250 * DIA
    ref, _, _ = integrar(r_p, v_p, T, 1e-12)
    res = {}
    for eps in [1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11]:
        r, _, _ = integrar(r_p, v_p, T, eps)
        res[f"{eps:.0e}"] = float(np.linalg.norm(r - ref))

    # Integrador independiente: DOP853 de scipy, con Sol y Júpiter en órbitas circulares
    # prescritas alrededor del baricentro (sistema inercial).
    al = MU2 / (MU1 + MU2)
    a = np.linalg.norm(jupiter_heliocentrico(0)[0])

    def pos(t):
        u = np.array([np.cos(N_J * t), np.sin(N_J * t), 0.0])
        return -al * a * u, (1 - al) * a * u

    def f(t, y):
        P1, P2 = pos(t)
        d1, d2 = y[:3] - P1, y[:3] - P2
        return np.concatenate([y[3:], -MU1 * d1 / np.linalg.norm(d1) ** 3 - MU2 * d2 / np.linalg.norm(d2) ** 3])
    P1, P2 = pos(0)
    vJbar = (1 - al) * a * N_J * np.array([0, 1.0, 0])
    y0 = np.concatenate([P2 + r_p, vJbar + v_p])
    sol = solve_ivp(f, (0, T), y0, method="DOP853", rtol=1e-13, atol=1e-6)
    P1T, _ = pos(T)
    r_dop = sol.y[:3, -1] - P1T
    return res, float(np.linalg.norm(r_dop - ref))


def parte_c():
    vinfs = [2.5, 3, 4, 5, 6, 7, 8, 10, 12, 15]
    radios = {"0.25 r_SOI": 0.25, "0.5 r_SOI": 0.5, "r_SOI (Laplace)": 1.0,
              "r_Hill": R_HILL_JUPITER / R_SOI_JUPITER, "2 r_SOI": 2.0, "3 r_SOI": 3.0}
    barrido_R = {k: [comparar(v, 10, f * R_SOI_JUPITER) for v in vinfs] for k, f in radios.items()}
    barrido_rp = {f"r_p = {rp} R_J": [comparar(v, rp, R_SOI_JUPITER) for v in vinfs] for rp in (3, 10, 30)}
    return vinfs, barrido_R, barrido_rp


def parte_d():
    out = {}
    for v in (5, 7.6, 12):
        sin = comparar(v, 10, R_SOI_JUPITER)
        con = comparar(v, 10, R_SOI_JUPITER, saturno=True)
        out[f"{v}"] = {"dE_sin_saturno": sin["dE_nbody"], "dE_con_saturno": con["dE_nbody"],
                       "cambio_rel": (con["dE_nbody"] - sin["dE_nbody"]) / sin["dE_nbody"]}
    return out


def parte_e():
    a, b = cobertura()
    t0, t1 = a + 0.01, b - 0.01
    sim = rebound.Simulation()
    sim.G = 1.0
    sim.integrator = "ias15"
    try:
        sim.integrator.epsilon = 1e-10
    except AttributeError:                 # REBOUND 4
        sim.ri_ias15.epsilon = 1e-10
    sim.add(m=GM_SOL)
    for n, gm in (("Jupiter", GM_JUPITER), ("Saturno", GM_SATURNO)):
        r, v = estado_planeta(n, t0, "de421")
        sim.add(m=gm, x=r[0], y=r[1], z=r[2], vx=v[0], vy=v[1], vz=v[2])
    sim.N_active = sim.N
    r, v = estado_heliocentrico(t0)
    sim.add(m=0.0, x=r[0], y=r[1], z=r[2], vx=v[0], vy=v[1], vz=v[2])
    sim.move_to_com()
    # Modelo 1 desde el mismo estado inicial: hipérbola kepleriana alrededor de Júpiter
    rho0, vrho0 = estado_relativo_jupiter(t0)
    ts = np.linspace(t0, t1, 1500)
    filas = []
    for t in ts:
        sim.integrate((t - t0) * DIA, exact_finish_time=1)
        pj, ps = sim.particles[1], sim.particles[3]
        rho_n = np.array([ps.x - pj.x, ps.y - pj.y, ps.z - pj.z])
        vrho_n = np.array([ps.vx - pj.vx, ps.vy - pj.vy, ps.vz - pj.vz])
        rho_s, vrho_s = estado_relativo_jupiter(t)
        rho_k, _ = propagar(rho0, vrho0, (t - t0) * DIA, GM_JUPITER)
        filas.append([t - t0, np.linalg.norm(rho_n - rho_s), np.linalg.norm(rho_k - rho_s),
                      np.sqrt(2 * energia_especifica(rho_n, vrho_n, GM_JUPITER)),
                      np.sqrt(2 * energia_especifica(rho_s, vrho_s, GM_JUPITER)),
                      np.linalg.norm(rho_s)])
    f = np.array(filas)
    i_ca = int(np.argmin(f[:, 5]))
    return f, i_ca, jd_a_fecha(t0, True), jd_a_fecha(t1, True)


def main():
    print("A. Jacobi ...")
    sA, EA = parte_a()
    dh = np.abs(sA[:, 1] - sA[len(sA) // 2, 1]) / abs(sA[len(sA) // 2, 1])
    print("B. Convergencia ...")
    conv, err_dop = parte_b()
    print("C. Barrido ...")
    vinfs, bR, brp = parte_c()
    print("D. Saturno ...")
    sat = parte_d()
    print("E. Voyager 2 ...")
    fv, i_ca, ini, fin = parte_e()

    def r3(x):
        return float(f"{x:.3g}")
    tabla_R = {k: {f"{d['vinf']:g}": {"error_dE_%": r3(100 * d["error_dE_rel"]),
                                      "jacobi_violada_conicas": r3(d["jacobi_violada_conicas_rel"])}
                   for d in lst} for k, lst in bR.items()}
    tabla_rp = {k: {f"{d['vinf']:g}": r3(100 * d["error_dE_rel"]) for d in lst} for k, lst in brp.items()}
    out = {
        "A_jacobi": {"max_variacion_relativa_h": r3(dh.max()),
                     "E_helio_inicio": r3(EA[0]), "E_helio_fin": r3(EA[-1])},
        "B_convergencia": {"posicion_final_vs_eps1e-12_km": {k: r3(v) for k, v in conv.items()},
                           "DOP853_vs_IAS15_km": r3(err_dop)},
        "C_error_conicas_vs_vinf_y_Re": tabla_R,
        "C_error_conicas_vs_rp_(Re=rSOI)_en_%": tabla_rp,
        "C_jacobi_nbody_max": r3(max(abs(d["jacobi_nbody_rel"]) for lst in bR.values() for d in lst)),
        "D_efecto_saturno": {k: {kk: r3(vv) for kk, vv in v.items()} for k, v in sat.items()},
        "E_voyager2": {"intervalo": [ini, fin],
                       "error_nbody_en_periapsis_km": r3(fv[i_ca, 1]),
                       "error_nbody_final_km": r3(fv[-1, 1]),
                       "error_conica_en_periapsis_km": r3(fv[i_ca, 2]),
                       "error_conica_final_km": r3(fv[-1, 2]),
                       "vinf_eq_nbody_inicio_fin": [r3(fv[0, 3]), r3(fv[-1, 3])],
                       "vinf_eq_spice_inicio_fin": [r3(fv[0, 4]), r3(fv[-1, 4])]},
    }
    (RES / "modelo2.json").write_text(json.dumps(out, indent=2, ensure_ascii=False))
    print(json.dumps(out, indent=2, ensure_ascii=False))

    # ----- fig06: Jacobi y convergencia -----
    fig, ax = figura(1, 2, figsize=(12, 4.6))
    t = sA[:, 0] / DIA
    ax[0].plot(t, EA, color=COL["sonda"], lw=2, label="energía heliocéntrica $E$")
    ax[0].plot(t, sA[:, 1], color=COL["tierra"], lw=2, label=r"$h = E - n\,\ell_z$")
    ax[0].set_xlabel("días desde el periapsis")
    ax[0].set_ylabel("km²/s²")
    ax[0].legend(loc="center right", fontsize=9.5)
    ax[0].set_title("(a) $E$ cambia, $h$ no", loc="left", fontsize=12)
    ins = ax[0].inset_axes([0.08, 0.55, 0.38, 0.33])
    ins.semilogy(t, np.maximum(dh, 1e-17), color=COL["tierra"], lw=1)
    ins.set_title(r"$|\Delta h/h|$", fontsize=8)
    ins.tick_params(labelsize=7)
    eps = [float(k) for k in conv]
    ax[1].loglog(eps, list(conv.values()), "o-", color=COL["sonda"], lw=2, label="IAS15 vs IAS15 ($\\epsilon=10^{-12}$)")
    ax[1].axhline(err_dop, color=COL["alterno"], ls="--", lw=1.5, label="DOP853 (integrador independiente)")
    ax[1].set_xlabel(r"tolerancia de IAS15 $\epsilon$")
    ax[1].set_ylabel("diferencia en la posición final (km)")
    ax[1].legend(fontsize=9.5)
    ax[1].set_title("(b) Convergencia (250 días)", loc="left", fontsize=12)
    titulo(fig, "Modelo 2: invariante y convergencia", x=0.01, ha="left", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "fig06_jacobi_convergencia.png")

    # ----- fig07: error de las cónicas -----
    fig, ax = figura(1, 2, figsize=(12, 4.8))
    paleta = [COL["alterno"], COL["tierra"], COL["sonda"], COL["jupiter"], COL["saturno"], COL["tinta"]]
    def etiqueta(k):
        import re
        k = k.replace("r_SOI", r"r_{SOI}").replace("r_Hill", r"r_{Hill}").replace("R_J", r"R_J")
        m = re.match(r"([\d.]+) (r_\{SOI\})(.*)", k)
        if m:
            return f"${m.group(1)}\\,{m.group(2)}${m.group(3)}"
        if k.startswith("r_p"):
            return "$" + k.replace(" R_J", r"\,R_J") + "$"
        return k.replace("r_{SOI} ", "$r_{SOI}$ ").replace("r_{Hill}", "$r_{Hill}$")

    def curva(a, lst, c, k):
        k = etiqueta(k)
        y = np.array([100 * d["error_dE_rel"] for d in lst])
        a.loglog(vinfs, np.abs(y), "-", color=c, lw=1.8, label=k)
        a.loglog(np.array(vinfs)[y > 0], y[y > 0], "o", color=c, ms=5)
        a.loglog(np.array(vinfs)[y < 0], -y[y < 0], "o", mfc="white", mec=c, ms=5)
    for (k, lst), c in zip(bR.items(), paleta):
        curva(ax[0], lst, c, k)
    ax[0].set_xlabel(r"$v_\infty$ (km/s)")
    ax[0].set_ylabel(r"error en $\Delta E$ (%)")
    ax[0].set_title(r"(a) Según el radio de empalme ($r_p = 10\,R_J$)", loc="left", fontsize=12)
    ax[0].set_ylim(0.015, 150)
    ax[0].legend(fontsize=9, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.27), columnspacing=1.4,
                 handlelength=1.6)
    for (k, lst), c in zip(brp.items(), [COL["alterno"], COL["sonda"], COL["saturno"]]):
        curva(ax[1], lst, c, k)
    ax[1].set_xlabel(r"$v_\infty$ (km/s)")
    ax[1].set_title(r"(b) Según el periapsis ($R_e = r_{SOI}$)", loc="left", fontsize=12)
    nota(ax[1], 0.02, 0.03, "punto lleno: las cónicas sobreestiman la ganancia\npunto vacío: la subestiman",
         transform=ax[1].transAxes, fontsize=8.5, color=COL["gris"])
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
    for a in ax:
        a.xaxis.set_major_locator(FixedLocator([2.5, 3, 4, 5, 6, 8, 10, 15]))
        a.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        a.xaxis.set_minor_formatter(NullFormatter())
        a.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax[1].set_ylim(0.05, 40)
    ax[1].legend(fontsize=9, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.27), columnspacing=1.4,
                 handlelength=1.6)
    titulo(fig, "Error de las cónicas empalmadas frente a N cuerpos", x=0.01, ha="left", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "fig07_error_conicas.png")

    # ----- fig08: Voyager 2 -----
    fig, ax = figura(1, 2, figsize=(12, 4.6))
    dias = fv[:, 0] - fv[i_ca, 0]
    ax[0].semilogy(dias, fv[:, 2], color=COL["alterno"], lw=2, label="Modelo 1 (hipérbola)")
    ax[0].semilogy(dias, fv[:, 1], color=COL["sonda"], lw=2, label="Modelo 2 (Sol+Júpiter+Saturno)")
    ax[0].set_xlabel("días desde el máximo acercamiento")
    ax[0].set_ylabel("distancia a la trayectoria real (km)")
    ax[0].legend(fontsize=9.5)
    ax[0].set_title("(a) Error frente a los datos SPICE", loc="left", fontsize=12)
    ax[1].plot(dias, fv[:, 4], color=COL["tinta"], lw=2, label="Voyager 2 real")
    ax[1].plot(dias, fv[:, 3], color=COL["sonda"], lw=2, ls="--", label="Modelo 2")
    ax[1].set_xlabel("días desde el máximo acercamiento")
    ax[1].set_ylabel(r"$\sqrt{2\varepsilon_J}$ (km/s)")
    ax[1].legend(fontsize=9.5)
    ax[1].set_title("(b) Energía respecto a Júpiter", loc="left", fontsize=12)
    titulo(fig, "Voyager 2: Modelo 2 desde el estado real del 4 jun 1979", x=0.01, ha="left",
           fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "fig08_voyager2_modelo2.png")


if __name__ == "__main__":
    main()

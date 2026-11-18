"""Figuras del Modelo 1. Requiere haber corrido 01_buscar_ventanas.py con ambas fuentes.

fig01: mapa de ventanas de lanzamiento (efemérides DE421).
fig02: tiempo mínimo hasta Saturno según el C3 permitido, circular frente a DE421
       (mide el impacto de la hipótesis H1).
fig03: trayectoria heliocéntrica de la misión óptima con C3 <= 100 (DE421).
fig04: el sobrevuelo en el sistema de Júpiter y el triángulo de velocidades.
"""
import json

import numpy as np
from estilo import COL, FIG, RES, figura, nota, titulo
from matplotlib.patches import Circle

from asistencia.constantes import AU, DIA, GM_JUPITER, GM_SOL, R_JUPITER, R_SOI_JUPITER
from asistencia.efemerides import estado_planeta, fecha_a_jd
from asistencia.kepler import propagar
from asistencia.modelo1 import evaluar_mision
from asistencia.sobrevuelo import estado_periapsis

CASO = "C3<=100"
VOY = {"Voyager 2": ("1977-08-20T14:29", "1979-07-09T22:29"),
       "Voyager 1": ("1977-09-05T12:56", "1979-03-05T12:05")}
T_REF = fecha_a_jd("1977-01-01")


def anio(jd):
    return 1977 + (jd - T_REF) / 365.25


def fig_ventanas():
    z = np.load(RES / "busqueda_de421.npz")
    jd0, tof1, C3, T = z["jd0"], z["tof1"], z["C3"], z["T_min_C3_105"]
    fig, ax = figura(1, 2, figsize=(12, 5), gridspec_kw={"width_ratios": [1, 1.25]})
    for k, a in enumerate(ax):
        m = a.pcolormesh(anio(jd0), tof1, T.T, cmap="Blues_r", shading="auto",
                         vmin=np.nanmin(T), vmax=np.nanpercentile(T, 90))
        cs = a.contour(anio(jd0), tof1, C3.T, levels=[90, 100, 110], colors=COL["gris"], linewidths=1)
        a.grid(False)
        a.set_xlabel("fecha de salida" + (" (1977)" if k else " (año)"))
        if k == 0:
            a.set_xticks(range(1976, 1981))
        if k == 1:
            a.clabel(cs, fmt=lambda v: f"C3={v:g}", fontsize=9)
            a.set_xlim(1977.55, 1977.8)
            a.set_xticks([anio(fecha_a_jd(f)) for f in ("1977-08-01", "1977-09-01", "1977-10-01")])
            a.set_xticklabels(["1 ago", "1 sep", "1 oct"])
            a.set_ylim(450, 800)
            for (n, (s, j)), mk in zip(VOY.items(), ["o", "s"]):
                x, y = anio(fecha_a_jd(s)), fecha_a_jd(j) - fecha_a_jd(s)
                a.plot(x, y, mk, ms=9, mfc="white", mec=COL["tinta"], mew=2)
                a.annotate(n, (x, y), xytext=(8, 6), textcoords="offset points", fontsize=10)
    ax[0].set_ylabel("tiempo de vuelo a Júpiter (días)")
    ax[0].add_patch(__import__("matplotlib").patches.Rectangle((1977.55, 450), 0.25, 350, fill=False,
                                                                ec=COL["tinta"], lw=1.2))
    ax[0].set_title("(a) 1976–1980", loc="left", fontsize=12)
    ax[1].set_title("(b) ventana de 1977", loc="left", fontsize=12)
    cb = fig.colorbar(m, ax=ax[1], pad=0.02)
    cb.set_label("tiempo mínimo hasta Saturno (años)")
    titulo(fig, "Ventanas Tierra → Júpiter → Saturno (efemérides DE421, C3 ≤ 105 km²/s²)",
           x=0.01, ha="left", fontsize=14, fontweight="bold")
    nota(fig, 0.01, -0.03, "Blanco: sin sobrevuelo sin motor (|Δ|v∞|| < 0.3 km/s, r_p ≥ 2 R_J). "
         "Líneas: C3 de salida. Las fechas de las Voyager no se usaron en la búsqueda.",
         fontsize=9, color=COL["gris"])
    fig.tight_layout()
    fig.savefig(FIG / "fig01_ventanas.png")


def fig_pareto():
    fig, ax = figura(figsize=(8, 4.8))
    for f, col, lab in [("circular", COL["alterno"], "Órbitas circulares (H1)"),
                        ("de421", COL["sonda"], "Efemérides reales DE421")]:
        z = np.load(RES / f"busqueda_{f}.npz")
        ax.step(z["limites_C3"], z["frente_T"], where="post", color=col, lw=2, label=lab)
        ref = json.loads((RES / f"mision_optima_{f}.json").read_text())["fechas_historicas_evaluadas_con_modelo1"]
        for (n, v), mk in zip(ref.items(), ["o", "s"]):
            ax.plot(v["C3_km2_s2"], v["tiempo_total_anios"], mk, ms=8, mfc="white", mec=col, mew=2)
            if f == "de421":
                ax.annotate(n, (v["C3_km2_s2"], v["tiempo_total_anios"]), xytext=(8, 4),
                            textcoords="offset points", fontsize=9)
    ax.set_xlabel("C3 máximo permitido (km²/s²)")
    ax.set_ylabel("tiempo hasta Saturno (años)")
    ax.legend(loc="upper right")
    titulo(ax, "Impacto de la hipótesis de órbitas circulares", loc="left")
    fig.savefig(FIG / "fig02_pareto_C3_tiempo.png")


def cargar(caso=CASO):
    d = json.loads((RES / "mision_optima_de421.json").read_text())["minimo_tiempo_con_limite_C3"][caso]
    x = d["exacto"]
    return evaluar_mision(x["jd_salida_tdb"], x["tof1_dias"], x["tof2_dias"], "de421"), d


def fig_trayectoria(m, d):
    fig, ax = figura(figsize=(7.4, 7.4))
    for n, col in [("Tierra", COL["tierra"]), ("Jupiter", COL["jupiter"]), ("Saturno", COL["saturno"])]:
        jds = np.linspace(m["jd_salida"], m["jd_salida"] + {"Tierra": 366, "Jupiter": 4333, "Saturno": 10760}[n], 500)
        p = np.array([estado_planeta(n, t, "de421")[0] for t in jds]) / AU
        ax.plot(p[:, 0], p[:, 1], color=col, lw=1, alpha=0.5)
    for r0, v0, T, ls, lab in [(m["rE"], m["v_salida"], m["tof1"], "-", "tramo 1"),
                               (m["rJ"], m["v_sale_J"], m["tof2"], (0, (6, 3)), "tramo 2")]:
        p = np.array([propagar(r0, v0, t, GM_SOL)[0] for t in np.linspace(0, T * DIA, 400)]) / AU
        ax.plot(p[:, 0], p[:, 1], color=COL["sonda"], lw=2.2, ls=ls, label=f"Sonda, {lab}")
    ax.plot(0, 0, "o", ms=12, color="#E4A900")
    for txt, r, col in [("Salida " + d["salida"], m["rE"], COL["tierra"]),
                        ("Júpiter " + d["llegada_jupiter"], m["rJ"], COL["jupiter"]),
                        ("Saturno " + d["llegada_saturno"], m["rS"], COL["saturno"])]:
        ax.plot(r[0] / AU, r[1] / AU, "o", ms=9, color=col, mec="white")
        ax.annotate(txt, (r[0] / AU, r[1] / AU), xytext=(8, 6), textcoords="offset points", fontsize=9.5)
    ax.set_aspect("equal")
    ax.set_xlim(-10.5, 10.5)
    ax.set_ylim(-10.5, 10.5)
    ax.set_xlabel("x eclíptico (UA)")
    ax.set_ylabel("y eclíptico (UA)")
    ax.legend(loc="lower left", fontsize=9.5)
    titulo(ax, f"Misión óptima del Modelo 1 con DE421 ({CASO.replace('<=', ' ≤ ')})", loc="left")
    fig.savefig(FIG / "fig03_trayectoria_modelo1.png")


def fig_sobrevuelo(m):
    fig, (ax, bx) = figura(1, 2, figsize=(13, 6.2))
    rp = m["rp"]
    r_p, v_p = estado_periapsis(m["vinf_in_vec"], rp, GM_JUPITER, m["normal"])
    # base del plano de la hipérbola, para dibujarla sin proyección
    e1 = r_p / np.linalg.norm(r_p)
    e2 = v_p / np.linalg.norm(v_p)
    pts = np.array([propagar(r_p, v_p, t, GM_JUPITER)[0] for t in np.linspace(-6, 6, 801) * DIA])
    pts = np.c_[pts @ e1, pts @ e2] / R_JUPITER
    L = 50
    pts = pts[np.linalg.norm(pts, axis=1) < 0.9 * L]
    ax.plot(pts[:, 0], pts[:, 1], color=COL["sonda"], lw=2.2)
    ax.add_patch(Circle((0, 0), 1, color=COL["jupiter"]))
    ax.add_patch(Circle((0, 0), rp / R_JUPITER, fill=False, ls=":", ec=COL["gris"]))
    # flechas de sentido a mitad de cada rama (la sonda entra por abajo y sale por arriba)
    dist = np.linalg.norm(pts, axis=1)
    i_p = int(np.argmin(dist))
    for i, lab in [(int(np.argmin(np.abs(dist[:i_p] - 28))), "entrada"),
                   (i_p + int(np.argmin(np.abs(dist[i_p:] - 28))), "salida")]:
        ax.annotate("", xy=pts[i + 6], xytext=pts[i - 6],
                    arrowprops=dict(arrowstyle="-|>", color=COL["sonda"], lw=2, mutation_scale=22))
        ax.annotate(lab, pts[i], xytext=(12, -4), textcoords="offset points", fontsize=10.5)
    ax.text(0.97, 0.97, f"$e$ = {m['e']:.3f}\n$\\delta$ = {np.degrees(m['delta']):.1f}°\n"
            f"$r_p$ = {rp / R_JUPITER:.2f} $R_J$\n$v_\\infty$ = {m['vinf_in']:.2f} km/s",
            transform=ax.transAxes, fontsize=10, ha="right", va="top",
            bbox=dict(fc="white", ec=COL["rejilla"], boxstyle="round,pad=0.4"))
    ax.set_aspect("equal")
    ax.set_xlim(-L, L)
    ax.set_ylim(-L, L)
    ax.set_xlabel("$R_J$ (plano de la hipérbola)")
    ax.set_ylabel("$R_J$")
    ax.set_title("(a) Sistema de Júpiter: una hipérbola", loc="left", fontsize=12)

    vJ, vi, vo = m["vJ"], m["v_llega_J"], m["v_sale_J"]
    u1 = vJ / np.linalg.norm(vJ)
    u2 = np.cross(m["normal"], u1)
    u2 /= np.linalg.norm(u2)
    P = lambda v: np.array([v @ u1, v @ u2])  # noqa: E731
    vJ2, vi2, vo2 = P(vJ), P(vi), P(vo)

    def flecha(a, b, col, lw=2.2, ls="-"):
        bx.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="-|>", color=col, lw=lw, ls=ls,
                                                       shrinkA=0, shrinkB=0))
    th = np.linspace(0, 2 * np.pi, 300)
    bx.plot(vJ2[0] + m["vinf_in"] * np.cos(th), vJ2[1] + m["vinf_in"] * np.sin(th), ":", color=COL["gris"])
    flecha((0, 0), vJ2, COL["jupiter"], 2.8)
    flecha((0, 0), vi2, COL["gris"], 2, "--")
    flecha((0, 0), vo2, COL["sonda"], 2.6)
    flecha(vJ2, vi2, COL["gris"], 1.6)
    flecha(vJ2, vo2, COL["sonda"], 1.6)
    from matplotlib.lines import Line2D
    bx.legend(handles=[
        Line2D([], [], color=COL["jupiter"], lw=2.8, label=f"Júpiter: {np.linalg.norm(vJ):.2f} km/s"),
        Line2D([], [], color=COL["gris"], lw=2, ls="--", label=f"sonda antes: {np.linalg.norm(vi):.2f} km/s"),
        Line2D([], [], color=COL["sonda"], lw=2.6, label=f"sonda después: {np.linalg.norm(vo):.2f} km/s"),
        Line2D([], [], color=COL["gris"], lw=1.2, ls=":", label=f"$|\\vec v_\\infty|$ = {m['vinf_in']:.2f} km/s"),
    ], loc="upper left", fontsize=9.5, frameon=True, facecolor="white", framealpha=0.95, edgecolor="none",
        borderaxespad=0.2)
    for v2, txt, dx in [(vi2, r"$\vec v_{\infty}$ entrada", (10, 2)), (vo2, r"$\vec v_{\infty}$ salida", (0, 9))]:
        bx.annotate(txt, 0.5 * (vJ2 + v2), xytext=dx, textcoords="offset points", fontsize=9.5,
                    color=COL["gris"] if dx[0] else COL["sonda"], ha="left" if dx[0] else "center")
    bx.set_aspect("equal")
    allp = np.array([[0, 0], vJ2, vi2, vo2])
    bx.set_xlim(allp[:, 0].min() - 2, max(allp[:, 0].max(), vJ2[0] + m["vinf_in"]) + 3)
    bx.set_ylim(vJ2[1] - m["vinf_in"] - 2, vJ2[1] + m["vinf_in"] + 9)
    bx.set_xlabel("componente a lo largo de $\\vec v_J$ (km/s)")
    bx.set_ylabel("componente perpendicular (km/s)")
    bx.set_title(r"(b) Sistema del Sol: $\vec v = \vec v_J + \vec v_\infty$", loc="left", fontsize=12)
    titulo(fig, "El sobrevuelo en dos sistemas de referencia", x=0.01, ha="left", fontsize=15, fontweight="bold")
    nota(ax, 0.03, 0.03, f"Esfera de influencia: {R_SOI_JUPITER / R_JUPITER:.0f} $R_J$ (fuera del cuadro)",
         transform=ax.transAxes, va="bottom", fontsize=9, color=COL["gris"])
    fig.tight_layout()
    fig.savefig(FIG / "fig04_sobrevuelo.png")


if __name__ == "__main__":
    fig_ventanas()
    fig_pareto()
    m, d = cargar()
    fig_trayectoria(m, d)
    fig_sobrevuelo(m)
    print("Figuras 1–4 listas")

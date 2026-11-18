"""Datos de la misión óptima (Tierra -> Júpiter -> Saturno) para los videos de docs/animaciones/video.

Posiciones diarias de la Tierra, Júpiter y Saturno (DE421) y de la sonda sobre las dos cónicas
heliocéntricas del Modelo 1 (Lambert), más su rapidez. Salida: docs/animaciones/video/mision.json
"""
import json
from pathlib import Path

import numpy as np

from asistencia.constantes import AU, DIA, GM_SOL
from asistencia.efemerides import estado_planeta
from asistencia.kepler import propagar
from asistencia.modelo1 import evaluar_mision

SALIDA = Path(__file__).resolve().parent.parent / "docs" / "animaciones" / "video" / "mision.json"


def main():
    jd0, t1, t2 = 2443389.776902212, 586.0884799157138, 661.6068238714446
    m = evaluar_mision(jd0, t1, t2)
    n = int(round(t1 + t2)) + 1
    dias = np.linspace(0, t1 + t2, n)
    out = {"dias": [round(float(d), 2) for d in dias], "t_jupiter": t1, "t_saturno": t1 + t2,
           "tierra": [], "jupiter": [], "saturno": [], "sonda": [], "rapidez": []}
    for d in dias:
        jd = jd0 + d
        for k, nom in (("tierra", "Tierra"), ("jupiter", "Jupiter"), ("saturno", "Saturno")):
            r, _ = estado_planeta(nom, jd)
            out[k].append([round(float(r[0] / AU), 4), round(float(r[1] / AU), 4)])
        if d <= t1:
            r, v = propagar(m["rE"], m["v_salida"], d * DIA, GM_SOL)
        else:
            r, v = propagar(m["rJ"], m["v_sale_J"], (d - t1) * DIA, GM_SOL)
        out["sonda"].append([round(float(r[0] / AU), 4), round(float(r[1] / AU), 4)])
        out["rapidez"].append(round(float(np.linalg.norm(v)), 3))
    out["resumen"] = {"C3": m["C3"], "vinf": m["vinf_in"], "rp_RJ": m["rp_RJ"], "v_ent": float(np.linalg.norm(m["v_llega_J"])),
                      "v_sal": float(np.linalg.norm(m["v_sale_J"])), "vJ": float(np.linalg.norm(m["vJ"])),
                      "dE": float(0.5 * (np.dot(m["v_sale_J"], m["v_sale_J"]) - np.dot(m["v_llega_J"], m["v_llega_J"])))}
    SALIDA.parent.mkdir(parents=True, exist_ok=True)
    SALIDA.write_text(json.dumps(out, separators=(",", ":")))
    print("escrito", SALIDA, SALIDA.stat().st_size // 1024, "KB", out["resumen"])


if __name__ == "__main__":
    main()

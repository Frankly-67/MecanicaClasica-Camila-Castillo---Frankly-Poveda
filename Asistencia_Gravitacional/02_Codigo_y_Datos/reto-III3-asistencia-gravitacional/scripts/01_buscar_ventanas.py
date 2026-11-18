"""Búsqueda de ventanas Tierra -> Júpiter -> Saturno con el Modelo 1 (cónicas empalmadas).

Uso:  python scripts/01_buscar_ventanas.py [de421|circular]     (por omisión: de421)

Barre fecha de salida, tiempo de vuelo a Júpiter (TOF1) y a Saturno (TOF2). Una
combinación es factible si el sobrevuelo puede hacerse sin motor:
    | |v_inf_out| - |v_inf_in| | < TOL_DESAJUSTE    y    r_p >= RP_MIN_RJ
Las mejores se refinan con optimización continua (desajuste exactamente cero).
Las fechas reales de las Voyager NO entran en la búsqueda: se evalúan aparte para comparar.

Salidas: resultados/busqueda_<fuente>.npz, resultados/mision_optima_<fuente>.json
"""
import json
import sys
import time

import numpy as np
from estilo import RES

from asistencia.constantes import DIA, GM_JUPITER, GM_SOL, R_JUPITER
from asistencia.efemerides import estado_planeta, fecha_a_jd, jd_a_fecha
from asistencia.lambert import lambert
from asistencia.modelo1 import evaluar_mision, refinar

TOL_DESAJUSTE = 0.3     # km/s, solo para la rejilla
RP_MIN_RJ = 2.0         # hipótesis H7
LIMITES_C3 = [100.0, 105.0, 115.0]   # hipótesis H8: capacidad del lanzador, km^2/s^2

PASO = 4                # días; TOF1 usa múltiplos para que las llegadas caigan en la rejilla
JD0 = np.arange(fecha_a_jd("1976-01-01"), fecha_a_jd("1980-01-01"), PASO)
TOF1 = np.arange(400, 1101, 2 * PASO)
TOF2 = np.arange(400, 3001, 20)

REFERENCIAS = {   # fechas históricas (PDS Ring-Moon Systems Node), solo para comparar
    "Voyager 2": ("1977-08-20T14:29", "1979-07-09T22:29", "1981-08-26T03:24"),
    "Voyager 1": ("1977-09-05T12:56", "1979-03-05T12:05", "1980-11-12T23:46"),
}


def resumen(m):
    return {"salida": jd_a_fecha(m["jd_salida"]), "llegada_jupiter": jd_a_fecha(m["jd_jupiter"]),
            "llegada_saturno": jd_a_fecha(m["jd_saturno"]),
            "tof1_dias": round(float(m["tof1"]), 2), "tof2_dias": round(float(m["tof2"]), 2),
            "tiempo_total_anios": round(float((m["tof1"] + m["tof2"]) / 365.25), 3),
            "C3_km2_s2": round(float(m["C3"]), 2),
            "vinf_jupiter_entrada_km_s": round(float(m["vinf_in"]), 3),
            "vinf_jupiter_salida_km_s": round(float(m["vinf_out"]), 3),
            "desajuste_km_s": round(float(m["desajuste"]), 4),
            "delta_grados": round(float(np.degrees(m["delta"])), 2),
            "rp_radios_jupiter": round(float(m["rp_RJ"]), 2),
            "vinf_saturno_km_s": round(float(m["vinf_saturno"]), 3),
            # sin redondear, para que otros scripts reconstruyan exactamente la misma misión
            "exacto": {"jd_salida_tdb": float(m["jd_salida"]), "tof1_dias": float(m["tof1"]),
                       "tof2_dias": float(m["tof2"])}}


def main(fuente):
    t_ini = time.time()
    # --- Tramo 1: Tierra -> Júpiter
    C3 = np.full((len(JD0), len(TOF1)), np.nan)
    VIN = np.full((len(JD0), len(TOF1), 3), np.nan)
    for i, t0 in enumerate(JD0):
        rE, vE = estado_planeta("Tierra", t0, fuente)
        for j, T1 in enumerate(TOF1):
            rJ, vJ = estado_planeta("Jupiter", t0 + T1, fuente)
            try:
                v1, v2 = lambert(rE, rJ, T1 * DIA, GM_SOL)
            except ValueError:
                continue
            C3[i, j] = np.dot(v1 - vE, v1 - vE)
            VIN[i, j] = v2 - vJ
    print(f"[{fuente}] tramo 1 listo ({time.time() - t_ini:.0f} s)")

    # --- Tramo 2: Júpiter -> Saturno, por fecha de llegada a Júpiter
    jdJ = np.arange(JD0[0] + TOF1[0], JD0[-1] + TOF1[-1] + 1, PASO)
    VOUT = np.full((len(jdJ), len(TOF2), 3), np.nan)
    for k, tJ in enumerate(jdJ):
        rJ, vJ = estado_planeta("Jupiter", tJ, fuente)
        for m, T2 in enumerate(TOF2):
            rS, _ = estado_planeta("Saturno", tJ + T2, fuente)
            try:
                v3, _ = lambert(rJ, rS, T2 * DIA, GM_SOL)
            except ValueError:
                continue
            VOUT[k, m] = v3 - vJ
    print(f"[{fuente}] tramo 2 listo ({time.time() - t_ini:.0f} s)")

    # --- Combinar: qué sobrevuelo haría falta en cada caso
    kJ = np.rint((JD0[:, None] + TOF1[None, :] - jdJ[0]) / PASO).astype(int)
    vin, vout = VIN[:, :, None, :], VOUT[kJ]
    nin, nout = np.linalg.norm(vin, axis=-1), np.linalg.norm(vout, axis=-1)
    delta = np.arccos(np.clip(np.sum(vin * vout, axis=-1) / (nin * nout), -1, 1))
    with np.errstate(divide="ignore", invalid="ignore"):
        rp = GM_JUPITER * (1 / np.sin(delta / 2) - 1) / (0.5 * (nin + nout)) ** 2 / R_JUPITER
    fact = (np.abs(nout - nin) < TOL_DESAJUSTE) & (rp >= RP_MIN_RJ)
    T_total = (TOF1[None, :, None] + TOF2[None, None, :]) / 365.25
    C3b = np.broadcast_to(C3[:, :, None], fact.shape)
    Tf = np.where(fact, np.broadcast_to(T_total, fact.shape), np.inf)
    print(f"[{fuente}] combinaciones: {fact.size:,}; factibles: {int(fact.sum()):,}")

    T_min = {}
    for L in LIMITES_C3:
        Tl = np.where(C3b <= L, Tf, np.inf).min(axis=2)
        T_min[L] = np.where(np.isfinite(Tl), Tl, np.nan)
    limites = np.arange(76, 131, 2.0)
    frente = np.array([float(np.min(np.where(C3b <= L, Tf, np.inf))) for L in limites])
    np.savez_compressed(RES / f"busqueda_{fuente}.npz", jd0=JD0, tof1=TOF1, tof2=TOF2, C3=C3,
                        **{f"T_min_C3_{int(L)}": T_min[L] for L in LIMITES_C3},
                        limites_C3=limites, frente_T=frente)

    # --- Refinar desde el mejor punto de la rejilla para cada límite de C3
    salida = {"fuente_efemerides": fuente,
              "criterios": {"tol_desajuste_rejilla_km_s": TOL_DESAJUSTE, "rp_min_RJ": RP_MIN_RJ},
              "minimo_tiempo_con_limite_C3": {}}
    for L in LIMITES_C3:
        Tl = np.where(C3b <= L, Tf, np.inf)
        i, j, k = np.unravel_index(np.argmin(Tl), Tl.shape)
        m, res = refinar(JD0[i], TOF1[j], TOF2[k], RP_MIN_RJ, c3_max=L, fuente=fuente)
        ok = (res.success and abs(m["desajuste"]) < 1e-6 and m["C3"] <= L + 1e-3
              and m["rp_RJ"] >= RP_MIN_RJ - 1e-3)
        salida["minimo_tiempo_con_limite_C3"][f"C3<={L:g}"] = resumen(m) if ok else None
        print(f"[{fuente}] C3 <= {L:g}: {'ok' if ok else 'no convergió'}")
    salida["fechas_historicas_evaluadas_con_modelo1"] = {
        n: resumen(evaluar_mision(fecha_a_jd(a), fecha_a_jd(b) - fecha_a_jd(a),
                                  fecha_a_jd(c) - fecha_a_jd(b), fuente))
        for n, (a, b, c) in REFERENCIAS.items()}
    (RES / f"mision_optima_{fuente}.json").write_text(json.dumps(salida, indent=2, ensure_ascii=False))
    print(json.dumps(salida, indent=2, ensure_ascii=False))
    print(f"[{fuente}] listo en {time.time() - t_ini:.0f} s")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "de421")

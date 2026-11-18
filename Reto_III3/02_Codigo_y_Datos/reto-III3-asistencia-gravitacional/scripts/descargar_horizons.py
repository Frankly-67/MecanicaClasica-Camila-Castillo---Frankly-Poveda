"""(Opcional) Descarga vectores de JPL Horizons para comparar con DE421 y con el kernel SPICE.

Correr en un computador con internet:
    pip install astroquery
    python scripts/descargar_horizons.py
Guarda CSV en datos/horizons/: posición (km) y velocidad (km/s) heliocéntricas, eclíptica J2000.
IDs: Júpiter = 599, Saturno = 699, Tierra = 399, Voyager 2 = -32, Voyager 1 = -31.
"""
from pathlib import Path

from astroquery.jplhorizons import Horizons

DATOS = Path(__file__).resolve().parents[1] / "datos" / "horizons"
DATOS.mkdir(parents=True, exist_ok=True)
AU_KM, DIA_S = 1.495978707e8, 86400.0
CUERPOS = {"jupiter": "599", "saturno": "699", "tierra": "399", "voyager2": "-32", "voyager1": "-31"}
EPOCAS = {"crucero": {"start": "1977-08-25", "stop": "1981-09-30", "step": "1d"},
          "encuentro": {"start": "1979-06-05", "stop": "1979-07-20", "step": "1h"}}

if __name__ == "__main__":
    for nombre, idh in CUERPOS.items():
        for tramo, ep in EPOCAS.items():
            if tramo == "encuentro" and nombre not in ("jupiter", "voyager2"):
                continue
            tab = Horizons(id=idh, location="@sun", epochs=ep, id_type=None).vectors(refplane="ecliptic")
            tab = tab[["datetime_jd", "datetime_str", "x", "y", "z", "vx", "vy", "vz"]]
            for c in "xyz":
                tab[c] = tab[c] * AU_KM
                tab["v" + c] = tab["v" + c] * AU_KM / DIA_S
            tab.write(DATOS / f"{nombre}_{tramo}.csv", format="csv", overwrite=True)
            print(nombre, tramo, len(tab))

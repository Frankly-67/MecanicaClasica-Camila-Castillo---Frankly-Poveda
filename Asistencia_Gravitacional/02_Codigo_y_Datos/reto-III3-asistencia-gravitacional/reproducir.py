"""Regenera TODOS los resultados y figuras desde cero.

    python reproducir.py            # figuras para pantalla (figuras/)
    python reproducir.py --reporte  # además, versión para el informe LaTeX (figuras/reporte/)

No hay números aleatorios en el proyecto, así que no se necesita semilla: el resultado
es determinista. Tolerancias declaradas: Lambert (Brent) xtol = rtol = 1e-12;
propagación kepleriana |Δχ| < 1e-10; optimización SLSQP ftol = 1e-10;
integrador de control DOP853 rtol = 1e-12 a 1e-13;
REBOUND/IAS15 epsilon = 1e-9 (barrido) y 1e-10 (Voyager 2);
Modelo 2 ampliado (satélites y J2): DOP853 rtol = 1e-12, atol = 1e-6 km.
"""
import os
import subprocess
import sys
import time
from pathlib import Path

RAIZ = Path(__file__).resolve().parent
PASOS = [
    [sys.executable, "-m", "pytest", "-q"],
    [sys.executable, "scripts/01_buscar_ventanas.py", "circular"],
    [sys.executable, "scripts/01_buscar_ventanas.py", "de421"],
    [sys.executable, "scripts/02_figuras_modelo1.py"],
    [sys.executable, "scripts/03_voyager2.py"],
    [sys.executable, "scripts/05_modelo2.py"],
    [sys.executable, "scripts/06_voyager2_satelites.py"],
    [sys.executable, "scripts/07_exportar_animaciones.py"],
    [sys.executable, "scripts/04_verificar_lagrangiano.py"],   # ≈ 1 min (sympy)
]


def correr(cmd, env=None):
    t = time.time()
    print(">>", " ".join(cmd[1:]), flush=True)
    r = subprocess.run(cmd, cwd=RAIZ, env=env, stdout=subprocess.DEVNULL if "pytest" not in cmd[2:3] else None)
    if r.returncode:
        sys.exit(f"Falló: {' '.join(cmd)}")
    print(f"   ok ({time.time() - t:.0f} s)")


if __name__ == "__main__":
    for cmd in PASOS:
        correr(cmd)
    if "--reporte" in sys.argv:
        env = {**os.environ, "REPORTE": "1"}
        for cmd in PASOS[3:7]:
            correr(cmd, env)
    print("Todo regenerado: resultados/ y figuras/")

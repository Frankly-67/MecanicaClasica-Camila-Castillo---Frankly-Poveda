# Bitácora de interacción con IA (metodología 4D)

Formato tomado de la presentación del curso "Las cuatro D de la colaboración humano–IA"
(L. A. Núñez, 2026, basada en Dakan y Feller 2025):

- **Delegación:** qué se conservó para el humano y qué se delegó. Regla: conservar lo que es
  objetivo de aprendizaje y no puede verificarse; delegar lo verificable.
- **Descripción:** el prompt, con sistema de referencia, supuestos y forma de la respuesta.
- **Discernimiento:** cada respuesta es una hipótesis. Mínimo una verificación dimensional y un
  caso límite (capas: superficial, estructural, profunda).
- **Diligencia:** verificar antes de confiar, declarar, registrar (prompt, respuesta,
  aceptado/rechazado/corregido, cómo se verificó) y re-derivar lo corregido.

**Declaración de uso.** Asistente: Claude (Anthropic). La IA escribió la mayor parte del
código, las figuras y borradores de texto (marco teórico, informes de avance). Los estudiantes
(Frankly Giovanni Poveda Pinzón y Camila Castillo López) eligieron el reto, decidieron el enfoque y revisan cada entrada. Esta versión de la bitácora se
reconstruyó el 7 nov 2026 a partir de la conversación; las entradas 1–7 reemplazan la
versión anterior. La sesión 5 (10 nov) es una auditoría completa pedida por los estudiantes y
hecha con otra versión del asistente: encontró ocho errores en el trabajo anterior del
propio asistente, cada uno con evidencia que no depende de su código (valores publicados,
un balance de energía, una fuente histórica).

> Pendiente para los estudiantes: revisar cada entrada, completar la columna "Revisado por"
> y pegar los prompts textuales que falten.

---

## Sesión 1 (25 sep 2026): plan, visualización y Modelo 1

**Delegación.** Se delegó: lectura del enunciado, plan de trabajo, simulación interactiva
para entender los dos modelos, código de Kepler/Lambert/sobrevuelo y sus pruebas. Se
conservó: elección del reto y del enfoque, interpretación.

| # | Prompt (resumen) | Respuesta / afirmación de la IA | Estado | Verificación | Revisado por |
|---|---|---|---|---|---|
| 1 | "Hazme un plan para el Reto III.3" | "En el sistema de Júpiter \|v∞ entrada\| = \|v∞ salida\|" | Aceptado | Integración numérica independiente (DOP853) de la hipérbola desde el periapsis: igualdad a 10⁻⁶ km/s (`test_hiperbola_integrada_confirma_giro_y_simetria`) | |
| 2 | idem | "sin(δ/2) = 1/e, e = 1 + r_p v∞²/μ" | Aceptado | Dimensional (`test_analisis_dimensional_excentricidad`); límites μ→0 y r_p→∞ (δ→0); ángulo integrado vs fórmula < 0.5°. Re-derivación en el marco teórico | |
| 3 | idem | "La esfera de influencia de Júpiter mide unos 48 millones de km" | Aceptado | Cálculo propio a(m/M)^(2/5) = 4.82×10⁷ km = 674 R_J | |
| 4 | "Comencemos con el trabajo" (código de Lambert) | Primera versión de `lambert.py` | **Corregido** | `test_lambert_llega_al_punto_pedido` falló para transferencias de 230°: la búsqueda del intervalo suponía F(z) creciente desde z negativo. Se cambió por una malla con cambio de signo + Brent. Evidencia: el caso de 230° no llegaba al punto pedido | |
| 5 | idem | Prueba "r_p = 10⁴ R_J no desvía (δ < 0.5°)" | **Rechazado** (la prueba, no el código) | Cálculo a mano: e = 204, δ = 0.56°. La tolerancia estaba mal elegida. Se reemplazó por prueba de monotonía y límite | |
| 6 | idem | "Las fechas óptimas del modelo coinciden con las de Voyager 2 con C3 ≤ 85" | **Corregido en la sesión 3** | Con órbitas circulares era cierto (24 ago 1977 / 8 jul 1979). Con efemérides reales, en las fechas de Voyager 2 se necesita C3 = 103 km²/s² y con C3 ≤ 85 no hay solución. La conclusión dependía de la hipótesis H1 | |

**Diligencia.** 17 pruebas al cierre de la sesión; los números del informe salen de
`resultados/`, no del chat.

## Sesión 2 (2 oct 2026): identificación del problema y marco teórico

**Delegación.** Se delegó la redacción del borrador del informe de avance y del marco
teórico; se conservó la decisión de qué entregar y la revisión.

| # | Prompt (resumen) | Respuesta / afirmación de la IA | Estado | Verificación | Revisado por |
|---|---|---|---|---|---|
| 7 | "Avance de identificación del problema en LaTeX" | Distancias de máximo acercamiento: Voyager 2 a 10.1 R_J, Voyager 1 a 4.9 R_J | Aceptado | PDS Ring-Moon Systems Node: 721 670 km y 348 890 km. En la sesión 3 se confirmó Voyager 2 con el kernel SPICE: 721 557 km | |
| 8 | idem | "En el borde de la esfera de influencia el Sol atrae a la sonda unas 4 veces más que Júpiter" | Aceptado | Cálculo propio: (GM☉/a_J²)/(GM_J/r_SOI²) = 4.0. Orden de magnitud y dimensiones ✓ | |
| 9 | "Marco teórico" | Δε☉ = v_J · (v∞,sal − v∞,ent) | Aceptado, falta re-derivar | Signo: si v∞ gira hacia v_J, Δε > 0 ✓. Límite δ→0: Δε→0 ✓. **Pendiente que los estudiantes lo deriven a mano** | |
| 10 | idem | Primer borrador: "la razón de masas es del orden de 10⁻²⁷" | **Corregido** | Orden de magnitud: m_sonda/m_J ≈ 800 kg / 1.9×10²⁷ kg ≈ 4×10⁻²⁵ | |

## Sesión 3 (7 nov 2026): reconstrucción con datos reales

**Delegación.** Se delegó buscar fuentes de datos y reescribir el código; se conservó la
decisión de usar datos reales y la interpretación de los resultados.

**Descripción.** "Reconstruye el repositorio; también toca usar datos reales", con el
enunciado completo del reto y las pautas de evaluación adjuntos.

| # | Prompt (resumen) | Respuesta / afirmación de la IA | Estado | Verificación | Revisado por |
|---|---|---|---|---|---|
| 11 | "Usar datos reales" | Proponer DE421 (PyPI) y el kernel SPICE `vgr2_jup230.bsp` porque JPL Horizons no era accesible | Aceptado | Kernel: máximo acercamiento 22:25 UTC a 721 557 km vs publicado 22:29 UTC y 721 670 km (`test_voyager2_maximo_acercamiento_coincide_con_lo_publicado`). **La diferencia de 4 min era un error numérico propio: ver la entrada 25.** DE421: Tierra a 0.9833 UA y 100.4° en J2000; distancias e inclinaciones de J y S dentro de lo publicado (`test_datos.py`) | |
| 12 | idem | "jplephem puede leer el kernel de Voyager" | **Rechazado** | Al ejecutarlo: el segmento es SPK tipo 1, que jplephem no implementa. Se usó CSPICE (spiceypy) | |
| 13 | idem | Optimización con DE421 | **Corregido** | Falló: el optimizador visitaba puntos donde Lambert no tiene solución. Se añadió un valor de penalización. Comprobación: los óptimos cumplen desajuste = 0, C3 ≤ límite, r_p ≥ 2 R_J | |
| 14 | idem | "Con efemérides reales el modelo predice mejor el sobrevuelo real" | Aceptado | Datos no usados en el ajuste: r_p predicho 9.72 R_J vs real 10.09 R_J (Voyager 2) y 4.74 vs 4.88 R_J (Voyager 1); con órbitas circulares 12.6 y 5.75 R_J | |
| 15 | idem | "El parámetro de Tisserand se conserva en el encuentro" | **Puesto en duda por los datos** | Con Voyager 2 real: T_J = 2.549 antes y 2.481 después (−2.7 %). Posibles causas: Júpiter no está en órbita circular (e = 0.048) y los extremos del registro están dentro de la esfera de influencia. Se verificará con el Modelo 2. **Verificado en la sesión 5: la causa es la excentricidad de Júpiter, no la cercanía (entrada 28)** | |

| 16 | "Haz tú los cálculos de la formulación lagrangiana y pásamelos en PDF" | Lagrangiano en el sistema rotante, ecuaciones de Euler–Lagrange, h = E − n ℓ_z (Jacobi), hamiltoniano, límite de Tisserand | Aceptado tras verificar. **Nota de delegación:** según la regla del curso, esta derivación es objetivo de aprendizaje; se delegó por tiempo y los estudiantes deben re-derivarla a mano antes de la defensa | (a) sympy, independiente del texto: E–L, dh/dt = 0, h = E − nℓ_z y límite de Tisserand, los cinco verdaderos (`04_verificar_lagrangiano.py`). (b) Numérica: h constante a 3.7×10⁻¹² mientras E pasa de −134 a +22 km²/s²; ΔE = nΔℓ_z = 155.75. (c) Dimensional y límites M₂→0, n→0, T_J = 3. (d) Datos reales de Voyager 2: ΔE = 172.42 vs nΔℓ_z = 173.38 (0.6 %). **La comprobación (d) estaba mal hecha: ver la entrada 23** | |

## Sesión 4 (7 nov 2026): Modelo 2

**Delegación.** Se delegó el código del Modelo 2 y del barrido; se conservó la elección de
la medida del error y la interpretación.

| # | Prompt (resumen) | Respuesta / afirmación de la IA | Estado | Verificación | Revisado por |
|---|---|---|---|---|---|
| 17 | "Siguiente" (Modelo 2) | `sim.ri_ias15.epsilon` para fijar la tolerancia de REBOUND | **Rechazado** | Error al ejecutar: en REBOUND 5 el atributo es `sim.integrator.epsilon`. El asistente usaba una API antigua | |
| 18 | idem | Propagación kepleriana con la estimación inicial de Curtis para hipérbolas | **Corregido** | El barrido hizo fallar a Newton para tiempos largos en hipérbolas. Se usó la estimación de Vallado; las 25 pruebas previas siguen pasando | |
| 19 | idem | Medir la violación de Jacobi con el Sol en el origen | **Corregido** | Con esa fórmula el propio Modelo 2 "violaba" Jacobi en 9×10⁻⁴, lo que es imposible (es su invariante exacto). Con la fórmula baricéntrica exacta: 2×10⁻¹⁴ | |
| 20 | idem | Función `conicas_empalmadas` | **Corregido** | La prueba de límite (cerca del periapsis los modelos deben coincidir) falló por 2.9×10⁶ km: la función no contemplaba tiempos menores que el empalme | |
| 21 | idem | "El error de las cónicas crece al bajar v∞" (predicción 3, escrita antes de simular) | Aceptado | Barrido: 17.5 % a 2.5 km/s, 0.13 % a 15 km/s (R_e = r_SOI). Comparación entre dos métodos independientes | |
| 22 | idem | "El salto de energía de Voyager 2 en el periapsis se debe a los satélites galileanos" | **Hipótesis, no verificada** | El Modelo 2 sin satélites no reproduce el salto, y el error que acumula (6.8×10⁴ km en 11 días ≈ 0.07 km/s) es del mismo orden que el salto (0.06 km/s). Falta una efeméride de los satélites para confirmarlo. **Resuelto en la sesión 5: los satélites explican la mitad; la otra mitad es una maniobra (entrada 24)** | |

## Sesión 5 (10 nov 2026): auditoría completa

**Delegación.** Se delegó una revisión de todo el proyecto (código, física, cifras de cada
documento, figuras) con la instrucción de no dar por bueno nada de las sesiones anteriores.
Se conservó la decisión de qué correcciones entran en la entrega y la re-derivación a mano.

**Descripción.** "Vuelve a analizar todo, y hazlo de nuevo".

**Discernimiento.** Regla usada: una cifra solo se acepta si coincide con algo que no sale del
mismo cálculo. Las que coincidían "demasiado bien" o "casi" con lo publicado se trataron como
sospechosas (entradas 23 y 25).

| # | Afirmación revisada (de sesiones anteriores) | Estado | Cómo se detectó y qué quedó | Revisado por |
|---|---|---|---|---|
| 23 | "Con Voyager 2, ΔE = 172.4 y nΔℓ_z = 173.4: Jacobi se cumple al 0.6 %" (entrada 16d) | **Rechazado** | Al rehacer el balance con la definición de la formulación: la E usada no incluía el potencial de Júpiter (cambia 8.1 km²/s² entre los extremos del registro) y se usó el n medio, cuando Júpiter, cerca del afelio (5.33 UA), giraba a 0.952 n. Los dos errores, de 5 % cada uno, se compensaban. Correcto: ΔE = 164.20; n medio × Δℓ_z = 173.21 (5.5 % de diferencia); ∫ω dℓ_z con la velocidad angular instantánea = 164.83 (0.4 %). La relación exacta solo vale para órbita circular; se añadió a la formulación la versión con ω(t) (sección 5.6) | |
| 24 | "El salto de energía de Voyager 2 se debe a los satélites galileanos" (entrada 22) | **Corregido: cierto solo a medias** | Se añadieron los satélites y J2 al Modelo 2 (`jupiter_real.py`), validados contra las distancias publicadas de los sobrevuelos (Ganímedes: 62 078 km frente a 62 130). Explican −0.029 km/s de la caída de v∞. La aceleración residual (real menos gravitacional) mostró el resto: 2.83 mm/s² constantes durante 76.6 min, desde 130 min después del máximo acercamiento (Δv = 13.0 m/s, −0.033 km/s en v∞). Fuente independiente: NASA SP-439 (Morrison y Samz 1980) describe un encendido de 76 minutos unas dos horas después del máximo acercamiento. Con satélites, J2 y maniobra: v∞ final 7.5780 frente a 7.5777 real; el tramo de 11 días posterior a la maniobra se reproduce con 6 km de error | |
| 25 | "El kernel da el máximo acercamiento a las 22:25 UTC a 721 557 km; lo publicado es 22:29" (entrada 11) | **Corregido** | La diferencia de 4 min con lo publicado se había aceptado como tolerable. Causa: el refinado del mínimo se hacía en días julianos, y la tolerancia del método acotado de scipy es relativa al valor de la variable (1.5×10⁻⁸·2.4×10⁶ días = 50 min). Rehecho en segundos: 22:28:52 UTC y 721 552 km respecto al baricentro joviano; 22:29:01 UTC y 721 673 km respecto al centro del planeta. Publicado: 22:29 UTC, 721 670 km | |
| 26 | "TDB = UTC + 49.18 s (1979)" (H11) | **Corregido** | En 1979 TAI − UTC = 18 s, es decir 50.184 s; 49.184 s es el valor de 1978. Se reemplazó por la tabla de segundos intercalares del IERS. Prueba: la época J2000 (11:58:55.816 UTC del 1 ene 2000) da JD 2451545.0 exacto. Efecto en los resultados: ninguno a la precisión citada | |
| 27 | "Saturno detrás de Júpiter (−60°), como en 1979" (prueba de la hipótesis H12) | **Rechazado** | DE421: en julio de 1979 Saturno iba 27.8° por delante de Júpiter. Con la posición correcta el efecto de Saturno sobre la ganancia de energía es 1–4×10⁻⁴, no 3–6×10⁻⁵. Sigue siendo despreciable | |
| 28 | "Tisserand cambia 2.7 % porque la sonda no está lejos de Júpiter" y "la energía pasa de −121 a +51 km²/s²" (entrada 15, informe) | **Corregido** | Se prolongó la trayectoria real 200 días a cada lado con el Modelo 2 ampliado, hasta más de 3 r_SOI. Lejos de Júpiter la energía heliocéntrica pasa de −121.8 a +26.3 km²/s² (los +51 se medían dentro de la esfera de influencia, donde falta el potencial de Júpiter) y Tisserand pasa de 2.576 a 2.670 (3.6 %): no se conserva ni siquiera lejos. La causa es la excentricidad de la órbita de Júpiter, la misma de la entrada 23. La conclusión de que la sonda queda en escape se mantiene (e = 1.30) | |
| 29 | Solucionador de Lambert (versión corregida en la entrada 4) | **Corregido** | Prueba nueva con transferencias muy rápidas (90° en 60 días): respondía "no hay solución" cuando sí la hay, porque exigía que los dos puntos vecinos de la malla fueran válidos. Fuera del intervalo usado en la búsqueda (≥ 400 días): los mapas de ventanas no cambian | |
| 30 | Figura del sobrevuelo óptimo | **Corregido** | El recuadro de la figura decía r_p = 5.9 R_J y δ = 98.0°, y el texto 5.96 y 97.9°. La figura reconstruía la misión con la fecha redondeada al día (C3 = 100.02, desajuste −0.004 km/s). Ahora se guardan y se usan los valores sin redondear | |
| 31 | "El Modelo 2 pasa a 1 300 km de la trayectoria real en el periapsis" (entrada 22) | Aceptado | Segundo método: DOP853 con el Sol, Júpiter y Saturno prescritos por DE421, en vez de integrarlos con REBOUND. Da 1 305 km en el periapsis y 67 582 km al final, frente a 1 310 y 67 600 km | |
| 32 | Reproducibilidad | Aceptado | Clon limpio del repositorio y `python reproducir.py --reporte`: los JSON y NPZ regenerados son idénticos byte a byte a los entregados | |

**Diligencia.** 40 pruebas (10 nuevas: escala de tiempo, sobrevuelos de los cuatro satélites,
centro de Júpiter, maniobra, tramo balístico, Lambert rápido). Todas las cifras de los
documentos se volvieron a leer de `resultados/*.json`.

## Errores detectados (resumen para el informe)

1. Lambert fallaba para transferencias > 180° (entrada 4). Evidencia: prueba de integración.
2. Una prueba escrita por la IA tenía una tolerancia imposible (entrada 5). Evidencia: cálculo a mano.
3. La conclusión "C3 ≤ 85 reproduce a Voyager 2" dependía de órbitas circulares (entrada 6). Evidencia: efemérides reales.
4. Orden de magnitud de la razón de masas (entrada 10).
5. La herramienta propuesta no leía el formato de los datos (entrada 12).
6. API desactualizada de REBOUND (entrada 17).
7. Medida aproximada que hacía parecer que el Modelo 2 violaba su propio invariante (entrada 19). Evidencia: un invariante exacto no puede variar en 10⁻⁴.
8. Función de cónicas incorrecta dentro de la esfera (entrada 20). Evidencia: prueba de límite.
9. La comprobación de Jacobi con datos reales omitía el potencial de Júpiter y usaba el n medio; la coincidencia del 0.6 % era fortuita (entrada 23). Evidencia: balance de energía con la definición de la formulación.
10. El salto de energía de Voyager 2 se atribuyó solo a los satélites; la mitad es una maniobra de 76 minutos (entrada 24). Evidencia: aceleración residual y NASA SP-439.
11. El mínimo de distancia se refinaba en días julianos y quedaba con minutos de error (entrada 25). Evidencia: hora y distancia publicadas.
12. Segundos intercalares de 1979 (entrada 26). Evidencia: tabla del IERS, época J2000.
13. Posición de Saturno en 1979 con el signo cambiado (entrada 27). Evidencia: DE421.
14. Explicación equivocada de por qué no se conserva Tisserand, y energía "de escape" medida dentro de la esfera de influencia (entrada 28). Evidencia: prolongación de la trayectoria.
15. Lambert perdía la solución de transferencias muy rápidas (entrada 29). Evidencia: prueba de integración.
16. Figura hecha con valores redondeados que no coincidía con el texto (entrada 30).

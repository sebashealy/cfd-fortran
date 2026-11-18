# cfd-fortran

Programas en Fortran de dinámica de fluidos y transferencia de calor computacional (desarrollados en el curso de **Análisis computacional en ingeniería**). Resuelven
ecuaciones de transporte (conducción, advección-difusión, Navier-Stokes y
convección natural) con el **método de volumen finito** sobre mallas
estructuradas, usando Jacobi y Gauss-TDMA línea por línea como solucionadores
lineales y SIMPLEC / proyección para el acoplamiento presión-velocidad.

## Programas

| # | Archivo | Qué resuelve |
|---|---------|--------------|
| 1 | `p01_calor1d.f90` | Conducción de calor 1D con generación de calor y temperaturas fijas en los extremos (Jacobi). Compara contra la solución analítica. |
| 2a | `p02a_adv_dif_central.f90` | Advección-difusión 1D, Pe = 100, esquema de diferencias centrales (Jacobi). |
| 2b | `p02b_adv_dif_upwind.f90` | Igual que 2a con esquema **upwind**, doble precisión y paro por tolerancia. Compara contra la solución analítica. |
| 3 | `p03_poisson2d.f90` | Ecuación de Poisson 2D en [-1, 1]² con Gauss-TDMA 2D. Escribe `phi`, `phia` y el error. |
| 4 | `p04_advdif2d.f90` | Advección-difusión 2D estacionaria con campo de velocidades prescrito; condiciones de Neumann en este y oeste. |
| 5 | `p05_advdif2d_temporal.f90` | Versión transitoria de 4 con upwind; genera `anim.gnp` para animar en gnuplot. |
| 6 | `p06_lid_driven_simplec.f90` | Cavidad con tapa deslizante (Re = 1000) con SIMPLEC en malla escalonada. |
| 7 | `p07_flujo_canal_simplec.f90` | Flujo en canal de 10 × 1 con entrada uniforme y salida tipo Neumann (SIMPLEC). |
| 8 | `p08_canal_obstaculo.f90` | Flujo en canal con un obstáculo sólido (frontera inmersa), Re = 100, upwind. |
| 9 | `p09_dipolo_magnetico.f90` | Flujo forzado por la fuerza de Lorentz de un imán permanente, con seguimiento lagrangiano de una interfaz de partículas. |
| 10 | `p10_proyeccion.f90` | Cavidad con tapa deslizante con el método de proyección (paso fraccionado). |
| 11 | `p11_conveccion_natural.f90` | Convección natural en cavidad cuadrada (Ra = 8·10³, Pr = 0.7): energía + Navier-Stokes con SIMPLEC. |
| 12 | `p12_difusion_2d_cilindrica.f90` | Difusión 2D en coordenadas cilíndricas (r, θ) en un sector de θ ∈ [π/4, 3π/4]. |

## Compilar

Requiere `gfortran` (probado con GNU Fortran 13).

```bash
make              # compila todo en bin/
make bin/p03_poisson2d
make clean
```

## Ejecutar

Cada programa escribe sus resultados (`.txt`, `anim.gnp`) en el directorio
actual, así que conviene correrlos cada uno en su propia carpeta:

```bash
mkdir -p salidas/p03 && cd salidas/p03
../../bin/p03_poisson2d
```

Los resultados se pueden graficar con gnuplot, por ejemplo
`splot 'phi0.txt' w pm3d` o `load 'anim.gnp'` para las animaciones.

Los programas 5 a 11 avanzan decenas de miles de pasos de tiempo, así que
pueden tardar bastante en terminar.

# Dinámica del Electrón en un Campo de Anillo Cargado

Este proyecto simula el movimiento de un electrón en el campo eléctrico creado por un anillo de cargas positivas, utilizando integración de Runge-Kutta de 4to orden (RK4). La simulación incluye visualización 3D de la trayectoria y análisis de energía.

## Análisis y Visualización

Los resultados de la simulación pueden explorarse en el siguiente notebook:

[![Ver Notebook de Análisis](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.org/github/isaultirado77/particle_on_charged_field/blob/main/simulation_analysis.ipynb)
```

## Estructura del Proyecto
- `electron_on_charged_field.py`: Script principal de la simulación (implementación de RK4)
- `plotting.py`: Visualización de trayectorias y dinámica de energía
- `animator.py`: Animación 3D del movimiento del electrón
- `simulation_analysis.ipynb`: Análisis teórico y discusión de resultados

Directorios de salida:
- `data/`: Archivos de datos de la simulación
- `plots/`: Gráficos y figuras generadas
- `animations/`: Animaciones de la trayectoria 3D

## Ejecución de la Simulación
Configura los parámetros desde la línea de comandos:

```bash
python electron_on_charged_field.py \
    --n_charges 50 \
    --ring_radius 1.0 \
    --charge_ring 50.0e-9 \
    --electron_x -1.0 \
    --electron_vy 6.0e4 \
    --dt 1e-9 \
    --t_max 7e-6 \
    --filename "mi_simulacion"
```

Parámetros clave:
- `n_charges`: Número de cargas puntuales en el anillo (se recomiendan ≥20)
- `ring_radius`: Radio del anillo en metros
- `charge_ring`: Carga total del anillo en Coulombs
- `electron_x/vy/vz`: Componentes de la posición/velocidad iniciales
- `dt`: Paso de tiempo (más pequeño para mayor precisión)
- `t_max`: Tiempo total de simulación

## Requisitos
- Python 3.7+
- Paquetes necesarios:
```bash
pip install numpy matplotlib scipy
```
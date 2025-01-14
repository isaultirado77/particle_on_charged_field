# Simulación de un Electrón en un Campo Eléctrico

Este proyecto simula el movimiento de un electrón bajo la influencia de un campo eléctrico generado por un anillo cargado. La simulación utiliza métodos numéricos para calcular la dinámica del sistema y generar resultados como las trayectorias y las energías del electrón.

## Descripción

Se considera un anillo cargado uniformemente con una carga total \( Q = 50 \times 10^{-9} \, \text{C} \). Un electrón es liberado en una posición inicial definida, y su movimiento es determinado únicamente por la interacción eléctrica con el anillo. La simulación permite calcular:
- La trayectoria del electrón.
- La energía cinética, potencial y total del sistema en función del tiempo.

### Características
- Representación modular del sistema (anillo y electrón).
- Uso del método de Euler para la integración numérica.
- Registro de datos de energía para análisis detallado.

## Requisitos

- Python 3.7 o superior.
- Librerías necesarias:
  - `numpy`

Instala las dependencias con:
```bash
pip install numpy
```

## Ejecución

1. Clona este repositorio:
```bash
git clone <https://github.com/isaultirado77/particle_on_charged_field>
```
2. Navega al directorio del proyecto:
```bash
cd electron_field_simulation
```
3. Ejecuta la simulación:
```bash
python electron_field_simulation.py
```

## Ejemplo de Configuración Inicial
- Posición inicial del electrón: \( (0.15, 0, 0) \, \text{m} \).
- Velocidad inicial del electrón: \( (0, 0, 0) \, \text{m/s} \).
- Número de partículas en el anillo: 100.
- Radio del anillo: \( 0.1 \, \text{m} \).
- Tiempo total de simulación: \( 1 \times 10^{-5} \, \text{s} \).
- Paso de tiempo: \( 1 \times 10^{-17} \, \text{s} \).

## Resultados

La simulación genera:
- **Trayectoria del electrón**: Información sobre su movimiento en el espacio.
- **Energías del sistema**: Datos de energía cinética, potencial y total.

Ejemplo de salida:
```
Time(s)    Kinetic Energy(J)    Potential Energy(J)    Total Energy(J)
1.00e-17    1.23e-20             -4.56e-19              -4.33e-19
...
```
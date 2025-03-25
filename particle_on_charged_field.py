# electron_field_simulation.py
# Simulation of the dynamics of an electron in an electric field generated by a charged ring

import numpy as np
import os

# Constantes
PI = np.pi
K_COULOMB = 8.9875517873681764e9  # Constante de Coulomb (N*m^2/C^2)
ELECTRON_MASS = 9.11e-31  # Masa del electrón (kg)
ELECTRON_CHARGE = -1.602176634e-19  # Carga del electrón (C)

class Ring:
    def __init__(self, n_charges, radius, total_charge):
        """
        Inicializa un anillo cargado.
        
        Parámetros:
        - n_charges: número de partículas con carga distribuidas en el anillo.
        - radius: radio del anillo.
        - total_charge: carga total en el anillo.
        """
        self.n_charges = n_charges
        self.radius = radius
        self.total_charge = total_charge
        self.charge_per_particle = total_charge / n_charges  # Carga por partícula en el anillo
        self.positions = self._compute_positions()  # Calcula las posiciones de las cargas en el anillo

    def _compute_positions(self):
        """Calcula las posiciones de las partículas en el anillo."""
        positions = []
        theta_step = 2 * PI / self.n_charges  # Paso angular para la distribución de carga
        for i in range(self.n_charges):
            theta = i * theta_step
            y = self.radius * np.cos(theta)  # Coordenada Y
            z = self.radius * np.sin(theta)  # Coordenada Z
            positions.append([0, y, z])  # Se asume que está en el plano YZ
        return np.array(positions)

class Electron:
    def __init__(self, position, velocity):
        """
        Inicializa un electrón con posición y velocidad.
        
        Parámetros:
        - position: posición inicial del electrón como lista o array [x, y, z].
        - velocity: velocidad inicial del electrón como lista o array [vx, vy, vz].
        """
        self.position = np.array(position, dtype=float)
        self.velocity = np.array(velocity, dtype=float)

class Simulation:
    def __init__(self, electron, ring, dt, total_time):
        """
        Inicializa la simulación.
        
        Parámetros:
        - electron: objeto Electron a simular.
        - ring: objeto Ring que genera el campo eléctrico.
        - dt: paso de tiempo de la simulación (en segundos).
        - total_time: tiempo total de la simulación (en segundos).
        """
        self.electron = electron
        self.ring = ring
        self.dt = dt
        self.total_time = total_time
        self.time = 0
        self.energy_data = []  # Almacena datos de energía en cada paso

    def _compute_electric_field(self):
        """
        Calcula el campo eléctrico total en la posición del electrón debido al anillo cargado.
        """
        total_field = np.zeros(3)
        for charge_position in self.ring.positions:
            r_vector = self.electron.position - charge_position
            distance = np.linalg.norm(r_vector)
            if distance == 0:
                continue  # Evita la singularidad si el electrón está en la misma posición que una carga
            field_magnitude = K_COULOMB * self.ring.charge_per_particle / distance**2
            field_vector = field_magnitude * (r_vector / distance)
            total_field += field_vector
        return total_field

    def _update_position_and_velocity(self):
        """
        Actualiza la posición y la velocidad del electrón usando el método de Euler.
        """
        electric_field = self._compute_electric_field()
        acceleration = (electric_field * ELECTRON_CHARGE) / ELECTRON_MASS
        self.electron.velocity += acceleration * self.dt
        self.electron.position += self.electron.velocity * self.dt

    def _compute_energies(self):
        """
        Calcula la energía cinética y potencial del electrón.
        """
        kinetic_energy = 0.5 * ELECTRON_MASS * np.linalg.norm(self.electron.velocity)**2
        potential_energy = 0
        for charge_position in self.ring.positions:
            r_vector = self.electron.position - charge_position
            distance = np.linalg.norm(r_vector)
            if distance != 0:
                potential_energy += ELECTRON_CHARGE * self.ring.charge_per_particle * K_COULOMB / distance
        total_energy = kinetic_energy + potential_energy
        return kinetic_energy, potential_energy, total_energy

    def run(self):
        """
        Ejecuta la simulación actualizando la posición, velocidad y calculando la energía en cada paso.
        """
        while self.time < self.total_time:
            self._update_position_and_velocity()
            energies = self._compute_energies()
            self.energy_data.append((self.time, *energies))
            self.time += self.dt
        self._save_results()

    def _save_results(self):
        """Guarda los datos de la simulación en el archivo 'dat/test_simulation.dat'."""
        os.makedirs("data", exist_ok=True)
        with open("data/test_simulation.dat", "w") as file:
            file.write("Tiempo(s)    Energía Cinética(J)    Energía Potencial(J)    Energía Total(J)\n")
            for data in self.energy_data:
                file.write(" ".join(f"{d:.5e}" for d in data) + "\n")

if __name__ == "__main__":
    # Configuración inicial
    electron = Electron(position=[0.15, 0.0, 0.0], velocity=[0.0, 0.0, 0.0])
    ring = Ring(n_charges=100, radius=0.1, total_charge=50e-9)
    simulation = Simulation(electron=electron, ring=ring, dt=1e-17, total_time=1e-5)

    # Ejecutar la simulación
    simulation.run()

    # Mensaje de finalización
    print("Simulación completada. Datos guardados en 'dat/test_simulation.dat'.")

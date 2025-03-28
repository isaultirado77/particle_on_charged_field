import os
import numpy as np
import argparse

class ChargedRingSystem:
    def __init__(self, n_charges, ring_radius, charge_ring, electron_position, electron_velocity, electron_mass=9.10938356e-31, electron_charge=-1.602176634e-19, k=8.9875517923e9):
        """
        Inicializa el sistema de anillo de cargas y electrón.
        
        :param n_charges: Número de cargas en el anillo.
        :param ring_radius: Radio del anillo (m).
        :param charge_ring: Carga de cada partícula en el anillo (C).
        :param electron_position: Posición inicial del electrón [x, y, z] (m).
        :param electron_velocity: Velocidad inicial del electrón [vx, vy, vz] (m/s).
        :param electron_mass: Masa del electrón (kg).
        :param electron_charge: Carga del electrón (C).
        :param k: Constante de Coulomb (N m²/C²).
        """
        self.n_charges = n_charges
        self.ring_radius = ring_radius
        self.charge_ring = charge_ring
        self.electron_mass = electron_mass
        self.electron_charge = electron_charge
        self.k = k
        
        # Inicializar estado del electrón: [x, y, z, vx, vy, vz]
        self.state = np.array([
            electron_position[0], electron_position[1], electron_position[2],
            electron_velocity[0], electron_velocity[1], electron_velocity[2]
        ])
        
        # Calcular posiciones de las cargas en el anillo (plano yz)
        self.ring_charges_pos = np.zeros((n_charges, 3))
        for i in range(n_charges):
            theta = 2 * np.pi * i / n_charges
            self.ring_charges_pos[i] = [0, ring_radius * np.cos(theta), ring_radius * np.sin(theta)]
    
    def electric_field(self, pos):
        """Calcula el campo eléctrico en la posición 'pos' debido al anillo de cargas."""
        Ex, Ey, Ez = 0.0, 0.0, 0.0
        x, y, z = pos
        
        for charge_pos in self.ring_charges_pos:
            dx = x - charge_pos[0]
            dy = y - charge_pos[1]
            dz = z - charge_pos[2]
            r_squared = dx**2 + dy**2 + dz**2
            r = np.sqrt(r_squared)
            
            if r_squared > 0:
                factor = self.k * self.charge_ring / (r_squared * r)
                Ex += factor * dx
                Ey += factor * dy
                Ez += factor * dz
        
        return np.array([Ex, Ey, Ez])
    
    def equations_of_motion(self, state):
        """Devuelve las derivadas del estado [dx/dt, dy/dt, dz/dt, dvx/dt, dvy/dt, dvz/dt]."""
        x, y, z, vx, vy, vz = state
        pos = np.array([x, y, z])
        
        # Campo eléctrico en la posición actual
        E = self.electric_field(pos)
        
        # Fuerza sobre el electrón: F = q_e * E
        ax = (self.electron_charge / self.electron_mass) * E[0]
        ay = (self.electron_charge / self.electron_mass) * E[1]
        az = (self.electron_charge / self.electron_mass) * E[2]
        
        return np.array([vx, vy, vz, ax, ay, az])
    
    def kinetic_energy(self):
        """Calcula la energía cinética del electrón."""
        vx, vy, vz = self.state[3], self.state[4], self.state[5]
        return 0.5 * self.electron_mass * (vx**2 + vy**2 + vz**2)
    
    def potential_energy(self):
        """Calcula la energía potencial del electrón."""
        x, y, z = self.state[0], self.state[1], self.state[2]
        pos = np.array([x, y, z])
        
        phi = 0.0  # Potencial eléctrico
        for charge_pos in self.ring_charges_pos:
            r = np.linalg.norm(pos - charge_pos)
            if r > 0:
                phi += self.k * self.charge_ring / r
        
        return self.electron_charge * phi
    
    def total_energy(self):
        """Calcula la energía total del electrón."""
        return self.kinetic_energy() + self.potential_energy()

class Simulator:
    def __init__(self, system, filename="electron_data"):
        """
        Inicializa el sistema de simulación.
        :param system: Instancia del sistema físico (e.g., ChargedRingSystem).
        :param filename: Nombre del archivo de salida.
        """
        self.system = system
        self.filename = filename
        os.makedirs("data", exist_ok=True)
    
    def runge_kutta_step(self, state, dt):
        """Realiza un paso de integración con RK4."""
        k1 = dt * self.system.equations_of_motion(state)
        k2 = dt * self.system.equations_of_motion(state + 0.5 * k1)
        k3 = dt * self.system.equations_of_motion(state + 0.5 * k2)
        k4 = dt * self.system.equations_of_motion(state + k3)
        return state + (k1 + 2*k2 + 2*k3 + k4) / 6
    
    def simulate(self, t_max, dt):
        """Ejecuta la simulación y guarda los datos."""
        steps = int(t_max / dt)
        state = self.system.state
        
        with open(f"data/{self.filename}.dat", "w") as file:
            file.write("# t x y z vx vy vz E_kin E_pot E_tot\n")
            time = 0.0
            for _ in range(steps):
                x, y, z, vx, vy, vz = state
                self.system.state = state  # Actualizar estado
                
                E_kin = self.system.kinetic_energy()
                E_pot = self.system.potential_energy()
                E_tot = self.system.total_energy()
                
                file.write(f"{time:.5f} {x:.5e} {y:.5e} {z:.5e} {vx:.5e} {vy:.5e} {vz:.5e} {E_kin:.5e} {E_pot:.5e} {E_tot:.5e}\n")
                state = self.runge_kutta_step(state, dt)
                time += dt
        
        print(f"Simulación completada. Datos guardados en data/{self.filename}.dat")

def main(n_charges=10, ring_radius=1.0, charge_ring=1.0e-9, 
         electron_position=[1.0, 0.0, 0.0], electron_velocity=[0.0, 1.0e4, 0.0], 
         dt=1.0e-9, t_max=1.0e-6, filename="electron_ring_simulation"):
    
    system = ChargedRingSystem(
        n_charges=n_charges,
        ring_radius=ring_radius,
        charge_ring=charge_ring,
        electron_position=electron_position,
        electron_velocity=electron_velocity
    )
    
    os.makedirs("data", exist_ok=True)
    with open(f"data/{filename}_params.txt", "w") as param_file:
        param_file.write("# Simulation parameters\n")
        param_file.write("{\n")
        param_file.write(f"    'n_charges': {n_charges},\n")
        param_file.write(f"    'ring_radius': {ring_radius},\n")
        param_file.write(f"    'charge_ring': {charge_ring},\n")
        param_file.write(f"    'electron_position': [{electron_position[0]}, {electron_position[1]}, {electron_position[2]}],\n")
        param_file.write(f"    'electron_velocity': [{electron_velocity[0]}, {electron_velocity[1]}, {electron_velocity[2]}],\n")
        param_file.write(f"    'dt': {dt},\n")
        param_file.write(f"    't_max': {t_max},\n")
        param_file.write(f"    'electron_mass': {system.electron_mass},\n")
        param_file.write(f"    'electron_charge': {system.electron_charge},\n")
        param_file.write(f"    'k': {system.k}\n")
        param_file.write("}\n")
    
    simulator = Simulator(system, filename)
    simulator.simulate(t_max, dt)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulador de electrón en un anillo de cargas usando RK4.")
    parser.add_argument("--n_charges", type=int, default=10, help="Número de cargas en el anillo")
    parser.add_argument("--ring_radius", type=float, default=1.0, help="Radio del anillo (m)")
    parser.add_argument("--charge_ring", type=float, default=1.0e-9, help="Carga de cada partícula en el anillo (C)")
    parser.add_argument("--electron_x", type=float, default=1.0, help="Posición x inicial del electrón (m)")
    parser.add_argument("--electron_vx", type=float, default=0.0, help="Velocidad x inicial del electrón (m/s)")
    parser.add_argument("--electron_vy", type=float, default=1.0e4, help="Velocidad y inicial del electrón (m/s)")
    parser.add_argument("--electron_vz", type=float, default=0.0, help="Velocidad z inicial del electrón (m/s)")
    parser.add_argument("--dt", type=float, default=1.0e-9, help="Paso de tiempo (s)")
    parser.add_argument("--t_max", type=float, default=1.0e-6, help="Tiempo total de simulación (s)")
    parser.add_argument("--filename", type=str, default="electron_ring_simulation", help="Nombre base para los archivos de salida")
    
    args = parser.parse_args()
    
    electron_position = [args.electron_x, 0.0, 0.0]
    electron_velocity = [args.electron_vx, args.electron_vy, args.electron_vz]
    
    main(
        n_charges=args.n_charges,
        ring_radius=args.ring_radius,
        charge_ring=args.charge_ring,
        electron_position=electron_position,
        electron_velocity=electron_velocity,
        dt=args.dt,
        t_max=args.t_max,
        filename=args.filename
    )
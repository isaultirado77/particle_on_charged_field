import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

def plot_3d_trajectory(n_charges, x, y, z, show=True, save=False, filename=""):
    """
    Grafica la trayectoria 3D del electrón y el anillo de protones
    Args:
        n_charges (int): Número de protones en el anillo
        x, y, z (arrays): Posiciones del electrón
        show (bool): Mostrar gráfico interactivo
        save (bool): Guardar imagen
        filename (str): Nombre base del archivo
    """
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    ax.plot(x, y, z, 'b-', linewidth=1.5, label='Trayectoria del electrón')
    ax.scatter(x[0], y[0], z[0], 'go', s=100, label='Posición inicial')
    ax.scatter(x[-1], y[-1], z[-1], 'ro', s=100, label='Posición final')
    
    theta = np.linspace(0, 2*np.pi, n_charges)
    ring_y = np.mean(y) * np.ones_like(theta)  # Centrar en el plano yz
    ring_z = np.mean(z) * np.ones_like(theta)
    
    ring_x = np.zeros_like(theta)
    ring_y = np.cos(theta) * np.max(np.sqrt(y**2 + z**2)) * 1.2  # Escalar para visualización
    ring_z = np.sin(theta) * np.max(np.sqrt(y**2 + z**2)) * 1.2
    
    ax.scatter(ring_x, ring_y, ring_z, 'co', s=50, label=f'Anillo ({n_charges} protones)')
    
    ax.set_xlabel('Posición X (m)')
    ax.set_ylabel('Posición Y (m)')
    ax.set_zlabel('Posición Z (m)')
    ax.set_title('Trayectoria del electrón')
    ax.legend()
    
    ax.view_init(elev=30, azim=45)

    if save:
        os.makedirs("plots", exist_ok=True)
        imname = f"_{filename}" if filename else ""
        plt.savefig(f"plots/positions{imname}.png", dpi=300, bbox_inches='tight')
        print(f"Gráfico guardado en plots/positions{imname}.png")
    
    if show:
        plt.show()
    
    plt.close()

def plot_energy(time, E_kin, E_pot, E_tot, show=True, save=False, filename=""):
    """
    Grafica las energías cinética, potencial y total del sistema.
    
    Args:
        time (array): Tiempo de la simulación en cada paso
        E_kin (array): Energía cinética en cada paso de tiempo
        E_pot (array): Energía potencial en cada paso de tiempo
        E_tot (array): Energía total en cada paso de tiempo
        show (bool): Mostrar gráfico interactivo
        save (bool): Guardar imagen
        filename (str): Nombre base del archivo
    """
    plt.figure(figsize=(10, 6))
    
    plt.plot(time, E_kin, 'b-', label='Energía Cinética', linewidth=1.5)
    plt.plot(time, E_pot, 'r-', label='Energía Potencial', linewidth=1.5)
    plt.plot(time, E_tot, 'g--', label='Energía Total', linewidth=1.5)
    
    plt.xlabel('Tiempo (s)')
    plt.ylabel('Energía (J)')
    plt.title('Evolución de las Energías del Sistema')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    
    y_min = min(min(E_kin), min(E_pot), min(E_tot)) * 1.1
    y_max = max(max(E_kin), max(E_pot), max(E_tot)) * 1.1
    plt.ylim(y_min, y_max)
    
    if save:
        os.makedirs("plots", exist_ok=True)
        imname = f"_{filename}" if filename else ""
        plt.savefig(f"plots/energias{imname}.png", dpi=300, bbox_inches='tight')
        print(f"Gráfico guardado en plots/energies{imname}.png")
    
    if show:
        plt.show()
    
    plt.close()
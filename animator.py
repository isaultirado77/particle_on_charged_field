import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import os

class ElectronTrajectoryAnimator:
    def __init__(self, data_path, n_charges=16, skip_frames=5):
        """
        Animador 3D de la trayectoria del electrón con cámara fija.
        
        Args:
            data_path (str): Ruta al archivo de datos de simulación
            n_charges (int): Número de cargas en el anillo
            skip_frames (int): Saltar frames para animación más ligera
        """
        self.data = np.loadtxt(data_path)
        self.t = self.data[:, 0]
        self.x = self.data[:, 1]
        self.y = self.data[:, 2]
        self.z = self.data[:, 3]
        
        self.n_charges = n_charges
        self.skip_frames = skip_frames
        
        self.fig = plt.figure(figsize=(12, 8))
        self.ax = self.fig.add_subplot(111, projection='3d')
        
        max_range = max(np.ptp(self.x), np.ptp(self.y), np.ptp(self.z)) * 0.5
        mid_x = (self.x.max() + self.x.min()) * 0.5
        mid_y = (self.y.max() + self.y.min()) * 0.5
        mid_z = (self.z.max() + self.z.min()) * 0.5
        
        self.ax.set_xlim(mid_x - max_range, mid_x + max_range)
        self.ax.set_ylim(mid_y - max_range, mid_y + max_range)
        self.ax.set_zlim(mid_z - max_range, mid_z + max_range)
        
        self.ax.view_init(elev=30, azim=45)  # 30° elevación, 45° azimut
        
        self.ax.set_xlabel('Posición X (m)')
        self.ax.set_ylabel('Posición Y (m)')
        self.ax.set_zlabel('Posición Z (m)')
        self.ax.set_title('Trayectoria del electrón en campo de anillo de cargas')
        
        # Dibujar anillo de cargas (plano YZ)
        theta = np.linspace(0, 2*np.pi, n_charges)
        ring_y = np.cos(theta) * max_range * 0.8
        ring_z = np.sin(theta) * max_range * 0.8
        self.ax.plot(np.zeros_like(theta), ring_y, ring_z, 'co-', 
                    markersize=5, alpha=0.5, label=f'Anillo ({n_charges} protones)')
        
        self.trajectory, = self.ax.plot([], [], [], 'b-', linewidth=1, alpha=0.7)
        self.electron, = self.ax.plot([], [], [], 'ro', markersize=8)
        self.time_text = self.ax.text2D(0.05, 0.95, '', transform=self.ax.transAxes)
        self.ax.legend(loc='upper right')
        
    def init_animation(self):
        """Inicializa los elementos de la animación."""
        self.trajectory.set_data([], [])
        self.trajectory.set_3d_properties([])
        self.electron.set_data([], [])
        self.electron.set_3d_properties([])
        return self.trajectory, self.electron
    
    def update(self, frame):
        """Actualiza la animación para cada frame sin cambiar la vista."""
        idx = frame * self.skip_frames
        if idx >= len(self.t):
            idx = len(self.t) - 1

        self.ax.view_init(elev=30, azim=45)

        self.trajectory.set_data(self.x[:idx], self.y[:idx])
        self.trajectory.set_3d_properties(self.z[:idx])

        self.electron.set_data([self.x[idx]], [self.y[idx]])
        self.electron.set_3d_properties([self.z[idx]])

        return self.trajectory, self.electron

    
    def animate(self, output_file="electron_trajectory.gif", fps=30):
        """Genera y guarda la animación."""
        os.makedirs("animations", exist_ok=True)
        
        total_frames = len(self.t) // self.skip_frames
        
        ani = animation.FuncAnimation(
            self.fig, self.update, frames=total_frames,
            init_func=self.init_animation, blit=True,
            interval=50, repeat=True
        )

        writer = animation.PillowWriter(fps=fps)
        ani.save(f"animations/{output_file}", writer=writer, dpi=100)
        print(f"Animación guardada en animations/{output_file}")
        
        plt.close()

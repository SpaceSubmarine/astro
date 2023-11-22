import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constantes
G = 6.67430e-11  # Constante gravitacional, m^3 kg^-1 s^-2
M_earth = 5.972e24  # Masa de la Tierra, kg
M_moon = 7.342e22   # Masa de la Luna, kg
R = 3.844e8  # Distancia media Tierra-Luna, m

# Crear una rejilla 3D
size = 5e8
num_points = 100
x = np.linspace(-size, size, num_points)
y = np.linspace(-size, size, num_points)
z = np.linspace(-size, size, num_points)
X, Y, Z = np.meshgrid(x, y, z)

# Calcular el campo gravitatorio debido a la Tierra y la Luna
R_earth = np.sqrt(X**2 + Y**2 + Z**2)
R_moon = np.sqrt((X-R)**2 + Y**2 + Z**2)

# Ley de gravitación universal
g_earth = G * M_earth / R_earth**3
g_moon = G * M_moon / R_moon**3

# Componentes del campo gravitatorio
gX_earth = g_earth * X / R_earth
gY_earth = g_earth * Y / R_earth
gZ_earth = g_earth * Z / R_earth

gX_moon = g_moon * (X-R) / R_moon
gY_moon = g_moon * Y / R_moon
gZ_moon = g_moon * Z / R_moon

# Superponer los campos gravitatorios
gX_total = gX_earth + gX_moon
gY_total = gY_earth + gY_moon
gZ_total = gZ_earth + gZ_moon

# Gráfica
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='3d')
ax.quiver(X, Y, Z, gX_total, gY_total, gZ_total, length=5e6, normalize=True, color='blue', linewidth=0.5)
ax.set_xlim([-size, size])
ax.set_ylim([-size, size])
ax.set_zlim([-size, size])
ax.set_title("Campo gravitatorio superpuesto de la Tierra y la Luna")
plt.show()

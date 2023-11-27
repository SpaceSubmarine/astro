import numpy as np
import matplotlib.pyplot as plt


# Función para calcular las coordenadas (x, y, z) a partir de los elementos orbitales
def orbit_coordinates(a, e, i, omega, w, theta):
    # Convertir ángulos de grados a radianes
    i = np.radians(i)
    omega = np.radians(omega)
    w = np.radians(w)
    theta = np.radians(theta)

    # Calcular el argumento del periastro
    argument_of_periapsis = omega + w

    # Calcular la distancia radial
    r = (a * (1 - e**2)) / (1 + e * np.cos(theta))

    # Calcular las coordenadas en el plano orbital
    x_orbital = r * np.cos(theta)
    y_orbital = r * np.sin(theta)

    # Transformar a coordenadas 3D
    x = (np.cos(argument_of_periapsis) * np.cos(omega) - np.sin(argument_of_periapsis) * np.sin(omega) * np.cos(i)) * x_orbital + (-np.sin(argument_of_periapsis) * np.cos(omega) - np.cos(argument_of_periapsis) * np.sin(omega) * np.cos(i)) * y_orbital
    y = (np.cos(argument_of_periapsis) * np.sin(omega) + np.sin(argument_of_periapsis) * np.cos(omega) * np.cos(i)) * x_orbital + (-np.sin(argument_of_periapsis) * np.sin(omega) + np.cos(argument_of_periapsis) * np.cos(omega) * np.cos(i)) * y_orbital
    z = (np.sin(argument_of_periapsis) * np.sin(i)) * x_orbital + (np.cos(argument_of_periapsis) * np.sin(i)) * y_orbital

    return x, y, z

# Definir elementos orbitales de una órbita Molniya
a = 26600  # Semi-major axis in km
e = 0.74  # Eccentricity
i = 63.4  # Inclination in degrees
omega = 0  # Right ascension of the ascending node in degrees
w = 270  # Argument of periapsis in degrees
theta_range = np.linspace(0, 360, 500)  # Range of true anomaly

# Calcular coordenadas de la órbita
x_orbit, y_orbit, z_orbit = [], [], []
for theta in theta_range:
    x, y, z = orbit_coordinates(a, e, i, omega, w, theta)
    x_orbit.append(x)
    y_orbit.append(y)
    z_orbit.append(z)

# Crear gráfico 3D
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
# Dibujar órbita
ax.plot(x_orbit, y_orbit, z_orbit, label='Molniya Orbit')
# Configuración del gráfico
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Z (km)')
ax.set_title('3D Visualization of a Molniya Orbit')
ax.legend()
plt.show()

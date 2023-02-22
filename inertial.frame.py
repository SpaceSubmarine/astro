import numpy as np
import matplotlib.pyplot as plt
import astro_constants

# INPUT SECTION
n = 5  # Number of points along each dimension of the 3D space
mars_d = astro_constants.distances['Mars']  # Distance from the Sun to Mars
mars_m = astro_constants.masses['Mars']  # Mass of Mars
earth_m = astro_constants.masses['Earth']  # Mass of Earth
earth_d = astro_constants.distances['Earth']  # Distance from the Sun to Earth
sun_m = astro_constants.masses['Sun']  # Mass of the Sun
max_distance = mars_d  # Maximum distance from the center of the 3D space
spacing = max_distance * 2 / (n - 1)  # Spacing between adjacent points along each dimension

# Meshgrid
x_grid = np.linspace(-max_distance, max_distance, n)  # X coordinates of the points in the 3D space
y_grid = np.linspace(-max_distance, max_distance, n)  # Y coordinates of the points in the 3D space
z_grid = np.linspace(-max_distance, max_distance, n)  # Z coordinates of the points in the 3D space

x, y, z = np.meshgrid(x_grid, y_grid, z_grid, indexing='ij')  # Create a 3D grid of points

space_3d = np.zeros((n, n, n))  # Initialize a 3D array to represent the 3D space
masses_3d = np.zeros_like(space_3d)  # Initialize a 3D array to represent the masses at each point in the 3D space
masses_3d[int(n/2), int(n/2), int(n/2)] = sun_m  # Set the mass of the Sun at the center of the 3D space
masses_3d[int((n/2)+(earth_d/spacing)), int(n/2), int(n/2)] = earth_m  # Set the mass of the Earth

# PLOT
plt.style.use('dark_background')
fig = plt.figure()  # Create a new figure
ax = fig.add_subplot(111, projection='3d')  # Add a 3D axis to the figure

# Plot the masses_3d using scatter
mask = masses_3d != 0  # Create a boolean mask for the points with non-zero mass
sc = ax.scatter(x, y, z, c=masses_3d)  # Create a scatter plot of the points with color mapped to their mass

# Apply mask to the plot
alpha = np.zeros_like(masses_3d)  # Create an array to represent the transparency of the points
alpha[mask] = 1  # Set the transparency to 1 for the points with non-zero mass
alpha[~mask] = 0  # Set the transparency to 0 for the points with zero mass
sc.set_alpha(alpha)  # Apply the transparency to the scatter plot

# Set the limits of the plot
ax.set_xlim(-max_distance, max_distance)  # Set the limits of the x-axis
ax.set_ylim(-max_distance, max_distance)  # Set the limits of the y-axis
ax.set_zlim(-max_distance, max_distance)  # Set the limits of the z-axis

# Show the colorbar
fig.colorbar(sc)  # Add a colorbar to the figure to show the mapping of mass to color

plt.show()  # Show the plot

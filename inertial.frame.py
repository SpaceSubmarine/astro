import numpy as np
import matplotlib.pyplot as plt
from astrod import astro_plots, astro_constants


# INPUT SECTION
n = 10
mars_d = astro_constants.distances['Mars']
mars_m = astro_constants.masses['Mars']
sun_m = astro_constants.masses['Sun']
max_distance = mars_d
spacing = max_distance * 2 / (n - 1)

# Meshgrid
x_grid = np.linspace(-max_distance, max_distance, n)
y_grid = np.linspace(-max_distance, max_distance, n)
z_grid = np.linspace(-max_distance, max_distance, n)
void = np.zeros_like(x_grid)
x, y, z = np.meshgrid(x_grid, y_grid, z_grid, indexing='ij')

space_3d = np.zeros((n, n, n))
x_void, y_void, z_void = np.meshgrid(void, void, z_grid, indexing='ij')

# PLOT
ax = astro_plots.create_3d_plot(n, max_distance)
plt.show()

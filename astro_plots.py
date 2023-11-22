import numpy as np
import matplotlib.pyplot as plt

def create_3d_plot(size, max_distance):
    # Calculate the spacing between points in the matrix
    spacing = max_distance * 2 / (size - 1)

    # Create a meshgrid of x, y, and z values
    x_vals = np.linspace(-max_distance, max_distance, size)
    y_vals = np.linspace(-max_distance, max_distance, size)
    z_vals = np.linspace(-max_distance, max_distance, size)
    x, y, z = np.meshgrid(x_vals, y_vals, z_vals, indexing='ij')

    # Create the 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Set the limits of the plot
    ax.set_xlim(-max_distance, max_distance)
    ax.set_ylim(-max_distance, max_distance)
    ax.set_zlim(-max_distance, max_distance)

    return ax
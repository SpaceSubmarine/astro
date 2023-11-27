import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

"""
Reaction wheels simulation 
"""

# Aluminum disc properties (reaction wheel)
disc_diameter = 0.05  # meters
disc_thickness = 0.003  # meters
aluminum_dens = 2700  # kg/m^3

# Volume and mass calculation
disc_volume = np.pi * (disc_diameter / 2) ** 2 * disc_thickness
disc_mass = aluminum_dens * disc_volume

# Disc Inertia (I = 1/2 *m r^2 for solid disc)
disc_radi = disc_diameter / 2
disc_inertia = 0.5 * disc_mass * disc_radi ** 2

# BLDC properties
acceleration_time = 2  # seconds
breaking_time = 2  # seconds

# Assuming constant acceleration during acceleration time and braking
# Assuming max angular velocity at the end of the acceleration
v_max_angular = 2200 * 2 * np.pi / 60
angular_acc = v_max_angular / acceleration_time
# Transient simulation
dt = 0.01  # delta T
total_time = 2 * acceleration_time + breaking_time
time_list = np.arange(0, total_time + dt, dt)

# Initializing variables
v_angular = np.zeros_like(time_list)
accelerations = np.zeros_like(time_list)
angular_momentum = np.zeros_like(time_list)

# simulation
for i, t in enumerate(time_list):
    print("Current time: ", t)
    if t < acceleration_time:
        # Accelerating
        accelerations[i] = angular_acc
        v_angular[i] = accelerations[i] * t
    elif t < acceleration_time + breaking_time:
        # Constant velocity
        accelerations[i] = 0
        v_angular[i] = v_max_angular
    else:
        accelerations[i] = -angular_acc
        v_angular[i] = v_max_angular + (accelerations[i] * (t - acceleration_time - breaking_time))
    angular_momentum[i] = disc_inertia * v_angular[i]


def plot_color_curve(time_list, v_angular, accelerations, angular_momentum):
    plt.style.use('dark_background')
    fig, axes = plt.subplots(3, 1, figsize=(12, 8))

    # Background color configuration
    for ax in axes:
        ax.set_facecolor((0.08, 0.08, 0.1))  # Custom background color for each subplot
        ax.grid(color='grey', linestyle='--', linewidth=0.5, alpha=0.6)

    fig.patch.set_facecolor((0.05, 0.05, 0.05))  # Custom background color for the figure

    # Subplot 1: Angular Velocity
    axes[0].plot(time_list, v_angular, label="Angular Velocity", color="cyan")
    axes[0].set_xlabel("Time (s)")
    axes[0].set_ylabel("Angular Velocity (rad/s)")
    axes[0].legend()

    # Subplot 2: Angular Acceleration
    axes[1].plot(time_list, accelerations, label="Angular Acceleration", color="#FF0000")
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("Angular Acceleration (rad/s²)")
    axes[1].legend()

    # Subplot 3: Angular Momentum
    axes[2].plot(time_list, angular_momentum, label="Angular Momentum", color="lime")
    axes[2].set_xlabel("Time (s)")
    axes[2].set_ylabel("Angular Momentum (kg·m²/s)")
    axes[2].legend()

    plt.tight_layout()
    plt.show()


# Running the function with the simulation data
# (These data should be defined in the previous cells)
plot_color_curve(time_list, v_angular, accelerations, angular_momentum)


##################################################################################################
# Re-defining the initial angular velocities of the CubeSat before simulation
omega_x = 0.0  # rad/s
omega_y = 0.0  # rad/s
omega_z = 0.0  # rad/s

# Initialization of the Euler angles and angular velocities of the CubeSat
pitch = 0.0  # rotation around the X-axis
roll = 0.0   # rotation around the Y-axis
yaw = 0.0    # rotation around the Z-axis

# CubeSat 1U data
cubesat_mass = 1.0  # kg
cubesat_side = 0.1  # meters (a 1U CubeSat has dimensions of 10cm x 10cm x 10cm)

# Calculation of the CubeSat's inertia tensor (modeled as a solid cube)
cubesat_inertia_x = 1/6 * cubesat_mass * (cubesat_side**2 + cubesat_side**2)
cubesat_inertia_y = cubesat_inertia_x  # Symmetry in a cube
cubesat_inertia_z = cubesat_inertia_x  # Symmetry in a cube

# Total inertia tensor (CubeSat + reaction wheel)
total_inertia_tensor = np.array([[cubesat_inertia_x, 0, 0],
                                 [0, cubesat_inertia_y + disc_inertia, 0],
                                 [0, 0, cubesat_inertia_z]])

# Adjusting the simulation to focus torque on the Y-axis and compare with other axes
# Resuming the simulation with defined angular velocities
pitch_list = []
roll_list = []
yaw_list = []

# Resuming the simulation taking into account the adjustments
for i in range(len(time_list)):
    t = time_list[i]
    dt = time_list[i] - time_list[i-1] if i > 0 else 0

    # Calculating the torque induced by the reaction wheel
    torque_y = disc_inertia * accelerations[i]

    # Updating the CubeSat's angular velocities
    # Only the Y-axis (roll) will be significantly affected by the torque
    omega_y += (torque_y / total_inertia_tensor[1, 1]) * dt

    # Updating the CubeSat's Euler angles (orientation)
    pitch += omega_x * dt
    roll += omega_y * dt
    yaw += omega_z * dt

    # Storing the results
    pitch_list.append(pitch)
    roll_list.append(roll)
    yaw_list.append(yaw)

# Displaying the results
plt.figure(figsize=(12, 8))
plt.subplot(3, 1, 1)
plt.plot(time_list, pitch_list, label="Pitch", color="cyan")
plt.xlabel("Time (s)")
plt.ylabel("Pitch (rad)")
plt.legend()

plt.subplot(3, 1, 2)
plt.plot(time_list, roll_list, label="Roll", color="lime")
plt.xlabel("Time (s)")
plt.ylabel("Roll (rad)")
plt.legend()

plt.subplot(3, 1, 3)
plt.plot(time_list, yaw_list, label="Yaw", color="yellow")
plt.xlabel("Time (s)")
plt.ylabel("Yaw (rad)")
plt.legend()

plt.tight_layout()
plt.show()








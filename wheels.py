import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

"""
GYRO-STATOR / REACTION WHEELS SIMULATION 

IMPORTANT:
X-AXIS IS DEFINED AS THE PITCH AXIS OF ROTATION OF THE CUBESAT
Y-AXIS IS DEFINED AS THE ROLL AXIS OF ROTATION OF THE CUBESAT
Z-AXIS IS DEFINED AS THE YAW AXIS OF ROTATION OF THE CUBESAT

REACTION WHEEL

"""

##################################################################################################
""" INPUT/VARIABLES SECTION """

# Re-defining the initial angular velocities of the CubeSat before simulation
omega_x = 0.0  # rad/s
omega_y = 0.0  # rad/s
omega_z = 0.0  # rad/s

# Initialization of the Euler angles and angular velocities of the CubeSat
pitch = 0.0  # rotation around the X-axis
roll = 0.0  # rotation around the Y-axis
yaw = 0.0  # rotation around the Z-axis

# Aluminum disc properties (reaction wheel)
disc_diameter = 0.05  # meters
disc_thickness = 0.003  # meters
aluminum_dens = 2700  # kg/m^3

# CubeSat 1U data
cubesat_mass = 1.0  # kg
cubesat_side = 0.1  # meters (a 1U CubeSat has dimensions of 10cm x 10cm x 10cm)

##################################################################################################
""" FUNCTION DEFINITION SECTION """


def plot_combined(time_list, omega_data, alpha_data, momentum_data, dynamics_data,labels, colors, ylabels, title):
    plt.style.use('dark_background')
    fig, axes = plt.subplots(4, 1, figsize=(12, 8))

    fig.suptitle(title)
    for ax in axes:
        ax.set_facecolor((0.08, 0.08, 0.1))
        ax.grid(color='grey', linestyle='--', linewidth=0.5, alpha=0.6)

    fig.patch.set_facecolor((0.05, 0.05, 0.05))

    # Velocidad angular en el primer subplot
    axes[0].plot(time_list, omega_data[0], label=labels[0][0], color=colors[0])
    axes[0].plot(time_list, omega_data[1], label=labels[0][1], color=colors[1])
    axes[0].plot(time_list, omega_data[2], label=labels[0][2], color=colors[2])
    axes[0].set_xlabel("Time (s)")
    axes[0].set_ylabel(ylabels[0])
    axes[0].legend()

    # Aceleración angular en el segundo subplot
    axes[1].plot(time_list, alpha_data[0], label=labels[1][0], color=colors[0])
    axes[1].plot(time_list, alpha_data[1], label=labels[1][1], color=colors[1])
    axes[1].plot(time_list, alpha_data[2], label=labels[1][2], color=colors[2])
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel(ylabels[1])
    axes[1].legend()

    # Momento angular en el tercer subplot
    axes[2].plot(time_list, momentum_data[0], label=labels[2][0], color=colors[0])
    axes[2].plot(time_list, momentum_data[1], label=labels[2][1], color=colors[1])
    axes[2].plot(time_list, momentum_data[2], label=labels[2][2], color=colors[2])
    axes[2].set_xlabel("Time (s)")
    axes[2].set_ylabel(ylabels[2])
    axes[2].legend()

    axes[3].plot(time_list, dynamics_data[0], label=labels[3][0], color=colors[0])
    axes[3].plot(time_list, dynamics_data[1], label=labels[3][1], color=colors[1])
    axes[3].plot(time_list, dynamics_data[2], label=labels[3][2], color=colors[2])
    axes[3].set_xlabel("Time (s)")
    axes[3].set_ylabel(ylabels[3])
    axes[3].legend()

    plt.tight_layout()



##################################################################################################
""" INITIAL CALCULATIONS """

# Volume and mass calculation
disc_volume = np.pi * (disc_diameter / 2) ** 2 * disc_thickness
disc_mass = aluminum_dens * disc_volume

# Disc Inertia (I = 1/2 *m r^2 for solid disc)
disc_radi = disc_diameter / 2
disc_inertia_x = (1/4) * disc_mass * disc_radi**2 + (1/12) * disc_mass * disc_thickness**2
disc_inertia_y = disc_inertia_x
disc_inertia_z = 0.5 * disc_mass * disc_radi ** 2

# BLDC properties
acceleration_time = 2  # seconds
breaking_time = 2  # seconds

# Assuming constant acceleration during acceleration time and braking
# Assuming max angular velocity at the end of the acceleration
v_max_angular = 2200 * 2 * np.pi / 60  # rpm to rad/s
angular_acc_x = 0
angular_acc_y = 0
angular_acc_z = v_max_angular / acceleration_time

# Transient simulation
dt = 0.01  # delta T
total_time = 2 * acceleration_time + breaking_time
time_list = np.arange(0, total_time + dt, dt)

# Initializing variables
omega_disc_x = np.zeros_like(time_list)
omega_disc_y = np.zeros_like(time_list)
omega_disc_z = np.zeros_like(time_list)
alpha_disc_x = np.zeros_like(time_list)
alpha_disc_y = np.zeros_like(time_list)
alpha_disc_z = np.zeros_like(time_list)
angular_momentum_x = np.zeros_like(time_list)
angular_momentum_y = np.zeros_like(time_list)
angular_momentum_z = np.zeros_like(time_list)

##################################################################################################
""" REACTIONS SIMULATION """

for i, t in enumerate(time_list):
    if t < acceleration_time:
        # Accelerating
        alpha_disc_x[i] = angular_acc_x
        alpha_disc_y[i] = angular_acc_y
        alpha_disc_z[i] = angular_acc_z
        omega_disc_x[i] = alpha_disc_x[i] * t
        omega_disc_y[i] = alpha_disc_y[i] * t
        omega_disc_z[i] = alpha_disc_z[i] * t

    elif t < acceleration_time + breaking_time:
        # Constant velocity
        alpha_disc_x[i] = 0
        alpha_disc_y[i] = 0
        alpha_disc_z[i] = 0
        omega_disc_x[i] = omega_disc_x[i - 1]
        omega_disc_y[i] = omega_disc_y[i - 1]
        omega_disc_z[i] = omega_disc_z[i - 1]
    else:
        # Decelerating
        time_since_braking = t - (acceleration_time + breaking_time)
        alpha_disc_x[i] = -angular_acc_x
        alpha_disc_y[i] = -angular_acc_y
        alpha_disc_z[i] = -angular_acc_z
        omega_disc_x[i] = omega_disc_x[i - 1] + alpha_disc_x[i] * dt
        omega_disc_y[i] = omega_disc_y[i - 1] + alpha_disc_y[i] * dt
        omega_disc_z[i] = omega_disc_z[i - 1] + alpha_disc_z[i] * dt

    angular_momentum_x[i] = disc_inertia_x * omega_disc_x[i]
    angular_momentum_y[i] = disc_inertia_y * omega_disc_y[i]
    angular_momentum_z[i] = disc_inertia_z * omega_disc_z[i]

##################################################################################################
""" GYRO-DYNAMICS SIMULATION"""

# Calculation of the CubeSat's inertia tensor (modeled as a solid cube)
cubesat_inertia_x = 1 / 6 * cubesat_mass * (cubesat_side ** 2 + cubesat_side ** 2)
cubesat_inertia_y = cubesat_inertia_x  # Symmetry in a cube
cubesat_inertia_z = cubesat_inertia_x  # Symmetry in a cube

# Total inertia tensor (CubeSat + reaction wheel)
total_inertia_tensor = np.array([[cubesat_inertia_x + disc_inertia_x, 0, 0],
                                 [0, cubesat_inertia_y + disc_inertia_y, 0],
                                 [0, 0, cubesat_inertia_z + disc_inertia_z]])


# Adjusting the simulation to focus torque on the Y-axis and compare with other axes
# Resuming the simulation with defined angular velocities
pitch_list = []
roll_list = []
yaw_list = []

# Resuming the simulation taking into account the adjustments
for i in range(len(time_list)):
    t = time_list[i]
    dt = time_list[i] - time_list[i - 1] if i > 0 else 0

    # Calculating the torque induced by the reaction wheel
    torque_x = disc_inertia_x * alpha_disc_x[i]
    torque_y = disc_inertia_y * alpha_disc_y[i]
    torque_z = disc_inertia_z * alpha_disc_z[i]

    # Updating the CubeSat's angular velocities
    # Only the Y-axis (roll) will be significantly affected by the torque
    omega_x += (torque_x / total_inertia_tensor[0, 0]) * dt
    omega_y += (torque_y / total_inertia_tensor[1, 1]) * dt
    omega_z += (torque_z / total_inertia_tensor[2, 2]) * dt

    # Updating the CubeSat's Euler angles (orientation)
    pitch += omega_x * dt
    roll += omega_y * dt
    yaw += omega_z * dt

    # Storing the results
    pitch_list.append(pitch)
    roll_list.append(roll)
    yaw_list.append(yaw)

##################################################################################################
""" PLOT SECTION """

# Plot Data
omega_data = [omega_disc_x, omega_disc_y, omega_disc_z]
alpha_data = [alpha_disc_x, alpha_disc_y, alpha_disc_z]
momentum_data = [angular_momentum_x, angular_momentum_y, angular_momentum_z]
gyro_data = [pitch_list, roll_list, yaw_list]

# Labels and Colors
labels_reaction = [["Omega X", "Omega Y", "Omega Z"],
                   ["Alpha X", "Alpha Y", "Alpha Z"],
                   ["Momentum X", "Momentum Y", "Momentum Z"],
                   ["Pitch", "Roll", "Yaw"]]


colors = ["red", "lime", "dodgerblue"]

# Plotting
plot_combined(time_list, omega_data, alpha_data, momentum_data, gyro_data, labels_reaction, colors,
              ["Angular Velocity (rad/s)", "Angular Acceleration (rad/s²)", "Angular Momentum (kg·m²/s)",
               "Rotation (rad)"],
              "Reaction Wheels Dynamics")

plt.show()

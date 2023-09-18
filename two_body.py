import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

earth_radius = 6378  # km
earth_mu = 398600  # km^3/s^2

'''
This code is extracted from the Youtube tutorial 
The Two Body Problem / ODE Solvers | Orbital Mechanics with Python 2
by Alfonso Gonzalez 
url = https://www.youtube.com/watch?v=7JY44m6eemo&list=PLOIRBaljOV8gn074rWFWYP1dCr2dJqWab&index=52&t=997s&ab_channel=AlfonsoGonzalez-Astrodynamics%26SEPodcast
'''


def plot():
    plt.style.use('dark_background')
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Draw the Earth as a sphere
    u, v = np.mgrid[0:2 * np.pi:40j, 0:np.pi:20j]
    x = earth_radius * np.cos(u) * np.sin(v)
    y = earth_radius * np.sin(u) * np.sin(v)
    z = earth_radius * np.cos(v)
    ax.plot_surface(x, y, z, color="blue", alpha=0.3)  # The color and alpha are adjustable
    # ax.grid(False)
    # ax.set_axis_off()  # This deletes axis
    # Draw the satellite's trajectory
    ax.plot(rs[:, 0], rs[:, 1], rs[:, 2], color="red", linewidth=1)
    plt.savefig("3D_orbit.png")
    plt.show()


def diff_eq(t, y, mu):  # our differential equation
    # unpack state
    # "t" is the time
    # "y" is the position and velocity
    # "mu" is the gravitational parameter
    rx, ry, rz, vx, vy, vz = y

    # for convenience, we define the position array as a position vector
    r = np.array([rx, ry, rz])
    # we also want the norm of the radius vector
    norm_r = np.linalg.norm(r)
    # Gravitational acceleration for the three components
    # with negative sign as the acceleration is pointing to the
    # center of the attracting body
    ax, ay, az = -r * mu / norm_r ** 3

    return [vx, vy, vz, ax, ay, az]


if __name__ == '__main__':
    # initial conditions of orbit parameters
    r_mag = earth_radius + 500  # km
    v_mag = np.sqrt(earth_mu / r_mag)

    # initial position and velocity vectors
    # for a circular orbit, the inclination is zero
    r0 = [r_mag, 0, 0]
    v0 = [0, v_mag, 0]

    # timespan
    tspan = 100 * 60  # minutes * seconds

    # timestep
    dt = 100  # seconds

    # total number of steps
    n_steps = int(np.ceil(tspan / dt))

    # initialize arrays (pre allocating memory)
    ys = np.zeros((n_steps, 6))
    ts = np.zeros((n_steps, 1))
    y0 = r0 + v0

    ys[0] = np.array(y0)  # initial condition
    step = 1

    """
    Here the solver 'ode' initializes with the function 'diff_eq', 
    which calculates the time derivatives from the current state 'y' and 'mu'
    """
    solver = ode(diff_eq)
    solver.set_integrator('lsoda')  # Integration algorithm
    solver.set_initial_value(y0, 0)  # Initial conditions for t=0
    solver.set_f_params(earth_mu)  # Introducing the gravitational parameter

    # Orbit propagation
    while solver.successful() and step < n_steps:
        solver.integrate(solver.t + dt)
        ts[step] = solver.t
        ys[step] = solver.y
        step += 1

    # extracting values of position from 'ys'
    rs = ys[:, :3]

    plot()

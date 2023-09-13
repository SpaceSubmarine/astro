import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

earth_radius = 6378  # km
earth_mu = 398600  # km^3/s^2


def diff_eq(t, y, mu):
    rx, ry, rz, vx, vy, vz = y
    r = np.array([rx, ry, rz])
    norm_r = np.linalg.norm(r)
    ax, ay, az = -r * mu / norm_r ** 3
    return [vx, vy, vz, ax, ay, az]


if __name__ == '__main__':
    r_mag = earth_radius + 500  # km
    v_mag = np.sqrt(earth_mu / r_mag)
    r0 = [r_mag, 0, 0]
    v0 = [0, v_mag, 0]
    tspan = 100 * 60  # seconds
    dt = 100  # seconds
    n_steps = int(np.ceil(tspan / dt))
    ys = np.zeros((n_steps, 6))
    ts = np.zeros((n_steps, 1))
    y0 = r0 + v0
    ys[0] = np.array(y0)
    step = 1

    solver = ode(diff_eq)
    solver.set_integrator('lsoda')
    solver.set_initial_value(y0, 0)
    solver.set_f_params(earth_mu)

    while solver.successful() and step < n_steps:
        solver.integrate(solver.t + dt)
        ts[step] = solver.t
        ys[step] = solver.y
        step += 1

    rs = ys[:, :3]

    plt.style.use('dark_background')
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Draw the Earth as a sphere
    u, v = np.mgrid[0:2 * np.pi:40j, 0:np.pi:20j]
    x = earth_radius * np.cos(u) * np.sin(v)
    y = earth_radius * np.sin(u) * np.sin(v)
    z = earth_radius * np.cos(v)
    ax.plot_surface(x, y, z, color="blue", alpha=0.5)  # The color and alpha are adjustable
    ax.grid(False)
    ax.set_axis_off()  # Esto quita los ejes
    # Draw the satellite's trajectory
    ax.plot(rs[:, 0], rs[:, 1], rs[:, 2], color="white", linewidth=2)

    plt.show()
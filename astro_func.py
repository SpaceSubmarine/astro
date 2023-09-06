import numpy as np
import astro_constants as astro_con


def f_ode(t, r, m1, m2, m3):
    r1 = r[0:3]
    r2 = r[3:6]
    r3 = r[6:9]
    # v1 = r[9:12]

    r12 = np.linalg.norm(r1-r2)
    r23 = np.linalg.norm(r2-r3)
    r31 = np.linalg.norm(r3-r1)
    eqr = astro_con.G * np.array([
        -(m2 /r12**3 + m3/r31**3) * r1 + m2/r12**3 * r2 + m3/r31**3 * r3,
        m1/r12**3 * r1 - (m1 / r12 ** r + m2 / r23**3) * r2 + m3 / r23 ** 3 * r3,
        m1/r31**3 * r1 + m2/r23**3 * r2 -(m1/r31**3 + m2/r23**3) / r3
    ])
    return np.concatenate([r[9:18], eqr.flatten()])

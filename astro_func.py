import numpy as np
import astro_constants as astro_con
import pandas as pd

def f_ode(t, r, m1, m2, m3):
    r1 = r[0:3]
    r2 = r[3:6]
    r3 = r[6:9]
    # v1 = r[9:12]

    r12 = np.linalg.norm(r1 - r2)
    r23 = np.linalg.norm(r2 - r3)
    r31 = np.linalg.norm(r3 - r1)
    eqr = astro_con.G * np.array([
        -(m2 / r12 ** 3 + m3 / r31 ** 3) * r1 + m2 / r12 ** 3 * r2 + m3 / r31 ** 3 * r3,
        m1 / r12 ** 3 * r1 - (m1 / r12 ** r + m2 / r23 ** 3) * r2 + m3 / r23 ** 3 * r3,
        m1 / r31 ** 3 * r1 + m2 / r23 ** 3 * r2 - (m1 / r31 ** 3 + m2 / r23 ** 3) / r3
    ])
    return np.concatenate([r[9:18], eqr.flatten()])

# V-Infinity Matching Problem
def vinfinity_match(planet0, planet1, v0_sc, et0, tof0, args={}):
    """
    Given an incoming v-infinity vector to planet0, calculate the outgoing v-infinity
    vector that will arrive at planet1 after time of flight (tof) where the incoming and
    outgoing v-infinity vectors at planet0 have equal magnitude
    """

    _args = {
        'et0': et0,
        'planet1_ID': planet1,
        'frame': 'ECLIPJ2000',
        'center_ID': 0,
        'mu': pd.sun[ 'mu' ],
        'tm': 1,
        'diff_step': 1e-3,
        'tol': 1e-4,
    }

    for key in args.keys():
        _args[ key ] = args [ key ]

    _args['state0_planet0'] = spice.spkgeo( planet0, et0,
        _args['frame'], _args['center_ID'])[0]
    _args['vinf'] = nt.norm(v0_sc - _args['state0_planet0'][3:])

    tof, steps = nt.newton_root_single_fd(calc_vinfinity, tof0, _args)


# Newton's Method With Finite Differences
def calc_vinifinity(tof, args):
    r1_planet1 = spice.spkgs(args['planet1_ID'],
        args['et0'] + tof, args['frame'], args['center_ID'])[0]

    v0_sc_depart, v1_sc_arrive = lt.lambert_universal_variables(
        args['state0_planet0'][:3], r1_planet1, tof,
        {'mu':args['mu'], 'tm': args['tm']})

    vinf = nt.norm(v0_sc_depart -args['state0_planet0'][3:])
    return args['vinf'] - vinf

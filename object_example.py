import astro_constants


class Planet():
    """
    mass (kg)
    radius (m)
    sun_distance (km)
    """

    # Initial state of the object:
    def __init__(self, planet_name, M, R, sun_D):
        self.name = str(planet_name)
        self.mass = M
        self.radius = R
        self.sun_distance = sun_D


# Radius of all bodies missing in astro_constants.py
Planet_4 = Planet('Mars', astro_constants.masses["Mars"], 1000, astro_constants.distances["Mars"])

print("Mass of " + str(Planet_4.name) + ": " + str(Planet_4.mass) + "kg")
print("Distance to the Sun of " + str(Planet_4.name) + ": " + str(Planet_4.sun_distance) + "km")

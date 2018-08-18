from amuse.lab import Particles, units, ph4, nbody_system, new_plummer_model
from amuse.units import constants
import matplotlib.pyplot as plt

"""
def cluster(mass, radius, N):
    particles = Particles(1)
    conv = nbody_system.nbody_to_si(radius, mass)
    particles[0] = new_plummer_model(N, convert_nbody=conv, do_scale=True)
    pm = particles[0]
    pm.positon = (100, 100, 100) | units.parsec
    pm.velocity = (-1000, -1000, -1000) | (units.km/units.s)
    particles.move_to_center()
    return particles, conv
def plummer(particles, conv):
    gravity = ph4(conv)
    
    gravity.particles.add_particles(particles)
    cluster = gravity.particles[0]
    print cluster.position
    while gravity.model_time < 10000000 | units.yr:
        gravity.evolve_model(gravity.model_time + (1000 | units.yr))
    print cluster.position
"""

mass = 1e7 | units.MSun
radius = 10 | units.parsec
N = 1000
#particles, conv = cluster(mass, radius, N)
#plummer(particles, conv)

def plummer(mass, radius, N):
    conv = nbody_system.nbody_to_si(radius, mass)
    pm = new_plummer_model(N, convert_nbody=conv, do_scale=True)
    gravity = ph4(conv)
    pm.positon = pm.position + ((100, 100, 100) | units.parsec)
    pm.velocity= pm.velocity + ((-1000, -1000, -1000) | (units.km/units.s))

    gravity.particles.add_particles(pm)

    print pm[0].position
    while gravity.model_time < 10000 | units.yr:
        gravity.evolve_model(gravity.model_time + (1000 | units.yr))
    print pm[0].position
    gravity.stop()
    print pm.center_of_mass_velocity()
    plt.scatter(pm.x.value_in(units.parsec), pm.y.value_in(units.parsec))
    plt.show()
plummer(mass, radius, N)

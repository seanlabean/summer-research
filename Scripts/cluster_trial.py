from amuse.lab import *
import numpy as np
from interaction import bbh_orbit
from interaction import impact_parameter
from interaction import bbh_hvgc
from interaction import integrate_system
import math
import matplotlib.pyplot as plt
#Trial 1#
#2336 & 0.25 & 3.0 & 4.5 & 4.61 & 4771 & 7768 \\
#3436 & 0.25 & 1.7 & 2.5 & 1.47 & 28584 & 67621 \\ 
#3000 & 0.25 & 1.7 & 2.5 & 1.68 & 26316 & 36027 \\
bbh_mass = 7e9 | units.MSun
mass_ratio = 0.25
bbh_separation = 1.7 | units.parsec
gc_closest = 2.5 | units.parsec
bbh_phase = 1.68
gc_vinf = 500 | (units.km/units.s)
gc_mass = 1e7 | units.MSun
gc_radius = 6 | units.parsec
def plummer(mass, radius, N, impact_parameter):
    conv = nbody_system.nbody_to_si(radius, mass)
    pm = new_plummer_model(N, convert_nbody=conv, do_scale=True)
    pm.position = pm.position + ((impact_parameter, -100, 0) | units.parsec)
    pm.velocity = (0|units.kms, gc_vinf, 0|units.kms)
    return pm
def particle_set(bbh_mass, mass_ratio, bbh_separation, bbh_phase):
    bh1_mass = bbh_mass/(1+mass_ratio)
    bh2_mass = mass_ratio*bh1_mass
    r1, r2, v1, v2 = bbh_orbit(bbh_mass, mass_ratio, bbh_separation)

    particles = Particles(3)
    bh1 = particles[0]
    bh1.mass = bh1_mass
    bh1.position = (r1*math.cos(bbh_phase), r1*math.sin(bbh_phase), 0|units.km)
    bh1.velocity = (-v1*math.sin(bbh_phase), v1*math.cos(bbh_phase), 0|units.kms)
    
    bh2 = particles[1]
    bh2.mass = bh2_mass
    bh2.position = (-r2*math.cos(bbh_phase), -r2*math.sin(bbh_phase), 0|units.km)
    bh2.velocity = (v2*math.sin(bbh_phase), -v2*math.cos(bbh_phase), 0|units.kms)
    

    particles.move_to_center()
    return particles

b = impact_parameter(bbh_mass, gc_closest, gc_vinf)
pm_particles = plummer(gc_mass, gc_radius, 1000, b)
bbh_particles = particle_set(bbh_mass, mass_ratio, bbh_separation, bbh_phase)
'''
Plot Testing to see if particles are placed correctly.

gravitating_bodies = ParticlesSuperset([pm_particles, bbh_particles])

plt.plot(gravitating_bodies.x.value_in(units.parsec), gravitating_bodies.y.value_in(units.parsec), '.')
plt.plot(pm_particles.x.value_in(units.parsec), pm_particles.y.value_in(units.parsec), '.')
plt.plot(bbh_particles.x.value_in(units.parsec), bbh_particles.y.value_in(units.parsec), '.')
plt.show()
'''
def gravity_integration(pm_particles, bbh_particles):
    gravitating_bodies = ParticlesSuperset([pm_particles, bbh_particles])
    convert_nbody = nbody_system.nbody_to_si(gravitating_bodies.mass.sum(), bbh_particles[0].position.length())
    gravity = ph4(convert_nbody)
    gravity.particles.add_particles(gravitating_bodies)
    bbh1 = gravity.particles[1001]
    bbh2 = gravity.particles[1002]
    pm = gravity.particles[0:1000]
    pm_com = []
    pm_Lag_rad1 = []
    pm_Lag_rad5 = []
    time = []
    pm_old_x = []
    pm_old_y = []
    pm_new_x = []
    pm_new_y = []
    pm_old_x.append(pm.x.value_in(units.parsec))
    pm_old_y.append(pm.y.value_in(units.parsec))
    #plt.subplot()
    #plt.scatter(pm.x.value_in(units.parsec),pm.y.value_in(units.parsec))
    epsilon = (0.1*bbh1.position.length())
    print epsilon
    gravity.parameters.epsilon_squared = epsilon**2
    print gravity.parameters
    while gravity.model_time < (2000000 | units.yr):
        gravity.evolve_model(gravity.model_time + (10000 | units.yr))
        pm_com.append(pm.center_of_mass().length().value_in(units.parsec))
        time.append(gravity.model_time.value_in(units.yr))
        lr, mf = pm.LagrangianRadii(unit_converter=convert_nbody)
        pm_Lag_rad1.append(lr[1].value_in(units.parsec))
        pm_Lag_rad5.append(lr[5].value_in(units.parsec))
    print pm.center_of_mass_velocity().in_(units.kms)
#    plt.subplot()
#    plt.scatter(pm.x.value_in(units.parsec),pm.y.value_in(units.parsec))
#    plt.show()
    pm_new_x.append(pm.x.value_in(units.parsec))
    pm_new_y.append(pm.y.value_in(units.parsec))
    gravity.stop()
    return pm_com, time, pm_old_x, pm_old_y, pm_new_x, pm_new_y, pm_Lag_rad1, pm_Lag_rad5

pm_com, time, old_x, old_y, new_x, new_y, lr1, lr5 = gravity_integration(pm_particles, bbh_particles)

plt.figure()
plt.scatter(old_x,old_y)
plt.scatter(new_x, new_y)
plt.show()
plt.figure()
plt.plot(time,lr1)
plt.plot(time,lr5)
plt.show()

#plt.figure()
#plt.plot(pm_com, time)
#plt.show()

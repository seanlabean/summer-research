from amuse.lab import *
import numpy as np
from interaction import bbh_orbit
from interaction import impact_parameter
from interaction import bbh_hvgc
from interaction import integrate_system
import math
import json


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


def gravity_integration(pm_particles, bbh_particles, gc_radius, runtime):
    gravitating_bodies = ParticlesSuperset([pm_particles, bbh_particles])
    convert_nbody = nbody_system.nbody_to_si(gravitating_bodies.mass.sum(),
                                             bbh_particles[0].position.length())
    gravity = ph4(convert_nbody)
    gravity.particles.add_particles(gravitating_bodies)
    bbh1 = gravity.particles[1001]
    bbh2 = gravity.particles[1002]
    pm = gravity.particles[0:1000]
    pm_com = []
    pm_Lag_rad1 = []
    pm_Lag_rad2 = []
    pm_Lag_rad4 = []
    pm_Lag_rad6 = []
    time = []
    
    epsilon = (0.1*gc_radius)
    print epsilon
    gravity.parameters.epsilon_squared = epsilon**2
    print gravity.parameters
    
    timestep = runtime.value_in(units.yr)/1000
    print timestep
    count = 0
    while gravity.model_time < (runtime):
        gravity.evolve_model(gravity.model_time + (timestep | units.yr))
        pm_com.append(pm.center_of_mass().length().value_in(units.parsec))
        time.append(gravity.model_time.value_in(units.yr))
        lr, mf = pm.LagrangianRadii(unit_converter=convert_nbody)
        pm_Lag_rad1.append(lr[1].value_in(units.parsec))
        pm_Lag_rad2.append(lr[2].value_in(units.parsec))
        pm_Lag_rad4.append(lr[4].value_in(units.parsec))
        pm_Lag_rad6.append(lr[6].value_in(units.parsec))
        count += 1
        print 'step'+str(count)
    
    pm_end = gravity.particles[0:1000].copy()
    gravity.stop()
    return pm_start, pm_end, time, pm_Lag_rad1, pm_Lag_rad2, pm_Lag_rad4, pm_Lag_rad6

bbh_mass = 7e9 | units.MSun
mass_ratio = 0.25
bbh_separation = 1.7 | units.parsec
gc_closest = 2.5 | units.parsec
bbh_phase = 4.61
gc_vinf = 500 | (units.km/units.s)
gc_mass = 1e7 | units.MSun
gc_radius = 0.1 | units.parsec

b = impact_parameter(bbh_mass, gc_closest, gc_vinf)
pm_particles = plummer(gc_mass, gc_radius, 1000, b)
pm_start = pm_particles.copy()
bbh_particles = particle_set(bbh_mass, mass_ratio, bbh_separation, bbh_phase)
'''
runtime = 1000000 | units.yr
pm_start, pm_end, time, lr1, lr2, lr4, lr6 = gravity_integration(pm_particles, bbh_particles, gc_radius, runtime)

file_location = './Text/Cluster_Trials/0-1pc/'
trial_name = '3000kms1Myr'
#write_set_to_file(pm_start, file_location+trial_name+'_start.hdf5',format='hdf5')
write_set_to_file(pm_end, file_location+trial_name+'_end.hdf5',format='hdf5')


with open(file_location+trial_name+'_time.txt',"w+") as time_file:
    json.dump(time, time_file)
with open(file_location+trial_name+'_lr1.txt',"w+") as lr1_file:
    json.dump(lr1, lr1_file)
with open(file_location+trial_name+'_lr2.txt',"w+") as lr2_file:
    json.dump(lr2, lr2_file)
with open(file_location+trial_name+'_lr4.txt',"w+") as lr4_file:
    json.dump(lr4, lr4_file)
with open(file_location+trial_name+'_lr6.txt',"w+") as lr6_file:
    json.dump(lr6, lr6_file)
'''

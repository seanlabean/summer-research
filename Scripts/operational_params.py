from amuse.lab import Particles, units
from amuse.units import constants
import numpy as np
from interaction import bbh_orbit
from interaction import impact_parameter
from integrator import integrate_system_final
from interaction import bbh_hvgc
import math
import time
#import matplotlib.pyplot as plt
#Maybe import *?

#import numpy as np
#bbh_mass_ratio = [0.01, 0.05, 0.1]
bbh_mass_ratio = [0.1]
bbh_separation = [1.7] | units.parsec
#bbh_separation = [1.7, 3, 5] | units.parsec
#gc_closest_ratio = [1.5, 2, 2.5]
gc_closest_ratio = [1.5, 2.5] 
phase_step = np.pi/4
#bbh_phase = np.arange(0, 2*np.pi, phase_step)
#bbh_phase = bbh_phase.tolist()
#particles = bbh_hvgc(bbh_mass, mass_ratio, bbh_separation, bbh_phase, gc_closest_ratio, gc_inf)

def operational_params(bbh_mass_ratio, bbh_separation,  gc_closest_ratio, bbh_phase_step):
    bbh_mass = 7e9 | units.MSun
    gc_inf = 500 | (units.km/units.s)
    end_position = 120 | units.parsec
    vel_list = np.zeros(40)
    phase_list = np.zeros_like(vel_list)
    bbh_phase = 0.
    index = 0
    gc_closest = bbh_separation*gc_closest_ratio
    while bbh_phase < 2*np.pi:
        particles = bbh_hvgc(bbh_mass, bbh_mass_ratio, bbh_separation, bbh_phase, gc_closest, gc_inf)
        final_vel, energy = integrate_system_final(particles, end_position)
        if final_vel.value_in(units.kms) > 1000:
            vel_list[index] = final_vel.value_in(units.kms)
        print bbh_phase
        phase_list[index] = bbh_phase
        index += 1
        bbh_phase += phase_step 
    #max_vel= max(vel_list)
    #max_vel_index = vel_list.index(max_vel)
    #phase_max = phase_list[max_vel_index]

    return vel_list, phase_list

start_time = time.time()
f = open("results1.txt", "w+")
count = 1
total = len(gc_closest_ratio)*len(bbh_mass_ratio)*len(bbh_separation)
for a in gc_closest_ratio:
    for i in bbh_mass_ratio:
        for j in bbh_separation:
            print "%d out of %d Calculating maximum velocity..." %(count, total)
            vel_list, phase_list = operational_params(i, j, a, phase_step)
            for v in enumerate(vel_list):
                vel = v[1]
                vel_index = v[0]
                phase = phase_list[vel_index]
                if vel != 0:
                    result_str = "{} {} {} {} {} {} {}\n".format(vel, i, j.value_in(units.parsec), a, phase, 0., 0.)
                    f.write(result_str)
                    print "Wrote"
                else:
                    continue
            count += 1
f.close()
print "Done."
print "Data written to results.txt"
print("--- %s seconds ---" % (time.time() - start_time))

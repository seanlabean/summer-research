from amuse.lab import *
import numpy as np
import matplotlib.pyplot as plt
def integrate_cluster_long(cluster_set, integration_time):
    convert_nbody = nbody_system.nbody_to_si(cluster_set.mass.sum(),
                                             cluster_set.center_of_mass().length()
                                             )
    gravity = ph4(convert_nbody)
    epsilon = 0.01 | units.parsec
    gravity.parameters.epsilon_squared = epsilon**2
    print gravity.parameters
    gravity.particles.add_particles(cluster_set)
    count = 10
    while gravity.model_time < integration_time:
        gravity.evolve_model(gravity.model_time + (10000000 | units.yr))
        print str(count)+"million"
        count += 10
    cluster_end = gravity.particles[:].copy()
    gravity.stop()
    return cluster_end

pc1 = './Text/Cluster_Trials/1pc/3000kms1Myr_end.hdf5'
cluster_begin = read_set_from_file(pc1, format='hdf5')
cluster_end = integrate_cluster_long(cluster_begin, 2e9 | units.yr)
write_set_to_file(cluster_end, './Text/Cluster_Trials/1pc/3000kms2Gyr_end.hdf5', format='hdf5')

print 'File saved to ' + pc1
'''
neighbour_dist = []

for particle in cluster_end:
    nearest = np.sqrt(particle.nearest_neighbour(cluster_end).x.value_in(units.parsec)**2 +
                      particle.nearest_neighbour(cluster_end).y.value_in(units.parsec)**2 +
                      particle.nearest_neighbour(cluster_end).z.value_in(units.parsec)**2)
    neighbour_dist.append(nearest)

print min(neighbour_dist)
'''
plt.figure(figsize=(7,7))
plt.scatter(cluster_end.x.value_in(units.parsec), 
            cluster_end.y.value_in(units.parsec),
            s=1
            )
plt.show()

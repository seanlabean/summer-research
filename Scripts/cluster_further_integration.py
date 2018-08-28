from amuse.lab import *
import numpy as np
import matplotlib.pyplot as plt
def integrate_cluster_long(cluster_set, integration_time):
    convert_nbody = nbody_system.nbody_to_si(cluster_set.mass.sum(),
                                             cluster_set.center_of_mass().length()
                                             )
    gravity = ph4(convert_nbody)
    print gravity.parameters
    gravity.particles.add_particles(cluster_set)
    while gravity.model_time < integration_time:
        gravity.evolve_model(gravity.model_time + (1000000 | units.yr))
    cluster_end = gravity.particles[:].copy()
#    print len(cluster_end.x.value_in(units.parsec))
    gravity.stop()
    return cluster_end

pc6 = './Text/Cluster_Trials/6pc/3000kms200kyr_end.hdf5'
cluster_begin = read_set_from_file(pc6, format='hdf5')
cluster_end = integrate_cluster_long(cluster_begin, 1e6 | units.yr)
write_set_to_file(cluster_end, './Text/Cluster_Trials/6pc/3000kms1Gyr_end.hdf5', format='hdf5')

print 'File saved to ' + pc6
neighbour_dist = []

for particle in cluster_end:
    nearest = np.sqrt(particle.nearest_neighbour(cluster_end).x.value_in(units.parsec)**2 +
                      particle.nearest_neighbour(cluster_end).y.value_in(units.parsec)**2 +
                      particle.nearest_neighbour(cluster_end).z.value_in(units.parsec)**2)
    neighbour_dist.append(nearest)

print min(neighbour_dist)
plt.figure(figsize=(7,7))
plt.scatter(cluster_end.x.value_in(units.parsec), 
            cluster_end.y.value_in(units.parsec),
            s=1
            )
plt.show()

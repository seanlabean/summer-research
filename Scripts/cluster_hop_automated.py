from amuse.lab import *
import matplotlib.pyplot as pyplot
import numpy as np

def execute_hop(cluster):
    cluster_total_E = cluster.kinetic_energy() + cluster.potential_energy()
    cluster_total_M = cluster.mass.sum()
    conv = nbody_system.nbody_to_si(cluster_total_M, cluster_total_E)
    hop = Hop(unit_converter=conv)
    
    hop.parameters.number_of_neighbors_for_hop = 4
    hop.parameters.number_of_particles_per_group_pair_boundary = 2
    hop.relative_saddle_density_threshold = True
    hop.set_relative_saddle_density_threshold(True)
    #hop.set_nHop(16)                                              
    hop.parameters.number_of_buckets = 16
    print hop.parameters
    hop.particles.add_particles(cluster)
    hop.calculate_densities()
    hop.do_hop()
    groups = [x.get_intersecting_subset_in(cluster) for x in hop.groups()]
    n_groups = hop.get_number_of_groups()
    print n_groups
    for g in groups:
        print g
    return groups

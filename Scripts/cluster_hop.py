from amuse.lab import *
import matplotlib.pyplot as pyplot
import numpy as np
cluster = read_set_from_file('./Text/Cluster/cluster_set.hdf5',format='hdf5')

cluster_total_E = cluster.kinetic_energy() + cluster.potential_energy()
cluster_total_M = cluster.mass.sum()
conv = nbody_system.nbody_to_si(cluster_total_M, cluster_total_E)
print cluster.LagrangianRadii(conv)
hop = Hop(unit_converter=conv)

hop.parameters.number_of_neighbors_for_hop = 10
hop.parameters.number_of_particles_per_group_pair_boundary = 2
hop.relative_saddle_density_threshold = True
hop.set_relative_saddle_density_threshold(True)
#hop.set_nHop(16)
hop.parameters.number_of_buckets = 16
print hop.parameters
hop.particles.add_particles(cluster)
hop.calculate_densities()
hop.do_hop()

n_groups = hop.get_number_of_groups()
print n_groups

def plot_clumps(groups, total_mass):
    number_of_particles_in_group = []
    fraction_of_mass_in_group =  []

    for group in groups:
        number_of_particles_in_group.append(len(group))
        fraction = (group.mass.sum()/total_mass)
        fraction_of_mass_in_group.append(fraction)

    figure = pyplot.figure(figsize= (12,6))


    subplot = figure.add_subplot(1, 2, 1)

    colormap = pyplot.cm.Dark2 #pyplot.cm.Paired
    for index, group in enumerate(groups):
        color = colormap(index) #colormap(1.0 * index / len(groups))
        subplot.scatter(
            group.x.value_in(units.parsec),
            group.y.value_in(units.parsec),
            s = 1,
            edgecolors = color,
            facecolors = color,
            label = "{:.2F}".format(fraction_of_mass_in_group[index])
        )

    #subplot.set_xlim(0,1)
    #subplot.set_ylim(0,1)
    subplot.set_xlabel('x (parsec)')
    subplot.set_ylabel('y (parsec)')
    subplot.legend()

    subplot = figure.add_subplot(1, 2, 2)

    subplot.scatter(
        number_of_particles_in_group,
        fraction_of_mass_in_group,
    )

    subplot.set_xscale('log')
    subplot.set_yscale('log')
    subplot.set_xlabel('N')
    subplot.set_ylabel('df/d(Log_10 N)')

    #figure.savefig('x.png')
    pyplot.show()
def vector_norm(vector):
    norm = np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
    return norm

groups = [x.get_intersecting_subset_in(cluster) for x in hop.groups()]
#plot_clumps(groups,cluster_total_M)
for group in groups:
    print group.mass.sum().in_(units.MSun)
    print vector_norm(group.center_of_mass_velocity().in_(units.kms))
    print group.LagrangianRadii(conv)[0][5].in_(units.parsec), group.LagrangianRadii(conv)[1][5]
hop.stop()

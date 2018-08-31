from amuse.lab import *
import matplotlib.pyplot as pyplot
import numpy as np
import matplotlib
from cluster_hop_automated import execute_hop
cluster = read_set_from_file('./Text/Cluster_Trials/6pc/3000kms2Gyr_end.hdf5',format='hdf5')

cluster_total_E = cluster.kinetic_energy() + cluster.potential_energy()
cluster_total_M = cluster.mass.sum()
conv = nbody_system.nbody_to_si(cluster_total_M, cluster_total_E)
print cluster.LagrangianRadii(conv)
hop = Hop(unit_converter=conv)

hop.parameters.number_of_neighbors_for_hop = 10
hop.parameters.number_of_particles_per_group_pair_boundary = 4
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

    figure = pyplot.figure(figsize= (6,6))
    params = {'axes.labelsize': 16,
          'axes.titlesize':18, 
          'font.size': 14, 
          'legend.fontsize': 14, 
          'xtick.labelsize': 14, 
          'ytick.labelsize': 14
         }
    matplotlib.rcParams.update(params)

    #subplot = figure.add_subplot(1, 2, 1)

    colormap = pyplot.cm.Dark2 #pyplot.cm.Paired
    for index, group in enumerate(groups):
        color = colormap(index) #colormap(1.0 * index / len(groups))
        pos,coreradius,coredens=group.densitycentre_coreradius_coredens(conv)
        pyplot.scatter(
            group.x.value_in(units.parsec)/1000000,
            group.y.value_in(units.parsec)/1000000,
            s = 1,
            edgecolors = color,
            facecolors = color,
            label = "mf = {:.3F}".format(fraction_of_mass_in_group[index])
        )
        density_string = ''#r"$\frac{M_\odot}{kpc^3}$"
        pyplot.scatter(pos.x.value_in(units.parsec)/1000000, 
                       pos.y.value_in(units.parsec)/1000000,
                       s=15,
                       color = 'black',
                       marker = 'x',
                       label = '{:.3F}'.format(coredens.value_in(units.MSun/(1000*units.parsec)**3))+density_string
                       )
        '''
        subplot.scatter(
            group.center_of_mass().x.value_in(units.parsec), 
            group.center_of_mass().y.value_in(units.parsec),
            marker='o',
            color='red',
            s=30,
            #edgecolors = color,
            #facecolors = color,
            label = "CoM"
        )
        '''
    #subplot.set_xlim(0,1)
    #subplot.set_ylim(0,1)
    #subplot.set_xlabel('x (parsec)')
    #subplot.set_ylabel('y (parsec)')
    pyplot.xlim(-10,10)
    pyplot.ylim(-10,10)
    pyplot.yticks(np.arange(-10,11,step=5))
    pyplot.title('Subgroups in 1pc Cluster After 2 Gyr')
    pyplot.minorticks_on()
    pyplot.tick_params(which='minor',direction='in')
    pyplot.tick_params(size = 6, width = 1, direction='in')
    pyplot.grid(alpha=0.25)
    pyplot.xlabel('x (Mpc)')
    pyplot.ylabel('y (Mpc)')
    pyplot.legend(markerscale=2)
    """
    subplot = figure.add_subplot(1, 2, 2)

    subplot.scatter(
        number_of_particles_in_group,
        fraction_of_mass_in_group,
    )

    subplot.set_xscale('log')
    subplot.set_yscale('log')
    subplot.set_xlabel('N')
    subplot.set_ylabel('df/d(Log_10 N)')
    """
    #figure.savefig('./Images/groups_1pc2Gyr')
    pyplot.show()
def vector_norm(vector):
    norm = np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
    return norm

groups = [x.get_intersecting_subset_in(cluster) for x in hop.groups()]
plot_clumps(groups,cluster_total_M)
for group in groups:
    pos,coreradius,coredens=group.densitycentre_coreradius_coredens(conv)
    print "core radius: ", coreradius.value_in(units.parsec), "core density", coredens.value_in(units.MSun/units.parsec**3)
    print "group mass: ", group.mass.sum().in_(units.MSun)
    print "group c_o_m velocity: ", vector_norm(group.center_of_mass_velocity().in_(units.kms))
    print "group LR 0.02 mf: ", group.LagrangianRadii(conv)[0][1].in_(units.parsec), group.LagrangianRadii(conv)[1][1]
hop.stop()
def evolve_group(group, end_time):
    group_conv = nbody_system.nbody_to_si(group.mass.sum(), group.center_of_mass().length())
    gravity = ph4(group_conv)
    gravity.particles.add_particles(group)

    lr1 = []
    lr2 = []
    time = []
    while gravity.model_time < end_time:
        gravity.evolve_model(gravity.model_time + (1000 | units.yr))
        lr, mf = gravity.particles.LagrangianRadii(group_conv)
        time.append(gravity.model_time.value_in(units.yr))
        lr1.append(lr[0].value_in(units.parsec))
        lr2.append(lr[1].value_in(units.parsec))
    pyplot.scatter(gravity.particles.x.value_in(units.parsec), gravity.particles.y.value_in(units.parsec),s=1)
    pyplot.show()
    pyplot.figure()
    pyplot.plot(time,lr1)
    gravity.stop()
#evolve_group(groups[4], 1000000 | units.yr)
#groups_new = execute_hop(evolved_group)
plot_clumps(cluster, cluster_total_M)

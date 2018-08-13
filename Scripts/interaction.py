from amuse.lab import Particles, units
from amuse.units import constants
import math


def bbh_orbit(bbh_mass, mass_ratio, separation):
    """Input: total BBH mass in Solar Masses, black hole mass ratio, separation in parsecs.
    Will generate orbital radius and velocity of each black hole that satisfies a e=0 planar orbit.
    Output: BH orbit radius 1 and 2 in parsecs, BH orbit velocities 1 and 2 in km/s.
    """
    G = constants.G
    a = separation
    m1 = bbh_mass/(1+mass_ratio)
    m2 = mass_ratio*m1
    r1 = (m2*a)/(bbh_mass)
    r2 = (m1*a)/(bbh_mass)
    T = ((4*(math.pi**2)*a**3)/(G*bbh_mass)).sqrt()
    v1 = (2*math.pi*r1)/T
    v2 = (2*math.pi*r2)/T

    return r1.in_(units.parsec), r2.in_(units.parsec), v1.in_(units.km/units.s), v2.in_(units.km/units.s)

def impact_parameter(bbh_mass, gc_closest, gc_vinf):
    """Input: total BBH mass in Solar Masses, GC closest approach to BBH CM, GC initial velocity at infinity.
    Calculates impact parameter necessary for given closest approach.
    Output: Impact parameter VALUE in units of parsecs.
    """
    G = constants.G.value_in((units.parsec*units.km**2)/(units.MSun*units.s**2))
    M = bbh_mass.value_in(units.MSun)
    r_min = gc_closest.value_in(units.parsec)
    v_inf = gc_vinf.value_in(units.km/units.s)
    b = r_min*math.sqrt(1 + (2*G*M)/(r_min*v_inf**2)) | units.parsec

    return b.value_in(units.parsec)

def bbh_hvgc(bbh_mass, mass_ratio, bbh_separation, bbh_phase, gc_closest, gc_vinf):
    """Creates particle list, defines circular-planar BBH orbit, places GC. """
    bh1_mass = bbh_mass/(1+mass_ratio)
    bh2_mass = mass_ratio*bh1_mass
    r1, r2, v1, v2 = bbh_orbit(bbh_mass, mass_ratio, bbh_separation)

    particles = Particles(3)
    bh1 = particles[0]
    bh1.mass = bh1_mass
    bh1.position = (r1*math.cos(bbh_phase), r1*math.sin(bbh_phase), 0|units.km) 
    bh1.velocity = (-v1*math.sin(bbh_phase), v1*math.cos(bbh_phase), 0|units.kms) 
    #bh1.position = (0,r1,0) 
    #bh1.velocity = (-v1,0,0) 

    bh2 = particles[1]
    bh2.mass = bh2_mass
    bh2.position = (-r2*math.cos(bbh_phase), -r2*math.sin(bbh_phase), 0|units.km)
    bh2.velocity = (v2*math.sin(bbh_phase), -v2*math.cos(bbh_phase), 0|units.kms)
#    bh2.position = (0,-r2,0) 
#    bh2.velocity = (v1,0,0) 

    b = impact_parameter(bbh_mass, gc_closest, gc_vinf)
    hvgc = particles[2]
    hvgc.mass = 1.0e7 | units.MSun
    hvgc.position = (b, -100, 0) | units.parsec
    hvgc.velocity = (0|units.kms, gc_vinf, 0|units.kms)  
    
    particles.move_to_center()

    return particles


def integrate_system(particles, end_position):
    """Given a set of particles, will integrate system until particle[2] reaches end_position. Returns the velocity of particle[2] at end_position."""
    from amuse.lab import ph4, nbody_system
    convert_nbody = nbody_system.nbody_to_si(particles.mass.sum(), particles[0].position.length())
    gravity = ph4(convert_nbody)
    gravity.particles.add_particles(particles)
    
    hvgc = gravity.particles[2]
    hvgc_pos = (hvgc.x**2 + hvgc.y**2 + hvgc.z**2).sqrt()
    captured = False
    while hvgc_pos < end_position:
        gravity.evolve_model(gravity.model_time + (100 | units.yr))
        hvgc_pos = (hvgc.x**2 + hvgc.y**2 + hvgc.z**2).sqrt()
        if gravity.model_time > (200000 | units.yr):
            hvgc_pos = end_position
            captured = True
            break
    if captured == True:
        hvgc_vel = 0 | units.kms
    else:
        hvgc_vel = (hvgc.vx**2 + hvgc.vy**2 + hvgc.vz**2).sqrt()
    gravity.stop()
    return hvgc_vel

def bbh_phase_loop(bbh_mass, mass_ratio, separation, bh_phase, gc_closest, gc_vinf, end_position):
    """Makes particles given initial properties, integrates until GC reaches end position. Repeats for different initial black hole phase. Saves final velocity of GC and corresponding BH phase. Plots phase vs. velocity to show phase dependence and maximum velocity kick."""
    from matplotlib import pyplot
    phase_list = []
    velocity_list = []
    phase = 0
    while phase <= 2*math.pi:
        particles = bbh_hvgc(bbh_mass, mass_ratio, separation, phase, gc_closest, gc_vinf)
        hvgc_vel = integrate_system(particles, end_position).value_in(units.kms)
        phase_list.append(phase)
        velocity_list.append(hvgc_vel)
        phase += math.pi/20
        print phase

    figure = pyplot.figure(figsize=(9,9))
    pyplot.rcParams.update({'font.size': 24})
    plot = figure.add_subplot(1,1,1)
    ax = pyplot.gca()
    ax.minorticks_on()
    ax.locator_params(nbins=3)
    x_label = 'phase [rad]'
    y_label = 'velocity at 120pc [km/s]'
    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)
    plot.plot(phase_list, velocity_list, ls = '-')
    save_file = 'phase_loop.png'
    pyplot.savefig(save_file)
    pyplot.show()


########################################################################################
#Plotting#
########################################################################################

def integrate_plotting(particles, end_position):
    """Input: Set of particles and end condition. Integrates until GC reaches end condition. Records positions of each particle throughout integration in steps of 100 years."""
    from amuse.lab import ph4, nbody_system
    convert_nbody = nbody_system.nbody_to_si(particles.mass.sum(), particles[0].position.length())
    gravity = ph4(convert_nbody)
    gravity.particles.add_particles(particles)
    
    bh1 = gravity.particles[0]
    bh2 = gravity.particles[1]
    hvgc = gravity.particles[2]
    
    x_hvgc = [] | units.parsec
    y_hvgc = [] | units.parsec
    x_bh1 = [] | units.parsec
    y_bh1 = [] | units.parsec
    x_bh2 = [] | units.parsec
    y_bh2 = [] | units.parsec
    gcbh1_dist = [] | units.parsec
    gcbh2_dist = [] | units.parsec
    hvgc_pos = (hvgc.x**2 + hvgc.y**2 + hvgc.z**2).sqrt()
    while hvgc_pos < end_position:
        gravity.evolve_model(gravity.model_time + (100 | units.yr))
        x_hvgc.append(hvgc.x)
        y_hvgc.append(hvgc.y)
        x_bh1.append(bh1.x)
        y_bh1.append(bh1.y)
        x_bh2.append(bh2.x)
        y_bh2.append(bh2.y)
        gcbh1_dist.append(hvgc.position.length()-bh1.position.length())
        hvgc_pos = (hvgc.x**2 + hvgc.y**2 + hvgc.z**2).sqrt()
    gcbh1_distmin = min(gcbh1_dist)
    hvgc_vel = (hvgc.vx**2 + hvgc.vy**2 + hvgc.vz**2).sqrt()
    gravity.stop()
    return hvgc_vel, x_hvgc, y_hvgc, x_bh1, y_bh1, x_bh2, y_bh2

def plot_track(xgc, ygc, xbh1, ybh1, xbh2, ybh2, output_filename):
    """Input: particle position lists. Output: visual of particle paths."""
    from matplotlib import pyplot
    figure = pyplot.figure(figsize=(7,7))
    pyplot.rcParams.update({'font.size': 16})
    plot = figure.add_subplot(1,1,1)
    ax = pyplot.gca()
    ax.minorticks_on()
    ax.locator_params(nbins=3)

    x_label = 'x [pc]'
    y_label = 'y [pc]'
    title_ = 'GC Capture by SMBBH'
    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)
    pyplot.title(title_)
#    pyplot.xlim(-5,5)                                                                                       
#    pyplot.ylim(-5,5)                                                                                       
    plot.plot(xgc.value_in(units.parsec), ygc.value_in(units.parsec), color = 'r', label='GC')
    plot.plot(xbh1.value_in(units.parsec), ybh1.value_in(units.parsec), color = 'k',label='BH1')
    plot.plot(xbh2.value_in(units.parsec), ybh2.value_in(units.parsec), color = 'b',label='BH2')
    plot.legend(loc='lower center')
    save_file = './Images/'+output_filename
    pyplot.savefig(save_file)
    print '\nSaved figure in file', save_file, '\n'
    pyplot.show()


########################################################################################                          
#MINIMUM DISTANCE PLOT#                                                                              
########################################################################################

def integrate_minimum_distance(particles, end_position):
    from amuse.lab import ph4, nbody_system
    convert_nbody = nbody_system.nbody_to_si(particles.mass.sum(), particles[0].position.length())
    gravity = ph4(convert_nbody)
    gravity.particles.add_particles(particles)

    bh1 = gravity.particles[0]
    bh2 = gravity.particles[1]
    hvgc = gravity.particles[2]
    hvgc_vel_list = [] | units.kms
    gcbh1_dist = [] | units.parsec
    gcbh2_dist = [] | units.parsec
    hvgc_pos = (hvgc.x**2 + hvgc.y**2 + hvgc.z**2).sqrt()
    while hvgc_pos < end_position:
        gravity.evolve_model(gravity.model_time + (100 | units.yr))
        gcbh1_dist.append((hvgc.position-bh1.position).length())
        gcbh2_dist.append((hvgc.position-bh2.position).length())
        hvgc_vel_list.append(hvgc.velocity.length())
        hvgc_pos = (hvgc.x**2 + hvgc.y**2 + hvgc.z**2).sqrt()
    gcbh1_distmin = min(gcbh1_dist)
    gcbh2_distmin = min(gcbh2_dist)
    hvgc_max_vel = max(hvgc_vel_list)
    print hvgc_max_vel
    print gcbh1_distmin
    print gcbh2_distmin
    hvgc_vel = (hvgc.vx**2 + hvgc.vy**2 + hvgc.vz**2).sqrt()
    gravity.stop()
    return hvgc_vel, gcbh1_distmin, gcbh2_distmin

def plot_min(bbh_mass, mass_ratio, separation, bh_phase, gc_closest, gc_vinf, end_position):
    from matplotlib import pyplot
    phase_list = []
    velocity_list = []
    bh1min_dist = []
    bh2min_dist = []
    phase = 0
    while phase <= 2*math.pi:
        particles = bbh_hvgc(bbh_mass, mass_ratio, separation, phase, gc_closest, gc_vinf)
        hvgc_vel, gcbh1_distmin, gcbh2_distmin = integrate_minimum_distance(particles, end_position)
        phase_list.append(phase)
        velocity_list.append(hvgc_vel.value_in(units.kms))
        bh1min_dist.append(gcbh1_distmin.value_in(units.parsec))
        bh2min_dist.append(gcbh2_distmin.value_in(units.parsec))
        phase += math.pi/4
        print phase

    figure = pyplot.figure(figsize=(9,9))
    pyplot.rcParams.update({'font.size': 30})
    plot = figure.add_subplot(1,1,1)
    ax = pyplot.gca()
    ax.minorticks_on()
    ax.locator_params(nbins=3)

    x_label = 'phase [rad]'
    y_label = 'minimum separation [parsec]'
    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)
    plot.plot(phase_list, bh1min_dist, ls = '-')
    plot.plot(phase_list, bh2min_dist, ls = '-')
    save_file = './Images/'+'bh_min_dist.png'
    pyplot.savefig(save_file)
    pyplot.show()
#plot_min(bbh_mass, mass_ratio, separation, bh_phase, gc_closest, gc_vinf, end_position)


########################################################################################    
#ENERGY PLOT#                                                    
########################################################################################
def integrate_energy(particles, end_condition):
    """Given particle set and end condition, integrate system and record energy of each particle. Return individual energies, total energy, time."""
    from amuse.lab import ph4, nbody_system
    convert_nbody = nbody_system.nbody_to_si(particles.mass.sum(), particles[0].position.length())
    gravity = ph4(convert_nbody)
    gravity.particles.add_particles(particles)

    bh1 = gravity.particles[0]
    bh2 = gravity.particles[1]
    hvgc = gravity.particles[2]
    bh1_energy = [] 
    bh2_energy = [] 
    bbh_energy = [] 
    hvgc_energy = [] 
    total_energy = [] 
    time = []
    hvgc_pos = (hvgc.x**2 + hvgc.y**2 + hvgc.z**2).sqrt()
    while hvgc_pos < end_condition:
        gravity.evolve_model(gravity.model_time + (1000 | units.yr))
        com_pos = particles.center_of_mass()
        com_vel = particles.center_of_mass_velocity()
        rel_hvgc_pos = hvgc.position - com_pos
        rel_hvgc_vel = hvgc.velocity - com_vel
        #bh1_total_energy = bh1.specific_kinetic_energy()*bh1.mass.in_(units.kg) - bh1.potential()*bh1.mass.in_(units.kg)
        #bh2_total_energy = bh2.specific_kinetic_energy()*bh2.mass.in_(units.kg) - bh2.potential()*bh2.mass.in_(units.kg)
        #hvgc_total_energy = hvgc.specific_kinetic_energy()*hvgc.mass.in_(units.kg) - hvgc.potential()*hvgc.mass.in_(units.kg)
        hvgc_total_energy = 0.5*hvgc.mass*(rel_hvgc_vel.x**2 + rel_hvgc_vel.y**2 + rel_hvgc_vel.z**2) - (constants.G*hvgc.mass*(bh1.mass+bh2.mass))/(rel_hvgc_pos.x**2 + rel_hvgc_pos.y**2 + rel_hvgc_pos.z**2).sqrt()

        time.append(gravity.model_time.value_in(units.yr))
        #bh1_energy.append(bh1_total_energy.value_in(units.J))
        #bh2_energy.append(bh2_total_energy.value_in(units.J))
        #bbh_energy.append(bh1_total_energy.value_in(units.J) + bh2_total_energy.value_in(units.J))
        hvgc_energy.append(hvgc_total_energy.value_in(units.J))
        #total_energy.append(-1*(particles.kinetic_energy() + particles.potential_energy()).value_in(units.J))
        
        hvgc_pos = (hvgc.x**2 + hvgc.y**2 + hvgc.z**2).sqrt()
#        if gravity.model_time > (600000 | units.yr):
#            hvgc_pos = end_condition
        
    hvgc_vel = (hvgc.vx**2 + hvgc.vy**2 + hvgc.vz**2).sqrt()
    gravity.stop()


    print hvgc_vel
    print "Initial cluster energy: {:.2E}".format(hvgc_energy[0])
    print "Final cluster energy: {:.2E}".format(hvgc_energy[-1])
    print (hvgc_energy[-1] - hvgc_energy[0])/hvgc_energy[-1]
    return bh1_energy, bh2_energy, bbh_energy, hvgc_energy, total_energy, time

def plot_energy(bh1_energy, bh2_energy, bbh_energy, hvgc_energy, total_energy, time, filename):
    """Given energy and time lists, plot time vs. energy."""
    from matplotlib import pyplot
    figure = pyplot.figure(figsize=(6,6))
    pyplot.rcParams.update({'font.size': 16})
    plot = figure.add_subplot(1,1,1)
    ax = pyplot.gca()
    ax.minorticks_on()
    ax.locator_params(nbins=3)
    pyplot.title('Relative Energy of GC')
    x_label = 'time [years]'
    y_label = 'energy [J]'
    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)
#    pyplot.xlim(-5,5)                                                                                                  
#    pyplot.ylim(-5,5)                                                                                                  
    plot.plot(time, hvgc_energy, color = 'r', label='GC')
#    plot.semilogy(time, bbh_energy, color = 'k', label='BBH')
#    plot.semilogy(bh2_energy, label='BH2')
#    plot.semilogy(total_energy, color = 'b',label='Total')
#    plot.legend()
    save_file = './Images/' + filename
    pyplot.savefig(save_file)
    print '\nSaved figure in file', save_file, '\n'
    pyplot.show()

########################################################################################   
#TIDAL ACCELERATION#                                           
########################################################################################

def tidal_ratio(gc_closest_ratio, bbh_mass_ratio, bbh_separation, phase):
    bbh_mass = 7e9 | units.MSun
    gc_inf = 500 | (units.km/units.s)
    gc_closest = bbh_separation*gc_closest_ratio
    end_position = 120 | units.parsec
    particles = bbh_hvgc(bbh_mass, bbh_mass_ratio, bbh_separation, phase, gc_closest, gc_inf)
    from amuse.lab import ph4, nbody_system
    convert_nbody = nbody_system.nbody_to_si(particles.mass.sum(), particles[0].position.length())
    gravity = ph4(convert_nbody)
    gravity.particles.add_particles(particles)

    bh1 = gravity.particles[0]
    bh2 = gravity.particles[1]
    hvgc = gravity.particles[2]
    hvgc_pos = (hvgc.x**2 + hvgc.y**2 + hvgc.z**2).sqrt()
    hvgc_int_acc = 1e7/10
    tide_ratio1 = []
    tide_ratio2 = []
    
    while hvgc_pos < end_position:
        gravity.evolve_model(gravity.model_time + (200 | units.yr))
        r_bh1 = (hvgc.position-bh1.position).length().value_in(units.parsec)
        r_bh2 = (hvgc.position-bh2.position).length().value_in(units.parsec)
        bh1_tide = bh1.mass.value_in(units.MSun)/r_bh1**3
        bh2_tide = bh2.mass.value_in(units.MSun)/r_bh2**3
        bh1_tide_ratio = bh1_tide/hvgc_int_acc
        bh2_tide_ratio = bh2_tide/hvgc_int_acc
        tide_ratio1.append(bh1_tide_ratio)
        tide_ratio2.append(bh2_tide_ratio)
        hvgc_pos = (hvgc.x**2 + hvgc.y**2 + hvgc.z**2).sqrt()
    max_tide1 = max(tide_ratio1)
    max_tide2 = max(tide_ratio2)
    return max_tide1, max_tide2

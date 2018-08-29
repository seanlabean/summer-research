########################################
###
### A module to find binaries in an
### AMUSE particle set. And also
### calculate some binary things.
###
### Written by J. Wall
### Drexel University
### Summer 2018
###
########################################

from amuse.lab import Particles
from scipy.spatial.distance import cdist, pdist
from amuse.lab import units as u
import numpy as np

### Feel free to ignore the first couple of functions
### if you don't use FLASH or yt. - JW

def get_ps_from_hdf5(partfile, particle_type='all'):
    
    """ Returns an AMUSE particle set from a FLASH
        hdf5 particle file (NOT a plot file) or checkpoint file.
        
        Keyword arguments:
        partfile      -- the FLASH particle file
        particle_type -- the type of particle that you
                         want from FLASH. Can be:
                         star -- active FLASH particle type
                         sink -- sink FLASH particle type
                         any other string or none returns all particle types.
                         The default is all.
                         
        Returns:
        stars -- An AMUSE particle set.
    """
    import h5py

    ds = h5py.File(partfile, 'r')
    
    part_data    = ds.get('tracer particles')
    part_names   = ds.get('particle names')
    part_data    = np.array(part_data)
    part_names   = np.array(part_names)
    real_scalars = ds.get('real scalars')
    part_time    = real_scalars[0][1]
    
    if (part_data.size > 1):
        
        mass_ind  = np.where(np.char.strip(part_names)=='mass')[0]
        tag_ind   = np.where(np.char.strip(part_names)=='tag')[0]
        type_ind  = np.where(np.char.strip(part_names)=='type')[0]
        ct_ind    = np.where(np.char.strip(part_names)=='creation_time')[0]
        posx_ind  = np.where(np.char.strip(part_names)=='posx')[0]
        posy_ind  = np.where(np.char.strip(part_names)=='posy')[0]
        posz_ind  = np.where(np.char.strip(part_names)=='posz')[0]
        velx_ind  = np.where(np.char.strip(part_names)=='velx')[0]
        vely_ind  = np.where(np.char.strip(part_names)=='vely')[0]
        velz_ind  = np.where(np.char.strip(part_names)=='velz')[0]
        
        # Get the particle types
        type_part = part_data[:,type_ind].flatten()
        
        if (particle_type == 'star' or particle_type == 'stars'):
            # Find only star particles (particle_type=1)
            star_ind = np.where(type_part==1)[0]
        elif (particle_type == 'sink' or particle_type == 'sinks'):
            # Find only sink particles (particle_type=2)
            star_ind = np.where(type_part==2)[0]
        else:
            star_ind = np.arange(len(type_part))
        
        num_parts = len(star_ind)
        # allocate the array
        stars = Particles(num_parts)
        
        # Masses of stars
        stars.mass = part_data[star_ind,mass_ind] | u.g
        # Tags of stars
        stars.tag = part_data[star_ind,tag_ind]
        # Creation times of stars
        stars.ct   = part_data[star_ind,ct_ind]   | u.s
        # Positions
        stars.x = part_data[star_ind,posx_ind]    | u.cm
        stars.y = part_data[star_ind,posy_ind]    | u.cm
        stars.z = part_data[star_ind,posz_ind]    | u.cm
        # Velocities
        stars.vx = part_data[star_ind,velx_ind]   | u.cm/u.s
        stars.vy = part_data[star_ind,vely_ind]   | u.cm/u.s
        stars.vz = part_data[star_ind,velz_ind]   | u.cm/u.s
        # Set an index attribute. Useful for some functions/codes in AMUSE.
        stars.index_of_the_particle = np.arange(num_parts)

    
    else:
        print "Error: No particles found in this file!"
        stars = None
    
    return stars


def get_nparrays_from_hdf5(partfile, particle_type="all"):
    
    """ Returns an set of numpy arrays from a FLASH
        hdf5 particle file (NOT a plot file) or checkpoint file.
        
        Keyword arguments:
        partfile      -- the FLASH particle file
        particle_type -- the type of particle that you
                         want from FLASH. Can be:
                         star -- active FLASH particle type
                         sink -- sink FLASH particle type
                         any other string or none returns 
                         all particle types.
                         The default is all.
                         
       Returns:
       mass             -- particle mass
       tags             -- particle ID tags
       ct               -- simulation time when particle formed
       pp               -- x,y,z positions 
       pv               -- x,y,z velocities
    """
    import h5py
    ds = h5py.File(partfile, 'r')
    
    part_data    = ds.get('tracer particles')
    part_names   = ds.get('particle names')
    part_data    = np.array(part_data)
    part_names   = np.array(part_names)
    real_scalars = ds.get('real scalars')
    part_time    = real_scalars[0][1]
    
    if (part_data.size > 1):
        
        mass_ind  = np.where(np.char.strip(part_names)=='mass')[0]
        tag_ind   = np.where(np.char.strip(part_names)=='tag')[0]
        type_ind  = np.where(np.char.strip(part_names)=='type')[0]
        ct_ind    = np.where(np.char.strip(part_names)=='creation_time')[0]
        posx_ind  = np.where(np.char.strip(part_names)=='posx')[0]
        posy_ind  = np.where(np.char.strip(part_names)=='posy')[0]
        posz_ind  = np.where(np.char.strip(part_names)=='posz')[0]
        velx_ind  = np.where(np.char.strip(part_names)=='velx')[0]
        vely_ind  = np.where(np.char.strip(part_names)=='vely')[0]
        velz_ind  = np.where(np.char.strip(part_names)=='velz')[0]
        
        # Get the particle types
        type_part = part_data[:,type_ind].flatten()
        
        if (particle_type == 'star' or particle_type == 'stars'):
            # Find only star particles (particle_type=1)
            star_ind = np.where(type_part==1)[0]
        elif (particle_type == 'sink' or particle_type == 'sinks'):
            # Find only sink particles (particle_type=2)
            star_ind = np.where(type_part==2)[0]
        else:
            star_ind = np.arange(len(type_part))
            
        # Masses of stars
        mass = part_data[star_ind,mass_ind]
        # Tags of stars
        tags = part_data[star_ind,tag_ind]
        # Creation times of stars
        ct   = part_data[star_ind,ct_ind]
        # Positions
        posx = part_data[star_ind,posx_ind]
        posy = part_data[star_ind,posy_ind]
        posz = part_data[star_ind,posz_ind]
        # Velocities
        velx = part_data[star_ind,velx_ind]
        vely = part_data[star_ind,vely_ind]
        velz = part_data[star_ind,velz_ind]
        pp = np.array([posx, posy, posz]).T
        pv = np.array([velx, vely, velz]).T
    else:
        print "Error: No particles found in this file!"
        mass = tags = ct = posx = posy = posz = velx = vely = velz = None
        
    return mass, tags, ct, pp, pv
    

def get_ps_from_yt(pltfile, particle_type="all"):
    
    """Take a FLASH plotfile and return an AMUSE
       particle set.
       
       Keyword arguments:
       pltfile      -- The FLASH hdf5 plot or checkpoint file.
                       Note a plotfile means there must also be
                       a particle file with the same ID number.
       particle_type -- the type of particle that you
                         want from FLASH. Can be:
                         star -- active FLASH particle type
                         sink -- sink FLASH particle type
                         any other string or none returns 
                         all particle types.
                         The default is all
       
       Returns:
       stars -- An AMUSE particle set.
    """
    import yt
    ds = yt.load(pltfile)
    dd = ds.all_data()
    
    all_parts_mass = dd['particle_mass'].in_units('Msun')
    all_parts_pos  = dd['particle_position'].in_units('cm')
    all_parts_vel  = dd['particle_velocity'].in_units('cm/s')
    
    if (particle_type == 'star' or particle_type == 'stars'):
        # Index for only the star particles, and not including sink particles.
        mass_ind = np.where(dd['particle_type'].v == 1)
    elif (particle_type == 'sink' or particle_type == 'sinks'):
        mass_ind = np.where(dd['particle_type'].v == 2)
    else:
        mass_ind = np.arange(len(all_parts_mass))
        
    all_parts_mass = all_parts_mass[mass_ind]
    all_parts_pos  = all_parts_pos[mass_ind]
    all_parts_vel  = all_parts_vel[mass_ind]
        
    num_parts = len(all_parts_mass.v)
    
    stars = Particles(num_parts)
    stars.mass     = all_parts_mass.v | u.MSun
    stars.position = all_parts_pos.v  | u.cm
    stars.velocity = all_parts_vel.v  | u.cm/u.s
    stars.index_of_the_particle = range(num_parts)
    
    return stars


def get_energies_from_particles(particles, pm=None, pp=None, pv=None):

    """Get the total energy, relative kinetic energy, and potential energy
    between every particle with respect to all other particles in the set.
    
    Generally takes an AMUSE particle set, but you can also pass directly
    the positions, masses and velocities of the particles if you pass
    particles == None.
    
    Returns matricies of the energies of each particle (on the 0 axis)
    with every other particle (on the 1 axis). The diagonial is zeroed.
    To get the total energy of the set, simply sum along one axis.
    
    Keywords:
    particles -- AMUSE particle set (set to None if using pm, pp, and pv)
    pm -- particle masses     (default = None)
    pp -- particle positions  (default = None)
    pv -- particle velocities (default = None)
    
    Returns:
    TE     -- total relative energies
    rel_ke -- relative kinectic energies
    pu     -- potential energies
    """
    
    gc  = u.constants.G.value_in(u.cm**3 * u.g**-1 * u.s**-2)
    
    if (particles is None and (pm is None 
                               or pp is None 
                               or pv is None)):
        print "Error: particles and pm/pp/pv cannot be None."
    elif (particles is not None
          and (pm is not None
               or pp is not None 
               or pv is not None)):
        print "Error: received both an AMUSE \
               particle set and values for pm/pp/pv. I'm not \
               sure who to believe here!"
    elif (particles is not None):
        pm = particles.mass.value_in(u.g)
        pp = particles.position.value_in(u.cm)
        pv = particles.velocity.value_in(u.cm/u.s)
        
    lp = len(pm)

    print "Number of stars = {}".format(lp)

    # Particle KE
    pke = 0.5*pm*np.sum(pv**2.0, axis=1)

    # Relative particle KE
    # center of mass velocities

    rel_pke = np.zeros((lp,lp))

    for i in range(lp):
        for j in range(i+1,lp):

            vcom_x_ij = (pm[i]*pv[i,0] + pm[j]*pv[j,0])/(pm[i]+pm[j])
            vcom_y_ij = (pm[i]*pv[i,1] + pm[j]*pv[j,1])/(pm[i]+pm[j])
            vcom_z_ij = (pm[i]*pv[i,2] + pm[j]*pv[j,2])/(pm[i]+pm[j])

            # rel vel for ith particle
            rel_vx = pv[i,0] - vcom_x_ij
            rel_vy = pv[i,1] - vcom_y_ij
            rel_vz = pv[i,2] - vcom_z_ij
            # rel kinetic energy from the center of velocity frame for ith
            rel_pke[i,j] = 0.5*pm[i]*np.sum(rel_vx**2.0+rel_vy**2.0+rel_vz**2.0)
            # rel vel for jth particle
            rel_vx = pv[j,0] - vcom_x_ij
            rel_vy = pv[j,1] - vcom_y_ij
            rel_vz = pv[j,2] - vcom_z_ij
            # rel kinetic energy from the center of velocity frame for jth
            rel_pke[j,i] = 0.5*pm[j]*np.sum(rel_vx**2.0
                                           +rel_vy**2.0
                                           +rel_vz**2.0)

    # Particle-Particle PE
    r = cdist(pp,pp)
    masses = np.outer(pm,pm)
    pu = gc*masses/r
    np.fill_diagonal(pu,0)
    #ppe = 0.5*yt.units.G.v*(pu).sum(axis=0)
    
    TE = rel_pke - pu

    print "Max total energy = {}".format(TE.max())
    print "Min total energy = {}".format(TE.min())
    
    return TE, rel_pke, pu
    

def find_bound_relations(TE):
    
    """Given a total relative energy matrix that has total energy
    TE[i,j] for the ith star's total energy relative to the jth
    star, find any jth star the ith star might be bound to and
    store in mul_im_bnd_to and any jth star bound to the ith
    star and store in mul_bnd_to_me.
    
    Keywords:
    
    TE  -- total relative energy matrix
    
    Returns:
    
    mul_im_bnd_to -- dictionary whose key is the index of the ith
                     star and whose values are the indices of all
                     stars the ith star is bound to.
    mul_bnd_to_me -- dictionary whose key is the index of the ith
                     star and whose values are the indices of all
                     stars that are bound to the ith star.
    mutally_bnd   -- dictionary whose key is the ith star and
                     whose (single) value is one star to which
                     it is both most bound and that is most bound
                     to it.
    """
    
    # Note TE[:,i] in this loop returns the stars this star index/key
    # is bound to while TE[i,:] in this loop returns all stars bound
    # this index/key.
    mul_im_bnd_to = {}
    mul_bnd_to_me = {}

    for i in range(len(TE[:,0])):

        bnd_ind = np.where(TE[:,i] < 0.0)[0]
        if len(bnd_ind > 0):
            mul_im_bnd_to[i] = bnd_ind[np.argsort(TE[bnd_ind,i])]

        bnd_ind = np.where(TE[i,:] < 0.0)[0]
        if len(bnd_ind > 0):
            mul_bnd_to_me[i] = bnd_ind[np.argsort(TE[i,bnd_ind])]
            
    # Now lets find any pairs that are mutually most bound.
    mutally_bnd = {}

    for key, val in mul_bnd_to_me.iteritems():
        # Do we love each other?
        other_key = val[0]
        #print other_key
        if (other_key in mul_bnd_to_me.keys()
            and other_key not in mutally_bnd.keys()):
            
            if (mul_bnd_to_me[other_key][0] == key):
                mutally_bnd[key]=val[0]
            
    return mul_im_bnd_to, mul_bnd_to_me, mutally_bnd
    

def check_perturbed_binary(star1, star2, pert):
    

    """Given a binary of star1 and star2 and a third perturber pert,
    determine the perturbation on the binary by comparing the
    difference in the accelerations of star1 by pert and star2
    by pert to the accelerations of star1 and star2 on each other.
    
    Keywords:
    star1 -- first member of binary
    star2 -- second member of binary
    pert  -- perturbing star
    
    Returns:
    gamma -- the ratio of perturb acceleration to member acceleration
    """
    
    gc  = u.constants.G.value_in(u.cm**3 * u.g**-1 * u.s**-2)
    
    s1m = star1.mass.value_in(u.g)
    s1p = star1.position.value_in(u.cm)
    
    s2m = star2.mass.value_in(u.g)
    s2p = star2.position.value_in(u.cm)
    
    pm = pert.mass.value_in(u.g)
    pp = pert.position.value_in(u.cm)
    # Individual accelerations of the binary members on
    # each other.
    r12vec = s1p-s2p
    r12mag = np.sqrt(np.sum((s1p-s2p)**2.0))
    r21vec = s2p-s1p
    r21mag = np.sqrt(np.sum((s2p-s1p)**2.0))
    
    r32vec = pp-s2p
    r32mag = np.sqrt(np.sum((pp-s2p)**2.0))
    r31vec = pp-s1p
    r31mag = np.sqrt(np.sum((pp-s1p)**2.0))
    
    accel_1on2 = gc * s1m * r12vec / r12mag**3.
    accel_2on1 = gc * s2m * r21vec / r21mag**3.
    # Average binary mutual acceleration.
    mutual_accel = gc * 0.5*(s1m+s2m) / np.sum((s2p-s1p)**2.0)
    
    # Direct acceleration caluclations. This makes no assumptions about
    # the orbits, it just calculates the tides at the current binary
    # positions.
    # Pertuber acceleration on star1
    ps1a = gc *  pm * r31vec / r31mag**3.
    # Perturber acceleration on star2
    ps2a = gc *  pm * r32vec / r32mag**3.
        
    # Gamma based on averaging the internal binary acceleration
    gamma = (  np.sqrt(np.sum((ps2a - ps1a)**2.0)) 
            /  np.sqrt(np.sum((mutual_accel)**2.0))  )
    
    return gamma
    
def find_all_binaries(stars, limiting_gamma=0.01):
    
    """Given a particle set, find all mutally bound binaries
    that are not perturbed by a gamma larger than limiting_gamma.
    Return the particle set with these binaries appended to the
    set as stars.bp, where bp is None if this star has no
    binary partner or equals the index of the binary partner.
    
    Keyword arguments:
    stars          -- AMUSE particle set
    limiting gamma -- the limiting ratio of perturber to binary
                      accelerations.
                      
    Returns:
    stars          -- AMUSE particle set with the binary
                      indicies appended
    """
    
    # First calculate all the relative energies of all stars
    # to each other.
    TE, rel_pke, pu = get_energies_from_particles(stars)
    
    # Now find all first pass mutually bound pairs.
    mul_im_bnd_to, mul_bnd_to_me, mutually_bnd = \
                            find_bound_relations(TE)
    
    # Now test each pair for perturbers.
    
    check_binaries_for_perturbers(stars,
                                  mutually_bnd,
                                  limiting_gamma=limiting_gamma)
                
    return stars, mutually_bnd

def check_binaries_for_perturbers(stars,
                                  mutually_bnd,
                                  limiting_gamma=0.01):
    """Given an AMUSE particle set and a dictionary with mutually
    bound key-value pairs, check all the mututally bound objects
    to see if any neighbors tidal accelerations on the binary pair
    are larger than limiting_gamma*binary mutual acceleration.

    If the tidal acceleration is too large, remove the binary key-
    value pair from the dictionary. Otherwise, place the index of
    the companion star in the primary star.bp and vice versa.

    Keyword arguments
    stars          -- AMUSE particle set
    mutually_bnd   -- dictionary of binary indices from stars as key-value pairs
    limiting gamma -- the limiting ratio of perturber to binary
                     accelerations.
    """
    stars.bp = None
    rm_keys = []
    for key, val in mutually_bnd.iteritems():
        
        star1 = stars[key]
        star2 = stars[val]
        
        other_stars = stars - star1 - star2
        
        perturbed_too_much = False

        for pert in other_stars:
            
            gamma = check_perturbed_binary(star1, star2, pert)
            if (gamma > limiting_gamma):
                perturbed_too_much = True
                #print "gamma = {}".format(gamma)
                break

        if (perturbed_too_much):
            # Add this key to a list to be removed
            # from the mutually_bnd dictionary.
            rm_keys.append(key)

        else:
            # Add this pair to the stars list.
            stars[key].bp = val
            stars[val].bp = key
        
    # Remove from mututally_bnd dict
    for k in rm_keys:
        mutually_bnd.pop(k, None)
        
def get_binary_number_fractions(stars, mutually_bnd, mass_bin_edges):

    """Given an AMUSE particle set, a dictionary of mutually bound
    pairs, and a list of the edges of the required mass bins,
    calculate:

    1. frac_in_bin:     the number fraction of stars in binaries per mass bin.
    2. bin_num_in_bin:  the number of binaries in each mass bin.
    3. num_in_bin:      total number of stars (including binaries) in each mass bin.
    4. bin_mass_in_bin: the summed binary mass in each mass bin.
    5. mass_in_bin:     total stellar mass (including binaries) in each mass bin.

    And return all those things.

    Keyword arguments:
    stars          -- AMUSE particle set
    mutually_bnd   -- dictionary of binary indices from stars as key-value pairs
    limiting gamma -- the limiting ratio of perturber to binary
                     accelerations.
    mass_bin_edges -- The values of the mass bin edges, including the last right edge (len=nbins+1).

    Returns:
    Items 1-5 above.
    """
    
    nbins           = len(mass_bin_edges)-1
    pm              = stars.mass.value_in(units.MSun)
    frac_in_bin     = np.zeros(nbins)
    mass_in_bin     = np.zeros(nbins) # total mass in bin
    num_in_mass_bin = np.zeros(nbins) # total # in bin
    bin_mass_in_bin = np.zeros(nbins) # binary mass in bin
    bin_num_in_bin  = np.zeros(nbins) # binary number in bin
    #mass_bin_edges  = np.zeros(nbins+1)

    mass_sorted_ind    = np.argsort(pm) #[::-1]
    sorted_m           = pm[mass_sorted_ind]
    tot_m              = sorted_m.sum()
    frac_m             = sorted_m/tot_m
    cul_m              = np.cumsum(frac_m)
    last_upper_bin_ind = 0

    #print cul_m

    #mass_bin_edges[:] = bin_edges[:]*tot_m/ms

    for i in range(nbins):

        #upper_bin_ind = np.where((bin_edges[i] < cul_m) &
        #                         (cul_m <= bin_edges[i+1]))[0]

        upper_bin_ind = np.where((mass_bin_edges[i] < sorted_m) &
                         (sorted_m <= mass_bin_edges[i+1]))[0]

        #print pm[mass_sorted_ind[upper_bin_ind]]
        
        # total number of all stars in this mass bin
        num_in_mass_bin[i] = float(len(upper_bin_ind))
        # the total mass in this bin
        mass_in_bin[i]     = sorted_m[upper_bin_ind].sum()
        #if (np.size(upper_bin_ind) > 0):
        #    mass_bin_edges[i+1] = sorted_m[upper_bin_ind[-1]]
        #else:
        #    mass_bin_edges[i+1] = mass_bin_edges[i]
        #last_upper_bin_ind = upper_bin_ind
        #print mass_sorted_ind[upper_bin_ind]

        cnt_m = 0.

        for key, val in mutually_bnd.iteritems():
            
            if (key in mass_sorted_ind[upper_bin_ind]):
                cnt_m += 1.
                bin_mass_in_bin[i] += pm[key]
            #if (val in mass_sorted_ind[upper_bin_ind]):
            #    cnt_m += 1.
            #    bin_mass_in_bin[i] += pm[val]

        if (cnt_m > 0):
            frac_in_bin[i] = cnt_m/num_in_mass_bin[i]
            bin_num_in_bin[i] = cnt_m
            
    return (frac_in_bin, bin_num_in_bin, num_in_mass_bin,
            bin_mass_in_bin, mass_in_bin)

def get_binary_mass_ratios(stars, mutually_bnd):
    """Given an AMUSE particle set and a dictionary with mutually
    bound key-value pairs, calculate the mass ratio q = m2/m1
    where m1 is the most massive (primary) star and m2 is the smaller
    (companion) star.

    Keyword arguments:
    stars          -- AMUSE particle set
    mutually_bnd   -- dictionary of binary indices from stars as key-value pairs

    Returns:
    q  -- Mass ratios of the binary
    s1 -- Masses of all the primaries
    s2 -- Masses of all the companions.
    """
    pm = stars.mass.value_in(units.MSun)
    nb = len(mutually_bnd)
    
    q  = np.zeros(nb)
    s1 = np.zeros(nb)
    s2 = np.zeros(nb)
    i  = 0

    for key, val in mutually_bnd.iteritems():

        m1 = pm[key]
        m2 = pm[val]

        # s1 is the primary, s2 is the companion.
        if (m1 > m2):
            s1[i] = m1
            s2[i] = m2
        else:         
            s2[i] = m1
            s1[i] = m2
        
        i = i+1    
    
    q = s2/s1
    return q, s1, s2

def get_binary_separations(stars, mutually_bnd):
    """Given an AMUSE particle set and a dictionary with mutually
    bound key-value pairs, calculate the binary separation distances
    and return that value in AU.

    Keyword arguments:
    stars          -- AMUSE particle set
    mutually_bnd   -- dictionary of binary indices from stars as key-value pairs

    Returns:
    bin_sep        -- Binary separation distances in AU (without the AMUSE units).
    """
    # Returning this in AU.
    pp      = stars.position.value_in(units.AU)
    r       = cdist(pp,pp)
    bin_sep = r[mutually_bnd.keys(),mutually_bnd.values()]
    
    return bin_sep

'''
### This is a test!

from amuse.lab import *

from amuse.community.fractalcluster.interface import MakeFractalCluster
from amuse.community.fractalcluster.interface import new_fractal_cluster_model

masses = new_kroupa_mass_distribution(1000)
conv = nbody_system.nbody_to_si(masses.sum(), 1.0 | units.parsec)
stars = new_fractal_cluster_model(convert_nbody=conv, 
                           masses=masses,
                           do_scale=False,
                          virial_ratio=0.01)

#stars = read_set_from_file('test_binary_stars.hdf5', format='hdf5')

stars, mutually_bnd = find_all_binaries(stars, limiting_gamma=0.01)

num_of_bin = len(mutually_bnd.keys())
print "Total number of binaries = {}".format(num_of_bin)

max_mass = np.max(stars.mass.value_in(units.MSun))
min_mass = np.min(stars.mass.value_in(units.MSun))
nbins          = 15
mass_bin_edges = np.logspace(np.log10(min_mass), np.log10(max_mass), nbins+1)

bin_num_in_bin   = np.zeros(nbins)
num_in_bin       = np.zeros(nbins)
bin_mass_in_bin  = np.zeros(nbins)
mass_in_bin      = np.zeros(nbins)
frac_in_bin      = np.zeros(nbins)

frac_in_bin, bin_num_in_bin, \
num_in_bin, bin_mass_in_bin, \
mass_in_bin = get_binary_number_fractions(
                     stars, mutually_bnd, mass_bin_edges)

import matplotlib

font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 24}

matplotlib.rc('font', **font)

import matplotlib.pyplot as plt

lw = 3
fs = 20 

f = plt.figure(figsize=(8,8))
ax = f.add_subplot(111)
ax.tick_params(axis="both", which="both", direction="in")
ax.tick_params(axis="both", which="both", width=1)
ax.tick_params(axis="both", which="major", length=10)
ax.tick_params(axis="both", which="minor", length=5)
        
ax.plot(mass_bin_edges[1:],frac_in_bin,lw=lw)
plt.ylabel("Fraction of Stars in Binaries")
plt.xlabel(r"Primary Mass ($\rm M_{\odot}$)")
plt.xscale('log')
ax.annotate(s='# of binaries = {}'.format(int(num_of_bin)),
             xy=(0.2,0.88), xycoords='figure fraction', color='k')
plt.show()

f = plt.figure(figsize=(8,8))
ax = f.add_subplot(111)
ax.tick_params(axis="both", which="both", direction="in")
ax.tick_params(axis="both", which="both", width=1)
ax.tick_params(axis="both", which="major", length=10)
ax.tick_params(axis="both", which="minor", length=5)

ax.plot(mass_bin_edges[1:],bin_num_in_bin, lw=lw)
plt.ylabel("Number of Stars in Binaries")
plt.xlabel(r"Primary Mass ($\rm M_{\odot}$)")
plt.xscale('log')
plt.show()
'''

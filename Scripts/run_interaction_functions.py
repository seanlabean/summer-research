from interaction import bbh_orbit
from interaction import impact_parameter
from interaction import bbh_hvgc
from interaction import integrate_system
from interaction import integrate_plotting
from interaction import plot_track
from interaction import bbh_phase_loop
from interaction import integrate_energy
from interaction import plot_energy
from amuse.lab import Particles, units
from amuse.units import constants
import math
bbh_mass = 7e9 | units.MSun
mass_ratio = 0.1
separation = 3 | units.parsec
bh_phase = 23*math.pi/15
gc_closest = 5.5 | units.parsec
gc_vinf = 500 | (units.km/units.s)

particles = bbh_hvgc(bbh_mass, mass_ratio, separation, bh_phase, gc_closest, gc_vinf)
end_position = 120 | units.parsec

#Final Velocity vs. Phase#
#bbh_phase_loop(bbh_mass, mass_ratio, separation, bh_phase, gc_closest, gc_vinf, end_position,'phase_loop.png') 

#GC Behavior#
hvgc_vel, xgc, ygc, xbh1, ybh1, xbh2, ybh2 = integrate_plotting(particles, end_position)
plot_track(xgc, ygc, xbh1, ybh1, xbh2, ybh2, "hvgc1.png")

#plot_min(bbh_mass, mass_ratio, separation, bh_phase, gc_closest, gc_vinf, end_position, 'bh_min_dist.png')


#Particle Energy#
#bh1_energy, bh2_energy, bbh_energy, hvgc_energy, total_energy, time  = integrate_energy(particles, 120 | units.parsec)
#plot_energy(bh1_energy, bh2_energy, bbh_energy, hvgc_energy, total_energy, time, 'unbound_energy.png')

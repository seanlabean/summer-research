from amuse.lab import Particles, units
from amuse.units import constants

def integrate_system_final(particles, end_position):
    from amuse.lab import ph4, nbody_system
    convert_nbody = nbody_system.nbody_to_si(particles.mass.sum(), particles[0].position.length())
    gravity = ph4(convert_nbody)
    gravity.particles.add_particles(particles)

    hvgc = gravity.particles[2]
    bh1 = gravity.particles[0]
    bh2 = gravity.particles[1]
    hvgc_energy = []
    hvgc_pos = (hvgc.x**2 + hvgc.y**2 + hvgc.z**2).sqrt()
    while hvgc_pos < end_position:
        captured = False
        gravity.evolve_model(gravity.model_time + (2000 | units.yr))
       # hvgc_v = (hvgc.vx**2 + hvgc.vy**2 + hvgc.vz**2).sqrt()
       # hvgc_pos = (hvgc.x**2 + hvgc.y**2 + hvgc.z**2).sqrt()
       # energy = 0.5*hvgc.mass*hvgc_v**2 - constants.G*(bh1.mass + bh2.mass)*hvgc.mass/hvgc_pos
       # hvgc_energy.append(energy.value_in(units.J))
        hvgc_pos = (hvgc.x**2 + hvgc.y**2 + hvgc.z**2).sqrt()
        if gravity.model_time > (300000 | units.yr):
            hvgc_pos = end_position
            captured = True
            break
    hvgc_vel = (hvgc.vx**2 + hvgc.vy**2 + hvgc.vz**2).sqrt()
    if captured == True:
        hvgc_vel = 0 | units.kms
    gravity.stop()
    return hvgc_vel, hvgc_energy

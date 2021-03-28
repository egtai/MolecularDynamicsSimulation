import numpy as np
from .LJ_forces import *
from .rescale import *
from .observables import *

def vervlet_evolution(pos, vel, time_step, steps, mass, epsilon, sigma, temp, desired_temp, size, relaxation_steps):
    """
    Execution of the Verlet simulation algorithm for a given system.
    Initially, the system is rescaled until it reaches the desired temperature.
    After equilibrium, observables are calculated.

    Parameters:
    -----------
    pos : np.ndarray
        The positions of the particles in cartesian space. pos.shape = (N,3)
    vel: np.darray
        The velocity vectors of each particle. vel.shape = (N,3)
    time_step: float
        Duration of a single time step
    steps: int
        Total number of time steps
    mass: float
        Mass of a particle
    epsilon: float
        Energy scale LJ
    sigma: float
        Length scale LJ
    temp : float
        The (unitless) initial temperature of the system.
    desired_temp : float
        The (unitless) final temperature of the system.
    size: float
        Dimension of the simulation domain
    relaxation_steps: int
        Number of steps at which temperature is evaluated to be rescaled

    Return:
    -------
    evolution_data: np.darray
        trajectories: np.ndarray
            Position of all particles during each time step of the simulation
        velocities: np.ndarray
            Velocity of all particles during each time step of the simulation

    energy_data:
        kinetic: np.ndarray
            Kinetic energy of each particle w.r.t. time
        potential: np.ndarray
            Potential energy of each particle w.r.t. time

    observables_data: np.darray

        pressure_data: np.ndarray
            pressure: np.darray
                Pressure w.r.t. time
            af_pressure: np.darray
                Pressure autocorrelation w.r.t. time
            error_pressure: float
                Error associated to pressure

        diffusion_data: np.ndarray
            diffusion: np.darray
                Diffusion constant w.r.t. time
            af_diffusion: np.darray
                Diffusion autocorrelation function w.r.t. time
            error_diffusion: float
                Error associated to diffusion

        pair_correlation_data: np.ndarray
            pair_correlation: np.ndarray
                Pair correlation function w.r.t. distance
            bins_pc: np.ndarray
                Distances associated to pair correlation function
    """
    #----------------------------EVOLUTION-------------------------------------#
    #Initialize arrays for evolution of position, velocity and energy.
    n_atoms = len(pos)
    trajectories = np.zeros((steps, n_atoms, 3))
    velocities = np.zeros((steps, n_atoms, 3))

    kinetic = np.zeros((steps, n_atoms),dtype=float)
    potential = np.zeros((steps, n_atoms),dtype=float)

    #--------------------------OBSERVABLES-------------------------------------#
    #Initialize arrays to save data of observables over time
    pressure = np.zeros((steps),dtype=float)
    diffusion = np.zeros((steps),dtype=float)

    #Number of pieces in which we divide the histogram of pair_correlation
    nbins = 200
    histogram_pc = np.zeros((nbins))
    bins_pc = np.zeros((nbins+1))

    #System starts out of equilibrum
    equilibrium = False
    #Save initial position for calculation of difussion constant
    initial_position = pos

    #--------------------------RUN-SIMULATION----------------------------------#
    for time in range(steps):

        #Save current state of the system
        trajectories[time,:,:] = pos
        velocities[time,:,:] = vel

        #Calculate energy
        kinetic[time,:], potential[time,:] = get_energies_particles(pos, vel, size)
        #Relaxation and rescaling
        std_dev = np.sqrt(temp)

        lamb, equilibrium = rescale(equilibrium=equilibrium,
                                    time=time,
                                    relaxation_steps=relaxation_steps,
                                    kinetic_energy=kinetic,
                                    T=temp,
                                    desired_T=desired_temp,
                                    N=n_atoms,
                                    std_dev=std_dev,
                                    vel=vel)

        #Calculate initial force
        F_old = lj_force(position=pos,
                         box_dim=size)
        #Update positions
        pos = vervlet_update_positions(init_velocities=vel,
                                       init_position=pos,
                                       time_step=time_step,
                                       mass=mass,
                                       sigma=sigma,
                                       epsilon=epsilon,
                                       F_old=F_old,
                                       box_dim=size)
        #Calculate force in updated positions
        F_new = lj_force(position=pos,
                         box_dim=size)
        #Update velocities with F_old and F_new
        vel = lamb * vervlet_update_velocities(init_velocities=vel,
                                               init_position=pos,
                                               time_step=time_step,
                                               mass=mass,
                                               sigma=sigma,
                                               epsilon=epsilon,
                                               F_old=F_old,
                                               F_new=F_new,
                                               box_dim=size)

        #Save observables once equilibrium is reached
        if equilibrium:
            #Calculate relative distances in the system
            _, dist = atomic_distances(pos=pos, box_dim=size)
            #Get pressure
            pressure[time] = get_pressure(T=desired_temp,
                                          sigma=sigma,
                                          distances=dist,
                                          box_dim=size)
            #Get diffusion constant
            diffusion[time] = get_difussion_constant(init_pos=initial_position,
                                                     current_pos=pos)
            #Update data of pair correlation function
            histogram_pc, bins_pc = update_pair_correlation(histogram=histogram_pc,
                                                            bins=bins_pc,
                                                            nbins=nbins,
                                                            distances=dist,
                                                            n_atoms=n_atoms)

    #Format data for output
    trajectories = order_evolution_data(data=trajectories,
                                        n_atoms=n_atoms)

    velocities =  order_evolution_data(data=velocities,
                                       n_atoms=n_atoms)

    #Pressure/pressure autocorrelation/pressure error
    pressure = np.ma.masked_equal(pressure,0).compressed()
    af_pressure = get_autocorrelation_function(pressure)
    error_pressure = get_error_observable(pressure)
    pressure_data= [pressure,af_pressure,error_pressure]

    #Diffusion/diffusion aurtocorrelation/diffusion error
    diffusion = np.ma.masked_equal(diffusion,0).compressed()
    af_diffusion = get_autocorrelation_function(diffusion)
    error_diffusion = get_error_observable(diffusion)
    diffusion_data = [diffusion,af_diffusion,error_diffusion]
    #Pair correlation function
    pair_correlation = get_pair_correlation(histogram=histogram_pc,
                                            bins=bins_pc,
                                            n_atoms=n_atoms,
                                            box_dim=size)
    pair_correlation_data = [pair_correlation,bins_pc]

    evolution_data = [trajectories, velocities]
    energy_data = [kinetic, potential]
    observables_data = [pressure_data, diffusion_data, pair_correlation_data]

    return evolution_data, energy_data, observables_data


def vervlet_update_positions(init_velocities, init_position, time_step, mass, sigma, epsilon, F_old, box_dim):
    """
    Update positions according to vervlet algorithm

    Parameters:
    -----------
        init_pos: np.darray
            Initial atomic configuration
        init_vel: np.darray
            Initial velocity distribution
        time_step: float
            Time interval

    Returns:
    --------
        np.darray
            Updated positions
    """
    new_positions = ( init_position + init_velocities*time_step + F_old*time_step**2/2 )%box_dim
    return new_positions

def vervlet_update_velocities(init_velocities, init_position, time_step, mass, sigma, epsilon, F_old, F_new, box_dim):
    """
    Update velocities according to vervlet algorithm

    Parameters:
    -----------
        init_pos: np.darray
            Initial atomic configuration
        init_vel: np.darray
            Initial velocity distribution
        time_step: float
            Time interval

    Returns:
    --------
        np.darray
            Updated velocities
    """
    new_velocities = init_velocities + time_step * (F_old+F_new)/2
    return np.nan_to_num(new_velocities)




def order_evolution_data(data, n_atoms):
    """
    Transforms input data from (time, atom, coordinate) to (atom, coordinate, time)

    Parameters:
    ----------
        data: np.darray
            Data from evolution in the format (time, atom, coordinate)
        n_atoms: np.int
            Number of atoms, correspond to dimension atom in data.

    Returns:
    --------
        np.darray
            Data ordered in format (atom, coordinate, time)
    """

    data_ordered = [[] for i in range(n_atoms)]

    for atom in range(n_atoms):
        data_ordered[atom] = [data.T[j][atom] for j in range(3)]

    return np.array(data_ordered)

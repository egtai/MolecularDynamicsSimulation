import numpy as np
from .LJ_forces import *
from .rescale import *

def euler_evolution(init_pos, init_vel, time_step, steps, mass, sigma, epsilon, temp, desired_temp, size, relaxation_steps=50):
    """
    Implements the Euler algorithm for the evolution of the positions and velocities of all the particles in the system.

    Parameters:
    -----------
        init_pos: np.darray
            Initial atomic configuration
        init_vel: np.darray
            Initial velocity distribution
        time_step: float
            Time interval
        steps: int
            Number of time steps

    Returns:
    --------
        np.darray
            trajectories: array of size (n_atoms, coordinates, steps) containing the evolution of position
        np.darray
            velocities: array of size (n_atoms, coordinates, steps) containing the evolution of velocity
    """

    n_atoms = len(init_pos)
    trajectories = np.zeros((steps, n_atoms, 3))
    velocities = np.zeros((steps, n_atoms, 3))

    kinetic = np.zeros((steps, n_atoms),dtype=float)
    potential = np.zeros((steps, n_atoms),dtype=float)

    equilibrium = False

    for time in range(steps):
        #Save state of the system
        trajectories[time,:,:] = init_pos
        velocities[time,:,:] = init_vel

        #Calculate total energy
        kinetic[time,:], potential[time,:] = get_energies_particles(init_pos, init_vel, size)
        #Relaxation and rescaling
        std_dev = np.sqrt(temp) # (without reduced units -> np.sqrt(temp*kB/mass))
        #lamb, equilibrium = relaxation_from_equilibrium(equilibrium, time, relaxation_steps, kinetic, temp, n_atoms, std_dev)
        lamb, equilibrium = rescale(equilibrium,
                                    time,
                                    relaxation_steps,
                                    kinetic,
                                    temp,
                                    desired_temp,
                                    n_atoms,
                                    std_dev,
                                    init_vel)

        #Update positions and velocities
        init_pos = euler_update_positions(init_vel, init_pos, time_step, mass, sigma, epsilon, 0, 0, size)
        init_vel = lamb * euler_update_velocities(init_vel, init_pos, time_step, mass, sigma, epsilon, 0, 0, size)

    #Reorder data
    trajectories = order_evolution_data(trajectories, n_atoms)
    velocities = order_evolution_data(velocities, n_atoms)

    return [trajectories, velocities], [kinetic, potential]


def euler_update_positions(init_velocities, init_position, time_step, mass, sigma, epsilon, F_old, F_new, box_dim):
    """
    Update positions according to Euler algorithm

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
    new_positions = (init_position + init_velocities*time_step)%box_dim
    return new_positions

def euler_update_velocities(init_velocities, init_position, time_step, mass, sigma, epsilon, F_old, F_new, box_dim):
    """
    Update velocities according to Euler algorithm

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
    new_velocities = init_velocities + lj_force(init_position, box_dim) * time_step
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

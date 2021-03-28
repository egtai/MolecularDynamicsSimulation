import numpy as np
from scipy.constants import k as kB

def unit_vectorise(rel_pos, rel_dist):
    return np.nan_to_num(rel_pos/rel_dist[:,:,None])

def lj_force(position, box_dim):
    """
    Calculates the net forces on each atom.

    Parameters
    ----------
    rel_pos : np.ndarray
        Relative particle positions as obtained from atomic_distances
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    np.ndarray
        The net force acting on particle i due to all other particles
    """
    rel_pos, rel_dist = atomic_distances(position, box_dim)
    uv_rel_pos = unit_vectorise(rel_pos, rel_dist)
    F = dlennard_jones(rel_dist)[:,:,None] * uv_rel_pos
    return np.sum( np.nan_to_num(F) , axis=1 )

def total_energy(pos, vel, rel_dist, box_dim):
    _, rel_dist = atomic_distances(pos, box_dim)
    return kinetic_energy(vel), potential_energy(rel_dist)

def kinetic_energy(vel):
    """
    Computes the kinetic energy of an atomic system.

    Parameters
    ----------
    vel: np.ndarray
        Velocity of particle

    Returns
    -------
    float
        The total kinetic energy of the system.
    """

    # E_kin_particle = np.sum(vel**2,axis=0)/2
    E_kin_particle = np.linalg.norm(vel, ord=2, axis=1)**2/2

    return E_kin_particle

def get_energies_particles(pos, vel, box_dim):
    _, rel_dist = atomic_distances(pos, box_dim)
    e_kin = kinetic_energy(vel)
    e_pot = potential_energy(rel_dist)
    return e_kin, np.sum(e_pot,axis=1)

def potential_energy(rel_dist):
    """
    Computes the potential energy of an atomic system.

    Parameters
    ----------
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    float
        The total potential energy of the system.
    """

    E_pot_particles = lennard_jones(rel_dist)

    return E_pot_particles

def lennard_jones(rel_dist):
    """
    Calculates the lennard jones potential for an array of distances.

    Parameters
    ----------
    rel_dist : np.ndarray
        Relative distances obtained from atomic_distances
    epsilon: float
        Energy scale LJ
    sigma: float
        Length scale LJ

    Returns
    -------
    np.darray
        potential energy of each atom
    """

    aux = (1/rel_dist)**6
    E_pot = 4*(aux**2-aux)
    E_pot[E_pot == np.inf] = 0
    return np.nan_to_num(E_pot)

def dlennard_jones(rel_dist):
    """
    Calculates the derivative of lennard jones potential w.r.t rᵢ such that r=|√(rᵢ-rⱼ)| for an array of distances.

    Parameters
    ----------
    rel_dist : np.ndarray
        Relative distances obtained from atomic_distances
    epsilon: float
        Energy scale LJ
    sigma: float
        Length scale LJ

    Returns
    -------
    np.darray
        derivative of potential energy of each atom w.r.t rᵢ such that r=|√(rᵢ-rⱼ)|.
    """
    rel_dist2 = rel_dist**2
    rel_dist6 = rel_dist**6
    rel_dist8 = rel_dist2 * rel_dist6
    rel_dist14 = rel_dist6**2 *rel_dist2

    dU_LJ = np.nan_to_num(2*24*rel_dist*(1/rel_dist8 - 2*1/rel_dist14))

    return dU_LJ

def atomic_distances(pos, box_dim):
    """
    Calculates relative positions and distances between particles.
    Assume that np.shape(pos) = (N, 3).

    parameters
    ----------
    pos : np.ndarray
        The positions of the particles in cartesian space. pos.shape = (N,3)
    box_dim : float
        The dimension of the simulation box

    returns
    -------
    rel_pos : np.ndarray
        Relative positions of particles w.r.t. other every other particle
    rel_dist : np.ndarray
        The distance between particles
    """

    """"""
    N = len(pos)
    rel_pos = np.zeros((N, N, 3),dtype=float)

    for i in range(N):
        rel_pos[i, :, :] = (pos - pos[i, :] + box_dim/2) % box_dim - box_dim/2

    rel_dist = np.linalg.norm(rel_pos, ord=2, axis=2)
    return rel_pos, rel_dist

def init_velocity(num_atoms, temp, mass):
    """
    Initializes the system with Gaussian distributed velocities.

    Parameters:
    -----------
    num_atoms : int
        The number of particles in the system.
    temp : float
        The (unitless) temperature of the system.

    Returns:
    --------
    vel_vec : np.ndarray
        Array of particle velocities
    """
    #std_dev = np.sqrt(temp*kB/mass)
    std_dev = np.sqrt(temp)
    print(f"Initial temperature: {temp}")
    vel_vec = np.random.normal(loc=0.0,scale=std_dev, size=(num_atoms, 3))

    return vel_vec

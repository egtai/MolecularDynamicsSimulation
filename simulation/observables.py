import numpy as np
from .LJ_forces import *
from scipy.optimize.minpack import curve_fit

def get_density(n_atoms, box_dim):
    """
    Parameters:
    -----------
    box_dim: float
        Dimension of the simulation domain
    n_atoms: int
        Number of atoms in the system

    Return:
    --------
    rho: float
        Density of the system
    """
    rho = (1/box_dim)**3 * n_atoms
    return rho

def get_pressure(T,sigma,distances,box_dim):
    """
    Obtain pressure at a given configuration of the system.

    Parameters:
    -----------
    T: float
        Temperature
    sigma: float
        Length scale LJ
    distances: np.ndarray
        Relative distances between all particles
    box_dim: float
        Dimension of the simulation domain

    Return:
    -------
    float:
        Pressure
    """
    n_atoms = len(distances[0])
    rho = get_density(n_atoms,box_dim)
    ave_factor = 1/2 * np.mean( distances * dlennard_jones(distances) )
    return T * rho * ( 1 - 1/(3*n_atoms*T) * ave_factor )

def get_difussion_constant(init_pos, current_pos):
    """
    Parameters:
    -----------
    init_pos: np.darray
        Initial position of every atom in the system
    current_pos: np.darray
        Current position of every atom in the system
    Return:
    --------
    _: float
        Average distance of all particles from origin
    """
    return np.mean((init_pos-current_pos)**2)

def get_pair_correlation(histogram, bins, n_atoms, box_dim):
    """
    Parameters:
    -----------
    histogram: np.darray
        Histogram of distances between particles summed over all times
    bins: np.darray
        Distances at which histogram is evaluated
    n_atoms: int
        Number of atoms in the system
    box_dim: float
        Dimension of the simulation domain
    Return:
    --------
    _: float
        pair correlation function
    """
    return (2*box_dim**3/(n_atoms*(n_atoms-1)))*histogram/(4*np.pi*(bins[-1]-bins[1])*bins[1:]**2)

def update_pair_correlation(histogram, bins, nbins, distances, n_atoms):
    """
    Parameters:
    -----------
    histogram: np.darray
        Current histogram of distances between particles
    bins: np.darray
        Distances at which histogram is evaluated
    nbins: int
        Number of pieces in which distances are going to be classified
    distances: np.darray
        Relative distances of all the particles in the system at a given time
    n_atoms: int
        Number of atoms in the system
    Return:
    --------
    histogram: np.darray
        updated histogram
    bins: np.darray
        updated bins
    """
    distances=distances[np.triu_indices(n_atoms,k=1)]
    h,bins=np.histogram(distances,bins=nbins,density=True)
    histogram += h
    bins += bins
    return histogram, bins

def autocorrelation_function(observable, time):
    """
    Parameters:
    -----------
    observable: np.ndarray
        data of an observable w.r.t. time
    time: float
        time at which autocorrelation function is computed

    Return:
    autocorrelation: float
        autocorrelation function at a given time
    --------
    """
    #Number of time steps
    N = len(observable)

    #Long time corrections tend to mess up data, this next line ignores corrections after half time
    #observable[N//:]=0

    A_nt = observable[time:N-time:time]
    A_n = observable[:(N-time)//time]
    A_nn = observable[:N-time]

    #Make sure that A_nt and A_n have the same length
    if len(A_nt)!=len(A_n):
        A_n = observable[:(N-time)//time-1]

    #Here the term sum_n A_{n} A_{n+t} was implemented as np.sum(A_nt * A_n)
    #Other terms like sum_n A_n were implemented as np.sum(A_nn)
    #Note that len(A_nn) > len(A_n)
    sigma_A_n = np.sqrt((N-time)*np.sum(A_nn**2)-(np.sum(A_nn)**2))
    sigma_A_nt = np.sqrt((N-time)*np.sum(A_nt**2)-(np.sum(A_nt)**2))
    autocorrelation = ((N-time)*np.sum(A_n*A_nt)-np.sum(A_n)*np.sum(A_nt))/(sigma_A_n*sigma_A_nt)

    return np.abs(autocorrelation)

def get_autocorrelation_function(observable):
    """
    Parameters:
    -----------
    observable: np.ndarray
        data of an observable w.r.t. time

    Return:
    af: np.darray
        autocorrelation function of an observable over all time
    --------
    """
    number_steps = len(observable)-1
    af = np.zeros(number_steps)
    for time in range(1,number_steps):
        af[time]=autocorrelation_function(observable, time)

    return af[1:]

def func(x,a,b,c):
    """
    Parameters:
    -----------
    x: float

    a: float
        magnitude
    b: float
        characteristic distance
    c: float
        shift from zero

    Return:
    --------
    _: float
        function of exponential decay evaluated at x,a,b,c.
    """
    return a*np.exp(-x/b)+c

def get_tau(autocorrelation):
    """
    Parameters:
    -----------
    autocorrelation: np.darray
        autocorrelation function of a given observable

    Return:
    -------
    params: np.darray
        fitting parameters of autocorrelation onto an exponential decay
    """
    N = len(autocorrelation)//2
    x = np.arange(N)
    y = np.nan_to_num(autocorrelation[:N])
    params, errors = curve_fit(func,x, y)
    return params

def get_error_observable(observable):
    """
    Parameters:
    -----------
    observable: np.ndarray
        data of an observable w.r.t. time

    Return:
    --------
    sigma: float
        error associated to a given observable
    """
    observable=np.nan_to_num(observable)
    tau = get_tau(np.nan_to_num(get_autocorrelation_function(observable)))[1]
    N = len(observable)
    sigma = np.sqrt(2*tau/N)*np.sqrt(np.mean(observable**2)-np.mean(observable)**2)
    return sigma

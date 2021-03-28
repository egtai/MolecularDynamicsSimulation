import numpy as np

def rescale(equilibrium,time,relaxation_steps,kinetic_energy,T,desired_T,N,std_dev,vel):
    """
    System described by the canonical ensemble. This function rescales the
    temperature until it reaches the desired value. This state is called equilibrum.

    Parameters:
    -----------
    equilibrium: bool
        Flag indicating equilibrum of the system before rescaling
    time: int
        Current time step
    relaxation_steps: int
        Duration of the relaxation interval
    kinetic_energy: np.ndarray
        Kinetic energy of the system at a given time
    T: float
        Initial temperature
    desired_T: float
        Desired temperature
    N: int
        Number of atoms
    std_dev: float
        Standart deviation of the temperature
    vel: np.ndarray
        Velocities of all particles at a given time

    Returns:
    --------
    lamb: float
        Rescaling parameters
    equilibrium: bool
        Flag indicating equilibrum of the system after rescaling
    """
    nvel = np.linalg.norm(vel, ord=2, axis=1)**2
    rescaling_factor = np.sqrt( 3*desired_T*(N-1)/np.sum(nvel) )
    if time%relaxation_steps==0 and time != 0:
        # average_temperature = np.sum(kinetic_energy[time-relaxation_steps:time,:],axis=1)/relaxation_steps
        average_temperature = 2/3 * np.sum(kinetic_energy[time-relaxation_steps:time,:]) / (relaxation_steps*(N-1))
        print(f"average_temperature: {average_temperature}")
        if not equilibrium and np.abs(desired_T-average_temperature)>0.01:
            print(f"Rescale velocity: {rescaling_factor}")
            lamb = rescaling_factor
        else:
            lamb = 1
            equilibrium = True
    else:
        lamb = 1

    return lamb, equilibrium

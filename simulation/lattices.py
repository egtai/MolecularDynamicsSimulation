import numpy as np

def fcc_lattice(len_x, len_y, len_z, lat_const):
    """
    Initializes a system of atoms on an fcc lattice.

    Parameters
    ----------
    num_atoms : int
        The number of particles in the system
    lattice_const : float
        The lattice constant for an fcc lattice

    Returns
    -------
    pos_vec : np.ndarray
        Array of particle coordinates
    """
    #First create the unit cell
    atom1=[0,0,0]
    atom2=[0.5,0.5,0]
    atom3=[0,0.5,0.5]
    atom4=[0.5,0,0.5]
    unit_cell = lat_const*np.array([atom1, atom2, atom3, atom4])

    #Create displacement arrays
    x_displacement = np.zeros((4,3))
    x_displacement[:,0] = 1
    y_displacement = np.zeros((4,3))
    y_displacement[:,1] = 1
    z_displacement = np.zeros((4,3))
    z_displacement[:,2] = 1

    #Generate the full lattice
    lattice = []

    for i in range(len_x):
        for j in range(len_y):
            for k in range(len_z):
                lattice.append(lat_const*(i*x_displacement+j*y_displacement+k*z_displacement)+unit_cell)

    lattice = np.array(lattice)
    lattice = lattice.reshape(len(lattice[0])*len(lattice),3)

    return lattice


def square_lattice():
    """
    Initializes a system of atoms on an fcc lattice.

    Parameters
    ----------
    side : int
        The number of atoms in one side of the square, total number of atoms will be side**2.
    lattice_const : float
        The lattice constant

    Returns
    -------
    pos_vec : np.ndarray
        Array of particle coordinates
    """
    x,y = np.meshgrid(np.arange(-1,2,1),np.arange(-1,2,1))
    lista = np.dstack([x,y])[:, :].reshape(3**2,2).T.tolist()
    l = len(lista[0])
    lista.append([0 for i in range(l)])

    return np.array(lista).T


def make_2_particles_2d(sigma):
    """
    Creates 2 particles at a particular position of a 2D lattice.

    Parameters:
        box_size: float
            Length of the 2D unit cell (square) side

    Returns:
        atoms: nd.array
            Vectors indicating the position in three dimensions.
    """
    atoms = np.array([[0.3,0.5,0],[0.9,0.5,0]])
    return atoms/sigma

def make_4_particles_2d(sigma):
    """
    Creates 2 particles at a particular position of a 2D lattice.

    Parameters:
    -----------
        box_size: float
            Length of the 2D unit cell (square) side

    Returns:
    --------
        atoms: nd.array
            Vectors indicating the position in three dimensions.
    """
    atoms = np.array([[0.2,0.2,0],[0.2,0.8,0],[0.8,0.2,0],[0.8,0.8,0]])
    return atoms/sigma


def make_n_particles_2d(n, sigma):
    """
    Creates n particles at a ramdom positions of a 2D lattice.

    Parameters:
    -----------
        n: int
            Number of particles
        box_size: float
            Length of the 2D unit cell (square) side

    Returns:
    --------
        atoms: nd.array
            Vectors indicating the position in three dimensions.
    """
    atoms = np.random.rand(n,3)
    atoms[:,2] = 0
    return atoms/sigma

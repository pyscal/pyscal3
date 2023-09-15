import pyscal3.csystem as pc

def calculate_centrosymmetry(system, nmax=12):
    """
    Calculate the centrosymmetry parameter

    Parameters
    ----------
    system: System object

    nmax : int, optional
        number of neighbors to be considered for centrosymmetry 
        parameters. Has to be a positive, even integer. Default 12

    Returns
    -------
    None
    
    Notes
    -----
    Calculate the centrosymmetry parameter for each atom which can be accessed by
    :attr:`~pyscal.catom.centrosymmetry` attribute. It calculates the degree of inversion
    symmetry of an atomic environment. Centrosymmetry recalculates the neighbor using
    the number method as specified in :func:`¬pyscal.core.System.find_neighbors` method. This
    is the ensure that the required number of neighbors are found for calculation of the parameter.

    The Greedy Edge Selection (GES) [1] as specified in [2] in used in this method. 
    GES algorithm is implemented in LAMMPS and Ovito. Please see [2] for
    a detailed description of the algorithms. 
    References
    ----------
    .. [1] Stukowski, A, Model Simul Mater SC 20, 2012
    .. [2] Larsen, arXiv:2003.08879v1, 2020

    """
    if not nmax>0:
        raise ValueError("nmax cannot be negative")

    if not nmax%2 == 0:
        raise ValueError("nmax has to even integer")

    system.atoms.create_attribute('centrosymmetry', fill_with = 0)
    system.find.neighbors(method='number', nmax=nmax, assign_neighbor=True)
    pc.calculate_centrosymmetry(system.atoms, nmax)

    return system.atoms.centrosymmetry

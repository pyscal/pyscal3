import numpy as np
import pyscal3.csystem as pc

def average_over_neighbors(system, key, include_self=True):
    """
    Perform a simple average over neighbor atoms

    Parameters
    ----------
    key: string
        atom property

    include_self: bool, optional
        If True, include the host atom in the calculation

    Returns
    -------

    """
    system._check_neighbors()

    if not key in system.atoms.keys():
        raise KeyError("required property not found!")

    test = system.atoms[key][0]

    if isinstance(test, list):
        raise TypeError("Averaging can only be done over 1D quantities")

    avgarr = []
    for i in range(len(system.atoms["positions"])):
        arr = []
        if include_self:
            arr.append(system.atoms[key][i])
        for j in system.atoms["neighbors"][i]:
            arr.append(system.atoms[key][j])
        avgarr.append(np.mean(arr))
    
    return avgarr 

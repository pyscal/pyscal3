import numpy as np
import pyscal3.csystem as pc


def calculate_entropy(system, rm, sigma=0.2,
    rstart=0.001, h=0.001, local=False, average=False):
    """
    Calculate the entropy parameter for each atom
    
    Parameters
    ----------
    rm : float
        cutoff distance for integration of entropy parameter in distance units

    sigma : float
        broadening parameter

    rstart : float, optional
        minimum limit for integration, default 0.00001

    h : float, optional
        width for trapezoidal integration, default 0.0001

    local : bool, optional
        if True, use the local density instead of global density
        default False

    average : bool, optional
        if True find the averaged entropy parameters
        default False
    Returns
    -------
    None
    """
    kb = 1
    if local:
        rho = 0
    else:
        rho = system.natoms/system.volume

    pc.calculate_entropy(system.atoms, sigma, rho, rstart, rm, h, kb)

    mapdict = {}
    mapdict["entropy"] = {}
    mapdict["entropy"]["norm"] = "entropy"

    if average:
        pc.calculate_average_entropy(system.atoms)
        mapdict["entropy"]["average"] = "average_entropy"
        system.atoms._add_attribute(mapdict)
        return system.atoms.entropy.average
    else:
        system.atoms._add_attribute(mapdict)
        return system.atoms.entropy.norm



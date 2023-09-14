import numpy as np
import pyscal3.csystem as pc


def get_distance(system, pos1, pos2, vector=False):
    """
    Get the distance between two atoms.

    Parameters
    ----------
    pos1 : list
            first atom position
    pos2 : list
            second atom position
    vector: bool, optional
        If True, return the vector between two atoms

    Returns
    -------
    distance : double
            distance between the first and second atom.

    Notes
    -----
    Periodic boundary conditions are assumed by default.
    """
    diff = pc.get_distance_vector(pos1, pos2, system.triclinic,
        system.rot, system.rotinv, system.boxdims)
    dist = np.linalg.norm(diff)
            
    if vector:
        return dist, diff
    else:
        return dist

def reset_neighbors(system):
    """
    Reset the neighbors of all atoms in the system.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    It is used automatically when neighbors are recalculated.

    """
    system.atoms["neighbors"] = []
    system.atoms["neighbordist"] = []
    system.atoms["temp_neighbors"] = []
    system.atoms["temp_neighbordist"] = []
    system.atoms["neighborweight"] = []
    system.atoms["diff"] = []
    system.atoms["r"] = []
    system.atoms["theta"] = []
    system.atoms["phi"] = []
    system.atoms["cutoff"] = []
    system.neighbors_found = False

    
    mapdict = {}
    mapdict["neighbors"] = {}
    mapdict["neighbors"]["index"] = "neighbors"
    mapdict["neighbors"]["distance"] = "neighbordist"
    mapdict["neighbors"]["weight"] = "neighborweight"
    mapdict["neighbors"]["displacement"] = "diff"
    mapdict["neighbors"]["cutoff"] = "cutoff"

    mapdict["neighbors"]["angle"] = {}
    mapdict["neighbors"]["angle"]["polar"] = "theta"
    mapdict["neighbors"]["angle"]["azimuthal"] = "phi"

    mapdict["neighbors"]["temporary"] = {}
    mapdict["neighbors"]["temporary"]["index"] = "temp_neighbors"
    mapdict["neighbors"]["temporary"]["distance"] = "temp_neighbordist"

    system.atoms._add_attribute(mapdict)


	
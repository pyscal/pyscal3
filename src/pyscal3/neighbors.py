"""
Neighbor finding for ASE Atoms using pyscal's C++ routines.

All functions take an ASE Atoms object and store neighbor data
back into it (via atoms.info for ragged arrays).

Example
-------
>>> from ase.build import bulk
>>> import pyscal3
>>> atoms = bulk("Fe", "bcc", cubic=True).repeat(4)
>>> pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
>>> atoms.info["pyscal_neighbors"]  # list of lists
"""
import warnings
import numpy as np
from ase import Atoms

import pyscal3.csystem as pc
from pyscal3._bridge import (
    get_box_params, atoms_to_dict, dict_to_atoms,
    pad_atoms_for_neighbor_finding,
)


def find_neighbors(atoms: Atoms, method='cutoff', cutoff=0, 
                   shell_thickness=0, threshold=2, voroexp=1,
                   padding=1.2, nlimit=6, cells=None, nmax=12,
                   assign_neighbor=True):
    """
    Find neighbors of all atoms.

    Parameters
    ----------
    atoms : ase.Atoms
        The atomic structure.
    method : {'cutoff', 'voronoi', 'number'}
        Neighbor finding algorithm.
    cutoff : float or str
        Cutoff distance. Use 'sann' or 'adaptive' for adaptive methods.
        0 defaults to adaptive.
    shell_thickness : float, optional
        If > 0, find neighbors in a shell [cutoff, cutoff+shell_thickness].
    threshold : float, optional
        Safety multiplier for adaptive methods. Default 2.
    voroexp : int, optional
        Voronoi face area weight exponent. Default 1.
    padding : float, optional
        Safety padding for adaptive/number methods. Default 1.2.
    nlimit : int, optional
        Number of atoms for adaptive cutoff estimation. Default 6.
    cells : bool or None, optional
        Force cell lists on/off. None = auto (>250 atoms).
    nmax : int, optional
        Number of neighbors for 'number' method. Default 12.
    assign_neighbor : bool, optional
        Whether to assign neighbors (for 'number' method). Default True.

    Returns
    -------
    None
        Results are stored in-place on the ``atoms`` object:

        - ``atoms.arrays["pyscal_neighbors"]`` — neighbor indices (natoms, max_neighbors)
        - ``atoms.arrays["pyscal_neighbordist"]`` — neighbor distances
        - ``atoms.arrays["pyscal_theta"]`` — polar angles to neighbors
        - ``atoms.arrays["pyscal_phi"]`` — azimuthal angles to neighbors
        - ``atoms.info["pyscal_neighbors_found"]`` — set to True
        - ``atoms.info["pyscal_neighbor_method"]`` — the method used

        For Voronoi, additional keys include ``face_vertices``, ``face_perimeters``,
        ``face_areas``, ``vertex_vectors``, and ``voronoivol``.
    """
    if threshold < 1:
        raise ValueError("threshold must be >= 1.0")

    # Use ghost padding for small cells
    d, (triclinic, rot, rotinv, boxdims), nreal = pad_atoms_for_neighbor_finding(atoms)
    natoms = len(d["positions"])
    
    if cells is None:
        cells = (natoms > 250)

    # Reset existing neighbor data
    _reset_neighbors(d)
    
    if method == 'cutoff':
        if cutoff == 'sann':
            finished = False
            for i in range(1, 10):
                finished = pc.get_all_neighbors_sann(
                    d, 0.0, triclinic, rot, rotinv, boxdims, threshold * i, cells)
                if finished:
                    if i > 1:
                        warnings.warn("Found neighbors with higher threshold than default/user input")
                    break
                warnings.warn("Could not find sann cutoff. Trying with higher threshold", RuntimeWarning)
            else:
                raise RuntimeError(
                    "SANN cutoff could not be converged. Try increasing threshold.")
        
        elif cutoff == 'adaptive' or (cutoff == 0 and shell_thickness == 0):
            finished = pc.get_all_neighbors_adaptive(
                d, 0.0, triclinic, rot, rotinv, boxdims,
                threshold, nlimit, padding, cells)
            if not bool(finished):
                raise RuntimeError("Could not find adaptive cutoff")
        
        else:
            if cutoff == 0 and shell_thickness > 0:
                cutoff = shell_thickness
                shell_thickness = 0
            if shell_thickness == 0:
                if cells:
                    pc.get_all_neighbors_cells(
                        d, cutoff, triclinic, rot, rotinv, boxdims)
                else:
                    pc.get_all_neighbors_normal(
                        d, cutoff, triclinic, rot, rotinv, boxdims)
            else:
                if cells:
                    pc.get_all_neighbors_shell_cells(
                        d, cutoff, cutoff + shell_thickness,
                        triclinic, rot, rotinv, boxdims)
                else:
                    pc.get_all_neighbors_shell_normal(
                        d, cutoff, cutoff + shell_thickness,
                        triclinic, rot, rotinv, boxdims)
    
    elif method == 'number':
        finished = pc.get_all_neighbors_bynumber(
            d, 0.0, triclinic, rot, rotinv, boxdims,
            threshold, nmax, cells, assign_neighbor)
        if not finished:
            raise RuntimeError("Could not find enough neighbors - try increasing threshold")
    
    elif method == 'voronoi':
        pc.get_all_neighbors_voronoi(
            d, 0.0, triclinic, rot, rotinv, boxdims, voroexp)
        
        if cutoff > 0:
            unique_vertices = pc.clean_voronoi_vertices(
                d, triclinic, rot, rotinv, boxdims, cutoff)
            atoms.info["pyscal_unique_vertices"] = unique_vertices
    
    else:
        raise ValueError(f"Unknown method: {method}. Use 'cutoff', 'voronoi', or 'number'.")
    
    # Write results back to ASE atoms
    dict_to_atoms(d, atoms, nreal=nreal)
    atoms.info["pyscal_neighbors_found"] = True
    atoms.info["pyscal_neighbor_method"] = method


def get_distance(atoms: Atoms, pos1, pos2, vector=False):
    """
    Get the distance between two positions respecting periodic boundaries.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Structure (used for box/PBC info).
    pos1, pos2 : array-like
        Positions.
    vector : bool, optional
        If True, also return the displacement vector.
        
    Returns
    -------
    float or (float, list)
        Distance, and optionally the displacement vector.
    """
    triclinic, rot, rotinv, boxdims = get_box_params(atoms)
    diff = pc.get_distance_vector(list(pos1), list(pos2), triclinic, rot, rotinv, boxdims)
    dist = np.linalg.norm(diff)
    if vector:
        return dist, diff
    return dist


def _reset_neighbors(d: dict):
    """Reset all neighbor data in the dict."""
    n = len(d["positions"])
    d["neighbors"] = [[] for _ in range(n)]
    d["neighbordist"] = [[] for _ in range(n)]
    d["temp_neighbors"] = [[] for _ in range(n)]
    d["temp_neighbordist"] = [[] for _ in range(n)]
    d["neighborweight"] = [[] for _ in range(n)]
    d["diff"] = [[] for _ in range(n)]
    d["r"] = [[] for _ in range(n)]
    d["theta"] = [[] for _ in range(n)]
    d["phi"] = [[] for _ in range(n)]
    d["cutoff"] = [0.0] * n

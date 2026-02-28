"""
pyscal3 — Fast structural descriptors for atomistic simulations.

Built on ASE (Atomic Simulation Environment) with C++ performance.

Quick Start
-----------
>>> from ase.build import bulk
>>> import pyscal3
>>>
>>> atoms = bulk("Cu", "fcc", cubic=True).repeat(4)
>>> pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0.9)
>>> q = pyscal3.steinhardt_parameter(atoms, l=[4, 6])
>>> cna = pyscal3.common_neighbor_analysis(atoms)
>>>
>>> # Create structures
>>> from pyscal3.structures import make_crystal, make_element
>>> fe = make_element("Fe", repetitions=(4, 4, 4))
"""

# Neighbor finding
from pyscal3.neighbors import find_neighbors, get_distance

# Descriptors (the core value)
from pyscal3.descriptors import (
    steinhardt_parameter,
    disorder,
    common_neighbor_analysis,
    diamond_structure,
    centrosymmetry,
    voronoi_vector,
    entropy,
    short_range_order,
    radial_distribution_function,
    angular_criteria,
    chi_params,
    atomic_strain,
    von_mises_strain,
    d2min,
    slip_vector,
    find_solids,
    find_clusters,
    average_over_neighbors,
)

# Trajectory
from pyscal3.trajectory import Trajectory

__version__ = "4.0.0"

__all__ = [
    # Neighbors
    "find_neighbors",
    "get_distance",
    # Descriptors
    "steinhardt_parameter",
    "disorder",
    "common_neighbor_analysis",
    "diamond_structure",
    "centrosymmetry",
    "voronoi_vector",
    "entropy",
    "short_range_order",
    "radial_distribution_function",
    "angular_criteria",
    "chi_params",
    "atomic_strain",
    "von_mises_strain",
    "d2min",
    "slip_vector",
    "find_solids",
    "find_clusters",
    "average_over_neighbors",
    # Trajectory
    "Trajectory",
]

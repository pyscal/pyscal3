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
    get_box_params,
    dict_to_atoms,
    pad_atoms_for_neighbor_finding,
    guess_cutoff,
    clear_neighbor_data,
    effective_periodic_cell,
    gc_paused,
    store_bond_arrays,
    _NONPERIODIC_PAD,
)


def find_neighbors(
    atoms: Atoms,
    method="cutoff",
    cutoff=0,
    shell_thickness=0,
    threshold=2,
    voroexp=1,
    padding=1.2,
    nlimit=6,
    cells=None,
    nmax=12,
    assign_neighbor=True,
    store_rows=True,
):
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
        0 defaults to adaptive. For ``method='voronoi'`` a positive value
        is the distance below which Voronoi vertices are merged into the
        unique interstitial sites stored in
        ``atoms.info["pyscal_unique_vertices"]``.
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
        Ignored. It selected the cell-list search in earlier versions and is
        kept so that existing calls keep working.
    nmax : int, optional
        Number of neighbors for 'number' method. Default 12.
    assign_neighbor : bool, optional
        Whether to assign neighbors (for 'number' method). Default True.
    store_rows : bool, optional
        Whether to store the per-atom row keys (``pyscal_neighbors``,
        ``pyscal_neighbordist``, ...). Default True. With False, only the flat
        ``pyscal_bond_*`` keys and ``pyscal_cutoff`` are stored; all
        descriptors work from these, and building the per-atom Python lists
        for atoms with different numbers of neighbors is skipped.

    Returns
    -------
    None
        Results are stored in-place on the ``atoms`` object with the
        ``pyscal_`` prefix. Per-atom data whose shape is the same for every
        atom is stored in ``atoms.arrays``; ragged data (atoms with
        different numbers of neighbors) and arrays with three or more
        dimensions are stored in ``atoms.info`` under the same key. Any
        neighbor-derived key from a previous search is removed first.

        - ``pyscal_neighbors`` — neighbor indices
        - ``pyscal_neighbordist`` — neighbor distances
        - ``pyscal_neighborweight`` — weights (1, or Voronoi face-area
          fractions)
        - ``pyscal_r``, ``pyscal_theta``, ``pyscal_phi`` — spherical
          coordinates of the neighbor vectors
        - ``pyscal_diff`` — neighbor vectors, shape (natoms, nn, 3), always in
          ``atoms.info``
        - ``pyscal_cutoff`` — per-atom cutoff that was used
        - ``pyscal_neighbors_found`` (True) and ``pyscal_neighbor_method`` in
          ``atoms.info``

        The same bonds are also stored flat in ``atoms.info``:
        ``pyscal_bond_offsets`` (natoms + 1 entries; the bonds of atom ``i``
        are ``offsets[i]:offsets[i + 1]``), ``pyscal_bond_neighbors``,
        ``pyscal_bond_distance``, ``pyscal_bond_weight``,
        ``pyscal_bond_vector`` (nbonds, 3), ``pyscal_bond_theta`` and
        ``pyscal_bond_phi``. The adaptive, SANN and number methods also store
        their candidates as ``pyscal_candidate_offsets``,
        ``pyscal_candidate_neighbors`` and ``pyscal_candidate_distance``.
        With ``store_rows=False`` only these flat keys and ``pyscal_cutoff``
        are stored.

        The Voronoi method additionally stores ``pyscal_voronoi_volume``,
        ``pyscal_face_vertices``, ``pyscal_face_perimeters``,
        ``pyscal_vertex_vectors``, ``pyscal_vertex_numbers`` and
        ``pyscal_vertex_positions``.

    Notes
    -----
    Except for Voronoi, the search uses the bundled matscipy-neighbours
    library. A fixed cutoff keeps pairs with ``d < cutoff``, a shell keeps
    ``cutoff <= d <= cutoff + shell_thickness``. The adaptive, SANN and number
    methods start from the candidates with ``d <= threshold * (V / N)**(1/3)``,
    sorted by distance and, for distances equal to within 1e-10, by atom
    index. ``pyscal_diff`` holds ``r_i - r_j``. When the cutoff exceeds half
    the cell width, each periodic image of a neighbor is a separate entry.
    """
    if threshold < 1:
        raise ValueError("threshold must be >= 1.0")

    # drop everything derived from a previous neighbor search
    clear_neighbor_data(atoms)

    if method == "cutoff":
        if cutoff == "sann":
            for i in range(1, 10):
                res = pc.nl_sann(*_geometry(atoms, threshold * i), threshold * i)
                if res["finished"]:
                    if i > 1:
                        warnings.warn(
                            "Found neighbors with higher threshold than default/user input"
                        )
                    break
                warnings.warn(
                    "Could not find sann cutoff. Trying with higher threshold",
                    RuntimeWarning,
                )
            else:
                raise RuntimeError(
                    "SANN cutoff could not be converged. Try increasing threshold."
                )

        elif cutoff == "adaptive" or (cutoff == 0 and shell_thickness == 0):
            res = pc.nl_adaptive(
                *_geometry(atoms, threshold), threshold, nlimit, padding
            )
            if not res["finished"]:
                raise RuntimeError("Could not find adaptive cutoff")

        else:
            if cutoff == 0 and shell_thickness > 0:
                cutoff = shell_thickness
                shell_thickness = 0
            pad = cutoff + shell_thickness
            positions, cell, pbc = _working_cell(atoms, pad if pad > 0 else None)
            if shell_thickness == 0:
                res = pc.nl_cutoff(positions, cell, pbc, cutoff)
            else:
                res = pc.nl_shell(positions, cell, pbc, cutoff, cutoff + shell_thickness)

    elif method == "number":
        res = pc.nl_number(
            *_geometry(atoms, threshold), threshold, nmax, bool(assign_neighbor)
        )
        if not res["finished"]:
            raise RuntimeError(
                "Could not find enough neighbors - try increasing threshold"
            )

    elif method == "voronoi":
        d, (triclinic, rot, rotinv, boxdims), nreal = pad_atoms_for_neighbor_finding(atoms)
        _reset_neighbors(d)
        pc.get_all_neighbors_voronoi(d, 0.0, triclinic, rot, rotinv, boxdims, voroexp)
        dict_to_atoms(d, atoms, nreal=nreal)
        store_bond_arrays(atoms)
        if not store_rows:
            _drop_rows(atoms)
        if isinstance(cutoff, (int, float)) and cutoff > 0:
            # merge Voronoi vertices closer than `cutoff` into unique sites
            atoms.info["pyscal_unique_vertices"] = _unique_voronoi_vertices(
                atoms, d["vertex_positions"][:nreal], cutoff
            )

    else:
        raise ValueError(
            f"Unknown method: {method}. Use 'cutoff', 'voronoi', or 'number'."
        )

    if method != "voronoi":
        _store_neighbors(atoms, res, store_rows)
    atoms.info["pyscal_neighbors_found"] = True
    atoms.info["pyscal_neighbor_method"] = method


def _working_cell(atoms: Atoms, pad):
    """Positions, cell and periodicity passed to the C++ search.

    Non-periodic directions and zero cell vectors get the vacuum cell vector
    of :func:`effective_periodic_cell` (which also validates the cell) and
    stay non-periodic. ``pad`` is the search radius, None if unknown.
    """
    if len(atoms) == 0:
        raise ValueError("Cannot find neighbors of an empty Atoms object.")
    cell, periodic = effective_periodic_cell(
        atoms, pad if pad is not None else _NONPERIODIC_PAD
    )
    positions = np.ascontiguousarray(atoms.positions, dtype=float)
    return positions, np.ascontiguousarray(cell, dtype=float), [bool(p) for p in periodic]


def _geometry(atoms: Atoms, prefactor):
    """Working cell for the candidate-based methods (adaptive, SANN, number).

    The C++ code takes the candidate radius as prefactor * (V / N)^(1/3) of
    this cell, which includes the vacuum added along non-periodic directions.
    """
    pad = guess_cutoff(atoms, prefactor) if len(atoms) else 0.0
    return _working_cell(atoms, pad if pad > 0 else None)


_PER_PAIR_KEYS = {
    "neighbors": "j",
    "neighbordist": "d",
    "neighborweight": "weight",
    "r": "r",
    "theta": "theta",
    "phi": "phi",
}


def _store_rows(atoms: Atoms, rows: dict, offsets, vectors=None):
    """Store per-pair values, one row per atom, under ``pyscal_<key>``.

    Rows of equal length go to ``atoms.arrays`` as (n, k) arrays, with
    ``vectors`` as an (n, k, 3) array in ``atoms.info``. Rows of different
    length go to ``atoms.info`` as lists of lists. Empty rows for every atom
    are stored as (n, 0) float arrays in ``atoms.arrays``.
    """
    n = len(atoms)
    counts = np.diff(offsets)
    k = int(counts[0]) if n else 0
    if (counts == k).all():
        for key, values in rows.items():
            if k == 0:
                atoms.arrays["pyscal_" + key] = np.zeros((n, 0))
            else:
                atoms.arrays["pyscal_" + key] = values.reshape(n, k)
        if vectors is not None:
            if k == 0:
                atoms.arrays["pyscal_diff"] = np.zeros((n, 0))
            else:
                atoms.info["pyscal_diff"] = vectors.reshape(n, k, 3)
        return
    with gc_paused():
        bounds = list(zip(offsets[:-1].tolist(), offsets[1:].tolist()))
        for key, values in rows.items():
            flat = values.tolist()
            atoms.info["pyscal_" + key] = [flat[a:b] for a, b in bounds]
        if vectors is not None:
            flat = vectors.tolist()
            atoms.info["pyscal_diff"] = [flat[a:b] for a, b in bounds]


_ROW_KEYS = ("neighbors", "neighbordist", "neighborweight", "diff", "r", "theta",
             "phi", "temp_neighbors", "temp_neighbordist")


def _drop_rows(atoms: Atoms):
    """Remove the per-atom neighbor row keys (the flat keys stay)."""
    for key in _ROW_KEYS:
        atoms.arrays.pop("pyscal_" + key, None)
        atoms.info.pop("pyscal_" + key, None)


def _store_neighbors(atoms: Atoms, res: dict, store_rows=True):
    """Write the result of a pc.nl_* search to ``atoms``."""
    info = atoms.info
    info["pyscal_bond_offsets"] = res["offsets"]
    info["pyscal_bond_neighbors"] = res["j"]
    info["pyscal_bond_distance"] = res["d"]
    info["pyscal_bond_weight"] = res["weight"]
    info["pyscal_bond_vector"] = res["diff"]
    info["pyscal_bond_theta"] = res["theta"]
    info["pyscal_bond_phi"] = res["phi"]
    if "temp_offsets" in res:
        info["pyscal_candidate_offsets"] = res["temp_offsets"]
        info["pyscal_candidate_neighbors"] = res["temp_j"]
        info["pyscal_candidate_distance"] = res["temp_d"]
    atoms.arrays["pyscal_cutoff"] = res["cutoff"]
    if not store_rows:
        return
    _store_rows(
        atoms,
        {key: res[src] for key, src in _PER_PAIR_KEYS.items()},
        res["offsets"],
        vectors=res["diff"],
    )
    if "temp_offsets" in res:
        _store_rows(
            atoms,
            {"temp_neighbors": res["temp_j"], "temp_neighbordist": res["temp_d"]},
            res["temp_offsets"],
        )
    else:
        _store_rows(
            atoms,
            {"temp_neighbors": np.zeros(0), "temp_neighbordist": np.zeros(0)},
            np.zeros(len(atoms) + 1, dtype=np.int64),
        )


def _unique_voronoi_vertices(atoms: Atoms, vertex_positions, cutoff):
    """Merge the Voronoi vertices of all atoms into unique sites.

    Vertices closer than ``cutoff`` (under the periodic boundary conditions
    of ``atoms``) are considered the same site; one representative per
    group is returned, wrapped into the cell. A Voronoi vertex is shared by
    every cell meeting there, which are not necessarily Voronoi neighbours
    of each other, so the merge has to be done over all vertices.
    """
    from ase.neighborlist import neighbor_list

    pts = [np.asarray(v, dtype=float).reshape(-1, 3) for v in vertex_positions if len(v)]
    if not pts:
        return np.zeros((0, 3))
    pts = np.concatenate(pts)

    work_cell, _ = effective_periodic_cell(atoms, cutoff)
    dummy = Atoms(positions=pts, cell=work_cell, pbc=True)
    dummy.wrap()
    i, j = neighbor_list("ij", dummy, cutoff)

    # union-find over the close pairs
    parent = np.arange(len(pts))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b in zip(i, j):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    roots = np.array([find(a) for a in range(len(pts))])
    return dummy.positions[np.unique(roots)]


def get_distance(atoms: Atoms, pos1, pos2, vector=False):
    """
    Get the distance between two positions respecting periodic boundaries
    (non-periodic directions of ``atoms`` are not wrapped).

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
    plain = float(np.linalg.norm(np.asarray(pos2, float) - np.asarray(pos1, float)))
    work_cell, periodic = effective_periodic_cell(atoms, max(plain, 1.0))
    if periodic.all():
        triclinic, rot, rotinv, boxdims = get_box_params(atoms)
    else:
        work = Atoms(cell=work_cell, pbc=True)
        triclinic, rot, rotinv, boxdims = get_box_params(work)
    diff = pc.get_distance_vector(
        list(pos1), list(pos2), triclinic, rot, rotinv, boxdims
    )
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

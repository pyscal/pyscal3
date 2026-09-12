"""Tests for trajectory module."""

from pathlib import Path
import numpy as np
from pyscal3.trajectory import Trajectory

ROOT = Path(__file__).resolve().parent.parent
TRAJ = str(ROOT / "examples" / "traj.light")


def test_traj_basics():
    traj = Trajectory(TRAJ)
    assert traj.nblocks == 10

    block = traj.get_block(0)
    assert len(block) == 509


def test_traj_load_unload():
    traj = Trajectory(TRAJ)
    traj.load(0)
    data = traj.data[0]
    assert data["box"][0][0] == -7.34762
    assert data["atoms"]["x"][0] == -4.72745
    traj.unload(0)
    assert traj.data[0] is None


def test_timeslice_to_atoms():
    traj = Trajectory(TRAJ)
    atoms_list = traj[0].to_atoms(species=["Au"])
    assert len(atoms_list) == 1
    # Check box length
    assert abs(atoms_list[0].cell[0][0] - 18.21922) < 0.001
    # Check first position
    assert abs(atoms_list[0].positions[0][0] - (-4.72745)) < 0.001


def test_timeslice_to_dict():
    traj = Trajectory(TRAJ)
    od = traj[0].to_dict()
    assert od[0]["box"][0][0] == -7.34762


def test_timeslice_to_file():
    import os

    traj = Trajectory(TRAJ)
    traj[0].to_file("test_traj_output.dump")
    assert os.path.exists("test_traj_output.dump")
    os.remove("test_traj_output.dump")


def test_timeslice_slice():
    traj = Trajectory(TRAJ)
    sl = traj[0:3]
    atoms_list = sl.to_atoms(species=["Au"])
    assert len(atoms_list) == 3


def _triclinic_dump(scaled, frac, cell, origin):
    """Write a one-frame LAMMPS dump for a triclinic box (bounding-box format)."""
    (ax, _, _), (xy, by, _), (xz, yz, cz) = cell
    xlo, ylo, zlo = origin
    xhi, yhi, zhi = xlo + ax, ylo + by, zlo + cz
    xlo_b = xlo + min(0.0, xy, xz, xy + xz)
    xhi_b = xhi + max(0.0, xy, xz, xy + xz)
    ylo_b = ylo + min(0.0, yz)
    yhi_b = yhi + max(0.0, yz)
    lines = ["ITEM: TIMESTEP", "0", "ITEM: NUMBER OF ATOMS", str(len(frac)),
             "ITEM: BOX BOUNDS xy xz yz pp pp pp",
             f"{xlo_b} {xhi_b} {xy}", f"{ylo_b} {yhi_b} {xz}", f"{zlo_b if False else zlo} {zhi} {yz}"]
    if scaled:
        lines.append("ITEM: ATOMS id type xs ys zs")
        coords = frac
    else:
        lines.append("ITEM: ATOMS id type x y z")
        coords = frac @ cell + origin
    for i, c in enumerate(coords):
        lines.append(f"{i + 1} 1 {c[0]!r} {c[1]!r} {c[2]!r}")
    return [l + "\n" for l in lines]


def test_parse_triclinic_scaled_and_unscaled_agree():
    from pyscal3.trajectory import _parse_lammps_lines_to_atoms

    rng = np.random.default_rng(0)
    cell = np.array([[10.0, 0.0, 0.0], [2.0, 9.0, 0.0], [1.5, -1.0, 8.0]])
    origin = np.array([-1.0, 0.5, 2.0])
    frac = rng.uniform(0, 1, size=(20, 3))
    expected = frac @ cell + origin

    for scaled in (False, True):
        atoms = _parse_lammps_lines_to_atoms(_triclinic_dump(scaled, frac, cell, origin), species=["Cu"])
        assert np.allclose(np.array(atoms.cell), cell)
        assert np.allclose(atoms.positions, expected, atol=1e-8), f"scaled={scaled}"
        assert np.allclose(atoms.get_celldisp().ravel(), origin)
        assert atoms.get_chemical_symbols() == ["Cu"] * 20


def test_parse_orthogonal_scaled():
    from pyscal3.trajectory import _parse_lammps_lines_to_atoms

    rng = np.random.default_rng(1)
    cell = np.diag([10.0, 9.0, 8.0])
    origin = np.array([-1.0, 0.5, 2.0])
    frac = rng.uniform(0, 1, size=(10, 3))
    a = _parse_lammps_lines_to_atoms(_triclinic_dump(False, frac, cell, origin), species=["Al"])
    b = _parse_lammps_lines_to_atoms(_triclinic_dump(True, frac, cell, origin), species=["Al"])
    assert np.allclose(a.positions, b.positions, atol=1e-8)
    assert np.allclose(a.positions, frac @ cell + origin, atol=1e-8)

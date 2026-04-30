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

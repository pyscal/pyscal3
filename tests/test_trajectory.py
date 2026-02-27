"""Tests for trajectory module."""
import numpy as np
from pyscal3.trajectory import Trajectory


def test_traj_basics():
    traj = Trajectory("examples/traj.light")
    assert traj.nblocks == 10

    block = traj.get_block(0)
    assert len(block) == 509


def test_traj_load_unload():
    traj = Trajectory("examples/traj.light")
    traj.load(0)
    data = traj.data[0]
    assert data["box"][0][0] == -7.34762
    assert data["atoms"]["x"][0] == -4.72745
    traj.unload(0)
    assert traj.data[0] is None


def test_timeslice_to_atoms():
    traj = Trajectory("examples/traj.light")
    atoms_list = traj[0].to_atoms(species=["Au"])
    assert len(atoms_list) == 1
    # Check box length
    assert abs(atoms_list[0].cell[0][0] - 18.21922) < 0.001
    # Check first position
    assert abs(atoms_list[0].positions[0][0] - (-4.72745)) < 0.001


def test_timeslice_to_dict():
    traj = Trajectory("examples/traj.light")
    od = traj[0].to_dict()
    assert od[0]["box"][0][0] == -7.34762


def test_timeslice_to_file():
    import os
    traj = Trajectory("examples/traj.light")
    traj[0].to_file("test_traj_output.dump")
    assert os.path.exists("test_traj_output.dump")
    os.remove("test_traj_output.dump")


def test_timeslice_slice():
    traj = Trajectory("examples/traj.light")
    sl = traj[0:3]
    atoms_list = sl.to_atoms(species=["Au"])
    assert len(atoms_list) == 3

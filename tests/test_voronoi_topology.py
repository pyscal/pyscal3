"""Tests for Voronoi topology descriptors."""

import numpy as np
import pytest
from ase.build import bulk

import pyscal3 as pc


def _voronoi_atoms(structure, cubic=True, repeat=3):
    """Build atoms with Voronoi tessellation."""
    atoms = bulk(structure, cubic=cubic)
    atoms = atoms.repeat(repeat)
    pc.find_neighbors(atoms, method="voronoi")
    return atoms


# ---------------------------------------------------------------------------
# Cell volume
# ---------------------------------------------------------------------------

class TestVoronoiCellVolume:

    def test_fcc_volume(self):
        """FCC Al: each Voronoi cell = a^3/4 ~ 16.6 A^3."""
        atoms = _voronoi_atoms("Al")
        vol = pc.voronoi_cell_volume(atoms)
        expected = atoms.cell.volume / len(atoms)
        assert np.allclose(vol, expected, rtol=0.01)

    def test_bcc_volume(self):
        """BCC Fe: each Voronoi cell = a^3/2 ~ 11.8 A^3."""
        atoms = _voronoi_atoms("Fe")
        vol = pc.voronoi_cell_volume(atoms)
        expected = atoms.cell.volume / len(atoms)
        assert np.allclose(vol, expected, rtol=0.01)

    def test_stored_in_arrays(self):
        atoms = _voronoi_atoms("Al")
        pc.voronoi_cell_volume(atoms)
        assert "pyscal_voronoi_volume" in atoms.arrays


# ---------------------------------------------------------------------------
# Sphericity
# ---------------------------------------------------------------------------

class TestVoronoiSphericity:

    def test_fcc_sphericity(self):
        """FCC rhombic dodecahedron IQ should be > 0.7."""
        atoms = _voronoi_atoms("Al")
        iq = pc.voronoi_sphericity(atoms)
        # Rhombic dodecahedron IQ ~ 0.740
        assert np.all(iq > 0.7), f"FCC IQ too low: {iq[0]}"
        assert np.all(iq < 1.0), "IQ should be < 1"

    def test_bcc_sphericity(self):
        """BCC truncated octahedron IQ should be > 0.7."""
        atoms = _voronoi_atoms("Fe")
        iq = pc.voronoi_sphericity(atoms)
        # Truncated octahedron IQ ~ 0.910
        assert np.all(iq > 0.7), f"BCC IQ too low: {iq[0]}"

    def test_stored_in_arrays(self):
        atoms = _voronoi_atoms("Al")
        pc.voronoi_sphericity(atoms)
        assert "pyscal_voronoi_sphericity" in atoms.arrays

    def test_uniform_for_perfect_crystal(self):
        """All atoms in a perfect crystal should have the same IQ."""
        atoms = _voronoi_atoms("Al")
        iq = pc.voronoi_sphericity(atoms)
        assert np.allclose(iq, iq[0], rtol=1e-6)


# ---------------------------------------------------------------------------
# Face analysis
# ---------------------------------------------------------------------------

class TestVoronoiFaceAnalysis:

    def test_fcc_nfaces(self):
        """FCC Voronoi cell has 12 faces (rhombic dodecahedron)."""
        atoms = _voronoi_atoms("Al")
        result = pc.voronoi_face_analysis(atoms)
        assert np.all(result["n_faces"] == 12)

    def test_bcc_nfaces(self):
        """BCC Voronoi cell has 14 faces (truncated octahedron)."""
        atoms = _voronoi_atoms("Fe")
        result = pc.voronoi_face_analysis(atoms)
        assert np.all(result["n_faces"] == 14)

    def test_face_area_positive(self):
        atoms = _voronoi_atoms("Al")
        result = pc.voronoi_face_analysis(atoms)
        assert np.all(result["mean_face_area"] > 0)

    def test_stored_in_arrays(self):
        atoms = _voronoi_atoms("Al")
        pc.voronoi_face_analysis(atoms)
        assert "pyscal_voronoi_nfaces" in atoms.arrays
        assert "pyscal_voronoi_mean_face_area" in atoms.arrays
        assert "pyscal_voronoi_std_face_area" in atoms.arrays

    def test_all_functions_together(self):
        """All Voronoi topology functions should work on the same atoms."""
        atoms = _voronoi_atoms("Cu")
        vol = pc.voronoi_cell_volume(atoms)
        iq = pc.voronoi_sphericity(atoms)
        fa = pc.voronoi_face_analysis(atoms)
        assert vol.shape == (len(atoms),)
        assert iq.shape == (len(atoms),)
        assert fa["n_faces"].shape == (len(atoms),)


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------

class TestVoronoiErrors:

    def test_no_voronoi_raises(self):
        """Should raise if Voronoi tessellation not computed."""
        atoms = bulk("Al", cubic=True).repeat(3)
        pc.find_neighbors(atoms, method="cutoff", cutoff=0)
        with pytest.raises(ValueError, match="Voronoi"):
            pc.voronoi_cell_volume(atoms)

    def test_no_voronoi_sphericity_raises(self):
        atoms = bulk("Al", cubic=True).repeat(3)
        pc.find_neighbors(atoms, method="cutoff", cutoff=0)
        with pytest.raises(ValueError, match="Voronoi"):
            pc.voronoi_sphericity(atoms)

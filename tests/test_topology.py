"""Tests for topological/graph-based descriptors."""

import numpy as np
import pytest
from ase.build import bulk, molecule
from ase import Atoms

import pyscal3


class TestCoordinationNumbers:
    """Tests for coordination_numbers function."""

    def test_fcc_coordination(self):
        """FCC should have 12 nearest neighbors."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        cn = pyscal3.coordination_numbers(atoms)
        
        # Interior atoms should have 12 neighbors
        assert cn.max() == 12
        assert len(cn) == len(atoms)

    def test_bcc_coordination(self):
        """BCC should have 8 nearest neighbors."""
        atoms = bulk("Fe", "bcc", cubic=True, a=2.87).repeat(3)
        # Use smaller cutoff to get only 1st shell (sqrt(3)/2 * a ≈ 2.49 Å)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=2.6)
        cn = pyscal3.coordination_numbers(atoms)
        
        # Should be 8 for 1st shell
        assert cn.max() == 8

    def test_diamond_coordination(self):
        """Diamond structure should have 4 nearest neighbors."""
        atoms = bulk("Si", "diamond", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=2.8)
        cn = pyscal3.coordination_numbers(atoms)
        
        # Diamond has tetrahedral coordination = 4
        assert np.all(cn == 4)

    def test_stored_in_arrays(self):
        """Coordination should be stored in atoms.arrays."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        cn = pyscal3.coordination_numbers(atoms)
        
        assert "pyscal_coordination" in atoms.arrays
        np.testing.assert_array_equal(atoms.arrays["pyscal_coordination"], cn)


class TestCoordinationStats:
    """Tests for coordination_stats function."""

    def test_perfect_crystal_stats(self):
        """Perfect crystal should have low variance in coordination."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        stats = pyscal3.coordination_stats(atoms)
        
        assert 'coordination' in stats
        assert 'mean' in stats
        assert 'std' in stats
        assert 'distribution' in stats
        
        # Mean should be close to 12 for FCC interior
        assert stats['mean'] > 10
        
    def test_distribution_sums_to_one(self):
        """Distribution fractions should sum to 1."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        stats = pyscal3.coordination_stats(atoms)
        
        total = sum(stats['distribution'].values())
        assert abs(total - 1.0) < 1e-10


class TestBondAngleDistribution:
    """Tests for bond_angle_distribution function."""

    def test_fcc_bond_angles(self):
        """FCC should have characteristic angles at 60, 90, 120, 180 degrees."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        bad = pyscal3.bond_angle_distribution(atoms, bins=90)
        
        assert 'angles' in bad
        assert 'bin_centers' in bad
        assert 'histogram' in bad
        assert 'mean' in bad
        assert 'std' in bad
        
        assert len(bad['angles']) > 0
        assert len(bad['bin_centers']) == 90
        assert len(bad['histogram']) == 90

    def test_diamond_tetrahedral_angle(self):
        """Diamond structure should have dominant angle near 109.47 degrees."""
        atoms = bulk("Si", "diamond", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=2.8)
        bad = pyscal3.bond_angle_distribution(atoms, bins=180)
        
        # Find peak in histogram
        peak_idx = np.argmax(bad['histogram'])
        peak_angle = bad['bin_centers'][peak_idx]
        
        # Tetrahedral angle is 109.47 degrees
        assert abs(peak_angle - 109.47) < 5.0

    def test_custom_bins_and_range(self):
        """Test custom binning parameters."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        bad = pyscal3.bond_angle_distribution(atoms, bins=36, range_deg=(30, 150))
        
        assert len(bad['bin_centers']) == 36
        assert bad['bin_centers'].min() > 30
        assert bad['bin_centers'].max() < 150

    def test_stored_in_info(self):
        """Statistics should be stored in atoms.info."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        bad = pyscal3.bond_angle_distribution(atoms)
        
        assert "pyscal_bond_angle_distribution" in atoms.info
        assert 'mean' in atoms.info["pyscal_bond_angle_distribution"]


class TestClusteringCoefficient:
    """Tests for clustering_coefficient function."""

    def test_clustering_shape(self):
        """Clustering coefficient should be computed for all atoms."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        cc = pyscal3.clustering_coefficient(atoms)
        
        assert 'local' in cc
        assert 'global' in cc
        assert 'transitivity' in cc
        
        assert len(cc['local']) == len(atoms)

    def test_clustering_bounds(self):
        """Clustering coefficients should be between 0 and 1."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        cc = pyscal3.clustering_coefficient(atoms)
        
        assert np.all(cc['local'] >= 0.0)
        assert np.all(cc['local'] <= 1.0)
        assert 0.0 <= cc['global'] <= 1.0
        assert 0.0 <= cc['transitivity'] <= 1.0

    def test_crystal_has_nonzero_clustering(self):
        """Crystal structures should have non-zero clustering."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        cc = pyscal3.clustering_coefficient(atoms)
        
        assert cc['global'] > 0.2  # FCC has significant clustering

    def test_stored_in_arrays(self):
        """Clustering should be stored in atoms.arrays."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        cc = pyscal3.clustering_coefficient(atoms)
        
        assert "pyscal_clustering_coefficient" in atoms.arrays
        np.testing.assert_array_equal(
            atoms.arrays["pyscal_clustering_coefficient"], cc['local']
        )


class TestRingStatistics:
    """Tests for ring_statistics function."""

    def test_ring_counts_keys(self):
        """Ring counts should have keys from 3 to max_ring_size."""
        atoms = bulk("Si", "diamond", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=2.8)
        rings = pyscal3.ring_statistics(atoms, max_ring_size=8)
        
        assert 'counts' in rings
        assert 'per_atom' in rings
        assert 'mean_size' in rings
        assert 'rings' in rings
        
        for n in range(3, 9):
            assert n in rings['counts']

    def test_diamond_has_6_rings(self):
        """Diamond structure should have 6-membered rings."""
        atoms = bulk("Si", "diamond", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=2.8)
        rings = pyscal3.ring_statistics(atoms, max_ring_size=8)
        
        # Diamond structure has 6-membered chair rings
        assert rings['counts'][6] > 0

    def test_max_ring_size(self):
        """Should respect max_ring_size parameter."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        
        rings6 = pyscal3.ring_statistics(atoms, max_ring_size=6)
        rings8 = pyscal3.ring_statistics(atoms, max_ring_size=8)
        
        assert max(rings6['counts'].keys()) == 6
        assert max(rings8['counts'].keys()) == 8

    def test_stored_in_info(self):
        """Ring stats should be stored in atoms.info."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        rings = pyscal3.ring_statistics(atoms)
        
        assert "pyscal_ring_statistics" in atoms.info


class TestTopologicalDescriptors:
    """Tests for combined topological_descriptors function."""

    def test_all_descriptors_computed(self):
        """Should compute all descriptor types."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        topo = pyscal3.topological_descriptors(atoms)
        
        assert 'coordination' in topo
        assert 'clustering' in topo
        assert 'bond_angles' in topo
        # rings not computed by default
        assert 'rings' not in topo

    def test_with_rings(self):
        """Should compute rings when requested."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        topo = pyscal3.topological_descriptors(atoms, compute_rings=True)
        
        assert 'rings' in topo
        assert 'counts' in topo['rings']

    def test_coordination_stats_valid(self):
        """Coordination stats should be valid."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        topo = pyscal3.topological_descriptors(atoms)
        
        assert isinstance(topo['coordination']['mean'], float)
        assert isinstance(topo['coordination']['std'], float)

    def test_clustering_valid(self):
        """Clustering stats should be valid."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        topo = pyscal3.topological_descriptors(atoms)
        
        assert isinstance(topo['clustering']['global'], float)
        assert 0.0 <= topo['clustering']['global'] <= 1.0

    def test_bond_angles_valid(self):
        """Bond angle stats should be valid."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        topo = pyscal3.topological_descriptors(atoms)
        
        assert isinstance(topo['bond_angles']['mean'], float)
        assert 0.0 <= topo['bond_angles']['mean'] <= 180.0


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_single_atom(self):
        """Handle single atom case - single atom has no neighbors."""
        atoms = Atoms('Cu', positions=[[0, 0, 0]], cell=[10, 10, 10], pbc=True)
        # Cannot use cutoff method for single atom (will find self)
        # Use voronoi which handles single atom edge case
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=2.0)
        
        cn = pyscal3.coordination_numbers(atoms)
        assert cn[0] == 0
        
        cc = pyscal3.clustering_coefficient(atoms)
        assert cc['local'][0] == 0.0
        assert cc['global'] == 0.0

    def test_no_neighbors(self):
        """Handle atoms with no neighbors."""
        atoms = Atoms('Cu2', positions=[[0, 0, 0], [10, 10, 10]], cell=[20, 20, 20], pbc=True)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=2.0)
        
        cn = pyscal3.coordination_numbers(atoms)
        assert np.all(cn == 0)
        
        bad = pyscal3.bond_angle_distribution(atoms)
        assert len(bad['angles']) == 0

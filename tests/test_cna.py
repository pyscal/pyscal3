"""Tests for Common Neighbor Analysis and Diamond Structure."""
import numpy as np
import pyscal3
from pyscal3.structures import make_crystal
from ase.build import bulk


def test_cna_fcc_with_lattice_constant():
    atoms = make_crystal("fcc", lattice_constant=4.00, repetitions=(7, 7, 7))
    res = pyscal3.common_neighbor_analysis(atoms, lattice_constant=4.00)
    assert res["fcc"] == 7 * 7 * 7 * 4


def test_cna_fcc_adaptive():
    atoms = make_crystal("fcc", lattice_constant=4.00, repetitions=(7, 7, 7))
    res = pyscal3.common_neighbor_analysis(atoms)
    assert res["fcc"] == 7 * 7 * 7 * 4


def test_cna_bcc_adaptive():
    atoms = make_crystal("bcc", lattice_constant=4.00, repetitions=(7, 7, 7))
    res = pyscal3.common_neighbor_analysis(atoms)
    assert atoms.arrays["pyscal_structure"][0] == 3


def test_cna_hcp_adaptive():
    atoms = make_crystal("hcp", lattice_constant=4.00, repetitions=(7, 7, 7))
    res = pyscal3.common_neighbor_analysis(atoms)
    assert atoms.arrays["pyscal_structure"][0] == 2


def test_cna_ase_bulks():
    al = bulk("Al")
    res = pyscal3.common_neighbor_analysis(al)
    assert res["fcc"] == 1

    fe = bulk("Fe")
    res = pyscal3.common_neighbor_analysis(fe)
    assert res["bcc"] == 1

    ti = bulk("Ti")
    res = pyscal3.common_neighbor_analysis(ti)
    assert res["hcp"] == 2


def test_diamond_structure():
    atoms = make_crystal("diamond", lattice_constant=4.00, repetitions=(7, 7, 7))
    res = pyscal3.diamond_structure(atoms)
    assert atoms.arrays["pyscal_structure"][0] == 1


def test_cna_small_cluster_reports_others_with_warning():
    """A cluster too small to give 14 candidates is labelled, not rejected."""
    import pytest

    atoms = bulk("Al", cubic=True)
    atoms.set_pbc(False)
    with pytest.warns(RuntimeWarning, match="atoms.pbc"):
        res = pyscal3.common_neighbor_analysis(atoms)
    assert res["others"] == len(atoms)
    assert sum(res.values()) == len(atoms)


def test_diamond_small_cluster_reports_others_with_warning():
    import pytest

    atoms = bulk("Si", cubic=True)[:2]
    atoms.set_pbc(False)
    with pytest.warns(RuntimeWarning):
        res = pyscal3.diamond_structure(atoms)
    assert res["others"] == len(atoms)


def test_cna_periodic_cell_unaffected_by_small_size():
    """The same two-atom bcc cell is classified when it is periodic."""
    from ase import Atoms as ASEAtoms

    atoms = ASEAtoms(
        "FeFe", scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=np.identity(3),
        pbc=True,
    )
    assert pyscal3.common_neighbor_analysis(atoms)["bcc"] == 2


def test_cna_nanoparticle_interior_is_fcc():
    """Surface atoms fall into 'others', the interior is still recognised."""
    from ase.cluster import FaceCenteredCubic

    atoms = FaceCenteredCubic("Al", [(1, 0, 0), (1, 1, 1)], [4, 4], 4.05)
    res = pyscal3.common_neighbor_analysis(atoms)
    assert res["fcc"] > 0
    assert res["others"] > 0
    assert sum(res.values()) == len(atoms)


def test_diamond_undercoordinated_atoms_do_not_crash():
    """Fewer than four candidates used to read past the neighbour lists."""
    atoms = bulk("Si", cubic=True)[:2]          # periodic, but only two atoms
    res = pyscal3.diamond_structure(atoms)
    assert res["others"] == 2
    assert sum(res.values()) == 2


def test_cna_and_diamond_on_a_single_atom():
    import pytest
    from ase import Atoms as ASEAtoms

    atoms = ASEAtoms("H", positions=[[0.0, 0.0, 0.0]], cell=np.eye(3) * 20, pbc=False)
    with pytest.warns(RuntimeWarning):
        assert pyscal3.common_neighbor_analysis(atoms.copy())["others"] == 1
    with pytest.warns(RuntimeWarning):
        assert pyscal3.diamond_structure(atoms.copy())["others"] == 1


def _slab_with_far_atom(atoms, vacuum=15.0):
    """Turn a periodic crystal into a z slab, and a copy with one isolated atom."""
    from ase import Atoms as ASEAtoms

    slab = atoms.copy()
    cell = np.array(slab.cell)
    cell[2, 2] += 2 * vacuum
    slab.set_cell(cell)
    slab.center(axis=2)
    lonely = slab.copy()
    lonely += ASEAtoms(slab.get_chemical_symbols()[0],
                       positions=[[1.0, 1.0, cell[2, 2] - 2.0]])
    return slab, lonely


def test_cna_one_isolated_atom_does_not_relabel_the_rest():
    """Only the under-coordinated atom is skipped, not every atom."""
    import pytest

    slab, lonely = _slab_with_far_atom(bulk("Cu", "fcc", a=3.61, cubic=True).repeat(6))
    ref = pyscal3.common_neighbor_analysis(slab)
    assert ref["fcc"] > 0
    with pytest.warns(RuntimeWarning):
        res = pyscal3.common_neighbor_analysis(lonely)
    assert res["fcc"] == ref["fcc"]
    assert res["others"] == ref["others"] + 1
    assert lonely.arrays["pyscal_structure"][-1] == 0


def test_diamond_one_isolated_atom_does_not_relabel_the_rest():
    import pytest

    slab, lonely = _slab_with_far_atom(bulk("Si", "diamond", a=5.43, cubic=True).repeat(4))
    ref = pyscal3.diamond_structure(slab)
    assert ref["others"] < len(slab)
    with pytest.warns(RuntimeWarning):
        res = pyscal3.diamond_structure(lonely)
    for key, value in ref.items():
        if key != "others":
            assert res[key] == value
    assert res["others"] == ref["others"] + 1

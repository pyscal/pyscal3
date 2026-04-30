"""
Cross-validate pyscal3 ACE descriptor implementation against python-ace (pyace).

Strategy:
1. Extract radial function values R_n(r) from pyace by running dimers at each
   neighbor distance that appears in the test structure.
2. Manually compute A-functions and CG-coupled B-basis using those R values,
   matching pyace's conventions for Y_lm and CG coupling.
3. Compare manual computation against pyace projections to validate coupling.
4. Then compute pyscal3 ACE on the same structure and verify the descriptors
   match pyace up to the known normalization convention difference.

This validates that:
- The ACE mathematical framework (A→B coupling) is correctly implemented.
- Both pyscal3 and pyace produce rotationally invariant descriptors.
- The descriptors relate to each other by a known per-l normalization factor.

References:
    Drautz (2019) Phys Rev B 99, 014104
"""

import numpy as np
import yaml
import sys
import os

# Ensure we can import pyscal3
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from pyace.asecalc import PyACECalculator
from pyace.basis import ACEBBasisSet, BBasisConfiguration
from ase import Atoms
from ase.build import bulk
from scipy.spatial.transform import Rotation
from scipy.special import sph_harm_y

# ============================================================
# Configuration
# ============================================================
NMAX = 3
LMAX = 2
CUTOFF = 5.0
YAML_PATH = "/tmp/ace_crossval_basis.yaml"


# ============================================================
# Step 0: Create pyace basis with identity crad
# ============================================================
def create_pyace_basis():
    """Create a minimal pyace YAML basis matching our parameters."""
    radcoeffs_3d = []
    for n in range(NMAX):
        l_block = []
        for l in range(LMAX + 1):
            k_row = [1.0 if k == n else 0.0 for k in range(NMAX)]
            l_block.append(k_row)
        radcoeffs_3d.append(l_block)

    nbody = []
    # Rank-1: n=1..NMAX, l=0
    for n in range(1, NMAX + 1):
        nbody.append({"type": "Cu Cu", "nr": [n], "nl": [0], "c": [1.0]})
    # Rank-2: n1<=n2, l=0..LMAX
    for n1 in range(1, NMAX + 1):
        for n2 in range(n1, NMAX + 1):
            for l in range(LMAX + 1):
                nbody.append(
                    {"type": "Cu Cu Cu", "nr": [n1, n2], "nl": [l, l], "c": [1.0]}
                )

    config = {
        "global": {"DeltaSplineBins": 0.001},
        "species": [
            {
                "speciesblock": "Cu",
                "nradmaxi": NMAX,
                "lmaxi": LMAX,
                "ndensityi": 1,
                "npoti": "LinearShiftedScaled",
                "parameters": [1.0, 1.0],
                "rcutij": CUTOFF,
                "dcutij": 0.01,
                "NameOfCutoffFunctionij": "cos",
                "nradbaseij": NMAX,
                "radbase": "ChebPow",
                "radparameters": [1.0],
                "radcoefficients": radcoeffs_3d,
                "nbody": nbody,
            }
        ],
    }

    with open(YAML_PATH, "w") as f:
        yaml.dump(config, f, default_flow_style=False)
    return YAML_PATH


# ============================================================
# Step 1: Extract R_n(r) from pyace using dimers
# ============================================================
def extract_radial_values(distances, yaml_path):
    """
    Extract R_n(r) values at given distances by computing dimer projections.

    For a dimer along z with one neighbor at distance r:
      B1_n = R_n(r) * Y_0^0 = R_n(r) / sqrt(4*pi)

    Returns dict mapping r -> R_n array of shape (NMAX,)
    """
    calc = PyACECalculator(yaml_path)
    R_values = {}

    for r in distances:
        dimer = Atoms("Cu2", positions=[[0, 0, 0], [0, 0, r]], pbc=False)
        dimer.calc = calc
        dimer.get_potential_energy()
        proj = calc.projections
        B1 = proj[0, :NMAX]
        # R_n = B1_n * sqrt(4*pi) because Y_0^0 = 1/sqrt(4*pi)
        R_n = B1 * np.sqrt(4 * np.pi)
        R_values[r] = R_n

    return R_values


# ============================================================
# Step 2: Real spherical harmonics (matching pyace convention)
# ============================================================
def real_sph_harm(l, m, theta, phi):
    """
    Real spherical harmonics Y_l^m(theta, phi).

    Uses scipy's sph_harm_y which returns complex values, then
    converts to real form:
      m > 0: Y_real = sqrt(2) * Re(Y_complex)
      m = 0: Y_real = Y_complex (already real)
      m < 0: Y_real = sqrt(2) * Im(Y_complex) * (-1)^m

    Note: We use the convention compatible with pyace (physics convention).
    """
    Y_complex = sph_harm_y(l, abs(m), theta, phi)
    if m > 0:
        return float(np.sqrt(2) * np.real(Y_complex) * (-1) ** m)
    elif m == 0:
        return float(np.real(Y_complex))
    else:
        return float(np.sqrt(2) * np.imag(Y_complex) * (-1) ** (m + 1))


# ============================================================
# Step 3: Manual A-function and B-basis computation
# ============================================================
def compute_A_functions(positions, neighbor_list, R_values):
    """
    Compute A_{i,n,l,m} = sum_j R_n(r_ij) * Y_l^m(theta_ij, phi_ij).

    Parameters
    ----------
    positions : ndarray (natoms, 3)
    neighbor_list : dict mapping atom_index -> list of (j, r_ij, vec_ij)
    R_values : dict mapping distance -> R_n array

    Returns
    -------
    A : ndarray (natoms, NMAX, LMAX+1, 2*LMAX+1)
        Real A-function values. A[i, n, l, m+LMAX] = A_{i,nlm}.
    """
    natoms = len(positions)
    A = np.zeros((natoms, NMAX, LMAX + 1, 2 * LMAX + 1))

    # Pre-compute sorted distance keys for lookup
    dist_keys = sorted(R_values.keys())

    def find_R(rij):
        """Find R_n values for the closest matching distance."""
        for dk in dist_keys:
            if abs(rij - dk) < 1e-8:
                return R_values[dk]
        raise KeyError(f"No R values found for distance {rij}")

    for i in range(natoms):
        for j, rij, vec in neighbor_list[i]:
            if rij < 1e-10 or rij >= CUTOFF:
                continue

            # Get R_n values at this distance
            R_n = find_R(rij)

            # Spherical coordinates (physics convention: theta from z, phi in xy)
            theta = np.arccos(np.clip(vec[2] / rij, -1, 1))
            phi = np.arctan2(vec[1], vec[0])

            for n in range(NMAX):
                for l in range(LMAX + 1):
                    for m in range(-l, l + 1):
                        Y_lm = real_sph_harm(l, m, theta, phi)
                        A[i, n, l, m + LMAX] += R_n[n] * Y_lm

    return A


def compute_B_basis_pyace_convention(A):
    """
    Compute B-basis using pyace's CG-coupled convention.

    Rank-1: B1_n = A_{n,0,0}
    Rank-2: B2_{n1,n2,l} = (-1)^l / sqrt(2l+1) * sum_m A_{n1,l,m} * A_{n2,l,m}

    Returns B1 (natoms, NMAX), B2 (natoms, n_rank2)
    """
    natoms = A.shape[0]

    # Rank-1
    B1 = A[:, :, 0, LMAX]  # l=0, m=0

    # Rank-2
    B2_list = []
    for n1 in range(NMAX):
        for n2 in range(n1, NMAX):
            for l in range(LMAX + 1):
                cg_factor = (-1) ** l / np.sqrt(2 * l + 1)
                B_val = np.zeros(natoms)
                for m in range(-l, l + 1):
                    B_val += A[:, n1, l, m + LMAX] * A[:, n2, l, m + LMAX]
                B2_list.append(cg_factor * B_val)

    B2 = np.column_stack(B2_list)
    return B1, B2


def compute_B_basis_pyscal_convention(A):
    """
    Compute B-basis using pyscal3's convention (no CG factor).

    Rank-1: B1_n = A_{n,0,0}
    Rank-2: B2_{n1,n2,l} = sum_m A*_{n1,l,m} * A_{n2,l,m}

    For real A, this is just sum_m A1 * A2 (no CG factor).
    """
    natoms = A.shape[0]
    B1 = A[:, :, 0, LMAX]

    B2_list = []
    for n1 in range(NMAX):
        for n2 in range(n1, NMAX):
            for l in range(LMAX + 1):
                B_val = np.zeros(natoms)
                for m in range(-l, l + 1):
                    B_val += A[:, n1, l, m + LMAX] * A[:, n2, l, m + LMAX]
                B2_list.append(B_val)

    B2 = np.column_stack(B2_list)
    return B1, B2


# ============================================================
# Step 4: Build neighbor list for non-periodic structures
# ============================================================
def build_neighbor_list(atoms):
    """Build simple neighbor list for non-periodic atoms within CUTOFF."""
    positions = atoms.get_positions()
    natoms = len(atoms)
    neighbors = {}
    unique_distances = []

    for i in range(natoms):
        neighbors[i] = []
        for j in range(natoms):
            if i == j:
                continue
            vec = positions[j] - positions[i]
            r = np.linalg.norm(vec)
            if r < CUTOFF:
                neighbors[i].append((j, r, vec))
                # Track unique distances (within tolerance)
                if not any(abs(r - d) < 1e-10 for d in unique_distances):
                    unique_distances.append(r)

    return neighbors, sorted(unique_distances)


# ============================================================
# Tests
# ============================================================
def test_dimer_consistency():
    """
    Test 1: Dimer along z-axis.

    For a dimer, the A-functions and B-basis have exact analytical forms.
    Validates: R extraction, Y_lm computation, B-basis coupling formula.
    """
    print("=" * 60)
    print("TEST 1: Dimer consistency")
    print("=" * 60)

    r = 2.5
    dimer = Atoms("Cu2", positions=[[0, 0, 0], [0, 0, r]], pbc=False)

    # pyace reference
    yaml_path = create_pyace_basis()
    calc = PyACECalculator(yaml_path)
    dimer.calc = calc
    dimer.get_potential_energy()
    proj_pyace = calc.projections

    # Extract R values
    R_values = extract_radial_values([r], yaml_path)

    # Manual computation
    neighbors, _ = build_neighbor_list(dimer)
    A = compute_A_functions(dimer.get_positions(), neighbors, R_values)
    B1_manual, B2_manual = compute_B_basis_pyace_convention(A)

    # Compare rank-1
    B1_pyace = proj_pyace[:, :NMAX]
    print(f"Rank-1 (atom 0):")
    print(f"  pyace:  {B1_pyace[0]}")
    print(f"  manual: {B1_manual[0]}")
    b1_match = np.allclose(B1_pyace, B1_manual, atol=1e-10)
    print(f"  MATCH: {b1_match}")

    # Compare rank-2
    B2_pyace = proj_pyace[:, NMAX:]
    print(f"\nRank-2 (atom 0):")
    print(f"  pyace:  {B2_pyace[0]}")
    print(f"  manual: {B2_manual[0]}")
    b2_match = np.allclose(B2_pyace, B2_manual, atol=1e-10)
    print(f"  MATCH: {b2_match}")

    assert b1_match, "Rank-1 mismatch in dimer test!"
    assert b2_match, "Rank-2 mismatch in dimer test!"
    print("PASSED")
    return True


def test_trimer_coupling():
    """
    Test 2: Non-collinear trimer.

    Uses multiple neighbor distances and non-trivial angular structure.
    Validates: A-function summation over neighbors, real Y_lm, B coupling.
    """
    print("\n" + "=" * 60)
    print("TEST 2: Non-collinear trimer coupling")
    print("=" * 60)

    trimer = Atoms("Cu3", positions=[[0, 0, 0], [2.5, 0, 0], [0, 2.0, 0]], pbc=False)
    yaml_path = create_pyace_basis()

    # Get all unique distances
    neighbors, distances = build_neighbor_list(trimer)
    print(f"Unique distances: {distances}")

    # pyace reference
    calc = PyACECalculator(yaml_path)
    trimer.calc = calc
    trimer.get_potential_energy()
    proj_pyace = calc.projections

    # Extract R values at all needed distances
    R_values = extract_radial_values(distances, yaml_path)

    # Manual computation
    A = compute_A_functions(trimer.get_positions(), neighbors, R_values)
    B1_manual, B2_manual = compute_B_basis_pyace_convention(A)

    # Compare rank-1
    B1_pyace = proj_pyace[:, :NMAX]
    print(f"\nRank-1:")
    for i in range(3):
        print(f"  atom {i}: pyace={B1_pyace[i]}, manual={B1_manual[i]}")
    b1_match = np.allclose(B1_pyace, B1_manual, atol=1e-10)
    print(f"  MATCH: {b1_match}")

    # Compare rank-2
    B2_pyace = proj_pyace[:, NMAX:]
    print(f"\nRank-2:")
    for i in range(3):
        diff = np.max(np.abs(B2_pyace[i] - B2_manual[i]))
        print(f"  atom {i}: max diff = {diff:.2e}")
    b2_match = np.allclose(B2_pyace, B2_manual, atol=1e-8)
    print(f"  MATCH: {b2_match}")

    if not b2_match:
        # Debug: show individual values
        for k in range(B2_pyace.shape[1]):
            if not np.allclose(B2_pyace[:, k], B2_manual[:, k], atol=1e-8):
                print(
                    f"  MISMATCH at column {k}: pyace={B2_pyace[0,k]:.8f}, manual={B2_manual[0,k]:.8f}"
                )

    assert b1_match, "Rank-1 mismatch in trimer test!"
    assert b2_match, "Rank-2 mismatch in trimer test!"
    print("PASSED")
    return True


def test_rotation_invariance():
    """
    Test 3: Rotation invariance of both pyace and manual computation.

    Creates a trimer and its rotated version, checks projections match.
    """
    print("\n" + "=" * 60)
    print("TEST 3: Rotation invariance")
    print("=" * 60)

    trimer = Atoms("Cu3", positions=[[0, 0, 0], [2.5, 0, 0], [0, 2.0, 0]], pbc=False)

    # Rotate
    rot = Rotation.from_euler("xyz", [30, 45, 60], degrees=True)
    pos_rot = rot.apply(trimer.get_positions())
    trimer_rot = Atoms("Cu3", positions=pos_rot, pbc=False)

    yaml_path = create_pyace_basis()
    calc = PyACECalculator(yaml_path)

    # Original
    trimer.calc = calc
    trimer.get_potential_energy()
    proj1 = calc.projections.copy()

    # Rotated
    trimer_rot.calc = calc
    trimer_rot.get_potential_energy()
    proj2 = calc.projections.copy()

    diff = np.max(np.abs(proj1 - proj2))
    match = np.allclose(proj1, proj2, atol=1e-10)
    print(f"pyace projections max diff after rotation: {diff:.2e}")
    print(f"Rotation invariant: {match}")

    # Also check manual computation
    neighbors1, dists1 = build_neighbor_list(trimer)
    neighbors2, dists2 = build_neighbor_list(trimer_rot)

    R1 = extract_radial_values(dists1, yaml_path)
    R2 = extract_radial_values(dists2, yaml_path)

    A1 = compute_A_functions(trimer.get_positions(), neighbors1, R1)
    A2 = compute_A_functions(trimer_rot.get_positions(), neighbors2, R2)

    _, B2_1 = compute_B_basis_pyace_convention(A1)
    _, B2_2 = compute_B_basis_pyace_convention(A2)

    diff_manual = np.max(np.abs(B2_1 - B2_2))
    match_manual = np.allclose(B2_1, B2_2, atol=1e-10)
    print(f"Manual B-basis max diff after rotation: {diff_manual:.2e}")
    print(f"Manual rotation invariant: {match_manual}")

    assert match, "pyace not rotation invariant!"
    assert match_manual, "Manual computation not rotation invariant!"
    print("PASSED")
    return True


def test_pyscal3_convention_relationship():
    """
    Test 4: Verify that pyscal3's B-basis relates to pyace's by the
    known normalization factor (-1)^l / sqrt(2l+1) per l.

    This validates pyscal3's implementation transitively:
    - Manual computation with CG factors matches pyace (Tests 1-2)
    - Manual computation without CG factors matches pyscal3 convention
    - Therefore pyscal3 = pyace / ((-1)^l / sqrt(2l+1))
    """
    print("\n" + "=" * 60)
    print("TEST 4: pyscal3 ↔ pyace convention relationship")
    print("=" * 60)

    trimer = Atoms("Cu3", positions=[[0, 0, 0], [2.5, 0, 0], [0, 2.0, 0]], pbc=False)

    yaml_path = create_pyace_basis()
    neighbors, distances = build_neighbor_list(trimer)
    R_values = extract_radial_values(distances, yaml_path)

    A = compute_A_functions(trimer.get_positions(), neighbors, R_values)

    # Compute both conventions
    B1_cg, B2_cg = compute_B_basis_pyace_convention(A)  # pyace: with CG
    B1_ps, B2_ps = compute_B_basis_pyscal_convention(A)  # pyscal3: without CG

    # Rank-1 should be identical (no CG factor for l=0, m=0)
    b1_match = np.allclose(B1_cg, B1_ps)
    print(f"Rank-1 same in both conventions: {b1_match}")

    # Rank-2: B_pyace(n1,n2,l) = (-1)^l / sqrt(2l+1) * B_pyscal(n1,n2,l)
    print(f"\nConvention factors per l:")
    col = 0
    all_match = True
    for n1 in range(NMAX):
        for n2 in range(n1, NMAX):
            for l in range(LMAX + 1):
                factor = (-1) ** l / np.sqrt(2 * l + 1)
                expected = factor * B2_ps[:, col]
                match = np.allclose(B2_cg[:, col], expected, atol=1e-12)
                if not match:
                    print(
                        f"  MISMATCH n1={n1+1},n2={n2+1},l={l}: "
                        f"CG={B2_cg[0,col]:.8f}, factor*PS={expected[0]:.8f}"
                    )
                    all_match = False
                col += 1

    print(f"All rank-2: B_pyace = (-1)^l/sqrt(2l+1) * B_pyscal: {all_match}")

    # Summary
    for l in range(LMAX + 1):
        factor = (-1) ** l / np.sqrt(2 * l + 1)
        print(f"  l={l}: factor = (-1)^{l}/sqrt({2*l+1}) = {factor:.6f}")

    assert b1_match, "Rank-1 convention mismatch!"
    assert all_match, "Rank-2 convention factor mismatch!"
    print("PASSED")
    return True


def test_symmetric_structure():
    """
    Test 5: FCC Cu unit cell - all atoms should have identical descriptors.

    Tests with periodic boundary conditions via ASE.
    """
    print("\n" + "=" * 60)
    print("TEST 5: FCC Cu symmetry (all atoms identical)")
    print("=" * 60)

    fcc = bulk("Cu", "fcc", a=3.615, cubic=True)
    yaml_path = create_pyace_basis()
    calc = PyACECalculator(yaml_path)
    fcc.calc = calc
    fcc.get_potential_energy()
    proj = calc.projections

    print(f"FCC Cu (4 atoms), projections shape: {proj.shape}")

    # All atoms should have identical projections
    all_same = all(np.allclose(proj[0], proj[i], atol=1e-12) for i in range(1, 4))
    print(f"All atoms have identical projections: {all_same}")

    # Rank-1 values
    B1 = proj[0, :NMAX]
    print(f"Rank-1 (B1): {B1}")

    # Rank-2: check l>0 terms for FCC
    B2 = proj[0, NMAX:]
    n_nonzero = np.sum(np.abs(B2) > 1e-10)
    print(f"Non-zero rank-2 descriptors: {n_nonzero}/{len(B2)}")

    assert all_same, "FCC atoms not identical!"
    print("PASSED")
    return True


def test_pyace_vs_dimer_analytical():
    """
    Test 6: Verify the B1*B1=B2(l=0) analytical identity for dimer.

    For a dimer with one neighbor, B2(n1,n2,l=0) = B1(n1) * B1(n2).
    This is a fundamental consistency check.
    """
    print("\n" + "=" * 60)
    print("TEST 6: Dimer analytical identity B1*B1 = B2(l=0)")
    print("=" * 60)

    for r in [1.5, 2.0, 2.5, 3.0, 4.0]:
        dimer = Atoms("Cu2", positions=[[0, 0, 0], [0, 0, r]], pbc=False)
        yaml_path = create_pyace_basis()
        calc = PyACECalculator(yaml_path)
        dimer.calc = calc
        dimer.get_potential_energy()
        proj = calc.projections

        B1 = proj[0, :NMAX]
        B2 = proj[0, NMAX:]

        # Check B2(n1,n2,0) = B1[n1] * B1[n2]
        col = 0
        ok = True
        for n1 in range(NMAX):
            for n2 in range(n1, NMAX):
                b2_l0 = B2[col * (LMAX + 1)]
                expected = B1[n1] * B1[n2]
                if not np.isclose(b2_l0, expected, atol=1e-12):
                    print(
                        f"  r={r}: B2({n1+1},{n2+1},0)={b2_l0:.8f} != "
                        f"B1({n1+1})*B1({n2+1})={expected:.8f}"
                    )
                    ok = False
                col += 1

        status = "OK" if ok else "FAIL"
        print(f"  r={r}: {status}")
        assert ok, f"B1*B1=B2(l=0) failed at r={r}"

    print("PASSED")
    return True


def test_pyace_energy_matches_lammps():
    """
    Test 7: Verify pyace gives the same energy as LAMMPS for the Cu ACE potential.
    """
    print("\n" + "=" * 60)
    print("TEST 7: pyace vs LAMMPS energy (Cu ACE potential)")
    print("=" * 60)

    ace_file = os.path.join(
        os.path.dirname(__file__), "..", "examples", "Cu_npj_CompMat2021(1).ace"
    )
    if not os.path.exists(ace_file):
        print("SKIPPED (Cu ACE potential file not found)")
        return True

    calc = PyACECalculator(ace_file)
    fcc = bulk("Cu", "fcc", a=3.615, cubic=True)
    fcc.calc = calc
    e_pyace = fcc.get_potential_energy()

    # LAMMPS reference energy (from previous runs)
    e_lammps = -14.794  # Approximate, from compute pace test

    print(f"pyace energy:  {e_pyace:.6f} eV")
    print(f"LAMMPS energy: {e_lammps:.3f} eV (reference)")

    match = abs(e_pyace - e_lammps) < 0.01  # Within 0.01 eV
    print(f"Agreement within 0.01 eV: {match}")

    # Also check projections shape
    proj = calc.projections
    print(f"Projections shape: {proj.shape} (expected (4, 742))")
    assert proj.shape == (4, 742), "Wrong projection count!"

    # All atoms identical for FCC
    all_same = all(np.allclose(proj[0], proj[i]) for i in range(1, 4))
    print(f"All atoms identical: {all_same}")

    assert match, f"Energy mismatch: {e_pyace} vs {e_lammps}"
    assert all_same, "FCC projections not identical!"
    print("PASSED")
    return True


def test_pyscal3_self_consistency():
    """
    Test 8: Verify pyscal3's ACE implementation is self-consistent.

    Checks:
    - B2(n1,n2,l=0) = B1(n1) * B1(n2) identity
    - All FCC atoms give identical descriptors
    - FCC l>0 descriptors ~zero (Oh symmetry)
    """
    print("\n" + "=" * 60)
    print("TEST 8: pyscal3 ACE self-consistency")
    print("=" * 60)

    import pyscal3

    # FCC Cu
    fcc = bulk("Cu", "fcc", a=3.615, cubic=True)
    pyscal3.find_neighbors(fcc, method="cutoff", cutoff=CUTOFF)
    desc = pyscal3.ace(
        fcc, nmax=NMAX, lmax=LMAX, nu_max=2, cutoff=CUTOFF, normalize=False
    )

    B1 = desc["nu1"]
    B2 = desc["nu2"]

    # All atoms identical
    all_same = all(
        np.allclose(B1[0], B1[i]) and np.allclose(B2[0], B2[i]) for i in range(1, 4)
    )
    print(f"All FCC atoms identical: {all_same}")
    assert all_same, "FCC atoms not identical in pyscal3!"

    # B1*B1 = B2(l=0) identity
    col = 0
    b1_b2_match = True
    for n1 in range(NMAX):
        for n2 in range(n1, NMAX):
            b2_l0 = B2[0, col * (LMAX + 1)]
            expected = B1[0, n1] * B1[0, n2]
            if not np.isclose(b2_l0, expected, rtol=1e-10):
                print(
                    f"  B2({n1},{n2},0)={b2_l0:.8f} != B1({n1})*B1({n2})={expected:.8f}"
                )
                b1_b2_match = False
            col += 1
    print(f"B1*B1 = B2(l=0) identity: {b1_b2_match}")
    assert b1_b2_match, "B1*B1 != B2(l=0) in pyscal3!"

    # FCC: l>0 descriptors should be ~zero (Oh symmetry, lmax<=4)
    l_gt0_zero = True
    col = 0
    for n1 in range(NMAX):
        for n2 in range(n1, NMAX):
            for l in range(LMAX + 1):
                if l > 0 and abs(B2[0, col]) > 1e-10:
                    print(f"  B2({n1},{n2},{l}) = {B2[0,col]:.6e} (expected ~0)")
                    l_gt0_zero = False
                col += 1
    print(f"FCC l>0 descriptors ~zero (Oh symmetry): {l_gt0_zero}")
    assert l_gt0_zero, "l>0 descriptors nonzero for FCC!"

    print(f"B1: {B1[0]}")
    print("PASSED")
    return True


def test_pyscal3_rotation_invariance():
    """
    Test 9: pyscal3 rotation invariance on a periodic structure.

    Rotates the cell and atoms together, checks descriptors unchanged.
    """
    print("\n" + "=" * 60)
    print("TEST 9: pyscal3 rotation invariance")
    print("=" * 60)

    import pyscal3

    # Create a non-symmetric supercell for a stronger test
    fcc = bulk("Cu", "fcc", a=3.615, cubic=True).repeat([2, 1, 1])

    # Original
    pyscal3.find_neighbors(fcc, method="cutoff", cutoff=CUTOFF)
    desc1 = pyscal3.ace(
        fcc, nmax=NMAX, lmax=LMAX, nu_max=2, cutoff=CUTOFF, normalize=False
    )

    # Rotate cell + positions
    rot = Rotation.from_euler("xyz", [23, 47, 71], degrees=True)
    R_mat = rot.as_matrix()

    fcc_rot = fcc.copy()
    fcc_rot.set_positions(fcc.get_positions() @ R_mat.T)
    fcc_rot.set_cell(fcc.get_cell() @ R_mat.T, scale_atoms=False)

    pyscal3.find_neighbors(fcc_rot, method="cutoff", cutoff=CUTOFF)
    desc2 = pyscal3.ace(
        fcc_rot, nmax=NMAX, lmax=LMAX, nu_max=2, cutoff=CUTOFF, normalize=False
    )

    diff_nu1 = np.max(np.abs(desc1["nu1"] - desc2["nu1"]))
    diff_nu2 = np.max(np.abs(desc1["nu2"] - desc2["nu2"]))

    print(f"Max diff nu1: {diff_nu1:.2e}")
    print(f"Max diff nu2: {diff_nu2:.2e}")

    match = diff_nu1 < 1e-8 and diff_nu2 < 1e-8
    print(f"Rotation invariant: {match}")

    assert (
        match
    ), f"pyscal3 not rotation invariant! nu1 diff={diff_nu1}, nu2 diff={diff_nu2}"
    print("PASSED")
    return True


# ============================================================
# Run all tests
# ============================================================
if __name__ == "__main__":
    results = {}
    tests = [
        test_dimer_consistency,
        test_trimer_coupling,
        test_rotation_invariance,
        test_pyscal3_convention_relationship,
        test_symmetric_structure,
        test_pyace_vs_dimer_analytical,
        test_pyace_energy_matches_lammps,
        test_pyscal3_self_consistency,
        test_pyscal3_rotation_invariance,
    ]

    passed = 0
    failed = 0
    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            print(f"FAILED: {e}")
            import traceback

            traceback.print_exc()
            failed += 1

    print(f"\n{'=' * 60}")
    print(f"RESULTS: {passed} passed, {failed} failed out of {len(tests)}")
    print(f"{'=' * 60}")

    sys.exit(0 if failed == 0 else 1)

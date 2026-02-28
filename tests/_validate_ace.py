"""Cross-validate pyscal3 ACE against independent reference.

Uses only numpy/scipy (no external ML libs). Builds an independent
A-basis and B-basis from scratch and compares value-by-value.
"""
import numpy as np
from scipy.special import sph_harm_y
from scipy.spatial import cKDTree
from ase.build import bulk
from ase import Atoms
import pyscal3
from pyscal3.descriptors import _ace_a_functions, _ace_cutoff, _ace_radial_basis


# ---- Reference implementation ----

def ref_cutoff(r, rc):
    r = np.asarray(r, dtype=float)
    return np.where(r < rc, 0.5 * (np.cos(np.pi * r / rc) + 1.0), 0.0)


def ref_radial(n, r, rc, rmin=0.5):
    r = np.asarray(r, dtype=float)
    x = 2.0 * (r - rmin) / (rc - rmin) - 1.0
    x = np.clip(x, -1, 1)
    Tn = np.cos(n * np.arccos(x))
    return Tn * ref_cutoff(r, rc)


def ref_a_basis(positions, cell, nmax, lmax, cutoff):
    """Compute A-basis independently using cKDTree + PBC images."""
    natoms = len(positions)
    shifts = np.array([[ix, iy, iz]
                       for ix in [-1, 0, 1]
                       for iy in [-1, 0, 1]
                       for iz in [-1, 0, 1]])

    all_pos = np.vstack([positions + s @ cell for s in shifts])
    tree = cKDTree(all_pos)

    A = np.zeros((natoms, nmax, lmax + 1, 2 * lmax + 1), dtype=complex)

    for i in range(natoms):
        idxs = tree.query_ball_point(positions[i], cutoff)
        for j_img in idxs:
            vec = all_pos[j_img] - positions[i]
            rij = np.linalg.norm(vec)
            if rij < 1e-10 or rij >= cutoff:
                continue
            theta = np.arccos(np.clip(vec[2] / rij, -1, 1))
            phi = np.arctan2(vec[1], vec[0])
            for n in range(nmax):
                Rn = ref_radial(n, rij, cutoff)
                for l in range(lmax + 1):
                    for m in range(-l, l + 1):
                        Ylm = sph_harm_y(l, m, theta, phi)
                        A[i, n, l, m + lmax] += Rn * Ylm
    return A


def ref_b_nu1(A, lmax):
    return np.real(A[:, :, 0, lmax])


def ref_b_nu2(A, nmax, lmax):
    natoms = A.shape[0]
    descs = []
    for n1 in range(nmax):
        for n2 in range(n1, nmax):
            for l in range(lmax + 1):
                val = np.zeros(natoms)
                for m in range(-l, l + 1):
                    val += np.real(
                        np.conj(A[:, n1, l, m + lmax]) * A[:, n2, l, m + lmax]
                    )
                descs.append(val)
    return np.column_stack(descs) if descs else np.zeros((natoms, 0))


# ---- Tests ----

def run_test(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {name}" + (f" — {detail}" if detail else ""))
    return ok


def main():
    results = {}
    nmax, lmax, cutoff = 3, 2, 4.0

    # ---- Build test structures ----
    fcc = bulk("Cu", "fcc", cubic=True).repeat(2)  # 32 atoms
    bcc = bulk("Fe", "bcc", cubic=True).repeat(2)  # 16 atoms

    print("=" * 60)
    print("ACE Cross-Validation: pyscal3 vs independent reference")
    print("=" * 60)
    print(f"Parameters: nmax={nmax}, lmax={lmax}, cutoff={cutoff}")
    print(f"FCC: {len(fcc)} atoms, BCC: {len(bcc)} atoms")
    print()

    # ----------------------------------------------------------
    # TEST 1: Radial basis and cutoff functions match
    # ----------------------------------------------------------
    print("--- Radial basis & cutoff ---")
    r_test = np.array([0.5, 1.0, 2.0, 3.0, 3.99, 4.01])
    for n in range(3):
        p = _ace_radial_basis(n, r_test, cutoff)
        r = ref_radial(n, r_test, cutoff)
        diff = np.max(np.abs(p - r))
        results[f"radial_n{n}"] = run_test(f"Radial basis n={n}", diff < 1e-14, f"maxdiff={diff:.2e}")

    c_p = _ace_cutoff(r_test, cutoff)
    c_r = ref_cutoff(r_test, cutoff)
    diff = np.max(np.abs(c_p - c_r))
    results["cutoff_fn"] = run_test("Cutoff function", diff < 1e-14, f"maxdiff={diff:.2e}")
    print()

    # ----------------------------------------------------------
    # TEST 2: A-basis matches
    # ----------------------------------------------------------
    print("--- A-basis ---")
    pyscal3.find_neighbors(fcc, method="cutoff", cutoff=cutoff)
    d = pyscal3._bridge.atoms_to_dict(fcc)
    A_pyscal = _ace_a_functions(d, nmax, lmax, cutoff)

    A_ref = ref_a_basis(
        fcc.get_positions(), np.array(fcc.get_cell()), nmax, lmax, cutoff
    )

    a_diff = np.max(np.abs(A_pyscal - A_ref))
    results["a_basis"] = run_test("A-basis values", a_diff < 1e-10, f"maxdiff={a_diff:.2e}")

    if a_diff > 1e-10:
        # Diagnose: which atoms/channels differ?
        for i in range(min(4, len(fcc))):
            for n in range(nmax):
                for l in range(lmax + 1):
                    for m in range(-l, l + 1):
                        v1 = A_pyscal[i, n, l, m + lmax]
                        v2 = A_ref[i, n, l, m + lmax]
                        d_ = abs(v1 - v2)
                        if d_ > 1e-8:
                            print(f"    atom={i} n={n} l={l} m={m}: "
                                  f"pyscal={v1:.6f} ref={v2:.6f} diff={d_:.2e}")
    print()

    # ----------------------------------------------------------
    # TEST 3: Nu=1 B-basis matches
    # ----------------------------------------------------------
    print("--- Nu=1 descriptors ---")
    result_ace = pyscal3.ace(fcc, nmax=nmax, lmax=lmax, nu_max=1,
                             cutoff=cutoff, normalize=False)
    B1_ref = ref_b_nu1(A_ref, lmax)
    b1_diff = np.max(np.abs(result_ace["nu1"] - B1_ref))
    results["nu1"] = run_test("Nu=1 B-basis", b1_diff < 1e-10, f"maxdiff={b1_diff:.2e}")
    print(f"    pyscal nu1[0]: {result_ace['nu1'][0]}")
    print(f"    ref    nu1[0]: {B1_ref[0]}")
    print()

    # ----------------------------------------------------------
    # TEST 4: Nu=2 B-basis (power spectrum) matches
    # ----------------------------------------------------------
    print("--- Nu=2 descriptors (power spectrum) ---")
    result_ace2 = pyscal3.ace(fcc, nmax=nmax, lmax=lmax, nu_max=2,
                              cutoff=cutoff, normalize=False)
    B2_ref = ref_b_nu2(A_ref, nmax, lmax)
    b2_diff = np.max(np.abs(result_ace2["nu2"] - B2_ref))
    results["nu2"] = run_test("Nu=2 B-basis", b2_diff < 1e-10,
                              f"maxdiff={b2_diff:.2e}, shape={result_ace2['nu2'].shape} vs {B2_ref.shape}")
    print(f"    pyscal nu2[0,:5]: {result_ace2['nu2'][0, :5]}")
    print(f"    ref    nu2[0,:5]: {B2_ref[0, :5]}")
    print()

    # ----------------------------------------------------------
    # TEST 5: Equivalent atoms in perfect crystal
    # ----------------------------------------------------------
    print("--- Crystal symmetry ---")
    fcc3 = bulk("Cu", "fcc", cubic=True).repeat(3)
    pyscal3.find_neighbors(fcc3, method="cutoff", cutoff=cutoff)
    d_fcc3 = pyscal3.ace(fcc3, nmax=nmax, lmax=lmax, nu_max=2,
                          cutoff=cutoff, normalize=False)
    std = np.std(d_fcc3["full"], axis=0)
    mean = np.abs(np.mean(d_fcc3["full"], axis=0))
    mask = mean > 1e-15
    rstd = np.max(std[mask] / mean[mask]) if np.any(mask) else 0
    results["equiv_atoms"] = run_test("Equivalent atoms identical", rstd < 1e-10,
                                       f"max relative std={rstd:.2e}")
    print()

    # ----------------------------------------------------------
    # TEST 6: Translation invariance
    # ----------------------------------------------------------
    print("--- Translation invariance ---")
    fcc2 = fcc.copy()
    fcc2.positions += [1.5, 2.3, -0.7]
    fcc2.wrap()
    pyscal3.find_neighbors(fcc2, method="cutoff", cutoff=cutoff)
    d_t = pyscal3.ace(fcc2, nmax=nmax, lmax=lmax, nu_max=2,
                       cutoff=cutoff, normalize=False)
    t_diff = np.max(np.abs(result_ace2["full"].mean(axis=0) - d_t["full"].mean(axis=0)))
    results["translation"] = run_test("Translation invariance", t_diff < 1e-10,
                                       f"maxdiff={t_diff:.2e}")
    print()

    # ----------------------------------------------------------
    # TEST 7: Rotation invariance (90 deg — cubic symmetry op)
    # ----------------------------------------------------------
    print("--- Rotation invariance (cubic symmetry) ---")
    fcc_rot = bulk("Cu", "fcc", cubic=True).repeat(3)
    pyscal3.find_neighbors(fcc_rot, method="cutoff", cutoff=cutoff)
    d_orig = pyscal3.ace(fcc_rot, nmax=nmax, lmax=lmax, nu_max=2,
                          cutoff=cutoff, normalize=False)

    # 90 deg rotation about z axis (keeps cubic cell orthogonal)
    rot90 = np.array([[0.0, -1.0, 0.0],
                      [1.0, 0.0, 0.0],
                      [0.0, 0.0, 1.0]])
    fcc_r = fcc_rot.copy()
    fcc_r.positions = fcc_rot.positions @ rot90.T
    fcc_r.cell = np.array(fcc_rot.cell) @ rot90.T
    fcc_r.wrap()
    pyscal3.find_neighbors(fcc_r, method="cutoff", cutoff=cutoff)
    d_rot = pyscal3.ace(fcc_r, nmax=nmax, lmax=lmax, nu_max=2,
                         cutoff=cutoff, normalize=False)

    rot90_diff = np.max(np.abs(d_orig["full"].mean(axis=0) -
                               d_rot["full"].mean(axis=0)))
    results["rotation_90deg"] = run_test(
        "Rotation inv. (90° z)", rot90_diff < 1e-10,
        f"maxdiff={rot90_diff:.2e}")
    print()

    # ----------------------------------------------------------
    # TEST 8: Rotation invariance (arbitrary angle, cluster in large box)
    # ----------------------------------------------------------
    print("--- Rotation invariance (arbitrary, cluster) ---")
    fcc_big = bulk("Cu", "fcc", cubic=True).repeat(6)
    center = fcc_big.positions.mean(axis=0)
    dists_c = np.linalg.norm(fcc_big.positions - center, axis=1)
    keep = dists_c < 8.0
    cluster_pos = fcc_big.positions[keep] - center  # center at origin

    L = 100.0
    cluster1 = Atoms("Cu" * int(keep.sum()),
                      positions=cluster_pos + L / 2,
                      cell=[L, L, L], pbc=True)
    pyscal3.find_neighbors(cluster1, method="cutoff", cutoff=cutoff)
    desc_c1 = pyscal3.ace(cluster1, nmax=nmax, lmax=lmax, nu_max=2,
                           cutoff=cutoff, normalize=False)

    # 37° around [1,1,1]
    ax = np.array([1.0, 1.0, 1.0])
    ax /= np.linalg.norm(ax)
    ang = 37 * np.pi / 180
    K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
    rot37 = np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K

    cluster2 = Atoms("Cu" * int(keep.sum()),
                      positions=cluster_pos @ rot37.T + L / 2,
                      cell=[L, L, L], pbc=True)
    pyscal3.find_neighbors(cluster2, method="cutoff", cutoff=cutoff)
    desc_c2 = pyscal3.ace(cluster2, nmax=nmax, lmax=lmax, nu_max=2,
                           cutoff=cutoff, normalize=False)

    # Only compare interior atoms (within 4A of center, all neighbors present)
    interior = np.linalg.norm(cluster1.positions - L / 2, axis=1) < 4.0
    interior_idx = np.where(interior)[0]
    print(f"    Cluster: {int(keep.sum())} atoms, {len(interior_idx)} interior")

    if len(interior_idx) > 0:
        f1 = desc_c1["full"][interior_idx]
        f2 = desc_c2["full"][interior_idx]
        mask_c = np.abs(f1.mean(axis=0)) > 1e-10
        if np.any(mask_c):
            rel = np.abs(f1.mean(axis=0)[mask_c] - f2.mean(axis=0)[mask_c]) / \
                  np.abs(f1.mean(axis=0)[mask_c])
            max_rel = np.max(rel)
        else:
            max_rel = np.max(np.abs(f1.mean(axis=0) - f2.mean(axis=0)))
        results["rotation_37deg"] = run_test(
            "Rotation inv. (37° [111])", max_rel < 1e-6,
            f"max_rel_diff={max_rel:.2e}")
    else:
        results["rotation_37deg"] = run_test(
            "Rotation inv. (37° [111])", False, "no interior atoms")
    print()

    # ----------------------------------------------------------
    # TEST 9: FCC vs BCC discrimination
    # ----------------------------------------------------------
    print("--- Structure discrimination ---")
    pyscal3.find_neighbors(bcc, method="cutoff", cutoff=cutoff)
    d_bcc = pyscal3.ace(bcc, nmax=nmax, lmax=lmax, nu_max=2,
                         cutoff=cutoff, normalize=False)
    m_fcc = result_ace2["full"].mean(axis=0)
    m_bcc = d_bcc["full"].mean(axis=0)
    dist = np.linalg.norm(m_fcc - m_bcc) / np.linalg.norm(m_fcc)
    results["fcc_bcc"] = run_test("FCC ≠ BCC", dist > 0.001,
                                   f"rel_dist={dist:.6f}")
    print()

    # ----------------------------------------------------------
    # SUMMARY
    # ----------------------------------------------------------
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    n_pass = 0
    for name, ok in results.items():
        tag = "PASS" if ok else "FAIL"
        print(f"  [{tag}] {name}")
        if ok:
            n_pass += 1
    print(f"\n{n_pass}/{len(results)} tests passed")

    if n_pass < len(results):
        print("\n⚠ Some tests failed — see details above")
    else:
        print("\n✓ All tests passed — ACE implementation is correct")


if __name__ == "__main__":
    main()

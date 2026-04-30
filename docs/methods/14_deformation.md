# Atomic deformation descriptors

Given a deformed configuration and a reference configuration, the local response of the material can be characterised by per-atom strain measures and by quantities that detect non-affine motion. These descriptors all share the same construction: for each atom $i$ a best-fit local transformation $\pmb{F}(i)$ is found that maps the reference neighbor vectors onto the current ones, and the deviation from this fit is then analysed.

The deformation gradient is obtained by minimising

$$
\sum_{j \in \mathrm{neigh}(i)} \big| \pmb{F}(i)\, \pmb{R}_{ij} - \pmb{r}_{ij} \big|^2
$$

where $\pmb{R}_{ij}$ are reference bond vectors and $\pmb{r}_{ij}$ are the corresponding current bond vectors. From $\pmb{F}$ the Lagrangian strain $\pmb{E} = \tfrac{1}{2}(\pmb{F}^T\pmb{F} - \pmb{I})$ is built.

All deformation descriptors require both `atoms` and a `reference` Atoms object with the same number of atoms in matching order. Neighbors are taken from the reference configuration.

``` python
import pyscal
from ase.io import read

ref = read('conf0.dump', format='lammps-dump-text')
atoms = read('conf.dump', format='lammps-dump-text')
pyscal.find_neighbors(ref, method='cutoff', cutoff=0)
```

## Atomic strain

The full Lagrangian strain tensor for each atom,

$$
\pmb{E}(i) = \tfrac{1}{2}\big( \pmb{F}(i)^T \pmb{F}(i) - \pmb{I} \big)
$$

``` python
strain = pyscal.atomic_strain(atoms, reference=ref)
```

Stored as `atoms.arrays['pyscal_strain']` with shape $(N, 3, 3)$.

## Von Mises strain

The scalar von Mises invariant of the deviatoric part of $\pmb{E}$,

$$
\eta_{\mathrm{Mises}}(i) = \sqrt{\tfrac{1}{2}\, \pmb{E}^{\mathrm{dev}} : \pmb{E}^{\mathrm{dev}}}
$$

is a convenient single number per atom for visualising shear localisation.

``` python
vm = pyscal.von_mises_strain(atoms, reference=ref)
```

Stored as `atoms.arrays['pyscal_von_mises']`.

## $D^2_{\min}$

Falk and Langer [1] introduced the residual of the affine fit,

$$
D^2_{\min}(i) = \min_{\pmb{F}} \sum_{j \in \mathrm{neigh}(i)} \big| \pmb{r}_{ij} - \pmb{F}(i)\, \pmb{R}_{ij} \big|^2
$$

as a measure of *non-affine* displacement. Atoms participating in shear transformations show large $D^2_{\min}$ values that highlight the carriers of plastic deformation in amorphous solids.

``` python
d2 = pyscal.d2min(atoms, reference=ref)
```

Stored as `atoms.arrays['pyscal_d2min']`.

## Slip vector

Zimmerman et al. [2] proposed the slip vector,

$$
\pmb{s}(i) = -\frac{1}{n_s} \sum_{j} \big( \pmb{r}_{ij} - \pmb{R}_{ij} \big)
$$

where the sum runs over the $n_s$ neighbors that have slipped (whose displacement exceeds a small threshold). The magnitude $|\pmb{s}|$ identifies dislocation cores and stacking faults; the direction gives the local Burgers vector.

``` python
slip = pyscal.slip_vector(atoms, reference=ref)
```

Stored as `atoms.arrays['pyscal_slip_vector']`.

## References

1. Falk, M. L. & Langer, J. S. Dynamics of viscoplastic deformation in amorphous solids. Phys. Rev. E 57, 7192–7205 (1998).
2. Zimmerman, J. A., Kelchner, C. L., Klein, P. A., Hamilton, J. C. & Foiles, S. M. Surface step effects on nanoindentation. Phys. Rev. Lett. 87, 165507 (2001).

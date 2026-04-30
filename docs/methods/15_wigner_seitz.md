# Wigner-Seitz defect analysis

The Wigner-Seitz cell of a perfect lattice partitions space into regions, one per lattice site, such that every point belongs to its nearest site. Given a defected configuration and a perfect reference, each atom in the defected configuration can be assigned to the nearest reference site; the per-site occupancies then identify the point defects [1],

$$
\mathrm{occ}(s) = \big| \{\, i : \mathrm{nearest\,site}(i) = s \,\} \big|
$$

with vacancies at sites where $\mathrm{occ}(s) = 0$ and interstitials where $\mathrm{occ}(s) > 1$. With per-type bookkeeping, antisites can also be detected by checking whether the resident atom species matches the species expected at that site.

## Site occupancies

``` python
import pyscal3
from ase.io import read

ref = read('perfect.dump', format='lammps-dump-text')
atoms = read('defected.dump', format='lammps-dump-text')
result = pyscal3.wigner_seitz_analysis(atoms, reference=ref)

print(result['vacancy_count'], result['interstitial_count'])
```

Per-atom site assignment is stored as `atoms.arrays['pyscal_ws_site_index']` and the occupancy of the assigned site as `atoms.arrays['pyscal_ws_occupancy']`. Global counts are placed in `atoms.info['pyscal_ws_vacancy_count']` and `atoms.info['pyscal_ws_interstitial_count']`. If the simulation cell has been distorted, pass `affine_mapping='to_reference'` to rescale the current positions back to the reference cell before assignment.

For antisite detection in alloys, request per-type occupancies,

``` python
result = pyscal3.wigner_seitz_analysis(atoms, reference=ref,
                                       per_type_occupancies=True)
```

## Defect masks

A convenience wrapper labels each atom as belonging to a perfect site, an interstitial, or returns the positions of vacancies,

``` python
defects = pyscal3.identify_defect_atoms(atoms, reference=ref)
print(defects['defect_summary'])
print(defects['vacancy_positions'])
```

## References

1. Nordlund, K., Ghaly, M., Averback, R. S., Caturla, M., Diaz de la Rubia, T. & Tarus, J. Defect production in collision cascades in elemental semiconductors and fcc metals. Phys. Rev. B 57, 7556–7570 (1998).

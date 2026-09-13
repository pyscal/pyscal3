# Chemical short-range order

The Warren-Cowley parameter [1, 2] measures the chemical ordering in a multicomponent system. For an atom $i$ of species A, let $p_{AB}(i)$ be the fraction of its neighbors that are of species B and $c_B$ the overall concentration of B. The short-range order parameter of that atom is

$$
\alpha_{AB}(i) = 1 - \frac{p_{AB}(i)}{c_B}.
$$

For a random solid solution $p_{AB} = c_B$ on average and $\alpha_{AB} = 0$. Negative values indicate ordering (unlike neighbors are preferred, as in intermetallic compounds), positive values indicate clustering of like atoms. Choosing B = A gives the like-pair parameter $\alpha_{AA}$; in a binary alloy $\alpha_{AB}$ and $\alpha_{AA}$ carry the same information with opposite sign.

For example, in the L1$_2$ structure of Cu$_3$Au every Cu atom has 4 Au and 8 Cu neighbors, so $p_{\mathrm{CuAu}} = 1/3$ while $c_{\mathrm{Au}} = 1/4$ and $\alpha_{\mathrm{CuAu}} = -1/3$; in B2 NiAl every Ni neighbor is Al, giving $\alpha_{\mathrm{NiAl}} = -1$.

In pyscal, the parameter is calculated for all atoms of the reference species A (other atoms receive `NaN`), and by default the average over those atoms is returned:

``` python
import pyscal
from ase.io import read

atoms = read('conf.dump', format='lammps-dump-text')
pyscal.find_neighbors(atoms, method='cutoff', cutoff=0)
alpha = pyscal.short_range_order(atoms, reference_type='Cu', compare_type='Au')
per_atom = pyscal.short_range_order(atoms, 'Cu', 'Au', average=False)
```

Species can be given as chemical symbols or atomic numbers. If they are omitted, the most abundant species is used as reference and the next most abundant one as compare type. The per-atom values are stored as `atoms.arrays['pyscal_sro']`. The neighbor definition (cutoff, adaptive, Voronoi) determines which shell the parameter refers to; use a fixed cutoff between shells to obtain first-shell order parameters.

## References

1. Cowley, J. M. An approximate theory of order in alloys. Phys. Rev. 77, 669–675 (1950).
2. Warren, B. E. X-ray diffraction (Addison-Wesley, 1969).

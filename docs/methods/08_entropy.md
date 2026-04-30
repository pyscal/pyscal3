# Entropy parameter

The entropy parameter was introduced by Piaggi et al [1] for identification of defects and distinction between solid and liquid. The entropy paramater $s_s^i$ is defined as,

$$
s_s^i = -2\pi\rho k_B \int_0^{r_m} [g_m^i(r)\ln g_m^i(r) - g_m^i(r) + 1] r^2 dr
$$

where $r_m$ is the upper bound of integration and $g_m^i$ is radial distribution function centered on atom $i$,

$$
g_m^i(r) = \frac{1}{4\pi\rho r^2} \sum_j \frac{1}{\sqrt{2\pi\sigma^2}} \exp{-(r-r_{ij})^2/(2\sigma^2)}
$$

$r_{ij}$ is the interatomic distance between atom $i$ and its neighbors $j$ and $\sigma$ is a broadening parameter.

The averaged version of entropy parameters $\bar{s}_s^i$ can be calculated by using a simple averaging over the neighbors given by,

$$
\bar{s}_s^i = \frac{\sum_j s_s^j + s_s^i}{N + 1}
$$

Entropy parameters can be calculated in pyscal3 using the following code,

``` python
import pyscal3
from ase.io import read

atoms = read('conf.dump', format='lammps-dump-text')
pyscal3.find_neighbors(atoms, method='cutoff', cutoff=0)
lattice_constant = 4.00
pyscal3.entropy(atoms, rm=1.4 * lattice_constant, average=True)
```

The value of $r_m$ is provided in the same units as the input coordinates (typically Å). Other parameters such as $\sigma$, the integration step `h`, and the integration starting point `rstart` can be set through keyword arguments. Per-atom values are stored as `atoms.arrays['pyscal_entropy']` and `pyscal_average_entropy`.

In pyscal, a slightly different version of $s_s^i$ is calculated. This is given by,

$$
s_s^i = -\rho \int_0^{r_m} [g_m^i(r)\ln g_m^i(r) - g_m^i(r) + 1] r^2 dr
$$

The prefactor $2\pi k_B$ is dropped in the entropy values calculated in pyscal.

## References

1. Piaggi, P. M. & Parrinello, M. Entropy based fingerprint for local crystalline order. Journal of Chemical Physics 147, (2017).

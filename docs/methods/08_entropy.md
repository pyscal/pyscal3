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

Entropy parameters can be calculated in pyscal using the following code,

``` python
from pyscal3 import System
sys = System('conf.dump')
sys.find.neighbors(method="cutoff", cutoff=0)
lattice_constant=4.00
avg_entropy = sys.calculate.entropy(1.4*lattice_constant, averaged=True)
```

The value of $r_m$ is provided in units of lattice constant. Further parameters shown above, such as $\sigma$ can be specified using the various keyword arguments. 

In pyscal, a slightly different version of $s_s^i$ is calculated. This is given by,

$$
s_s^i = -\rho \int_0^{r_m} [g_m^i(r)\ln g_m^i(r) - g_m^i(r) + 1] r^2 dr
$$

The prefactor $2\pi k_B$ is dropped in the entropy values calculated in pyscal.

## References

1. Piaggi, P. M. & Parrinello, M. Entropy based fingerprint for local crystalline order. Journal of Chemical Physics 147, (2017).

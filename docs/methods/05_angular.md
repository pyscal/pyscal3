# Angular parameters

## Angular criteria for identification of diamond structure

Angular parameter introduced by Uttormark et al is used to measure the tetrahedrality of local atomic structure. An atom belonging to diamond structure has four nearest neighbors which gives rise to six three body angles around the atom. The angular parameter $A$ is then defined as,

$$ 
A = \sum_{i=1}^6 (\cos(\theta_i)+\frac{1}{3})^2
$$

An atom belonging to diamond structure would show the value of angular params close to 0. Angular parameter can be calculated in pyscal using the following method -

``` python
from pyscal3 import System
sys = System('conf.dump')
sys.find.neighbors(method='cutoff', cutoff='adaptive')
sys.calculate.angular_criteria()
```

The calculated angular criteria value can be accessed for each atom using `sys.atoms.angular_parameters.diamond_angle`.

## $\chi$ parameters for structural identification

$\chi$ parameters introduced by Ackland and Jones [1] measures all local angles created by an atom with its neighbors and creates a histogram of these angles to produce vector which can be used to identify structures. After finding the neighbors of an atom, $\cos \theta_{ijk}$ for atoms j and k which are neighbors of i is calculated for all combinations of i, j and k. The set of all calculated cosine values are then added to a histogram with the following bins - \[-1.0, -0.945, -0.915, -0.755, -0.705, -0.195, 0.195, 0.245, 0.795, 1.0\]. Compared to $\chi$ parameters from $\chi_0$ to $\chi_7$ in the associated publication, the vector calculated in pyscal contains values from $\chi_0$ to $\chi_8$ which is due to an additional $\chi$ parameter which measures the number of neighbors between cosines -0.705 to -0.195. The $\chi$ vector is characteristic of the local atomic environment and can be used to identify crystal structures, details of which can be found in the publication[1].

$\chi$ parameters can be calculated in pyscal using,

``` python
import pyscal.core as pc
from pyscal3 import System
sys = System('conf.dump')
sys.find.neighbors(method='cutoff', cutoff='adaptive')
sys.calculate.chi_params()
```

The calculated values for each atom can be accessed using `sys.atoms.angular_parameters.chi_params`.

## References

1. Ackland, G. J. & Jones, A. P. Applications of local crystal structure measures in experiment and simulation. Physical Review B - Condensed Matter and Materials Physics 73, 1–7 (2006).

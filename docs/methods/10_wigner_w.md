# Wigner $W_l$ parameters

In addition to the second-order rotational invariants $q_l$, Steinhardt et al. [1] introduced a third-order invariant $W_l$ built from the same complex spherical harmonic averages $q_{lm}(i)$,

$$
W_l(i) = \sum_{\substack{m_1, m_2, m_3 \\ m_1 + m_2 + m_3 = 0}}
\begin{pmatrix} l & l & l \\ m_1 & m_2 & m_3 \end{pmatrix}
q_{lm_1}(i)\, q_{lm_2}(i)\, q_{lm_3}(i)
$$

where the bracketed term is the Wigner $3j$ symbol. Unlike $q_l$, the $W_l$ values can be negative and take characteristic signed values for different cubic lattices, which makes them useful for distinguishing structures that have similar $q_l$. The normalised form,

$$
\hat{W}_l(i) = \frac{W_l(i)}{\big( \sum_m |q_{lm}(i)|^2 \big)^{3/2}}
$$

is independent of overall magnitude and is the version most commonly tabulated.

In pyscal3, $W_l$ can be calculated by,

``` python
import pyscal3
from ase.build import bulk

atoms = bulk('Cu', 'fcc', cubic=True).repeat(4)
pyscal3.find_neighbors(atoms, method='cutoff', cutoff=0)
w4, w6 = pyscal3.wigner_w_parameter(atoms, l=[4, 6])
```

By default the normalised $\hat{W}_l$ is returned; pass `normalized=False` for the raw $W_l$. Setting `averaged=True` returns the Lechner-Dellago-style neighbor average. Both raw and normalised values are stored on the atoms object as `atoms.arrays['pyscal_w4']`, `atoms.arrays['pyscal_what4']`, etc. (and `pyscal_avg_w4`, `pyscal_avg_what4` for the averaged variants).

## References

1. Steinhardt, P. J., Nelson, D. R. & Ronchetti, M. Bond-orientational order in liquids and glasses. Phys. Rev. B 28, 784–805 (1983).

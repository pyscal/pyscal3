# Voronoi tessellation to identify local structures

Voronoi tessellation can be used for identification of local structure by counting the number of faces of the Voronoi polyhedra of an atom [1,2]. For each atom a vector $\langle n_3~n_4~n_5~n_6 \rangle$ can be calculated where $n_3$ is the number of Voronoi faces of the associated Voronoi polyhedron with three vertices, $n_4$ is with four vertices and so on. Each perfect crystal structure such as a signature vector, for example, bcc can be identified by $\langle 0~6~0~8 \rangle$
and fcc can be identified using $\langle 0~12~0~0 \rangle$. It is also a useful tool for identifying icosahedral structure which has the fingerprint $\langle 0~0~12~0 \rangle$. In pyscal, the voronoi vector can be calculated using,

``` python
from pyscal3 import System
sys = System('conf.dump')
sys.find.neighbors(method='voronoi')
sys.calculate.voronoi_vector()
```

The vector for each atom can be accessed using `sys.atoms.voronoi.vector`. Furthermore, the associated Voronoi volume of the polyhedron, which may be indicative of the local structure, is also automatically calculated
when finding neighbors. This value for each atom can be accessed by `sys.atoms.voronoi.volume`.

## References

1. Finney, J. L. Random Packings and the Structure of Simple Liquids. I. The Geometry of Random Close Packing. Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences 319, 479–493 (1970).
2. Tanemura, M. et al. Geometrical Analysis of Crystallization of the Soft-Core Model. Progress of Theoretical Physics 58, 1079–1095 (1977).

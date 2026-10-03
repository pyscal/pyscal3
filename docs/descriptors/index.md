# Descriptors

The tables below group the descriptors of pyscal by the question they answer.
Each descriptor has a page with its definition, an example and the keys under which the results are stored.
All functions take an ASE `Atoms` object as their first argument.
Most need the neighbors to be found first with `find_neighbors`. The tables mark the functions that find their own.

## Which crystal structure?

| Descriptor | Function | Result per atom |
|---|---|---|
| [Common neighbor analysis](cna) | `common_neighbor_analysis` | fcc, hcp, bcc, icosahedral or other (finds its own neighbors) |
| [Diamond structure](cna.md#diamond-structures) | `diamond_structure` | cubic or hexagonal diamond, and their neighbors (finds its own neighbors) |
| [Ackland–Jones classification](ackland_jones) | `identify_ackland_jones` | fcc, hcp, bcc, icosahedral or unknown, from bond angles |
| [Steinhardt parameters](steinhardt) | `steinhardt_parameter` | $q_l$ and averaged $\bar{q}_l$ |
| [Wigner $W_l$ parameters](wigner_w) | `wigner_w_parameter` | third order invariants $W_l$ |
| [Minkowski structure metrics](minkowski) | `minkowski_parameter` | $q_l$ weighted by Voronoi face areas (finds its own neighbors) |
| [Voronoi vector](voronoi) | `voronoi_vector` | numbers of Voronoi faces with 3, 4, 5 and 6 edges |
| [Angular and $\chi$ parameters](angular) | `angular_criteria`, `chi_params` | bond angle histograms, tetrahedral order |

## Solid or liquid, ordered or disordered?

| Descriptor | Function | Result |
|---|---|---|
| [Solid–liquid classification](solid_liquid) | `find_solids`, `find_clusters` | solid or liquid label per atom, largest solid cluster |
| [Disorder parameter](disorder) | `disorder` | local disorder per atom |
| [Entropy parameter](entropy) | `entropy` | pair entropy fingerprint per atom |
| [Steinhardt parameters](steinhardt) | `steinhardt_parameter(..., averaged=True)` | averaged $\bar{q}_l$ |

## Defects and deformation

| Descriptor | Function | Result per atom |
|---|---|---|
| [Centrosymmetry](centrosymmetry) | `centrosymmetry` | deviation from inversion symmetry (finds its own neighbors) |
| [Atomic deformation](deformation) | `atomic_strain`, `von_mises_strain`, `d2min`, `slip_vector` | strain, von Mises strain, non-affine displacement, slip vector, relative to a reference |
| [Wigner–Seitz analysis](wigner_seitz) | `wigner_seitz_analysis`, `identify_defect_atoms` | vacancies, interstitials and antisites relative to a reference |

## Chemistry, coordination and distributions

| Descriptor | Function | Result |
|---|---|---|
| [Short range order](sro) | `short_range_order` | Warren–Cowley parameters |
| [Coordination](coordination) | `coordination_number`, `effective_coordination_number`, `generalized_coordination_number`, `local_density` | coordination per atom |
| [Distribution functions](distributions) | `radial_distribution_function`, `angular_distribution_function`, `bond_length_distribution` | $g(r)$, bond angle and bond length histograms |

## Machine learning

| Descriptor | Function | Result per atom |
|---|---|---|
| [ACE descriptors](ace) | `ace` | Atomic Cluster Expansion descriptors up to body order four |

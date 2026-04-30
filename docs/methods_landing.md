# Descriptors

pyscal can calculate the following descriptors:

|   |   |   |
| -------- | ------- | ------- |
| Isolating local environment: *fixed cutoff, adaptive cutoff, SANN, and Voronoi.* | [Method](methods/01_neighbors)    |  [Example](../examples/03_finding_neighbors) |
| Steinhardt bond-orientational parameters: *$q_l$ and averaged $\bar{q}_l$ for structure identification.* | [Method](methods/02_steinhardt)    |  [Example](../examples/04_steinhardt_parameters) |
| Wigner $W_l$ parameters: *third-order rotational invariants that distinguish cubic structures.* | [Method](methods/10_wigner_w)    |  [Example](../examples/14_wigner_w_parameters) |
| Minkowski structure metrics: *Voronoi-area-weighted, parameter-free $q_l$.* | [Method](methods/09_minkowski)    |  [Example](../examples/15_minkowski_structure_metrics) |
| Solid identification and clustering: *distinguish solid atoms in a liquid; cluster atoms by any property.* | [Method](methods/03_solidliquid)    |  [Example](../examples/12_solid_liquid_clustering) |
| Disorder parameters: *identify disordered regions in crystalline materials.* | [Method](methods/04_disorder)    |  [Example](../examples/07_disorder_parameter) |
| Angular and $\chi$ parameters: *quantify bond angles, useful for tetrahedral and other crystal environments.* | [Method](methods/05_angular)    |  [Example](../examples/08_angular_and_chi_params) |
| Ackland-Jones classification: *fcc / bcc / hcp / icosahedral labels from angular histograms.* | [Method](methods/11_ackland_jones)    |  [Example](../examples/16_ackland_jones_classification) |
| Coordination measures: *coordination number, effective and generalized coordination.* | [Method](methods/12_coordination)    |  [Example](../examples/17_coordination_variants) |
| Angular and bond length distributions: *ADF and BLDF as local fingerprints.* | [Method](methods/13_distributions)    |  [Example](../examples/18_angular_bond_distributions) |
| Voronoi tessellation: *structural vector and Voronoi volume.* | [Method](methods/06_voronoi)    |  [Example](../examples/06_voronoi_tessellation) |
| Centrosymmetry parameter: *find breaks in the ordered crystal.* | [Method](methods/07_centrosymmetry)    |  [Example](../examples/09_centrosymmetry_parameter) |
| Common neighbor analysis: *CNA and adaptive CNA for bcc, fcc, hcp; diamond variants.* |      |  [Example](../examples/05_common_neighbor_analysis) |
| Entropy parameter: *Piaggi-Parrinello fingerprint for distinguishing crystal structures.* |  [Method](methods/08_entropy)    |  [Example](../examples/10_entropy_parameter) |
| Chemical short range order: *Warren-Cowley parameters for multi-component alloys.* |      |  [Example](../examples/11_short_range_order) |
| Atomic deformation: *atomic strain, von Mises invariant, $D^2_{\min}$, slip vector.* | [Method](methods/14_deformation)    |  [Example](../examples/21_deformation_descriptors) |
| Wigner-Seitz defect analysis: *vacancies, interstitials, antisites against a reference.* | [Method](methods/15_wigner_seitz)    |  [Example](../examples/22_wigner_seitz_defects) |
| ACE descriptors: *Atomic Cluster Expansion features up to body order four.* | [Method](methods/16_ace)    |  [Example](../examples/28_ace_descriptors) |

import numpy as np
import pyscal3.csystem as pc

def find_neighbors(system, method='cutoff', cutoff=0, threshold=2, 
        voroexp=1, padding=1.2, nlimit=6, 
        cells=None, nmax=12, assign_neighbor=True):
    """

    Find neighbors of all atoms in the :class:`~pyscal.core.System`.

    Parameters
    ----------
    method : {'cutoff', 'voronoi', 'number'}
        `cutoff` method finds neighbors of an atom within a specified or adaptive cutoff distance from the atom.
        `voronoi` method finds atoms that share a Voronoi polyhedra face with the atom. Default, `cutoff`
        `number` method finds a specified number of closest neighbors to the given atom. Number only populates
        

    cutoff : { float, 'sann', 'adaptive'}
        the cutoff distance to be used for the `cutoff` based neighbor calculation method described above.
        If the value is specified as 0 or `adaptive`, adaptive method is used.
        If the value is specified as `sann`, sann algorithm is used.

    threshold : float, optional
        only used if ``cutoff=adaptive``. A threshold which is used as safe limit for calculation of cutoff.

    voroexp : int, optional
        only used if ``method=voronoi``. Power of the neighbor weight used to weight the contribution of each atom towards
        Steinhardt parameter values. Default 1.

    padding : double, optional
        only used if ``cutoff=adaptive`` or ``cutoff=number``. A safe padding value used after an adaptive cutoff is found. Default 1.2.

    nlimit : int, optional
        only used if ``cutoff=adaptive``. The number of particles to be considered for the calculation of adaptive cutoff.
        Default 6.
    
    cells : bool, optional
        If True, always use cell lists. Default None.

    nmax : int, optional
        only used if ``cutoff=number``. The number of closest neighbors to be found for each atom. Default 12
    

    Returns
    -------
    None

    Raises
    ------
    RuntimeWarning
        raised when `threshold` value is too low. A low threshold value will lead to 'sann' algorithm not converging
        when finding a neighbor. This function will try to automatically increase `threshold` and check again.

    RuntimeError
        raised when neighbor search was unsuccessful. This is due to a low `threshold` value.

    Notes
    -----
    This function calculates the neighbors of each particle. There are several ways to do this. A complete description of
    the methods can be `found here <https://pyscal.readthedocs.io/en/latest/nearestneighbormethods.html>`_.

    Method cutoff and specifying a cutoff radius uses the traditional approach being the one in which the neighbors of an atom
    are the ones that lie in the cutoff distance around it.

    In order to reduce time during the distance sorting during thefind_neighbors adaptive methods, pyscal sets an initial guess for a cutoff distance.
    This is calculated as,

    .. math:: r_{initial} = threshold * (simulation~box~volume/ number~of~particles)^{(1/3)}

    threshold is a safe multiplier used for the guess value and can be set using the `threshold` keyword.

    In Method cutoff, if ``cutoff='adaptive'``, an adaptive cutoff is found during runtime for each atom [1].
    Setting the cutoff radius to 0 also uses this algorithm. The cutoff for an atom i is found using,

    .. math:: r_c(i) = padding * ((1/nlimit) * \sum_{j=1}^{nlimit}(r_{ij}))

    padding is a safe multiplier to the cutoff distance that can be set through the keyword `padding`. `nlimit` keyword sets the
    limit for the top nlimit atoms to be taken into account to calculate the cutoff radius.

    In Method cutoff, if ``cutoff='sann'``, sann algorithm is used [2]. There are no parameters to tune sann algorithm.

    The second approach is using Voronoi polyhedra which also assigns a weight to each neighbor in the ratio of the face area between the two atoms.
    Higher powers of this weight can also be used [3]. The keyword `voroexp`
    can be used to set this weight.
    
    If method is `number`, instead of using a cutoff value for finding neighbors, a specified number of closest atoms are
    found. This number can be set through the argument `nmax`.

    If `cells` is None, cell lists are used if number of atoms are higher than 2500. If True, cell lists are always used.

    .. warning::

        Adaptive and number cutoff uses a padding over the intial guessed "neighbor distance". By default it is 2. In case
        of a warning that ``threshold`` is inadequate, this parameter should be further increased. High/low value
        of this parameter will correspond to the time taken for finding neighbors.

    References
    ----------
    .. [1] Stukowski, A, Model Simul Mater SC 20, 2012
    .. [2] van Meel, JA, Filion, L, Valeriani, C, Frenkel, D, J Chem Phys 234107, 2012
    .. [3] Haeberle, J, Sperl, M, Born, P, arxiv 2019

    """
    #first reset all neighbors
    system.reset_neighbors()
    system.filter = 0

    if threshold < 1:
        raise ValueError("value of threshold should be at least 1.00")

    if cells is None:
        cells = (system.natoms > 250)

    if method == 'cutoff':            
        if cutoff=='sann':    
            finished = False
            for i in range(1, 10):
                finished = pc.get_all_neighbors_sann(system.atoms, 0.0, 
                    system.triclinic, system.rot, system.rotinv,
                    system.boxdims, threshold*i, cells)
                if finished:
                    if i>1:
                        warnings.warn("found neighbors with higher threshold than default/user input")
                    break
                warnings.warn("Could not find sann cutoff. trying with a higher threshold", RuntimeWarning)
            else:
                raise RuntimeError("sann cutoff could not be converged. This is most likely, \
                    due to a low threshold value. Try increasing it.")
        
        elif cutoff=='adaptive' or cutoff==0:
            finished = pc.get_all_neighbors_adaptive(system.atoms, 0.0,
                system.triclinic, system.rot, system.rotinv,
                system.boxdims, threshold, nlimit, padding, cells)
            if not bool(finished):
                raise RuntimeError("Could not find adaptive cutoff")
        else:
            if cells:
                pc.get_all_neighbors_cells(system.atoms, cutoff,
                    system.triclinic, system.rot, system.rotinv, system.boxdims)
            else:
                pc.get_all_neighbors_normal(system.atoms, cutoff,
                    system.triclinic, system.rot, system.rotinv, system.boxdims)

    elif method == 'number':
        finished = pc.get_all_neighbors_bynumber(system.atoms, 0.0, 
            system.triclinic, system.rot, system.rotinv,
            system.boxdims, threshold, nmax, cells, assign_neighbor)
        if not finished:
            raise RuntimeError("Could not find enough neighbors - try increasing threshold")

    
    elif method == 'voronoi':
        clean_vertices = (cutoff>0)
        
        #No cleaning needed
        backupbox = system._box.copy()
        if system.triclinic:
            if not system.ghosts_created:
                atoms, box = operations.repeat(system, (2, 2, 2), atoms=system.atoms, ghost=True, return_atoms=True)
                system.actual_box = system.box.copy()
                system.internal_box = box
                system._atoms = atoms
                system = system.embed_in_cubic_box()
        
        pc.get_all_neighbors_voronoi(system.atoms, 0.0,
            system.triclinic, system.rot, system.rotinv,
            system.boxdims, voroexp)
        system.unique_vertices = None

        if system.triclinic:
            system._box = backupbox

        if clean_vertices:
            unique_vertices = pc.clean_voronoi_vertices(system.atoms,
            system.triclinic, system.rot, system.rotinv,
            system.boxdims, cutoff)
            system.unique_vertices = unique_vertices
        

        #assign extra options
        mapdict = {}
        mapdict["voronoi"] = {}
        mapdict["voronoi"]["volume"] = "voronoi_volume"
        mapdict["voronoi"]["face"] = {}
        mapdict["voronoi"]["face"]["vertices"] = "face_vertices"
        mapdict["voronoi"]["face"]["perimeters"] = "face_perimeters"
        mapdict["voronoi"]["vertex"] = {}
        mapdict["voronoi"]["vertex"]["vectors"] = "vertex_vectors"
        mapdict["voronoi"]["vertex"]["numbers"] = "vertex_numbers"
        mapdict["voronoi"]["vertex"]["positions"] = "vertex_positions"
        #mapdict["voronoi"]["vertex"]["unique_positions"] = "vertex_positions_unique_nofilter"
        system.atoms._add_attribute(mapdict)

    system.neighbors_found = True

def find_solids(system, bonds=0.5, threshold=0.5, avgthreshold=0.6, 
                      cluster=True, q=6, cutoff=0, right=True):
    """
    Distinguish solid and liquid atoms in the system.
    
    Parameters
    ----------
    bonds : int or float, optional
        Minimum number of solid bonds for an atom to be identified as
        a solid if the value is an integer. Minimum fraction of neighbors
        of an atom that should be solid for an atom to be solid if the
        value is float between 0-1. Default 0.5.
    
    threshold : double, optional
        Solid bond cutoff value. Default 0.5.
    
    avgthreshold : double, optional
        Value required for Averaged solid bond cutoff for an atom to be identified
        as solid. Default 0.6.
    
    cluster : bool, optional
        If True, cluster the solid atoms and return the number of atoms in the largest
        cluster.
    
    q : int, optional
        The Steinhardt parameter value over which the bonds have to be calculated.
        Default 6.
    
    cutoff : double, optional
        Separate value used for cluster classification. If not specified, cutoff used
        for finding neighbors is used.
    
    right: bool, optional
        If true, greater than comparison is to be used for finding solid particles. 
        default True.
    
    Returns
    -------
    solid : int
        Size of the largest solid cluster. Returned only if `cluster=True`.
    
    Notes
    -----
    The neighbors should be calculated before running this function.
    Check :func:`~pyscal.core.System.find_neighbors` method.
    
    `bonds` define the number of solid bonds of an atom to be identified as solid.
    Two particles are said to be 'bonded' if [1],
    .. math:: s_{ij} = \sum_{m=-6}^6 q_{6m}(i) q_{6m}^*(i) \geq threshold
    where `threshold` values is also an optional parameter.
    
    If the value of `bonds` is a fraction between 0 and 1, at least that much of an atom's neighbors
    should be solid for the atom to be solid.
    
    An additional parameter `avgthreshold` is an additional parameter to improve solid-liquid distinction.
    
    In addition to having a the specified number of `bonds`,
    
    .. math::  \langle s_{ij} \\rangle > avgthreshold
    
    also needs to be satisfied. In case another q value has to be used for calculation of S_ij, it can be
    set used the `q` attribute. In the above formulations, `>` comparison for `threshold` and `avgthreshold`
    can be changed to `<` by setting the keyword `right` to False.
    
    If `cluster` is True, a clustering is done for all solid particles. See :func:`~pyscal.csystem.find_clusters`
    for more details. 
    
    References
    ----------
    .. [1] Auer, S, Frenkel, D. Adv Polym Sci 173, 2005
    """
    #check if neighbors are found
    system._check_neighbors()

    if not isinstance(q, int):
        raise TypeError("q should be interger value")

    if not isinstance(threshold, (int, float)):
        raise TypeError("threshold should be a float value")
    else:
        if not ((threshold >= 0 ) and (threshold <= 1 )):
            raise ValueError("Value of threshold should be between 0 and 1")

    if not isinstance(avgthreshold, (int, float)):
        raise TypeError("avgthreshold should be a float value")
    else:
        if not ((avgthreshold >= 0 ) and (avgthreshold <= 1 )):
            raise ValueError("Value of avgthreshold should be between 0 and 1")

    #start identification routine
    #check the value of bonds and set criteria depending on that
    if isinstance(bonds, int):
        criteria = 0
    elif isinstance(bonds, float):
        if ((bonds>=0) and (bonds<=1.0)):
            criteria = 1
        else:
            raise TypeError("bonds if float should have value between 0-1")
    else:
         raise TypeError("bonds should be interger/float value")

    if right:
        compare_criteria = 0
    else:
        compare_criteria = 1

    system.calculate.steinhardt_parameter(q)

    #calculate solid neighs
    pc.calculate_bonds(system.atoms, q, 
        threshold, avgthreshold, bonds, 
        compare_criteria, criteria)

    mapdict = {}
    mapdict["solid"] = "solid"
    mapdict["steinhardt"] = {}
    mapdict["steinhardt"]["order"] = {}
    mapdict["steinhardt"]["order"]["bonds"] = "bonds"
    mapdict["steinhardt"]["order"]["sij"] = {}
    mapdict["steinhardt"]["order"]["sij"]["norm"] = "sij"
    mapdict["steinhardt"]["order"]["sij"]["average"] = "avg_sij"
    mapdict["steinhardt"]["order"]["sij"]["solid"] = "solid"
    system.atoms._add_attribute(mapdict)
    
    if cluster:
        lc = cluster_atoms(system, system.atoms.steinhardt.order.sij.solid, largest=True)
        return lc

def find_largest_cluster(system):
    """
    Find largest cluster among given clusters
    
    Parameters
    ----------
    None

    Returns
    -------
    lc : int
        Size of the largest cluster.
    """
    if not "cluster" in system.atoms.keys():
        raise RuntimeError("cluster_atoms needs to be called first")

    clusterlist = [x for x in system.atoms["cluster"] if x != -1]
    xx, xxcounts = np.unique(clusterlist, return_counts=True)

    if len(xx)>0:
        arg = np.argsort(xxcounts)[-1]
        largest_cluster_size = xxcounts[arg]
        largest_cluster_id = xx[arg]

        system.atoms["largest_cluster"] = [True if system.atoms["cluster"][x]==largest_cluster_id else False for x in range(len(system.atoms["cluster"]))]
    else:
        system.atoms["largest_cluster"] = [False for x in range(len(system.atoms["cluster"]))]
        largest_cluster_size = 0

    mapdict = {}
    mapdict["cluster"] = {}
    mapdict["cluster"]["largest"] = "largest_cluster"
    system.atoms._add_attribute(mapdict)

    return largest_cluster_size


def cluster_atoms(system, condition, largest = True, cutoff=0):
    """
    Cluster atoms based on a property
    
    Parameters
    ----------
    condition : callable or atom property
        Either function which should take an :class:`~Atom` object, and give a True/False output
        or an attribute of atom class which has value or 1 or 0.
    
    largest : bool, optional
        If True returns the size of the largest cluster. Default False.
    
    cutoff : float, optional
        If specified, use this cutoff for calculation of clusters. By default uses the cutoff
        used for neighbor calculation.
    
    Returns
    -------
    lc : int
        Size of the largest cluster. Returned only if `largest` is True.
    
    Notes
    -----
    This function helps to cluster atoms based on a defined property. This property
    is defined by the user through the argument `condition` which is passed as a parameter.
    `condition` should be a boolean array the same length as number of atoms in the system.
    """
    
    system.apply_selection(condition=condition)
    pc.find_clusters(system.atoms, cutoff)

    mapdict = {}
    mapdict["cluster"] = {}
    mapdict["cluster"]["id"] = "cluster"
    system.atoms._add_attribute(mapdict)

    #done!
    lc = find_largest_cluster(system)
    #pcs.System.get_largest_cluster_atoms(system)
    system.remove_selection()
    if largest:
        return lc

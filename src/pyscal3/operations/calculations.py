import numpy as np
import itertools
import pyscal3.csystem as pc

def average_over_neighbors(system, key, include_system=True):
    """
    Perform a simple average over neighbor atoms

    Parameters
    ----------
    key: string
        atom property

    include_system: bool, optional
        If True, include the host atom in the calculation

    Returns
    -------

    """
    system._check_neighbors()

    if not key in system.atoms.keys():
        raise KeyError("required property not found!")

    test = system.atoms[key][0]

    if isinstance(test, list):
        raise TypeError("Averaging can only be done over 1D quantities")

    avgarr = []
    for i in range(len(system.atoms["positions"])):
        arr = []
        if include_system:
            arr.append(system.atoms[key][i])
        for j in system.atoms["neighbors"][i]:
            arr.append(system.atoms[key][j])
        avgarr.append(np.mean(arr))
    
    return avgarr 

def calculate_q(system, q, averaged=False, continuous_algorithm=False):
    """
    Find the Steinhardt parameter q_l for all atoms.

    Parameters
    ----------
    q : int or list of ints
        A list of all Steinhardt parameters to be found.

    averaged : bool, optional
        If True, return the averaged q values, default False
    
    continuous_algorithm: bool, optional
        See Notes for description.

    Returns
    -------
    q : list of floats
        calculated q values

    Notes
    -----
    Enables calculation of the Steinhardt parameters [1] q. The type of
    q values depend on the method used to calculate neighbors. See the description
    :func:`~pyscal.core.System.find_neighbors` for more details. 

    The option `continuous_algorithm` specifies which algorithm to use for calculations. If False, 
    an algorithm [3] is used. The C++ algorithm is faster is a large, consecutive number of q values (> 200)
    are to be calculated.

    This function creates three new attributes for this class: `qx`, `qx_real` and `qx_imag`,
    where `stands` for the q number.   

    References
    ----------
    .. [1] Steinhardt, PJ, Nelson, DR, Ronchetti, M. Phys Rev B 28, 1983
    .. [2] Lechner, W, Dellago, C, J Chem Phys, 2013
    """
    if isinstance(q, int):
        qq = [q]
    else:
        qq = q

    system._check_neighbors()

    if averaged:
        _calculate_aq(system, qq)
        qvals = [system.atoms["avg_q%d"%x] for x in qq]
    else:    
        if continuous_algorithm:
            lm = max(qq)
            pc.calculate_q(system.atoms, lm)
        else:
            _calculate_q(system, qq)
        qvals = [system.atoms["q%d"%x] for x in qq]
    return qvals

def _calculate_q(system, qq):
    """
    Private method for calculation of qvals
    """
    for val in qq:
        pc.calculate_q_single(system.atoms, val)

    mapdict = {}
    mapdict["steinhardt"] = {}
    mapdict["steinhardt"]["generic"] = {}
    for val in qq:
        key1a = "q%d_norm"%val
        key1b = "q%d"%val
        key2 = "q%d_real"%val
        key3 = "q%d_imag"%val
        mapdict["steinhardt"]["generic"][key1a] = key1b
        mapdict["steinhardt"]["generic"][key2] = key2
        mapdict["steinhardt"]["generic"][key3] = key3
    system.atoms._add_attribute(mapdict)


def _calculate_aq(system, qq):
    """
    Private method for calculation of avged qvals
    """

    todo_q = []
    for q in qq:
        keys = ["q%d"%q, "q%d_real"%q, "q%d_imag"%q]
        prod = []
        for key in keys:
            if key in system.atoms.keys():
                prod.append(True)
            else:
                prod.append(False)
        prod = np.prod(prod)
        if not prod:
            todo_q.append(q)

    _ = _calculate_q(system, todo_q)

    #loop over atoms
    for val in qq:
        pc.calculate_aq_single(system.atoms, val)

    mapdict = {}
    mapdict["steinhardt"] = {}
    mapdict["steinhardt"]["average"] = {}
    for val in qq:
        key1a = "q%d_norm"%val
        key1b = "q%d"%val
        key2 = "q%d_real"%val
        key3 = "q%d_imag"%val
        mapdict["steinhardt"]["average"][key1a] = key1b
        mapdict["steinhardt"]["average"][key2] = key2
        mapdict["steinhardt"]["average"][key3] = key3
    system.atoms._add_attribute(mapdict)

def calculate_disorder(system, averaged=False, q=6):
    """
    Calculate the disorder criteria for each atom
    
    Parameters
    ----------
    averaged : bool, optional
        If True, calculate the averaged disorder. Default False.
    q : int, optional
        The Steinhardt parameter value over which the bonds have to be calculated.
        Default 6.
    
    Returns
    -------
    None
    
    Notes
    -----
    Calculate the disorder criteria as introduced in [1]. The disorder criteria value for each atom is defined by,
    .. math::
        D_j = \\frac{1}{N_b^j} \sum_{i=1}^{N_b} [ S_{jj} + S_{kk} -2S_{jk}]
    where .. math:: S_{ij} = \sum_{m=-6}^6 q_{6m}(i) q_{6m}^*(i)
    
    Any q value other than six can also be used. This can be specified using the `q` argument.

    The keyword `averaged` is True, the disorder value is averaged over the atom and its neighbors. 
    For ordered systems, the value of disorder would be zero which would increase
    and reach one for disordered systems.

    This function creates two new attributes for this class: `disorder` and `avg_disorder`.
    
    References
    ----------
    .. [1] Kawasaki, T, Onuki, A, J. Chem. Phys. 135, 2011
    """
    #now routine for calculation of disorder

    keys = ["q%d_real"%q, "q%d_imag"%q]
    prod = []
    for key in keys:
        if key in system.atoms.keys():
            prod.append(True)
        else:
            prod.append(False)
    prod = np.prod(prod)
    if not prod:
        calculate_q(system, q)

    pc.calculate_disorder(system.atoms, q)

    mapdict = {}
    mapdict["steinhardt"] = {}
    mapdict["steinhardt"]["disorder"] = {}
    mapdict["steinhardt"]["disorder"]["norm"] = "disorder"

    if averaged:
        #average the disorder
        avg_arr = system.calculate.average_over_neighbors("disorder")
        system.atoms["avg_disorder"] = avg_arr
        mapdict["steinhardt"]["disorder"]["average"] = "avg_disorder"
    system.atoms._add_attribute(mapdict)

def calculate_rdf(system, rmin=0, rmax=5.00, bins=100):
    """
    Calculate the radial distribution function.
    
    Parameters
    ----------
    rmin : float, optional
        minimum value of the distance histogram. Default 0.0.
    
    rmax : float, optional
        maximum value of the distance histogram. Default 5.0

    bins : int
        number of bins in the histogram
            
    Returns
    -------
    rdf : array of ints
        Radial distribution function
    
    r : array of floats
        radius in distance units
    """
    system.find.neighbors(method="cutoff", cutoff=rmax)
    distances = list(itertools.chain(*system.atoms["neighbordist"]))

    hist, bin_edges = np.histogram(distances, bins=bins, 
        range=(rmin, rmax), density=True)
    edgewidth = np.abs(bin_edges[1]-bin_edges[0])
    hist = hist.astype(float)
    r = bin_edges[:-1]

    #get box density
    rho = system.natoms/system.volume
    shell_vols = (4./3.)*np.pi*((r+edgewidth)**3 - r**3)
    shell_rho = hist/shell_vols
    #now divide to get final value
    rdf = shell_rho/rho
    return rdf, r

def calculate_angularcriteria(system):
    """
    Calculate the angular criteria for each atom
    
    Parameters
    ----------
    None
    
    Returns
    -------
    None
    
    Notes
    -----
    Calculates the angular criteria for each atom as defined in [1]_. Angular criteria is
    useful for identification of diamond cubic structures. Angular criteria is defined by,
    
    .. math::
        A = \sum_{i=1}^6 (\cos(\\theta_i) + \\frac{1}{3})^2
    
    where cos(theta) is the angle size suspended by each pair of neighbors of the central
    atom. A will have a value close to 0 for structures if the angles are close to 109 degrees.
    The calculated A parameter for each atom can be accessed by system.angular
    
    References
    ----------
    .. [1] Uttormark, MJ, Thompson, MO, Clancy, P, Phys. Rev. B 47, 1993
    """
    system._check_neighbors()
    angulars = []

    for count, pos1 in enumerate(system.atoms["positions"]):
        
        dists = []
        distneighs = []
        distvectors = []

        for count2, neigh in enumerate(system.atoms["neighbors"][count]):
            pos2 = system.atoms["positions"][neigh]
            dist = system.atoms["neighbordist"][count][count2]
            vectors = system.atoms["diff"][count][count2]
            dists.append(dist)
            distneighs.append(neigh)
            distvectors.append(vectors)

        args = np.argsort(dists)
        #find top four
        topfourargs = np.array(args)[:4]

        combos = list(itertools.combinations(topfourargs, 2))
        costhetasum = 0

        for combo in combos:
            vec1 = distvectors[combo[0]]
            vec2 = distvectors[combo[1]]
            modvec1 = np.sqrt(np.sum([x**2 for x in vec1]))
            modvec2 = np.sqrt(np.sum([x**2 for x in vec2]))
            costheta = np.dot(vec1, vec2)/(modvec1*modvec2)
            costhetasum += (costheta +(1./3.))**2
        angulars.append(costhetasum)

    system.atoms["angular"] = angulars
    
    mapdict = {}
    mapdict["angular_parameters"] = {}
    mapdict["angular_parameters"]["diamond_angle"] = "angular"
    system.atoms._add_attribute(mapdict)

def calculate_chiparams(system, angles=False):
    """
    Calculate the chi param vector for each atom
    
    Parameters
    ----------
    angles : bool, optional
        If True, return the list of cosines of all neighbor pairs
    
    Returns
    -------
    angles : array of floats
        list of all cosine values, returned only if `angles` is True.
    
    Notes
    -----
    This method tries to distinguish between crystal structures by finding the cosines of angles
    formed by an atom with its neighbors. These cosines are then historgrammed with bins
    `[-1.0, -0.945, -0.915, -0.755, -0.705, -0.195, 0.195, 0.245, 0.795, 1.0]` to find a vector for
    each atom that is indicative of its local coordination. Compared to chi parameters from chi_0 to
    chi_7 in the associated publication, the vector here is from chi_0 to chi_8. This is due to an additional
    chi parameter which measures the number of neighbors between cosines -0.705 to -0.195.
    Parameter `nlimit` specifies the number of nearest neighbors to be included in the analysis to find the cutoff.
    If parameter `angles` is true, an array of all cosine values is returned. The publication further provides
    combinations of chi parameters for structural identification which is not implemented here. The calculated
    chi params can be accessed using :attr:`~pyscal.catom.chiparams`.
    
    References
    ----------
    .. [1] Ackland, Jones, Phys. Rev. B 73, 2006
    """

    system._check_neighbors()

    bins = [-1.0, -0.945, -0.915, -0.755, -0.705, -0.195, 0.195, 0.245, 0.795, 1.0]
    chiparams = []
    cosines = []

    for count, pos in enumerate(system.atoms["positions"]):

        dists = system.atoms["neighbordist"][count]
        neighs = system.atoms["neighbors"][count]

        args = range(len(dists))
        combos = list(itertools.combinations(args, 2))
        costhetas = []
        
        for combo in combos:
            vec1 = system.atoms["diff"][count][combo[0]]
            vec2 = system.atoms["diff"][count][combo[1]]
            modvec1 = np.linalg.norm(vec1)
            modvec2 = np.linalg.norm(vec2)
            costheta = np.dot(vec1, vec2)/(modvec1*modvec2)
            #found costheta
            costhetas.append(costheta)


        #now add according to classification in paper
        chivector = np.histogram(costhetas, bins=bins)
        chiparams.append(chivector[0])
        if angles:
            cosines.append(costhetas)
    
    system.atoms["chiparams"] = chiparams
    
    mapdict = {}
    mapdict["angular_parameters"] = {}
    mapdict["angular_parameters"]["chi_params"] = "chiparams"
    
    if angles:
        system.atoms["cosines"] = cosines
        mapdict["angular_parameters"]["cosines"] = "cosines"

    system.atoms._add_attribute(mapdict)

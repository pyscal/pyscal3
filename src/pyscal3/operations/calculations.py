import numpy as np
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

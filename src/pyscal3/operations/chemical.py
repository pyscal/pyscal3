import numpy as np

def calculate_sro(system, reference_type=1, compare_type=2, average=True):
    """
    Calculate short range order

    Parameters
    ----------
    reference_type: int, optional
        type of the atom to be used a reference. default 1

    compare_type: int, optional
        type of the atom to be compared to. default 2

    average: bool, optional
        if True, average over all atoms of the reference type in the system.
        default True.

    Returns
    -------
    vec: list of float
        The short range order averaged over the whole system for atom of
        the reference type. Only returned if `average` is True. 

    Notes
    -----
    Calculates the short range order for an AB alloy using the approach by
    Cowley [1]. Short range order is calculated as,

    .. math::

        \\alpha_i = 1 - \\frac{n_i}{m_A c_i}

    where n_i is the number of atoms of the non reference type among the c_i atoms
    in the ith shell. m_A is the concentration of the non reference atom. Please
    note that the value is calculated for shells 1 and 2 by default. In order for
    this to be possible, neighbors have to be found first using the :func:`~pyscal.core.System.find_neighbors`
    method. The selected neighbor method should include the second shell as well. For this
    purpose `method=cutoff` can be chosen with a cutoff long enough to include the second
    shell. In order to estimate this cutoff, one can use the :func:`~pyscal.core.System.calculate_rdf`
    method.

    The method also works for alloying elements greater than 2 [2]. For this, you can choose the reference type.

    References
    ----------
    .. [1] Cowley J. M., PR 77(5), 1950.
    .. [2] de Fountaine D., J. Appl. Cryst. 4(15), 1971.

    """

    def _get_dict(val, cdict):
        cx, cxc = np.unique(val, return_counts=True)
        cx = list(cx)
        cxc = list(cxc)
        for key, val in cdict.items():
            if key not in cx:
                cx.append(key)
                cxc.append(0)
        d = dict([(x, cxc[c]/np.sum(cxc)) for c, x in enumerate(cx)])
        return d

    system._check_neighbors()
    
    cdict = system.atoms.composition_ints

    neighbortypes = system.atoms["types"][system.atoms['neighbors']]
    comp_dicts = [_get_dict(neighbortype, cdict) for neighbortype in neighbortypes]
    local_comp = np.array([d[reference_type] for d in comp_dicts])
    global_comp = cdict[reference_type]

    if reference_type == compare_type:
        sro = (local_comp - global_comp)/(1 - global_comp)
    else:
        sro = 1 - (local_comp/global_comp)

    system.atoms["sro"]= sro
    mapdict = {}
    mapdict['chemical'] = {}
    mapdict['chemical']['short_range_order'] = "sro"
    system.atoms._add_attribute(mapdict)

    if average:
        np.mean(sro)
        
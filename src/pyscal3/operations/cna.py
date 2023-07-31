
import pyscal3.csystem as pc

def calculate_cna(system, lattice_constant=None):
    """
    Calculate the Common Neighbor Analysis indices

    Parameters
    ----------
    system: System object

    lattice_constant : float, optional
        lattice constant to calculate CNA. If not specified,
        adaptive CNA will be used

    Returns
    -------
    resdict: dict
        dictionary containing the calculated results

    Notes
    -----
    Performs the common neighbor analysis [1][2] or the adaptive common neighbor
    analysis [2] and assigns a structure to each atom.
    
    If `lattice_constant` is specified, a convential common neighbor analysis is
    used. If `lattice_constant` is not specified, adaptive common neighbor analysis is used. 
    The assigned structures can be accessed by :attr:`~pyscal.catom.Atom.structure`.
    The values assigned for stucture are 0 Unknown, 1 fcc, 2 hcp, 3 bcc, 4 icosahedral.

    References
    ----------
    .. [1] Faken, Jonsson, CMS 2, 1994
    .. [2] Honeycutt, Andersen, JPC 91, 1987
    .. [3] Stukowski, A, Model Simul Mater SC 20, 2012

    """
    system.atoms.create_attribute('structure', fill_with = 0)
    #run to calculate temp neighbors
    system.find_neighbors(method='number', nmax=14, assign_neighbor=False)

    if lattice_constant is None:
        #we need adaptive calculation
        pc.get_acna_neighbors_cn12(system.atoms, system.triclinic, 
            system.rot, system.rotinv, system.boxdims)
        pc.identify_cn12(system.atoms, system.triclinic, 
            system.rot, system.rotinv, system.boxdims)

        pc.get_acna_neighbors_cn14(system.atoms, system.triclinic, 
            system.rot, system.rotinv, system.boxdims)
        pc.identify_cn14(system.atoms, system.triclinic, 
            system.rot, system.rotinv, system.boxdims)

    else:
        self.lattice_constant = lattice_constant
        pc.get_cna_neighbors(system.atoms, system.triclinic, 
            system.rot, system.rotinv, system.boxdims,
            lattice_constant, 1)
        pc.identify_cn12(system.atoms, system.triclinic, 
            system.rot, system.rotinv, system.boxdims)

        pc.get_cna_neighbors(system.atoms, system.triclinic, 
            system.rot, system.rotinv, system.boxdims,
            lattice_constant, 2)
        pc.identify_cn14(system.atoms, system.triclinic, 
            system.rot, system.rotinv, system.boxdims)


    res_dict = {
    "others": len([x for x in system.atoms.structure if x==0]),
    "fcc": len([x for x in system.atoms.structure if x==1]),
    "hcp": len([x for x in system.atoms.structure if x==2]),
    "bcc": len([x for x in system.atoms.structure if x==3]),
    "ico": len([x for x in system.atoms.structure if x==4]),
    }
    return res_dict
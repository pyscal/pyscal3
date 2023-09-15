
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
    system.find.neighbors(method='number', nmax=14, assign_neighbor=False)

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
        system.lattice_constant = lattice_constant
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

def identify_diamond(system):
    """
    Identify diamond structure

    Parameters
    ----------
    system: System object

    Returns
    -------
    diamondstructure : dict
        dict of structure signature

    Notes
    -----
    Identify diamond structure using the algorithm mentioned in [1]. It is an
    extended CNA method. The integers 5, 6, 7, 8, 9 and 10 are assigned to the
    structure variable of the atom. 5 stands for cubic diamond, 6 stands for first
    nearest neighbors of cubic diamond and 7 stands for second nearest neighbors
    of cubic diamond. 8 signifies hexagonal diamond, the first nearest neighbors
    are marked with 9 and second nearest neighbors with 10.

    References
    ----------
    .. [1] Maras et al, CPC 205, 2016
    """
    system.atoms.create_attribute('structure', fill_with = 0)
    #run to calculate temp neighbors
    system.find.neighbors(method='number', nmax=4, assign_neighbor=False)

    pc.identify_diamond_cna(system.atoms, system.triclinic, 
            system.rot, system.rotinv, system.boxdims)

    res_dict = {
    "others": len([x for x in system.atoms.structure if x==0]),
    "cubic diamond": len([x for x in system.atoms.structure if x==5]),
    "cubic diamond 1NN": len([x for x in system.atoms.structure if x==6]),
    "cubic diamond 2NN": len([x for x in system.atoms.structure if x==7]),
    "hex diamond": len([x for x in system.atoms.structure if x==8]),
    "hex diamond 1NN": len([x for x in system.atoms.structure if x==9]),
    "hex diamond 2NN": len([x for x in system.atoms.structure if x==10]),
    }
    return res_dict    
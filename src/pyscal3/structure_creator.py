
import numpy as np
import warnings
import os
import copy

from pyscal3.atoms import Atoms
from pyscal3.attributes import read_yaml
structures = read_yaml(os.path.join(os.path.dirname(__file__), "data/structure_data.yaml"))

def make_crystal(structure, 
    lattice_constant = 1.00, 
    repetitions = None, 
    ca_ratio = 1.633, 
    noise = 0, 
    element=None,
    return_structure_dict=False, 
    structures=structures,
    primitive=False):
    """
    Create a basic crystal structure and return it as a list of `Atom` objects
    and box dimensions.

    Parameters
    ----------
    structure : {'sc', 'bcc', 'fcc', 'hcp', 'diamond', 'a15' or 'l12'}
        type of the crystal structure

    lattice_constant : float, optional
        lattice constant of the crystal structure, default 1

    repetitions : list of ints of len 3, optional
        of type `[nx, ny, nz]`, repetions of the unit cell in x, y and z directions.
        default `[1, 1, 1]`.

    ca_ratio : float, optional
        ratio of c/a for hcp structures, default 1.633

    noise : float, optional
        If provided add normally distributed noise with standard deviation `noise` to the atomic positions.

    Returns
    -------
    atoms : list of `Atom` objects
        list of all atoms as created by user input

    box : list of list of floats
        list of the type `[[xlow, xhigh], [ylow, yhigh], [zlow, zhigh]]` where each of them are the lower
        and upper limits of the simulation box in x, y and z directions respectively.

    Examples
    --------
    >>> atoms, box = make_crystal('bcc', lattice_constant=3.48, repetitions=[2,2,2])
    >>> sys = System()
    >>> sys.assign_atoms(atoms, box)

    """
    if repetitions == None:
        repetitions = [1, 1, 1]
    elif isinstance(repetitions, int):
        repetitions = [repetitions, repetitions, repetitions]

    if structure in structures.keys():
        if primitive:
            if 'primitive' not in structures[structure].keys():
                raise ValueError('primitive not found, try setting primitive=False')
            sdict = copy.copy(structures[structure]['primitive']) 
        else:
            if 'conventional' not in structures[structure].keys():
                raise ValueError('conventional not found, try setting primitive=True')
            sdict = copy.copy(structures[structure]['conventional']) 
    else:
        raise ValueError("Unknown crystal structure")
    
    unique_types = np.unique(sdict["species"])
    #print(element)
    if element is not None:
        if isinstance(element, str):
            element = [element]
        if not (len(element) == len(unique_types)):
            raise ValueError("Elements should equal number of species in the system")
        
        element_dict = dict([x for x in zip(unique_types, element)])
    else:
        element = [None for x in range(len(unique_types))]
    element_dict = dict([x for x in zip(unique_types, element)])
    

    #fix ca ration
    if structure in ['hcp', 'dhcp']:
        sdict['box'][2][2] = ca_ratio 

    #from here, the creation starts
    box = lattice_constant*np.array(sdict["box"])        

    pos = np.array([_unfold_positions(p, box) for p in sdict["positions"]])
    pos = np.array([_generate_noise(x, noise) for x in pos])
    nop = len(pos)
    types = sdict['species']
    species = [element_dict[typ] for typ in types] 

    for d in range(3):        
        pos_list = []
        for i in range(1, repetitions[d]):
            npos = copy.copy(pos)
            npos = npos + i*box[d]
            pos_list.append([_generate_noise(n, noise) for n in npos])

        pos = np.concatenate((pos, *pos_list))
        types = np.concatenate((types, np.tile(types, len(pos_list))))
        species = np.concatenate((species, np.tile(species, len(pos_list))))

                    
    atoms = {}
    atoms['positions'] = pos
    atoms['ids'] = [x+1 for x in range(len(pos))]
    atoms['types'] = types
    atoms['species'] = species
    atoms['ghost'] = [False for x in range(len(types))]

    patoms = Atoms()
    patoms.from_dict(atoms)

    #scale box 
    box[0] = repetitions[0]*np.array(box[0])
    box[1] = repetitions[1]*np.array(box[1])
    box[2] = repetitions[2]*np.array(box[2])

    if return_structure_dict:
        return patoms, box, sdict

    return patoms, box

def general_lattice(positions,
    types, 
    box,
    lattice_constant = 1.00, 
    repetitions = None, 
    noise = 0,
    element=None,
    return_structure_dict=False):
    """
    Create a general lattice structure.

    species: list
        list of per-atom species

    positions:
        list of relative positions positions of reach atom (between 0-1)

    scaling_fractors:
        factors with which the unit cell should be scaled, for example hcp could
        have [1,1.73, 1.63]. Default [1,1,1]

    lattice_constant : float, optional
        lattice constant of the crystal structure, default 1

    repetitions : list of ints of len 3, optional
        of type `[nx, ny, nz]`, repetions of the unit cell in x, y and z directions.
        default `[1, 1, 1]`.

    noise : float, optional
        If provided add normally distributed noise with standard deviation `noise` to the atomic positions.
    
    element : string, optional
        The chemical element
    """
    if not (len(types) == len(positions)):
        raise ValueError("types and positions should have same length!")

    sdict = {"custom":
                {"conventional":
                    {
                     "species": types,
                     "positions": positions,
                     "box": box,
                    }
                }
            }

    atoms, box = make_crystal("custom", lattice_constant=lattice_constant,
        repetitions=repetitions, noise=noise, element=element,
        structures=sdict)

    if return_structure_dict:
        return atoms, box, sdict

    return atoms, box


def _generate_noise(p, noise):
    if noise > 0:
        p[0] = np.random.normal(loc=p[0], scale=noise)
        p[1] = np.random.normal(loc=p[1], scale=noise)
        p[2] = np.random.normal(loc=p[2], scale=noise)
    return p

def _update_list_of_elements():
    """
    Only run when needed to update database
    """
    import mendeleev
    from mendeleev import element
    el_list = dir(mendeleev)
    el_dict = {}

    for el in el_list:
        if len(el) == 2:
            if el=="db":
                break
            chem = element(el)
            struct = chem.lattice_structure
            lc = chem.lattice_constant
            if (struct is not None) and (lc is not None):
                if struct.lower() in ['sc', 'bcc', 'fcc', 'hex', 'dia']:
                    el_dict[el] = {}
                    if struct.lower() == "dia":
                        struct = "diamond"
                    elif struct.lower() == "hex":
                        struct = "hcp"
                    el_dict[el]["structure"] = str(struct.lower())
                    el_dict[el]["lattice_constant"] = float(lc)
    return el_dict

def _unfold_positions(p, box):
    return p[0]*box[0] + p[1]*box[1] + p[2]*box[2]
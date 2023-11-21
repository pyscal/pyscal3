import numpy as np
import copy
from scipy.spatial import distance
import pyscal3.csystem as pc

#from pyscal3.attributes import DocumentedKeywords
def pad_repeat(atoms, box, repetitions, ghost=True):
    """
    This method creates a padded layer around the atoms
    """
    box = np.array(copy.copy(box))        
    idstart = len(atoms) + 1

    pos = atoms.positions
    nop = len(pos)
    ids = atoms.ids
    head = [x for x in range(len(pos))]
    ghosts = [False for x in range(nop)]


    #all the other keys
    datadict = {key: atoms[key][:atoms.nreal] for key in atoms.keys()}
    del datadict['positions']
    del datadict['ids']
    del datadict['head']
    del datadict['ghost']

    #prepopulate
    for d in range(3):
        #first do to the right
        pos_list = []
        new_id_list = []
        ghost_list = []

        for i in range(1, repetitions[d]):
            npos = copy.copy(pos)
            npos = npos + i*box[d]
            pos_list.append(npos)
            new_ids = [idstart+i for i in range(len(pos))]
            new_id_list.append(new_ids)
            ghost_list.append([ghost for x in range(len(pos))])
            idstart = idstart + len(new_ids)

        for i in range(-repetitions[d]+1, 0):
            npos = copy.copy(pos)
            npos = npos + i*box[d]
            pos_list.append(npos)
            new_ids = [idstart+i for i in range(len(pos))]
            new_id_list.append(new_ids)
            ghost_list.append([ghost for x in range(len(pos))])
            #change id start
            idstart = idstart + len(new_ids)

        pos = np.concatenate((pos, *pos_list))
        ids = np.concatenate((ids, *new_id_list))
        ghosts = np.concatenate((ghosts, *ghost_list))
        
        head = np.concatenate((head, np.tile(head, len(pos_list))))

        for key in datadict.keys():
            datadict[key] = np.concatenate((datadict[key], np.tile(datadict[key], len(pos_list))))

    atoms["positions"] = pos
    atoms["ids"] = ids
    atoms["head"] = head
    atoms["ghost"] = ghosts

    for key in datadict.keys():
        atoms[key] = datadict[key]


    box[0] = (2*repetitions[0]-1)*np.array(box[0])
    box[1] = (2*repetitions[1]-1)*np.array(box[1])
    box[2] = (2*repetitions[2]-1)*np.array(box[2])

    return atoms, box


def repeat(system, repetitions, 
    ghost = False, scale_box = True, 
    atoms = None, return_atoms = False):
    """
    Repeat the given system

    Parameters
    ----------
    system: pyscal System object
        the input system to be repeated

    repetitions: tuple of ints
        number of times the system is to be rotated in each direction

    ghost: bool, optional
        if True, make the new atoms ghost, default False

    scale_box: bool, optional
        if True, scale the simulation box, default True

    atoms: None, optional
        if provided use the given atoms, and not the atoms from the system

    return_atoms: bool, optional
        if True, return atoms instead of adding them to the system. Default False
    
    Returns
    -------
    system: pyscal System object
        the system with repetitions. Only returned if `return_atoms` is False.
    
    atoms: Atoms object
        only returned if `return_atoms` is True.


    """
    box = np.array(copy.copy(system.box))        
    system.actual_box = box.copy()

    if atoms is None:
        atoms = system.atoms

    idstart = len(atoms) + 1

    pos = atoms.positions
    nop = len(pos)
    ids = atoms.ids
    head = [x for x in range(len(pos))]
    ghosts = [False for x in range(nop)]

    #all the other keys
    datadict = {key: atoms[key][:atoms.nreal] for key in atoms.keys()}
    del datadict['positions']
    del datadict['ids']
    del datadict['head']
    del datadict['ghost']

    #prepopulate
    for d in range(3):
        pos_list = []
        new_id_list = []
        ghost_list = []
        for i in range(1, repetitions[d]):
            npos = copy.copy(pos)
            npos = npos + i*box[d]
            pos_list.append(npos)
            new_ids = [idstart+i for i in range(len(pos))]
            new_id_list.append(new_ids)
            ghost_list.append([ghost for x in range(len(pos))])
            #change id start
            idstart = idstart + len(new_ids)

        pos = np.concatenate((pos, *pos_list))
        ids = np.concatenate((ids, *new_id_list))
        ghosts = np.concatenate((ghosts, *ghost_list))
        #generate new ids
        
        head = np.concatenate((head, np.tile(head, len(pos_list))))

        for key in datadict.keys():
            datadict[key] = np.concatenate((datadict[key], np.tile(datadict[key], len(pos_list))))

    atoms["positions"] = pos
    atoms["ids"] = ids
    atoms["head"] = head
    atoms["ghost"] = ghosts

    for key in datadict.keys():
        atoms[key] = datadict[key]

    if scale_box:
        box[0] = repetitions[0]*np.array(box[0])
        box[1] = repetitions[1]*np.array(box[1])
        box[2] = repetitions[2]*np.array(box[2])
    if ghost:
        system.ghosts_created = True

    if return_atoms:
        return atoms, box

    else:
        system.box = box
        system.atoms = atoms
        return system


def embed_in_cubic_box(system, input_box=None, 
    return_box=False):
    """
    Embedded the triclinic box in a cubic box
    
    Parameters
    ----------
    system : pyscal System
        input system

    input_box: optional
        if specified, this box is used instead of the
        system simulation box

    return_box: bool, optional
        if specified, return the box instead of assigning
        to the system. Default False

    Returns
    -------
    system: 
        system object with the modified box
        only returned if `return_box` is False

    box:
        the cubic simulation box
        only returned if `return_box` is True
    """
    if input_box is None:
        box = copy.copy(system.box)
        backupbox = copy.copy(box)
    else:
        box = copy.copy(input_box)
        backupbox = copy.copy(input_box)

    a = np.array(box[0])
    b = np.array(box[1])
    c = np.array(box[2])

    cosa = np.dot(b, c)/(np.linalg.norm(b)*np.linalg.norm(c))
    cosb = np.dot(c, a)/(np.linalg.norm(c)*np.linalg.norm(a))
    cosc = np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b))

    lx = np.linalg.norm(a)
    xy = np.linalg.norm(b)*cosc
    xz = np.linalg.norm(c)*cosb
    ly = np.sqrt(np.linalg.norm(b)*np.linalg.norm(b) - xy*xy)
    yz = (np.linalg.norm(b)*np.linalg.norm(c)*cosa - xy*xz)/ly
    lz = np.sqrt(np.linalg.norm(c)*np.linalg.norm(c) - xz*xz - yz*yz)

    xlo = ylo = zlo = 0
    xhi = lx
    yhi = ly
    zhi = lz

    xlo_bound = xlo + min(0.0,xy,xz,xy+xz)
    xhi_bound = xhi + max(0.0,xy,xz,xy+xz)
    ylo_bound = ylo + min(0.0,yz)
    yhi_bound = yhi + max(0.0,yz)
    zlo_bound = zlo
    zhi_bound = zhi

    newbox = np.array([[xhi_bound-xlo_bound, 0, 0], [0, yhi_bound-ylo_bound, 0], [0, 0, zhi_bound-zlo_bound]])
    
    if not return_box:
        system.newbox = newbox
        system.box = newbox
        system.box_backup = backupbox
        system.actual_box = None
        return system
    else:
        return newbox


def extract_cubic_representation(system, repetitions = (3,3,3), return_atoms = False):
    """
    Extract a cubic representation of a given box.

    Parameters
    ----------
    system: pyscal System object
        the input system to be repeated

    repetitions: tuple of ints
        number of times the system is to be rotated in each direction

    return_atoms: bool, optional
        if True, return atoms instead of adding them to the system. Default False
    
    Returns
    -------
    system: pyscal System object
        the system with repetitions. Only returned if `return_atoms` is False.
    
    atoms: Atoms object
        only returned if `return_atoms` is True.
    

    Notes
    -----
    This algorithm is not guaranteed to work. Althought it can extract a cubic box,
    there is no guarantee that the extracted box is actually a unit cell.
    """
    atoms, box = pad_repeat(system.atoms, system.box, repetitions, ghost=False)
    
    def _is_in_bound(pos, lo, hi):
        c1 = np.where(pos >= lo, True, False)
        c2 = np.where(pos < hi, True, False)
        return c1*c2

    bounds = []
    seed = 0

    for search_dir, plane in enumerate([[1, 2], [0, 2], [0, 1]]):
        anchor_type = atoms.types[seed]
        anchor = atoms.positions[seed]

        plane_distance_1 = np.abs(atoms.positions[:,plane[0]] - anchor[plane[0]])
        plane_distance_2 = np.abs(atoms.positions[:,plane[1]] - anchor[plane[1]])
        indices_1 = [count for count, p in enumerate(plane_distance_1) if p < 1e-5]
        indices_2 = [count for count, p in enumerate(plane_distance_2) if p < 1e-5]
        indices = list(set(indices_1).intersection(indices_2))
        indices = [index for index in indices if atoms.types[index] == anchor_type]
        distances = distance.cdist([anchor], atoms.positions[indices])[0]
        sorted_args = np.argsort(distances)
        distances = np.array(distances)[sorted_args][1:]
        indices = np.array(indices)[sorted_args][1:]

        if len(indices) == 0:
            raise ValueError(f'Could not find cubic representation, please increase repetitions and try!')
        closest_atom = atoms.positions[indices][0]
        bounds.append([anchor[search_dir], closest_atom[search_dir]])

    #now filter out all atoms which fall in the bounds
    x_atoms = _is_in_bound(atoms.positions[:,0], min(bounds[0]), max(bounds[0]))
    y_atoms = _is_in_bound(atoms.positions[:,1], min(bounds[1]), max(bounds[1]))
    z_atoms = _is_in_bound(atoms.positions[:,2], min(bounds[2]), max(bounds[2]))
    is_needed = x_atoms*y_atoms*z_atoms

    #now get those atoms
    #delete unneeded ones
    unneeded = [count for count in range(len(atoms.positions)) if not is_needed[count]]
    atoms.delete(indices=unneeded)
    box = []
    box.append([max(bounds[0])-min(bounds[0]), 0, 0])
    box.append([0, max(bounds[1])-min(bounds[1]), 0])
    box.append([0, 0, max(bounds[2])-min(bounds[2])])
    
    if return_atoms:
        return atoms, box
    else:
        system.box = box
        system.atoms = atoms
        system.modify.remap_to_box()
        return system


def remap_to_box(system, ghosts=True):
    """
    Remap the atom to back inside the box.
    
    Parameters
    ----------
    None

    Returns
    -------
    system
    """
    #rot = np.array(system._box).T
    #rotinv = np.linalg.inv(rot)
    if ghosts:
        for x in range(system.atoms.ntotal):
            pos = pc.remap_atom_into_box(system.atoms["positions"][x], 
                system.triclinic,
                system.rot, 
                system.rotinv, 
                system.boxdims)
            #print(f'{system.atoms["positions"][x]} changed to {pos} ')
            system.atoms["positions"][x] = pos
    else:
        box = system.box
        rot = np.array(box).T
        rotinv = np.linalg.inv(rot)
        boxdims = [0,0,0]
        boxdims[0] = np.sum(np.array(box[0])**2)**0.5
        boxdims[1] = np.sum(np.array(box[1])**2)**0.5
        boxdims[2] = np.sum(np.array(box[2])**2)**0.5

        for x in range(system.atoms.natoms):
            pos = pc.remap_atom_into_box(system.atoms["positions"][x], 
                system.triclinic,
                rot, 
                rotinv, 
                boxdims)
            #print(f'{system.atoms["positions"][x]} changed to {pos} ')
            system.atoms["positions"][x] = pos

    return system

def remap_position_to_box(system, pos):
    pos = pc.remap_atom_into_box(pos, 
        system.triclinic,
        system.rot, 
        system.rotinv, 
        system.boxdims)
    #print(f'{system.atoms["positions"][x]} changed to {pos} ')
    return pos


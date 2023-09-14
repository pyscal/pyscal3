
import numpy as np
import warnings
import os


def delete(atoms, ids=None, indices=None, condition=None, selection=False):
    """
    Delete atoms by a given property

    Parameters
    ----------
    atoms: Atoms object

    ids: list of ints
        list of ids to be deleted

    indices: list of ints
        list of indices to be deleted

    condition: callable or list of bool

    selection: bool, optional
        can be used to delete selected atoms
    """
    #generate a list of atoms to be deleted based on the needed conditions
    masks = _generate_bool_list(atoms, 
        ids=ids, 
        indices=indices, 
        condition=condition, 
        selection=selection)

    #get ids of atoms to be deleted
    delete_list = [masks[atoms["head"][x]] for x in range(atoms.ntotal)]
    delete_ids = [x for x in range(atoms.ntotal) if masks[x]]
    
    del_real = np.sum([1 for x in indices if x < atoms._nreal])
    del_ghost = np.sum([1 for x in indices if x >= atoms._nreal])

    #find atoms that need not be deleted; and trim all properties to meet the length
    #this is equivalent to deletion
    reverse_indices = [x for x in range(len(atoms['positions'])) if x not in indices]
    for key in atoms.keys():
        atoms[key] = atoms[key][reverse_indices]

    #update lengths
    td = len(indices)
    atoms._nreal = int(atoms.nreal - del_real)
    atoms._nghost = int(atoms.nghost - del_ghost)
    #atoms is dict, no need to return

def _generate_bool_list(atoms, ids=None, indices=None, condition=None, selection=False):
    #necessary checks
    non_nones = sum(x is not None for x in [ids, indices, condition])
    if non_nones > 1:
        raise ValueError("Only one of ids, indices or condition should be provided")
    elif ((non_nones == 0) and (selection==False)):
        warnings.warn("No conditions provided, all atoms will be included")
    #generate a list of indices
    if selection:
        indices = [x for x in range(atoms.nreal) if atoms["condition"][x]]
    elif ids is not None:
        if np.isscalar(ids):
            ids = [ids]
        indices = [x for x in range(len(atoms["ids"])) if atoms["ids"][x] in ids]
        if len(indices) == 0:
            raise ValueError("No ids found to delete")
        if len(indices) != len(ids):
            warnings.warn("Not all ids were found")
    elif condition is not None:
        indices = [x for x in range(atoms.nreal) if condition(atoms._get_atoms(x))]
    elif indices is None:
        indices = [x for x in range(atoms.nreal)]
    
    if np.isscalar(indices):
        indices = [indices]

    bool_list = [ True if x in indices else False for x in range(self.nreal)]
    return bool_list

def _validate_condition(atoms, condition):
    if not (len(condition)==atoms.nreal):
        raise ValueError("condition should have same length as atoms")
    for c, x in enumerate(condition):
        try:
            x = bool(x)
            condition[c] = x
        except:
            pass 
        if not isinstance(x, bool):
            raise ValueError("Condition elements should be boolean")
    return condition

def apply_selection(atoms, ids=None, indices=None, condition=None):
    if isinstance(condition, (list, np.ndarray)):
        masks = atoms._validate_condition(condition)
    else:
        masks = atoms._generate_bool_list(ids=ids, indices=indices, condition=condition)
    for i in range(atoms.ntotal):
        atoms["condition"][i] = masks[atoms["head"][i]]

def remove_selection(atoms, ids=None, indices=None, condition=None):
    if isinstance(condition, (list, np.ndarray)):
        masks = atoms._validate_condition(condition)
    else:
        masks = atoms._generate_bool_list(ids=ids, indices=indices, condition=condition)
    masks = [not x for x in masks]
    for i in range(atoms.ntotal):
        atoms["condition"][i] = masks[atoms["head"][i]]

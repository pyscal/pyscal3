"""
This file includes the definition of an Atoms class which can be used with
System

TODO
----
- Iterators for masked/unmasked atoms
- Iterators for selected/unselected atoms
"""

import numpy as np
import warnings
import os
from pyscal3.attributes import AttrSetter, read_yaml, MyList

class Atoms(dict, AttrSetter):
    def __init__(self, atoms=None):
        #self.update(atoms=atoms)
        self._nreal = 0
        self._nghost = 0
        self._lattice_constant = None
        self._lattice = None
        AttrSetter.__init__(self)
        if atoms is not None:
            self.from_dict(atoms)
    
    def __dir__(self):
        attrs = ["natoms", "nreal", "nghost", 
        "ntotal", "from_dict", "iter_atoms", "apply_mask", "remove_mask",
        "apply_selection", "remove_selection", "delete", "composition"]
        return attrs + list(self._map_dict.keys())

    def __getitem__(self, key):
        if isinstance(key, slice):
            return self._get_atoms(key)
        elif isinstance(key, int):
            return self._get_atoms(key)
        else:
            val = dict.__getitem__(self, key)
            return val

    def __setitem__(self, key, val):
        dict.__setitem__(self, key, np.array(val))

    def __len__(self):
        return self.nreal

    def __repr__(self):
        dictrepr = dict.__repr__(self)
        return '%s(%s)' % (type(self).__name__, dictrepr)

    def _repr_json(self):
        #convert to atom base dict
        disp_atoms = {f"atom {x}": self._get_atoms(x) for x in range(self.natoms)}
        return disp_atoms
        
    def __add__(self, atoms):
        if not 'positions' in atoms.keys():
            raise ValueError('positions is a necessary key in atoms')
        nop = len(atoms["positions"])
        val_length_check = np.prod([len(val)==nop for key, val in atoms.items()])
        if not val_length_check:
            raise ValueError("All times in the atoms dict should have same length as positions")
        
        #now add necessary keys-ids, types, ghost
        maxid = max(self["ids"])
        if not 'ids' in atoms.keys():
            atoms['ids'] = [maxid+x+1 for x in range(nop)]
            #print(self["ids"], atoms["ids"])
        else:
            if len(set(atoms['ids']).intersection(set(self['ids']))):
                raise ValueError("Atom id already exists, unique ID is required")
        
        atoms['ghost'] = [False for x in range(nop)]
        if not 'types' in atoms.keys():
            atoms['types'] = [1 for x in range(nop)]
        if not 'species' in atoms.keys():
            atoms['species'] = [None for x in range(nop)]
        if not 'mask_1' in atoms.keys():
            atoms['mask_1'] = [False for x in range(nop)]
        if not 'mask_2' in atoms.keys():
            atoms['mask_2'] = [False for x in range(nop)]
        if not 'condition' in atoms.keys():
            atoms['condition'] = [True for x in range(nop)]
        if not 'head' in atoms.keys():
            atoms['head'] = [self.natoms+x for x in range(nop)]
        
        common_keys = list(set(self.keys()).intersection(set(atoms.keys())))
        for key in common_keys:
            self[key] = np.concatenate((self[key], atoms[key]))
        extra_keys_add = len(atoms.keys()) - len(common_keys)
        extra_keys_exist = len(self.keys()) - len(common_keys)
        
        if extra_keys_add > 0:
            warnings.warn("Some keys present in the atoms are add are not present in existing atoms, they were ignored")
        if extra_keys_exist > 0:
            warnings.warn("Some keys present in the existing Atoms were not present in the atoms to add, please recalculate")
        self._nreal += nop
        return self


    def __radd__(self, atoms):
        """
        Reverse add method
        """
        if ntraj == 0:
            return self
        else:
            return self.__add__(atoms)

    #def _add_attribute()
      
    @property
    def natoms(self):
        return len([1 for x in self['ghost'] if x==False])

    @property
    def nreal(self):
        print(self.keys())
        return len([1 for x in self['ghost'] if x==False])
    
    @property
    def nghost(self):
        return len([1 for x in self['ghost'] if x==True])
    
    @property
    def ntotal(self):
        return len(self['positions'])

    def create_attribute(self, key, fill_with=None, alias=None):
        """
        Create an attribute in atoms, and will with given value
        """
        if alias is None:
            alias = key
        self[key] = np.array([fill_with for x in range(self.ntotal)])
        mapdict = {alias: key}
        self._add_attribute(mapdict)

    def from_dict(self, atoms):
        if not 'positions' in atoms.keys():
            raise ValueError('positions is a necessary key in atoms')
        nop = len(atoms["positions"])
        
        if not 'ids' in atoms.keys():
            atoms['ids'] = np.array([x+1 for x in range(nop)])
        
        atoms['ghost'] = np.array([False for x in range(nop)])
        if not 'types' in atoms.keys():
            atoms['types'] = np.array([1 for x in range(nop)])
        if not 'species' in atoms.keys():
            atoms['species'] = np.array([None for x in range(nop)])
        if not 'mask_1' in atoms.keys():
            atoms['mask_1'] = np.array([False for x in range(nop)])
        if not 'mask_2' in atoms.keys():
            atoms['mask_2'] = np.array([False for x in range(nop)])
        if not 'condition' in atoms.keys():
            atoms['condition'] = np.array([True for x in range(nop)])
        if not 'head' in atoms.keys():
            atoms['head'] = np.array([x for x in range(nop)])
        
        for key, val in atoms.items():
            self[key] = np.array(val)
        self._nreal = len(val)

        #add attributes
        mapdict = {"positions": "positions",
        "ids": "ids",
        "types": "types",
        "species": "species",
        "mask": {"primary": "mask_1", "secondary": "mask_2"},
        "selection": "condition",
        "condition": "condition",
        "head": "head"}

        #add extra keys that might be needed; non-standard ones
        for key, val in atoms.items():
            if key not in ["positions", "ids", "types", "species", "mask_1", "mask_2", "condition", "head"]:
                mapdict[key] = key

        self._add_attribute(mapdict)

    def _convert_to_list(self, data):
        """
        Check if the given item is a list, if not convert to a single item list
        """
        if np.isscalar(data):
            data = np.array([data])
        return data
       
    def _get_atoms(self, index):
        atom_dict = {}
        for key in self.keys():
            if index < len(self[key]):
                atom_dict[key] = self._convert_to_list(self[key][index]) 
        #atom_dict = {key: self._convert_to_list(self[key][index]) for key in self.keys()}
        return Atoms(atom_dict)

    def iter_atoms(self):
        for index in range(self.nreal):
            atom_dict = {}
            for key in self.keys():
                if index < len(self[key]):
                    atom_dict[key] = self._convert_to_list(self[key][index]) 
            yield Atoms(atom_dict)

    #def _apply_mask(self, masks, mask_type):
    #    if (mask_type == 'primary') or (mask_type == 'all'):
    #        for i in range(self.ntotal):
    #            self["mask_1"][i] = masks[self["head"][i]]
    #    if (mask_type == 'secondary') or (mask_type == 'all'):
    #        for i in range(self.ntotal):
    #            self["mask_2"][i] = masks[self["head"][i]]

    #def apply_mask(self, mask_type="primary", ids=None, indices=None, condition=None, selection=False):
    #    masks = self._generate_bool_list(ids=ids, indices=indices, condition=condition, selection=selection)
    #    self._apply_mask(masks, mask_type)


    #def remove_mask(self, mask_type="primary", ids=None, indices=None, condition=None, selection=False):
    #    masks = self._generate_bool_list(ids=ids, indices=indices, condition=condition, selection=selection)
    #    masks = [not x for x in masks]
    #    self._apply_mask(masks, mask_type)


    
    @property
    def _type_dict(self):
        sp = []
        types, typecounts = np.unique(self["types"][:self.nreal], return_counts=True)
        for t in types:
            for count, tx in enumerate(self["types"][:self.nreal]):
                if t==tx:
                    sp.append(self["species"][count])
                    break
        return dict([x for x in zip(types, sp)])


    @property
    def composition(self):
        if self["species"][0] is None:
            typelist = self["types"][:self.nreal]
            types, typecounts = np.unique(typelist, return_counts=True)
            concdict = dict([(t, typecounts[c]/np.sum(typecounts)) for c, t in enumerate(types)])
        else:
            typelist = self["species"][:self.nreal]
            types, typecounts = np.unique(typelist, return_counts=True)
            concdict = {str(t): typecounts[c]/np.sum(typecounts) for c, t in enumerate(types)}
        return concdict




    
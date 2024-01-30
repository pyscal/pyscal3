"""


"""


import os
import numpy as np
import warnings
import itertools
from ase.io import write
import uuid
import gzip
import io
from scipy.special import sph_harm
import copy
from functools import partial, update_wrapper

from pyscal3.attributes import read_yaml
from pyscal3.atoms import Atoms, AttrSetter
import pyscal3.csystem as pc
import pyscal3.traj_process as ptp
from pyscal3.formats.ase import convert_snap
import pyscal3.structure_creator as pcs
from pyscal3.grain_boundary import GrainBoundary

import pyscal3.operations.operations as operations
import pyscal3.operations.cna as cna
import pyscal3.operations.centrosymmetry as centrosymmetry
import pyscal3.operations.neighbor as neighbor
import pyscal3.operations.input as inputmethods
import pyscal3.operations.calculations as calculations
import pyscal3.operations.identify as identify
import pyscal3.operations.voronoi as voronoi
import pyscal3.operations.chemical as chemical
import pyscal3.operations.visualize as visualize
import pyscal3.operations.serialize as serialize

#import pyscal.routines as routines
#import pyscal.visualization as pv

structure_dict = read_yaml(os.path.join(os.path.dirname(__file__), "data/structure_data.yaml"))
element_dict = read_yaml(os.path.join(os.path.dirname(__file__), "data/element_data.yaml"))

def _make_crystal(structure, 
    lattice_constant = 1.00, 
    repetitions = None, 
    ca_ratio = 1.633, 
    noise = 0, 
    element=None,
    primitive=False):
    
    atoms, box, sdict = pcs.make_crystal(structure, 
        lattice_constant=lattice_constant,
        repetitions=repetitions, 
        ca_ratio=ca_ratio,
        noise=noise, 
        element=element, 
        return_structure_dict=True,
        primitive=primitive)
    
    s = System()
    s.box = box
    s.atoms = atoms
    s.atoms._lattice = structure
    s.atoms._lattice_constant = lattice_constant
    s._structure_dict = sdict
    return s

def _make_general_lattice(positions,
    types, 
    box,
    lattice_constant = 1.00, 
    repetitions = None, 
    noise = 0,
    element=None):

    atoms, box, sdict = pcs.general_lattice(positions,
        types,
        box,
        lattice_constant=lattice_constant,
        repetitions=repetitions,
        noise=noise,
        element=element,
        return_structure_dict=True)
    s = System()
    s.box = box
    s.atoms = atoms
    s.atoms._lattice = 'custom'
    s.atoms._lattice_constant = lattice_constant
    s._structure_dict = sdict
    return s

def _make_grain_boundary(axis, 
    sigma, gb_plane,
    structure = None,
    element = None, 
    lattice_constant = 1,
    repetitions = (1,1,1),
    overlap=0.0):

    gb = GrainBoundary()
    gb.create_grain_boundary(axis=axis, sigma=sigma, 
                             gb_plane=gb_plane)

    if structure is not None:
        atoms, box, sdict = gb.populate_grain_boundary(structure, 
                                        repetitions = repetitions,
                                        lattice_parameter = lattice_constant,
                                        overlap=overlap)
    elif element is not None:
        atoms, box, sdict = gb.populate_grain_boundary(element, 
                                        repetitions=repetitions,
                                        overlap=overlap)
    s = System()
    s.box = box
    s.atoms = atoms
    s.atoms._lattice = structure
    s.atoms._lattice_constant = lattice_constant
    s._structure_dict = sdict
    #s = operations.remap_to_box(s)
    return s


class System:
    """
    Python class for holding the properties of an atomic configuration 
    """ 
    #system wide things available before structure creation
    create = AttrSetter()
    #create.head = pcs
    mapdict = {}
    mapdict["lattice"] = {}
    for key in structure_dict.keys():
        mapdict["lattice"][key] = update_wrapper(partial(_make_crystal, key), 
            _make_crystal)
    mapdict["lattice"]["custom"] = _make_general_lattice

    mapdict["element"] = {}
    for key in element_dict.keys():
        mapdict["element"][key] = update_wrapper(partial(_make_crystal,
            element_dict[key]['structure'],
            lattice_constant=element_dict[key]['lattice_constant'],
            element = key), pcs.make_crystal)

    mapdict["defect"] = {}
    mapdict["defect"]["grain_boundary"] = _make_grain_boundary

    create._add_attribute(mapdict)

    def __init__(self, filename=None, format="lammps-dump", 
                                            compressed = False, customkeys=None):
        self.initialized = True
        self.neighbors_found = False
        self.neighbor_method = None
        self.ghosts_created = False
        self.actual_box = None

        #box parameters
        self.rot = [[0,0,0], [0,0,0], [0,0,0]]
        self.rotinv = [[0,0,0], [0,0,0], [0,0,0]]
        self.boxdims = [0,0,0]
        self.triclinic = 0
        self._structure_dict = None
        self._atoms = Atoms()
        
        if filename is not None:
            inputmethods.read_inputfile(self, filename, 
                format=format, 
                compressed = compressed, 
                customkeys=customkeys)

        #customised methods for the class
        self.modify = AttrSetter()
        self.modify.head = operations
        mapdict = {}
        #repeat methid
        mapdict["repeat"] = update_wrapper(partial(operations.repeat, self), operations.repeat)
        mapdict["transform_to_cubic_cell"] = update_wrapper(partial(operations.extract_cubic_representation, self), operations.extract_cubic_representation)
        mapdict["remap_to_box"] = update_wrapper(partial(operations.remap_to_box, self), operations.remap_to_box)
        mapdict["remap_position_to_box"] = update_wrapper(partial(operations.remap_position_to_box, self), operations.remap_position_to_box)
        mapdict["embed_in_cubic_box"] = update_wrapper(partial(operations.embed_in_cubic_box, self), operations.embed_in_cubic_box)
        self.modify._add_attribute(mapdict)

        self.find = AttrSetter()
        mapdict = {}
        mapdict['neighbors'] = update_wrapper(partial(identify.find_neighbors, self), identify.find_neighbors)
        mapdict['solids'] = update_wrapper(partial(identify.find_solids, self), identify.find_solids)
        mapdict['clusters'] = update_wrapper(partial(identify.cluster_atoms, self), identify.cluster_atoms)
        mapdict['largest_cluster'] = update_wrapper(partial(identify.find_largest_cluster, self), identify.find_largest_cluster)
        self.find._add_attribute(mapdict)

        self.calculate = AttrSetter()
        mapdict = {}
        mapdict['distance'] = update_wrapper(partial(neighbor.get_distance, self), neighbor.get_distance)
        mapdict['average_over_neighbors'] = update_wrapper(partial(calculations.average_over_neighbors, self), calculations.average_over_neighbors)
        mapdict['steinhardt_parameter'] = update_wrapper(partial(calculations.calculate_q, self), calculations.calculate_q)
        mapdict['disorder'] = update_wrapper(partial(calculations.calculate_disorder, self), calculations.calculate_disorder)        
        mapdict['radial_distribution_function'] = update_wrapper(partial(calculations.calculate_rdf, self), calculations.calculate_rdf)        
        mapdict['angular_criteria'] = update_wrapper(partial(calculations.calculate_angularcriteria, self), calculations.calculate_angularcriteria)        
        mapdict['chi_params'] = update_wrapper(partial(calculations.calculate_chiparams, self), calculations.calculate_chiparams)
        mapdict['common_neighbor_analysis'] = update_wrapper(partial(cna.calculate_cna, self), cna.calculate_cna)
        mapdict['diamond_structure'] = update_wrapper(partial(cna.identify_diamond, self), cna.identify_diamond)
        mapdict['centrosymmetry'] = update_wrapper(partial(centrosymmetry.calculate_centrosymmetry, self), centrosymmetry.calculate_centrosymmetry)
        mapdict['voronoi_vector'] = update_wrapper(partial(voronoi.calculate_vorovector, self), voronoi.calculate_vorovector)
        self.calculate._add_attribute(mapdict)

        self.analyze = AttrSetter()
        mapdict = {}
        mapdict['common_neighbor_analysis'] = update_wrapper(partial(cna.calculate_cna, self), cna.calculate_cna)
        mapdict['diamond_structure'] = update_wrapper(partial(cna.identify_diamond, self), cna.identify_diamond)
        mapdict['short_range_order'] = update_wrapper(partial(chemical.calculate_sro, self), chemical.calculate_sro)
        self.analyze._add_attribute(mapdict)

        self.chemical = AttrSetter()
        mapdict = {}
        mapdict['short_range_order'] = update_wrapper(partial(chemical.calculate_sro, self), chemical.calculate_sro)
        self.chemical._add_attribute(mapdict)

        self.read = AttrSetter()
        mapdict = {}
        mapdict['file'] = update_wrapper(partial(inputmethods.read_inputfile, self), inputmethods.read_inputfile)
        mapdict['ase'] = update_wrapper(partial(inputmethods.read_inputfile, self, format='ase'), inputmethods.read_inputfile)
        self.read._add_attribute(mapdict)

        self.convert_to = AttrSetter()
        mapdict = {}
        mapdict['ase'] = update_wrapper(partial(convert_snap, self), convert_snap)
        self.convert_to._add_attribute(mapdict)

        self.write = AttrSetter()
        mapdict = {}
        mapdict['ase'] = update_wrapper(partial(convert_snap, self), convert_snap)
        mapdict['file'] = update_wrapper(partial(inputmethods.to_file, self), inputmethods.to_file)
        mapdict['dict'] = update_wrapper(partial(serialize.serialize, self, return_type='dict'), serialize.serialize)
        mapdict['json'] = update_wrapper(partial(serialize.serialize, self, return_type='json'), serialize.serialize)
        mapdict['pydantic'] = update_wrapper(partial(serialize.serialize, self, return_type='model'), serialize.serialize)
        mapdict['json_file'] = update_wrapper(partial(serialize.serialize, self, return_type='file'), serialize.serialize)
        self.write._add_attribute(mapdict)

        self.show = AttrSetter()
        mapdict = {}
        mapdict['all'] = update_wrapper(partial(visualize.plot_simple, self), visualize.plot_simple)
        mapdict['continuous_property'] = update_wrapper(partial(visualize.plot_by_property, self), visualize.plot_by_property)
        mapdict['boolean_property'] = update_wrapper(partial(visualize.plot_by_boolean, self), visualize.plot_by_boolean)
        mapdict['selection'] = update_wrapper(partial(visualize.plot_by_selection, self), visualize.plot_by_selection)
        self.show._add_attribute(mapdict)


    def iter_atoms(self):
        return self.atoms.iter_atoms()

    @property
    def natoms(self):
        return self.atoms.natoms

    @property
    def concentration(self):
        return self.atoms.composition

    @property
    def composition(self):
        return self.atoms.composition

    @property
    def box(self):
        """
        Wrap for inbuilt box
        """
        if self.actual_box is not None:
            return self.actual_box
        else:
            return self._box

    @box.setter
    def box(self, userbox):
        """
        Box setter 
        """
        #we should automatically check for triclinic cells here
        self.internal_box = userbox
        self.actual_box = userbox

    @property
    def internal_box(self):
        return self._box

    @internal_box.setter
    def internal_box(self, userbox):
        """
        Box setter 
        """
        #we should automatically check for triclinic cells here
        summ = 0
        for i in range(3):
            box1 = np.array(userbox[i-1])
            box2 = np.array(userbox[i])
            summ += np.dot(box1, box2)/(np.linalg.norm(box1)*np.linalg.norm(box2))

        #check if the summ is zero
        if (np.abs(summ) > 1E-6):
            #this is a triclinic box
            rot = np.array(userbox).T
            rotinv = np.linalg.inv(rot)
            self.triclinic = 1
            self.rot = rot
            self.rotinv = rotinv

        self.boxdims[0] = np.sum(np.array(userbox[0])**2)**0.5
        self.boxdims[1] = np.sum(np.array(userbox[1])**2)**0.5
        self.boxdims[2] = np.sum(np.array(userbox[2])**2)**0.5

        #and we reset the original memory of the box
        self._box = userbox

    @property
    def box_dimensions(self):
        return [np.linalg.norm(self.box[x]) for x in range(3)]

    @property
    def direct_coordinates(self):
        dims = self.box_dimensions
        xfrac = np.array(self.atoms.positions)[:,0]/dims[0]
        yfrac = np.array(self.atoms.positions)[:,1]/dims[1]
        zfrac = np.array(self.atoms.positions)[:,2]/dims[2]
        coords = np.column_stack((xfrac, yfrac, zfrac))
        return coords    
    
    @property
    def volume(self):
        """
        Volume of box
        """
        vol = np.dot(np.cross(self._box[0], self._box[1]), self._box[2])
        return vol

    @property
    def lattice_properties(self):
        """
        Return lattice properties
        """
        if self._structure_dict is not None:
            ldict = copy.copy(self._structure_dict)
        else:
            ldict = {}
        if self.atoms._lattice is not None:
            ldict['lattice'] = self.atoms._lattice 
        if self.atoms._lattice_constant is not None:
            ldict['lattice_constant'] = self.atoms._lattice_constant 
        return ldict

    @lattice_properties.setter
    def lattice_properties(self, ldict):
        if self._structure_dict is None:
            self._structure_dict = {}

        for key, val in ldict.items():
            if key in ['species', 'box', 'positions']:
                self._structure_dict[key] = val
            elif key == 'lattice':
                self.atoms._lattice = val
            elif key == 'lattice_constant':
                self.atoms._lattice_constant = val

    @property
    def atoms(self):
        return self._atoms
    
    @atoms.setter
    def atoms(self, atoms):
        """
        Set atoms
        """

        if(len(atoms['positions']) < 200):
            #we need to estimate a rough idea
            needed_atoms = 200 - len(atoms['positions'])
            #get a rough cell
            #print(needed_atoms)
            needed_cells = np.ceil(needed_atoms/len(atoms['positions']))
            nx = int(needed_cells**(1/3))
            if nx < 2:
                nx = 2
            if np.sum(self.box) == 0:
                raise ValueError("Simulation box should be initialized before atoms")
            atoms, box = operations.repeat(self, (nx, nx, nx), atoms=atoms, ghost=True, return_atoms=True)
            self.actual_box = self.box.copy()
            self.internal_box = box

        self._atoms = atoms
        


    def add_atoms(self, atoms):
        """
        Cleanly add a given list of atoms

        Parameters
        ----------
        atoms : dict

        Returns
        -------
        None
        """ 
        ## MOVE TO ATOMS
        self._atoms.add_atoms(atoms)
     
    def apply_selection(self, ids=None, indices=None, condition=None):
        self._atoms.apply_selection(ids=ids, indices=indices, condition=condition)    
    
    def remove_selection(self, ids=None, indices=None, condition=None):
        self._atoms.remove_selection(ids=ids, indices=indices, condition=condition)
    
    def delete(self, ids=None, indices=None, condition=None, selection=False):
        self._atoms.delete(ids=ids, indices=indices, condition=condition, selection=selection)



    def reset_neighbors(self):
        """
        Reset the neighbors of all atoms in the system.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        It is used automatically when neighbors are recalculated.

        """
        neighbor.reset_neighbors(self)


    def _check_neighbors(self):
        """
        Check if neighbors are calculated
        """
        if not self.neighbors_found:
            raise ValueError("This calculation needs neighbors to be calculated")


    
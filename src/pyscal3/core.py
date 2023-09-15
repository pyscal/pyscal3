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

import pyscal3.operations.operations as operations
import pyscal3.operations.cna as cna
import pyscal3.operations.centrosymmetry
import pyscal3.operations.neighbor as neighbor
import pyscal3.operations.input as inputmethods
import pyscal3.operations.calculations as calculations
import pyscal3.operations.identify as identify
#import pyscal.routines as routines
#import pyscal.visualization as pv

structure_dict = read_yaml(os.path.join(os.path.dirname(__file__), "data/structure_data.yaml"))
element_dict = read_yaml(os.path.join(os.path.dirname(__file__), "data/element_data.yaml"))

def _make_crystal(structure, 
    lattice_constant = 1.00, 
    repetitions = None, 
    ca_ratio = 1.633, 
    noise = 0, 
    element=None):
    
    atoms, box, sdict = pcs.make_crystal(structure, 
        lattice_constant=lattice_constant,
        repetitions=repetitions, 
        ca_ratio=ca_ratio,
        noise=noise, 
        element=element, 
        return_structure_dict=True)
    s = System()
    s.box = box
    s.atoms = atoms
    s.atoms._lattice = structure
    s.atoms._lattice_constant = lattice_constant
    s._structure_dict = sdict
    return s

def _make_general_lattice(positions,
    types, 
    scaling_factors=[1.0, 1.0, 1.0],
    lattice_constant = 1.00, 
    repetitions = None, 
    noise = 0,
    element=None):

    atoms, box, sdict = pcs.general_lattice(positions,
        types,
        scaling_factors=scaling_factors,
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
        self.calculate._add_attribute(mapdict)

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
        self.write._add_attribute(mapdict)

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


    def calculate_rdf(self, rmin=0, rmax=5.00, bins=100):
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
        self.find.neighbors(method="cutoff", cutoff=rmax)
        distances = list(itertools.chain(*self.atoms["neighbordist"]))

        hist, bin_edges = np.histogram(distances, bins=bins, 
            range=(rmin, rmax), density=True)
        edgewidth = np.abs(bin_edges[1]-bin_edges[0])
        hist = hist.astype(float)
        r = bin_edges[:-1]

        #get box density
        rho = self.natoms/self.volume
        shell_vols = (4./3.)*np.pi*((r+edgewidth)**3 - r**3)
        shell_rho = hist/shell_vols
        #now divide to get final value
        rdf = shell_rho/rho
        return rdf, r

    def calculate_angularcriteria(self):
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
        self._check_neighbors()
        angulars = []

        for count, pos1 in enumerate(self.atoms["positions"]):
            
            dists = []
            distneighs = []
            distvectors = []

            for count2, neigh in enumerate(self.atoms["neighbors"][count]):
                pos2 = self.atoms["positions"][neigh]
                dist = self.atoms["neighbordist"][count][count2]
                vectors = self.atoms["diff"][count][count2]
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

        self.atoms["angular"] = angulars
        
        mapdict = {}
        mapdict["angular_parameters"] = {}
        mapdict["angular_parameters"]["diamond_angle"] = "angular"
        self.atoms._add_attribute(mapdict)

    def calculate_chiparams(self, angles=False):
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

        self._check_neighbors()

        bins = [-1.0, -0.945, -0.915, -0.755, -0.705, -0.195, 0.195, 0.245, 0.795, 1.0]
        chiparams = []
        cosines = []

        for count, pos in enumerate(self.atoms["positions"]):

            dists = self.atoms["neighbordist"][count]
            neighs = self.atoms["neighbors"][count]

            args = range(len(dists))
            combos = list(itertools.combinations(args, 2))
            costhetas = []
            
            for combo in combos:
                vec1 = self.atoms["diff"][count][combo[0]]
                vec2 = self.atoms["diff"][count][combo[1]]
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
        
        self.atoms["chiparams"] = chiparams
        
        mapdict = {}
        mapdict["angular_parameters"] = {}
        mapdict["angular_parameters"]["chi_params"] = "chiparams"
        
        if angles:
            self.atoms["cosines"] = cosines
            mapdict["angular_parameters"]["cosines"] = "cosines"

        self.atoms._add_attribute(mapdict)


    def calculate_cna(self, lattice_constant=None):
        """
        Calculate the Common Neighbor Analysis indices

        Parameters
        ----------

        lattice_constant : float, optional
            lattice constant to calculate CNA. If not specified,
            adaptive CNA will be used

        Returns
        -------
        resdict: dict
            dictionary of calculated structure

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
        resdict = cna.calculate_cna(self)
        return resdict

    def identify_diamond(self):
        """
        Identify diamond structure

        Parameters
        ----------
        None

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
        resdict = cna.identify_diamond(self)
        return resdict 

    def calculate_centrosymmetry(self, nmax=12):
        """
        Calculate the centrosymmetry parameter

        Parameters
        ----------        
        nmax : int, optional
            number of neighbors to be considered for centrosymmetry 
            parameters. Has to be a positive, even integer. Default 12

        Returns
        -------
        None
        
        Notes
        -----
        Calculate the centrosymmetry parameter for each atom which can be accessed by
        :attr:`~pyscal.catom.centrosymmetry` attribute. It calculates the degree of inversion
        symmetry of an atomic environment. Centrosymmetry recalculates the neighbor using
        the number method as specified in :func:`¬pyscal.core.System.find_neighbors` method. This
        is the ensure that the required number of neighbors are found for calculation of the parameter.

        The Greedy Edge Selection (GES) [1] as specified in [2] in used in this method. 
        GES algorithm is implemented in LAMMPS and Ovito. Please see [2] for
        a detailed description of the algorithms. 
        References
        ----------
        .. [1] Stukowski, A, Model Simul Mater SC 20, 2012
        .. [2] Larsen, arXiv:2003.08879v1, 2020

        """
        return pyscal3.operations.centrosymmetry.calculate_centrosymmetry(self, nmax)


          
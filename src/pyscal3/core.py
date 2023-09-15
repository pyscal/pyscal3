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
        self.find._add_attribute(mapdict)

        self.calculate = AttrSetter()
        mapdict = {}
        mapdict['distance'] = update_wrapper(partial(neighbor.get_distance, self), neighbor.get_distance)
        mapdict['average_over_neighbors'] = update_wrapper(partial(calculations.average_over_neighbors, self), calculations.average_over_neighbors)
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


    def calculate_q(self, q, averaged=False, continuous_algorithm=False):
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

        self._check_neighbors()

        if averaged:
            self._calculate_aq(qq)
            qvals = [self.atoms["avg_q%d"%x] for x in qq]
        else:    
            if continuous_algorithm:
                lm = max(qq)
                pc.calculate_q(self.atoms, lm)
            else:
                self._calculate_q(qq)
            qvals = [self.atoms["q%d"%x] for x in qq]
        return qvals

    def _calculate_q(self, qq):
        """
        Private method for calculation of qvals
        """
        for val in qq:
            pc.calculate_q_single(self.atoms, val)
  
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
        self.atoms._add_attribute(mapdict)


    def _calculate_aq(self, qq):
        """
        Private method for calculation of avged qvals
        """

        todo_q = []
        for q in qq:
            keys = ["q%d"%q, "q%d_real"%q, "q%d_imag"%q]
            prod = []
            for key in keys:
                if key in self.atoms.keys():
                    prod.append(True)
                else:
                    prod.append(False)
            prod = np.prod(prod)
            if not prod:
                todo_q.append(q)

        _ = self._calculate_q(todo_q)

        #loop over atoms
        for val in qq:
            pc.calculate_aq_single(self.atoms, val)

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
        self.atoms._add_attribute(mapdict)

    def calculate_disorder(self, averaged=False, q=6):
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
            if key in self.atoms.keys():
                prod.append(True)
            else:
                prod.append(False)
        prod = np.prod(prod)
        if not prod:
            self.calculate_q(q)

        pc.calculate_disorder(self.atoms, q)

        mapdict = {}
        mapdict["steinhardt"] = {}
        mapdict["steinhardt"]["disorder"] = {}
        mapdict["steinhardt"]["disorder"]["norm"] = "disorder"

        if averaged:
            #average the disorder
            avg_arr = self.calculate.average_over_neighbors("disorder")
            self.atoms["avg_disorder"] = avg_arr
            mapdict["steinhardt"]["disorder"]["average"] = "avg_disorder"
        self.atoms._add_attribute(mapdict)


    def find_solids(self, bonds=0.5, threshold=0.5, avgthreshold=0.6, 
                          cluster=True, q=6, cutoff=0, right=True):
        """
        Distinguish solid and liquid atoms in the system.
        Parameters
        ----------
        bonds : int or float, optional
            Minimum number of solid bonds for an atom to be identified as
            a solid if the value is an integer. Minimum fraction of neighbors
            of an atom that should be solid for an atom to be solid if the
            value is float between 0-1. Default 0.5.
        
        threshold : double, optional
            Solid bond cutoff value. Default 0.5.
        
        avgthreshold : double, optional
            Value required for Averaged solid bond cutoff for an atom to be identified
            as solid. Default 0.6.
        
        cluster : bool, optional
            If True, cluster the solid atoms and return the number of atoms in the largest
            cluster.
        
        q : int, optional
            The Steinhardt parameter value over which the bonds have to be calculated.
            Default 6.
        
        cutoff : double, optional
            Separate value used for cluster classification. If not specified, cutoff used
            for finding neighbors is used.
        
        right: bool, optional
            If true, greater than comparison is to be used for finding solid particles. 
            default True.
        
        Returns
        -------
        solid : int
            Size of the largest solid cluster. Returned only if `cluster=True`.
        
        Notes
        -----
        The neighbors should be calculated before running this function.
        Check :func:`~pyscal.core.System.find_neighbors` method.
        
        `bonds` define the number of solid bonds of an atom to be identified as solid.
        Two particles are said to be 'bonded' if [1],
        .. math:: s_{ij} = \sum_{m=-6}^6 q_{6m}(i) q_{6m}^*(i) \geq threshold
        where `threshold` values is also an optional parameter.
        
        If the value of `bonds` is a fraction between 0 and 1, at least that much of an atom's neighbors
        should be solid for the atom to be solid.
        
        An additional parameter `avgthreshold` is an additional parameter to improve solid-liquid distinction.
        
        In addition to having a the specified number of `bonds`,
        
        .. math::  \langle s_{ij} \\rangle > avgthreshold
        
        also needs to be satisfied. In case another q value has to be used for calculation of S_ij, it can be
        set used the `q` attribute. In the above formulations, `>` comparison for `threshold` and `avgthreshold`
        can be changed to `<` by setting the keyword `right` to False.
        
        If `cluster` is True, a clustering is done for all solid particles. See :func:`~pyscal.csystem.find_clusters`
        for more details. 
        
        References
        ----------
        .. [1] Auer, S, Frenkel, D. Adv Polym Sci 173, 2005
        """
        #check if neighbors are found
        self._check_neighbors()

        if not isinstance(q, int):
            raise TypeError("q should be interger value")

        if not isinstance(threshold, (int, float)):
            raise TypeError("threshold should be a float value")
        else:
            if not ((threshold >= 0 ) and (threshold <= 1 )):
                raise ValueError("Value of threshold should be between 0 and 1")

        if not isinstance(avgthreshold, (int, float)):
            raise TypeError("avgthreshold should be a float value")
        else:
            if not ((avgthreshold >= 0 ) and (avgthreshold <= 1 )):
                raise ValueError("Value of avgthreshold should be between 0 and 1")

        #start identification routine
        #check the value of bonds and set criteria depending on that
        if isinstance(bonds, int):
            criteria = 0
        elif isinstance(bonds, float):
            if ((bonds>=0) and (bonds<=1.0)):
                criteria = 1
            else:
                raise TypeError("bonds if float should have value between 0-1")
        else:
             raise TypeError("bonds should be interger/float value")

        if right:
            compare_criteria = 0
        else:
            compare_criteria = 1

        self.calculate_q(q)

        #calculate solid neighs
        pc.calculate_bonds(self.atoms, q, 
            threshold, avgthreshold, bonds, 
            compare_criteria, criteria)

        mapdict = {}
        mapdict["steinhardt"] = {}
        mapdict["steinhardt"]["order"] = {}
        mapdict["steinhardt"]["order"]["bonds"] = "bonds"
        mapdict["steinhardt"]["order"]["sij"] = {}
        mapdict["steinhardt"]["order"]["sij"]["norm"] = "sij"
        mapdict["steinhardt"]["order"]["sij"]["average"] = "avg_sij"
        mapdict["steinhardt"]["order"]["sij"]["solid"] = "solid"
        self.atoms._add_attribute(mapdict)
        
        if cluster:
            lc = self.cluster_atoms(self.atoms.steinhardt.order.sij.solid, largest=True)
            return lc

    def find_largest_cluster(self):
        """
        Find largest cluster among given clusters
        
        Parameters
        ----------
        None

        Returns
        -------
        lc : int
            Size of the largest cluster.
        """
        if not "cluster" in self.atoms.keys():
            raise RuntimeError("cluster_atoms needs to be called first")

        clusterlist = [x for x in self.atoms["cluster"] if x != -1]
        xx, xxcounts = np.unique(clusterlist, return_counts=True)
        arg = np.argsort(xxcounts)[-1]
        largest_cluster_size = xxcounts[arg]
        largest_cluster_id = xx[arg]


        self.atoms["largest_cluster"] = [True if self.atoms["cluster"][x]==largest_cluster_id else False for x in range(len(self.atoms["cluster"]))]
        
        mapdict = {}
        mapdict["cluster"] = {}
        mapdict["cluster"]["largest"] = "largest_cluster"
        self.atoms._add_attribute(mapdict)

        return largest_cluster_size


    def cluster_atoms(self, condition, largest = True, cutoff=0):
        """
        Cluster atoms based on a property
        
        Parameters
        ----------
        condition : callable or atom property
            Either function which should take an :class:`~Atom` object, and give a True/False output
            or an attribute of atom class which has value or 1 or 0.
        
        largest : bool, optional
            If True returns the size of the largest cluster. Default False.
        
        cutoff : float, optional
            If specified, use this cutoff for calculation of clusters. By default uses the cutoff
            used for neighbor calculation.
        
        Returns
        -------
        lc : int
            Size of the largest cluster. Returned only if `largest` is True.
        
        Notes
        -----
        This function helps to cluster atoms based on a defined property. This property
        is defined by the user through the argument `condition` which is passed as a parameter.
        `condition` should be a boolean array the same length as number of atoms in the system.
        """
        
        self.apply_selection(condition=condition)
        pc.find_clusters(self.atoms, cutoff)

        mapdict = {}
        mapdict["cluster"] = {}
        mapdict["cluster"]["id"] = "cluster"
        self.atoms._add_attribute(mapdict)

        #done!
        lc = self.find_largest_cluster()
        #pcs.System.get_largest_cluster_atoms(self)
        self.remove_selection()
        if largest:
            return lc


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


          
import pyscal3.traj_process as ptp

def read_inputfile(system, filename, format="lammps-dump", 
                                        compressed = False, customkeys=None):
    """

    Read input file that contains the information of system configuration.

    Parameters
    ----------
    filename : string
        name of the input file.

    format : {'lammps-dump', 'poscar', 'ase', 'mdtraj'}
        format of the input file, in case of `ase` the ASE Atoms object

    compressed : bool, optional
        If True, force to read a `gz` compressed format, default False.

    customkeys : list
        A list containing names of headers of extra data that needs to be read in from the
        input file.

    Returns
    -------
    None

    Notes
    -----
    `format` keyword specifies the format of the input file. Currently only
    a `lammps-dump` and `poscar` files are supported.  Additionaly, the widely
    use Atomic Simulation environment (https://wiki.fysik.dtu.dk/ase/ase/ase.html).
    mdtraj objects (http://mdtraj.org/1.9.3/) are also supported by using the keyword
    `'mdtraj'` for format. Please note that triclinic boxes are not yet supported for
    mdtraj format.
    Atoms object can also be used directly. This function uses the
    :func:`~pyscal.traj_process` module to process a file which is then assigned to system.

    `compressed` keyword is not required if a file ends with `.gz` extension, it is
    automatically treated as a compressed file.

    Triclinic simulation boxes can also be read in.

    If `custom_keys` are provided, this extra information is read in from input files if
    available. This information is can be accessed directly as `self.atoms['customkey']`


    """
    system.neighbors_found = False
    system.neighbor_method = None
    system.ghosts_created = False
    system.actual_box = None

    #box parameters
    system.rot = [[0,0,0], [0,0,0], [0,0,0]]
    system.rotinv = [[0,0,0], [0,0,0], [0,0,0]]
    system.boxdims = [0,0,0]
    system.triclinic = 0

    atoms, box = ptp.read_file(filename, format=format, 
                                compressed=compressed, customkeys=customkeys,)
    system.box = box
    system.atoms = atoms
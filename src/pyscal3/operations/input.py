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

def to_file(system, outfile, format='lammps-dump', customkeys=None, customvals=None,
            compressed=False, timestep=0, species=None):
    """
    Save the system instance to a trajectory file.

    Parameters
    ----------
    outfile : string
        name of the output file

    format : string, {'lammps-dump', 'lammps-data', 'poscar'}
        format of the output file, default `lammps-dump`
        Currently only `lammps-dump` format is supported.

    customkeys : list of strings, optional
        a list of extra atom wise values to be written in the output file.

    customvals : list or list of lists, optional
        If `customkey` is specified, `customvals` take an array of the same length
        as number of atoms, which contains the values to be written out.

    compressed : bool, optional
        If true, the output is written as a compressed file.

    timestep : int, optional
        timestep to be written to file. default 0

    species : None, optional
        species of the atoms. Required if any format other than 'lammps-dump' is used. Required
        for convertion to ase object.

    Returns
    -------
    None

    Notes
    -----
    `to_file` method can handle a number of file formats. The most customizable format is the
    `lammps-dump` which can take a custom options using customkeys and customvals. customkeys
    will be the header written to the dump file. It can be any Atom attribute, any property
    stored in custom variable of the Atom, or calculated q values which can be given by `q4`,
    `aq4` etc. External values can also be provided using `customvals` option. `customvals` array
    should be of the same length as the number of atoms in the system.

    For all other formats, ASE is used to write out the file, and hence the `species` keyword
    needs to be specified. If initially, an ASE object was used to create the System, `species`
    keyword will already be saved, and need not be specified. In other cases, `species` should
    be a list of atomic species in the System. For example `["Cu"]` or `["Cu", "Al"]`, depending
    on the number of species in the System. In the above case, atoms of type 1 will be mapped to
    Cu and of type 2 will be mapped to Al. For a complete list of formats that ASE can handle,
    see `here <https://wiki.fysik.dtu.dk/ase/ase/io/io.html>`_ . 
    """

    ptp.write_file(system, outfile, format = format,
        compressed = compressed, customkeys = customkeys, customvals = customvals,
        timestep = timestep, species = species)
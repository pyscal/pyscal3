"""
Trajectory module for lazy reading of LAMMPS dump files.

Provides efficient block-level access to large trajectory files
with support for slicing, lazy loading, and conversion to ASE Atoms.
"""

import os
import numpy as np
from ase import Atoms as ASEAtoms
import warnings


def _parse_lammps_lines_to_atoms(lines, species=None, customkeys=None):
    """Parse LAMMPS dump lines (list of str) into an ASE Atoms object.

    Parameters
    ----------
    lines : list of str
        Raw lines from one dump block (header + atom data).
    species : list of str, optional
        Element symbols mapped to LAMMPS type indices, e.g. ``['Cu', 'Zr']``.
    customkeys : list of str, optional
        Extra per-atom columns to store in ``atoms.arrays``.

    Returns
    -------
    ase.Atoms
    """
    if customkeys is None:
        customkeys = []

    # --- parse header ---
    natoms = int(lines[3].strip())

    # box bounds (lines 5, 6, 7)
    box_raw = [line.strip().split() for line in lines[5:8]]
    triclinic = len(box_raw[0]) == 3

    boxx = [float(box_raw[0][0]), float(box_raw[0][1])]
    boxy = [float(box_raw[1][0]), float(box_raw[1][1])]
    boxz = [float(box_raw[2][0]), float(box_raw[2][1])]

    if triclinic:
        xy = float(box_raw[0][2])
        xz = float(box_raw[1][2])
        yz = float(box_raw[2][2])

    # column header (line 8)
    header = lines[8].strip().split()
    headerdict = {header[i]: i - 2 for i in range(len(header))}

    # detect scaled vs unscaled coordinates
    if "x" in headerdict:
        scaled = False
    elif "xs" in headerdict:
        scaled = True
        headerdict["x"] = headerdict.pop("xs")
        headerdict["y"] = headerdict.pop("ys")
        headerdict["z"] = headerdict.pop("zs")
    else:
        raise ValueError("Only x/xs, y/ys, z/zs keys are supported in LAMMPS dump")

    # --- parse atom data ---
    ids = np.empty(natoms, dtype=int)
    types = np.empty(natoms, dtype=int)
    positions = np.empty((natoms, 3), dtype=float)
    custom_data = {k: [] for k in customkeys}

    for i, line in enumerate(lines[9 : 9 + natoms]):
        raw = line.strip().split()
        ids[i] = int(raw[headerdict["id"]])
        types[i] = int(raw[headerdict["type"]])
        positions[i] = [
            float(raw[headerdict["x"]]),
            float(raw[headerdict["y"]]),
            float(raw[headerdict["z"]]),
        ]
        for k in customkeys:
            if k in headerdict:
                custom_data[k].append(raw[headerdict[k]])

    # --- build cell ---
    # LAMMPS dumps the *bounding box* of a triclinic cell; recover the
    # actual box limits (xlo, xhi, ...) from the bounds and the tilts.
    if triclinic:
        amin = min(0.0, xy, xz, xy + xz)
        amax = max(0.0, xy, xz, xy + xz)
        bmin = min(0.0, yz)
        bmax = max(0.0, yz)
        xlo = boxx[0] - amin
        xhi = boxx[1] - amax
        ylo = boxy[0] - bmin
        yhi = boxy[1] - bmax
        zlo = boxz[0]
        zhi = boxz[1]

        a = np.array([xhi - xlo, 0.0, 0.0])
        b = np.array([xy, yhi - ylo, 0.0])
        c = np.array([xz, yz, zhi - zlo])
        cell = np.array([a, b, c])
        origin = np.array([xlo, ylo, zlo])
    else:
        cell = np.diag([boxx[1] - boxx[0], boxy[1] - boxy[0], boxz[1] - boxz[0]])
        origin = np.array([boxx[0], boxy[0], boxz[0]])

    # handle scaled coordinates: fractional w.r.t. the cell, measured from
    # the box origin (as in LAMMPS), so convert and add the origin back to
    # obtain absolute coordinates consistent with the unscaled case
    if scaled:
        positions = positions @ cell + origin

    # --- determine species ---
    if species is not None:
        species = np.atleast_1d(species)
        symbols = [species[t - 1] for t in types]
    else:
        # default to element symbol from type index
        symbols = ["X"] * natoms

    # --- create ASE Atoms ---
    # positions are absolute (as written by LAMMPS); the box origin is kept
    # in celldisp like ASE's own LAMMPS dump reader does
    atoms = ASEAtoms(symbols=symbols, positions=positions, cell=cell, pbc=True)
    atoms.set_celldisp(origin)

    # store LAMMPS metadata
    atoms.arrays["lammps_ids"] = ids
    atoms.arrays["lammps_types"] = types

    for k in customkeys:
        if k in custom_data and len(custom_data[k]) == natoms:
            try:
                atoms.arrays[k] = np.array(custom_data[k], dtype=float)
            except ValueError:
                atoms.info[k] = custom_data[k]

    return atoms


class Timeslice:
    """
    A slice of a trajectory containing one or more consecutive timesteps.

    Timeslices are created by indexing a :class:`Trajectory` object. Multiple
    slices can be concatenated with ``+``.

    Attributes
    ----------
    trajectory : Trajectory
        Parent trajectory.
    blocklist : range or list of int
        Indices of the blocks in this slice.
    """

    def __init__(self, trajectory, blocklist):
        """
        Parameters
        ----------
        trajectory : Trajectory
            The source trajectory.
        blocklist : range or list of int
            Block indices for this slice.
        """
        self.trajectory = trajectory
        self.blocklist = blocklist
        self.trajectories = [trajectory]
        self.blocklists = [blocklist]

    def __repr__(self):
        """
        String of the class
        """
        blockstring = ["%d-%d" % (x[0], x[-1]) for x in self.blocklists]
        blockstring = "/".join(blockstring)

        data = "Trajectory slice\n %s\n natoms=%d\n" % (
            blockstring,
            self.trajectory.natoms,
        )
        return data

    def __add__(self, ntraj):
        """
        Add a method for summing
        """
        for traj in ntraj.trajectories:
            self.trajectories.append(traj)

        for block in ntraj.blocklists:
            self.blocklists.append(block)

        return self

    def __radd__(self, ntraj):
        """
        Reverse add method
        """
        if ntraj == 0:
            return self
        else:
            return self.__add__(ntraj)

    def to_atoms(self, species=None, customkeys=None):
        """
        Get blocks as ASE Atoms objects.

        Parameters
        ----------
        species : list of str, optional
            Species names mapped to LAMMPS types, e.g. ['Cu', 'Zr'].
            Required if the dump file does not contain element information.
        customkeys : list of str, optional
            Extra per-atom columns to read from the dump file.

        Returns
        -------
        list of ase.Atoms
        """
        atoms_list = []
        for count, traj in enumerate(self.trajectories):
            for x in self.blocklists[count]:
                a = self.trajectories[count]._get_block_as_atoms(
                    x, species=species, customkeys=customkeys
                )
                atoms_list.append(a)
        return atoms_list

    # Backward-compatible aliases
    def to_ase(self, species=None, customkeys=None):
        """Alias for :meth:`to_atoms`."""
        return self.to_atoms(species=species, customkeys=customkeys)

    def to_system(self, customkeys=None, species=None):
        """Alias for :meth:`to_atoms` (legacy name)."""
        return self.to_atoms(species=species, customkeys=customkeys)

    def to_dict(self):
        """
        Get raw lines for each block in this slice.

        Returns
        -------
        list of list of str
            Raw line data for each block.
        """
        data = []
        for count, traj in enumerate(self.trajectories):
            for x in self.blocklists[count]:
                self.trajectories[count].load(x)
                data.append(self.trajectories[count].data[x])
                self.trajectories[count].unload(x)
        return data

    def to_file(self, outfile, mode="w"):
        """
        Get block as outputfile

        Parameters
        ----------
        outfile : string
            name of output file

        mode : string
            write mode to be used, optional
            default "w" write
            also can be "a" to append.

        Returns
        -------
        None

        """
        # the lines carry their own newlines (copied from the source file or
        # generated with "\n"); newline="" stops the text layer from turning
        # them into "\r\n" on Windows
        with open(outfile, mode, encoding="utf-8", newline="") as fout:
            for count, traj in enumerate(self.trajectories):
                self.trajectories[count]._get_blocks_to_file(fout, self.blocklists[count])


class Trajectory:
    """
    Lazy reader for LAMMPS dump trajectory files.

    Supports block-level random access via indexing and slicing.
    Blocks are loaded on demand and can be converted to ASE Atoms.

    Parameters
    ----------
    filename : str
        Path to a LAMMPS dump file.

    Attributes
    ----------
    nblocks : int
        Number of timestep blocks in the file.
    natoms : int
        Number of atoms per block.

    Examples
    --------
    >>> traj = pyscal3.Trajectory("dump.lammpstrj")
    >>> ts = traj[0]                    # single block → Timeslice
    >>> atoms_list = ts.to_atoms(species=["Cu"])
    >>> atoms_list = traj[0:10].to_atoms(species=["Cu"])  # slice
    """

    def __init__(self, filename):
        """
        Parameters
        ----------
        filename : str
            Path to a LAMMPS dump file.
        """
        if os.path.exists(filename):
            self.filename = filename
        else:
            raise FileNotFoundError("%s file not found" % filename)
        self.natoms = 0
        self.blocksize = 0
        self.nblocks = 0
        self.loadlist = None
        self.data = None

        self._get_natoms()
        self._get_nblocks()

    def __repr__(self):
        """
        String of the class
        """
        return "Trajectory of %d slices with %d atoms" % (self.nblocks, self.natoms)

    def __getitem__(self, blockno):
        """
        Allow for slice indexing
        """
        if isinstance(blockno, slice):
            blocklist = range(*blockno.indices(self.nblocks))
            timeslice = Timeslice(self, blocklist)
            return timeslice
        else:
            blocklist = [blockno]
            timeslice = Timeslice(self, blocklist)
            return timeslice

    def _get_natoms(self):
        """
        Get number of atoms in the system

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        with open(self.filename, "rb") as fout:
            data = [next(fout) for x in range(0, 4)]
        self.natoms = int(data[-1])

    def _get_nlines(self):
        """
        Get total number of lines in the file

        Parameters
        ----------
        None

        Returns
        -------
        nlines : int
            number of lines
        """
        line_offset = []
        offset = 0
        nlines = 0
        with open(self.filename, "rb") as fin:
            for line in fin:
                line_offset.append(offset)
                offset += len(line)
                nlines += 1

        self.nlines = nlines
        self.line_offset = line_offset
        return nlines

    def _get_nblocks(self):
        """
        Get number of blocks in the trajectory file

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self._get_natoms()
        nlines = self._get_nlines()
        self.blocksize = self.natoms + 9
        self.nblocks = nlines // self.blocksize
        self.straylines = nlines - self.nblocks * self.blocksize
        # set load list to False
        self.loadlist = [False for x in range(self.nblocks)]
        self.data = [None for x in range(self.nblocks)]

    def get_block(self, blockno):
        """
        Get a block from the file as raw data

        Parameters
        ----------
        blockno : int
            number of the block to be read, starts from 0

        Returns
        -------
        data : list
            list of strings containing data
        """
        start = blockno * self.blocksize

        data = []
        with open(self.filename, "rb") as fin:
            fin.seek(self.line_offset[start])
            for i in range(self.blocksize):
                data.append(fin.readline().decode("utf-8"))
        return data

    def load(self, blockno):
        """
        Load the data of a block into memory as a dictionary
        of numpy arrays

        Parameters
        ----------
        blockno : int
            number of the block to be read, starts from 0

        Returns
        -------
        None

        Notes
        -----
        When the data of a block is loaded, it is accessible
        through `Trajectory.data[x]`. This data can then be
        modified. When the block is written out, the modified
        data is written instead of existing one. But, loaded
        data is kept in memory until unloaded using `unload`
        method.
        """
        data = self.get_block(blockno)
        box = np.loadtxt(data[5:8])
        columns = np.loadtxt(data[9:])
        header = np.loadtxt(data[8:9], dtype=str)[2:]
        outdict = {}
        outdict["box"] = box
        outdict["atoms"] = {}
        for count, h in enumerate(header):
            outdict["atoms"][h] = columns[:, count]

        self.data[blockno] = outdict
        self.loadlist[blockno] = True

    def unload(self, blockno):
        """
        Unload the data that is loaded to memory using
        `load` method

        Parameters
        ----------
        blockno : int
            number of the block to be read, starts from 0

        Returns
        -------
        None
        """
        self.data[blockno] = None
        self.loadlist[blockno] = False

    def _convert_data_to_lines(self, blockno):
        """
        Create lines from loaded data

        Parameters
        ----------
        blockno : int
            number of the block to be read, starts from 0

        Returns
        -------
        data : list of strs
            list of lines
        """
        dd = self.data[blockno]

        data = []
        data.append("ITEM: TIMESTEP\n")
        data.append("".join([str(0), "\n"]))
        data.append("ITEM: NUMBER OF ATOMS\n")
        data.append("".join([str(self.natoms), "\n"]))
        data.append("ITEM: BOX BOUNDS pp pp pp\n")
        for b in dd["box"]:
            dstr = " ".join(b.astype(str))
            data.append("".join([dstr, "\n"]))

        xf = []
        xd = []
        xfkeys = []
        xdkeys = []
        for key, val in dd["atoms"].items():
            if key in ["id", "type"]:
                val.astype(int)
                xd.append(val)
                xdkeys.append(key)
            else:
                xf.append(val)
                xfkeys.append(key)

        xdstrs = []
        if len(xd) > 0:
            for i in range(len(xd[0])):
                substr = []
                for j in range(len(xdkeys)):
                    substr.append("%d" % xd[j][i])
                xdstrs.append(" ".join(substr))

        xdheader = " ".join(xdkeys)
        xdheader = " ".join(["ITEM: ATOMS", xdheader])

        xfstrs = []
        xf = np.array(xf)
        if len(xf) > 0:
            for i in range(len(xf[0])):
                dstr = " ".join((xf[:, i]).astype(str))
                xfstrs.append("".join([dstr, "\n"]))

        xfheader = " ".join(xfkeys)
        mainheader = " ".join([xdheader, xfheader])
        mainheader = "".join([mainheader, "\n"])

        data.append(mainheader)

        for i in range(len(xfstrs)):
            valstr = " ".join([xdstrs[i], xfstrs[i]])
            data.append(valstr)

        return data

    def _get_block_as_atoms(self, blockno, species=None, customkeys=None):
        """
        Parse a block from the dump file directly into an ASE Atoms object.

        Parameters
        ----------
        blockno : int
            Block index (0-based).
        species : list of str, optional
            Element names mapped to LAMMPS types.
        customkeys : list of str, optional
            Extra per-atom columns to read.

        Returns
        -------
        ase.Atoms
        """
        data = self.get_block(blockno)
        return _parse_lammps_lines_to_atoms(
            data, species=species, customkeys=customkeys
        )

    def _get_blocks_to_file(self, fout, blocklist):
        """
        Write a series of blocks to an open file handle.

        Parameters
        ----------
        fout : file object
            open, writable file handle
        blocklist : iterable of int
            block numbers to write, starting from 0

        Returns
        -------
        None
        """
        xl = [x for x in blocklist]
        xl = np.array(xl)

        # open file
        # convert lines to start from end
        for x in xl:
            if self.loadlist[x]:
                data = self._convert_data_to_lines(x)
            else:
                data = self.get_block(x)
            for line in data:
                fout.write(line)

"""module for handling plain text 1D profiles sorted in
columns and with the following structure:
# col names
length of data(rows int)
<data block>
# comments or anything
"""

import linecache
import os
import numpy as np
import operator
from flashy.utils import by_mass


class DataMatrix(object):
    """Object for interacting with columned data in a plain text file.

    Format for file is strict at the start since it uses np.genfromtxt
    to handle the data block. After row 3 any comment is allowed (# marked).
    Essential columns: radius and density.
    Any tag with a number becomes a species.

    File format:
        # col names
        length of data(rows int)
        <data block>
        # comments or anything

    Attributes:
        radius (np.ndarray): Radial coordinates.
        density (np.ndarray): Density values.
        masses (np.ndarray): Cumulative mass enclosed.
        species (list): Species names found in the data.
        bulkprops (list): Bulk property names found in the data.
        meta (dict): Metadata including mass, central density, radius, resolution.
    """

    _rnames = ["radius", "r", "rad"]
    _dnames = ["density", "rho", "dens"]
    _mnames = ["masses", "m", "mass"]
    _tnames = ["temperature", "temp"]
    _pnames = ["pressure", "pres"]

    _classn = ["species", "data", "filekeys", "bulkprops", "meta"]

    def __init__(self, filename, comment=""):
        if isinstance(filename, list):
            self.filekeys, self.data = filename[0], filename[1]
        else:
            self.filekeys, self.data = chopFile(filename)
        self.species = []
        self.bulkprops = []
        self.meta = {}

        # sift data
        skip = []
        _, num = multiIndex(["radius", "r", "rad"], self.filekeys)
        self.radius = self.data[:, num]
        skip.append(num)
        _, num = multiIndex(["rho", "dens", "density"], self.filekeys)
        self.density = self.data[:, num]
        skip.append(num)
        for i, k in enumerate(self.filekeys):
            if i not in skip:
                setattr(self, k, self.data[:, i])
            else:
                continue
        self.masses = by_mass(self.radius, self.density)
        self.setMeta(comment=comment)

    def __setattr__(self, name, value):
        # run name through the 'dictionary' of names
        # bulk = True
        if name in self._dnames:
            rname = "density"
        elif name in self._rnames:
            rname = "radius"
        elif name in self._mnames:
            rname = "masses"
        elif name in self._tnames:
            rname = "temperature"
        elif name in self._pnames:
            rname = "pressure"
        else:
            rname = name
        if rname in self._classn:  # skip init settings
            return super().__setattr__(rname, value)
        if rname not in self.__dict__:
            if any(char.isdigit() for char in rname) or len(rname) == 1:
                self.__dict__["species"].append(rname)
            else:
                self.__dict__["bulkprops"].append(rname)
            return super().__setattr__(rname, value)
        else:
            return super().__setattr__(rname, value)

    def __getattr__(self, name):
        # run name through the 'dictionary' of names
        if name in self._dnames:
            return self.__getattribute__("density")
        elif name in self._rnames:
            return self.__getattribute__("radius")
        elif name in self._mnames:
            return self.__getattribute__("masses")
        elif name in self._tnames:
            return self.__getattribute__("temperature")
        elif name in self._pnames:
            return self.__getattribute__("pressure")
        else:
            return self.__getattribute__(name)

    def setMeta(self, comment=""):
        self.meta["mass"] = self.masses[-1]
        self.meta["central dens"] = self.density[0]
        self.meta["radius"] = self.radius[-1]
        distance = self.radius[-1] - self.radius[0]
        self.meta["resolution"] = distance / len(self.radius)
        self.meta["points"] = len(self.radius)
        if comment:
            self.meta["comment"] = comment
        elif "comment" not in self.meta:
            self.meta["comment"] = ""

    def printMeta(self):
        for k, v in self.meta.items():
            pformat = "{}: {:}" if isinstance(v, str) else "{}: {:E}"
            print(pformat.format(k.capitalize(), v))

    def writeProf(self, output, subset=[], autotag=False, Xthresh=1e-15):
        """Write profile to file
        bug: ndarray ninja breaks formatting (np.array([float]))

        Args:
            output(str): otp filename.
            subset(list): write a subset of keys to file.
            autotag(bool): setup name based on the profile (WDs)

        """
        basedir = os.path.dirname(output)
        if autotag:
            rhoc = "{:2.4e}".format(self.density[0])
            rhoc = rhoc.replace("+", "").replace(".", "p")
            c12 = "{:2.2f}".format(self.c12[0])
            c12 = c12.replace("+", "").replace(".", "p")
            o16 = "{:2.2f}".format(self.o16[0])
            o16 = o16.replace("+", "").replace(".", "p")
            filename = "wd_fermi_helm_{}_{}_{}.dat".format(rhoc, c12, o16)
        else:
            filename = os.path.basename(output)
        otp = os.path.join(basedir, filename)
        if subset:
            keys = DataMatrix.checkbulkprops(subset)
        else:
            keys = self.bulkprops + self.species
        print("Writing: {}".format(" ".join(keys)))
        missing = set()
        with open(otp, "w") as f:
            header = " ".join(keys)
            f.write("# {}\n".format(header))
            f.write("{}\n".format(self.meta["points"]))
            for i in range(self.meta["points"]):
                line = []
                for k in keys:
                    try:
                        # clip stupid numbers such as 1e-99
                        value = getattr(self, k)[i]
                        if value < Xthresh or np.isnan(value):
                            value = 0.0
                        line.append("{:15.8e}".format(value))
                    except AttributeError:
                        missing.add(k)
                        line.append("{:15.8e}".format(0.0))
                f.write("{}\n".format(" ".join(line)))
            if self.meta["comment"]:
                f.write("# {}\n".format(self.meta["comment"]))
            f.write("# Total Mass {} Msun".format(self.meta["mass"]))
        if missing:
            print("Missing keys: " "{} were set to zero.".format(" ".join(missing)))
        print("Wrote: {}".format(otp))

    @staticmethod
    def checkbulkprops(klist):
        nklist = []
        for k in klist:
            # hard bulk props are always named density and radius
            if k in DataMatrix._rnames:
                # key requested is a valid radius name
                nklist.append("radius")
            elif k in DataMatrix._dnames:
                # key requested is a valid density name
                nklist.append("density")
            else:
                # key is something else.
                nklist.append(k)
        return nklist


def spliceProfs(left, right):
    """Joins profiles at ends, 'right' takes precedence on overlap.

    Note: bulkprops should stay ordered unless a new one is added to
    DataMatrix. Now species mix so one might need to do a sortnuclides
    down the road.

    Args:
        left(DataMatrix): innermost profile.
        right(DataMatrix): outermost profile.

    Returns:
        (DataMatrix): spliced profile.

    """
    offs = right.radius[0]
    left = snipProf(left, offs)
    hsr = np.hstack((left.radius, right.radius))
    hsd = np.hstack((left.density, right.density))
    dblock = np.column_stack((hsr, hsd))
    lprops = left.bulkprops + left.species
    rprops = right.bulkprops + right.species
    if lprops == rprops:
        keys = left.filekeys
    else:
        # set() yields unordered members
        # simply check repeated since they are already ordered from the object
        keys = []
        for k in lprops + rprops:
            if k not in keys:
                keys.append(k)
    for k in keys[2:]:
        try:
            latt = getattr(left, k)
        except AttributeError:
            latt = np.zeros(len(left.radius))
        try:
            ratt = getattr(right, k)
        except AttributeError:
            ratt = np.zeros(len(right.radius))
        hst = np.hstack((latt, ratt))
        dblock = np.column_stack((dblock, hst))
    return DataMatrix([keys, dblock])


def snipProf(orig, cut, byM=False, left=True, edgecell=False):
    """Cuts a profile, returning a new profile object.

    Conv: center is 0, edge is -1.

    Args:
        orig(DataMatrix): DataMatrix object to cut.
        cut(float): cut coordinate.
        byM(bool): specify cut is by mass coordinate.
        left(bool): return data at the left/right of the cut.

    Returns:
        (DataMatrix): new DataMatrix object.

    """
    abscissa = orig.masses if byM else orig.radius
    npabs = np.array(abscissa)
    flow = operator.le(npabs, cut) if left else operator.ge(npabs, cut)
    if np.any(flow) is False:
        print("Cut outside the profile, returning whole profile")
        return orig
    allcells = np.where(flow)
    # remove edge cell
    # toggle was added late, it broke some plots so it was removed by doing the same for both cases
    # if edgecell:
    #     cells = (allcells[0][:],) if left else (allcells[0][0:],)
    # else:
    #     cells = (allcells[0][:],) if left else (allcells[0][0:],)
    cells = (allcells[0][:],) if left else (allcells[0][0:],)
    # start block with essential properties (all dmat have these 2).
    nra = orig.radius[cells]
    nde = orig.density[cells]
    dblock = np.column_stack([nra, nde])
    keys = orig.filekeys
    rem = []
    for i, k in enumerate(keys):
        if k in orig._rnames or k in orig._dnames:
            continue
        elif k in orig._mnames:
            rem.append(i)
            continue
        dblock = np.column_stack((dblock, getattr(orig, k)[cells]))
    for index in rem:
        del keys[index]
    return DataMatrix([keys, dblock])


def multiIndex(names, keys):
    """returns the index of a key with multiple probable names in a
    list. (assumes such key only happens once)"""
    exists = False
    for n in names:
        try:
            return n, keys.index(n)
        except ValueError:
            continue
    return names[-1], -1


def chopFile(filename):
    """Returns header names and a data matrix from a file.

    Args:
        filename(str): file path.

    Returns:
        (list): lowercase header names.
        (np.array): data matrix of shape (coords, properties)

    """
    # exception loop to hack this with StringIO
    try:
        header = linecache.getline(filename, 1).strip(" #\n").split()
        data = np.genfromtxt(filename, skip_header=2, comments="#")
    except AttributeError:
        header = filename.readline()
        header = header.strip(" #\n").split()
        data = np.genfromtxt(filename, skip_header=1, comments="#")
    return [h.lower() for h in header], data

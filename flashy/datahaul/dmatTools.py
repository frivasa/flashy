"""various profile(dmat) 1D methods"""
import flashy.nuclear as nuc
import flashy.post as post
from flashy.utils import np, kb


def getSphericalVolumes(prof):
    """returns cell volumes times density for spherical cells.
    flash returns cell centers so there's a half-shift to get
    a correct volume.
    """
    fac = 4.0/3.0*np.pi
    rdiff = np.diff(prof.radius)
    rdiff = np.append(rdiff, rdiff[-1])
    skewed = 0.5*rdiff + prof.radius
    skewed = np.insert(skewed, 0, 0.0)
    rcub = skewed**3
    vols = fac*np.diff(rcub)
    return prof.density*vols


def getMassEnergy(prof, trueA=False):
    """calculate total binding energy per cell
    for a dmat(lineout) in erg/g.
    """
    # molar in mol/g
    factors = nuc.getBinding(prof.species, trueA=trueA)
    binding = []
    sumy = []
    for cell in range(len(prof.density)):
        # for each cell get all the mass fractions and convert to molar
        xis = prof.data[cell, len(prof.bulkprops)-1:]
        molar, abar, zbar = nuc.convXmass2Abun(prof.species, xis)
        sumy.append(np.sum(molar))
        binding.append(np.dot(factors, molar)*nuc.Avogadro)
    return sumy, binding  # binding in erg/g


def getProfDegen(prof, relativistic=False):
    """calculates degeneracy parameters for a dmat profile,
    namely Y_e and \eta = T / T_fermi.
    \eta near zero implies degeneracy,
    while \eta >> 1 or negative implies non-degenerate matter.

    Args:
        prof (dataMatrix): dataMatrix obj.

    Returns:
        (np.arrays): fermi temps, electron fractions, number fractions

    """
    # Ye and Yion
    yes, yis = [], []
    npnts = len(getattr(prof, prof.species[0]))
    for i in range(npnts):
        xis = []
        for s in prof.species:
            xis.append(getattr(prof, s)[i])
        invyi, invye = nuc.getMus(prof.species, xis)
        yes.append(1.0/invye)
        yis.append(1.0/invyi)
    # get fermi temperatures through the lineout
    if relativistic:
        fermT = [post.extRelFermi(d)/kb for d in prof.density]
    else:
        fermT = [post.nonRelFermi(d, ye=y)/kb
                 for d, y in zip(prof.density, yes)]
    return fermT, yes, yis

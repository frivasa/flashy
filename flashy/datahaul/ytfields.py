"""custom fields for yt.
method() describe operations over the ds
_method fills the metadata for yt:
units, take_log, sampling_type
"""
from flashy.utils import np
from flashy.nuclear import splitSpecies
_speed = {'units': 'cm/s', 'take_log': False, 'sampling_type': 'cell'}
_fermiDeg = {'units': 'auto', 'take_log': False,
             'sampling_type': 'cell', 'dimensions': 'dimensionless'}
# soundspeed already in derived field list: sound_speed
_machNumber = {'units': 'auto', 'take_log': False,
                'sampling_type': 'cell', 'dimensions': 'dimensionless'}
_burnLimiter = {'units': 'auto', 'take_log': False,
                'sampling_type': 'cell', 'dimensions': 'dimensionless'}
_pressure_scale_height = {'units': 'auto', 'take_log': False,
                          'sampling_type': 'cell', 'dimensions': 'dimensionless'}
_mean_molec_weight = {'units': 'auto', 'take_log': False,
                      'sampling_type': 'cell', 'dimensions': 'dimensionless'}
_cylgrav = {'units': 'cm/s/s', 'take_log': False, 'sampling_type': 'cell'}
_cylmass = {'units': 'g', 'take_log': False, 'sampling_type': 'cell'}
_cylvol = {'units': 'cm*cm*cm', 'take_log': False, 'sampling_type': 'cell'}
_alphas = ['he4 ', 'c12 ', 'o16 ', 'ne20',
           'mg24', 'si28', 's32 ', 'ar36',
           'ca40', 'ti44', 'cr48', 'fe52', 'ni56']
_massHe4 = {'units': 'g', 'take_log': False, 'sampling_type': 'cell'}
_massC12 = {'units': 'g', 'take_log': False, 'sampling_type': 'cell'}
_massO16 = {'units': 'g', 'take_log': False, 'sampling_type': 'cell'}
_massNe20 = {'units': 'g', 'take_log': False, 'sampling_type': 'cell'}
_massMg24 = {'units': 'g', 'take_log': False, 'sampling_type': 'cell'}
_massSi28 = {'units': 'g', 'take_log': False, 'sampling_type': 'cell'}
_massNi56 = {'units': 'g', 'take_log': False, 'sampling_type': 'cell'}


def speed(field, data):
    """velocity on grid magnitude"""
    vx = data['flash', 'velx']
    vy = data['flash', 'vely']
    vz = data['flash', 'velz']
    spd = np.sqrt(vx*vx + vy*vy, vz*vz)
    return spd


def fermiDeg(field, data):
    """ratio of the temperature on grid vs
    the fermi temperature (non-relativistic)
    """
    temp = data['flash', 'temp'].v
    dens = data['flash', 'dens'].v
    yes = data['flash', 'ye  '].v
    # local version for speed
    prefac1 = 6.1042642146064e-28  # 0.5*((0.5*h/pi)^2)/m_e
    prefac2 = 9.5707800006273  # (3pi^2)^(2/3)
    # prefac3: (avogadro*rho*Y_e) ^(2/3)
    prefac3 = np.power(6.02214076e23*dens*yes, 2.0/3.0)
    fermT = prefac1*prefac2*prefac3/1.380649e-16  # /k_b
    return temp/fermT


def machNumber(field, data):
    """gridspeed/soundspeed"""
    gc = data['flash', 'gamc'].v
    pr = data['flash', 'pres'].v
    de = data['flash', 'dens'].v
    cs = np.sqrt(gc*pr/de)
    vx = data['flash', 'velx'].v
    vy = data['flash', 'vely'].v
    vz = data['flash', 'velz'].v
    spd = np.sqrt(vx*vx + vy*vy, vz*vz)
    return spd/cs


def burnLimiter(field, data):
    """Kushnir, 2019 energy limiter factor.
    WARN: using dx as in first coordinate, this is not
    sound for non-square cells.
    # f_old = (ei/abs(enuc))/(dx/c_s)
    # enuc > 0 and
    # f_old < bn_burnLimiter triggers limiter

    """
    gc = data['flash', 'gamc'].v
    pr = data['flash', 'pres'].v
    de = data['flash', 'dens'].v
    cs = np.sqrt(gc*pr/de)
    ei = data['flash', 'eint'].v
    ec = data['flash', 'enuc'].v
    dx = data['flash', 'dx'].v
    return (ei/np.abs(ec))/(dx/cs)


def pressure_scale_height(field, data):
    """Scale height (H) for cylindrical coordinates.
    WARN: cylindrical_r is not sound for other coord systems.
    approximates gravity by U/r
    H = kT/mg
    m = mean molec weight = rho * N_A/(n_I+n_e)
    n_I = rho*N_A \Sigma_i x_i / A_i
    n_e = rho*N_A \Sigma_i x_i * Z_i / A_i
    """
    _, ps, ns, ms = splitSpecies(_alphas)
#     ms = np.array(len(_alphas)*[np.array(ms)])
#     ps = np.array(len(_alphas)*[np.array(ps)])
    ms = np.array(ms)
    ps = np.array(ps)
    xms = [data['flash', s].v for s in _alphas]
    xms = np.array(xms)
    # 13 cell*cell*cell blocks
    nions = [xms[i]/ms[i] for i in range(len(_alphas))]
    nions = np.array(nions)
    neles = [xms[i]/ms[i]*ps[i] for i in range(len(_alphas))]
    neles = np.array(neles)
    totalnele = np.sum(neles)
    totalnion = np.sum(nions)
    mu_ele = 1/(totalnele+totalnion)
    num = data['kT'].v
    g = -data['gpot'].v/data['cylindrical_r'].v
    return num/g/mu_ele


def mean_molec_weight(field, data):
    """Mean molecular weight for alpha mixture
    m = mean molec weight = rho * N_A/(n_I+n_e)
    n_I = rho*N_A \Sigma_i x_i / A_i
    n_e = rho*N_A \Sigma_i x_i * Z_i / A_i
    """
    _, ps, ns, ms = splitSpecies(_alphas)
    ms = np.array(ms)
    ps = np.array(ps)
    xms = [data['flash', s].v for s in _alphas]
    xms = np.array(xms)
    # 13 cell*cell*cell blocks
    nions = [xms[i]/ms[i] for i in range(len(_alphas))]
    nions = np.array(nions)
    neles = [xms[i]/ms[i]*ps[i] for i in range(len(_alphas))]
    neles = np.array(neles)
    totalnele = np.sum(neles)
    totalnion = np.sum(nions)
    return 1/(totalnele+totalnion)


def cylgrav(field, data):
    """gravity approximation in cylindrical
    g = -gpot/cylindrical_r
    """
    return -data['gpot']/data['cylindrical_r']


def cylmass(field, data):
    """mass of species"""
    dr = data['flash', 'path_element_x']
    dz = data['flash', 'path_element_y']
    r  = data['flash', 'x']
    de = data['flash', 'dens']
    return 2.0*np.pi*r*dr*dz*de


def cylvol(field, data):
    """mass of species"""
    dr = data['flash', 'path_element_x']
    dz = data['flash', 'path_element_y']
    r  = data['flash', 'x']
    return 2.0*np.pi*r*dr*dz


def massHe4(field, data):
    """mass of species"""
    xi = data['flash', 'he4 ']
    de = data['flash', 'dens']
    vl = data['flash', 'cell_volume']
    return xi*de*vl


def massC12(field, data):
    """mass of species"""
    xi = data['flash', 'c12 ']
    de = data['flash', 'dens']
    vl = data['flash', 'cell_volume']
    return xi*de*vl


def massO16(field, data):
    """mass of species"""
    xi = data['flash', 'o16 ']
    de = data['flash', 'dens']
    vl = data['flash', 'cell_volume']
    return xi*de*vl


def massNe20(field, data):
    """mass of species"""
    xi = data['flash', 'ne20']
    de = data['flash', 'dens']
    vl = data['flash', 'cell_volume']
    return xi*de*vl


def massMg24(field, data):
    """mass of species"""
    xi = data['flash', 'mg24']
    de = data['flash', 'dens']
    vl = data['flash', 'cell_volume']
    return xi*de*vl


def massSi28(field, data):
    """mass of species"""
    xi = data['flash', 'si28']
    de = data['flash', 'dens']
    vl = data['flash', 'cell_volume']
    return xi*de*vl


def massNi56(field, data):
    """mass of species"""
    xi = data['flash', 'ni56']
    de = data['flash', 'dens']
    vl = data['flash', 'cell_volume']
    return xi*de*vl

from _globals import *
import collections
import linecache

# custom fields
def _speed(field, data):
    vx = data['flash', 'velx']
    vy = data['flash', 'vely']
    spd = np.sqrt(vx*vx + vy*vy)
    return spd
#yt.add_field(("flash","speed"), function=_speed, units="cm/s", take_log=False)

# yt plotting functions
def fileProfile(fname, species=ap13, thresh=1e-6, mrange=[0.0, 0.0], radius=5e9, filetag='prof', show=False):
    """plots a checkpoint file ray through the domain.
    
    Args:
        fname (str): filename of checkpoint.
        species (list of str): list of species names to plot.
        thresh (float): ymin for species fraction plot.
        radius (float): reach of ray.
        mrange (list of float): if set, change the mass range of the abundance plot.
        filetag (str): change prefix of png output files.
        show (bool): return figure instead of saving it to file.
    
    """
    if 'chk' not in fname:
        print "Not a checkpoint, skipping abundance plots"
        plotsp = False
    else:
        plotsp = True
    ds = yt.load(fname)
    # spherical (r, theta, phi)
    cname = ds.coordinates.axis_name[0]
    ray = ds.ray([0e5, 0e5, 0] , [radius, 0, 0])
    ray_sort = np.argsort(ray['t'])
    #fig = plt.figure(figsize=(8, 7), dpi=80)
    fig = plt.figure(figsize=(12, 8))

    if plotsp:
        layout = (3, 2)
        ax4 = plt.subplot2grid(layout, (0, 1), aspect="auto", adjustable='box-forced', rowspan=2)
    else:
        layout = (3,1)
    ax1 = plt.subplot2grid(layout, (0, 0), aspect='auto')
    ax2 = plt.subplot2grid(layout, (1, 0), aspect="auto", adjustable='box-forced')
    ax3 = plt.subplot2grid(layout, (2, 0), aspect="auto", adjustable='box-forced')

    ax1.loglog(ray[cname][ray_sort], ray['density'][ray_sort], color='black')
    ax1.set_ylabel('Density($g/cm^3$)')
    ax1.set_xlabel('Radius ($cm$)')
    ax1.yaxis.set_major_formatter(StrMethodFormatter('{x:.2e}'))
    ax1.annotate("{:.5f} s".format(float(ds.current_time)),
                 xy=(0.0, 0.0), xytext=(0.05, 0.15), size=12,
                 textcoords='axes fraction', xycoords='axes fraction')

    ax2.loglog(ray[cname][ray_sort], ray['temperature'][ray_sort], color='red')
    ax2.set_ylabel('Temperature($K$)')
    ax2.set_xlabel('Radius ($cm$)')
    ax2.yaxis.set_major_formatter(StrMethodFormatter('{x:.2e}'))
    ax2.yaxis.set_minor_formatter(StrMethodFormatter(''))
    
    ax3.loglog(ray[cname][ray_sort], ray['pressure'][ray_sort], color='blue')
    ax3.set_ylabel('Pressure($dyne/cm^2$)')
    ax3.set_xlabel('Radius ($cm$)')
    ax3.yaxis.set_major_formatter(StrMethodFormatter('{x:.2e}'))
    
    if plotsp:
        # sort by atomic weight
        aws = []
        for s in species:
            if s=='n   ':
                aws.append(0)
            elif s=='p   ':
                aws.append(1)
            elif s=='d   ':
                aws.append(2)
            else:
                aws.append(int(''.join([c for c in s if c.isdigit()])))
        species = sorted(zip(aws, species))
        styleIter = colIter()
        # don't skip any plot to ensure colors stick to species, and legend doesn't 
        # shapeshift.
        xs = byMass(ray[cname][ray_sort], ray['density'][ray_sort])
        for s in [sp[1].strip() for sp in species]:
            tag = '$^{{{}}}{}$'.format(*elemSplit(s))
            c, ls = styleIter.next()
            ax4.semilogy(xs, ray[s][ray_sort], label=tag, color=c, linestyle=ls, alpha=0.7)
        lgd = ax4.legend(ncol=5, loc='upper left', bbox_to_anchor=(1.0, 1.0), 
          columnspacing=0.5, labelspacing=0.5, markerfirst=False, 
          numpoints=4)
        if np.sum(mrange)!=0.0:
            ax4.set_xlim(mrange)
        ax4.axhline(1e0, linewidth=1, linestyle=':', color='black')
        ax4.set_ylim(thresh, 2.0)
        ax4.set_xlabel('Mass ($M_{\odot}$)')
        ax4.set_ylabel('$X_{frac}$')
    plt.tight_layout(pad=1.0, h_pad=0.0, w_pad=0.5, rect=(0,0,0.67,1))
    if show:
        return
    else:
        num = ds.parameter_filename[-5:]
        otpf, _ = os.path.split(ds.fullpath)
        tag = filetag
        savn = '{}{}.png'.format(tag, num)
        savf = os.path.join(otpf, "png")
        savp = os.path.join(otpf, "png", tag, savn)
        # build filetree and show or save the figure
        if not os.path.exists(savf):
            os.mkdir(savf)
            os.mkdir(os.path.join(savf, tag))
        elif not os.path.exists(os.path.join(savf, tag)):
            os.mkdir(os.path.join(savf, tag))
        plt.savefig(savp,  bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.close(fig)
        print "Wrote: {}".format(savp)

# utility functions
def getVelocities(filename):
    """Returns positions of the shock, and both in an outer cj 
    velocities for a file, plus the starting x_match.
    (xin, xout, cjin, cjout, float(ray.ds.current_time), ray.ds.parameters['x_match'])
    """
    ray, ray_sort = getRay(filename)
    shockin, shockout = locateShock(ray, vvv=False)
    cname = ray.ds.coordinates.axis_name[0]
    xin, xout = ray[cname][ray_sort][shockin], ray[cname][ray_sort][shockout]
    cjin = roughCJ(ray, shockin)
    cjout = roughCJ(ray, shockout)
    return xin, xout, cjin, cjout, float(ray.ds.current_time), ray.ds.parameters['x_match']


def byMass(rads, dens):
    """Returns a mass abscissa for plots."""
    xs = len(rads)
    dr = rads[0]
    vol = dr**3 *4.0*np.pi/3.0
    mass = vol*dens[0]/Ms
    absc = []
    absc.append(mass)
    for i in range(1, xs):
        dr = rads[i] - rads[i-1]
        dvol = dr * ( 3.0*rads[i-1]*rads[i] + dr*dr ) * 4.0*np.pi/3.0
        mass = mass + dvol*dens[i]/Ms
        absc.append(mass)
    return absc


def roughCJ(ray, point):
    ray_sort = np.argsort(ray['t']) # not time, this is a grid variable for ordering points
    dp = float(ray['pressure'][ray_sort][point-1]-ray['pressure'][ray_sort][point+1])
    nu1 = 1.0/float(ray['density'][ray_sort][point+1])
    nu2 = 1.0/float(ray['density'][ray_sort][point-1])
    dnu = (nu1-nu2)/(nu1**2)
    # print dp, nu1, nu2, dnu
    if dnu==0.0:
        dnu=1.0
    svel = np.sqrt(dp/dnu)
    return svel


def locateShock(ray, vvv=True):
    """returns index within ray of detected shock."""
    ray_sort = np.argsort(ray['t'])
    cname = ray.ds.coordinates.axis_name[0]
    
    filt, offs1 = split(ray[cname][ray_sort], ray.ds.parameters['x_match'], True)
    shockin = shock1D(ray[cname][ray_sort][filt], ray['sound_speed'][ray_sort][filt][:-1], True)
    filt, offs2 = split(ray[cname][ray_sort], ray.ds.parameters['x_match'], False)
    shockout = shock1D(ray[cname][ray_sort][filt], ray['sound_speed'][ray_sort][filt][:-1], False)
    if vvv:
        print 'Ignition Center: ', ray.ds.parameters['x_match']
        print 'Inward Shock at: {:E}'.format(float(ray[cname][ray_sort][shockin+offs1]))
        print 'Outward Shock at: {:E}'.format(float(ray[cname][ray_sort][shockout+offs2]))
    return shockin+offs1, shockout+offs2


def shock1D(rad, soundspeeds, inward=True):
    """finds a shock in an array by detecting the last 
    large variation within it that is larger than the mean of deltas.
    """
    dr = np.diff(rad)[:-1]
    ds = np.diff(soundspeeds)
    div = np.nan_to_num(ds/dr)
    accel = np.abs(np.nan_to_num(ds/dr))
    mean = np.mean(accel)
    pertp = np.where(accel>mean)
    if inward:
        return pertp[0][0]
    else:
        return pertp[0][-1]


def getRay(fname, radius=5e9):
    """returns a ray trace through the domain."""
    ds = yt.load(fname)
    ray = ds.ray([0e5, 0e5, 0] , [radius, 0, 0])
    ray_sort = np.argsort(ray['t'])
    return ray, ray_sort


def customFormatter(factor, prec=1, width=2):
    """specialized format for plotting labels:
    width(prec) x 10^factor.
    """
    fstr = '{:{width}.{prec}f}'
    exp = 10.0**factor
    return FuncFormatter(lambda x, pos:fstr.format(x/exp, 
                                                   width=width, 
                                                   prec=prec))


def getExtrema(filename, flist=['density', 'temperature', 'pressure']):
    """returns a list of tuples with extrema for given fields(flist)."""
    ds = yt.load(filename)
    ad = ds.all_data()
    return ad.quantities.extrema(flist)


def gridInfo(filename, silent=True):
    """prints checkpoint file grid stats 
    through yt.print_stats() method.
    """
    ds = yt.load(filename)
    if not silent:
        ds.print_stats()
    return ds.parameters, float(ds.current_time)


def writeProf(filename, profD):
    """writes a 1D profile dict to a file."""
    points = len(profD['Radius'])
    header = " ".join(profD.keys())
    with open(filename, 'w') as f:
        f.write('# {}\n'.format(header))
        f.write('{}\n'.format(points))
        for i in range(points):
            line = []
            for k in profD.keys():
                line.append('{:15.8e}'.format(profD[k][i]))
            f.write("{}\n".format(" ".join(line)))
    print 'Wrote: {}'.format(filename)


def getProps(profD):
    """returns central dens, temp, rmin, rmax, volume and mass of a 1D profile."""
    pts = len(profD['Radius'])
    mass, vol = 0 , 0
    rads = np.insert(np.copy(profD['Radius']), 0, 0.0)
    for i in range(1, pts):
        r2, r1 = profD['Radius'][i], profD['Radius'][i-1]
        dr = r2-r1
        dvol = dr * ( 3.0*r1*r2 + dr*dr ) * 4.0*np.pi/3.0
        vol += dvol
        mass += dvol*profD['dens'][i-1]
    cd, ct = profD['dens'][0], profD['temp'][0]
    msuns = mass/(1.989e33)
    res = (r2 - profD['Radius'][0])/pts
    print "Points/resolution: {} / {:E}".format(pts, res)
    print "Central Dens: {:E}".format(cd)
    print "Central Temp: {:E}".format(ct)
    print "Radius: {:E}".format(r2)
    print "Volume: {:E}".format(vol)
    print "Mass (Msun): {:E}".format(msuns)
    return pts, res, cd, ct, r2, vol, msuns

# auxiliary functions
def getProfData(filename):
    """Returns arrays for radius, dens, temp and a matrix with abundances, 
    as an ordered dictionary.
    
    """
    data = np.genfromtxt(filename, skip_header=2, comments='#')
    header = linecache.getline(filename, 1).strip(' #\n').split()
    varsinFile = collections.OrderedDict(zip(header, range(len(header))))
    outDict = collections.OrderedDict()
    for k, v in varsinFile.items():
        var = data[:, v]
        outDict[k] = var
    return outDict


def plotter(ax, filename, key, xkey='Radius', label='', log=1, **kwargs):
    """plots 'key' values from a data dictionary"""
    dataDict = getProfData(filename)
    ax.set_ylabel(key)
    ax.set_xlabel(xkey)
    if log:
        ax.loglog(dataDict[xkey], dataDict[key], label=label, **kwargs)
    else:
        ax.plot(dataDict[xkey], dataDict[key], label=label, **kwargs)


def percentDiff(ax, file1, file2, diffkey, xkey='Radius', label='', **kw):
    """plots the percentage difference for diffkey for two filenames."""
    d1 = getProfData(file1)
    d2 = getProfData(file2)
    jake = np.interp(d1[xkey], d2[xkey], d2[diffkey])
    norm = np.max(jake)
    ax.plot(d1[xkey], 100*(d1[diffkey]-jake)/norm, label=label, **kw)


def getDset(filename):
    """Returns the yt dataset (not compatible with batch/parallel scripts)"""
    return yt.load(filename)
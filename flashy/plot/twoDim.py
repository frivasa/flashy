import os
import numpy as np
import yt
# avoid yt warnings
from yt.funcs import mylog
mylog.setLevel(50)
# mpl config
import matplotlib as mpl
nonguis = [u'agg', u'cairo', u'gdk', u'pdf', u'pgf', u'ps', u'svg', u'template']
for nongui in nonguis:
    try:
        #print "[flash_plot]: Testing", nongui
        mpl.use(nongui, warn=False, force=True)
        from matplotlib import pyplot as plt
        break
    except:
        continue
print "[flash.plot.twoD]: Using",mpl.get_backend()
from mpl_toolkits.axes_grid1 import AxesGrid
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from matplotlib.ticker import FuncFormatter, StrMethodFormatter
from matplotlib.colors import LogNorm, NoNorm
mpl.rc('lines', linewidth=2, linestyle='-', marker=None)
mpl.rc('font', family='monospace', size=12.0)
mpl.rc('text', color='000000')
mpl.rc('axes', linewidth=2, grid=False, titlepad=10.0, labelsize='large', 
       axisbelow=False, autolimit_mode='data')  # round_numbers
mpl.rc('axes.formatter', limits=(-1,1))
mpl.rc('xtick', top=True, direction='in')
mpl.rc('xtick.major', size=9, width=1.4, pad=7)
mpl.rc('xtick.minor', size=5, width=1.0, pad=7)
mpl.rc('ytick', right=True, direction='in')
mpl.rc('ytick.major', size=9, width=1.4, pad=7)
mpl.rc('ytick.minor', size=5, width=1.0, pad=7)
mpl.rc('savefig', facecolor='ffffff', dpi=100, bbox='tight')
# COLORMAPS
# https://jiffyclub.github.io/palettable/colorbrewer/
import palettable  # cmocean colors are included in palettable
# palettable schemes are not interpolated so make a new one for yt
# yt call: (name, type, num): num<=11 type={diverging, qualitative, sequential}
#_cScheme = ('GnBu', 'sequential', 9)
# import cmocean fails for ipyparallel. so hack it through palettable
_cmap = palettable.cmocean.sequential.Deep_20.mpl_colormap
# import cmocean fails for ipyparallel. 
# so hack it through palettable
#_ytcmap = palettable.cmocean.sequential.Ice_14_r.mpl_colors
_ytcmap = palettable.matplotlib.Inferno_20_r.mpl_colors
cols = [tuple([list(x),1]) for x in _ytcmap]
cols.append(([0.0, 0.0, 0.0], 0))
setcmap = yt.make_colormap(cols, name='custom')
# CONSTANTS AND SPECIES
_xnet = [
    'n   ', 'p   ', 'd   ', 'he3 ', 'he4 ', 'li6 ', 'li7 ', 'be7 ', 'be9 ', 
    'b8  ', 'b10 ', 'b11 ', 'c12 ', 'c13 ', 'c14 ', 'n13 ', 'n14 ', 'n15 ', 
    'o14 ', 'o15 ', 'o16 ', 'o17 ', 'o18 ', 'f17 ', 'f18 ', 'f19 ', 'ne18', 
    'ne19', 'ne20', 'ne21', 'ne22', 'na21', 'na22', 'na23', 'mg23', 'mg24', 
    'mg25', 'mg26', 'al24', 'al25', 'al26', 'al27', 'si28', 'si29', 'si30', 
    'si31', 'si32', 'p29 ', 'p30 ', 'p31 ', 'p32 ', 'p33 ', 's32 ', 's33 ', 
    's34 ', 's35 ', 's36 ', 'cl33', 'cl34', 'cl35', 'cl36', 'cl37', 'ar36', 
    'ar37', 'ar38', 'ar39', 'ar40', 'k37 ', 'k38 ', 'k39 ', 'k40 ', 'k41 ', 
    'ca40', 'ca41', 'ca42', 'ca43', 'ca44', 'ca45', 'ca46', 'ca47', 'ca48', 
    'sc43', 'sc44', 'sc45', 'sc46', 'sc47', 'sc48', 'sc49', 'ti44', 'ti45', 
    'ti46', 'ti47', 'ti48', 'ti49', 'ti50', 'v46 ', 'v47 ', 'v48 ', 'v49 ', 
    'v50 ', 'v51 ', 'cr48', 'cr49', 'cr50', 'cr51', 'cr52', 'cr53', 'cr54', 
    'mn50', 'mn51', 'mn52', 'mn53', 'mn54', 'mn55', 'fe52', 'fe53', 'fe54', 
    'fe55', 'fe56', 'fe57', 'fe58', 'co53', 'co54', 'co55', 'co56', 'co57', 
    'co58', 'co59', 'ni56', 'ni57', 'ni58', 'ni59', 'ni60', 'ni61', 'ni62', 
    'cu57', 'cu58', 'cu59', 'cu60', 'cu61', 'cu62', 'cu63', 'zn59', 'zn60',
    'zn61', 'zn62', 'zn63', 'zn64', 'zn65', 'zn66' 
]
_xnetReduced = [
    'neut', 'h1  ', 'h2  ', 'he3 ', 'he4 ', 'li7 ', 'b8  ', 'c12 ', 'n14 ',
    'o16 ', 'f18 ', 'ne20', 'na22', 'na23', 'mg24', 'al26', 'al27', 'si28', 
    'p30 ', 's32 ', 'cl34', 'ar36', 'ar40', 'k38 ', 'ca40', 'ca44', 'sc44',
    'ti44', 'ti45', 'ti46', 'v46 ', 'cr48', 'cr50', 'mn50', 'fe52', 'fe53', 
    'fe54', 'fe55', 'fe56', 'fe57', 'fe58', 'co54', 'co55', 'co56', 
    'co54', 'ni56', 'ni58', 'cu58', 'zn60'
]
_ap13 = [
    'he4 ', 'c12 ', 'o16 ', 'ne20',
    'mg24', 'si28', 's32 ', 'ar36',
    'ca40', 'ti44', 'cr48', 'fe52', 'ni56'
]
_Rs = 695700e5 # cm
_Ms = 1.988e33 # g
_plR = 1.1e9 # cm, sphere radius cutoff

# custom fields
def _speed(field, data):
    vx = data['flash', 'velx']
    vy = data['flash', 'vely']
    spd = np.sqrt(vx*vx + vy*vy)
    return spd
#yt.add_field(("flash","speed"), function=_speed, units="cm/s", take_log=False)

# yt plotting functions
def planeSlice(fname, field, lims, zcut=0.0, linear=False, show=False, width=_plR):
    """Makes a slice at zcut"""
    ds = yt.load(fname)
    p = yt.SlicePlot(ds, 'z', field, center=[0, 0, zcut])
    p.set_width((width, width))
    p.set_cmap(field, 'custom')
    p.set_axes_unit('cm')
    fig = plt.figure(figsize=(10,10))
    grid = AxesGrid(fig, 111, nrows_ncols = (1, 1),
                    axes_pad = 0.0, label_mode = "L",
                    share_all = False, cbar_location="right",
                    cbar_mode="each", cbar_size="5%", cbar_pad="0%")
    im, cax = plotFRB(grid[0], grid.cbar_axes[0], 
                      np.transpose(p.frb[field]), 
                      lims, top=False)
    im.axes.tick_params(axis='x', top='on')
    im.axes.annotate("Time: {:.5f} s".format(float(ds.current_time)),
                     xy=(0.75, 0.05), textcoords='figure fraction')
    im.axes.annotate("Z: {:.4e} cm".format(zcut),
                     xy=(0.10, 0.05), textcoords='figure fraction')
    if show:
        return fig
    else:
        num = ds.parameter_filename[-5:]
        otpf, _ = os.path.split(ds.fullpath)
        tag = 'slices{}'.format(num)
        savn = 'slice_{:+13.0f}.png'.format(zcut)
        savf = os.path.join(otpf, "png")
        savp = os.path.join(otpf, "png", tag, savn)
        # build filetree and show or save the figure
        if not os.path.exists(savf):
            os.mkdir(savf)
            os.mkdir(os.path.join(savf, tag))
        elif not os.path.exists(os.path.join(savf, tag)):
            os.mkdir(os.path.join(savf, tag))
        plt.savefig(savp)
        plt.close(fig)
        print "Wrote: {}".format(savp)


def bifold(fname, mhead=True, topfield='density', tlims=[1e-1, 4e7], width=_plR,
           botfield='temperature', blims=[3e7, 2e9], lintop=False,  show=False,
           linbot=False, name=''):
    """Plots 2 properties as a joined sphere."""
    ds = yt.load(fname)
    p = yt.SlicePlot(ds, 'z', [topfield, botfield])
    p.set_width((width, 2*width))
    p.set_center((width*0.5, 0.0))
    p.set_origin(("center", "left", "domain"))
    p.set_axes_unit('cm')
    fig = plt.figure(figsize=(12,9))
    #grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
    grid = AxesGrid(fig, 121,
                    nrows_ncols = (2, 1),
                    axes_pad = 0.0, label_mode = "L",
                    share_all = False, cbar_location="right",
                    cbar_mode="each", cbar_size="5%", cbar_pad="0%")
    plotFRB(grid[0], grid.cbar_axes[0], 
                  np.flipud(np.transpose(p.frb[topfield])), 
                  tlims, linear=lintop)
    plotFRB(grid[1], grid.cbar_axes[1], 
                  np.fliplr(np.transpose(p.frb[botfield])), 
                  blims, 
                  top=False, linear=linbot)
    if mhead:
        # we've transposed the matrix, so flip the _match positions
        x_match = ds.parameters['y_match']
        y_match = ds.parameters['x_match']
        grid[0].axes.plot(x_match, y_match, 'o', color='black')
        grid[1].axes.plot(x_match, y_match, 'o', color='black')
    grid[1].axes.annotate("Time: {:.5f} s".format(float(ds.current_time)),
                          xy=(0.75, -0.17), textcoords='axes fraction')
    fig.tight_layout()
    if show:
        return fig
    else:
        num = ds.parameter_filename[-5:]
        otpf, _ = os.path.split(ds.fullpath)
        if name:
            tag = name
        else:
            tag = 'bifold'
        savn = '{}{}.png'.format(tag, num)
        savf = os.path.join(otpf, "png")
        savp = os.path.join(otpf, "png", tag, savn)
        # build filetree and show or save the figure
        if not os.path.exists(savf):
            os.mkdir(savf)
            os.mkdir(os.path.join(savf, tag))
        elif not os.path.exists(os.path.join(savf, tag)):
            os.mkdir(os.path.join(savf, tag))
        plt.savefig(savp)
        plt.close(fig)
        print "Wrote: {}".format(savp)


def plotFRB(gridAx, cbgAx, imArr, lims, top=True, linear=False):
    unit = imArr.units
    field = imArr.info['field']
    yl, yr = imArr.info['ylim']
    xl, xr = imArr.info['xlim']
    vmin, vmax = lims
    if linear:
        norm = NoNorm(vmin=vmin, vmax=vmax)
    else:
        norm = LogNorm(vmin=vmin, vmax=vmax)
    if top:
        im = gridAx.imshow(imArr, 
                           cmap=_cmap, aspect='auto', 
                           extent=[yl, yr, xl, xr], # horizontal mi/ma vertical mi/ma
                           norm=norm)
        yticks = im.axes.yaxis.get_major_ticks()
        yticks[0].label.set_visible(False)
        im.axes.tick_params(axis='x', bottom='off')
    else:
        im = gridAx.imshow(imArr, 
                           cmap=_cmap, aspect='auto', 
                           extent=[yr, yl, xr, xl], norm=norm)
        im.axes.tick_params(axis='x', top='off')
    
    # get order of mag of highest y value and rescale the ticks
    #print im.get_extent()
    scaling = int(np.log10(abs(xr.value)))
    #scaling = int(np.log10(abs(im.get_extent()[3].value)))
    im.axes.yaxis.set_major_formatter(customFormatter(scaling))
    # yt flips axes, plus imarray is transposed, this must be an hdf5 saving issue
    im.axes.set_ylabel(r'x ($10^{{{}}}$ cm)'.format(scaling))
    
    # get order of mag of highest x value and rescale the ticks
    #scaling = int(np.log10(abs(im.get_extent()[1].value)))
    scaling = int(np.log10(abs(yr.value)))
    im.axes.xaxis.set_major_formatter(customFormatter(scaling))
    # yt flips axes, plus imarray is transposed, this must be an hdf5 saving issue
    im.axes.set_xlabel(r'y ($10^{{{}}}$ cm)'.format(scaling))
    
    scaling = int(np.log10(max(im.get_clim())))
    cax = cbgAx.colorbar(im, format=customFormatter(scaling, width=0, prec=2))
    cax.set_label_text('{} ($10^{{{}}} {}$)'.format(field.capitalize(),
                                                scaling,  
                                                str(unit).replace('**', '^')),
                        fontsize=11.0)
    yticks = cax.ax.get_yticklabels()
    yticks[-1].set_visible(False)
    yticks[0].set_visible(False)
    return im, cax


def mainProps(fname, mhead=True, grids=False, show=False,
              fields=['density', 'pressure', 'temperature'], linear=False, 
              mins=[0.1, 1e+16, 3e7], maxs=[4e7, 6e+24, 2e9]):
    """Plots the list of fields specified in yt.
    
    Args:
        fname (str): filename to plot.
        mhead (bool): mark the position of the matchhead.
        grids (bool): overplot the grid structure.
        show (bool): return figure (true) or save to file (false)
        fields (list of str): list of named fields to plot.
        linear (bool): set linear or log scale(false).
        mins (list of float): minima of scale for each field.
        maxs (list of float): maxima of scale for each field.
    
    """
    ds = yt.load(fname)
    size = len(fields)
    fig = plt.figure(figsize=(5*size, 10))
    grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                    nrows_ncols = (1, size),
                    axes_pad = 1.2, label_mode = "L",
                    share_all = True, cbar_location="right",
                    cbar_mode="each", cbar_size="10%", cbar_pad="0%")
    p = yt.SlicePlot(ds, 'z', list(fields))
    p.set_width((_plR, 2*_plR))
    p.set_center((_plR*0.5, 0.0))
    p.set_origin(("center", "left", "domain"))
    p.set_axes_unit('cm')
    header = '{:17.6f} s'
    if mhead:
        x_match = ds.parameters['x_match']
        y_match = ds.parameters['y_match']
        p.annotate_marker((x_match, y_match), coord_system='plot',
                          plot_args={'color':'black', 's': 30})
    if grids:
        p.annotate_grids()
    pvars = zip(fields, mins, maxs)
    for i, (f, mi, mx) in enumerate(pvars):
        if not i:
            p.annotate_title(header.format(float(ds.current_time)))
        p.set_cmap(f, 'custom')
        p.set_zlim(f, mi, mx)
        if linear:
            p.set_log(field, False)
        plot = p.plots[f]
       # plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
    p._setup_plots()
    if show:
        return fig
    else:
        num = ds.parameter_filename[-5:]
        otpf, _ = os.path.split(ds.fullpath)
        tag = 'ytprops'
        savn = '{}{}.png'.format(tag, num)
        savf = os.path.join(otpf, "png")
        savp = os.path.join(otpf, "png", tag, savn)
        # build filetree and show or save the figure
        if not os.path.exists(savf):
            os.mkdir(savf)
            os.mkdir(os.path.join(savf, tag))
        elif not os.path.exists(os.path.join(savf, tag)):
            os.mkdir(os.path.join(savf, tag))
        plt.savefig(savp)
        plt.close(fig)
        print "Wrote: {}".format(savp)


def fileProfile(fname, species=_ap13, thresh=1e-6, angle=45, radius=5e9, show=False):
    """plots a checkpoint file ray through the domain.
    
    Args:
        fname (str): filename of checkpoint.
        species (list of str): list of species names to plot.
        thresh (float): ymin for species fraction plot.
        angle (float): angle of ray in degrees (measured from +y towards +x).
        radius (float): reach of ray.
        show (bool): return figure instead of saving it to file.
    
    """
    if 'chk' not in fname:
        print "Not a checkpoint, skipping abundance plots"
        plotsp = False
    else:
        plotsp = True
    ds = yt.load(fname)
    ray = ds.ray([0e5, 0e5, 0] , [radius*np.sin(np.radians(angle)), radius*np.cos(np.radians(angle)), 0])
    ray_sort = np.argsort(ray['t'])
    fig = plt.figure(figsize=(16, 8))
    if plotsp:
        layout = (3, 3) # extra plot space for legend spacing
        ax4 = plt.subplot2grid(layout, (0, 1), aspect="auto", adjustable='box-forced', rowspan=2)
    else:
        layout = (3,1)
    ax1 = plt.subplot2grid(layout, (0, 0), aspect='auto')
    ax2 = plt.subplot2grid(layout, (1, 0), aspect="auto", adjustable='box-forced')
    ax3 = plt.subplot2grid(layout, (2, 0), aspect="auto", adjustable='box-forced')

    ax1.loglog(ray['radius'][ray_sort], ray['density'][ray_sort], color='black')
    ax1.set_ylabel('Density($g/cm^3$)')
    ax1.set_xlabel('Radius ($cm$)')
    ax1.yaxis.set_major_formatter(StrMethodFormatter('{x:.2e}'))
    ax1.annotate("{:.5f} s".format(float(ds.current_time)),
                 xy=(0.0, 0.0), xytext=(0.05, 0.15), size=12,
                 textcoords='axes fraction', xycoords='axes fraction')

    ax2.loglog(ray['radius'][ray_sort], ray['temperature'][ray_sort], color='red')
    ax2.set_ylabel('Temperature($K$)')
    ax2.set_xlabel('Radius ($cm$)')
    ax2.yaxis.set_major_formatter(StrMethodFormatter('{x:.2e}'))
    ax2.yaxis.set_minor_formatter(StrMethodFormatter(''))
    
    ax3.loglog(ray['radius'][ray_sort], ray['pressure'][ray_sort], color='blue')
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
        for s in [sp[1] for sp in species]:
            tag = '$^{{{}}}{}$'.format(*elemSplit(s))
            c, ls = styleIter.next()
            ax4.loglog(ray['radius'], ray[s][ray_sort], label=tag, color=c, linestyle=ls, alpha=0.7)
        ax4.legend(ncol=5, loc='upper left', bbox_to_anchor=(1.0, 1.0), 
          columnspacing=0.5, labelspacing=0.5, markerfirst=False, 
          numpoints=4)
        ax4.axhline(1e0, linewidth=1, linestyle=':', color='black')
        #ax4.legend(bbox_to_anchor=(1.0, 1.0), ncol=2)
        ax4.set_ylim(thresh, 2.0)
        ax4.set_xlabel('Radius ($cm$)')
        ax4.set_ylabel('$X_{frac}$')
    plt.tight_layout()
    if show:
        return fig
    else:
        num = ds.parameter_filename[-5:]
        otpf, _ = os.path.split(ds.fullpath)
        tag = 'prof'
        savn = '{}{}.png'.format(tag, num)
        savf = os.path.join(otpf, "png")
        savp = os.path.join(otpf, "png", tag, savn)
        # build filetree and show or save the figure
        if not os.path.exists(savf):
            os.mkdir(savf)
            os.mkdir(os.path.join(savf, tag))
        elif not os.path.exists(os.path.join(savf, tag)):
            os.mkdir(os.path.join(savf, tag))
        plt.savefig(savp)
        plt.close(fig)
        print "Wrote: {}".format(savp)

# utility functions
def roughCJ(ray, point):
    ray_sort = np.argsort(ray['t']) # not time, this is a grid variable for ordering points
    dp = float(ray['pressure'][ray_sort][point-1]-ray['pressure'][ray_sort][point+1])
    nu1 = 1.0/float(ray['density'][ray_sort][point+1])
    nu2 = 1.0/float(ray['density'][ray_sort][point-1])
    dnu = (nu1-nu2)/(nu1**2)
    print dp, nu1, nu2, dnu
    svel = np.sqrt(dp/dnu)
    return svel


def locateShock(ray, vvv=True):
    """returns index within ray of detected shock."""
    ray_sort = np.argsort(ray['t'])
    # build reference vectors
    igr = np.array((ray.ds.parameters['x_match'], ray.ds.parameters['y_match'], ray.ds.parameters['z_match']))
    r = np.array((float(ray['x'][ray_sort][-1]), float(ray['y'][ray_sort][-1]), float(ray['z'][ray_sort][-1])))
    cross = np.linalg.norm(np.cross(igr,r))
    dot = np.linalg.norm(np.dot(igr, r))
    angle = np.arctan(cross/dot)

    spx = np.cos(angle)*np.linalg.norm(igr)
    filt, offs1 = split(ray['radius'][ray_sort], spx, True)
    shockin = shock1D(ray['radius'][ray_sort][filt], ray['sound_speed'][ray_sort][filt][:-1], True)
    filt, offs2 = split(ray['radius'][ray_sort], spx, False)
    shockout = shock1D(ray['radius'][ray_sort][filt], ray['sound_speed'][ray_sort][filt][:-1], False)
    if vvv:
        print 'Ray: ', r
        print 'Ignition Center: ', igr
        print 'Angle: ', angle
        print 'Inward Shock at: {:E}'.format(float(ray['radius'][ray_sort][shockin+offs1]))
        print 'Outward Shock at: {:E}'.format(float(ray['radius'][ray_sort][shockout+offs2]))
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


def getRay(fname, direction=[0.0, 1.0, 0.0], radius=5e9):
    """returns a ray trace through the domain."""
    ds = yt.load(fname)
    ray = ds.ray([0e5, 0e5, 0] , np.array(direction)*radius)
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

# auxiliary functions
def getDset(filename):
    """Returns the yt dataset (not compatible with batch/parallel scripts)"""
    return yt.load(filename)


def split(x, xsplit, inward=True):
    """returns indices below or above xsplit and offset
    [0,1,2,3,4,5,6]
    inward True, xsplit 3: [0,1,2], 0
    inward False, xsplit 3: [4,5,6], 4
    """
    if inward:
        return np.where(x<xsplit), 0
    else:
        filt = np.where(x>xsplit)
        return filt, filt[0][0]


def colIter2():
    """Simple color iterator. Colors selected from Sasha Trubetskoy's 
    simple 20 color list (based on metro lines).
    https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
    """
    cols = ['#e6194b', '#3cb44b', '#0082c8', '#000000', '#f58231', '#911eb4', 
            '#008080', '#e6beff', '#aa6e28', '#fffac8', '#800000', '#aaffc3', 
            '#808000', '#ffd8b1', '#000080', '#808080', '#ffe119', '#f032e6', 
            '#46f0f0', '#d2f53c', '#fabebe']
    i=-1
    while(True):
        i+=1
        if i==len(cols):
            i=0
        yield cols[i]


def colIter():
    """Simple color/linestyle iterator. Colors selected from Sasha 
    Trubetskoy's simple 20 color list (based on metro lines)
    https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
    """
    cols = ['#e6194b', '#3cb44b', '#0082c8', '#000000', '#f58231', '#911eb4', 
            '#008080', '#e6beff', '#aa6e28', '#999678', '#800000', '#88cc9c', 
            '#808000', '#ffb265', '#000080', '#fabebe', '#e5c700']
            #, '#f032e6', '#46f0f0', '#d2f53c', '#808080']ffe119
    styles = [(0, ()),
              # dotted loose/normal/dense
              (0, (1, 10)),
              #(0, (1, 5)),
              (0, (1, 1)),
              # dashed loose/normal/dense
              (0, (5, 10)), 
              #(0, (5, 5)),
              (0, (5, 1)),
              # dash-dot loose/normal/dense
              (0, (3, 10, 1, 10)), 
              #(0, (3, 5, 1, 5)),
              (0, (3, 1, 1, 1)),
              # dash-dot-dot loose/normal/dense
              (0, (3, 10, 1, 10, 1, 10)),
              #(0, (3, 5, 1, 5, 1, 5))] 
              (0, (3, 1, 1, 1, 1, 1))]
    lstyles = len(styles)
    lcols = len(cols)
    alphas = np.linspace(0.0, 1.0, num=lstyles)
    i, j = -1, 0
    while(True):
        i+=1
        if i==lcols:
            i=0
            j+=1
        yield cols[i], styles[j]#, alphas[i]
        #if i==lcols:
        #    i=0
        #yield cols[i], styles[i%lstyles]


def elemSplit(s):
    """Standalone element name spliter. 
    he4 -> (4, He)
    
    Args:
        s(str): element string of the form "numberSpecies".
    
    Returns:
        tuple of str: (element number, capitalized name)
    
    """
    sym = s.rstrip('0123456789 ')
    A = s[len(sym):].strip()
    return A, sym.title()


# batch functions
def writeBatch(filenames, runf, procs, method, args):
    sp = os.path.dirname(os.path.realpath(__file__))
    runf = os.path.abspath(runf)
    with open(os.path.join(runf, "batch_{}".format(method)), 'w') as f:
        for i, file in enumerate(filenames, 1):
            _, name = os.path.split(file)
            line = ['python {}'.format(os.path.join(sp, "flash_plot.py"))]
            line.append(method)
            line.append(os.path.join(runf, "otp", name))
            for a in args:
                line.append(str(a))
            line.append("&\n")
            f.write(" ".join(line))
            if i%procs==0:
                f.write("wait\n")


if __name__=="__main__":
    import sys
    import os
    methods = {'planeSlice': planeSlice, 
               'bifold': bifold, 
               'fileProfile': fileProfile,
               'mainProps': mainProps}
    method = sys.argv[1]
    if method in methods:
        if method=='fileProfile':
            fname, thresh = sys.argv[2:]
            folder = 'xProfile'
            fig = methods[method](fname, thresh=thresh)
        elif method=='bifold':
            fname, tfield, mi, ma, lin = sys.argv[2:]
            if len(tfield)<4:
                tfield = tfield.ljust(4)
            folder = tfield
            fig = methods[method](fname, topfield=tfield, 
                                  tlims=[float(mi), float(ma)],
                                  lintop=int(lin))
        else:
            print "Fill with more methods"
    else:
        raise Exception("Method {} is not implemented".format(method))
    otp, name = os.path.split(fname)
    runf, _ = os.path.split(otp)
    if not os.path.exists(os.path.join(runf, "png", folder)):
        os.mkdir(os.path.join(runf, "png", folder))
    plt.savefig('{}.png'.format(os.path.join(runf, "png", folder, name)))
    print "Wrote: {}.png".format(os.path.join(runf, "png", folder, name))
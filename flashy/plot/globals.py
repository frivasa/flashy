from ..IOutils import os
from ..utils import np
import pkg_resources
import matplotlib as mpl
# backends for script/parallel ploting
nonguis = [u'agg', u'cairo', u'gdk', u'pdf', u'pgf', u'ps', u'svg', u'template']
for nongui in nonguis:
    try:
        #print "[flash_plot]: Testing", nongui
        mpl.use(nongui, warn=False, force=True)
        from matplotlib import pyplot as plt
        break
    except:
        continue
# print("[flash.plot]: Using",mpl.get_backend())

# misc imports
from mpl_toolkits.axes_grid1 import AxesGrid
# from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from matplotlib.ticker import FuncFormatter, StrMethodFormatter, ScalarFormatter
from matplotlib.colors import LogNorm, NoNorm

# x = {1, 5, 10}
# (0, (1, x)) dotted
# (0, (5, x)) dashed
# (0, (3, x, 1, x)) dash-dot
# (0, (1, x, 1, x, 1, x)) dash-dot-dot
from cycler import cycler
lines = [(0, ()), 
         (0, (1, 5)),
         (0, (5, 5)),
         (0, (3, 5, 1, 5)),
         (0, (3, 5, 1, 5, 1, 5)),
         (0, (1, 1)),
         (0, (5, 1)),
         (0, (3, 1, 1, 1)),
         (0, (3, 1, 1, 1, 1, 1))]

# Colors modified from Sasha Trubetskoy's 
# simple 20 color list (based on metro lines).
# https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
colors = ['#e6194b', '#3cb44b', '#0082c8', '#000000', '#f58231', 
          '#911eb4', '#008080', '#e6beff', '#bddc36', '#ccc8a0', 
          '#800000', '#808080', '#808000', '#46f0f0', '#000080', 
          '#ffe119', '#aa6e28', '#f032e6', '#ffa64d', '#7fbf92', '#f68888']
cc = (cycler('linestyle', lines)*cycler('color', colors))

# styling
style = pkg_resources.resource_filename('flashy', '../data/mplstyles/flashydef.mplsty')
# mpl.rc('lines', linewidth=2, linestyle='-', marker=None)
# mpl.rc('font', family='monospace', size=12.0)
# mpl.rc('text', color='000000')
# mpl.rc('axes', linewidth=2, grid=False, labelsize='large', axisbelow=False )
#        # titlepad=10.0,  # these breaks RTD webhook
#        # autolimit_mode='data')  # round_numbers
# mpl.rc('axes.formatter', limits=(-1,1))
# mpl.rc('xtick', top=True, direction='in')
# mpl.rc('xtick.major', size=9, width=1.4, pad=7)
# mpl.rc('xtick.minor', size=5, width=1.0, pad=7)
# mpl.rc('ytick', right=True, direction='in')
# mpl.rc('ytick.major', size=9, width=1.4, pad=7)
# mpl.rc('ytick.minor', size=5, width=1.0, pad=7)
# mpl.rc('savefig', facecolor='ffffff', dpi=100, bbox='tight')
mpl.rc_file(style)
mpl.rc('axes', prop_cycle=cc)

# mpl.rc('axes', facecolor='ffffff')
# mpl.rc('figure', facecolor='E7E0D6')

# axes Formatter
def customFormatter(factor, prec=1, width=2):
    """create a mpl formatter whith specialized format 
    for plotting labels:
    width(prec) x 10^factor.
    """
    fstr = '{:{width}.{prec}f}'
    exp = 10.0**factor
    return FuncFormatter(lambda x, pos:fstr.format(x/exp, 
                         width=width, prec=prec))


def writeFig(fig, paths, filetag):
    """writes figure to file according to folders in path.
    
    Args:
        fig(mpl.figure): matplotlib object to store.
        paths(str list): output paths.
        filetag(str): preffix for output file.
        
    Returns:
        (str): destination path of the file.
        (str): file suffix number.
    
    """
    num = paths[1][-5:]  # checkpoint number 'flash_hdf5_chk_0001'
    dest = os.path.join(os.path.dirname(paths[0]), filetag)
    name = os.path.join(dest, '{}{}.png'.format(filetag, num))
    os.makedirs(dest, exist_ok=True)  # bless you, p3
    plt.savefig(name, format='png')
    plt.close(fig)
    print("Wrote: {}".format(name))
    return dest, num
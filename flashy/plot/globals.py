from ..IOutils import os
from ..utils import np
import pkg_resources
from cycler import cycler
from matplotlib.ticker import (FuncFormatter,
                               StrMethodFormatter, ScalarFormatter)
from matplotlib.colors import LogNorm, NoNorm
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import AxesGrid
# backends for script/parallel ploting
nonguis = [u'agg', u'cairo', u'gdk', u'pdf',
           u'pgf', u'ps', u'svg', u'template']
for nongui in nonguis:
    try:
        # print "[flash_plot]: Testing", nongui
        mpl.use(nongui, warn=False, force=True)
        from matplotlib import pyplot as plt
        break
    except:
        continue
# print("[flash.plot]: Using",mpl.get_backend())
# x = {1, 5, 10}
# (0, (1, x)) dotted
# (0, (5, x)) dashed
# (0, (3, x, 1, x)) dash-dot
# (0, (1, x, 1, x, 1, x)) dash-dot-dot
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
style = pkg_resources.resource_filename('flashy',
                                        '../data/mplstyles/flashydef.mplsty')
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
    return FuncFormatter(lambda x, pos: fstr.format(x/exp,
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


# hex/rgb handler from Ben Southgate:
# https://bsou.io/posts/color-gradients-with-python
def hex_to_RGB(hex):
    """"#FFFFFF" -> [255,255,255] """
    # Pass 16 to the integer function for change of base
    return [int(hex[i:i+2], 16) for i in range(1, 6, 2)]


def RGB_to_hex(RGB):
    """[255,255,255] -> "#FFFFFF" """
    # Components need to be integers for hex to make sense
    RGB = [int(x) for x in RGB]
    s = ["0{0:x}".format(v) if v < 16 else "{0:x}".format(v) for v in RGB]
    return "#"+"".join(s)


def color_dict(gradient):
    """Takes in a list of RGB sub-lists and returns dictionary of
    colors in RGB and hex form for use in a graphing function
    defined later on
    """
    return {"hex": [RGB_to_hex(RGB) for RGB in gradient],
            "r": [RGB[0] for RGB in gradient],
            "g": [RGB[1] for RGB in gradient],
            "b": [RGB[2] for RGB in gradient]}


def linear_gradient(start_hex, finish_hex="#FFFFFF", n=10):
    """returns a gradient list of (n) colors between
    two hex colors. start_hex and finish_hex
    should be the full six-digit color string,
    inlcuding the number sign ("#FFFFFF")
    """
    # Starting and ending colors in RGB form
    s = hex_to_RGB(start_hex)
    f = hex_to_RGB(finish_hex)
    # Initilize a list of the output colors with the starting color
    RGB_list = [s]
    # Calcuate a color at each evenly spaced value of t from 1 to n
    for t in range(1, n):
        # Interpolate RGB vector for color at the current value of t
        curr_vector = [
            int(s[j] + (float(t)/(n-1))*(f[j]-s[j]))
            for j in range(3)
        ]
        # Add it to our list of output colors
        RGB_list.append(curr_vector)

    return color_dict(RGB_list)

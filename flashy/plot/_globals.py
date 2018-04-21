from flashy.IOutils import os
from flashy.utils import np
from flashy.nuclear import sortNuclides, elemSplit
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
print "[flash.plot]: Using",mpl.get_backend()
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

# plot customizations
def customFormatter(factor, prec=1, width=2):
    """specialized format for plotting labels:
    width(prec) x 10^factor.
    """
    fstr = '{:{width}.{prec}f}'
    exp = 10.0**factor
    return FuncFormatter(lambda x, pos:fstr.format(x/exp, 
                         width=width, prec=prec))


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

# SPECIES just in case...
xnet = [
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
xnetReduced = [
    'neut', 'h1  ', 'h2  ', 'he3 ', 'he4 ', 'li7 ', 'b8  ', 'c12 ', 'n14 ',
    'o16 ', 'f18 ', 'ne20', 'na22', 'na23', 'mg24', 'al26', 'al27', 'si28', 
    'p30 ', 's32 ', 'cl34', 'ar36', 'ar40', 'k38 ', 'ca40', 'ca44', 'sc44',
    'ti44', 'ti45', 'ti46', 'v46 ', 'cr48', 'cr50', 'mn50', 'fe52', 'fe53', 
    'fe54', 'fe55', 'fe56', 'fe57', 'fe58', 'co54', 'co55', 'co56', 
    'co54', 'ni56', 'ni58', 'cu58', 'zn60'
]
ap13 = [
    'he4 ', 'c12 ', 'o16 ', 'ne20',
    'mg24', 'si28', 's32 ', 'ar36',
    'ca40', 'ti44', 'cr48', 'fe52', 'ni56'
]
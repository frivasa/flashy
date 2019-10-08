import qgrid
import ipywidgets as widgets
from .IOutils import os, sys, io, _otpfolder, _cdxfolder, _cdxpfol
from traitlets import traitlets
import flashy.meta as meta
from IPython.display import clear_output
_gopts = {'enableColumnReorder': True, 'editable': False}


def getPreRunBox(path):
    codever = widgets.Dropdown(options=os.listdir(path))
    prepath = '/'.join([path, codever.value, _cdxpfol])
    profopts = sorted([x for x in os.listdir(prepath)])
    profs = widgets.Dropdown(options=profopts)

    def update_profs(*args):
        prepath = '/'.join([path, codever.value, _cdxpfol])
        profopts = sorted([x for x in os.listdir(prepath)])
        profs.options = profopts

    codever.observe(update_profs, 'value')
    box = widgets.HBox([codever, profs])
    return box


def getPickBox(path):
    codever = widgets.Dropdown(options=os.listdir(path))
    opts = sorted([x for x in
                   os.listdir(os.path.join(path, codever.value))
                   if _cdxfolder not in x])
    runs = widgets.Dropdown(options=opts)
#             sorted([x for x in os.listdir(os.path.join(path, codever.value))
#                     if _cdxfolder not in x]))

    def update_runs(*args):
        runs.options = \
            sorted([x for x in os.listdir(os.path.join(path, codever.value))
                    if _cdxfolder not in x])

    codever.observe(update_runs, 'value')
    box = widgets.HBox([codever, runs])
    return box


def getChkPickBox(path):
    codever = widgets.Dropdown(options=os.listdir(path))
    ropts = sorted([x for x in
                    os.listdir(os.path.join(path, codever.value))
                    if 'cdx' not in x])
    runs = widgets.Dropdown(options=ropts)
    chkfs = sorted(os.listdir(os.path.join(path, codever.value,
                                           runs.value, _otpfolder)))
    chkopts = list(zip(chkfs, ['/'.join([_otpfolder, v]) for v in chkfs]))
    chks = widgets.Dropdown(options=chkopts)

    def update_runs(*args):
        ropts = sorted([x for x in
                        os.listdir(os.path.join(path, codever.value))
                        if 'cdx' not in x])
        runs.options = ropts

    def update_chks(*args):
        chkfs = sorted(os.listdir(os.path.join(path, codever.value,
                                               runs.value, _otpfolder)))
        chkopts = list(zip(chkfs, ['/'.join([_otpfolder, v]) for v in chkfs]))
        chks.options = chkopts

    codever.observe(update_runs, 'value')
    runs.observe(update_chks, 'value')
    box = widgets.VBox([codever, runs, chks])
    return box


def gridPandas(DataFrame):
    qwidget = qgrid.show_grid(DataFrame, show_toolbar=True,
                              grid_options=_gopts)
    return qwidget


class fatButton(widgets.Button):
    def __init__(self, dframe=None, *args, **kwargs):
        super(fatButton, self).__init__(*args, **kwargs)
        self.add_traits(dframe=traitlets.Any(dframe))


def metaNavigator(path, **kwargs):
    """layout is [[dropdown, dropdown], [button,button], [terminal]]"""
    box = getPickBox(path)
    btnR = widgets.Button(description='Check Run')
    btnF = fatButton(description='Check Folder', dframe=None)
    buttons = widgets.HBox([btnR, btnF])
    miniterm = widgets.Output()
    disp = widgets.VBox([box, buttons, miniterm])

    def getNameFromBox(path):
        # read the picked run
        for ch in box.children:
            path = os.path.join(path, ch.value)
        return path

    def runHandler(btn):
        # meta for a single picked run
        name = getNameFromBox(path)
        stdout = sys.stdout
        sys.stdout = io.StringIO()
        print('\n\t', name)
        tags, vals = meta.getRunMeta(name)
        # get output and restore sys.stdout
        output = sys.stdout.getvalue()
        sys.stdout = stdout
        with miniterm:
            clear_output()
            print(output)
            for t, v in zip(tags, vals):
                print('{:20} {:20}'.format(t, v))

    def folderHandler(fbtn):
        # meta for a whole folder
        name = getNameFromBox(path)
        folder = os.path.dirname(name)
        stdout = sys.stdout
        sys.stdout = io.StringIO()
        fbtn.dframe = meta.checkCodeFolder(folder, **kwargs)
        output = sys.stdout.getvalue()
        sys.stdout = stdout
        with miniterm:
            clear_output()
            print(output)
    btnR.on_click(runHandler)
    btnF.on_click(folderHandler)
    return disp

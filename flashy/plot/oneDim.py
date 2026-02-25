from .globals import (np, os, plt, setFolders,
                      StrMethodFormatter, writeFig)
from ..datahaul.hdf5yt import getLineout, probeDomain
from ..datahaul.hdfdirect import directMeta
import flashy.utils as ut
import flashy.paraMan as pman
from ..datahaul.plainText import DataMatrix
from ..nuclear import sort_nuclides, elemSplit, decayYield, getMus
from ..simulation import simulation
import flashy.post as fp
from scipy.integrate import trapz
from flashy.datahaul.parData import glue2dWedges
from ..datahaul.ytfields import _alphas


def plot2DtauRatio(bview, fname, factor=0.1, wedges=5, batch=False):
    """Plots sound crossing timescale to burning
    timescale for a file.

    Args:
        bview(ipp.LoadBalancedView): ipp setup workhorse.
        fname(str): file name.
        wedges(int): parallel extraction slices.
        batch(bool): write to file and skip returning figure.

    Returns:
        (mpl.figure or None)

    """
    time, pars, _, _, paths = directMeta(fname)
    dx, tauC, tauE = fp.get2Dtaus(bview, fname, wedges=wedges)
    f, ax = plt.subplots()
    ax.scatter(dx/1e5, tauC/(factor*tauE), marker='.',
               label=u'$\\tau_s/\\tau_e$')
    ax.set_xscale('log')
    tag = os.path.basename(os.path.dirname(os.path.dirname(fname)))
    ax.set_title('{} ({:.3f} s)'.format(tag, time))
    ax.set_xlabel('Scale (km)')
    ax.set_ylabel(u'$\\tau_s/\\tau_e$')
    lgd = ax.legend()
    if batch:
        filetag = 'timescales'
        writeFig(f, fname, filetag)
    else:
        return f


def plot2Dtaus(bview, fname, wedges=5, batch=False):
    """Plots sound crossing timescale to burning
    timescale for a file.

    Args:
        bview(ipp.LoadBalancedView): ipp setup workhorse.
        fname(str): file name.
        wedges(int): parallel extraction slices.
        batch(bool): write to file and skip returning figure.

    Returns:
        (mpl.figure or None)

    """
    time, pars, _, _, paths = directMeta(fname)
    dx, tauC, tauE = fp.get2Dtaus(bview, fname, wedges=wedges)
    f, ax = plt.subplots()
    ax.scatter(dx, tauC, marker='.', label=u'$\\tau_s$')
    ax.scatter(dx, tauE, marker='.', label=u'$\\tau_e$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    tag = os.path.basename(os.path.dirname(os.path.dirname(fname)))
    ax.set_title('{} ({:.3f} s)'.format(tag, time))
    ax.set_xlabel('Scale (cm)')
    ax.set_ylabel('Tau (s)')
    lgd = ax.legend()
    if batch:
        filetag = 'timescales'
        writeFig(f, fname, filetag)
    else:
        return f


def plotRadialSpeeds(bview, fname, slices=5, dimension=2,
                     ref='x', antipode=False, avoid=0.0,
                     batch=False):
    """Plots Radial speed vs radius for a file.

    Args:
        bview(ipp.LoadBalancedView): ipp setup workhorse.
        fname(str): file name.
        slices(int): parallel extraction slices.
        dimension(int): specify file dimension.
        ref(str): reference axis (3D only, see datahaul.hdf5yt.wedge3d).
        antipode(bool): antipodal wedge (see datahaul.hdf5yt.wedge3d).
        avoid(float): polar angle to exclude (from both poles).
        batch(bool): write to file and skip returning figure.

    Returns:
        (mpl.figure or None)

    """
    kwargs = {'fname': os.path.abspath(fname), 'dimension': dimension,
              'wedges': slices, 'ref': ref, 'avoid': avoid, 
              'antipode': antipode}
    res = pman.throwHammer(bview, slices, fp.par_radialSpeeds, **kwargs)
    dat = glue2dWedges(res.get())
    f, ax = plt.subplots()
    ax.scatter(*dat, marker='.', s=0.05)
    line = np.logspace(7, 11, num=30)
    ax.plot(line, line, ls='--', color='black', alpha=0.6)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([2e7, 2e10])
    t, dt, _, _ = probeDomain(fname)
    tag = os.path.basename(os.path.dirname(os.path.dirname(fname)))
    if avoid:
        ttl = '{} ({:.3f} s, -{:d}°)'.format(tag, t, avoid)
    else:
        ttl = '{} ({:.3f} s)'.format(tag, t)
    ax.set_title(ttl)
    ax.set_xlabel('Radius (cm)')
    ax.set_ylabel('Speed (cm/s)')
    if batch:
        filetag = 'radialSp'
        writeFig(f, fname, filetag)
    else:
        return f


def adhoc_pLoSSpeedHisto2d(fname, radius, sphereCapR=1e9,
                           geom='cartesian', grouping=False,
                           groups={'GrpLabel': ['he4', 'ni56']},
                           resolution=4e7, ref='y', antipode=False,
                           ylims=[1e-4, 0.5], xlims=[0.1, 4e9], maxR=1.0e9,
                           species=_alphas, massfrac=True, cylvel=False,
                           velrange=[1e7, 3e9],
                           showtotal=False, batch=False, filetag='sphisto'):
    """adhoc involves grouping species according to a dict.
    dict tags are plot labels and each is the list of species to sum.

    """
    databuckets, bins = fp.speedHisto2d(fname, radius, geom=geom,
                                        sphereCapR=sphereCapR,
                                        velrange=velrange, ref=ref,
                                        resolution=resolution, species=species)
    print('speedHisto done')
    # only 1 wedge to distribute
    histobuckets = []
    for s in species:
        histobuckets.append(np.zeros(len(bins)))
    tmass = np.zeros(len(bins))
    # wedge data already in mass units
    for s in range(len(species)):
        histobuckets[s] += databuckets[s]
        tmass += databuckets[s]
    tmass[tmass == 0.0] = 1.0
    print('histo buckets summed, starting plot')
    t, dt, _, _ = probeDomain(fname)
    # build the plot
    f, ax = plt.subplots()
    lss = ['-', '--']
    if massfrac:
        if grouping:
            for i, (tag, slist) in enumerate(groups.items()):
                bbucket = np.zeros(len(bins))
                for s in slist:
                    ll = species.index(s)
                    bbucket += histobuckets[ll]
                weights = np.nan_to_num(bbucket)/tmass
                mpln, mplbs, pts = ax.hist(bins, bins=len(bins), log=True,
                                           weights=weights, # ls=lss[i%len(lss)],
                                           histtype='step', label=tag)
        else:
            for i, s in enumerate(species):
                weights = np.nan_to_num(histobuckets[i])/tmass
                mpln, mplbs, pts = ax.hist(bins, bins=len(bins), log=True,
                                           weights=weights, # ls=lss[i%len(lss)],
                                           histtype='step', label=s)
        ax.set_ylabel(u'$X_i$')
        ax.set_ylim(ylims)
        ax.axhline(1, ls='--', alpha=0.6, c='k')
    else:
        print('massfrac false')
        # for i, (tag, slist) in enumerate(groups.items()):
        #     bbucket = np.zeros(len(bins))
        #     for s in slist:
        #         ll = species.index(s)
        #         bbucket += histobuckets[ll]
        #     mpln, mplbs, pts = ax.hist(bins, bins=len(bins), log=True,
        #                                weights=bbucket/ut.msol,
        #                                histtype='step', label=s)
        for i, s in enumerate(species):
            weight = np.nan_to_num(histobuckets[i])
            mpln, mplbs, pts = ax.hist(bins, bins=len(bins), log=True,
                                       weights=weight/ut.msol,
                                       histtype='step', label=s)
        if showtotal:
            mpln, mplbins, patches = ax.hist(bins,
                                             bins=len(bins), histtype='step',
                                             weights=tmass/ut.msol, log=True,
                                             label='Total', color='#000000')
        ax.set_ylabel(u'$M_{\\odot}$')
        ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    if cylvel:
        ax.set_xlabel('Axial Velocity (km/s)')    
    else:
        ax.set_xlabel('Velocity (km/s)')
    tag = os.path.basename(os.path.dirname(os.path.dirname(fname)))
    ax.set_title('{} ({:.3f} s)\n(radius:{:.0f}km,cap:{:.0f}km)'.format(tag, t, radius/1e5, sphereCapR/1e5))
    lgd = ax.legend(ncol=1, loc='upper left', bbox_to_anchor=(1.01, 1.0),
                    markerfirst=False, numpoints=3, frameon=False)
    if batch:
        writeFig(f, fname, filetag)
    else:
        return f



def adhoc_pLoSSpeedHisto3d(fname, radius, sphereCapR=1e9,
                           geom='cartesian', dimension=3, grouping=False,
                           groups={'GrpLabel': ['he4', 'ni56']},
                           resolution=4e7, ref='z', antipode=False,
                           ylims=[1e-4, 0.5], xlims=[0.1, 4e9], maxR=1.0e9,
                           species=_alphas, massfrac=True, cylvel=False,
                           velrange=[1e7, 3e9],
                           showtotal=False, batch=False, filetag='sphisto'):
    """adhoc involves grouping species according to a dict.
    dict tags are plot labels and each is the list of species to sum.

    Args:
        elevation(float): equator-north pole degree.
        depth(float): equator-south pole degree.
        *others analogous to plotSpeedHisto*
        
    Returns:
        (mpl.figure or None)

    """
    databuckets, bins = fp.speedHisto3d(fname, radius, geom=geom,
                                        sphereCapR=sphereCapR,
                                        dimension=dimension,
                                        velrange=velrange, ref=ref,
                                        resolution=resolution, species=species)
    print('speedHisto done')
    # only 1 wedge to distribute
    histobuckets = []
    for s in species:
        histobuckets.append(np.zeros(len(bins)))
    tmass = np.zeros(len(bins))
    # wedge data already in mass units
    for s in range(len(species)):
        histobuckets[s] += databuckets[s]
        tmass += databuckets[s]
    tmass[tmass == 0.0] = 1.0
    print('histo buckets summed, starting plot')
    t, dt, _, _ = probeDomain(fname)
    # build the plot
    f, ax = plt.subplots()
    lss = ['-', '--']
    if massfrac:
        if grouping:
            for i, (tag, slist) in enumerate(groups.items()):
                bbucket = np.zeros(len(bins))
                for s in slist:
                    ll = species.index(s)
                    bbucket += histobuckets[ll]
                weights = np.nan_to_num(bbucket)/tmass
                mpln, mplbs, pts = ax.hist(bins, bins=len(bins), log=True,
                                           weights=weights, # ls=lss[i%len(lss)],
                                           histtype='step', label=tag)
        else:
            for i, s in enumerate(species):
                weights = np.nan_to_num(histobuckets[i])/tmass
                mpln, mplbs, pts = ax.hist(bins, bins=len(bins), log=True,
                                           weights=weights, # ls=lss[i%len(lss)],
                                           histtype='step', label=s)
        ax.set_ylabel(u'$X_i$')
        ax.set_ylim(ylims)
        ax.axhline(1, ls='--', alpha=0.6, c='k')
    else:
        print('massfrac false')
        if grouping:
            print('grouping species')
            for i, (tag, slist) in enumerate(groups.items()):
                bbucket = np.zeros(len(bins))
                for s in slist:
                    ll = species.index(s)
                    bbucket += histobuckets[ll]
                mpln, mplbs, pts = ax.hist(bins, bins=len(bins), log=True,
                                           weights=bbucket/ut.msol,
                                           histtype='step', label=s)
        else:
            for i, s in enumerate(species):
                weight = np.nan_to_num(histobuckets[i])
                mpln, mplbs, pts = ax.hist(bins, bins=len(bins), log=True,
                                           weights=weight/ut.msol,
                                           histtype='step', label=s)
        if showtotal:
            mpln, mplbins, patches = ax.hist(bins,
                                             bins=len(bins), histtype='step',
                                             weights=tmass/ut.msol, log=True,
                                             label='Total', color='#000000')
        ax.set_ylabel(u'$M_{\\odot}$')
        ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    if cylvel:
        ax.set_xlabel('Axial Velocity (km/s)')    
    else:
        ax.set_xlabel('Velocity (km/s)')
    tag = os.path.basename(os.path.dirname(os.path.dirname(fname)))
    ax.set_title('{} ({:.3f}s)\n(radius:{:.0f}km,cap:{:.0f}km)'.format(tag, t, radius/1e5, sphereCapR/1e5))
    lgd = ax.legend(ncol=1, loc='upper left', bbox_to_anchor=(1.01, 1.0),
                    markerfirst=False, numpoints=3, frameon=False)
    if batch:
        writeFig(f, fname, filetag)
    else:
        return f


def adhoc_pLoSSpeedHisto(fname, elevation=89.9, depth=-40, 
                         geom='cartesian', dimension=2,
                         groups={'GrpLabel': ['he4', 'ni56']},
                         resolution=4e7, ref='x', antipode=False,
                         ylims=[1e-4, 0.5], xlims=[0.1, 4e9], maxR=1.0e9,
                         species=_alphas, massfrac=False, cylvel=False,
                         velrange=[1e7, 3e9],
                         showtotal=False, batch=False, filetag='sphisto'):
    """plots speed vs mass fraction for an angled wedge in the domain.
    adhoc involves grouping species according to a dict.
    dict tags are plot labels and each is the list of species to sum.

    Args:
        elevation(float): equator-north pole degree.
        depth(float): equator-south pole degree.
        *others analogous to plotSpeedHisto*
        
    Returns:
        (mpl.figure or None)

    """
    databuckets, bins = fp.speedHisto(fname, geom=geom, dimension=dimension,
                                      velrange=velrange,
                                      ref=ref, resolution=resolution, maxR=maxR,
                                      antipode=antipode, depth=depth, cylvel=cylvel,
                                      species=species, elevation=elevation)
    print('speedHisto done')
    # only 1 wedge to distribute
    histobuckets = []
    for s in species:
        histobuckets.append(np.zeros(len(bins)))
    tmass = np.zeros(len(bins))
    # wedge data already in mass units
    for s in range(len(species)):
        histobuckets[s] += databuckets[s]
        tmass += databuckets[s]
    tmass[tmass == 0.0] = 1.0
    print('histo buckets summed, starting plot')
    t, dt, _, _ = probeDomain(fname)
    # build the plot
    f, ax = plt.subplots()
    lss = ['-', '--']
    if massfrac:
        for i, (tag, slist) in enumerate(groups.items()):
            bbucket = np.zeros(len(bins))
            for s in slist:
                ll = species.index(s)
                bbucket += histobuckets[ll]
            weights = np.nan_to_num(bbucket)/tmass
            mpln, mplbs, pts = ax.hist(bins, bins=len(bins), log=True,
                                       weights=weights, # ls=lss[i%len(lss)],
                                       histtype='step', label=tag) 
#         for i, s in enumerate(species):
#             weights = np.nan_to_num(histobuckets[i])/tmass
#             mpln, mplbs, pts = ax.hist(bins, bins=len(bins), log=True,
#                                        weights=weights, # ls=lss[i%len(lss)],
#                                        histtype='step', label=s)
        ax.set_ylabel(u'$X_i$')
        ax.set_ylim(ylims)
        ax.axhline(1, ls='--', alpha=0.6, c='k')
    else:
        for i, (tag, slist) in enumerate(groups.items()):
            bbucket = np.zeros(len(bins))
            for s in slist:
                ll = species.index(s)
                bbucket += histobuckets[ll]
            mpln, mplbs, pts = ax.hist(bins, bins=len(bins), log=True,
                                       weights=bbucket/ut.msol,
                                       histtype='step', label=s)
        if showtotal:
            mpln, mplbins, patches = ax.hist(bins,
                                             bins=len(bins), histtype='step',
                                             weights=tmass/ut.msol, log=True,
                                             label='Total', color='#000000')
        ax.set_ylabel(u'$M_{\\odot}$')
        ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    if cylvel:
        ax.set_xlabel('Axial Velocity (cm/s)')    
    else:
        ax.set_xlabel('Velocity (cm/s)')
    tag = os.path.basename(os.path.dirname(os.path.dirname(fname)))
    ax.set_title('{} ({:.3f} s) ({},{})'.format(tag, t, elevation, depth))
    lgd = ax.legend(ncol=1, loc='upper left', bbox_to_anchor=(1.01, 1.0),
                    markerfirst=False, numpoints=3, frameon=False)
    if batch:
        writeFig(f, fname, filetag)
    else:
        return f


def pLoSSpeedHisto(fname, elevation=89.9, depth=-40, 
                   geom='cartesian', dimension=2,
                   resolution=4e7, ref='x', antipode=False,
                   ylims=[1e-4, 0.5], xlims=[0.1, 4e9], kms=False,
                   species=_alphas, massfrac=False, cylvel=False,
                   showtotal=False, batch=False, filetag='sphisto'):
    """plots speed vs mass fraction for an angled wedge in the domain.

    Args:
        elevation(float): equator-north pole degree.
        depth(float): equator-south pole degree.
        kms(bool): factor bins by 1e5 to plot km/s.
        *others analogous to plotSpeedHisto*
        
    Returns:
        (mpl.figure or None)

    """
    databuckets, bins = fp.speedHisto(fname, geom=geom, dimension=dimension,
                                      ref=ref, resolution=resolution,
                                      antipode=antipode, depth=depth, cylvel=cylvel,
                                      species=species, elevation=elevation)
    # only 1 wedge to distribute
    histobuckets = []
    for s in species:
        histobuckets.append(np.zeros(len(bins)))
    if kms:
        bins = bins/1e5
    tmass = np.zeros(len(bins))
    # wedge data already in mass units
    for s in range(len(species)):
        histobuckets[s] += databuckets[s]
        tmass += databuckets[s]
    tmass[tmass == 0.0] = 1.0

    t, dt, _, _ = probeDomain(fname)
    # build the plot
    f, ax = plt.subplots()
    lss = ['-', '--']
    if massfrac:
        for i, s in enumerate(species):
            weights = np.nan_to_num(histobuckets[i])/tmass
            mpln, mplbs, pts = ax.hist(bins, bins=len(bins), log=True,
                                       weights=weights, # ls=lss[i%len(lss)],
                                       histtype='step', label=s)
        ax.set_ylabel(u'$X_i$')
        ax.set_ylim(ylims)
        ax.axhline(1, ls='--', alpha=0.6, c='k')
    else:
        for i, s in enumerate(species):
            mpln, mplbs, pts = ax.hist(bins, bins=len(bins), log=True,
                                       weights=histobuckets[i]/ut.msol,
                                       histtype='step', label=s)
        if showtotal:
            mpln, mplbins, patches = ax.hist(bins,
                                             bins=len(bins), histtype='step',
                                             weights=tmass/ut.msol, log=True,
                                             label='Sum', color='#9F9F9F')
            print('Mass sum from all bins: {:.4E}'.format(np.sum(tmass/ut.msol)))
        ax.set_ylabel(u'$M_{\\odot}$')
        ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    unit = 'km/s' if kms else 'cm/s'
    if cylvel:
        ax.set_xlabel('Axial Velocity ({})'.format(unit))    
    else:
        ax.set_xlabel('Velocity ({})'.format(unit))
    tag = os.path.basename(os.path.dirname(os.path.dirname(fname)))
    ax.set_title('{} ({:.3f} s) ({},{})'.format(tag, t, elevation, depth))
    lgd = ax.legend(ncol=1, loc='upper left', bbox_to_anchor=(1.01, 1.0),
                    markerfirst=False, numpoints=3, frameon=False)
    if batch:
        writeFig(f, fname, filetag)
    else:
        return f


def plotSpeedHisto(bview, fname, geom='cartesian', dimension=2,
                   resolution=4e7, ref='x', antipode=False, slices=5,
                   ylims=[1e-4, 0.5], xlims=[0.1, 4e9],
                   species=_alphas, massfrac=False, cylvel=False,
                   showtotal=False, batch=False, filetag='sphisto'):
    """Plots figure with speed vs mass fraction histogram for a file.

    Args:
        bview(ipp.LoadBalancedView): ipp setup workhorse.
        fname(str): file name.
        geom(str): specify geometry for 1d file.
        dimension(int): specify file dimension.
        resolution(float): bin size in cm/s
        ref(str): reference axis (3D only, see datahaul.hdf5yt.wedge3d).
        antipode(bool): antipodal wedge (see datahaul.hdf5yt.wedge3d).
        slices(int): parallel extraction slices.
        ylims(float list): ordinate range.
        xlims(float list): abscissa range.
        species(str list): list of species to plot.
        massfrac(bool): draw X_i instead of mass.
        cylvel(bool): get velx as speed (cylindrical radial speed).
        showtotal(bool): draw summed mass histogram.
        batch(bool): write to file and skip returning figure.
        filetag(str): batch file nametag.

    Returns:
        (mpl.figure or None)

    """
    kwargs = {'fname': os.path.abspath(fname), 'dimension': dimension,
              'resolution': resolution, 'wedges': slices, 'cylvel': cylvel,
              'species': species, 'ref': ref, 'antipode': antipode}
    res = pman.throwHammer(bview, slices, fp.par_speedHisto, **kwargs)
    # retrieve the results and sum the wedges
    cont = res.get()
    databuckets, bins = cont[0]
    histobuckets = []
    for s in species:
        histobuckets.append(np.zeros(len(bins)))
    tmass = np.zeros(len(bins))
    # wedge data already in mass units
    for wedge in cont[0]:
        for s in range(len(species)):
            histobuckets[s] += wedge[s]
            tmass += wedge[s]
    tmass[tmass == 0.0] = 1.0
    t, dt, _, _ = probeDomain(fname)
    # build the plot
    f, ax = plt.subplots()
    lss = ['-', '--']
    if massfrac:
        for i, s in enumerate(species):
            weights = np.nan_to_num(histobuckets[i])/tmass
            mpln, mplbs, pts = ax.hist(bins, bins=len(bins), log=True,
                                       weights=weights, ls=lss[i%len(lss)],
                                       histtype='step', label=s)
        ax.set_ylabel(u'$X_i$')
        ax.set_ylim(ylims)
        ax.axhline(1, ls='--', alpha=0.6, c='k')
    else:
        for i, s in enumerate(species):
            mpln, mplbs, pts = ax.hist(bins, bins=len(bins), log=True,
                                       weights=histobuckets[i]/ut.msol,
                                       histtype='step', label=s)
        if showtotal:
            mpln, mplbins, patches = ax.hist(bins,
                                             bins=len(bins), histtype='step',
                                             weights=tmass/ut.msol, log=True,
                                             label='Total', color='#000000')
        ax.set_ylabel(u'$M_{\\odot}$')
        ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    if cylvel:
        ax.set_xlabel('Axial Velocity (cm/s)')    
    else:
        ax.set_xlabel('Velocity (cm/s)')
    tag = os.path.basename(os.path.dirname(os.path.dirname(fname)))
    ax.set_title('{} ({:.3f} s)'.format(tag, t))
    lgd = ax.legend()
    if batch:
        writeFig(f, fname, filetag)
    else:
        return f


def fetchData(fname, direction, origin=[0.0, 0.0, 0.0], fields=[]):
    """builds a profile DataMatrix and retrieves metadata from a file.

    Args:
        fname(str): filename.
        direction(float list): empty for 1d.
            angle from equator for 2d.
            spherical angles for 3d (polar angle and azimuth).
        origin(float list): override origin of ray to build profile.
        fields(str list): specify field to extract (optional).

    Returns:
        (float): timestamp for the file.
        (dict): parameter dictionary for the checkpoint.
        (str list): filesystem paths for the file.
        (DataMatrix): block of data for plotting.

    """
    time, pars, _, _, paths = directMeta(fname)
    if len(direction) > (pars['dimensionality']-1):
        print("WARN: Direction doesn't match "
              "dimensionality: {}".format(pars['dimensionality']))
    keys = ['radius']
    if fields:
        ad, allsp = getLineout(fname, geom=pars['geometry'],
                               fields=fields, origin=origin,
                               direction=direction, srcnames=False)
        keys += fields
    else:
        ad, allsp = getLineout(fname, geom=pars['geometry'],
                               origin=origin,
                               direction=direction, srcnames=False)
        keys += ['density', 'temperature', 'pressure']
    keys += allsp
    # add time as a column so it can be read elsewhere easily
    ad = np.vstack([ad, [time]*ad.shape[1]])
    keys += ['simtime']
    return time, pars, paths, DataMatrix([keys, ad.transpose()])


def shockFollow(fname, simfolder='', thresh=1e-4, batch=False, byM=False,
                direction=[], wakesize=1e8, inward=False, grid=False):
    """plot shock wake for a checkpoint.
    WARN: slowest pos will not match for non x-axis directions.

    Args:
        fname(str): path to checkpoint file.
        simfolder(str): simulation run metadata folder.
        thresh(float): species threshold.
        batch(bool): write file instead of returning figure.
        byM(bool): plot by mass coordinate instead of radius.
        direction(str list): polar angle and azimuth of profile ray.
        wakesize(float): extent of plot in fraction of shock position.
        inward(bool): shock to follow.
        grid(bool): add radial positions as a line grid.

    Return:
        (mpl figure) or (None)

    """
    fig, meta = PARsimProfile(fname, simfolder=simfolder, thresh=thresh,
                              filetag='', batch=False, byM=byM,
                              direction=direction, meta=True, grid=grid)
    # readjust the plot to center the shock
    ax = plt.gca()
    outsh, inwsh, slowx, paths = meta
    if inward:
        left, right = inwsh*0.95, inwsh*(1.0 + wakesize)
        filetag = 'inward_shock'
    else:
        left, right = outsh*(1.0 - wakesize), outsh*1.05
        filetag = 'outward_shock'
    ax.set_xlim([left, right])
    # I/O
    if not batch:
        return fig
    else:
        writeFig(fig, paths[1], filetag)


def PARsimProfile(fname, simfolder='', thresh=1e-4, xrange=[0.0, 0.0],
                  filetag='metaprof', batch=False, byM=False, direction=[],
                  meta=False, grid=False):
    """bloated method for IPP.
    WARN: slowest pos will not match for non x-axis directions.

    Args:
        fname(str): path to checkpoint file.
        simfolder(str): simulation run metadata folder.
        thresh(float): species threshold.
        xrange(float list): abscissa range for all plots.
        filetag(str): prefix for output files (if any).
        batch(bool): write file instead of returning figure.
        byM(bool): plot by mass coordinate instead of radius.
        direction(str list): polar angle and azimuth of profile ray.
        grid(bool): add radial positions as a line grid.

    Return:
        (mpl figure) or (None)

    """

    time, pars, paths, prof = fetchData(fname, direction)
    sim = simulation(simfolder)
    n, t = sim.time2step(time)
    step = sim.steps[int(n)]
    xm, ym, zm = step.slowx, step.slowy, step.slowz
    rm = np.sqrt(xm**2+ym**2+zm**2)
    fig = plotDMatMerged(prof, byM=byM,
                         thresh=thresh, xrange=xrange)
    # add extra marks from simulation obj
    ax = plt.gca()
    a = ax.annotate("{:.5f} s".format(time),
                    xy=(0.0, 0.0), xytext=(0.84, 0.10), size=12,
                    textcoords='figure fraction', xycoords='figure fraction',
                    bbox=dict(boxstyle='round', fc='w', ec='k'))
    axlist = fig.axes
    for ax in axlist[1:]:
        p = rm if not byM else prof.masses[np.where(prof.radius > rm)[0][0]]
        ax.axvline(p, ls=':', color='brown', alpha=0.5, label='slowest')
    if sim.CJ:
        times, xins, cjins, xouts, cjouts = np.genfromtxt(sim.CJ[0],
                                                          unpack=True)
        if byM:
            xo = xouts[pars['checkpointfilenumber']]
            xi = xins[pars['checkpointfilenumber']]
            no = np.where(prof.radius > xo)[0][0]
            ni = np.where(prof.radius > xi)[0][0]
            po, pi = prof.masses[no], prof.masses[ni]
        else:
            po = xouts[pars['checkpointfilenumber']]
            pi = xins[pars['checkpointfilenumber']]
        ax.axvline(po, ls=':', color='green', alpha=0.5, label='shock')
        ax.axvline(pi, ls=':', color='red', alpha=0.5, label='shock')
    lgd = ax.legend(ncol=1, loc='upper left', bbox_to_anchor=(1.00, 0.50),
                    columnspacing=0.0, labelspacing=0.0, markerfirst=False,
                    numpoints=3, handletextpad=0.0, edgecolor='k')
    if grid:
        drawGrid(ax, prof.radius)
    if meta:
        markings = [po, pi, p, paths]
        return fig, markings
    # I/O
    if not batch:
        return fig
    else:
        writeFig(fig, paths[1], filetag)


def flashProfile(fname, thresh=1e-6, xrange=[0.0, 0.0],
                 yrange=[0.0, 0.0], filetag='prof', batch=False,
                 byM=True, direction=[], grid=False, points=False,
                 origin=[0.0, 0.0, 0.0], extgrid='',
                 fields=['pressure', 'temperature']):
    """Plot bulk properties and species in a chekpoint through a ray.

    Args:
        fname(str): path of file.
        thresh(float): threshold for species fraction.
        xrange(float list): abscissa range.
        yrange(float list): bottom axes yrange.
        filetag(str): prefix for batch mode.
        batch(bool): skips returning figure,
        saving it to a structured directory instead.
        byM(bool): plot by mass instead of radius.
        direction(float list): list of spherical angles (alt, azimuth),
        empty for 1D.
        grid(bool): add radial positions as a line grid.
        points(bool): draw y-points in graphs.
        origin(float list): override origin of ray to build profile.
        extgrid(str): filename for external grid points.
        fields(str list): specify fields for lower plot.

    Returns:
        (mpl figure) or (None)

    """
    time, pars, paths, prof = fetchData(fname, direction, 
                                        origin=origin,
                                        fields=['density']+fields)
    fig = plotDMatMerged(prof, byM=byM, thresh=thresh,
                         xrange=xrange, marker=points)
    ax = plt.gca()
    a = ax.annotate("{:.5f} s".format(time),
                    xy=(0.0, 0.0), xytext=(0.82, 0.10), size=12,
                    textcoords='figure fraction', xycoords='figure fraction',
                    bbox=dict(boxstyle='round', fc='w', ec='k'))
    if direction:
        direc = [str(i) for i in direction]
        b = ax.annotate("Angle(s): {}".format(','.join(direc)),
                        xy=(0.0, 0.0), xytext=(0.82, 0.05), size=12,
                        textcoords='figure fraction',
                        xycoords='figure fraction',
                        bbox=dict(boxstyle='round', fc='w', ec='k'))
    if sum(yrange) != 0.0:
        ax.set_ylim(yrange)
    if grid:
        for ax in fig.axes:
            drawGrid(ax, prof.radius)
    if extgrid:
        exgrdm = DataMatrix(extgrid)
        for ax in fig.axes:
            drawGrid(ax, exgrdm.radius)
    if not batch:
        return fig
    else:
        dest, num, name = setFolders(paths[1], filetag)
        prof.writeProf(name + '.prof')
        writeFig(fig, paths[1], filetag)


def flashDegeneracy(fname, thresh=1e-6, filetag='deg', batch=False,
                    byM=True, direction=[],  xrange=[0.0, 0.0]):
    """compare fermi temperature to found temperature in a lineout."""
    time, pars, paths, prof = fetchData(fname, direction)
    fig = plotDegen(prof, byM=byM, thresh=thresh, xrange=xrange)
    ax = plt.gca()
    a = ax.annotate("{:.5f} s".format(time),
                    xy=(0.0, 0.0), xytext=(0.88, 0.10), size=12,
                    textcoords='figure fraction', xycoords='figure fraction',
                    bbox=dict(boxstyle='round', fc='w', ec='k'))
    if direction:
        direc = [str(i) for i in direction]
        b = ax.annotate("Angle(s): {}".format(','.join(direc)),
                        xy=(0.0, 0.0), xytext=(0.82, 0.05),
                        size=12, textcoords='figure fraction',
                        xycoords='figure fraction',
                        bbox=dict(boxstyle='round', fc='w', ec='k'))
    if not batch:
        return fig
    else:
        writeFig(fig, paths[1], filetag)


def flashSpecies(fname, thresh=1e-6, filetag='spec', batch=False,
                 xrange=[0.0, 0.0], byM=True, direction=[],
                 vmax=4e9, plotall=False):
    """Plot species and aggregated masses in a chekpoint through a ray.
    Masses are a trapezoid integrated ray for each species.

    Args:
        fname(str): path of file.
        thresh(float): threshold for species fraction and mass yields.
        filetag(str): prefix for batch mode.
        batch(bool): skips returning figure,
        saving it to a structured directory instead.
        byM(bool): plot by mass instead of radius.
        direction(float list): list of spherical angles (alt, azimuth),
        empty for 1D.
        vmax(floar): limit for species velocities.
        plotall(bool): force plotting every species found.

    Returns:
        (mpl figure) or (None)

    """
    fields = ['density', 'temperature', 'pressure', 'speed']
    time, pars, paths, prof = fetchData(fname, direction, fields=fields)
    fig = plt.figure(figsize=(10, 8))
    layout = (3, 3)
    # plot species
    ax1 = plt.subplot2grid(layout, (0, 0), colspan=2)
    skip = plotSpecies(ax1, prof, byMass=byM, thresh=thresh, plotall=plotall)
    if byM:
        ax1.set_xlabel('Mass ($M_{\odot}$)')
    else:
        ax1.set_xlabel('Radius (cm)')
    if sum(xrange) != 0.0:
        ax1.set_xlim(xrange)
    # timestamp and legend
    a = ax1.annotate("{:.5f} s".format(time),
                     xy=(0.0, 0.0), xytext=(0.65, 0.1), size=12,
                     textcoords='figure fraction', xycoords='figure fraction',
                     bbox=dict(boxstyle='round', fc='w', ec='k'))
    lgd = ax1.legend(ncol=4, loc='upper left', bbox_to_anchor=(1.05, 1.0),
                     columnspacing=0.0, labelspacing=0.0, markerfirst=False,
                     numpoints=3, handletextpad=0.0)
    lgd.get_frame().set_edgecolor('k')

    # write otp files for yields
    massfiletuples = []

    # plot masses
    ax2 = plt.subplot2grid(layout, (1, 0), colspan=2)
    for i, sp in enumerate(prof.species):
        props = next(ax2._get_lines.prop_cycler)
        spmass = trapz(getattr(prof, sp), x=prof.masses)
        massfiletuples.append((sp, spmass))
        if i in skip:
            continue
        a = elemSplit(sp)[1]
        spmass = trapz(getattr(prof, sp), x=prof.masses)
        ax2.scatter(a, spmass, label=sp, color=props['color'])
    ax2.set_ylabel('Total Mass ($M_\odot$)', rotation=90, labelpad=5)
    ax2.set_xlabel('Mass Number (A)')
    ax2.set_ylim([thresh, 1.5])
    ax2.set_yscale("log", nonposy='clip')
    ax2.yaxis.set_minor_formatter(StrMethodFormatter(''))
    ax2.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
    ax2.xaxis.set_minor_formatter(StrMethodFormatter(''))

    # decay yield
    names, decmasses = decayYield(*zip(*massfiletuples))
    # plot mass vs speed
    ax3 = plt.subplot2grid(layout, (2, 0), colspan=2)
    nbins = 60
    fac = 1e5
    for i, sp in enumerate(prof.species):
        props = next(ax2._get_lines.prop_cycler)
        if i in skip:
            continue
        counts, bins = np.histogram(abs(prof.speed), bins=nbins)
        joined = sorted(zip(abs(prof.speed), getattr(prof, sp)))
        speeds, masses = zip(*joined)
        start = 0
        weights = []
        for c in counts:
            massf = sum(masses[start:start+c])
            weights.append(massf/c)
            start += c
        weights = np.array(weights)
        bins = bins[1:]
        mpln, mplbins, patches = ax3.hist(bins/fac, bins=len(bins),
                                          weights=weights, histtype='step',
                                          log=True, label=sp)
    ax3.set_ylim([thresh, 2])
    ax3.set_xlim([0.0, vmax/fac])
    ax3.yaxis.set_minor_formatter(StrMethodFormatter(''))
    ax3.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
    ax3.xaxis.set_minor_formatter(StrMethodFormatter(''))
    ax3.set_xlabel('Speed (km/s)')
    fig.subplots_adjust(hspace=0.4)

    if not batch:
        return fig
    else:
        dest, num = writeFig(fig, paths[1], filetag)
        yieldfile = os.path.join(dest, '{}{}.yield'.format(filetag, num))
        with open(yieldfile, 'w') as f:
            f.write('# {}\n'.format(paths[1]))
            f.write('# time: {}\n'.format(time))
            f.write('\n'.join(['{} {}'.format(*a) for a in massfiletuples]))
            f.write('\n')
        print("Wrote: {}".format(yieldfile))
        decayfile = os.path.join(dest, '{}{}.decayed'.format(filetag, num))
        with open(decayfile, 'w') as f:
            f.write('# {}\n'.format(paths[1]))
            f.write('# time: {}\n'.format(time))
            nmass = ['{} {}'.format(*a) for a in zip(names, decmasses)]
            f.write('\n'.join(nmass))
            f.write('\n')
        print("Wrote: {}".format(decayfile))


def plainTprofile(fname, thresh=1e-4, xrange=[0.0, 0.0],
                  byM=True, merged=False):
    """plots main properties of a plain text file.

    Args:
        fname (str): filename of checkpoint.
        thresh (float): ymin for species fraction plot.
        xrange (list of float): if set, change the range of the plots.
        byM (bool): abundance plot xaxis (by Mass or by Radius).

    """
    prof = DataMatrix(fname)
    if merged:
        return plotDMatMerged(prof, thresh=thresh, xrange=xrange, byM=byM)
    else:
        return plotDMat(prof, thresh=thresh, xrange=xrange, byM=byM)


def plotDMat(prof, thresh=1e-4, xrange=[0.0, 0.0], byM=True):
    """plots main properties of a profile object.

    Args:
        prof (DataMatrix): DataMatrix obj.
        thresh (float): ymin for species fraction plot.
        xrange (list of float): if set, change the range of the plots.
        byM (bool): abundance plot xaxis (by Mass or by Radius).

    """
    fig = plt.figure()
    skip = ['radius', 'masses', 'density']
    plotp = [x for x in prof.bulkprops if x not in skip]
    keys = sort_nuclides(prof.species)
    ncol = 4
    labelspace = -0.15

    if byM:
        xs = prof.masses
        xlab = 'Mass ($M_{\odot}$)'
        log = False
    else:
        xs = prof.radius
        xlab = 'Radius ($cm$)'
        log = True

    layout = (len(plotp)+2, 3)
    # Species
    ax1 = plt.subplot2grid(layout, (0, 0), colspan=2)
    skip = plotSpecies(ax1, prof, byMass=byM, thresh=thresh, plotall=False)
    # remove last(lowest) yticklabel to avoid overlap
    ax1.get_yticklabels()[1].set_visible(False)
    lgd = ax1.legend(ncol=ncol, loc='upper left', bbox_to_anchor=(1.00, 1.05),
                     columnspacing=0.0, labelspacing=0.0, markerfirst=False,
                     numpoints=3, handletextpad=0.0)
    ax1.axhline(1e0, linewidth=1, linestyle=':', color='black')
    ax1.set_ylim(thresh, 2.0)
    ax1.set_ylabel('Mass Frac.($X_{i}$)', size=13, rotation=90, labelpad=0)
    ax1.yaxis.set_label_coords(labelspace, 0.5)
    ax1.yaxis.set_minor_formatter(StrMethodFormatter(''))
    ax1.tick_params(labelbottom=False)
    if log:
        ax1.set_yscale('log')
        ax1.set_xscale('log')

    ax2 = plt.subplot2grid(layout, (1, 0), sharex=ax1, colspan=2)
    ax2.semilogy(xs, prof.density)
    ax2.set_ylabel('Density($\\frac{g}{cm^3}$)',
                   size=13, rotation=90, labelpad=0)
    ax2.yaxis.set_label_coords(labelspace, 0.5)
    ax2.yaxis.set_minor_formatter(StrMethodFormatter(''))
    ax2.tick_params(labelbottom=False)

    for i, p in enumerate(plotp):
        ax3 = plt.subplot2grid(layout, (i+2, 0), sharex=ax1, colspan=2)
        for j in range(i+1):
            ax3.plot([], [])
        ax3.semilogy(xs, getattr(prof, p))  # , color='k')
        u = ut.get_unit(p)
        spacer = labelspace if '\\' in u else labelspace - 0.02
        ax3.set_ylabel('{}({})'.format(p.capitalize(), u),
                       rotation=90, size=13, labelpad=0)
        ax3.yaxis.set_minor_formatter(StrMethodFormatter(''))
        ax3.yaxis.set_label_coords(spacer, 0.5)
        ax3.tick_params(labelbottom=False)
    ax3.set_xlabel(xlab)
    ax3.tick_params(labelbottom=True)
    if sum(xrange) != 0.0:
        ax1.set_xlim(xrange)
    plt.subplots_adjust(hspace=0.001, wspace=0.0)
    plt.subplots_adjust(left=0.13, right=0.80)
    plt.subplots_adjust(top=0.99, bottom=0.10)
    fig.set_size_inches(8.5, 7.5, forward=True)
    return fig


def plotDMatMerged(prof, thresh=1e-4, xrange=[0.0, 0.0],
                   byM=True, alpha=1.0, marker=False, botplotshift=5):
    """plots main properties of a profile in only two axes, merging
    thermodynamic properties.

    Args:
        prof (DataMatrix): DataMatrix obj.
        thresh (float): ymin for species fraction plot.
        xrange (list of float): if set, change the range of the plots.
        byM (bool): abundance plot xaxis (by Mass or by Radius).
        alpha(float): species plot alpha.
        marker(bool): enable plot cell marker.
        botplotshift(int): order of magnitude shift for bottom plot.

    """
    fig = plt.figure()
    skip = ['radius', 'masses', 'density']
    plotp = [x for x in prof.bulkprops if x not in skip]
    keys = sort_nuclides(prof.species)
    ncol = 4
    labelspace = -0.1
    if byM:
        xs = prof.masses
        xlab = 'Mass ($M_{\odot}$)'
        log = False
    else:
        xs = prof.radius
        xlab = 'Radius ($cm$)'
        log = True

    layout = (2, 2)
    # Species
    spax = plt.subplot2grid(layout, (0, 0), colspan=2)
    skip = plotSpecies(spax, prof, byMass=byM, thresh=thresh,
                       plotall=False, alpha=alpha, marker=marker)
    # remove last(lowest) yticklabel to avoid overlap
    spax.get_yticklabels()[1].set_visible(False)
    # remove xticklabels to avoid peeking to the sides of plot
    for lab in spax.get_xticklabels():
        lab.set_visible(False)
    # for logscale clearing xticklabels is not enough
    for lab in spax.get_xminorticklabels():
        lab.set_visible(False)
    for lab in spax.get_xmajorticklabels():
        lab.set_visible(False)
    lgd = spax.legend(ncol=ncol, loc='upper left', bbox_to_anchor=(1.00, 1.02),
                      columnspacing=0.0, labelspacing=0.0, markerfirst=False,
                      numpoints=3, handletextpad=0.0, edgecolor='k')
    spax.axhline(1e0, linewidth=1, linestyle=':', color='black')
    spax.set_ylim(thresh, 2.0)
    spax.set_ylabel('Mass Frac.($X_{i}$)', size=13, rotation=90, labelpad=0)
    spax.yaxis.set_label_coords(labelspace, 0.5)
    spax.yaxis.set_minor_formatter(StrMethodFormatter(''))
    spax.tick_params(labelbottom=False)
    spax.set_xlabel('')
    if log:
        spax.set_yscale('log')
        spax.set_xscale('log')
    # Thermodynamical variables (reference is, as always, density)
    if marker:
        mark = '.'
    else:
        mark = None
    tdax = plt.subplot2grid(layout, (1, 0), colspan=2, sharex=spax)
    tdax.semilogy(xs, prof.density, label='Density', marker=mark)
    ylabels = ['$\\frac{g}{cm^3}$']
    tdax.yaxis.set_minor_formatter(StrMethodFormatter(''))
    tdax.tick_params(labelbottom=False)
    for i, p in enumerate(plotp):
        u = ut.get_unit(p)
        if p == 'pressure':
            labl = p.capitalize()+'$\cdot 10^{{-18}}$'
            tdax.plot(xs, getattr(prof, p)/1e18, label=labl, marker=mark)
            ylabels.append(u)
        else:
            attv = getattr(prof, p)
            if p == 'gpot' or p == 'gpol':
                attv = -attv
                print('inverting potential energy')
            if np.max(attv) <= 0:
                print(p + ' skipped due to <=0 everywhere.')
                continue
            expo = np.floor(np.log10(np.max(attv)))-botplotshift
            if expo < 0:
                printfactor = '$\cdot 10^{{{}}}$'
            else:
                printfactor = '$\cdot 10^{{-{}}}$'
            labl = p.capitalize()+printfactor.format(int(abs(expo)))
            tdax.plot(xs, attv/np.power(10, expo),
                      label=labl, marker=mark)
            ylabels.append(u)
    tdax.yaxis.set_minor_formatter(StrMethodFormatter(''))
    tdax.tick_params(labelbottom=False)
    tdax.set_xlabel(xlab)
    tdax.tick_params(labelbottom=True)
    tdax.set_ylabel('({})'.format(','.join(ylabels)),
                    size=13, rotation=90, labelpad=0)
    tdax.yaxis.set_label_coords(labelspace, 0.5)
    lgd = tdax.legend(ncol=1, loc='upper left', bbox_to_anchor=(1.00, 0.50),
                      columnspacing=0.0, labelspacing=0.0, markerfirst=False,
                      numpoints=3, handletextpad=0.0, edgecolor='k')
    if sum(xrange) != 0.0:
        spax.set_xlim(xrange)
    plt.subplots_adjust(hspace=0.001, wspace=0.01)
    plt.subplots_adjust(left=0.13, right=0.80)
    plt.subplots_adjust(top=0.99, bottom=0.10)
    fig.set_size_inches(8.5, 7.5, forward=True)
    return fig


def plotDegen(prof, thresh=1e-4, xrange=[0.0, 0.0], byM=True):
    """plots composition and degeneracy parameters for a lineout,
    namely Y_e and \eta = T / T_fermi.
    \eta near zero implies degeneracy,
    while \eta >> 1 or negative implies non-degenerate matter.

    Args:
        prof (DataMatrix): DataMatrix obj.
        thresh (float): ymin for species fraction plot.
        xrange (list of float): if set, change the range of the plots.
        byM (bool): abundance plot xaxis (by Mass or by Radius).

    Returns:
        (mpl figure) or (None)

    """
    fig = plt.figure()
    skip = ['radius', 'masses', 'density']
    plotp = [x for x in prof.bulkprops if x not in skip]
    keys = sort_nuclides(prof.species)
    ncol = 4
    labelspace = -0.1
    if byM:
        xs = prof.masses
        xlab = 'Mass ($M_{\odot}$)'
        log = False
    else:
        xs = prof.radius
        xlab = 'Radius ($cm$)'
        log = True

    layout = (2, 2)
    # Species
    spax = plt.subplot2grid(layout, (0, 0), colspan=2)
    skip = plotSpecies(spax, prof, byMass=byM, thresh=thresh, plotall=False)
    # remove last(lowest) yticklabel to avoid overlap
    spax.get_yticklabels()[1].set_visible(False)
    lgd = spax.legend(ncol=ncol, loc='upper left', bbox_to_anchor=(1.00, 1.02),
                      columnspacing=0.0, labelspacing=0.0, markerfirst=False,
                      numpoints=3, handletextpad=0.0, edgecolor='k')
    spax.axhline(1e0, linewidth=1, linestyle=':', color='black')
    spax.set_ylim(thresh, 2.0)
    spax.set_ylabel('Mass Frac.($X_{i}$)', size=13, rotation=90, labelpad=0)
    spax.yaxis.set_label_coords(labelspace, 0.5)
    spax.yaxis.set_minor_formatter(StrMethodFormatter(''))
    spax.tick_params(labelbottom=False)
    if log:
        spax.set_yscale('log')
        spax.set_xscale('log')
    # Ye and Yion
    yes, yis = [], []
    npnts = len(getattr(prof, prof.species[0]))
    for i in range(npnts):
        xis = []
        for s in prof.species:
            xis.append(getattr(prof, s)[i])
        invyi, invye = getMus(prof.species, xis)
        yes.append(1.0/invye)
        yis.append(1.0/invyi)

    yeax = plt.subplot2grid(layout, (1, 0), colspan=2, sharex=spax)
    pl1 = yeax.semilogy(xs, yes, label='$Y_e$', color='#0082c8', alpha=0.8)
    yeax.set_ylabel('$Y_e$', size=13, rotation=0, labelpad=10)
    yeax.yaxis.set_major_formatter(StrMethodFormatter('{x:.3f}'))
    offset = yeax.get_yaxis().get_offset_text()
    offset.set_visible(False)
    yeax.set_ylim(0.1, 1.2)
    # yeax.get_yticklabels()[1].set_visible(False)
    # print(yeax.get_yticklabels())
    # yeax.get_yticklabels()[0].set_visible(False)
    # temps
    tdax = yeax.twinx()
    # get fermi temperatures through the lineout
    # fermT = [ fp.extRelFermi(d)/ut.kb for d in prof.density ]
    # tdax.plot(xs, prof.temperature/fermT, label='extreme fermT')#, color='k')
    fermT = [fp.nonRelFermi(d, ye=y)/ut.kb for d, y in zip(prof.density, yes)]
    pl2 = tdax.semilogy(xs, prof.temperature/fermT,
                        label=r'$\eta=\frac{T}{T_{f}}$',
                        color='#3cb44b', alpha=0.8)
    tdax.yaxis.set_minor_formatter(StrMethodFormatter(''))
    tdax.set_ylabel('$\eta$', size=13, rotation=0, labelpad=8)
    plts = pl1 + pl2
    labels = [l.get_label() for l in plts]
    lgd = tdax.legend(plts, labels, ncol=1, loc='upper left',
                      labelspacing=0.0, markerfirst=False,
                      numpoints=3, handletextpad=0.0, frameon=False)
    yeax.set_xlabel(xlab)
    if sum(xrange) != 0.0:
        spax.set_xlim(xrange)
    plt.subplots_adjust(hspace=0.001, wspace=0.0)
    plt.subplots_adjust(left=0.13, right=0.80)
    plt.subplots_adjust(top=0.99, bottom=0.10)
    fig.set_size_inches(8.5, 7.5, forward=True)
    return fig


def plotSpecies(ax, dmatr, byMass=True, thresh=1e-4, log=True,
                plotall=False, alpha=1.0, marker=False, ls='-',
                offset=0.0):
    """draws species from a profile object.

    Args:
        ax(mpl.axes): axes instance to draw on.
        dmatr(DataMatrix): profile object to extract data.
        byMass(bool): plot by Mass or radius (toggle).
        thresh(float): lower bound for plot.
        log(bool): x-axis log toggle.
        plotall(bool): plot species even when under 'thresh'.
        alpha(float): transparency for all lines plotted.
        marker(bool): mark points in graph.
        ls(str): force linestyle for 'plotall'.
        offset(float): abscissa offset.

    """
    if byMass:
        absc = 'mass'
    else:
        absc = 'radius'
    if marker:
        mark = '.'
    else:
        mark = None
    skip = []
    if plotall:
        for s in dmatr.species:
            line = simplePlot(ax, dmatr, absc, s, marker=mark,
                              delt=offset, log=log,
                              alpha=alpha, ls=ls)
            line[0].set_label(s)
    else:
        for i, s in enumerate(dmatr.species):
            props = next(ax._get_lines.prop_cycler)
            if np.max(getattr(dmatr, s)) < thresh:
                skip.append(i)
                continue
            else:
                tag = '$^{{{}}}{}$'.format(*elemSplit(s, invert=True))
                line = simplePlot(ax, dmatr, absc, s,
                                  marker=mark, delt=offset,
                                  log=log, color=props['color'],
                                  ls=props['linestyle'], alpha=alpha)
                line[0].set_label(tag)
    l = ax.set_ylabel('$X_i$')
    l = ax.set_ylim([thresh, 2.0])
    return skip


def simplePlot(ax, dmatr, absc, attr, delt=0.0, log=True, **kwargs):
    """draws a pair of properties from a profile object.

    Args:
        ax(mpl.axes): axes instance to draw on.
        dmatr(DataMatrix): profile object to extract data.
        absc(str): abscissa for the plot.
        attr(str): attribute to plot.

    """
    x = np.array(getattr(dmatr, absc))
    y = np.array(getattr(dmatr, attr))
    if log:
        # print(absc, delt)
        # print(getattr(dmatr, absc))
        # print(type(getattr(dmatr, absc)), type(delt))
        line = plt.loglog(x + delt, y, **kwargs)
    else:
        line = plt.plot(getattr(dmatr, absc) + delt,
                        getattr(dmatr, attr), **kwargs)
    l = plt.xlabel(absc.capitalize())
    l = plt.ylabel(attr.capitalize())
    return line


def drawGrid(ax, gridpoints, alpha=0.6, color='salmon', lw=2.0):
    """draws gridpoints as a line grid on an axes."""
    print('Gridpoints: {}'.format(len(gridpoints)))
    for r in gridpoints:
        ax.axvline(r, alpha=alpha, color=color, lw=lw)

# This file is part of meas_extensions_psfex.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os
import re

import numpy as np
from astropy.io import fits
try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

import lsst.geom as geom
import lsst.afw.geom as afwGeom
from lsst.afw.fits import readMetadata
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.display as afwDisplay
from . import psfexLib

afwDisplay.setDefaultMaskTransparency(75)


def compute_fwhmrange(fwhm, maxvar, minin, maxin, plot=dict(fwhmHistogram=False)):
    """Compute the FWHM range associated to a series of FWHM measurements.
        AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
        VERSION 20/03/2008

    Parameters
    ----------
    fwhm: iterable of `float`
        Iterable of full width half-maxima.
    maxvar: `float`
        Maximum allowed FWHM variation.
    minin: `float`
        Minimum allowed FWHM.
    maxin: `float`
        Maximum allowed FWHM.
    plot: `dict`, optional
        Dict of plotting options.

    Returns
    -------
    fmin: `float`
        FWHM mode.
    minout: `float`
        Lower FWHM range.
    maxout: `float`
        Upper FWHM range.
    """

    nfwhm = len(fwhm)
    fwhm.sort()

    # Find the mode
    nw = nfwhm//4
    if nw < 4:
        nw = 1
    dfmin = psfexLib.BIG
    fmin = 0.0
    for i in range(nfwhm - nw):
        df = fwhm[i + nw] - fwhm[i]
        if df < dfmin:
            dfmin = df
            fmin = (fwhm[i + nw] + fwhm[i])/2.0

    if nfwhm < 2:
        fmin = fwhm[0]

    dfmin = (maxvar + 1.0)**0.3333333
    minout = fmin/dfmin if dfmin > 0.0 else 0.0
    if minout < minin:
        minout = minin

    maxout = fmin*dfmin**2
    if maxout > maxin:
        maxout = maxin

    if plt and plot.get("fwhmHistogram"):
        plt.clf()
        plt.hist(fwhm, nfwhm//10 + 1, normed=1, facecolor='g', alpha=0.75)
        plt.xlabel("FWHM")
        plt.axvline(fmin, color='red')
        [plt.axvline(_, color='blue') for _ in (minout, maxout)]

        input("Continue? ")

    return fmin, minout, maxout


def select_candidates(set, prefs, frmin, frmax,
                      flags, flux, fluxerr, rmsSize, elong, vignet,
                      plot=dict(), title=""):
    maxbad = prefs.getBadpixNmax()
    maxbadflag = prefs.getBadpixFlag()
    maxelong = (prefs.getMaxellip() + 1.0)/(1.0 - prefs.getMaxellip()) if prefs.getMaxellip() < 1.0 else 100.0
    minsn = prefs.getMinsn()

    sn = flux/np.where(fluxerr > 0, fluxerr, 1)
    sn[fluxerr <= 0] = -psfexLib.BIG
    # ---- Apply some selection over flags, fluxes...
    plotFlags = plot.get("showFlags") if plt else False
    plotRejection = plot.get("showRejection") if plt else False

    bad = flags & prefs.getFlagMask() != 0
    set.setBadFlags(int(sum(bad != 0)))

    if plotRejection:
        selectionVectors = []
        selectionVectors.append((bad, "flags %d" % sum(bad != 0)))

    dbad = sn < minsn
    set.setBadSN(int(sum(dbad)))
    bad = np.logical_or(bad, dbad)
    if plotRejection:
        selectionVectors.append((dbad, "S/N %d" % sum(dbad)))

    dbad = rmsSize < frmin
    set.setBadFrmin(int(sum(dbad)))
    bad = np.logical_or(bad, dbad)
    if plotRejection:
        selectionVectors.append((dbad, "frmin %d" % sum(dbad)))

    dbad = rmsSize > frmax
    set.setBadFrmax(int(sum(dbad)))
    bad = np.logical_or(bad, dbad)
    if plotRejection:
        selectionVectors.append((dbad, "frmax %d" % sum(dbad)))

    dbad = elong > maxelong
    set.setBadElong(int(sum(dbad)))
    bad = np.logical_or(bad, dbad)
    if plotRejection:
        selectionVectors.append((dbad, "elong %d" % sum(dbad)))

    # -- ... and check the integrity of the sample
    if maxbadflag:
        nbad = np.array([(v <= -psfexLib.BIG).sum() for v in vignet])
        dbad = nbad > maxbad
        set.setBadPix(int(sum(dbad)))
        bad = np.logical_or(bad, dbad)
        if plotRejection:
            selectionVectors.append((dbad, "badpix %d" % sum(dbad)))

    good = np.logical_not(bad)
    if plotFlags or plotRejection:
        imag = -2.5*np.log10(flux)
        plt.clf()

        alpha = 0.5
        if plotFlags:
            labels = getFlags()

            isSet = np.where(flags == 0x0)[0]
            plt.plot(imag[isSet], rmsSize[isSet], 'o', alpha=alpha, label="good")

            for i in range(16):
                mask = 1 << i
                if mask & prefs.getFlagMask():
                    isSet = np.where(np.bitwise_and(flags, mask))[0]
                    if isSet.any():
                        plt.plot(imag[isSet], rmsSize[isSet], 'o', alpha=alpha, label=labels[mask])
        else:
            for bad, label in selectionVectors:
                plt.plot(imag[bad], rmsSize[bad], 'o', alpha=alpha, label=label)

        plt.plot(imag[good], rmsSize[good], 'o', color="black", label="selected")
        [plt.axhline(_, color='red') for _ in [frmin, frmax]]
        plt.xlim(np.median(imag[good]) + 5*np.array([-1, 1]))
        plt.ylim(-0.1, 2*frmax)
        plt.legend(loc=2)
        plt.xlabel("Instrumental Magnitude")
        plt.ylabel("rmsSize")
        plt.title("%s %d selected" % (title, sum(good)))

        input("Continue? ")

    return good


def showPsf(psf, set, ext=None, wcsData=None, trim=0, nspot=5,
            diagnostics=False, outDir="", frame=None, title=None):
    """Show a PSF on display (e.g., ds9)
    """

    if ext is not None:
        psf = psf[ext]

    if wcsData:
        if ext is not None:
            wcsData = wcsData[ext]
        wcs, naxis1, naxis2 = wcsData
    else:
        wcs, naxis1, naxis2 = None, None, None

    naxis = [naxis1, naxis2]
    for i in range(2):
        if naxis[i] is None:
            # cmin, cmax are the range of input star positions
            cmin, cmax = [set.getContextOffset(i) + d*set.getContextScale(i) for d in (-0.5, 0.5)]
            naxis[i] = cmax + cmin          # a decent guess

    if naxis[0] > naxis[1]:
        nx, ny = int(nspot*naxis[0]/float(naxis[1]) + 0.5), nspot
    else:
        nx, ny = nspot, int(nspot*naxis[1]/float(naxis[0]) + 0.5)

    mos = afwDisplay.utils.Mosaic(gutter=2, background=0.02)

    xpos, ypos = np.linspace(0, naxis[0], nx), np.linspace(0, naxis[1], ny)
    for y in ypos:
        for x in xpos:
            psf.build(x, y)

            im = afwImage.ImageF(*psf.getLoc().shape)
            im.getArray()[:] = psf.getLoc()
            im /= float(im.getArray().max())
            if trim:
                if trim > im.getHeight()//2:
                    trim = im.getHeight()//2

                im = im[trim:-trim, trim:-trim]

            mos.append(im)

    mosaic = mos.makeMosaic(mode=nx)
    if frame is not None:
        afwDisplay.Display(frame=frame).mtv(mosaic, title=title)
    #
    # Figure out the WCS for the mosaic
    #
    pos = []
    pos.append([geom.PointD(0, 0), wcs.pixelToSky(geom.PointD(0, 0))])
    pos.append([geom.PointD(*mosaic.getDimensions()), wcs.pixelToSky(geom.PointD(naxis1, naxis2))])

    CD = []
    for i in range(2):
        delta = pos[1][1][i].asDegrees() - pos[0][1][i].asDegrees()
        CD.append(delta/(pos[1][0][i] - pos[0][0][i]))
    CD = np.array(CD)
    CD.shape = (2, 2)
    mosWcs = afwGeom.makeSkyWcs(crval=pos[0][0], crpix=pos[0][1], cdMatrix=CD)

    if ext is not None:
        title = "%s-%d" % (title, ext)

    if frame is not None:
        afwDisplay.Display(frame=frame).mtv(mosaic, title=title, wcs=mosWcs)

    if diagnostics:
        outFile = "%s-mod.fits" % title
        if outDir:
            outFile = os.path.join(outDir, outFile)
        mosaic.writeFits(outFile, mosWcs.getFitsMetadata())

    mos = afwDisplay.utils.Mosaic(gutter=4, background=0.002)
    for i in range(set.getNsample()):
        s = set.getSample(i)
        if ext is not None and s.getExtindex() != ext:
            continue

        smos = afwDisplay.utils.Mosaic(gutter=2, background=-0.003)
        for func in [s.getVig, s.getVigResi]:
            arr = func()
            if func == s.getVig:
                norm = float(arr.max()) if True else s.getNorm()

            arr /= norm
            im = afwImage.ImageF(*arr.shape)
            im.getArray()[:] = arr
            smos.append(im)

        mos.append(smos.makeMosaic(mode="x"))

    mosaic = mos.makeMosaic(title=title)

    if frame is not None:
        afwDisplay.Display(frame=frame + 1).mtv(mosaic, title=title)

    if diagnostics:
        outFile = "%s-psfstars.fits" % title
        if outDir:
            outFile = os.path.join(outDir, outFile)

        mosaic.writeFits(outFile)


def getFlags(tab=None):
    flagKeys = [
        "base_PixelFlags_flag_edge",
        # "base_PixelFlags_flag_interpolated",
        # "base_PixelFlags_flag_interpolatedCenter",
        # "base_PixelFlags_flag_saturated",
        "base_PixelFlags_flag_saturatedCenter",
        # "base_PixelFlags_flag_cr",
        "base_PixelFlags_flag_crCenter",
        "base_PixelFlags_flag_bad",
        "base_PsfFlux_flag",
        "parent",
    ]

    if tab is None:
        flags = {}
        for i, k in enumerate(flagKeys):
            flags[1 << i] = re.sub(r"\_flag", "",
                                   re.sub(r"^base\_", "", re.sub(r"^base\_PixelFlags\_flag\_", "", k)))
    else:
        flags = 0
        for i, k in enumerate(flagKeys):
            if k == "parent":
                try:
                    isSet = tab.get("deblend_nChild") > 0
                except KeyError:
                    isSet = 0
            else:
                isSet = tab.get(k)
            flags = np.bitwise_or(flags, np.where(isSet, 1 << i, 0))

    return flags


def guessCalexp(fileName):
    for guess in [
        re.sub("/src", r"", fileName),
        re.sub("(SRC([^.]+))", r"CORR\2-exp", fileName),
    ]:
        if guess != fileName and os.path.exists(guess):
            return guess

    raise RuntimeError("Unable to find a calexp to go with %s" % fileName)


def makeitLsst(prefs, context, saveWcs=False, plot=dict()):
    """This is the python wrapper that reads lsst tables
    """
    # Create an array of PSFs (one PSF for each extension)
    if prefs.getVerboseType() != prefs.QUIET:
        print("----- %d input catalogues:" % prefs.getNcat())

    if saveWcs:                         # only needed for making plots
        wcssList = []

    fields = psfexLib.vectorField()
    for cat in prefs.getCatalogs():
        field = psfexLib.Field(cat)
        wcss = []
        wcssList.append(wcss)
        with fits.open(cat):
            # Hack: I want the WCS so I'll guess where the calexp is to be
            # found
            calexpFile = guessCalexp(cat)
            md = readMetadata(calexpFile)
            wcs = afwGeom.makeSkyWcs(md)

            if not wcs:
                cdMatrix = np.array([1.0, 0.0, 0.0, 1.0])
                cdMatrix.shape = (2, 2)
                wcs = afwGeom.makeSkyWcs(crpix=geom.PointD(0, 0),
                                         crval=geom.SpherePoint(0.0, 0.0, geom.degrees),
                                         cdMatrix=cdMatrix)

            naxis1, naxis2 = md.getScalar("NAXIS1"), md.getScalar("NAXIS2")
            # Find how many rows there are in the catalogue
            md = readMetadata(cat)

            field.addExt(wcs, naxis1, naxis2, md.getScalar("NAXIS2"))
            if saveWcs:
                wcss.append((wcs, naxis1, naxis2))

        field.finalize()
        fields.append(field)

    fields[0].getNext()          # number of extensions

    prefs.getPsfStep()

    sets = psfexLib.vectorSet()
    for set in load_samplesLsst(prefs, context, plot=plot):
        sets.append(set)

    psfexLib.makeit(fields, sets)

    ret = [[f.getPsfs() for f in fields], sets]
    if saveWcs:
        ret.append(wcssList)

    return ret


def read_samplesLsst(prefs, set, filename, frmin, frmax, ext, next, catindex, context, pcval, nobj,
                     plot=dict(showFlags=False, showRejection=False)):
    # allocate a new set iff set is None
    if not set:
        set = psfexLib.Set(context)

    cmin, cmax = None, None
    if set.getNcontext():
        cmin = np.empty(set.getNcontext())
        cmax = np.empty(set.getNcontext())
        for i in range(set.getNcontext()):
            if set.getNsample():
                cmin[i] = set.getContextOffset(i) - set.getContextScale(i)/2.0
                cmax[i] = cmin[i] + set.getContextScale(i)
            else:
                cmin[i] = psfexLib.BIG
                cmax[i] = -psfexLib.BIG
    #
    # Read data
    #
    tab = afwTable.SourceCatalog.readFits(filename)

    centroid = tab.getCentroidDefinition()
    xm = tab.get("%s.x" % centroid)
    ym = tab.get("%s.y" % centroid)

    shape = tab.getShapeDefinition()
    ixx = tab.get("%s.xx" % shape)
    iyy = tab.get("%s.yy" % shape)

    rmsSize = np.sqrt(0.5*(ixx + iyy))
    elong = 0.5*(ixx - iyy)/(ixx + iyy)

    flux = tab.get(prefs.getPhotfluxRkey())
    fluxErr = tab.get(prefs.getPhotfluxerrRkey())
    flags = getFlags(tab)

    #
    # Now the VIGNET data
    #
    vigw, vigh = 35, 35  # [prefs.getPsfsize()[i] for i in range(2)]
    if set.empty():
        set.setVigSize(vigw, vigh)

    vignet = np.empty(nobj*vigw*vigh, "float32").reshape(nobj, vigw, vigh)

    # Hack: I want the WCS so I'll guess where the calexp is to be found
    calexpFile = guessCalexp(filename)
    mi = afwImage.MaskedImageF(calexpFile)
    backnoise2 = np.median(mi.getVariance().getArray())
    gain = 1.0

    edgeBit = [k for k, v in getFlags().items() if v == "edge"][0]

    for i, xc, yc in zip(range(nobj), xm, ym):
        try:
            x, y = int(xc), int(yc)
        except ValueError:
            flags[i] |= edgeBit         # mark star as bad

        try:
            pstamp = mi[x - vigw//2:x + vigw//2 + 1, y - vigh//2:y + vigh//2 + 1]
            vignet[i] = pstamp.getImage().getArray().transpose()
        except Exception:
            flags[i] |= edgeBit         # mark star as bad

    # Try to load the set of context keys
    pc = 0
    contextvalp = []
    for i, key in enumerate(context.getName()):
        if context.getPcflag(i):
            contextvalp.append(pcval[pc])
            pc += 1
        elif key[0] == ':':
            try:
                contextvalp.append(tab.header[key[1:]])
            except KeyError:
                raise RuntimeError("*Error*: %s parameter not found in the header of %s" %
                                   (key[1:], filename))
        else:
            try:
                contextvalp.append(tab.get(key))
            except KeyError:
                raise RuntimeError("*Error*: %s parameter not found in the header of %s" %
                                   (key, filename))
            set.setContextname(i, key)

    # Now examine each vector of the shipment
    good = select_candidates(set, prefs, frmin, frmax,
                             flags, flux, fluxErr, rmsSize, elong, vignet,
                             plot=plot, title="%s[%d]" % (filename, ext + 1) if next > 1 else filename)
    #
    # Insert the good candidates into the set
    #
    if not vignet.dtype.isnative:
        # without the swap setVig fails with
        # "ValueError: 'unaligned arrays cannot be converted to C++'"
        vignet = vignet.byteswap()

    for i in np.where(good)[0]:
        sample = set.newSample()
        sample.setCatindex(catindex)
        sample.setExtindex(ext)

        sample.setVig(vignet[i])
        sample.setNorm(float(flux[i]))
        sample.setBacknoise2(backnoise2)
        sample.setGain(gain)
        sample.setX(float(xm[i]))
        sample.setY(float(ym[i]))
        sample.setFluxrad(float(rmsSize[i]))

        for j in range(set.getNcontext()):
            sample.setContext(j, float(contextvalp[j][i]))

        set.finiSample(sample, prefs.getProfAccuracy())

    # ---- Update min and max
    for j in range(set.getNcontext()):
        cmin[j] = contextvalp[j][good].min()
        cmax[j] = contextvalp[j][good].max()

    # Update the scaling
    if set.getNsample():
        for i in range(set.getNcontext()):
            set.setContextScale(i, cmax[i] - cmin[i])
            set.setContextOffset(i, (cmin[i] + cmax[i])/2.0)

    # Don't waste memory!
    set.trimMemory()

    return set


def load_samplesLsst(prefs, context, ext=psfexLib.Prefs.ALL_EXTENSIONS, next=1, plot=dict()):
    minsn = prefs.getMinsn()
    maxelong = (prefs.getMaxellip() + 1.0)/(1.0 - prefs.getMaxellip()) if prefs.getMaxellip() < 1.0 else 100
    min = prefs.getFwhmrange()[0]
    max = prefs.getFwhmrange()[1]

    filenames = prefs.getCatalogs()

    ncat = len(filenames)
    fwhmmin = np.empty(ncat)
    fwhmmax = np.empty(ncat)

    if not prefs.getAutoselectFlag():
        fwhmmin = np.zeros(ncat) + prefs.getFwhmrange()[0]
        fwhmmax = np.zeros(ncat) + prefs.getFwhmrange()[1]
        fwhmmode = (fwhmmin + fwhmmax)/2.0
    else:
        fwhms = {}

        # -- Try to estimate the most appropriate Half-light Radius range
        # -- Get the Half-light radii
        nobj = 0
        for i, fileName in enumerate(filenames):
            fwhms[i] = []

            if prefs.getVerboseType() != prefs.QUIET:
                print("Examining Catalog #%d" % (i+1))

            # ---- Read input catalog
            tab = afwTable.SourceCatalog.readFits(fileName)

            # -------- Fill the FWHM array
            shape = tab.getShapeDefinition()
            ixx = tab.get("%s.xx" % shape)
            iyy = tab.get("%s.yy" % shape)

            rmsSize = np.sqrt(0.5*(ixx + iyy))
            elong = 0.5*(ixx - iyy)/(ixx + iyy)

            flux = tab.get(prefs.getPhotfluxRkey())
            fluxErr = tab.get(prefs.getPhotfluxerrRkey())

            flags = getFlags(tab)

            good = np.logical_and(flux/fluxErr > minsn,
                                  np.logical_not(flags & prefs.getFlagMask()))
            good = np.logical_and(good, elong < maxelong)
            fwhm = 2.0*rmsSize
            good = np.logical_and(good, fwhm >= min)
            good = np.logical_and(good, fwhm < max)
            fwhms[i] = fwhm[good]

        if prefs.getVarType() == prefs.VAR_NONE:
            if nobj:
                fwhms_all = np.empty(sum([len(f) for f in fwhms.values()]))
                i = 0
                for f in fwhms.values():
                    fwhms_all[i:len(f)] = f
                    i += len(f)
                mode, min, max = compute_fwhmrange(fwhms_all, prefs.getMaxvar(),
                                                   prefs.getFwhmrange()[0], prefs.getFwhmrange()[1],
                                                   plot=plot)
            else:
                raise RuntimeError("No source with appropriate FWHM found!!")
                mode = min = max = 2.35/(1.0 - 1.0/psfexLib.cvar.INTERPFAC)

                fwhmmin = np.zeros(ncat) + min
                fwhmmax = np.zeros(ncat) + max
                fwhmmode = np.zeros(ncat) + mode
        else:
            fwhmmode = np.empty(ncat)
            fwhmmin = np.empty(ncat)
            fwhmmax = np.empty(ncat)

            for i in range(ncat):
                nobj = len(fwhms[i])
                if (nobj):
                    fwhmmode[i], fwhmmin[i], fwhmmax[i] = \
                        compute_fwhmrange(fwhms[i], prefs.getMaxvar(),
                                          prefs.getFwhmrange()[0], prefs.getFwhmrange()[1], plot=plot)
                else:
                    raise RuntimeError("No source with appropriate FWHM found!!")
                    fwhmmode[i] = fwhmmin[i] = fwhmmax[i] = 2.35/(1.0 - 1.0/psfexLib.cvar.INTERPFAC)

    # Read the samples
    mode = psfexLib.BIG               # mode of FWHM distribution

    sets = []
    for i, fileName in enumerate(filenames):
        set = None
        for ext in range(next):
            set = read_samplesLsst(prefs, set, fileName, fwhmmin[i]/2.0, fwhmmax[i]/2.0,
                                   ext, next, i, context,
                                   context.getPc(i) if context.getNpc() else None, nobj, plot=plot)

        if fwhmmode[i] < mode:
            mode = fwhmmode[i]

        set.setFwhm(mode)

        if prefs.getVerboseType() != prefs.QUIET:
            if set.getNsample():
                print("%d samples loaded." % set.getNsample())
            else:
                raise RuntimeError("No appropriate source found!!")

        sets.append(set)

    return sets

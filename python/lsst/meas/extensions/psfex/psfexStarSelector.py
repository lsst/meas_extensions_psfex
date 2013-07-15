# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import re
import sys

import numpy as np
try:
    import matplotlib.pyplot as plt
    fig = None
except ImportError:
    plt = None

import lsst.pex.config as pexConfig
import lsst.pex.logging as pexLogging
import lsst.afw.detection as afwDetection
import lsst.afw.display.ds9 as ds9
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as geomEllip
import lsst.afw.cameraGeom as cameraGeom
import lsst.meas.algorithms as measAlg
from lsst.meas.algorithms.starSelectorRegistry import starSelectorRegistry
import lsst.meas.extensions.psfex as psfex

class PsfexStarSelectorConfig(pexConfig.Config):
    badFlags = pexConfig.ListField(
        doc="List of flags which cause a source to be rejected as bad",
        dtype=str,
        default=["initial.flags.pixel.edge",
                 "initial.flags.pixel.saturated.center",
                 "initial.flags.pixel.cr.center",
                 "initial.flags.pixel.bad",
                 "initial.flags.pixel.suspect.center",
                 "initial.flux.psf.flags",
                 #"parent",            # actually this is a test on deblend.nchild
                 ],
        )
    fluxName = pexConfig.Field(
        dtype=str,
        doc="Name of photometric flux key ",
        default="initial.flux.psf",
        )
    fluxErrName = pexConfig.Field(
        dtype=str,
        doc="Name of phot. flux err. key",
        default="",
        )
    minFwhm = pexConfig.Field(
        dtype=float,
        doc="Maximum allowed FWHM ",
        default=2,
        )
    maxFwhm = pexConfig.Field(
        dtype=float,
        doc="Minimum allowed FWHM ",
        default=10,
        )
    maxFwhmVariability = pexConfig.Field(
        dtype=float,
        doc="Allowed FWHM variability (1.0 = 100%)",
        default=0.2,
        )
    maxbad = pexConfig.Field(
        dtype=int,
        doc="Max number of bad pixels ",
        default=0,
        check = lambda x: x >= 0,
        )
    maxbadflag = pexConfig.Field(
        dtype=bool,
        doc="Filter bad pixels? ",
        default=True
        )
    maxellip = pexConfig.Field(
        dtype=float,
        doc="Maximum (A-B)/(A+B) ",
        default=0.3,
        check = lambda x: x >= 0.0,
        )
    minsn = pexConfig.Field(
        dtype=float,
        doc="Minimum S/N for candidates",
        default=100,
        check = lambda x: x >= 0.0,
        )
    kernelSize = pexConfig.Field(
        dtype=int,
        doc = "size of the Psf kernel to create",
        default = 21,
        )
    borderWidth = pexConfig.Field(
        doc = "number of pixels to ignore around the edge of PSF candidate postage stamps",
        dtype = int,
        default = 0,
    )

    def validate(self):
        pexConfig.Config.validate(self)

        if self.fluxErrName == "":
            self.fluxErrName = self.fluxName + ".err"
        elif self.fluxErrName != self.fluxName + ".err":
            raise pexConfig.FieldValidationError("fluxErrName (%s) doesn't correspond to fluxName (%s)"
                                                 % (self.fluxErrName, self.fluxName))

        if self.minFwhm > self.maxFwhm:
            raise pexConfig.FieldValidationError("minFwhm (%f) > maxFwhm (%f)" % (self.minFwhm, self.maxFwhm))

class EventHandler(object):
    """A class to handle key strokes with matplotlib displays"""
    def __init__(self, axes, xs, ys, x, y, frames=[0]):
        self.axes = axes
        self.xs = xs
        self.ys = ys
        self.x = x
        self.y = y
        self.frames = frames

        self.cid = self.axes.figure.canvas.mpl_connect('key_press_event', self)

    def __call__(self, ev):
        if ev.inaxes != self.axes:
            return
        
        if ev.key and ev.key in ("p"):
            dist = np.hypot(self.xs - ev.xdata, self.ys - ev.ydata)
            dist[np.where(np.isnan(dist))] = 1e30
            dmin = min(dist)

            which = np.where(dist == min(dist))

            x = self.x[which][0]
            y = self.y[which][0]
            for frame in self.frames:
                ds9.pan(x, y, frame=frame)
            ds9.cmdBuffer.flush()
        else:
            pass

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def plot(mag, width, centers, clusterId, marker="o", markersize=2, markeredgewidth=0, ltype='-',
         clear=True):

    global fig
    if not fig:
        fig = plt.figure()
        newFig = True
    else:
        newFig = False
        if clear:
            fig.clf()

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));

    xmin = sorted(mag)[int(0.05*len(mag))]
    xmax = sorted(mag)[int(0.95*len(mag))]

    axes.set_xlim(-17.5, -13)
    axes.set_xlim(xmin - 0.1*(xmax - xmin), xmax + 0.1*(xmax - xmin))
    axes.set_ylim(0, 10)

    colors = ["r", "g", "b", "c", "m", "k",]
    for k, mean in enumerate(centers):
        if k == 0:
            axes.plot(axes.get_xlim(), (mean, mean,), "k%s" % ltype)

        l = (clusterId == k)
        axes.plot(mag[l], width[l], marker, markersize=markersize, markeredgewidth=markeredgewidth,
                  color=colors[k%len(colors)])

    l = (clusterId == -1)
    axes.plot(mag[l], width[l], marker, markersize=markersize, markeredgewidth=markeredgewidth,
              color='k')

    if newFig:
        axes.set_xlabel("model")
        axes.set_ylabel(r"$\sqrt{I_{xx} + I_{yy}}$")

    return fig
        
class PsfexStarSelector(object):
    ConfigClass = PsfexStarSelectorConfig

    def __init__(self, config):
        """Construct a star selector using psfex's algorithm
        
        @param[in] config: An instance of PsfexStarSelectorConfig
        """
        self.config = config
            
    def selectStars(self, exposure, catalog, matches=None):
        """Return a list of PSF candidates that represent likely stars
        
        A list of PSF candidates may be used by a PSF fitter to construct a PSF.
        
        @param[in] exposure: the exposure containing the sources
        @param[in] catalog: a SourceCatalog containing sources that may be stars
        @param[in] matches: astrometric matches; ignored by this star selector
        
        @return psfCandidateList: a list of PSF candidates.
        """
        import lsstDebug
        display = lsstDebug.Info(__name__).display

        displayExposure = display and \
            lsstDebug.Info(__name__).displayExposure # display the Exposure + spatialCells
        plotFwhmHistogram = display and plt and \
            lsstDebug.Info(__name__).plotFwhmHistogram # Plot histogram of FWHM
        plotFlags = display and plt and \
            lsstDebug.Info(__name__).plotFlags # Plot the sources coloured by their flags
        plotRejection = display and plt and \
            lsstDebug.Info(__name__).plotRejection # Plot why sources are rejected
        # create a log for my application
        logger = pexLogging.Log(pexLogging.getDefaultLog(), "meas.extensions.psfex.psfexStarSelector")

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        #
        fluxName = self.config.fluxName
        fluxErrName = self.config.fluxErrName
        minFwhm = self.config.minFwhm
        maxFwhm = self.config.maxFwhm
        maxFwhmVariability = self.config.maxFwhmVariability
        maxbad = self.config.maxbad
        maxbadflag = self.config.maxbadflag
        maxellip = self.config.maxellip
        minsn = self.config.minsn

        maxelong = (maxellip + 1.0)/(1.0 - maxellip) if maxellip < 1.0 else 100

        # Unpack the catalogue
        shape = catalog.getShapeDefinition()
        ixx = catalog.get("%s.xx" % shape)
        iyy = catalog.get("%s.yy" % shape)
        ixy = catalog.get("%s.xy" % shape)

        fwhm = 2*np.sqrt(2*np.log(2))*np.sqrt(0.5*(ixx + iyy))
        elong = 0.5*(ixx - iyy)/(ixx + iyy)

        flux = catalog.get(fluxName)
        fluxErr = catalog.get(fluxErrName)
        sn = flux/np.where(fluxErr > 0, fluxErr, 1)
        sn[fluxErr <= 0] = -psfex.psfex.cvar.BIG

        flags = 0x0
        for i, f in enumerate(self.config.badFlags):
            flags = np.bitwise_or(flags, np.where(catalog.get(f), 1 << i, 0))
        #
        # Estimate the acceptable range of source widths
        #
        good = np.logical_and(sn > minsn, np.logical_not(flags))
        good = np.logical_and(good, elong < maxelong)
        good = np.logical_and(good, fwhm >= minFwhm)
        good = np.logical_and(good, fwhm <  maxFwhm)

        fwhmMode, fwhmMin, fwhmMax = psfex.compute_fwhmrange(fwhm[good], maxFwhmVariability, minFwhm, maxFwhm,
                                                             plot=dict(fwhmHistogram=plotFwhmHistogram))

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        #
        # Here's select_candidates
        #
        #---- Apply some selection over flags, fluxes...

        bad = (flags != 0)
        #set.setBadFlags(int(sum(bad)))

        if plotRejection:
            selectionVectors = []
            selectionVectors.append((bad, "flags %d" % sum(bad)))

        dbad = sn < minsn
        #set.setBadSN(int(sum(dbad)))
        bad = np.logical_or(bad, dbad)
        if plotRejection:
            selectionVectors.append((dbad, "S/N %d" % sum(dbad)))

        dbad = fwhm < fwhmMin
        #set.setBadFrmin(int(sum(dbad)))
        bad = np.logical_or(bad, dbad)
        if plotRejection:
            selectionVectors.append((dbad, "fwhmMin %d" % sum(dbad)))

        dbad = fwhm > fwhmMax
        #set.setBadFrmax(int(sum(dbad)))
        bad = np.logical_or(bad, dbad)
        if plotRejection:
            selectionVectors.append((dbad, "fwhmMax %d" % sum(dbad)))

        dbad = elong > maxelong
        #set.setBadElong(int(sum(dbad)))
        bad = np.logical_or(bad, dbad)
        if plotRejection:
            selectionVectors.append((dbad, "elong %d" % sum(dbad)))

        #-- ... and check the integrity of the sample
        if maxbadflag:
            nbad = np.array([(v <= -psfex.psfex.cvar.BIG).sum() for v in vignet])
            dbad = nbad > maxbad
            #set.setBadPix(int(sum(dbad)))
            bad = np.logical_or(bad, dbad)
            if plotRejection:
                selectionVectors.append((dbad, "badpix %d" % sum(dbad)))

        good = np.logical_not(bad)
        #
        # We know enough to plot, if so requested
        #
        frame = 0
        if displayExposure:
            mi = exposure.getMaskedImage()
    
            ds9.mtv(mi, frame=frame, title="PSF candidates")

            with ds9.Buffering():
                for i, source in enumerate(catalog):
                    if good[i]:
                        ctype = ds9.GREEN # star candidate
                    else:
                        ctype = ds9.RED # not star

                    ds9.dot("+", source.getX() - mi.getX0(), source.getY() - mi.getY0(),
                            frame=frame, ctype=ctype)

        if plotFlags or plotRejection:
            imag = -2.5*np.log10(flux)
            plt.clf()

            alpha = 0.5
            if plotFlags:
                isSet = np.where(flags == 0x0)[0]
                plt.plot(imag[isSet], fwhm[isSet], 'o', alpha=alpha, label="good")

                for i, f in enumerate(self.config.badFlags):
                    mask = 1 << i
                    isSet = np.where(np.bitwise_and(flags, mask))[0]
                    if isSet.any():
                        if np.isfinite(imag[isSet] + fwhm[isSet]).any():
                            plt.plot(imag[isSet], fwhm[isSet], 'o', alpha=alpha,
                                     label=re.sub(r"^.*(flags\.pixel|flux)\.", "", f))
            else:
                for bad, label in selectionVectors:
                    plt.plot(imag[bad], fwhm[bad], 'o', alpha=alpha, label=label)

            plt.plot(imag[good], fwhm[good], 'o', color="black", label="selected")
            [plt.axhline(_, color='red') for _ in [fwhmMin, fwhmMax]]
            plt.xlim(np.median(imag[good]) + 5*np.array([-1, 1]))
            plt.ylim(fwhm[np.where(np.isfinite(fwhm + imag))].min(), 2*fwhmMax)
            plt.legend(loc=2)
            plt.xlabel("Instrumental %s Magnitude" % fluxName.split(".")[-1].title())
            plt.ylabel("fwhm")
            title = "PSFEX Star Selection"
            plt.title("%s %d selected" % (title, sum(good)))

        if displayExposure:
            global eventHandler
            eventHandler = EventHandler(plt.axes(), imag, fwhm, catalog.getX(), catalog.getY(), frames=[frame])

        if plotFlags or plotRejection:
            while True:
                try:
                    reply = raw_input("continue? [y[es] h(elp) p(db) q(uit)] ").strip()
                except EOFError:
                    reply = "y"

                if not reply:
                    reply = "y"

                if reply[0] == "h":
                    print """\
At this prompt, you can continue with almost any key; 'p' enters pdb,
                                                      'q' returns to the shell, and
                                                      'h' prints this text
""",

                    if displayExposure:
                        print """
If you put the cursor on a point in the matplotlib scatter plot and hit 'p' you'll see it in ds9."""
                elif reply[0] == "p":
                    import pdb; pdb.set_trace()
                elif reply[0] == 'q':
                    sys.exit(1)
                else:
                    break

        #
        # Time to use that stellar classification to generate psfCandidateList
        #
        with ds9.Buffering():
            psfCandidateList = []
            if True:
                catalog = [s for s,g in zip(catalog, good) if g]
            else:
                catalog = catalog[good]

            for source in catalog:
                try:
                    psfCandidate = measAlg.makePsfCandidate(source, exposure)
                    # The setXXX methods are class static, but it's convenient to call them on
                    # an instance as we don't know Exposure's pixel type
                    # (and hence psfCandidate's exact type)
                    if psfCandidate.getWidth() == 0:
                        psfCandidate.setBorderWidth(self.config.borderWidth)
                        psfCandidate.setWidth(self.config.kernelSize + 2*self.config.borderWidth)
                        psfCandidate.setHeight(self.config.kernelSize + 2*self.config.borderWidth)

                    im = psfCandidate.getMaskedImage().getImage()
                    vmax = afwMath.makeStatistics(im, afwMath.MAX).getValue()
                    if not np.isfinite(vmax):
                        continue
                    psfCandidateList.append(psfCandidate)

                    if display and displayExposure:
                        ds9.dot("o", source.getX() - mi.getX0(), source.getY() - mi.getY0(),
                                size=4, frame=frame, ctype=ds9.CYAN)
                except Exception as err:
                    logger.logdebug("Failed to make a psfCandidate from source %d: %s" % (source.getId(), err))

        return psfCandidateList

starSelectorRegistry.register("psfex", PsfexStarSelector)

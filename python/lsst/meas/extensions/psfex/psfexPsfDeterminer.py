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
import os
import sys
import numpy as np

import lsst.daf.base as dafBase
import lsst.pex.config as pexConfig
import lsst.pex.exceptions as pexExceptions
import lsst.pex.logging as pexLog
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as afwEll
import lsst.afw.detection as afwDetection
import lsst.afw.display.ds9 as ds9
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
import lsst.meas.algorithms.utils as maUtils
import lsst.meas.algorithms.psfDeterminerRegistry as psfDeterminerRegistry
import lsst.meas.extensions.psfex as psfex

class PsfexPsfDeterminerConfig(pexConfig.Config):
    __nEigenComponents = pexConfig.Field(
        doc = "number of eigen components for PSF kernel creation",
        dtype = int,
        default = 4,
    )
    spatialOrder = pexConfig.Field(
        doc = "specify spatial order for PSF kernel creation",
        dtype = int,
        default = 2,
        check = lambda x: x >= 0,
    )
    sizeCellX = pexConfig.Field(
        doc = "size of cell used to determine PSF (pixels, column direction)",
        dtype = int,
        default = 256,
#        minValue = 10,
        check = lambda x: x >= 10,
    )
    sizeCellY = pexConfig.Field(
        doc = "size of cell used to determine PSF (pixels, row direction)",
        dtype = int,
        default = sizeCellX.default,
#        minValue = 10,
        check = lambda x: x >= 10,
    )
    __nStarPerCell = pexConfig.Field(
        doc = "number of stars per psf cell for PSF kernel creation",
        dtype = int,
        default = 3,
    )
    kernelSize = pexConfig.Field(
        doc = "radius of the kernel to create, relative to the square root of the stellar quadrupole moments",
        dtype = float,
        default = 10.0,
    )
    kernelSizeMin = pexConfig.Field(
        doc = "Minimum radius of the kernel",
        dtype = int,
        default = 25,
    )
    kernelSizeMax = pexConfig.Field(
        doc = "Maximum radius of the kernel",
        dtype = int,
        default = 45,
    )
    __borderWidth = pexConfig.Field(
        doc = "Number of pixels to ignore around the edge of PSF candidate postage stamps",
        dtype = int,
        default = 0,
    )
    __nStarPerCellSpatialFit = pexConfig.Field(
        doc = "number of stars per psf Cell for spatial fitting",
        dtype = int,
        default = 5,
    )
    __constantWeight = pexConfig.Field(
        doc = "Should each PSF candidate be given the same weight, independent of magnitude?",
        dtype = bool,
        default = True,
    )
    __nIterForPsf = pexConfig.Field(
        doc = "number of iterations of PSF candidate star list",
        dtype = int,
        default = 3,
    )
    tolerance = pexConfig.Field(
        doc = "tolerance of spatial fitting",
        dtype = float,
        default = 1e-2,
    )
    lam = pexConfig.Field(
        doc = "floor for variance is lam*data",
        dtype = float,
        default = 0.05,
    )
    reducedChi2ForPsfCandidates = pexConfig.Field(
        doc = "for psf candidate evaluation",
        dtype = float,
        default = 2.0,
    )
    spatialReject = pexConfig.Field(
        doc = "Rejection threshold (stdev) for candidates based on spatial fit",
        dtype = float,
        default = 3.0,
    )

class PsfexPsfDeterminer(object):
    ConfigClass = PsfexPsfDeterminerConfig

    def __init__(self, config):
        """Construct a PSFEX PSF Fitter

        @param[in] config: instance of PsfexPsfDeterminerConfig
        """
        self.config = config
        # N.b. name of component is meas.algorithms.psfDeterminer so you can turn on psf debugging
        # independent of which determiner is active
        self.debugLog = pexLog.Debug("meas.algorithms.psfDeterminer")
        self.warnLog = pexLog.Log(pexLog.getDefaultLog(), "meas.algorithms.psfDeterminer")

    def determinePsf(self, exposure, psfCandidateList, metadata=None, flagKey=None):
        """Determine a PSFEX PSF model for an exposure given a list of PSF candidates
        
        @param[in] exposure: exposure containing the psf candidates (lsst.afw.image.Exposure)
        @param[in] psfCandidateList: a sequence of PSF candidates (each an lsst.meas.algorithms.PsfCandidate);
            typically obtained by detecting sources and then running them through a star selector
        @param[in,out] metadata  a home for interesting tidbits of information
        @param[in] flagKey: schema key used to mark sources actually used in PSF determination
    
        @return psf: a meas.extensions.psfex.PsfexPsf
        """
        import lsstDebug
        display = lsstDebug.Info(__name__).display 
        displayExposure = display and \
            lsstDebug.Info(__name__).displayExposure     # display the Exposure + spatialCells 
        displayPsfCandidates = display and \
            lsstDebug.Info(__name__).displayPsfCandidates # show the viable candidates 
        displayPsfComponents = display and \
            lsstDebug.Info(__name__).displayPsfComponents # show the basis functions
        showBadCandidates = display and \
            lsstDebug.Info(__name__).showBadCandidates # Include bad candidates (meaningless, methinks)
        displayResiduals = display and \
            lsstDebug.Info(__name__).displayResiduals         # show residuals
        displayPsfMosaic = display and \
            lsstDebug.Info(__name__).displayPsfMosaic   # show mosaic of reconstructed PSF(x,y)
        matchKernelAmplitudes = lsstDebug.Info(__name__).matchKernelAmplitudes # match Kernel amplitudes for spatial plots
        normalizeResiduals = lsstDebug.Info(__name__).normalizeResiduals # Normalise residuals by object amplitude 

        mi = exposure.getMaskedImage()
        
        nCand = len(psfCandidateList)
        if nCand == 0:
            raise RuntimeError("No PSF candidates supplied.")
        #
        # How big should our PSF models be?
        #
        if display:                     # only needed for debug plots
            # construct and populate a spatial cell set
            bbox = mi.getBBox(afwImage.PARENT)
            psfCellSet = afwMath.SpatialCellSet(bbox, self.config.sizeCellX, self.config.sizeCellY)
        
        sizes = np.empty(nCand)
        for i, psfCandidate in enumerate(psfCandidateList):
            try:
                if psfCellSet:
                    psfCellSet.insertCandidate(psfCandidate)
            except Exception, e:
                self.debugLog.debug(2, "Skipping PSF candidate %d of %d: %s" % (i, len(psfCandidateList), e))
                continue

            source = psfCandidate.getSource()
            quad = afwEll.Quadrupole(source.getIxx(), source.getIyy(), source.getIxy())
            rmsSize = quad.getTraceRadius()
            sizes[i] = rmsSize

        if self.config.kernelSize >= 15:
            self.debugLog.debug(1, \
                "WARNING: NOT scaling kernelSize by stellar quadrupole moment, but using absolute value")
            actualKernelSize = int(self.config.kernelSize)
        else:
            actualKernelSize = 2 * int(self.config.kernelSize * np.sqrt(np.median(sizes)) + 0.5) + 1
            if actualKernelSize < self.config.kernelSizeMin:
                actualKernelSize = self.config.kernelSizeMin
            if actualKernelSize > self.config.kernelSizeMax:
                actualKernelSize = self.config.kernelSizeMax
            if display:
                rms = np.median(sizes)
                print "Median PSF RMS size=%.2f pixels (\"FWHM\"=%.2f)" % (rms, 2*np.sqrt(2*np.log(2))*rms)
        self.debugLog.debug(3, "Kernel size=%s" % (actualKernelSize,))

        # Set size of image returned around candidate
        psfCandidateList[0].setHeight(actualKernelSize)
        psfCandidateList[0].setWidth(actualKernelSize)

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- BEGIN PSFEX
        #
        # Insert the good candidates into the set
        #
        defaultsFile = os.path.join(os.environ["MEAS_EXTENSIONS_PSFEX_DIR"], "config", "default-lsst.psfex")
        args_md = dafBase.PropertySet()
        args_md.set("PSFVAR_DEGREES", str(self.config.spatialOrder))
        prefs = psfex.Prefs(defaultsFile, args_md)
        prefs.setCommandLine([])
        prefs.addCatalog("psfexPsfDeterminer")

        prefs.use()

        context = psfex.Context(prefs.getContextName(), prefs.getContextGroup(),
                                prefs.getGroupDeg(),
                                psfex.Context.REMOVEHIDDEN if False else psfex.Context.KEEPHIDDEN)
	set = psfex.Set(context)
        set.setVigSize(actualKernelSize, actualKernelSize)
        set.setFwhm(2*np.sqrt(2*np.log(2))*np.median(sizes))

        catindex, ext = 0, 0
        backnoise2 = afwMath.makeStatistics(mi.getImage(), afwMath.VARIANCECLIP).getValue()
        gain = 1.0
        self.warnLog.log(pexLog.Log.WARN, "Setting gain to %g" % gain)

        contextvalp = []
        for i, key in enumerate(context.getName()):
            if context.getPcflag(i):
                contextvalp.append(pcval[pc])
                pc += 1
            elif key[0] == ':':
                try:
                    contextvalp.append(exposure.getMetadata().get(key[1:]))
                except KeyError:
                    raise RuntimeError("*Error*: %s parameter not found in the header of %s" %
                                       (key[1:], filename))
            else:
                try:
                    contextvalp.append(np.array([psfCandidateList[_].getSource().get(key)
                                                    for _ in range(nCand)]))
                except KeyError:
                    raise RuntimeError("*Error*: %s parameter not found" % (key))
                set.setContextname(i, key)

        if display:
            frame = 0
            if displayExposure:
                ds9.mtv(exposure, frame=frame, title="psf determination")
            
        with ds9.Buffering():
            xpos, ypos = [], []
            for i, psfCandidate in enumerate(psfCandidateList):
                source = psfCandidate.getSource()

                sample = set.newSample()
                sample.setCatindex(catindex)
                sample.setExtindex(ext)

                xc, yc = source.getX(), source.getY()
                try:
                    x, y = int(xc), int(yc)
                except ValueError:
                    continue

                try:
                    pstamp = mi[x - actualKernelSize//2:x + actualKernelSize//2 + 1,
                                y - actualKernelSize//2:y + actualKernelSize//2 + 1]
                    sample.setVig(pstamp.getImage().getArray().transpose().copy())
                except Exception, e:
                    print e
                    continue

                if False:
                    pstamp.getImage().getArray()[:] = sample.getVig()

                    ds9.mtv(pstamp)

                sample.setNorm(source.get(prefs.getPhotfluxRkey()))
                sample.setBacknoise2(backnoise2)
                sample.setGain(gain)
                sample.setX(xc)
                sample.setY(float(yc))
                sample.setFluxrad(sizes[i])

                for j in range(set.getNcontext()):
                    sample.setContext(j, float(contextvalp[j][i]))

                set.finiSample(sample)

                if flagKey is not None:
                    source.set(flagKey, True)

                xpos.append(xc); ypos.append(yc) # for QA

            if displayExposure:
                ds9.dot("o", xc, yc, ctype=ds9.CYAN, size=4, frame=frame)

        #---- Update min and max and then the scaling
        for i in range(set.getNcontext()):
            cmin = contextvalp[i].min()
            cmax = contextvalp[i].max()
            set.setContextScale(i, cmax - cmin)
            set.setContextOffset(i, (cmin + cmax)/2.0)

        # Don't waste memory!
        set.trimMemory()

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- END PSFEX 
        #
        # Do a PSFEX decomposition of those PSF candidates
        #
        fields = psfex.vectorField()
        field = psfex.Field("Unknown")
        field.addExt(exposure.getWcs(), exposure.getWidth(), exposure.getHeight(), set.getNsample())
        field.finalize()

        fields.append(field)

        sets = psfex.vectorSet()
        sets.append(set)

        psfex.makeit(fields, sets)
        psfs = field.getPsfs()

        xpos = np.array(xpos); ypos = np.array(ypos)
        numGoodStars = len(xpos)
        avgX, avgY = np.mean(xpos), np.mean(ypos)

        psf = psfex.PsfexPsf(psfs[0], afwGeom.Point2D(avgX, avgY))

        if False and (displayResiduals or displayPsfMosaic):
            ext = 0
            frame = 1
            diagnostics = True
            catDir = "."
            title = "psfexPsfDeterminer"
            psfex.psfex.showPsf(psfs, set, ext,
                                [(exposure.getWcs(), exposure.getWidth(), exposure.getHeight())],
                                nspot=3, trim=5, frame=frame, diagnostics=diagnostics, outDir=catDir,
                                title=title)
        #
        # Display code for debugging
        #
        if display:
            assert psfCellSet is not None

            if displayExposure:
                maUtils.showPsfSpatialCells(exposure, psfCellSet, showChi2=True,
                                            symb="o", ctype=ds9.YELLOW, ctypeBad=ds9.RED, size=8, frame=frame)
            if displayResiduals:
                maUtils.showPsfCandidates(exposure, psfCellSet, psf=psf, frame=4,
                                          normalize=normalizeResiduals,
                                          showBadCandidates=showBadCandidates)
            if displayPsfComponents:
                maUtils.showPsf(psf, frame=6)
            if displayPsfMosaic:
                maUtils.showPsfMosaic(exposure, psf, frame=7, showFwhm=True)
                ds9.ds9Cmd(ds9.selectFrame(frame=7) + " ;scale limits 0 1")
        #
        # Generate some QA information
        #
        # Count PSF stars
        #
        if metadata != None:
            metadata.set("spatialFitChi2", np.nan)
            metadata.set("numAvailStars", nCand)
            metadata.set("numGoodStars", numGoodStars)
            metadata.set("avgX", avgX)
            metadata.set("avgY", avgY)

        psfCellSet = None
        return psf, psfCellSet

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-    

psfDeterminerRegistry.register("psfex", PsfexPsfDeterminer)

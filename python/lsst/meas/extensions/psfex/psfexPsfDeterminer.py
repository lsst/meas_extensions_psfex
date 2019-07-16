#
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#
import os
import numpy as np

import lsst.daf.base as dafBase
import lsst.pex.config as pexConfig
import lsst.geom as geom
import lsst.afw.geom.ellipses as afwEll
import lsst.afw.display as afwDisplay
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
import lsst.meas.algorithms.utils as maUtils
import lsst.meas.extensions.psfex as psfex


class PsfexPsfDeterminerConfig(measAlg.BasePsfDeterminerConfig):
    __nEigenComponents = pexConfig.Field(
        doc="number of eigen components for PSF kernel creation",
        dtype=int,
        default=4,
    )
    spatialOrder = pexConfig.Field(
        doc="specify spatial order for PSF kernel creation",
        dtype=int,
        default=2,
        check=lambda x: x >= 0,
    )
    sizeCellX = pexConfig.Field(
        doc="size of cell used to determine PSF (pixels, column direction)",
        dtype=int,
        default=256,
        #        minValue = 10,
        check=lambda x: x >= 10,
    )
    sizeCellY = pexConfig.Field(
        doc="size of cell used to determine PSF (pixels, row direction)",
        dtype=int,
        default=sizeCellX.default,
        #        minValue = 10,
        check=lambda x: x >= 10,
    )
    __nStarPerCell = pexConfig.Field(
        doc="number of stars per psf cell for PSF kernel creation",
        dtype=int,
        default=3,
    )
    samplingSize = pexConfig.Field(
        doc="Resolution of the internal PSF model relative to the pixel size; "
        "e.g. 0.5 is equal to 2x oversampling",
        dtype=float,
        default=1,
    )
    badMaskBits = pexConfig.ListField(
        doc="List of mask bits which cause a source to be rejected as bad "
        "N.b. INTRP is used specially in PsfCandidateSet; it means \"Contaminated by neighbour\"",
        dtype=str,
        default=["INTRP", "SAT"],
    )
    psfexBasis = pexConfig.ChoiceField(
        doc="BASIS value given to psfex.  PIXEL_AUTO will use the requested samplingSize only if "
        "the FWHM < 3 pixels.  Otherwise, it will use samplingSize=1.  PIXEL will always use the "
        "requested samplingSize",
        dtype=str,
        allowed={
            "PIXEL": "Always use requested samplingSize",
            "PIXEL_AUTO": "Only use requested samplingSize when FWHM < 3",
        },
        default='PIXEL',
        optional=False,
    )
    __borderWidth = pexConfig.Field(
        doc="Number of pixels to ignore around the edge of PSF candidate postage stamps",
        dtype=int,
        default=0,
    )
    __nStarPerCellSpatialFit = pexConfig.Field(
        doc="number of stars per psf Cell for spatial fitting",
        dtype=int,
        default=5,
    )
    __constantWeight = pexConfig.Field(
        doc="Should each PSF candidate be given the same weight, independent of magnitude?",
        dtype=bool,
        default=True,
    )
    __nIterForPsf = pexConfig.Field(
        doc="number of iterations of PSF candidate star list",
        dtype=int,
        default=3,
    )
    tolerance = pexConfig.Field(
        doc="tolerance of spatial fitting",
        dtype=float,
        default=1e-2,
    )
    lam = pexConfig.Field(
        doc="floor for variance is lam*data",
        dtype=float,
        default=0.05,
    )
    reducedChi2ForPsfCandidates = pexConfig.Field(
        doc="for psf candidate evaluation",
        dtype=float,
        default=2.0,
    )
    spatialReject = pexConfig.Field(
        doc="Rejection threshold (stdev) for candidates based on spatial fit",
        dtype=float,
        default=3.0,
    )
    recentroid = pexConfig.Field(
        doc="Should PSFEX be permitted to recentroid PSF candidates?",
        dtype=bool,
        default=False,
    )

    def setDefaults(self):
        self.kernelSize = 41


class PsfexPsfDeterminerTask(measAlg.BasePsfDeterminerTask):
    ConfigClass = PsfexPsfDeterminerConfig

    def determinePsf(self, exposure, psfCandidateList, metadata=None, flagKey=None):
        """Determine a PSFEX PSF model for an exposure given a list of PSF
        candidates.

        Parameters
        ----------
        exposure: `lsst.afw.image.Exposure`
            Exposure containing the PSF candidates.
        psfCandidateList: iterable of `lsst.meas.algorithms.PsfCandidate`
            Sequence of PSF candidates typically obtained by detecting sources
            and then running them through a star selector.
        metadata: metadata, optional
            A home for interesting tidbits of information.
        flagKey: `lsst.afw.table.Key`, optional
            Schema key used to mark sources actually used in PSF determination.

        Returns
        -------
        psf: `lsst.meas.extensions.psfex.PsfexPsf`
            The determined PSF.
        """

        import lsstDebug
        display = lsstDebug.Info(__name__).display
        displayExposure = display and \
            lsstDebug.Info(__name__).displayExposure      # display the Exposure + spatialCells
        displayPsfComponents = display and \
            lsstDebug.Info(__name__).displayPsfComponents  # show the basis functions
        showBadCandidates = display and \
            lsstDebug.Info(__name__).showBadCandidates    # Include bad candidates (meaningless, methinks)
        displayResiduals = display and \
            lsstDebug.Info(__name__).displayResiduals     # show residuals
        displayPsfMosaic = display and \
            lsstDebug.Info(__name__).displayPsfMosaic     # show mosaic of reconstructed PSF(x,y)
        normalizeResiduals = lsstDebug.Info(__name__).normalizeResiduals
        afwDisplay.setDefaultMaskTransparency(75)
        # Normalise residuals by object amplitude

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
        else:
            psfCellSet = None

        sizes = np.empty(nCand)
        for i, psfCandidate in enumerate(psfCandidateList):
            try:
                if psfCellSet:
                    psfCellSet.insertCandidate(psfCandidate)
            except Exception as e:
                self.log.debug("Skipping PSF candidate %d of %d: %s", i, len(psfCandidateList), e)
                continue

            source = psfCandidate.getSource()
            quad = afwEll.Quadrupole(source.getIxx(), source.getIyy(), source.getIxy())
            rmsSize = quad.getTraceRadius()
            sizes[i] = rmsSize

        if self.config.kernelSize >= 15:
            self.log.warn("NOT scaling kernelSize by stellar quadrupole moment, but using absolute value")
            actualKernelSize = int(self.config.kernelSize)
        else:
            actualKernelSize = 2 * int(self.config.kernelSize * np.sqrt(np.median(sizes)) + 0.5) + 1
            if actualKernelSize < self.config.kernelSizeMin:
                actualKernelSize = self.config.kernelSizeMin
            if actualKernelSize > self.config.kernelSizeMax:
                actualKernelSize = self.config.kernelSizeMax
            if display:
                rms = np.median(sizes)
                print("Median PSF RMS size=%.2f pixels (\"FWHM\"=%.2f)" % (rms, 2*np.sqrt(2*np.log(2))*rms))

        # If we manually set the resolution then we need the size in pixel
        # units
        pixKernelSize = actualKernelSize
        if self.config.samplingSize > 0:
            pixKernelSize = int(actualKernelSize*self.config.samplingSize)
            if pixKernelSize % 2 == 0:
                pixKernelSize += 1
        self.log.trace("Psfex Kernel size=%.2f, Image Kernel Size=%.2f", actualKernelSize, pixKernelSize)
        psfCandidateList[0].setHeight(pixKernelSize)
        psfCandidateList[0].setWidth(pixKernelSize)

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- BEGIN PSFEX
        #
        # Insert the good candidates into the set
        #
        defaultsFile = os.path.join(os.environ["MEAS_EXTENSIONS_PSFEX_DIR"], "config", "default-lsst.psfex")
        args_md = dafBase.PropertySet()
        args_md.set("BASIS_TYPE", str(self.config.psfexBasis))
        args_md.set("PSFVAR_DEGREES", str(self.config.spatialOrder))
        args_md.set("PSF_SIZE", str(actualKernelSize))
        args_md.set("PSF_SAMPLING", str(self.config.samplingSize))
        prefs = psfex.Prefs(defaultsFile, args_md)
        prefs.setCommandLine([])
        prefs.addCatalog("psfexPsfDeterminer")

        prefs.use()
        principalComponentExclusionFlag = bool(bool(psfex.Context.REMOVEHIDDEN)
                                               if False else psfex.Context.KEEPHIDDEN)
        context = psfex.Context(prefs.getContextName(), prefs.getContextGroup(),
                                prefs.getGroupDeg(), principalComponentExclusionFlag)
        set = psfex.Set(context)
        set.setVigSize(pixKernelSize, pixKernelSize)
        set.setFwhm(2*np.sqrt(2*np.log(2))*np.median(sizes))
        set.setRecentroid(self.config.recentroid)

        catindex, ext = 0, 0
        backnoise2 = afwMath.makeStatistics(mi.getImage(), afwMath.VARIANCECLIP).getValue()
        ccd = exposure.getDetector()
        if ccd:
            gain = np.mean(np.array([a.getGain() for a in ccd]))
        else:
            gain = 1.0
            self.log.warn("Setting gain to %g" % (gain,))

        pc = 0
        contextvalp = []
        for i, key in enumerate(context.getName()):
            if context.getPcflag(i):
                raise RuntimeError("Principal Components can not be accessed")
                contextvalp.append(pcval[pc])  # noqa: F821
                pc += 1
            elif key[0] == ':':
                try:
                    contextvalp.append(exposure.getMetadata().getScalar(key[1:]))
                except KeyError:
                    raise RuntimeError("*Error*: %s parameter not found in the header of %s" %
                                       (key[1:], prefs.getContextName()))
            else:
                try:
                    contextvalp.append(np.array([psfCandidateList[_].getSource().get(key)
                                                 for _ in range(nCand)]))
                except KeyError:
                    raise RuntimeError("*Error*: %s parameter not found" % (key,))
                set.setContextname(i, key)

        if display:
            frame = 0
            if displayExposure:
                disp = afwDisplay.Display(frame=frame)
                disp.mtv(exposure, title="psf determination")

        badBits = mi.getMask().getPlaneBitMask(self.config.badMaskBits)
        fluxName = prefs.getPhotfluxRkey()
        fluxFlagName = "base_" + fluxName + "_flag"

        xpos, ypos = [], []
        for i, psfCandidate in enumerate(psfCandidateList):
            source = psfCandidate.getSource()
            xc, yc = source.getX(), source.getY()
            try:
                int(xc), int(yc)
            except ValueError:
                continue

            try:
                pstamp = psfCandidate.getMaskedImage().clone()
            except Exception:
                continue

            if fluxFlagName in source.schema and source.get(fluxFlagName):
                continue

            flux = source.get(fluxName)
            if flux < 0 or np.isnan(flux):
                continue

            # From this point, we're configuring the "sample" (PSFEx's version
            # of a PSF candidate).
            # Having created the sample, we must proceed to configure it, and
            # then fini (finalize), or it will be malformed.
            try:
                sample = set.newSample()
                sample.setCatindex(catindex)
                sample.setExtindex(ext)
                sample.setObjindex(i)

                imArray = pstamp.getImage().getArray()
                imArray[np.where(np.bitwise_and(pstamp.getMask().getArray(), badBits))] = \
                    -2*psfex.BIG
                sample.setVig(imArray)

                sample.setNorm(flux)
                sample.setBacknoise2(backnoise2)
                sample.setGain(gain)
                sample.setX(xc)
                sample.setY(yc)
                sample.setFluxrad(sizes[i])

                for j in range(set.getNcontext()):
                    sample.setContext(j, float(contextvalp[j][i]))
            except Exception as e:
                self.log.debug("Exception when processing sample at (%f,%f): %s", xc, yc, e)
                continue
            else:
                set.finiSample(sample)

            xpos.append(xc)  # for QA
            ypos.append(yc)

        if displayExposure:
            with disp.Buffering():
                disp.dot("o", xc, yc, ctype=afwDisplay.CYAN, size=4)

        if set.getNsample() == 0:
            raise RuntimeError("No good PSF candidates to pass to PSFEx")

        # ---- Update min and max and then the scaling
        for i in range(set.getNcontext()):
            cmin = contextvalp[i].min()
            cmax = contextvalp[i].max()
            set.setContextScale(i, cmax - cmin)
            set.setContextOffset(i, (cmin + cmax)/2.0)

        # Don't waste memory!
        set.trimMemory()

        # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- END PSFEX
        #
        # Do a PSFEX decomposition of those PSF candidates
        #
        fields = []
        field = psfex.Field("Unknown")
        field.addExt(exposure.getWcs(), exposure.getWidth(), exposure.getHeight(), set.getNsample())
        field.finalize()

        fields.append(field)

        sets = []
        sets.append(set)

        psfex.makeit(fields, sets)
        psfs = field.getPsfs()

        # Flag which objects were actually used in psfex by
        good_indices = []
        for i in range(sets[0].getNsample()):
            index = sets[0].getSample(i).getObjindex()
            if index > -1:
                good_indices.append(index)

        if flagKey is not None:
            for i, psfCandidate in enumerate(psfCandidateList):
                source = psfCandidate.getSource()
                if i in good_indices:
                    source.set(flagKey, True)

        xpos = np.array(xpos)
        ypos = np.array(ypos)
        numGoodStars = len(good_indices)
        avgX, avgY = np.mean(xpos), np.mean(ypos)

        psf = psfex.PsfexPsf(psfs[0], geom.Point2D(avgX, avgY))

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
                                            symb="o", ctype=afwDisplay.YELLOW, ctypeBad=afwDisplay.RED,
                                            size=8, display=disp)
            if displayResiduals:
                disp4 = afwDisplay.Display(frame=4)
                maUtils.showPsfCandidates(exposure, psfCellSet, psf=psf, display=disp4,
                                          normalize=normalizeResiduals,
                                          showBadCandidates=showBadCandidates)
            if displayPsfComponents:
                disp6 = afwDisplay.Display(frame=6)
                maUtils.showPsf(psf, display=disp6)
            if displayPsfMosaic:
                disp7 = afwDisplay.Display(frame=7)
                maUtils.showPsfMosaic(exposure, psf, display=disp7, showFwhm=True)
                disp.scale('linear', 0, 1)
        #
        # Generate some QA information
        #
        # Count PSF stars
        #
        if metadata is not None:
            metadata.set("spatialFitChi2", np.nan)
            metadata.set("numAvailStars", nCand)
            metadata.set("numGoodStars", numGoodStars)
            metadata.set("avgX", avgX)
            metadata.set("avgY", avgY)

        return psf, psfCellSet


measAlg.psfDeterminerRegistry.register("psfex", PsfexPsfDeterminerTask)

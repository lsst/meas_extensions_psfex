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

from lsst.pipe.base import Struct
import lsst.pex.config as pexConfig
import lsst.afw.display as afwDisplay
from lsst.meas.algorithms import BaseSourceSelectorTask, sourceSelectorRegistry
from . import psfexLib
from .psfex import compute_fwhmrange

__all__ = ["PsfexStarSelectorConfig", "PsfexStarSelectorTask"]


class PsfexStarSelectorConfig(BaseSourceSelectorTask.ConfigClass):
    fluxName = pexConfig.Field(
        dtype=str,
        doc="Name of photometric flux key ",
        default="base_PsfFlux",
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
        check=lambda x: x >= 0,
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
        check=lambda x: x >= 0.0,
    )
    minsn = pexConfig.Field(
        dtype=float,
        doc="Minimum S/N for candidates",
        default=100,
        check=lambda x: x >= 0.0,
    )

    def validate(self):
        pexConfig.Config.validate(self)

        if self.fluxErrName == "":
            self.fluxErrName = self.fluxName + ".err"
        elif self.fluxErrName != self.fluxName + ".err":
            msg = f"fluxErrName ({self.fluxErrName}) doesn't correspond to fluxName ({self.fluxName})"
            raise pexConfig.FieldValidationError(PsfexStarSelectorConfig.fluxName, self, msg)

        if self.minFwhm > self.maxFwhm:
            raise pexConfig.FieldValidationError(PsfexStarSelectorConfig.minFwhm, self,
                                                 f"minFwhm ({self.minFwhm}) > maxFwhm ({self.maxFwhm})")

    def setDefaults(self):
        self.badFlags = [
            "base_PixelFlags_flag_edge",
            "base_PixelFlags_flag_saturatedCenter",
            "base_PixelFlags_flag_crCenter",
            "base_PixelFlags_flag_bad",
            "base_PixelFlags_flag_suspectCenter",
            "base_PsfFlux_flag",
            # "parent",            # actually this is a test on deblend_nChild
        ]


class EventHandler():
    """A class to handle key strokes with matplotlib displays
    """

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

            which = np.where(dist == min(dist))

            x = self.x[which][0]
            y = self.y[which][0]
            for frame in self.frames:
                disp = afwDisplay.Display(frame=frame)
                disp.pan(x, y)
            disp.flush()
        else:
            pass


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

    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80))

    xmin = sorted(mag)[int(0.05*len(mag))]
    xmax = sorted(mag)[int(0.95*len(mag))]

    axes.set_xlim(-17.5, -13)
    axes.set_xlim(xmin - 0.1*(xmax - xmin), xmax + 0.1*(xmax - xmin))
    axes.set_ylim(0, 10)

    colors = ["r", "g", "b", "c", "m", "k", ]
    for k, mean in enumerate(centers):
        if k == 0:
            axes.plot(axes.get_xlim(), (mean, mean,), "k%s" % ltype)

        ll = (clusterId == k)
        axes.plot(mag[ll], width[ll], marker, markersize=markersize, markeredgewidth=markeredgewidth,
                  color=colors[k%len(colors)])

    ll = (clusterId == -1)
    axes.plot(mag[ll], width[ll], marker, markersize=markersize, markeredgewidth=markeredgewidth,
              color='k')

    if newFig:
        axes.set_xlabel("model")
        axes.set_ylabel(r"$\sqrt{I_{xx} + I_{yy}}$")

    return fig

## \addtogroup LSST_task_documentation
## \{
## \page PsfexStarSelectorTask
## \ref PsfexStarSelectorTask_ "PsfexStarSelectorTask"
## \copybrief PsfexStarSelectorTask
## \}


@pexConfig.registerConfigurable("psfex", sourceSelectorRegistry)
class PsfexStarSelectorTask(BaseSourceSelectorTask):
    r"""A star selector whose algorithm is not yet documented.

    @anchor PsfexStarSelectorTask_

    @section meas_extensions_psfex_psfexStarSelectorStarSelector_Contents  Contents

     - @ref meas_extensions_psfex_psfexStarSelectorStarSelector_Purpose
     - @ref meas_extensions_psfex_psfexStarSelectorStarSelector_Initialize
     - @ref meas_extensions_psfex_psfexStarSelectorStarSelector_IO
     - @ref meas_extensions_psfex_psfexStarSelectorStarSelector_Config
     - @ref meas_extensions_psfex_psfexStarSelectorStarSelector_Debug

    @section meas_extensions_psfex_psfexStarSelectorStarSelector_Purpose  Description

    A star selector whose algorithm is not yet documented

    @section meas_extensions_psfex_psfexStarSelectorStarSelector_Initialize  Task initialisation

    @copydoc \_\_init\_\_

    @section meas_extensions_psfex_psfexStarSelectorStarSelector_IO  Invoking the Task

    Like all star selectors, the main method is `run`.

    @section meas_extensions_psfex_psfexStarSelectorStarSelector_Config  Configuration parameters

    See @ref PsfexStarSelectorConfig

    @section meas_extensions_psfex_psfexStarSelectorStarSelector_Debug  Debug variables

    PsfexStarSelectorTask has a debug dictionary with the following keys:
    <dl>
    <dt>display
    <dd>bool; if True display debug information
    <dt>displayExposure
    <dd>bool; if True display the exposure and spatial cells
    <dt>plotFwhmHistogram
    <dd>bool; if True plot histogram of FWHM
    <dt>plotFlags
    <dd>bool: if True plot the sources coloured by their flags
    <dt>plotRejection
    <dd>bool; if True plot why sources are rejected
    </dl>

    For example, put something like:
    @code{.py}
        import lsstDebug
        def DebugInfo(name):
            di = lsstDebug.getInfo(name)  # N.b. lsstDebug.Info(name) would call us recursively
            if name.endswith("objectSizeStarSelector"):
                di.display = True
                di.displayExposure = True
                di.plotFwhmHistogram = True

            return di

        lsstDebug.Info = DebugInfo
    @endcode
    into your `debug.py` file and run your task with the `--debug` flag.
    """  # noqa: W505
    ConfigClass = PsfexStarSelectorConfig
    usesMatches = False  # selectStars does not use its matches argument

    def selectSources(self, sourceCat, matches=None, exposure=None):
        """Return a selection of psf-like objects.

        Parameters
        ----------
        sourceCat : `lsst.afw.table.SourceCatalog`
            Catalog of sources to select from.
            This catalog must be contiguous in memory.
        matches : `list` of `lsst.afw.table.ReferenceMatch` or None
            Ignored by this source selector.
        exposure : `lsst.afw.image.Exposure` or None
            The exposure the catalog was built from; used for debug display.

        Return
        ------
        struct : `lsst.pipe.base.Struct`
            The struct contains the following data:

            - selected : `numpy.ndarray` of `bool``
                Boolean array of sources that were selected, same length as
                sourceCat.
        """
        import lsstDebug
        display = lsstDebug.Info(__name__).display

        displayExposure = display and \
            lsstDebug.Info(__name__).displayExposure  # display the Exposure + spatialCells
        plotFwhmHistogram = display and plt and \
            lsstDebug.Info(__name__).plotFwhmHistogram  # Plot histogram of FWHM
        plotFlags = display and plt and \
            lsstDebug.Info(__name__).plotFlags  # Plot the sources coloured by their flags
        plotRejection = display and plt and \
            lsstDebug.Info(__name__).plotRejection  # Plot why sources are rejected
        afwDisplay.setDefaultMaskTransparency(75)

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
        shape = sourceCat.getShapeDefinition()
        ixx = sourceCat.get("%s.xx" % shape)
        iyy = sourceCat.get("%s.yy" % shape)

        fwhm = 2*np.sqrt(2*np.log(2))*np.sqrt(0.5*(ixx + iyy))
        elong = 0.5*(ixx - iyy)/(ixx + iyy)

        flux = sourceCat.get(fluxName)
        fluxErr = sourceCat.get(fluxErrName)
        sn = flux/np.where(fluxErr > 0, fluxErr, 1)
        sn[fluxErr <= 0] = -psfexLib.BIG

        flags = 0x0
        for i, f in enumerate(self.config.badFlags):
            flags = np.bitwise_or(flags, np.where(sourceCat.get(f), 1 << i, 0))
        #
        # Estimate the acceptable range of source widths
        #
        good = np.logical_and(sn > minsn, np.logical_not(flags))
        good = np.logical_and(good, elong < maxelong)
        good = np.logical_and(good, fwhm >= minFwhm)
        good = np.logical_and(good, fwhm < maxFwhm)

        fwhmMode, fwhmMin, fwhmMax = compute_fwhmrange(fwhm[good], maxFwhmVariability, minFwhm, maxFwhm,
                                                       plot=dict(fwhmHistogram=plotFwhmHistogram))

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        #
        # Here's select_candidates
        #
        # ---- Apply some selection over flags, fluxes...

        bad = (flags != 0)
        # set.setBadFlags(int(sum(bad)))

        if plotRejection:
            selectionVectors = []
            selectionVectors.append((bad, "flags %d" % sum(bad)))

        dbad = sn < minsn
        # set.setBadSN(int(sum(dbad)))
        bad = np.logical_or(bad, dbad)
        if plotRejection:
            selectionVectors.append((dbad, "S/N %d" % sum(dbad)))

        dbad = fwhm < fwhmMin
        # set.setBadFrmin(int(sum(dbad)))
        bad = np.logical_or(bad, dbad)
        if plotRejection:
            selectionVectors.append((dbad, "fwhmMin %d" % sum(dbad)))

        dbad = fwhm > fwhmMax
        # set.setBadFrmax(int(sum(dbad)))
        bad = np.logical_or(bad, dbad)
        if plotRejection:
            selectionVectors.append((dbad, "fwhmMax %d" % sum(dbad)))

        dbad = elong > maxelong
        # set.setBadElong(int(sum(dbad)))
        bad = np.logical_or(bad, dbad)
        if plotRejection:
            selectionVectors.append((dbad, "elong %d" % sum(dbad)))

        # -- ... and check the integrity of the sample
        if maxbadflag:
            raise RuntimeError("vignet variable not defined. Code is broken.")
            nbad = np.array([(v <= -psfexLib.BIG).sum() for v in vignet])  # noqa: F821
            dbad = nbad > maxbad
            # set.setBadPix(int(sum(dbad)))
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
            disp = afwDisplay.Display(frame=frame)
            disp.mtv(mi, title="PSF candidates")

            with disp.Buffering():
                for i, source in enumerate(sourceCat):
                    if good[i]:
                        ctype = afwDisplay.GREEN  # star candidate
                    else:
                        ctype = afwDisplay.RED  # not star

                    disp.dot("+", source.getX() - mi.getX0(), source.getY() - mi.getY0(),
                             ctype=ctype)

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
                            label = re.sub(r"\_flag", "",
                                           re.sub(r"^base\_", "",
                                                  re.sub(r"^.*base\_PixelFlags\_flag\_", "", f)))
                            plt.plot(imag[isSet], fwhm[isSet], 'o', alpha=alpha, label=label)
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
            eventHandler = EventHandler(plt.axes(), imag, fwhm, sourceCat.getX(), sourceCat.getY(),
                                        frames=[frame])

        if plotFlags or plotRejection:
            while True:
                try:
                    reply = input("continue? [y[es] h(elp) p(db) q(uit)] ").strip()
                except EOFError:
                    reply = "y"

                if not reply:
                    reply = "y"

                if reply[0] == "h":
                    print("""\
At this prompt, you can continue with almost any key; 'p' enters pdb,
                                                      'q' returns to the shell, and
                                                      'h' prints this text
""", end=' ')

                    if displayExposure:
                        print("""
If you put the cursor on a point in the matplotlib scatter plot and hit 'p' you'll see it in ds9.""")
                elif reply[0] == "p":
                    import pdb
                    pdb.set_trace()
                elif reply[0] == 'q':
                    sys.exit(1)
                else:
                    break

        return Struct(selected=good)

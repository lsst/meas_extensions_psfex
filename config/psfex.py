import sys
import lsst.meas.extensions.psfex.psfexStarSelector
import lsst.meas.extensions.psfex.psfexPsfDeterminer

root.calibrate.measurePsf.starSelector.name  = "psfex" if False else "objectSize"
root.calibrate.measurePsf.psfDeterminer.name = "psfex" if True  else "pca"

root.calibrate.measurePsf.starSelector["objectSize"].sourceFluxField = "initial.flux.psf"
try:
    root.calibrate.measurePsf.starSelector["objectSize"].nSigmaClip = 1.75
except AttributeError, e:
    print >> sys.stderr, "Warning: %s" % (e,) # needs a post-June-2013 meas_algorithms

root.calibrate.measurePsf.starSelector["psfex"].maxFwhmVariability = 0.1
root.calibrate.measurePsf.starSelector["psfex"].maxbadflag = False

root.calibrate.measurePsf.psfDeterminer["psfex"].spatialOrder = 1

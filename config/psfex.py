import sys
import lsst.meas.extensions.psfex.psfexStarSelector
import lsst.meas.extensions.psfex.psfexPsfDeterminer

config.calibrate.measurePsf.starSelector.name  = "psfex" if False else "objectSize"
config.calibrate.measurePsf.psfDeterminer.name = "psfex" if True  else "pca"

config.calibrate.measurePsf.starSelector["objectSize"].sourceFluxField = "initial.flux.psf"
try:
    config.calibrate.measurePsf.starSelector["objectSize"].nSigmaClip = 1.75
except AttributeError, e:
    print >> sys.stderr, "Warning: %s" % (e,) # needs a post-June-2013 meas_algorithms

config.calibrate.measurePsf.starSelector["psfex"].maxFwhmVariability = 0.1
config.calibrate.measurePsf.starSelector["psfex"].maxbadflag = False

config.calibrate.measurePsf.psfDeterminer["psfex"].spatialOrder = 1 if False else 2

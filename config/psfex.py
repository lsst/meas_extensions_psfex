import sys
from lsst.meas.extensions.psfex import PsfExStarSelectorTask

config.calibrate.measurePsf.starSelector.retarget(PsfExStarSelectorTask)
config.calibrate.measurePsf.psfDeterminer.name = "psfex" if True else "pca"

config.calibrate.measurePsf.starSelector.sourceFluxField = "initial.flux.psf"
try:
    config.calibrate.measurePsf.starSelector.nSigmaClip = 1.75
except AttributeError, e:
    print >> sys.stderr, "Warning: %s" % (e,)  # needs a post-June-2013 meas_algorithms

config.calibrate.measurePsf.starSelector.maxFwhmVariability = 0.1
config.calibrate.measurePsf.starSelector.maxbadflag = False

config.calibrate.measurePsf.psfDeterminer["psfex"].spatialOrder = 1 if False else 2

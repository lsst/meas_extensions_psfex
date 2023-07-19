import lsst.meas.algorithms.makePsfCandidates
assert type(config)==lsst.meas.algorithms.makePsfCandidates.MakePsfCandidatesConfig, 'config is of type %s.%s instead of lsst.meas.algorithms.makePsfCandidates.MakePsfCandidatesConfig' % (type(config).__module__, type(config).__name__)
# Size of the postage stamp in pixels (excluding the border) around each star that is extracted for fitting. Should be odd and preferably at least 25.
config.kernelSize=25

# Number of pixels to ignore around the edge of PSF candidate postage stamps
config.borderWidth=0


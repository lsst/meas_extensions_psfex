import lsst.meas.extensions.psfex.psfexPsfDeterminer
assert type(config)==lsst.meas.extensions.psfex.psfexPsfDeterminer.PsfexPsfDeterminerConfig, 'config is of type %s.%s instead of lsst.meas.extensions.psfex.psfexPsfDeterminer.PsfexPsfDeterminerConfig' % (type(config).__module__, type(config).__name__)
# Size of the postage stamp (in native pixels) to render the PSF model. Should be odd.
config.stampSize=41

# specify spatial order for PSF kernel creation
config.spatialOrder=2

# size of cell used to determine PSF (pixels, column direction)
config.sizeCellX=256

# size of cell used to determine PSF (pixels, row direction)
config.sizeCellY=256

# Resolution of the internal PSF model relative to the pixel size; e.g. 0.5 is equal to 2x oversampling
config.samplingSize=0.5

# List of mask bits which cause a source to be rejected as bad N.b. INTRP is used specially in PsfCandidateSet; it means "Contaminated by neighbour"
config.badMaskBits=['INTRP', 'SAT']

# BASIS value given to psfex.  PIXEL_AUTO will use the requested samplingSize only if the FWHM < 3 pixels.  Otherwise, it will use samplingSize=1.  PIXEL will always use the requested samplingSize
config.psfexBasis='PIXEL_AUTO'

# tolerance of spatial fitting
config.tolerance=0.01

# floor for variance is lam*data
config.lam=0.05

# for psf candidate evaluation
config.reducedChi2ForPsfCandidates=2.0

# Rejection threshold (stdev) for candidates based on spatial fit
config.spatialReject=3.0

# Should PSFEX be permitted to recentroid PSF candidates?
config.recentroid=False


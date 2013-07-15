import re
import sys
import numpy as np
import pyfits
import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.afw.display.ds9 as ds9
import lsst.meas.extensions.psfex as psfex

def readSExtractor(filename):
    with pyfits.open(filename) as pf:
        for hdu in pf:
            if hdu.name == "PRIMARY":
                pass
            elif hdu.name == "LDAC_IMHEAD":
                hdr = hdu.data[0][0]    # the fits header from the original fits image
                print hdr[3]
            elif hdu.name == "LDAC_OBJECTS":
                print "%d objects" % (len(hdu.data))
                # Find the VIGNET column
                ttype = [k for k, v in hdu.header.items() if v == "VIGNET"]
                if not ttype:
                    raise RuntimeError("Unable to find a VIGNET column")
                vignetCol = int(re.search(r"^TTYPE(\d+)$", ttype[0]).group(1)) - 1

                for row in range(len(hdu.data)):
                    pixelData = hdu.data[row][vignetCol]
                    bad = np.where(pixelData < -1e29)
                    sat = np.where(pixelData > 99e3)
                    pixelData[bad] = 0.0
                    mi = afwImage.MaskedImageF(*hdu.data[row][vignetCol].shape)
                    im = mi.getImage()
                    im.getArray()[:] = pixelData
                    msk = mi.getMask().getArray()
                    msk[bad] = afwImage.MaskU.getPlaneBitMask("BAD")
                    msk[sat] = afwImage.MaskU.getPlaneBitMask("SAT")
                    ds9.mtv(mi, title=row)
                    raw_input("Next ")
            
def readPrefs(filename, md=None):
    return psfex.Prefs(filename, md)

def foo():
    pp = psfex.Psf(psfex.Context(["A", "B"], [1, 1], [2], 1, True), [1], 0.1, [1], 1)

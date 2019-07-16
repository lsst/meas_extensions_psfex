// -*- lsst-C++ -*-
#include <cstring>
#include "lsst/meas/extensions/psfex/Field.hh"
#undef PI
#include "lsst/geom/Point.h"
#include "lsst/afw/geom/SkyWcs.h"

extern "C" {
    #include "globals.h"
}

namespace lsst { namespace meas { namespace extensions { namespace psfex {

Field::Field(std::string const& ident) :
    impl(NULL), _isInitialized(false)
{
    QCALLOC(impl, fieldstruct, 1);
    impl->next = 0;

    strcpy(impl->catname, ident.c_str());
    impl->rcatname = impl->catname;
#if 0
    strncpy(impl->rtcatname, impl->rcatname, sizeof(impl->rtcatname) - 1);
    strncpy(impl->ident, "??", sizeof(impl->ident) - 1);
#elif 1
    if (!(impl->rcatname = strrchr(impl->catname, '/'))) {
        impl->rcatname = impl->catname;
    } else {
        ++impl->rcatname;
    }

    strncpy(impl->rtcatname, impl->rcatname, sizeof(impl->rtcatname) - 1);
    {
        char *pstr=strrchr(impl->rtcatname, '.');
        if (pstr) {
            *pstr = '\0';
        }
    }

    strncpy(impl->ident, "??", sizeof(impl->ident) - 1);
#endif

    impl->ndet = 0;
    impl->psf = NULL;
    impl->wcs = NULL;

    _finalize();
}

Field::~Field()
{
    for (int i = 0; i != impl->next; ++i) {
        free(impl->wcs[i]);             // psfex's wcs isn't quite the same as ours ...
        impl->wcs[i] = NULL;            // ... so don't let psfex free it
    }
    field_end(impl);
    impl = NULL;
}

/************************************************************************************************************/

void
Field::_finalize(bool force)
{
    if (force || !_isInitialized) {
        field_init_finalize(impl);
        _isInitialized = true;
    }
}

/************************************************************************************************************/

std::vector<Psf>
Field::getPsfs() const
{
    if (_psfs.empty()) {
        _psfs.reserve(impl->next);
        for (int i = 0; i != impl->next; ++i) {
            _psfs.push_back(Psf(impl->psf[i]));
        }
    }

    return _psfs;
}

void
Field::addExt(lsst::afw::geom::SkyWcs const& wcs_,
              int const naxis1, int const naxis2,
              int const nobj)
{
    QREALLOC(impl->psf, psfstruct *, impl->next + 1);
    impl->psf[impl->next] = 0;
    QREALLOC(impl->wcs, wcsstruct *, impl->next + 1);
    impl->wcs[impl->next] = 0;
    /*
     * We're going to fake psfex's wcsstruct object.  We only need enough of it for field_locate
     */
    QMALLOC(impl->wcs[impl->next], wcsstruct, 1);
    wcsstruct *wcs = impl->wcs[impl->next];

    wcs->naxis = 2;
    wcs->naxisn[0] = naxis1;
    wcs->naxisn[1] = naxis2;

    auto const crval = wcs_.getSkyOrigin();
    // crpix using the FITS standard = pixel origin using the LSST standard + 1
    auto const crpix = wcs_.getPixelOrigin() + geom::Extent2D(1, 1);
    auto const cdMatrix = wcs_.getCdMatrix();
    std::string const cunit("DEG");
    auto metadata = wcs_.getFitsMetadata();
    for (int i = 0; i != wcs->naxis; ++i) {
        auto ifits = i + 1;
        auto ctype = metadata->getAsString("CTYPE" + std::to_string(ifits));
        strncpy(wcs->ctype[i], ctype.c_str(), ctype.size() + 1);
        strncpy(wcs->cunit[i], cunit.c_str(), cunit.size() + 1);
        wcs->crpix[i] = crpix[i];
        wcs->crval[i] = crval[i].asDegrees();
        wcs->cdelt[i] = 1.0;  // scale is in the CD matrix (is this even needed?)
        wcs->crder[i] = 0;
        wcs->csyer[i] = 0;
    }
    for (int i = 0, k = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j, ++k) {
            wcs->cd[k] = cdMatrix(i, j);
        }
    }
    wcs->lng = 0;
    wcs->lat = 1;
    wcs->equinox = 2000;

    auto center = wcs_.pixelToSky(geom::Point2D(0.5*naxis1, 0.5*naxis2));
    wcs->wcsscalepos[0] = center.getLongitude().asDegrees();
    wcs->wcsscalepos[1] = center.getLatitude().asDegrees();

    double maxradius = 0.0;             // Maximum distance to wcsscalepos
    for (int x = 0; x <= 1; ++x) {
        for (int y = 0; y <= 1; ++y) {
            geom::Point2D point(x*naxis1, y*naxis2); // Corner
            double const radius = center.separation(wcs_.pixelToSky(point)).asDegrees();
            if (radius > maxradius) {
                maxradius = radius;
            }
        }
    }
    wcs->wcsmaxradius = maxradius;

    impl->ndet += nobj;

    ++impl->next;
}

}}}}

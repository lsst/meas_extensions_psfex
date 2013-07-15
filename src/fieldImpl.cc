// -*- lsst-C++ -*-
#include <cstring>
#include "lsst/meas/extensions/psfex/Field.hh"
#include "wcslib/wcs.h"
#undef PI
#include "lsst/afw/image/Wcs.h"

extern "C" {
    #include "globals.h"
}

namespace astromatic { namespace psfex {

Field::Field(std::string const& ident) :
    impl(new fieldstruct, field_end), _isInitialized(false)
{
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
    //field_end(impl);
}
        
/************************************************************************************************************/

void
Field::_finalize(bool force)
{
    if (force || !_isInitialized) {
        field_init_finalize(impl.get());
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
            _psfs.push_back(Psf(impl->psf[i], impl));
        }
    }

    return _psfs;
}

/************************************************************************************************************/
//
// This class exists solely so that I can recover the protected data member _wcsInfo
//
namespace {
struct PsfUnpack : private lsst::afw::image::Wcs {
    PsfUnpack(lsst::afw::image::Wcs const& wcs) : Wcs(wcs) { }
    const struct wcsprm* getWcsInfo() { return _wcsInfo; }
};
}

void
Field::addExt(lsst::afw::image::Wcs const& wcs_,
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
    PsfUnpack wcsUnpacked(wcs_);
    struct wcsprm const* wcsPrm = wcsUnpacked.getWcsInfo();
    QMALLOC(impl->wcs[impl->next], wcsstruct, 1);
    wcsstruct *wcs = impl->wcs[impl->next];
    
    wcs->naxis = wcsPrm->naxis;
    wcs->naxisn[0] = naxis1;
    wcs->naxisn[1] = naxis2;
    
    for (int i = 0; i != wcs->naxis; ++i) {
        strncpy(wcs->ctype[i], wcsPrm->ctype[i], sizeof(wcs->ctype[i]) - 1);
        strncpy(wcs->cunit[i], wcsPrm->cunit[i], sizeof(wcs->cunit[i]) - 1);
        wcs->crval[i] = wcsPrm->crval[i];
        
        wcs->cdelt[i] = wcsPrm->cdelt[i];
        wcs->crpix[i] = wcsPrm->crpix[i];
        wcs->crder[i] = wcsPrm->crder[i];
        wcs->csyer[i] = wcsPrm->csyer[i];
        wcs->crval[i] = wcsPrm->crval[i];
    }
    for (int i = 0; i != wcs->naxis*wcs->naxis; ++i) {
        wcs->cd[i] = wcsPrm->cd[i];
    }
    wcs->longpole = wcsPrm->lonpole;
    wcs->latpole = wcsPrm->latpole;
    wcs->lat = wcsPrm->lat;
    wcs->lng = wcsPrm->lng;
    wcs->equinox = wcsPrm->equinox;

    impl->ndet += nobj;
    
    ++impl->next;
}

}}

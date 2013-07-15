// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
/*!
 * @brief Represent a PSF using the representation from Emmanuel's PSFEX code
 *
 * @file
 */
#include <cmath>
#include <cassert>

#include "boost/make_shared.hpp"

#include "lsst/base.h"
#include "lsst/utils/ieee.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/math/Statistics.h"
extern "C" {
#include "vignet.h"
}
#include "lsst/meas/extensions/psfex/PsfexPsf.h"
#include "lsst/meas/algorithms/KernelPsfFactory.h"

namespace astromatic {
namespace psfex {

namespace afw = lsst::afw;

PsfexPsf::PsfexPsf(
    astromatic::psfex::Psf const& psf,
    afw::geom::Point2D const & averagePosition
                  ) : ImagePsf(), _averagePosition(averagePosition),
                      _size(psf.impl->dim),
                      _comp(psf.impl->npix),
                      _context(psf.impl->poly->ndim)
                      
{
    _poly = poly_copy(psf.impl->poly);

    _pixstep = psf.impl->pixstep;

    std::copy(psf.impl->size, psf.impl->size + psf.impl->dim, _size.begin());

    std::copy(psf.impl->comp, psf.impl->comp + psf.impl->npix, _comp.begin());

    for (int i = 0; i != psf.impl->poly->ndim; ++i) {
        _context[i].first = psf.impl->contextoffset[i];
        _context[i].second = psf.impl->contextscale[i];
    }
}

    PsfexPsf::PsfexPsf() : ImagePsf(),
                           _averagePosition(afw::geom::Point2I(0, 0)),
                           _poly(0),
                           _pixstep(0.0),
                           _size(),
                           _comp(),
                           _context()
{
    ;
}

PsfexPsf::~PsfexPsf()
{
    poly_end(_poly);
}

PTR(afw::detection::Psf)
PsfexPsf::clone() const {
    return boost::make_shared<PsfexPsf>(*this);
}

PTR(afw::math::LinearCombinationKernel const)
PsfexPsf::getKernel(afw::geom::Point2D position) const
{
    double pos[MAXCONTEXT];
    int const ndim = _context.size();
    if (ndim != 2) {                    // we're only handling spatial variation for now
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          str(boost::format("Only spatial variation (ndim == 2) is supported; saw %d")
                              % ndim));

    }
    // where we want to evaluate the basis function's weights
    if (!lsst::utils::isfinite(position[0])) {
        position = _averagePosition;
    }

    for (int i = 0; i < ndim; ++i) {
        pos[i] = (position[i] - _context[i].first)/_context[i].second;
    }

    poly_func(_poly, pos);              // evaluate polynomial

    int const w = _size[0], h = _size[1];
    std::vector<float> fullresIm(w*h);       // accumulate full-resolution image into this buffer
    /*
     * Create a fixed Kernel out of each component, and then create a LinearCombinationKernel from them
     */
    const int nbasis = _size.size() > 2 ? _size[2] : 1; // number of basis functions
    afw::math::KernelList kernels; kernels.reserve(nbasis); // the insides of the LinearCombinationKernel
    std::vector<double>   weights; weights.reserve(nbasis);

    float const vigstep = 1/_pixstep;
    float const dx = 0.0, dy = 0.0;
    std::vector<float> sampledBasis(w*h);

    afw::detection::Psf::Image kim(w, h); // a basis function image, to be copied into a FixedKernel
    kim.setXY0(-w/2, -h/2);

    for (int i = 0; i != nbasis; ++i) {
        /*
         * Resample the basis function onto the output resolution (and potentially subpixel offset)
         */
        vignet_resample(const_cast<float *>(&_comp[i*w*h]), w, h,
                        &sampledBasis[0],                   w, h,
                        -dy*vigstep, -dx*vigstep, vigstep, 1.0); // n.b. x and y are transposed
        //
        // And copy it into place
        //
        {
            float *pl = &sampledBasis[0];
            for (int y = 0; y != h; ++y) {
                for (int x = 0; x != w; ++x) {
                    kim(y, x) = *pl++;  // N.b.: (y, x) --- we're transposing the data
                }
            }
        }

        kernels.push_back(boost::make_shared<afw::math::FixedKernel>(kim));
        weights.push_back(_poly->basis[i]);
    }

    _kernel = boost::make_shared<afw::math::LinearCombinationKernel>(kernels, weights);

    return _kernel;
}

PTR(afw::detection::Psf::Image)
PsfexPsf::doComputeImage(afw::geom::Point2D const & position,
                         afw::image::Color const & color) const {
    return _doComputeImage(position, color, position);
}
    
PTR(afw::detection::Psf::Image)
PsfexPsf::doComputeKernelImage(afw::geom::Point2D const& position,
                               afw::image::Color const& color) const
{
    return _doComputeImage(position, color, afw::geom::Point2D(0, 0));
}

PTR(afw::detection::Psf::Image)
PsfexPsf::_doComputeImage(afw::geom::Point2D const& position,
                          afw::image::Color const&,
                          afw::geom::Point2D const& center
        ) const
{
    double pos[MAXCONTEXT];
    int const ndim = _context.size();
    if (ndim != 2) {                    // we're only handling spatial variation for now
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          str(boost::format("Only spatial variation (ndim == 2) is supported; saw %d")
                              % ndim));

    }

    for (int i = 0; i < ndim; ++i) {
        pos[i] = (position[i] - _context[i].first)/_context[i].second;
    }

    poly_func(_poly, pos);              // evaluate polynomial

    int const w = _size[0], h = _size[1];
    std::vector<float> fullresIm(w*h);       // accumulate full-resolution image into this buffer
    const int nbasis = _size.size() > 2 ? _size[2] : 1; // number of basis functions

    /* Sum each component */
    int const npix = w*h;
    for (int i = 0; i != nbasis; ++i) {
        float *pl = &fullresIm[0];
        float const fac = _poly->basis[i];
        float const *ppc = &_comp[i*w*h];

        for (int j = 0; j != npix; ++j) {
            pl[j] +=  fac*ppc[j];
        }
    }
    /*
     * We now have the image reconstructed at internal resolution; resample it onto the output resolution
     * and subpixel offset
     */
    float const vigstep = 1/_pixstep;
    float dx = center[0] - static_cast<int>(center[0]);
    float dy = center[1] - static_cast<int>(center[1]);
    if (dx > 0.5) dx -= 1.0;
    if (dy > 0.5) dy -= 1.0;
    
    std::vector<float> sampledIm(w*h);
    vignet_resample(&fullresIm[0], w, h,
                    &sampledIm[0], w, h,
                    -dy*vigstep, -dx*vigstep, vigstep, 1.0); // n.b. x and y are transposed
    //
    // And copy it into place
    //
    PTR(afw::detection::Psf::Image) im = boost::make_shared<afw::detection::Psf::Image>(w, h);
    // N.b. center[0] - dx == (int)center[] if we hadn't reduced dx to (-0.5, 0.5].
    // The + 0.5 is to handle floating point imprecision in this calculation
    im->setXY0(static_cast<int>(center[0] - dx - w/2 + 0.5),
               static_cast<int>(center[1] - dy - h/2 + 0.5));
    {
        float *pl = &sampledIm[0];
        for (int y = 0; y != h; ++y) {
            for (int x = 0; x != w; ++x) {
                (*im)(y, x) = *pl++;       // N.b.: (y, x) --- we're transposing the data
            }
        }
    }

    return im;
}

/************************************************************************************************************/
/*
 * All the rest of this file handles persistence to FITS files
 */
namespace table = afw::table;

namespace {

class PsfexPsfSchema1 {
public:
    PsfexPsfSchema1() :
        schema(),
        ndim(schema.addField<int>("ndim", "Number of elements in group")),
        ngroup(schema.addField<int>("ngroup", "Number of elements in degree")),
        ncoeff(schema.addField<int>("ncoeff", "Number of coefficients")),

        _size_size(schema.addField<int>("_size_size", "Size of _size array")),
        _comp_size(schema.addField<int>("_comp_size", "Size of _comp array")),
        _context_size(schema.addField<int>("_context_size", "Size of _context array")),
    // Other scalars
        averagePosition(schema.addField< afw::table::Point<double> >(
                        "averagePosition", "average position of stars used to make the PSF")),
        _pixstep(schema.addField<float>("_pixstep", "oversampling", "pixel"))
    {
        ;
    }

    table::Schema schema;

    // Sizes in _poly
    table::Key<int> ndim;
    table::Key<int> ngroup;
    table::Key<int> ncoeff;
    // Sizes of vectors
    table::Key<int> _size_size;
    table::Key<int> _comp_size;
    table::Key<int> _context_size;
    // Other scalars
    table::Key<table::Point<double> > averagePosition;
    table::Key<float> _pixstep;
};

class PsfexPsfSchema2 {
public:
    PsfexPsfSchema2(int const ndim, int const ngroup, int const ncoeff,
                    int size_size, int comp_size, int context_size) :
                    
        schema(),
        group(schema.addField<table::Array<int> >("group", "Groups (of coefficients?)", ndim)),
        degree(schema.addField<table::Array<int> >("degree", "Degree in each group", ngroup)),
        basis(schema.addField<table::Array<double> >("basis", "Values of the basis functions", ncoeff)),
        coeff(schema.addField<table::Array<double> >("coeff", "Polynomial coefficients", ncoeff)),
        _size(schema.addField<table::Array<int> >("_size", "PSF dimensions", size_size)),
        _comp(schema.addField<table::Array<float> >("_comp", "Complete pixel data", comp_size)),
        _context_first( schema.addField<table::Array<double> >("_context_first",
                                                   "Offset to apply to context data", context_size)),
        _context_second(schema.addField<table::Array<double> >("_context_second",
                                                   "Scale to apply to context data", context_size))
    {
        ;
    }

    table::Schema schema;
    // _poly
    table::Key<table::Array<int> > group;              // len(group) == ndim
    table::Key<table::Array<int> > degree;             // len(degree) == ngroup
    table::Key<table::Array<double> > basis;           // len(basis) == ncoeff
    table::Key<table::Array<double> > coeff;           // len(coeff) == ncoeff
    // vectors
    table::Key<table::Array<int> > _size;
    table::Key<table::Array<float> > _comp;
    table::Key<table::Array<double> > _context_first;
    table::Key<table::Array<double> > _context_second;
};

std::string getPsfexPsfPersistenceName() { return "PsfexPsf"; }

} // anonymous

/************************************************************************************************************/

namespace detail {                      // PsfexPsfFactory needs to be a friend of PsfexPsf
class PsfexPsfFactory : public table::io::PersistableFactory {
public:

virtual PTR(table::io::Persistable)
read(InputArchive const & archive, CatalogVector const & catalogs) const {
    LSST_ARCHIVE_ASSERT(catalogs.size() == 2u);

    PTR(PsfexPsf) result(new PsfexPsf());

    int ndim, ngroup, ncoeff;
    int size_size, comp_size, context_size;
    {
        PsfexPsfSchema1 const & keys = PsfexPsfSchema1();
        LSST_ARCHIVE_ASSERT(catalogs[0].size() == 1u);
        LSST_ARCHIVE_ASSERT(catalogs[0].getSchema() == keys.schema);
        table::BaseRecord const & record = catalogs[0].front();
        
        // fields in _poly
        ndim = record.get(keys.ndim);
        ngroup = record.get(keys.ngroup);
        ncoeff = record.get(keys.ncoeff);
        // Other scalars
        result->_averagePosition = record.get(keys.averagePosition);
        result->_pixstep = record.get(keys._pixstep);
        // sizes of vectors
        size_size = record.get(keys._size_size);
        comp_size = record.get(keys._comp_size);
        context_size = record.get(keys._context_size);
    }
    // Now we can read the data
    {
        PsfexPsfSchema2 const keys(ndim, ngroup, ncoeff,
                                   size_size, comp_size, context_size);

        LSST_ARCHIVE_ASSERT(catalogs[1].size() == 1u);
        LSST_ARCHIVE_ASSERT(catalogs[1].getSchema() == keys.schema);
        table::BaseRecord const & record = catalogs[1].front();

        // _poly
        std::vector<int> group;
        {
            int const *begin = record.getElement(keys.group);
            group.assign(begin, begin + ndim);

            for (int i = 0; i != ndim; ++i) {
                ++group[i];             // poly_init subtracts 1 from each element.  Sigh.
            }
        }
        std::vector<int> degree;
        {
            int const *begin = record.getElement(keys.degree);
            degree.assign(begin, begin + ngroup);
        }
        result->_poly = poly_init(&group[0], group.size(), &degree[0], degree.size());
        LSST_ARCHIVE_ASSERT(result->_poly->ncoeff == ncoeff);

        {
            double const *begin = record.getElement(keys.basis);
            std::copy(begin, begin + ncoeff, result->_poly->basis);
        }
        {
            double const *begin = record.getElement(keys.coeff);
            std::copy(begin, begin + ncoeff, result->_poly->coeff);
        }
        // vectors
        {
            int const *begin = record.getElement(keys._size);
            result->_size.assign(begin, begin + size_size);
        }
        {
            float const *begin = record.getElement(keys._comp);
            result->_comp.assign(begin, begin + comp_size);
        }
        {
            double const *begin1 = record.getElement(keys._context_first);
            double const *begin2 = record.getElement(keys._context_second);
            result->_context.resize(context_size);
            for (int i = 0; i != context_size; ++i) {
                result->_context[i].first  = begin1[i];
                result->_context[i].second = begin2[i];
            }
        }
    }

    return result;
}

explicit PsfexPsfFactory(std::string const & name) : table::io::PersistableFactory(name) {}

};
}

namespace {
    detail::PsfexPsfFactory registration(getPsfexPsfPersistenceName());
}

/************************************************************************************************************/

std::string PsfexPsf::getPythonModule() const { return "lsst.meas.extensions.psfex"; }

std::string PsfexPsf::getPersistenceName() const { return getPsfexPsfPersistenceName(); }

void PsfexPsf::write(afw::table::io::OutputArchiveHandle & handle) const {
    // First write the dimensions of the various arrays to an HDU, so we can construct the second
    // HDU using them when we read the Psf back
    {
        PsfexPsfSchema1 const keys;
        afw::table::BaseCatalog cat = handle.makeCatalog(keys.schema);
        PTR(afw::table::BaseRecord) record = cat.addNew();
        
        // Sizes in _poly
        record->set(keys.ndim, _poly->ndim);
        record->set(keys.ngroup, _poly->ngroup);
        record->set(keys.ncoeff, _poly->ncoeff);
        // Other scalars
        record->set(keys.averagePosition, _averagePosition);
        record->set(keys._pixstep, _pixstep);
        record->set(keys._size_size, _size.size());
        record->set(keys._comp_size, _comp.size());
        record->set(keys._context_size, _context.size());

        handle.saveCatalog(cat);
    }
    // Now we can write the data
    {
        PsfexPsfSchema2 const keys(_poly->ndim, _poly->ngroup, _poly->ncoeff,
                                   _size.size(), _comp.size(), _context.size());
        afw::table::BaseCatalog cat = handle.makeCatalog(keys.schema);
        // _poly
        PTR(afw::table::BaseRecord) record = cat.addNew();
        {
            int *begin = record->getElement(keys.group);
            std::copy(_poly->group, _poly->group + _poly->ndim, begin);
        }
        {
            int *begin = record->getElement(keys.degree);
            std::copy(_poly->degree, _poly->degree + _poly->ngroup, begin);
        }
        {
            double *begin = record->getElement(keys.basis);
            std::copy(_poly->basis, _poly->basis + _poly->ncoeff, begin);
        }
        {
            double *begin = record->getElement(keys.coeff);
            std::copy(_poly->coeff, _poly->coeff + _poly->ncoeff, begin);
        }
        // vectors
        {
            int *begin = record->getElement(keys._size);
            std::copy(_size.begin(), _size.end(), begin);
        }
        {
            float *begin = record->getElement(keys._comp);
            std::copy(_comp.begin(), _comp.end(), begin);
        }
        {
            double *begin1 = record->getElement(keys._context_first);
            double *begin2 = record->getElement(keys._context_second);
            for (int i = 0; i != _context.size(); ++i) {
                begin1[i] = _context[i].first;
                begin2[i] = _context[i].second;
            }
        }

        handle.saveCatalog(cat);
    }
}

}} // namespace astromatic::psfex

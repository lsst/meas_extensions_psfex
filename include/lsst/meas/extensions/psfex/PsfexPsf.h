// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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

#ifndef ASTROMATIC_PSFEX_PSFEX_H
#define ASTROMATIC_PSFEX_PSFEX_H

#include "lsst/meas/algorithms/ImagePsf.h"
#include "lsst/meas/extensions/psfex/psf.hh"

namespace astromatic { namespace psfex {
    namespace detail {
        class PsfexPsfFactory;
    }

/**
 * @brief Represent a PSF as a linear combination of PSFEX (== Karhunen-Loeve) basis functions
 */
class PsfexPsf : public lsst::afw::table::io::PersistableFacade<PsfexPsf>,
                 public lsst::meas::algorithms::ImagePsf {
    friend class detail::PsfexPsfFactory;
public:
    /**
     *  @brief Constructor for a PsfexPsf
     */
    explicit PsfexPsf(
        astromatic::psfex::Psf const& psf, ///< [in] Psfex PSF model that we want to wrap into an LSST Psf
        lsst::afw::geom::Point2D const & averagePosition=lsst::afw::geom::Point2D()
                                        ///< [in] Average position of stars used to construct the Psf.
    );
    virtual ~PsfexPsf();

    /// Polymorphic deep copy; should usually be unnecessary as Psfs are immutable.x
    virtual PTR(lsst::afw::detection::Psf) clone() const;

    /// Return average position of stars; used as default position.
    virtual lsst::afw::geom::Point2D getAveragePosition() const { return _averagePosition; }
    
    /// Return the PSF's basis functions as a spatially-invariant LinearCombinationKernel
    /// with unit weights
    PTR(lsst::afw::math::LinearCombinationKernel const)
    getKernel(lsst::afw::geom::Point2D =
              lsst::afw::geom::Point2D(std::numeric_limits<double>::quiet_NaN())) const;

    /// Is this object persistable?
    virtual bool isPersistable() const { return true; }
    
    void write(lsst::afw::table::io::OutputArchiveHandle & handle) const;
private:
    lsst::afw::geom::Point2D _averagePosition;
    // Here are the unpacked fields from the psfex psf struct
    struct poly *_poly;                 // Polynom describing the PSF variations
    float _pixstep;                     // Mask oversampling (pixel)
    std::vector<int> _size;             // PSF dimensions
    std::vector<float> _comp;           // Complete pix. data (PSF components)
    std::vector<std::pair<double, double> > _context; // Offset/scale to apply to context data
    mutable PTR(lsst::afw::math::LinearCombinationKernel) _kernel; // keep a reference to getKernel()'s kernel

    /// default ctor; needed for persistence
    explicit PsfexPsf();

    /// Compute an image of the Psf at the specified position/colour, at pixel position in the output image
    virtual PTR(lsst::afw::detection::Psf::Image) _doComputeImage(
        lsst::afw::geom::Point2D const & position, ///< position within the chip
        lsst::afw::image::Color const& color,      ///< colour of object
        lsst::afw::geom::Point2D const& center     ///< position of center of image in the output image
    ) const;

    /// Compute an image of the Psf at the specified position/colour, at pixel (0.0, 0.0) in the output image
    virtual PTR(lsst::afw::detection::Psf::Image) doComputeKernelImage(
        lsst::afw::geom::Point2D const & position, ///< position within the chip
        lsst::afw::image::Color const& color       ///< colour of object
    ) const;

    /// Compute an image of the Psf at the specified position/colour, at pixel position in the output image
    virtual PTR(lsst::afw::detection::Psf::Image) doComputeImage(
        lsst::afw::geom::Point2D const & position, ///< position within the chip
        lsst::afw::image::Color const& color       ///< colour of object
    ) const;

    /// Name used in table persistence
    virtual std::string getPersistenceName() const;
    /// The python module name (for use in table persistence)
    virtual std::string getPythonModule() const;
};

}}

#endif // !ASTROMATIC_PSFEX_PSFEX_H

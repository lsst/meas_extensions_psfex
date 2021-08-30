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

# if !defined(LSST_MEAS_EXTENSIONS_PSFEX_PSFEX_H)
#define LSST_MEAS_EXTENSIONS_PSFEX_PSFEX_H 1

#include "lsst/geom/Box.h"
#include "lsst/meas/algorithms/ImagePsf.h"
#include "lsst/meas/extensions/psfex/psf.hh"

namespace lsst { namespace meas { namespace extensions { namespace psfex {
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
        lsst::meas::extensions::psfex::Psf const& psf, ///< [in] Psfex PSF model to be wrapped into an LSST Psf
        lsst::geom::Point2D const & averagePosition=lsst::geom::Point2D()
                                        ///< [in] Average position of stars used to construct the Psf.
    );
    virtual ~PsfexPsf();

    /// Polymorphic deep copy; should usually be unnecessary as Psfs are immutable.x
    virtual std::shared_ptr<lsst::afw::detection::Psf> clone() const;

    /// Return a clone with specified kernel dimensions
    virtual std::shared_ptr<afw::detection::Psf> resized(int width, int height) const;

    /// Return average position of stars; used as default position.
    virtual lsst::geom::Point2D getAveragePosition() const { return _averagePosition; }

    /// Return the PSF's basis functions as a spatially-invariant LinearCombinationKernel
    /// with unit weights
    std::shared_ptr<lsst::afw::math::LinearCombinationKernel const>
    getKernel(lsst::geom::Point2D =
              lsst::geom::Point2D(std::numeric_limits<double>::quiet_NaN())) const;

    /// Is this object persistable?
    virtual bool isPersistable() const noexcept override { return true; }

    void write(lsst::afw::table::io::OutputArchiveHandle & handle) const;
private:
    lsst::geom::Point2D _averagePosition;
    // Here are the unpacked fields from the psfex psf struct
    struct poly *_poly;                 // Polynom describing the PSF variations
    float _pixstep;                     // Mask oversampling (pixel)
    std::vector<int> _size;             // PSF dimensions
    std::vector<float> _comp;           // Complete pix. data (PSF components)
    std::vector<std::pair<double, double> > _context; // Offset/scale to apply to context data
    mutable std::shared_ptr<lsst::afw::math::LinearCombinationKernel> _kernel; // keep a reference to getKernel()'s kernel

    /// default ctor; needed for persistence
    explicit PsfexPsf();

    /// Compute an image of the Psf at the specified position/colour, at pixel position in the output image
    virtual std::shared_ptr<lsst::afw::detection::Psf::Image> _doComputeImage(
        lsst::geom::Point2D const & position, ///< position within the chip
        lsst::afw::image::Color const& color,      ///< colour of object
        lsst::geom::Point2D const& center     ///< position of center of image in the output image
    ) const;

    /// Compute an image of the Psf at the specified position/colour, at pixel (0.0, 0.0) in the output image
    virtual std::shared_ptr<lsst::afw::detection::Psf::Image> doComputeKernelImage(
        lsst::geom::Point2D const & position, ///< position within the chip
        lsst::afw::image::Color const& color       ///< colour of object
    ) const;

    /// Compute an image of the Psf at the specified position/colour, at pixel position in the output image
    virtual std::shared_ptr<lsst::afw::detection::Psf::Image> doComputeImage(
        lsst::geom::Point2D const & position, ///< position within the chip
        lsst::afw::image::Color const& color       ///< colour of object
    ) const;

    /// Compute the bbox of the kernel image at the specified position/color
    virtual lsst::geom::Box2I doComputeBBox(
        lsst::geom::Point2D const & position,
        lsst::afw::image::Color const & color
    ) const;

    /// Compute bbox of either image or kernel image, depending on provided center
    /// Does not depend on color, which is left out of parameter list to permit reuse
    /// by doComputeBBox, doCompute[Kernel]Image, and getKernel.
    lsst::geom::Box2I _doComputeBBox(
        lsst::geom::Point2D const & position, ///< position within the chip
        lsst::geom::Point2D const & center ///< center of output image. Use (0., 0.) for a kernel image
    ) const;

    /// Name used in table persistence
    virtual std::string getPersistenceName() const;
    /// The python module name (for use in table persistence)
    virtual std::string getPythonModule() const;
};

}}}}

#endif // !LSST_MEAS_EXTENSIONS_PSFEX_PSFEX_H

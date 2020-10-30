/*
 * LSST Data Management System
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 * See the COPYRIGHT file
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
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */
#include "pybind11/pybind11.h"

#include "lsst/afw/table/io/python.h"  // for addPersistableMethods

#include "define.h"
#include "vignet.h"
static double PSFEX_SAVE_BIG = BIG;	// we'll #undef BIG and define a variable called BIG
static double PSFEX_SAVE_INTERPFAC = INTERPFAC;

#undef PI

#include "lsst/meas/extensions/psfex/PsfexPsf.h"

#include "lsst/afw/table/io/OutputArchive.h"

#undef BIG
#undef INTERPFAC
double BIG = PSFEX_SAVE_BIG;
double INTERPFAC = PSFEX_SAVE_INTERPFAC;

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace extensions {
namespace psfex {

PYBIND11_MODULE(psfexPsf, mod) {
    py::module::import("lsst.afw.table");
    py::module::import("lsst.meas.algorithms");

    mod.attr("BIG") = py::cast(BIG);
    mod.attr("INTERPFAC") = py::cast(INTERPFAC);

    py::class_<PsfexPsf, std::shared_ptr<PsfexPsf>, lsst::meas::algorithms::ImagePsf> clsPsfexPsf(mod,
                                                                                                  "PsfexPsf");
    lsst::afw::table::io::python::addPersistableMethods<PsfexPsf>(clsPsfexPsf);


    clsPsfexPsf.def(py::init<lsst::meas::extensions::psfex::Psf const&, lsst::geom::Point2D const &>(),
            "psf"_a, "averagePosition"_a=lsst::geom::Point2D());

    clsPsfexPsf.def("clone", &PsfexPsf::clone);
    clsPsfexPsf.def("getAveragePosition", &PsfexPsf::getAveragePosition);
    clsPsfexPsf.def("getKernel", &PsfexPsf::getKernel,
            "position"_a=lsst::geom::Point2D(std::numeric_limits<double>::quiet_NaN()));
    clsPsfexPsf.def("isPersistable", &PsfexPsf::isPersistable);
    clsPsfexPsf.def("write", &PsfexPsf::write);
}

}  // psfex
}  // extensions
}  // meas
}  // lsst

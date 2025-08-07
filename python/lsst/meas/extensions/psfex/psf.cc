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
#include "pybind11/stl.h"
#include "lsst/cpputils/python.h"


#include "ndarray/pybind11.h"

#include "lsst/meas/extensions/psfex/psf.hh"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace extensions {
namespace psfex {

void wrapPsf(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyContext = py::class_<Context>;

    wrappers.wrapType(PyContext(wrappers.module, "Context"), [](auto &mod, auto &clsContext) {

        clsContext.attr("KEEPHIDDEN") = py::cast(static_cast<int>(Context::KEEPHIDDEN));
        clsContext.attr("REMOVEHIDDEN") = py::cast(static_cast<int>(Context::REMOVEHIDDEN));

        clsContext.def(
                py::init<std::vector<std::string> const &, std::vector<int> const &, std::vector<int> const &, bool>(),
                "names"_a, "group"_a, "degree"_a, "pcexflag"_a);

        clsContext.def("getName", &Context::getName);
        clsContext.def("getNpc", &Context::getNpc);
        clsContext.def("getPcflag", &Context::getPcflag);
        clsContext.def("getPc", &Context::getPc);
    });

    using PySample = py::class_<Sample>;
    wrappers.wrapType(PySample(wrappers.module, "Sample"), [](auto &mod, auto &clsSample) {

        clsSample.def("getCatindex", &Sample::getCatindex);
        clsSample.def("setCatindex", &Sample::setCatindex);
        clsSample.def("getObjindex", &Sample::getObjindex);
        clsSample.def("setObjindex", &Sample::setObjindex);
        clsSample.def("getExtindex", &Sample::getExtindex);
        clsSample.def("setExtindex", &Sample::setExtindex);
        clsSample.def("setVig", &Sample::setVig);
        clsSample.def("setNorm", &Sample::setNorm);
        clsSample.def("setBacknoise2", &Sample::setBacknoise2);
        clsSample.def("setGain", &Sample::setGain);
        clsSample.def("setX", &Sample::setX);
        clsSample.def("setY", &Sample::setY);
        clsSample.def("setContext", &Sample::setContext);
        clsSample.def("setFluxrad", &Sample::setFluxrad);
        clsSample.def("getVig", &Sample::getVig);
        clsSample.def("getVigResi", &Sample::getVigResi);
        clsSample.def("getVigChi", &Sample::getVigChi);
        clsSample.def("getVigWeight", &Sample::getVigWeight);
        clsSample.def("getXY", &Sample::getXY);
        clsSample.def("getNorm", &Sample::getNorm);
    });

    using PySet = py::classh<Set>;
    wrappers.wrapType(PySet(wrappers.module, "Set"), [](auto &mod, auto &clsSet) {
        clsSet.def(py::init<Context &>());

        clsSet.def("newSample", &Set::newSample);
        clsSet.def("trimMemory", &Set::trimMemory);
        clsSet.def("getFwhm", &Set::getFwhm);
        clsSet.def("setFwhm", &Set::setFwhm);
        clsSet.def("getNcontext", &Set::getNcontext);
        clsSet.def("getNsample", &Set::getNsample);
        clsSet.def("getContextOffset", &Set::getContextOffset);
        clsSet.def("setContextOffset", &Set::setContextOffset);
        clsSet.def("getContextScale", &Set::getContextScale);
        clsSet.def("setContextScale", &Set::setContextScale);
        clsSet.def("getRecentroid", &Set::getRecentroid);
        clsSet.def("setRecentroid", &Set::setRecentroid);
        clsSet.def("setVigSize", &Set::setVigSize);
        clsSet.def("finiSample", &Set::finiSample);
        clsSet.def("empty", &Set::empty);
        clsSet.def("getContextNames", &Set::getContextNames);
        clsSet.def("setContextname", &Set::setContextname);
        clsSet.def("setBadFlags", &Set::setBadFlags);
        clsSet.def("getBadFlags", &Set::getBadFlags);
        clsSet.def("setBadSN", &Set::setBadSN);
        clsSet.def("getBadSN", &Set::getBadSN);
        clsSet.def("setBadFrmin", &Set::setBadFrmin);
        clsSet.def("getBadFrmin", &Set::getBadFrmin);
        clsSet.def("setBadFrmax", &Set::setBadFrmax);
        clsSet.def("getBadFrmax", &Set::getBadFrmax);
        clsSet.def("setBadElong", &Set::setBadElong);
        clsSet.def("getBadElong", &Set::getBadElong);
        clsSet.def("setBadPix", &Set::setBadPix);
        clsSet.def("getBadPix", &Set::getBadPix);
        clsSet.def("getSample", &Set::getSample);
    });

    using PyPsf =     py::class_<Psf>;
    wrappers.wrapType(PyPsf(wrappers.module, "Psf"), [](auto &mod, auto &clsPsf) {
        clsPsf.def(py::init<>());

        clsPsf.def("getLoc", &Psf::getLoc);
        clsPsf.def("getResi", &Psf::getResi);
        clsPsf.def("build", &Psf::build,
                   "x"_a, "y"_a, "other"_a = std::vector<double>());
        clsPsf.def("clip", &Psf::clip);
    });
}

}  // psfex
}  // extensions
}  // meas
}  // lsst

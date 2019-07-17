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

#include "lsst/afw/table/io/python.h"  // for declarePersistableFacade

#include "lsst/meas/extensions/psfex/prefs.hh"
#include "lsst/daf/base/PropertySet.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace extensions {
namespace psfex {

PYBIND11_MODULE(prefs, mod) {
    py::class_<Prefs> clsPrefs(mod, "Prefs");

    clsPrefs.attr("ALL_EXTENSIONS") = py::cast(static_cast<int>(Prefs::ALL_EXTENSIONS));

    clsPrefs.def(py::init<std::string const&, lsst::daf::base::PropertySet const*>(),
            "filename"_a, "values"_a=nullptr);

    clsPrefs.def("use", &Prefs::use);
    clsPrefs.def("setCommandLine", &Prefs::setCommandLine);
    clsPrefs.def("getNcat", &Prefs::getNcat);
    clsPrefs.def("getPsfStep", &Prefs::getPsfStep);
    clsPrefs.def("getMinsn", &Prefs::getMinsn);
    clsPrefs.def("getMaxellip", &Prefs::getMaxellip);
    clsPrefs.def("getFwhmrange", &Prefs::getFwhmrange);
    clsPrefs.def("getPsfsize", &Prefs::getPsfsize);
    clsPrefs.def("getAutoselectFlag", &Prefs::getAutoselectFlag);
    clsPrefs.def("getFlagMask", &Prefs::getFlagMask);
    clsPrefs.def("getMaxvar", &Prefs::getMaxvar);
    clsPrefs.def("getVarType", &Prefs::getVarType);
    clsPrefs.def("getBadpixNmax", &Prefs::getBadpixNmax);
    clsPrefs.def("getBadpixFlag", &Prefs::getBadpixFlag);
    clsPrefs.def("getCenterKey", &Prefs::getCenterKey);
    clsPrefs.def("getPhotfluxRkey", &Prefs::getPhotfluxRkey);
    clsPrefs.def("getPhotfluxNum", &Prefs::getPhotfluxNum);
    clsPrefs.def("getPhotfluxerrRkey", &Prefs::getPhotfluxerrRkey);
    clsPrefs.def("getPhotfluxerrNum", &Prefs::getPhotfluxerrNum);
    clsPrefs.def("getProfAccuracy", &Prefs::getProfAccuracy);
    clsPrefs.def("getVerboseType", &Prefs::getVerboseType);
    clsPrefs.def("getContextName", &Prefs::getContextName);
    clsPrefs.def("getContextGroup", &Prefs::getContextGroup);
    clsPrefs.def("getGroupDeg", &Prefs::getGroupDeg);
    clsPrefs.def("addCatalog", &Prefs::addCatalog);
    clsPrefs.def("getCatalogs", &Prefs::getCatalogs);
}

}  // psfex
}  // extensions
}  // meas
}  // lsst

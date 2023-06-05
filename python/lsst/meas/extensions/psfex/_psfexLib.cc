/*
 * This file is part of meas_extensions_psfex.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
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
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "pybind11/pybind11.h"
#include "lsst/cpputils/python.h"

namespace lsst {
namespace meas {
namespace extensions {
namespace psfex {

using lsst::cpputils::python::WrapperCollection;
void wrapField(WrapperCollection &wrappers);
void wrapPrefs(WrapperCollection &wrappers);
void wrapPsfexPsf(WrapperCollection &wrappers);
void wrapPsf(WrapperCollection &wrappers);

PYBIND11_MODULE(_psfexLib, mod) {
    lsst::cpputils::python::WrapperCollection wrappers(mod, "lsst.meas.extensions.psfex");
    wrappers.addInheritanceDependency("lsst.meas.algorithms");
    wrappers.addSignatureDependency("lsst.afw.table");
    wrapField(wrappers);
    wrapPrefs(wrappers);
    wrapPsfexPsf(wrappers);
    wrapPsf(wrappers);
    wrappers.finish();
}

}
}
}
}

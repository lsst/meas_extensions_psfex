# This file is part of meas_extensions_psfex.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import unittest
import os

from lsst.utils import getPackageDir

from lsst.afw.table import SourceCatalog
from lsst.afw.image import ExposureF
from lsst.pipe.base import Task, TaskMetadata
from lsst.meas.algorithms import MakePsfCandidatesTask, ReserveSourcesTask
from lsst.meas.extensions.psfex.psfexPsfDeterminer import PsfexPsfDeterminerTask


CONFIG_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), "config")

try:
    DATA_DIR = os.path.join(getPackageDir("afwdata"), "psfex")
except LookupError:
    DATA_DIR = None


@unittest.skipUnless(DATA_DIR is not None, "afwdata is not setup")
class TestMklBadness(unittest.TestCase):
    """Test a particular visit/detector combination for which a bad PSF was
    generated when PSFEx was run with MKL in the environment.

    See DM-40066 for details.
    """

    def setUp(self) -> None:
        self.stars = SourceCatalog.readFits(os.path.join(DATA_DIR, "stars.fits"))
        self.exposure = ExposureF(os.path.join(DATA_DIR, "exposure.fits"))
        makePsfCandidatesConfig = MakePsfCandidatesTask.ConfigClass()
        makePsfCandidatesConfig.load(os.path.join(CONFIG_DIR, "make_psf_candidates.py"))
        self.makePsfCandidates = MakePsfCandidatesTask(config=makePsfCandidatesConfig)
        reserveSourcesConfig = ReserveSourcesTask.ConfigClass()
        reserveSourcesConfig.load(os.path.join(CONFIG_DIR, "reserve_sources.py"))
        # ResourceSourceTask.__init__ insists on adding a new column to the
        # catalog, but this catalog already has that column; hack it.
        self.reserve = ReserveSourcesTask.__new__(ReserveSourcesTask)
        Task.__init__(self.reserve, config=reserveSourcesConfig)
        self.reserve.columnName = "calib_psf"
        self.reserve.key = self.stars.schema.find("calib_psf_reserved").key
        psfDeterminerConfig = PsfexPsfDeterminerTask.ConfigClass()
        psfDeterminerConfig.load(os.path.join(CONFIG_DIR, "psf_determiner.py"))
        self.psfDeterminer = PsfexPsfDeterminerTask(config=psfDeterminerConfig, schema=self.stars.schema)
        self.metadata = TaskMetadata()

    def test_mkl_badness(self) -> None:
        selectionResult = self.makePsfCandidates.run(self.stars, exposure=self.exposure)
        reserveResult = self.reserve.run(selectionResult.goodStarCat, expId=3581242)
        psfDeterminerList = [
            cand
            for cand, use in zip(selectionResult.psfCandidates, reserveResult.use)
            if use
        ]
        self.psfDeterminer.determinePsf(
            self.exposure, psfDeterminerList, self.metadata,
        )
        self.assertEqual(self.metadata.getScalar("numAvailStars"), 84)
        self.assertEqual(self.metadata.getScalar("numGoodStars"), 77)


if __name__ == "__main__":
    unittest.main()

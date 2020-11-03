.. _lsst.meas.extensions.psfex:

##########################
lsst.meas.extensions.psfex
##########################

``lsst.meas.extensions.psfex`` provides algorithms for PSF estimation based on the PSFEx code.


.. _lsst.meas.extensions.psfex-using:

Using lsst.meas.extensions.psfex
================================

.. toctree::
   :maxdepth: 1

``lsst.meas.extensions.psfex`` is developed at https://github.com/lsst/meas_extensions_psfex.
You can find Jira issues for this module under the `meas_extensions_psfex <https://jira.lsstcorp.org/issues/?jql=project%20%3D%20DM%20AND%20component%20%3D%20meas_extensions_psfex>`_ component.


Python API reference
====================

.. automodapi:: lsst.meas.extensions.psfex
   :no-main-docstr:
   :no-inheritance-diagram:

.. This fails because the module can't be imported; see DM-16899 and DM-9840.
.. .. automodapi:: lsst.meas.extensions.psfex.psfex
..   :no-main-docstr:
..   :no-inheritance-diagram:

.. automodapi:: lsst.meas.extensions.psfex.psfexPsfDeterminer
   :no-main-docstr:
   :no-inheritance-diagram:

.. Also fails; see comment on psfex.psfex above.
.. .. automodapi:: lsst.meas.extensions.psfex.psfexStarSelector
..   :no-main-docstr:
..   :no-inheritance-diagram:

import lsst.meas.algorithms.reserveSourcesTask
assert type(config)==lsst.meas.algorithms.reserveSourcesTask.ReserveSourcesConfig, 'config is of type %s.%s instead of lsst.meas.algorithms.reserveSourcesTask.ReserveSourcesConfig' % (type(config).__module__, type(config).__name__)
# Fraction of candidates to reserve from fitting; none if <= 0
config.fraction=0.2

# This number will be added to the exposure ID to set the random seed for reserving candidates
config.seed=1


// -*- lsst-C++ -*-
#include "lsst/meas/extensions/psfex/psf.hh"
#include "ndarray.h"

#define RETURN_IMAGE_FIELD(NAME, CNAME, SIZE) \
     ndarray::Array<float,2,2> \
     NAME() const         \
     { \
        ndarray::Array<float,2,2>::Index shape = ndarray::makeVector(SIZE[0], SIZE[1]); \
        ndarray::Array<float,2,2>::Index strides = ndarray::makeVector(1, SIZE[0]); \
        \
        return ndarray::external(impl->CNAME, shape, strides); \
    }

namespace astromatic { namespace psfex {

Context::Context(std::vector<std::string> const& names,
                 std::vector<int> const& group,
                 std::vector<int> const& degree,
                 bool pcexflag
                )
{
    _names = names;
    int const ndim = names.size();
    int const ngroup = group.size();

    std::vector<const char *> cnames(ndim);
    for (int i = 0; i != ndim; ++i) {
        cnames[i] = names[i].c_str();
    }

    impl = context_init(const_cast<char **>(&cnames[0]),
                        const_cast<int *>(&group[0]),
                        ndim,
                        const_cast<int *>(&degree[0]),
                        ngroup, (pcexflag ? 1 : 0));
    _pc_vectors.resize(impl->npc);
}

Context::~Context() {
    context_end(impl);
}

std::vector<double> & Context::getPc(int const i) const
{
    if (i < 0 || i >= impl->npc) {
        std::ostringstream s1;
        s1 << "Index " << i << " is out of range 0.." << impl->npc - 1;
        throw std::out_of_range(s1.str());
    }
    if (_pc_vectors[i].empty()) {
        _pc_vectors[i].reserve(impl->npc);
        for (int j = 0; j != impl->npc; ++j) {
            _pc_vectors[i][j] = impl->pc[i*impl->npc + j];
        }
    }
    return _pc_vectors[i];
}

/************************************************************************************************************/

void Sample::setVig(ndarray::Array<float,2,2> const& img)
{
    memcpy(impl->vig, &img(0,0), img.getSize<0>()*img.getSize<1>()*sizeof(float));
}

RETURN_IMAGE_FIELD(Sample::getVig,       vig,       _vigsize)
RETURN_IMAGE_FIELD(Sample::getVigResi,   vigresi,   _vigsize)
RETURN_IMAGE_FIELD(Sample::getVigChi,    vigchi,    _vigsize)
RETURN_IMAGE_FIELD(Sample::getVigWeight, vigweight, _vigsize)
        
/************************************************************************************************************/

Set::Set(Context &c) {
    impl = init_set(c.impl);
    impl->nsample = impl->nsamplemax = 0;
}
Set::~Set() {
    end_set(impl);
}

Sample
Set::newSample()
{
    if (impl->nsample >= impl->nsamplemax) {
        if (impl->nsamplemax == 0) {
            malloc_samples(impl, LSAMPLE_DEFSIZE);
        } else {
            realloc_samples(impl, (int)(1.62*impl->nsamplemax + 1));
        }
    }
    return Sample(impl->sample[impl->nsample++], impl->vigsize);
}
        
void
Set::trimMemory() const
{
    setstruct *mutable_impl = const_cast<setstruct *>(impl);
    realloc_samples(mutable_impl, impl->nsample);
}

void
Set::finiSample(Sample const& sample)
{
    make_weights(impl, sample.impl);
    recenter_sample(sample.impl, impl, sample._fluxrad);
}

/************************************************************************************************************/

Psf::~Psf() {
#if 0                                   // we don't own impl
    psf_end(impl);
#endif
}

void
Psf::build(double x, double y,
           std::vector<double> const& other)
{
    if (!impl->contextoffset) {         // we never set context{offset,scale}
        throw std::out_of_range("psf.contextoffset not set (probably no valid stars) so I can't index it");
    }

    std::vector<double> pos(2);
    pos[0] = x; pos[1] = y;

    pos.insert(pos.end(), other.begin(), other.end());
    
    for (int i = 0; i != pos.size(); ++i) {
        pos[i] = (pos[i] - impl->contextoffset[i])/impl->contextscale[i];
    }

    psf_build(impl, &pos[0]);
}		

RETURN_IMAGE_FIELD(Psf::getLoc,  loc,  impl->size)
RETURN_IMAGE_FIELD(Psf::getResi, resi, impl->size)

}}

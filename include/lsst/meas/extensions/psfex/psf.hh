// -*- lsst-C++ -*-
#if !defined(ASTROMATIC_PSFEX_PSF_HH)
#define ASTROMATIC_PSFEX_PSF_HH

#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>
#include "boost/shared_ptr.hpp"

extern "C" {
#include "context.h"
#include "psf.h"
#include "sample.h"
typedef struct field fieldstruct;
}

namespace ndarray {
    template<typename T, int n1, int n2> class Array;
}

namespace astromatic { namespace psfex {

class Field;
class Psf;
class Set;

/**
 * \brief The parameters that describe the PSF's variability
 */
class Context {
    friend class Psf;
    friend class Set;
public:
    enum { KEEPHIDDEN=CONTEXT_KEEPHIDDEN, REMOVEHIDDEN=CONTEXT_REMOVEHIDDEN };
        
    Context(std::vector<std::string> const& names, ///< names of fields to use
            std::vector<int> const& group,         ///< tags for each member of names
            std::vector<int> const& degree,        ///< polynomial degree for each group
            bool pcexflag                          ///< exclude PCA components?
           );
    ~Context();
    std::vector<std::string> const& getName() const { return _names; }
    int getNpc() const { return impl->npc; }
    int getPcflag(int i) const { return impl->pcflag[i]; }
    std::vector<double> & getPc(int const i) const;
private:
    contextstruct *impl;
    mutable std::vector<std::vector<double> > _pc_vectors;
    std::vector<std::string> _names;
};

class Sample {
    friend class Set;
public:
    ~Sample() { }

    int getCatindex() const { return impl->catindex; }
    void setCatindex(int val) { impl->catindex = val; }
    int getExtindex() const { return impl->extindex; }
    void setExtindex(int val) { impl->extindex = val; }
    void setVig(ndarray::Array<float,2,2> const& img);
    void setNorm(float val) { impl->norm = val; }
    void setBacknoise2(float val) { impl->backnoise2 = val; }
    void setGain(float val) { impl->gain = val; }
    void setX(double val) { impl->x = val; impl->dx = impl->x - (int)(impl->x + 0.49999); }
    void setY(double val) { impl->y = val; impl->dy = impl->y - (int)(impl->y + 0.49999); }
    void setContext(int i, double val) { impl->context[i] = val; }
    void setFluxrad(float val) { _fluxrad = val; }

    ndarray::Array<float,2,2> getVig() const;
    ndarray::Array<float,2,2> getVigResi() const;
    ndarray::Array<float,2,2> getVigChi() const;
    ndarray::Array<float,2,2> getVigWeight() const;
    std::pair<double, double> getXY() const { return std::pair<double,double>(impl->x, impl->y); }
    float getNorm() const { return impl->norm; }
    
private:
    Sample(samplestruct *s, int *vigsize) : impl(s), _vigsize(vigsize) { }

    samplestruct *impl;
    float _fluxrad;                     // needed by recenter_sample
    int *_vigsize;                      // needed to make sense of the arrays
};

/************************************************************************************************************/
    
class Set {
    friend class Psf;
    friend void makeit(std::vector<boost::shared_ptr<Field> > &fields,
                       std::vector<boost::shared_ptr<Set> > const& sets);
public:
    Set(Context &c);
    ~Set();
    Sample newSample();
    void trimMemory() const;

    int getFwhm() const { return impl->fwhm; }
    void setFwhm(double fwhm) { impl->fwhm = fwhm; }
    int getNcontext() const { return impl->ncontext; }
    int getNsample() const { return impl->nsample; }
    double getContextOffset(int i) const { return impl->contextoffset[i]; }
    void   setContextOffset(int i, double val) { impl->contextoffset[i] = val; }
    double getContextScale(int i) const { return impl->contextscale[i]; }
    void   setContextScale(int i, double val) { impl->contextscale[i] = val; }
    void setVigSize(int w, int h) {
        impl->vigsize[0] = w;
        impl->vigsize[1] = h;
        impl->nvig = w*h;
    }
    void finiSample(Sample const& sample);
    bool empty() const {
        return impl->nsample == 0;
    }
    std::vector<const char *> const& getContextNames() const {
        return _contextNames;
    }

    void setContextname(int i, std::string const& s) {
        strcpy(impl->contextname[i], s.c_str());
        _contextNames.push_back(impl->contextname[i]);
    }
    void setBadFlags(int n) { impl->badflags = n; }
    int  getBadFlags() const { return impl->badflags; }
    void setBadSN(int n) { impl->badsn = n; }
    int  getBadSN() const { return impl->badsn; }
    void setBadFrmin(int n) { impl->badfrmin = n; }
    int  getBadFrmin() const { return impl->badfrmin; }
    void setBadFrmax(int n) { impl->badfrmax = n; }
    int  getBadFrmax() const { return impl->badfrmax; }
    void setBadElong(int n) { impl->badelong = n; }
    int  getBadElong() const { return impl->badelong; }
    void setBadPix(int n) { impl->badpix = n; }
    int  getBadPix() const { return impl->badpix; }
    Sample getSample(int i) const {
        if (i < 0 || i > impl->nsample) {
            std::ostringstream s1;
            s1 << "Index " << i << " is out of range 0.." << impl->nsample - 1;
            throw std::out_of_range(s1.str());
        }

        return Sample(impl->sample[i], impl->vigsize);
    }

private:
    setstruct *impl;
    std::vector<const char *> _contextNames;
};

/** \brief PSF
 */
class Psf {
    friend class PsfexPsf;
public:
    Psf() : impl(0), _owner(boost::shared_ptr<fieldstruct>()) {}
    Psf(psfstruct *psf, boost::shared_ptr<fieldstruct> owner) : impl(psf), _owner(owner) {}
    ~Psf();
    
    ndarray::Array<float,2,2> getLoc() const;
    ndarray::Array<float,2,2> getResi() const;
#if 0
    void make(Set &s, double prof_accuracy) {
        psf_make(impl, s.impl, prof_accuracy);
    }
    void makeresi(Set &s, double prof_accuracy, int centflag=0) {
        psf_makeresi(impl, s.impl, centflag, prof_accuracy);
    }
#endif
    void build(double x, double y, std::vector<double> const& other=std::vector<double>());

    void clip() {
	psf_clip(impl);
    }
    
protected:
    psfstruct *impl;
private:
    boost::shared_ptr<fieldstruct> _owner;
};

}}

#endif

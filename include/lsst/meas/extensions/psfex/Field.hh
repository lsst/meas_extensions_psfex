// -*- lsst-C++ -*-
#if !defined(ASTROMATIC_PSFEX_FIELD_HH)
#define ASTROMATIC_PSFEX_FIELD_HH

#include <string>
#include <vector>
#include "psf.hh"

namespace lsst {
   namespace daf { namespace base {
      class PropertySet;
   }}
   namespace afw { namespace image {
      class Wcs;
   }}
}

extern "C" {
struct structcat; typedef struct structcat catstruct;
struct structtab; typedef struct structtab tabstruct;

#include "define.h"
#include "field.h"
}

namespace astromatic { namespace psfex {
class Psf;
class Set;

/**
 * \brief Store all we know about for a visit to a field (maybe including multiple chips)
 */
class Field {
    friend void makeit(std::vector<boost::shared_ptr<Field> > &fields,
                       std::vector<boost::shared_ptr<Set> > const& sets);
public:
    Field(std::string const& ident="unknown" ///< Name of Field
         );
    ~Field();
    //
    void finalize() { _finalize(true); }        
    
    void addExt(lsst::afw::image::Wcs const& wcs, int const naxis1, int const naxis2, int const nobj=0);
    
    /// Return the number of extensions
    int getNext() const { return impl->next; }
    /// Return the Psfs
    std::vector<Psf> getPsfs() const;
    
private:
   boost::shared_ptr<fieldstruct> impl;
   mutable std::vector<Psf> _psfs;
   mutable bool _isInitialized;

   void _finalize(const bool force=false);
};

void makeit(std::vector<boost::shared_ptr<Field> > &fields,
            std::vector<boost::shared_ptr<Set> > const& sets);

}}

#endif

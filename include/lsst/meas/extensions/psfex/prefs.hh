// -*- lsst-C++ -*-
#if !defined(ASTROMATIC_PSFEX_PREFS_HH)
#define ASTROMATIC_PSFEX_PREFS_HH

#include <string>
#include <vector>

namespace lsst { namespace daf { namespace base {
    class PropertySet;
}}}

extern "C" {
struct structcat; typedef struct structcat catstruct;
struct structtab; typedef struct structtab tabstruct;

#include "define.h"
#include "prefs.h"
}

namespace astromatic { namespace psfex {
/**
 * \brief Tuning parameters
 */
class Prefs {
public:
#if 0
   enum {NEWBASIS_NONE = prefstruct::NEWBASIS_NONE,
	 NEWBASIS_PCAINDEPENDENT = prefstruct::NEWBASIS_PCAINDEPENDENT,
	 NEWBASIS_PCACOMMON = prefstruct::NEWBASIS_PCACOMMON};
   enum	{HIDDEN_MEF_INDEPENDENT = prefstruct::HIDDEN_MEF_INDEPENDENT,
	 HIDDEN_MEF_COMMON = prefstruct::HIDDEN_MEF_COMMON };
   enum	{STABILITY_EXPOSURE = prefstruct::STABILITY_EXPOSURE,
	 STABILITY_SEQUENCE = prefstruct::STABILITY_SEQUENCE};
   enum	{PSF_MEF_INDEPENDENT = prefstruct::PSF_MEF_INDEPENDENT,
	 PSF_MEF_COMMON = prefstruct::PSF_MEF_COMMON};
   enum	{HOMOBASIS_NONE = prefstruct::HOMOBASIS_NONE,
	 HOMOBASIS_GAUSSLAGUERRE = prefstruct::HOMOBASIS_GAUSSLAGUERRE};
#endif
   enum {QUIET = prefstruct::QUIET, NORM = prefstruct::NORM,
	 LOG = prefstruct::LOG, FULL = prefstruct::FULL};
   enum {VAR_NONE =  prefstruct::VAR_NONE, VAR_SEEING = prefstruct::VAR_SEEING};

   enum {_ALL_EXTENSIONS = ALL_EXTENSIONS,
         #undef ALL_EXTENSIONS
	 ALL_EXTENSIONS = _ALL_EXTENSIONS};

   Prefs(std::string const& filename,	///< Filename
	 lsst::daf::base::PropertySet const* values=NULL ///< overrides
      );
   ~Prefs();
   //
   void use() {
      useprefs();
   }

   void setCommandLine(std::vector<std::string> const& argv);
   int getNcat() const { return prefs.ncat; }
   double getPsfStep() const { return prefs.psf_step; }
#if 0
   int getNgroupDeg() const { return prefs.ngroup_deg; }
   int getNewbasisType() const { return prefs.newbasis_type; }
   int getStabilityType() const { return prefs.stability_type; }
   int getHiddenMefType() const { return prefs.hidden_mef_type; }
   int getPsfMefType() const { return prefs.psf_mef_type; }
   int getHomobasisType() const { return prefs.homobasis_type; }
#endif
   double getMinsn() const { return prefs.minsn; }
   double getMaxellip() const { return prefs.maxellip; }
   std::pair<double, double> getFwhmrange() const {
      return std::make_pair(prefs.fwhmrange[0], prefs.fwhmrange[1]);
   }
    std::pair<int, int> getPsfsize() const {
        return std::make_pair(prefs.psf_size[0], prefs.psf_size[1]);
   }
   int getAutoselectFlag() const { return prefs.autoselect_flag; }
   int getFlagMask() const { return prefs.flag_mask; }
   double getMaxvar() const { return prefs.maxvar; }
   int getVarType() const { return prefs.var_type; }
   int getBadpixNmax() const { return prefs.badpix_nmax; }
   int getBadpixFlag() const { return prefs.badpix_flag; }
   char *getCenterKey(int i) const { return prefs.center_key[i]; }
   char *getPhotfluxRkey() const { return prefs.photflux_rkey; }
   int   getPhotfluxNum() const { return prefs.photflux_num; }
   char *getPhotfluxerrRkey() const { return prefs.photfluxerr_rkey; }
   int   getPhotfluxerrNum() const { return prefs.photfluxerr_num; }
   double getProfAccuracy() const { return prefs.prof_accuracy; }
   int   getVerboseType() const { return prefs.verbose_type; }
   
   std::vector<std::string> const& getContextName() const { return _context_names; }
   std::vector<int> const& getContextGroup() const { return _context_groups; }
   std::vector<int> const& getGroupDeg() const { return _group_degs; }

   void addCatalog(std::string const& filename);
   std::vector<std::string> const& getCatalogs() const { return _catalogs; }

private:
   char const** _command_line;		     // argv passed from the unix command line
   std::vector<std::string> _catalogs;	     // names of catalogs
   std::vector<std::string> _context_names;  // names of context-keys
   std::vector<int> _context_groups;	     // Context group
   std::vector<int> _group_degs;	     // Degree for each group
};


}}

#endif

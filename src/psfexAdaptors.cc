#include "boost/format.hpp"
#include "lsst/pex/exceptions.h"
#include "lsst/meas/extensions/psfex/Field.hh"

extern "C" {
#include "globals.h"
#include "context.h"
#include "prefs.h"
#include "sample.h"

/************************************************************************************************************/

setstruct *
load_samples(char **filenames, int catindex, int ncat, int ext,
             int next, contextstruct *context)
{
    /*
     * The C version of this is called two ways:
     *   catindex == 0, ncat == ncat            Read all catalogues
     *   catindex == c, ncat == 1               Read only catalogue c
     */
    setstruct *completeSet = reinterpret_cast<setstruct *>(filenames[catindex + 0]);
    /*
     * Make a new set, which may be a subset of the completeSet
     */
    setstruct *set = init_set(context);
    set->fwhm = completeSet->fwhm;
    for (int i = 0; i != completeSet->vigdim; ++i) {
        set->vigsize[i] = completeSet->vigsize[i];
    }
    for (int i = 0; i != completeSet->ncontext; ++i) {
        strcpy(set->contextname[i], completeSet->contextname[i]);
        set->contextoffset[i] = completeSet->contextoffset[i];
        set->contextscale[i] = completeSet->contextscale[i];
    }
    /*
     * Count how many samples we'll be including
     */
    int nsample_keep = 0;
    for (int i = 0; i != ncat; ++i) {
        setstruct *s = reinterpret_cast<setstruct *>(filenames[catindex + i]);
        for (int j = 0; j != completeSet->nsample; ++j) {
            samplestruct const *samp = s->sample[j];
            if (ext == ALL_EXTENSIONS || ext == samp->extindex) {
                ++nsample_keep;
            }
        }
    }

    set->samples_owner = 0;
    malloc_samples(set, nsample_keep);
    for (int i = 0; i != ncat; ++i) {
        setstruct *s = reinterpret_cast<setstruct *>(filenames[catindex + i]);
        for (int j = 0; j != completeSet->nsample; ++j) {
            samplestruct *samp = s->sample[j];
            if (ext == ALL_EXTENSIONS || ext == samp->extindex) {
                set->sample[set->nsample++] = samp;
            }
        }
    }

    return set;
}

}

/************************************************************************************************************/

namespace lsst { namespace meas { namespace extensions { namespace psfex {

void
makeit(std::vector<boost::shared_ptr<Field> > &fields_,
       std::vector<boost::shared_ptr<Set> > const& sets
      )
{
    if (sets.size() > MAXFILE) {
        throw LSST_EXCEPT(lsst::pex::exceptions::LengthError,
                          (boost::format("Too many sets: %d") % sets.size()).str());
    }

    std::vector<fieldstruct *> fields(fields_.size());
    for (unsigned int i = 0; i != fields.size(); ++i) {
        fields[i] = fields_[i]->impl;
    }
    /*
     * We are going to scribble on prefs.incat_name to replace the array of (char*) with
     * an array of data
     */

    /// RAII for the scribbling
    struct ScribbleRaii {
        int const ncat;                 // Original number
        std::vector<char *> incat_name; // Original data
        size_t const setsSize;          // New size
        ScribbleRaii(std::vector<boost::shared_ptr<Set> > const& sets) :
            ncat(prefs.ncat), incat_name(ncat), setsSize(sets.size())
        {
            for (int i = 0; i != prefs.ncat; ++i) {
                incat_name[i] = prefs.incat_name[i];
            }
            for (unsigned int i = 0; i != sets.size(); ++i) {
                prefs.incat_name[i] = reinterpret_cast<char *>(sets[i]->impl);
            }
            prefs.ncat = setsSize;
        }
        ~ScribbleRaii() {
            for (int i = 0; i != prefs.ncat; ++i) {
                prefs.incat_name[i] = incat_name[i];
            }
            for (size_t i = prefs.ncat; i < setsSize; ++i) {
                prefs.incat_name[i] = NULL;
            }
            prefs.ncat = ncat;
        }
    };


    contextstruct *context = NULL, *fullcontext = NULL;
    {
        ScribbleRaii scribble(sets);
        makeit_body(&fields[0], &context, &fullcontext, false);
    }

    if (context->npc) {
        context_end(fullcontext);
    }
    context_end(context);   
}

}}}}

#!/usr/bin/env python
import argparse
import os
import sys
from lsst.meas.extensions.psfex import psfex, readPrefs, makeitLsst, makeit, showPsf, dafBase

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Wrapper for Point Spread Function Extractor (PSFEX)")

    parser.add_argument('catalogs', type=str, nargs="+", help="Input catalogues from SExtractor")
    parser.add_argument('-c', type=str, dest="defaultsFile",
                        help="File containing default parameters", default="default-lsst.psfex")
    parser.add_argument('--lsst', action="store_true", help="Read LSST data")
    parser.add_argument('--overrides', type=str, nargs="+",
                        help="Overrides for default parameters", default=[])
    parser.add_argument('--plot', type=str, nargs="+",
                        help="Desired plots", default=[])
    parser.add_argument('--ds9', type=int,
                        help="Show the PSF on ds9", default=None)
    parser.add_argument('--diagnostics', action="store_true",
                        help="Write output diagnostic plots to outDir")
    parser.add_argument('--verbose', action="store_true", help="How chatty should I be?", default=False)

    argv = sys.argv[:]                  # argparse will mess with sys.argv
    args = parser.parse_args()

    args_md = dafBase.PropertySet()
    for x in args.overrides:
        try:
            k, v = x.split('=')
        except ValueError:
            print >> sys.stderr, "Overrides must be of the form key=value, saw %s" % x
            continue
        args_md.set(k, v)

    plotKeys = ["fwhmHistogram", "showFlags", "showRejection"]
    if "help" in args.plot:
        print "Valid plot types are %s" % " ".join(["none"] + plotKeys)
        sys.exit(0)
    plot = {}
    if "none" not in args.plot:
        for k in args.plot:
            if k not in plotKeys:
                print >> sys.stderr, "Unknown plot type %s (Valid types are %s)" % (k, " ".join(plotKeys))
                sys.exit(1)
            plot[k] = True

    if args.lsst:
        psfex.psfex.setDataType("LSST")
    #
    # To work
    #
    prefs = readPrefs(args.defaultsFile, args_md)
    prefs.setCommandLine(argv)

    for f in args.catalogs:
        prefs.addCatalog(f)

    prefs.use()

    context = psfex.Context(prefs.getContextName(), prefs.getContextGroup(),
                            prefs.getGroupDeg(),
                            psfex.Context.REMOVEHIDDEN if False else psfex.Context.KEEPHIDDEN)

    psfs, sets, wcss = makeitLsst(prefs, context, saveWcs=True, plot=plot) if args.lsst else \
        makeit(prefs, context, saveWcs=True, plot=plot)

    if args.diagnostics or args.ds9 is not None:
        ds9Frame = args.ds9
        for i in range(len(sets)):
            for ext in range(len(psfs[i])):
                catDir, catFile = os.path.split(prefs.getCatalogs()[i])

                showPsf(psfs[i], sets[i], ext, wcss[i], nspot=3, trim=5,
                        frame=ds9Frame, diagnostics=args.diagnostics,
                        outDir=catDir, title=os.path.splitext(catFile)[0])

                if ds9Frame is not None:
                    ds9Frame += 2

#!/usr/bin/env python
import argparse
import os
import sys
from lsst.meas.extensions.psfex import psfex, makeitLsst, makeit, showPsf, dafBase


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Wrapper for Point Spread Function Extractor (PSFEX)")

    parser.add_argument('catalogs', type=str, nargs="+", help="Input catalogues from SExtractor")
    parser.add_argument('-c', type=str, dest="defaultsFile",
                        help="File containing default parameters", default="default-lsst.psfex")
    parser.add_argument('--overrides', type=str, nargs="+",
                        help="Overrides for default parameters", default=[])
    parser.add_argument('--plot', type=str, nargs="+",
                        help="Desired plots", default=[])
    parser.add_argument('--doDisplay', type=int,
                        help="Show the PSF on the display", default=None)
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
            print("Overrides must be of the form key=value, saw %s" % x, file=sys.stderr)
            continue
        args_md.set(k, v)

    plotKeys = ["fwhmHistogram", "showFlags", "showRejection"]
    if "help" in args.plot:
        print("Valid plot types are %s" % " ".join(["none"] + plotKeys))
        sys.exit(0)
    plot = {}
    if "none" not in args.plot:
        for k in args.plot:
            if k not in plotKeys:
                print("Unknown plot type %s (Valid types are %s)" % (k, " ".join(plotKeys)), file=sys.stderr)
                sys.exit(1)
            plot[k] = True

    #
    # To work
    #
    prefs = psfex.Prefs(args.defaultsFile, args_md)
    prefs.setCommandLine(argv)

    for f in args.catalogs:
        prefs.addCatalog(f)

    prefs.use()

    context = psfex.Context(prefs.getContextName(), prefs.getContextGroup(),
                            prefs.getGroupDeg(),
                            psfex.Context.REMOVEHIDDEN if False else psfex.Context.KEEPHIDDEN)

    psfs, sets, wcss = makeitLsst(prefs, context, saveWcs=True, plot=plot) if args.lsst else \
        makeit(prefs, context, saveWcs=True, plot=plot)

    if args.diagnostics or args.doDisplay is not None:
        dispFrame = args.doDisplay
        for i in range(len(sets)):
            for ext in range(len(psfs[i])):
                catDir, catFile = os.path.split(prefs.getCatalogs()[i])

                showPsf(psfs[i], sets[i], ext, wcss[i], nspot=3, trim=5,
                        frame=dispFrame, diagnostics=args.diagnostics,
                        outDir=catDir, title=os.path.splitext(catFile)[0])

                if dispFrame is not None:
                    dispFrame += 2

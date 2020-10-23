Installing's a bit of a pain for now.

# rpm install atlas-devel, the version seems to be 3.8.4
mkdir -p /home/astro/hsc/products/Linux64/external/atlas/3.8.4/lib
cd       /home/astro/hsc/products/Linux64/external/atlas/3.8.4/lib
ln -s /usr/lib64/atlas/lib* .
cd ..
eups declare -m none -r . atlas 3.8.4

export LSST_CFG_PATH=~/LSST/devenv/buildFiles
# or:
#   mkdir ups
#   cp ~/LSST/devenv/buildFiles/atlas/atlas.cfg ups

# N.b. you need a version of fftw that has the single as well as double precision libraries

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

mkdir ~/Astromatic
cd ~/Astromatic
svn co https://astromatic.net/pubsvn/software/psfex/branches/rhl psfex
cd psfex

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

cd ~/LSST/meas/extensions
git clone git@github.com:LSST-nonproject/meas_extensions_psfex.git psfex
cd psfex
cp -r files_for_vanilla_psfex/{SConstruct,lib,ups} ~/Astromatic/psfex

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

cd ~/Astromatic/psfex

configure CC='cc -Wno-empty-body -Wno-format -Wno-format-security' \
    --prefix=/home/astro/hsc/products/Linux64/external/psfex/rhl \
    --with-fftw-incdir=$FFTW_DIR/include --with-fftw-libdir=$FFTW_DIR/lib \
    --with-atlas-incdir=$ATLAS_DIR/include --with-atlas-libdir=$ATLAS_DIR/lib
make
scons -Q opt=3
eups declare -r . psfex svn -c

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

cd ~/LSST/meas/extensions/psfex
setup -r .
setup -j -r ~/LSST/devenv/sconsUtils
scons -Q opt=3

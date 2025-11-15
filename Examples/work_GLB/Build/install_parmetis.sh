#!/bin/bash
# Install ParMetis in user directory
PARMETIS_INSTALL=/work/home/ac6vfo7a3a/soft/ParMetis

# Dependencies
#   * cmake
#   * Fortran90 compiler
#   * MPI
export CFLAGS=-fPIC
source ~/load_modules.sh
CC=mpicc
CXX=mpicxx

# Check install dir
mkdir -p $PARMETIS_INSTALL
cd $PARMETIS_INSTALL && echo "Current directory: $(pwd)"
[ -d "parmetis-4.0.3" ] && rm -rf parmetis-4.0.3 bin lib include
echo "Start to install ParMetis..."
if [ ! -f "parmetis_4.0.3.orig.tar.gz" ]; then
	wget https://launchpad.net/ubuntu/+archive/primary/+sourcefiles/parmetis/4.0.3-4/parmetis_4.0.3.orig.tar.gz
fi
tar -xvzf parmetis_4.0.3.orig.tar.gz

# Build: Metis
cd parmetis-4.0.3/metis
make config cc=$CC cxx=$CXX prefix=$PARMETIS_INSTALL |& tee metis-make-config.out
make VERBOSE=1 |& tee metis-make.out
make install

# Build: ParMetis
cd ..
make config cc=$CC cxx=$CXX prefix=$PARMETIS_INSTALL |& tee parmetis-make-config.out
make VERBOSE=1 |& tee parmetis-make.out
make install

# Check lib
cd ..
if [ -d "$PARMETIS_INSTALL/lib" ]; then
	echo "ParMetis installed successfully."
else
	echo "ParMetis installation failed."
fi

#!/bin/bash
# Build WW3 for UNST grid with ParMETIS support
if [ "$#" -eq 0 ]; then
	switch=$PWD/switch_PDLIB
	build=build_PDLIB
elif [ "$#" -eq 1 ]; then
	switch=$(realpath "$1")
	bn=$(basename "$switch")
	if [[ "$bn" == switch_* ]]; then
		build="build_${bn#switch_}"
	else
		echo "Error: Invalid switch file name: $bn" >&2
		exit 1
	fi
else
	echo "Error: Usage: $0 [<switch_file>]" >&2
	exit 1
fi
ww3_root=$HOME/caolang/WW3-7.14

if [ ! -f "$switch" ]; then
	echo "Error: Switch file not found: $switch"
	exit 1
fi
if [ ! -d "$ww3_root" ]; then
	echo "Error: WW3 root directory not found: $ww3_root"
	exit 1
fi

# Module load
source ~/load_modules.sh || { echo "ERROR: Failed to load modules"; exit 1; }

# Optionally set compiler and env vars to locate libraries
export CC=mpicc
export FC=mpifort
export NetCDF_ROOT=/public/software/mathlib/netcdf/4.4.1/intel
export ParMETIS_ROOT=$HOME/soft/ParMetis

# build
[ -d "$build" ] && rm -rf "$build"
mkdir "$build" && cd "$build"
cmake "$ww3_root" -DSWITCH="$switch" -DCMAKE_INSTALL_PREFIX=install
make
make install
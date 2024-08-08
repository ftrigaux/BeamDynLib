# How to install standAlone BeamDyn

## Install from this version
You will need
- GCC (gcc and gfortran)
- A version of Lapack and Blas (or OpenBlas)

Simply type
```
make
```
And the executable will be built. Use `make lib` to build the library version, used for Python and C. 
You can then run the simulation in the run directory
If you are using OpenBlas (usually available on cluster as a module to load), you can change this line in the MakeFile
```
#LAPACK_LINK  = -llapack -lblas # Must become ->
LAPACK_LINK  = -lopenblas
```
The specific makefile for the Lucia Tier-0 supercomputer is provided as Makefile.lucia, and uses OpenBLAS, as well as the standard installation directory.

An installation script `install.sh` is included to guide you to the whole installing process on Lucia.

## Install from updated OpenFast sources
### Get the sources, or update the source from the OpenFast repository
To obtain the sources from the openfast repository, copy the following directories:
```
cp    openfast/modules/beamdyn/src/*  BeamDynLib/src/
cp -r openfast/modules/nwtc-libraries BeamDynLib/src/
cp -r openfast/modules/version        BeamDynLib/src/
```
This will update the "base" version of BeamDyn. 
The interface consists of several types and routines included in the files starting with `BeamDynLib`. The C interface is provided in the `CBeamDyn.c` file. 

### Build the project
On Linux and MacOS, if you have gfortran installed, you can simply run
```
make -j 8
```
in your terminal. This will create a `obj` directory with all the compiled objects, and an executable `BeamDyn`.
To make a shared object, that can be used to access the library from C or Python, use 
```
make lib -j 8
```
To build the fortran only program version, use 
```
make fprog
```
To use intel compilers 'ifort' instead, rename `Makefile.intel` to `Makefile` and run `make`. Typically, intel compilers with MKL run faster than with gfortran (~20% in double precision)
### Run the simulation
In a directory (for example: `run/nrel5mw_dynamic`), copy the input files from openfast:
```
scp /openfast/docs/source/user/beamdyn/examples/*.inp .
```
Adapt the parameters, copy the BeamDyn executable, and run the simulation using 
```
./Beamdyn bd_driver_dynamic_nrel_5mw.inp
```
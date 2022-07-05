# How to install standAlone BeamDyn

## Install from this version
Simply type
```
make
```
And the software will be built provided you have gfortran.
You can then run the simulation in the run directory

## Install from updated OpenFast sources
### Get the sources
Fist, obtain the sources from the openfast repository. Copy the following directories:
```
modules/beamdyn/src
modules/nwtc-libraries
modules/version
```
In this framework, they are all included in the `src` directory.

### Build the project
On Linux and MacOS, if you have gfortran installed, you can simply run
```
make
```
in your terminal. This will create a `obj` directory with all the compiled objects, and an executable `BeamDyn`.

### Run the simulation
In a directory (for example: `run/nrel5mw_dybamic`), copy the input files from openfast:
```
scp /openfast/docs/source/user/beamdyn/examples/*.inp .
```
Adapt the parameters, copy the BeamDyn executable, and run the simulation using 
```
./Beamdyn bd_driver_dynamic_nrel_5mw.inp
```
# pHPitz

Computational model developed by Baard Kaasa as part of our CCS research on CO2 Impurities

A thermodynamic model that can calculate whether sulfuric acid will form needs to contain the following main elements: 

- A model to calculate fugacities in the CO2 rich phase 
- A model to calculate reactions in the CO2 rich phase 
- A model to calculate activities in the sulfuric acid 
- Model condensation of H2O and H2SO4 and equilibrium between aqueous phase and CO2 rich phase 

## Build on Windows

Here is one way to build for Windows:

- Install Visual Studio
- Install Intel oneAPI Fortran compiler with MKL libraries
- From Intel oneAPI command prompt

```
mkdir build
cd build
cmake -G "NMake Makefiles" -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_Fortran_FLAGS="/Qmkl" ..
nmake
```

## Build on Linux with gfortran

Here is one way to build it for Linux using gfortran:

```
mkdir build
cd build
cmake -DCMAKE_EXE_LINKER_FLAGS="-llapack -lblas" -DCMAKE_Fortran_FLAGS="-std=legacy -ffree-line-length-none" ..
make
```

for now there might be some manual patches needed in the source code:
- CO2phaseReactions.f90 line 60 - add .NE.0
- Pitzer_Init.f90 in Function PitzTfunc (line 263-315): change signature and ipmlementation from Par to ParArray

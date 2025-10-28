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
mkdir build_windows
cd build_windows
cmake -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_Fortran_FLAGS="/Qmkl" ..
nmake
```

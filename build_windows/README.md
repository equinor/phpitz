Here is one way to build for Windows:

- Install Visual Studio
- Install Intel oneAPI Fortran compiler with MKL libraries

```
cd ./build_windows
cmake -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_Fortran_FLAGS="/Qmkl" ..
nmake
````

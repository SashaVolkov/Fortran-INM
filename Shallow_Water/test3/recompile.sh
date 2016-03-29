#!/bin/bash
rm -rf *.o *.mod *.out *~ 2>/dev/null

~/Fortran/Comands/compc modnet.f90 modfunc.f90 method.f90 Schemes.f90 Shallow_water.f90
~/Fortran/Comands/compo modnet.o modfunc.o method.o Schemes.o Shallow_water.o
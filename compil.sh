#!/bin/sh
mpif90 -c -w spline.f90
mpif90 -c -w mpsi.f90
#mpif90 -Mpreprocess -Droutin -c -w orbit.f90
#mpif90 -cpp -Droutin -c -w orbit.f90
mpif90 -w $1.f90 spline.o mpsi.o -L$ARMNLIB/lib/Linux -lrmn_012 -o $1.exe
#mpif90 $1.f90 spline.o mpsi.o orbit.o -L$ARMNLIB/lib/Linux -lrmnbeta -L/opt/intel/fce/10.1.018/lib -limf -o $1.exe
rm spline.o mpsi.o
if [ -e $1.o ] ; then
   rm $1.o
fi

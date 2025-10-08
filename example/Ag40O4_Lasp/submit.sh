#!/bin/bash
CURDIR=$PWD

export OMP_NUM_THREADS=1
source /home_gsx/compiler/intel/oneapi/setvars.sh
cd $CURDIR
LASP=/home_gsx/users/nsgsx_cdj7/led_ws/llasp/soft/NN_ac3.7.3_intel18_double/Src/lasp


echo -n "start time  " > time
date >> time
#mpirun -np $NP $LASP 
$LASP
echo -n "end   time  " >> time ; date >> time
exit


#!/bin/bash
export AFEPACK_PATH=$HOME/local/src/AFEPack:$HOME/AFEPack
export AFEPACK_TEMPLATE_PATH=$AFEPACK_PATH/template/triangle:$AFEPACK_PATH/template/tetrahedron:$AFEPACK_PATH/template/twin_tetrahedron:$AFEPACK_PATH/template/four_tetrahedron:$AFEPACK_PATH/template/twin_triangle
export LD_LIBRARY_PATH=$HOME/local/lib

#export OMP_NUM_THREADS=$1
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OMP_NESTED=false
./main $1 $2

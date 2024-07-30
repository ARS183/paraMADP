#!/bin/bash

nprocs=$1
filename="MP3BC2_nx289f40_MADP-v3_np${nprocs}.log"
nohup mpirun -np $nprocs \
             --oversubscribe \
             ./helmholtz_2d_solver.out Input/Helmholtz.in </dev/null \
             1>./$filename \
             2>./error.log &

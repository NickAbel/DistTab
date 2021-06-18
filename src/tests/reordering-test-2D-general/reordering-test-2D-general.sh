#!/bin/bash

for PROCS in 1
do
  for N1 in 1 2 4
  do
    for N2 in 1 2 4
    do
      for Q1 in 1 2 4
      do
        for Q2 in 1 2 4
        do
          if [ "$Q1" -gt "$N1" -o "$Q2" -gt "$N2" ]
          then
            continue
          fi
          mpirun -np $PROCS --oversubscribe ../../lib/disttab 3 $N1 $N2 1 2 $Q1 $Q2
        done
      done
    done
  done
done

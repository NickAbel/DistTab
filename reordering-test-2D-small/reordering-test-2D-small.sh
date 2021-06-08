#!/bin/bash

for PROCS in 1 2
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
          mpirun -np $PROCS --oversubscribe ../disttab 3 $N1 $N2 1 2 $Q1 $Q2 2 | sort -g >| ./np-$PROCS\_n-$N1-$N2\_q-$Q1-$Q2.out
          diff --unified=0 ./np-$PROCS\_n-$N1-$N2\_q-$Q1-$Q2.out ./gold/np-$PROCS\_n-$N1-$N2\_q-$Q1-$Q2.gold
          rm ./np-$PROCS\_n-$N1-$N2\_q-$Q1-$Q2.out 
        done
      done
    done
  done
done

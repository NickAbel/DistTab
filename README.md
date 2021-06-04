# DistTab
Distributed Tabulation of Spray Flamelet Lookup Tables

## Instructions
* Compile
```shell
mpif90 disttab.f90 -o disttab
```

* Run
```shell
mpirun -np 2 ./disttab
```

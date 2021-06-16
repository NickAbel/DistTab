# DistTab
Distributed Tabulation of Spray Flamelet Lookup Tables

## Instructions
* Dependencies
MPI (I use OpenMPI)
* Compile
```shell
cd src/lib
mpifort disttab_table.f90 disttab.f90 -o disttab
```

* Run Test (Command Runstring To Be Explained Later)
```shell
cd src/tests/reordering-test-2D-small/
./reordering-test-2D-small.sh
```

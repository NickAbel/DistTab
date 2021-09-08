<img src="disttab-logo.png" align="left" />

# DistTab

Distributed Tabulation of Flamelet Tables

## Table of Contents

* [About](#about)
* [Functionality](#functionality)
* [Prerequisites](#prerequisites)
  * [Install](#install)
  * [Tree](#tree)
* [Documentation](#documentation)
* [Todo](#todo)
* [Alya](#alya)
* [Authors](#authors)

## About

DistTab is a library being developed for use with combustion codes that use spray flame lookup tables. Currently, simulations across multiple NUMA nodes store the entire flame lookup table on each node, which limits the size of the lookup table. DistTab distributes the lookup table across nodes, partitions the table in a logical manner and uses MPI one-sided communication to retrieve lookup table values and stores the partitions containing the lookup values on the local node using a FIFO queue.

## Functionality

* Lookup table storage
  * Loading/reading in table
  * Reorganizing table into partitions based on locality in thermochemical state space
* MPI get values from other nodes
  * Based on partition, stored in FIFO queue locally
* And more...

## Prerequisites

* Fortran compiler, such as gfortran
* MPI library, such as OpenMPI

### Install

To compile and run the tests

```sh
$ git clone https://github.com/NickAbel/DistTab.git
$ cd DistTab
$ mkdir build
$ cd build
$ cmake ..
$ make
$ cd ../bin
$ ./disttab
```


### Tree

```text
├── cmake
│   └── Modules
├── scripts
├── src
│   ├── old
│   ├── table
│   ├── tests
│   └── third_party
│       └── (...)
└── tables
```

## Documentation

Documentation may be generated using Doxygen with the (`Doxyfile`) in the base directory.

One can generate the documentation by running the following in the base directory:

```sh
$ doxygen Doxyfile
```

The documentation is output in the directory (`doc`).

## Todo

- [x] Relevant documentation/comments for the completed items below
- [x] Simple test for two-dimensional partitioning, divisible partition sizes
- [x] Two-dimensional partitioning, partition sizes perfectly divisible by table sizes
- [x] Fulsome test for two-dimensional partitioning, divisible partition sizes
- [x] Simple test for two-dimensional partitioning, nondivisible partition sizes
- [x] Two-dimensional partitioning, with zero-padding for nondivisible partition sizes
- [x] Fulsome test for two-dimensional partitioning, nondivisible partition sizes
- [x] Simple test for n-dimensional partitioning, divisible partition sizes
- [x] n-dimensional partitioning, partition sizes perfectly divisible by table sizes
- [x] Fulsome test for n-dimensional partitioning, divisible partition sizes
- [x] Simple test for n-dimensional partitioning, nondivisible partition sizes
- [x] n-dimensional partitioning, with zero-padding for nondivisible partition sizes
- [x] Fulsome test for n-dimensional partitioning, nondivisible partition sizes
- [x] When number of state variables > 1, the state variables per-value must be contiguous in memory
- [x] Given spatial coordinates in integer, local, or global form, return state variables from partitioned table
- [x] Given spatial coordinates in integer, local, or global form, return state variable value cloud from partitioned table
- [x] Given control variable coordinates in [0,1]x...x[0,1], return state variable value cloud from partitioned table
- [x] Repartitioning; accept 'old' partition dimensions as option instead of assuming 'Alya format'
- [x] Don't rely on the sort command to generate input for tests
- [ ] Make output verbiage from tests easier to understand, should see at a glance all PASS/FAIL/DIFFs.
- [ ] Re-introduce MPI calls and split table among ranks in simple manner
- [ ] Begin working on local pile
- [ ] Explain new functionalities in README.md

## Alya

Eventually, want to use library with combustion codes such as Alya, OpenFOAM, etc.

## Authors

* [**Nicholas Abel**](https://github.com/NickAbel) 

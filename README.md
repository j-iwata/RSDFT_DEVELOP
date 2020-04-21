# RSDFT_DEVELOP

## Installation
Download or clone the branch SX-AuroraTSUBASA, and
```
$ cd RSDFT_DEVELOP/src
$ make
```
then you get the execution rsdft.x

## Usage
You need some input files (fort.1, fort.970, and pseudopotentials). Please see the documents in RSDFT repogitry for details
(https://github.com/j-iwata/RSDFT/tree/master/doc). Example of 1-core execution is
```
$ mpirun -np 1 rsdft.x
```
Execution with 8 cores of a specific VE (VE #3 for example),
```
$ mpirun -np 8 -ve 3 rsdft.x
```
Example of the execution with multiple VEs (VE #0 and #1) is
```
$ mpirun -np 16 -ve 0-1 rsdft.x
```


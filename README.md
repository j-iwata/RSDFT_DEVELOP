# RS-CPMD

## Installation
Download or clone the branch RSCPMD, and
```
$ cd RSDFT_DEVELOP/src
$ cp make.inc.org make.inc
```
then you should edit make.inc to suite your compiler, libraries, etc.

The default FFT library used in RS-CPMD is FFTE (http://www.ffte.jp). Unfortunately the library has a restriction that the number of grid points in each direction (it is specified by NGRID keyword in the RSDFT input file) has to be (2^p)\*(3^p)\*(5^r).

FFTW is also possible to use. If you don't want to avoid the above restriction, we recommend to compile RS-CPMD with FFTW.

After editing the make.inc, just type following command.
```
$ make
```
then you get the execution file rsdft.x

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

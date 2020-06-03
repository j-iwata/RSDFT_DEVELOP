# RS-CPMD

Gamma-point sampling mode is only available for CPMD calculations.

## Installation
Download or clone the branch RSCPMD, and type the following commands
```
$ cd RSDFT_DEVELOP/src
$ cp make.inc.org make.inc
```
Next you should edit make.inc to suite your compiler and libraries.

The default FFT library used in RS-CPMD is FFTE (http://www.ffte.jp). Unfortunately the library has a restriction that the number of grid points in each direction (it is specified by NGRID keyword in the RSDFT input file) has to be (2^p)\*(3^p)\*(5^r).

FFTW is also possible to use. If you want to avoid the above restriction, we recommend to compile RS-CPMD with FFTW.

After editing the make.inc, type the following command
```
$ make
```
then you get the execution file rsdft.x

## Usage
You need several input files (fort.1, fort.970, cpmd_var.dat and pseudopotentials). Please see the documents in RSDFT repogitry for details
(https://github.com/j-iwata/RSDFT/tree/master/doc). Typical example of 1-core execution is
```
$ mpirun -np 1 rsdft.x
```
Execution with 8 cores is
```
$ mpirun -np 8 rsdft.x
```

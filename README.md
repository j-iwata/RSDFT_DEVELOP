# RSDFT

## NOTICE
The format of the input parameters for atomic structure optimization have been changed. If you use the old one, for example,

```
ATOMOPT1  50  6  5                     / ncycle  most  nrfr
ATOMOPT2  0.5d0  1.d-10  5.d-4  1.d-1  / okatom  eeps  feps  decr
ATOMOPT3  100                          / diter_opt
```

please change as follows

```
ATOMOPT1  50  5.d-4  .f.        / ncycle  feps  continue_flag
ATOMOPT2  6  5  100             / most  nrfr  diter_opt
ATOMOPT3  0.5d0  1.d-1  1.d-10  / diter_opt
```

If you use CG optimization (*SWOPT = 1* or *2*), *continue_flag* is omitted.

You can also omit the other aprameters, for example,

```
ATOMOPT1  50  5.d-4  / ncycle  feps
```

and in this case the omitted parameters are set by default values.
 

## ITO (Kyusyu university)

### Compile

```
module load intel fftw
cd RSDFT/src
cp make.inc.ITO.complex16 make.inc
make
```

after a few minutes, you get *__rsdft.x__* which is the execution file for the calculation with double-complex wavefunctions. The double-complex rsdft.x is necessary for the calculaiotns with BZ sampling for k/=(0,0,0) and band-structure calculations.

If your k-point sampling is only &Gamma; point, you can use *rsdft.x* for the calculation with double-real wavefunctions. The compilation procedure is

```
module load intel fftw
cd RSDFT/src
cp make.inc.ITO.real8 make.inc
make
```

With the doube-real *rsdft.x*, you can also perform RS-CPMD calculations.

### Samples
Please see *samples_for_supercomputers/ITO_KYUSHU-U/* .



## FUGAKU (RIKEN)



## FLOW-FX (Nagoya university)

The complex16 version of *rsdft.x* can be compiled as follows:

```
module load fftw
cd RSDFT/src
cp make.inc.FLOW-FX.complex16 make.inc
make
```

and the real8 version is

```
module load fftw
cd RSDFT/src
cp make.inc.FLOW-FX.real8 make.inc
make
```

### Samples
Please see *samples_for_supercomputers/FLOW-FX_NAGOYA-U/* .



## REEDBUSH-H (The university of Tokyo)

The complex16 version of *rsdft.x* can be compiled as follows:

```
cd RSDFT/src
cp make.inc.REEDBUSH-H.complex16 make.inc
make
```

and the real8 version is

```
cd RSDFT/src
cp make.inc.REEDBUSH-H.real8 make.inc
make
```

### Samples
Please see *samples_for_supercomputers/REEDBUSH-H_U-TOKYO/* .

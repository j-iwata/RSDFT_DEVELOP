# RSDFT

## ITO (Kyusyu university)

### Compile

```
module load intel fftw
cd RSDFT/src
cp make.inc.ITO.complex16 make.inc
make
```

after a few minutes, you get *__rsdft.x__* which is the execution file for the calculation with double-complex wavefunctions. The double-complex rsdft.x is necessary for the calculaiotns with BZ sampling for k/=(0,0,0) and band-structure calculations.

If your k-point sampling is only &Gamma; k point, you can use *rsdft.x* for the calculation with double-real wavefunctions. The compilation procedure is

```
module load intel fftw
cd RSDFT/src
cp make.inc.ITO.real8 make.inc
make
```

With the doube-real *rsdft.x*, you can also perform RS-CPMD calculations.

### Samples
Please see *samples_for_supercomputers/ITO_KYUSHU-U/* .


## Fugaku (RIKEN)


## Furou (Nagoya university)

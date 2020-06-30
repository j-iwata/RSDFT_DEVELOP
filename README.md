# RSDFT

## ITO (kyushu-u)

### Compile

```
module load intel fftw
cd RSDFT/src
cp make.inc.ITO.complex16 make.inc
make
```

after a few minutes, you get *__rsdft.x__* which is the execution file for the calculation with double-complex wavefunctions. The double-complex rsdft.x is necessary for the calculaiotns with BZ sampling for k/=(0,0,0) and band-structure calculations.

If you sample only \gamma k point, you can use *__rsdft.x__* for the calculation with double-real wavefunctions. The compilation procedure is

```
module load intel fftw
cd RSDFT/src
cp make.inc.ITO.real8 make.inc
make
```

With the doube-real rsdft.x, you can also perform RS-CPMD calculations.

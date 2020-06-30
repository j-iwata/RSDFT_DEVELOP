# RSDFT

## ITO (kyushu-u)

### Compile

```
module load intel fftw
cd RSDFT/src
cp make.inc.ITO.complex16 make.inc
make
```

after a few minutes, you get __rsdft.x__ which is the execution file for the calculations with double-complex wavefunctions. The double-complex rsdft.x is necessary to perform band-structure calculations.

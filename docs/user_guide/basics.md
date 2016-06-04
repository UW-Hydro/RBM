# Basics

## Fortran 
RBM is run with Fortran. As such it is run from the terminal window. If you have not used a computer through the terminal, there are a plethora of tutorials online for an introduction. Since Fortran is a compiled language, it must be compiled before it is run. Code is developed as ".f90" (i.e. Begin.f90) files.  Those files are then "compiled" with a fortran compiler (see the [support page] (INSERT) for details on the compiler). And before code is compiled, previously compiled code must be "cleaned" from the directory. In RBM, we have developed a "makefile" that allows these steps to be done easily by first entering 

```shell
make clean
```
 which cleans out previously compiled code, then enter

```shell
make
```
which compiles all .f90 files in that directory.  

## Other


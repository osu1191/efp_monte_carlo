# C compiler
CC= gcc

# Fortran compiler
FC= gfortran

# install prefix
PREFIX= /depot/lslipche/data/skp/tstgit/efp_monte_carlo

# fragment library path
FRAGLIB= $(PREFIX)/share/libefp/fraglib

# additional link libraries
MYLIBS= -lblas -llapack -lgfortran

# additional linker flags
MYLDFLAGS=

# additional C flags
MYCFLAGS= -std=c99 -O3 -g -fopenmp

# additional Fortran flags
MYFFLAGS= -g
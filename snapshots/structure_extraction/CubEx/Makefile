# -*- Makefile -*-
# =======================================================================
#
#				Makefile for CubEx-x.x
#
# Author: Sebastiano Cantalupo
#
# =======================================================================


#--- CHANGE THE FOLLOWING IF NECESSARY:
FC=gfortran
CFITSIO_LIBPATH= -L/usr/lib/x86_64-linux-gnu -lcfitsio
BINDIR=/scratch/jgacon/CubEx-1.7/bin
#
#-------------------------------------

#
# You should not change anything below this line, unless you want to change compiler and architecture specific flags
# 
ifeq ($(FC),ifort)
 ifeq ($(shell uname),Darwin)
  FFLAGS=-O2 -assume bscc -vec-report0 -heap-arrays 1600 -Wl,-stack_size,0x10000000
 else
  FFLAGS=-O2   # ifort has been only tested on Mac, you may need to add arch specific flags here
 endif
else 
 FFLAGS=-O2 -ffree-line-length-none -static-libgfortran #-fcheck=all -fbacktrace
endif 
DFLAGS= -O0 -g -debug all -check all -traceback -assume bscc -vec-report0 -heap-arrays 1600 -Wl,-stack_size,0x10000000
CFLAGS= -I. $(CFITSIO_LIBPATH) 
LINKER= $(FC)
#

# sources:
EXEC= CubEx
MOD=Globalmodule.f90 StatLib.f90
SRC=main.f90 ReadParameterFile.f90 CreateCube.f90 Extract.f90 WriteCheckCube.f90 \
    WriteCatalogue.f90 Photometry.f90 ReadInputFile.f90 Filter.f90 Associate.f90 \
    MergeAssoc.f90 Deblend.f90 UseIdCube.f90
#
MODO= $(MOD:.f90=.o)
MODB= $(MOD:.f90=.mod)
OBJ= $(SRC:.f90=.o) 

#
#


default: $(BINDIR)/$(EXEC)

$(BINDIR)/$(EXEC): $(MODO) $(OBJ) Makefile
	$(LINKER) $(FFLAGS) $(CFLAGS) $(MODO) $(OBJ) -o $(BINDIR)/$(EXEC) $(CFITSIO_LIBPATH)

debug:  
	make -f Makefile FFLAGS="$(DFLAGS)"

show: 
	echo $(MODO) $(OBJ) 

clean:
	rm -f -v *.o *.mod

tidy:	clean
	$(RM) $(BINDIR)/$(EXEC)

new:	clean default

%.o: %.f90 
	$(FC) $(CFLAGS) $(FFLAGS) -c $<


 




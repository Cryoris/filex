# -*- Makefile -*-


#--- CHANGE THE FOLLOWING IF NECESSARY:
#
FC=gfortran
CFITSIO_LIBPATH=-L/usr/lib/x86_64-linux-gnu -lcfitsio
BINDIR=$(HOME)/bin
#-------------------------------------

#
#
# You should not change anything below this line, unless you want to change compiler and architecture specific flags
#
ifeq ($(FC),ifort)
 ifeq ($(shell uname),Darwin)
  FFLAGS=-O2 -assume bscc -vec-report0 -heap-arrays 1200 -Wl,-stack_size,0x10000000 -openmp
 else
  FFLAGS=-O2 -openmp # ifort has been only tested on Mac, you may need to add arch specific flags here
 endif
else 
 FFLAGS=-O2 -ffree-line-length-none -static-libgfortran -fopenmp #-fcheck=all -fbacktrace
endif 
DFLAGS= -O0 -g -debug all -check all -traceback -assume bscc -vec-report0 
CFLAGS= -I. -I../ $(CFITSIO_LIBPATH) #$(DFLAGS)
LINKER= $(FC)
#

# sources:
EXEC= CubeSel Cube2Im pix2wcs Cube2Spc CubeBKGSub
#
MODO= $(MOD:.f90=.o)
MODB= $(MOD:.f90=.mod)
OBJ= $(SRC:.f90=.o) 

#
#
all:    ../StatLib.o CubeLib.o CubeSel.o Cube2Im.o Cube2Spc.o CubeBKGSub.o pix2wcs.o CubeCombine.o CubeArit.o CubePSFSub.o CubeStat.o CubeFix.o CubeSharp.o PixTableMask.o CubeAdd2Mask.o GetImageOffsets.o 
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o CubeSel.o -o $(BINDIR)/CubeSel $(CFITSIO_LIBPATH)
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o Cube2Im.o -o $(BINDIR)/Cube2Im $(CFITSIO_LIBPATH)
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o Cube2Spc.o -o $(BINDIR)/Cube2Spc $(CFITSIO_LIBPATH)
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o CubeBKGSub.o -o $(BINDIR)/CubeBKGSub $(CFITSIO_LIBPATH)
	$(LINKER) $(FFLAGS) $(CFLAGS) pix2wcs.o -o $(BINDIR)/pix2wcs $(CFITSIO_LIBPATH)
	$(LINKER) $(FFLAGS) $(CFLAGS) -fopenmp ../StatLib.o CubeLib.o CubeCombine.o -o $(BINDIR)/CubeCombine $(CFITSIO_LIBPATH)
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o CubeArit.o -o $(BINDIR)/CubeArit $(CFITSIO_LIBPATH)
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o CubePSFSub.o -o $(BINDIR)/CubePSFSub $(CFITSIO_LIBPATH)
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o CubeStat.o -o $(BINDIR)/CubeStat $(CFITSIO_LIBPATH)
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o CubeFix.o -o $(BINDIR)/CubeFix $(CFITSIO_LIBPATH)
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o PixTableMask.o -o $(BINDIR)/PixTableMask $(CFITSIO_LIBPATH)
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o CubeAdd2Mask.o -o $(BINDIR)/CubeAdd2Mask $(CFITSIO_LIBPATH)
	$(LINKER) $(FFLAGS) $(CFLAGS) -openmp ../StatLib.o CubeLib.o CubeSharp.o -o $(BINDIR)/CubeSharp $(CFITSIO_LIBPATH)
	$(LINKER) $(FFLAGS) $(CFLAGS) GetImageOffsets.o -o $(BINDIR)/GetImageOffsets.x $(CFITSIO_LIBPATH)

cubesel: ../StatLib.o CubeLib.o CubeSel.o
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o CubeSel.o -o $(BINDIR)/CubeSel $(CFITSIO_LIBPATH)

cube2im: ../StatLib.o CubeLib.o Cube2Im.o
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o Cube2Im.o -o $(BINDIR)/Cube2Im $(CFITSIO_LIBPATH)

cube2spc: ../StatLib.o CubeLib.o Cube2Spc.o
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o Cube2Spc.o -o $(BINDIR)/Cube2Spc $(CFITSIO_LIBPATH)

cubebkgsub: ../StatLib.o CubeLib.o CubeBKGSub.o
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o CubeBKGSub.o -o $(BINDIR)/CubeBKGSub $(CFITSIO_LIBPATH)

pix2wcs: pix2wcs.o
	$(LINKER) $(FFLAGS) $(CFLAGS) pix2wcs.o -o $(BINDIR)/pix2wcs $(CFITSIO_LIBPATH)

cubeheader: CubeHeader.o
	$(LINKER) $(FFLAGS) $(CFLAGS) CubeHeader.o -o $(BINDIR)/CubeHeader $(CFITSIO_LIBPATH)

getimageoffsets: GetImageOffsets.o
	$(LINKER) $(FFLAGS) $(CFLAGS) GetImageOffsets.o -o $(BINDIR)/GetImageOffsets.x $(CFITSIO_LIBPATH)

cubecomb: ../StatLib.o CubeLib.o CubeCombine.o
	$(LINKER) $(FFLAGS) $(CFLAGS) -fopenmp ../StatLib.o CubeLib.o CubeCombine.o -o $(BINDIR)/CubeCombine $(CFITSIO_LIBPATH)

cubestat: ../StatLib.o CubeLib.o CubeStat.o
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o CubeStat.o -o $(BINDIR)/CubeStat $(CFITSIO_LIBPATH)

cubearit: ../StatLib.o CubeLib.o CubeArit.o
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o CubeArit.o -o $(BINDIR)/CubeArit $(CFITSIO_LIBPATH)

cubepsfsub:  ../StatLib.o CubeLib.o CubePSFSub.o
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o CubePSFSub.o -o $(BINDIR)/CubePSFSub $(CFITSIO_LIBPATH)

cubefix: ../StatLib.o CubeLib.o CubeFix.o
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o CubeFix.o -o $(BINDIR)/CubeFix $(CFITSIO_LIBPATH)

cubesharp: ../StatLib.o CubeLib.o CubeSharp.o
	$(LINKER) $(FFLAGS) $(CFLAGS) -openmp ../StatLib.o CubeLib.o CubeSharp.o -o $(BINDIR)/CubeSharp $(CFITSIO_LIBPATH)

cubereplace: ../StatLib.o CubeLib.o CubeReplace.o
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o CubeReplace.o -o $(BINDIR)/CubeReplace $(CFITSIO_LIBPATH)

cuberesample: ../StatLib.o CubeLib.o CubeResample.o
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o CubeResample.o -o $(BINDIR)/CubeResample $(CFITSIO_LIBPATH)

cubereshape: ../StatLib.o CubeLib.o CubeReshape.o
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o CubeReshape.o -o $(BINDIR)/CubeReshape $(CFITSIO_LIBPATH)

cubetrace: ../StatLib.o CubeLib.o CubeTrace.o
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o CubeTrace.o -o $(BINDIR)/CubeTrace $(CFITSIO_LIBPATH)

cubetest: ../StatLib.o CubeLib.o CubeTest.o
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o CubeTest.o -o $(BINDIR)/CubeTest $(CFITSIO_LIBPATH)

pixtablemask: ../StatLib.o CubeLib.o PixTableMask.o
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o PixTableMask.o -o $(BINDIR)/PixTableMask $(CFITSIO_LIBPATH)

cubeadd2mask: ../StatLib.o CubeLib.o CubeAdd2Mask.o
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o CubeAdd2Mask.o -o $(BINDIR)/CubeAdd2Mask $(CFITSIO_LIBPATH)

cubequickspc: ../StatLib.o CubeLib.o CubeQuickSpc.o
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o CubeQuickSpc.o -o $(BINDIR)/CubeQuickSpc $(CFITSIO_LIBPATH)

sbprof: ../StatLib.o CubeLib.o SBProf.o
	$(LINKER) $(FFLAGS) $(CFLAGS) ../StatLib.o CubeLib.o SBProf.o -o $(BINDIR)/SBProf $(CFITSIO_LIBPATH)

clean:
	rm -f -v *.o *.mod

tidy:	clean
	$(RM) $(BINDIR)/$(EXEC)

new:	clean default

%.o: %.f90 
	$(FC) $(CFLAGS) $(FFLAGS) -c $<


 




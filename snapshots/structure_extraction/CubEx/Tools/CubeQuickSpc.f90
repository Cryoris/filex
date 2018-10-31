!-----------------------------------------
!
! performs 1d spectrum extraction (and plotting) of datacube spectra selected 
! interactively from a white-light image displayed by ds9
!
! Author: SC
! Last Mod: Aug 27, 2014
!
!------------------------------------------


PROGRAM main

  USE StatLib
  USE CubeLib
  IMPLICIT NONE
  CHARACTER(len=350) :: InpFile, IdCube, OutputImage, VarCube, ImType, Output1d, masktype, outEPS
  INTEGER(kind=4)    :: SelId, k, proj, i, j, CrossProj, zstart, zend, nend, nint, smooth, xmin, xmax, ymin, ymax, nsel
  LOGICAL :: plot, interactive
  REAL :: xlist(1000), ylist(1000), apsize

  Verbosity=2

  apsize=5

  !.. collect input parameters
  CALL ReadCommandLine

  !..read datacube
  CALL ReadCube(InpFile)

  !..get coords from DS9
  CALL GetCoordinatesFromDS9

  DO i=1,nsel

     xmin=NINT(xlist(i)-apsize*0.5)
     xmax=NINT(xmin+apsize)
     ymin=NINT(ylist(i)-apsize*0.5)
     ymax=NINT(ymin+apsize)

     OPEN(1,file=Output1d,action="write")
     DO k=1,SIZE(Cube,DIM=3)
        WRITE(1,*) k, (k-WCS(3))*WCS(11)+WCS(6), &
             SUM(Cube(xmin:xmax,ymin:ymax,k),MASK=Cube(xmin:xmax,ymin:ymax,k)/=UNDEF)
     END DO
     CLOSE(1)

     IF(plot) CALL showplot(TRIM(Output1d))
     IF(TRIM(outEPS)/="??") CALL saveEPS(TRIM(Output1d))

  END DO

CONTAINS

SUBROUTINE ReadCommandLine

  IMPLICIT NONE
  CHARACTER(len=300) :: ExeName, fname, string, ds, opt, arg
  INTEGER :: narg, iarg, i


!..default
  InpFile="??"
  IdCube="??"
  Output1d="tmp.out"
  VarCube="??"
  SelId=0
  CrossProj=0
  ImType="Flux"
  plot=.true.
  masktype="2d"
  interactive=.false.
  smooth=0
  outEPS="??"

  CALL GetArg(0,ExeName)
  narg=iargc()
  IF(narg>0) THEN
     CALL GetArg(1,fname)
  ELSE
     print *, " "
     WRITE(*,'(2a)')"        Cube2Spc (part of CubEx package)   "
     WRITE(*,'(a)') " "
     WRITE(*,'(a)') " by Sebastiano Cantalupo (cantal@ucolick.org)"
     print *, " "
     WRITE(*,'(a)') "usage: Cube2Spc.x -cube <name> -idcube <name> -id <val> [-option <val>]"
     WRITE(*,'(a)') " "
     WRITE(*,'(a)') "  options:"
     WRITE(*,'(a)') "  -cube          <name>             : datacube name (NO DEFAULT)"
     WRITE(*,'(a)') "  -idcube        <name>             : name of the 3d mask with object ids (NO DEFAULT)"
     WRITE(*,'(a)') "  -id            <name>             : object id to plot (NO DEFAULT)"
     WRITE(*,'(a)') "  -outname        <string>          : if a name is provided, save the 1d spectrum for the selected object in a ASCII file"
     WRITE(*,'(a)') "  -outEPS        <string>           : if a name is provided, save the 1d spectrum for the selected object in a EPS file (requires SM)"
     WRITE(*,'(a)') "  -show           <bol>             : if .true. (default), shows the 1d spectrum using python"
     WRITE(*,'(a)') "  -masktype       '2d'/'3d'         : if '2d' project the mask on XY before applying it (default)"
     WRITE(*,'(a)') "  -interactive    <bol>             : if .true. interactively asks for new spectra to plot without re-reading cube (default=.false.)"
     WRITE(*,'(a)') "  -smoothr        <bol>             : (gaussian) smoothing radius for the plotted spectrum (default=0, no smoothing)"
     STOP
  END IF

  !..read from command line
  DO i=1,narg,2
     CALL getarg(i,opt)
     CALL getarg(i+1,arg)
     SELECT CASE(TRIM(opt))
     CASE('-cube')       ; READ(arg,'(a)') InpFile
     CASE('-masktype')   ; READ(arg,'(a)') masktype
     CASE('-outname')       ; READ(arg,'(a)') Output1d
     CASE('-idcube')        ; READ(arg,'(a)') IdCube
     CASE('-varcube')       ; READ(arg,'(a)') VarCube
     CASE('-id')            ; READ(arg,*) SelId
     CASE('-interactive')   ; READ(arg,*) interactive
     CASE('-show')          ; READ(arg,*) plot
     CASE('-smoothr')       ; READ(arg,*) smooth
     CASE('-outEPS')        ; READ(arg,'(a)') outEPS
     CASE default
        print *, "command line argument ",TRIM(opt), " not recognized!"
        STOP
     END SELECT
  END DO

!..perform few checks
  IF(TRIM(InpFile)=="??") THEN
     print *, "please provide the input datacube with the -cube option!"
     STOP
  END IF

END SUBROUTINE ReadCommandLine

!-----------------------------------------------------

SUBROUTINE GetCoordinatesFromDS9

    IMPLICIT NONE
    CHARACTER(len=500) :: cmd, imagename, cubebase, datacube
    INTEGER :: is, ierr
    REAL    :: x, y
    LOGICAL :: ex

    datacube=InpFile

    !..creates a ds9 analysis file
    OPEN(11,file=".analysis.ds9",action="write")
    WRITE(11,*) "Print Coordinates"
    WRITE(11,*) "*.fits"
    WRITE(11,*) "bind s"
    WRITE(11,*) 'echo "$x $y" >> .coords.txt | $null'
    CLOSE(11)

    !..remove previous coords.txt file
    cmd="rm -f .coords.txt"
    CALL SYSTEM(cmd)

    !..get a image file name
    is=INDEX(TRIM(datacube),".fits")
    IF(is==0) is=LEN_TRIM(datacube)+1
    cubebase=TRIM(datacube(1:is-1))
    imagename=TRIM(cubebase)//".IM.fits"  
    INQUIRE(File=imagename, Exist=ex)
    IF(.not.ex) THEN
       !..create an image with Cube2Im
       cmd="Cube2Im "//TRIM(datacube)
       CALL SYSTEM(cmd)
    END IF

    print *, "type 's' to select the central position of the apertures on ds9 window. Close the ds9 to exit"

    !..run ds9
    cmd='ds9 '//TRIM(imagename)//" -analysis load .analysis.ds9 -scale mode zscale -cmap b -zoom to fit"
    CALL SYSTEM(cmd)

    !..retrieve coordinates 
    OPEN(11,file=".coords.txt",action="read")
    nsel=0
    DO
       READ(11,*,iostat=ierr) x, y
       IF(ierr/=0) EXIT
       nsel=nsel+1
       xlist(nsel)=x; ylist(nsel)=y
    END DO
    CLOSE(12)
 

  END SUBROUTINE GetCoordinatesFromDS9


SUBROUTINE showplot(fname)

  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN) :: fname

  !..write a python script:

  OPEN(99,file="plot.py",action="write")
  WRITE(99,'(a)') "import matplotlib.pyplot as plt"
  WRITE(99,'(a)') "import numpy as np"
  WRITE(99,'(a)') "from scipy import signal"
  WRITE(99,'(3a)') "f2 = open('",TRIM(fname),"', 'r')"
  WRITE(99,'(3a)') "lines = f2.readlines()"
  WRITE(99,'(3a)') "x1 = []"
  WRITE(99,'(3a)') "y1 = []"  
  WRITE(99,'(a)') "for line in lines:"
  WRITE(99,'(a)') "   p = line.split()"
  WRITE(99,'(a)') "   x1.append(float(p[1]))"
  WRITE(99,'(a)') "   y1.append(float(p[2]))"
  WRITE(99,'(a)') " "  
  WRITE(99,'(a)') "xv = np.array(x1)"
  WRITE(99,'(a)') "yv = np.array(y1)"
  IF(smooth>0) THEN
     WRITE(99,'(a,i2,a)') "w = signal.gaussian(20,std=",smooth,")"
     WRITE(99,'(a)') "yv = signal.convolve(yv, w, mode='same')"
  END IF
  WRITE(99,'(a)') "plt.plot(xv, yv)"
  WRITE(99,'(a)') "plt.show()"
  CLOSE(99)

  CALL SYSTEM("python plot.py")

END SUBROUTINE showplot

SUBROUTINE saveEPS(fname)

  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN) :: fname

  !..write a SM script:
  OPEN(99,file="plot.sm",action="write")
  WRITE(99,'(a)') "define TeX_Strings 1"
  WRITE(99,'(2a)') "dev postencap ",TRIM(outEPS)
  WRITE(99,'(2a)') "data ",TRIM(fname)
  WRITE(99,'(a)') "read {x 2 y 3}"
  IF(smooth>0) THEN
     WRITE(99,'(a,i2)') "smooth y yy",2*smooth
     WRITE(99,'(a)') "set y=yy"
  END IF
  WRITE(99,'(a)') "lim x yy"
  WRITE(99,'(a)') "expand 1.3"
  WRITE(99,'(a)') "lweight 3"
  WRITE(99,'(a)') "box"
  WRITE(99,'(a)') "xlabel \AA"
  WRITE(99,'(a)') "ylabel flux"
  WRITE(99,'(a)') "connect x yy"
  WRITE(99,'(a)') "dev x11"
  CLOSE(99)

  !..run the script
  CALL SYSTEM('sm < plot.sm')

  print *, "1d spectrum saved to: ", TRIM(outEPS)

END SUBROUTINE saveEPS


END PROGRAM main

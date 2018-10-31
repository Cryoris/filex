!-----------------------------------------
!
! performs 1d spectrum extraction (and plotting) of datacube extracted objects
! defined by the 3d mask provided by CubEx (CheckCubeType = "Object_Id" or "Object_Id_Assoc")
!
! Author: SC
! Last Mod: Nov 3, 2014
!
!------------------------------------------


PROGRAM main

  USE StatLib
  USE CubeLib
  IMPLICIT NONE
  CHARACTER(len=350) :: InpFile, IdCube, OutputImage, VarCube, ImType, Output1d, masktype, outEPS
  INTEGER(kind=4)    :: SelId, k, proj, i, j, CrossProj, zstart, zend, nend, nint, smooth
  REAL(kind=4), ALLOCATABLE :: CheckCube(:,:,:)
  INTEGER(kind=4), ALLOCATABLE :: mask2d(:,:)
  REAL(kind=8) :: CheckCubeWCS(SIZE(WCS))
  LOGICAL :: plot, interactive

  Verbosity=2

  !.. collect input parameters
  CALL ReadCommandLine

  !..read checkcube if requested
  IF(TRIM(IdCube)/="??") THEN
     CALL ReadCube(IdCube)
     ALLOCATE(CheckCube(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2),SIZE(Cube,DIM=3)))
     CheckCube=Cube
  END IF

  !..save WCS values for later
  CheckCubeWCS=WCS
  
  !..read checkcube if requested
  IF(TRIM(VarCube)/="??") THEN
     CALL ReadCube(VarCube)
     ALLOCATE(Var(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2),SIZE(Cube,DIM=3))) !..NB: Var is a global array
     Var=Cube
  END IF

  !..read datacube
  CALL ReadCube(InpFile)

  !..if Cube and CheckCube dimensions are different and if WCS are present, get initial pixel value of CheckCube in Cube
  IF(SIZE(Cube,DIM=3)>SIZE(CheckCube,DIM=3)) THEN
     IF(ALL(CheckCubeWCS==0.d0).and.ALL(WCS==0.d0)) STOP "CheckCube and Cube have different z-dimension and WCS were not found to correct for that!"
     zstart=INT((CheckCubeWCS(6)-WCS(6))/WCS(11)+1)
     zend=zstart+SIZE(CheckCube,DIM=3)-1
  ELSEIF(SIZE(Cube,DIM=3)<SIZE(CheckCube,DIM=3)) THEN
     STOP "CheckCube is larger than datacube!"
  ELSE
     zstart=1
     zend=SIZE(Cube,DIM=3)
  END IF

  IF(TRIM(masktype)=="2d") ALLOCATE(mask2d(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2)))


  IF(interactive) THEN
     nend=1000
  ELSE
     nend=1
  END IF

  nint=0

  int_loop:DO 

     nint=nint+1

     OPEN(1,file=Output1d,action="write")

     IF(TRIM(masktype)=="3d") THEN

        DO k=1,SIZE(CheckCube,DIM=3)
           WRITE(1,*) k, (k-CheckCubeWCS(3))*CheckCubeWCS(11)+CheckCubeWCS(6), &
                SUM(Cube(:,:,k+zstart-1),MASK=CheckCube(:,:,k)==REAL(SelId).and.Cube(:,:,k+zstart-1)/=UNDEF)
        END DO
        CLOSE(1)

     ELSEIF(TRIM(masktype)=="2d") THEN

        !..project mask to 2d
        mask2d=0
        IF(SelId/=0) THEN  !..object spectrum
           mask2d(:,:)=SUM(CheckCube(:,:,:),DIM=3,MASK=(CheckCube(:,:,:)==REAL(SelId).and.CheckCube(:,:,:)/=UNDEF))
           DO k=1,SIZE(Cube,DIM=3)
              WRITE(1,*) k, (k-WCS(3))*WCS(11)+WCS(6), &
                   SUM(Cube(:,:,k),MASK=mask2d(:,:)/=0.and.Cube(:,:,k)/=UNDEF)
           END DO
        ELSE !..do a sky spectrum, excluding all objects in the mask
           mask2d(:,:)=SUM(CheckCube(:,:,:),DIM=3,MASK=CheckCube(:,:,:)/=UNDEF)
          DO k=1,SIZE(Cube,DIM=3)
              WRITE(1,*) k, (k-WCS(3))*WCS(11)+WCS(6), &
                   SUM(Cube(:,:,k),MASK=mask2d(:,:)==0.and.Cube(:,:,k)/=UNDEF)
           END DO           
        END IF
        CLOSE(1)

     ELSE

        STOP "wrong masktype requested (options are: '2d' or '3d')!"

     END IF

     IF(plot) CALL showplot(TRIM(Output1d))
     IF(TRIM(outEPS)/="??") CALL saveEPS(TRIM(Output1d))

     IF(nint==nend) EXIT int_loop

     print *, "next Id to plot [use -1 to exit]"
     READ(*,*) SelId
     IF(SelId<0) EXIT int_loop
     !Output1d="temp.dat"

  END DO int_loop

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
  SelId=-1
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
     WRITE(*,'(a)') "  -id            <name>             : object id to plot (NO DEFAULT), use 0 for sky-spectrum"
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
  ELSEIF(TRIM(IdCube)=="??") THEN
     print *, "pleased provide the id cube with the -idcube option!"
     STOP
  ELSEIF(SelId==-1) THEN
     print *, "please provide the id number with the -id option!"
     STOP
  END IF

END SUBROUTINE ReadCommandLine


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
  WRITE(99,'(a)') "lim x y"
  WRITE(99,'(a)') "expand 1.3"
  WRITE(99,'(a)') "lweight 3"
  WRITE(99,'(a)') "box"
  WRITE(99,'(a)') "xlabel \AA"
  WRITE(99,'(a)') "ylabel flux"
  WRITE(99,'(a)') "connect x y"
  WRITE(99,'(a)') "dev x11"
  CLOSE(99)

  !..run the script
  CALL SYSTEM('sm < plot.sm')

  print *, "1d spectrum saved to: ", TRIM(outEPS)

END SUBROUTINE saveEPS


END PROGRAM main

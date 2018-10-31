PROGRAM SBProf

  USE CubeLib
  USE StatLib
  IMPLICIT NONE
  CHARACTER(len=500) :: InpFile, outfile, Output1d, VarCube, IdCube, rmbkg
  INTEGER(kind=4) :: SelId, smooth, cen(2), i, j, DimX, DImY, nbins, thisbin
  LOGICAL :: plot, logbin, logplot
  REAL(kind=4) :: binsize, maxdist, dist, logmin, meanclip, medianclip, sigma
  REAL(kind=4), ALLOCATABLE :: Prof(:)
  INTEGER(kind=4), ALLOCATABLE :: n(:)
  
  Verbosity=2

  !.. collect input parameters
  CALL ReadCommandLine

  CALL ReadCube(InpFile)
  DimX=SIZE(Cube,DIM=1); DimY=SIZE(Cube,DIM=2)
  ALLOCATE(Image(DimX,DImY))
  Image=Cube(:,:,1)

  IF(ALL(cen==-1)) cen=MAXLOC(Image)
  print *, "central pixel=", cen

  maxdist=MAX(DimX-cen(1),DimY-cen(2),cen(1),cen(2))
  IF(logbin) maxdist=log10(maxdist)
  nbins=INT(maxdist/binsize)
  print *, "maxdist=",maxdist
  print *, "nbins=",nbins
  ALLOCATE(Prof(nbins), n(nbins))

  print *, "MINMAX image=", MINVAL(Image), MAXVAL(image)
  CALL SigmaClip(PACK(Image,MASK=Image/=UNDEF),meanclip, medianclip, sigma)
  print *, "mean, median, sigma=", meanclip, medianclip, sigma
  IF(TRIM(rmbkg)=="mean") THEN
     print *, "removing bkg using clipped mean"
     WHERE(Image/=UNDEF) Image=Image-meanclip
  ELSEIF(TRIM(rmbkg)=="median") THEN
     print *, "removing bkg using clipped median"
     WHERE(Image/=UNDEF) Image=Image-medianclip
  END IF

  Prof=0.
  n=0
  DO j=1,DimY
     DO i=1,DimX
        IF(Image(i,j)==UNDEF) CYCLE
        IF(i==cen(1).and.j==cen(2)) CYCLE
        IF(logbin) THEN
           dist=log10(SQRT(REAL(i-cen(1))**2+REAL(j-cen(2))**2))
        ELSE
           dist=SQRT(REAL(i-cen(1))**2+REAL(j-cen(2))**2)
        END IF
        thisbin=INT((dist/maxdist)*nbins)+1
        IF(thisbin>nbins) CYCLE
        Prof(thisbin)=Prof(thisbin)+Image(i,j)
        n(thisbin)=n(thisbin)+1
     END DO
  END DO
  WHERE(n>0) Prof=Prof/n

  OPEN(1,file=Output1d,action="write")
  WRITE(1,*) "# bin          x (pixels)    Prof (image units)    nbins     std"
  DO i=1,nbins
     IF(logbin) THEN
        IF(logplot) THEN
           WRITE(1,*) i,REAL(i-0.5)*binsize, log10(MAX(Prof(i),logmin)), n(i), sigma/sqrt(MAX(1.,REAL(n(i))))
        ELSE
           WRITE(1,*) i,10**(REAL(i-0.5)*binsize),Prof(i), n(i), sigma/sqrt(MAX(1.,REAL(n(i))))
        END IF
     ELSE
        IF(logplot) THEN
           WRITE(1,*) i,log10((i-0.5)*binsize), log10(MAX(Prof(i),logmin)), n(i), sigma/sqrt(MAX(1.,REAL(n(i))))
        ELSE
           WRITE(1,*) i,(i-0.5)*binsize, Prof(i), n(i), sigma/sqrt(MAX(1.,REAL(n(i))))
        END IF
     END IF
  END DO
  CLOSE(1)

  IF(plot) CALL showplot(Output1d)



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
  plot=.true.
  smooth=0
  binsize=0.3
  cen=-1
  logbin=.true.
  logplot=.true.
  logmin=-22
  rmbkg='none'
  
  CALL GetArg(0,ExeName)
  narg=iargc()
  IF(narg>0) THEN
     CALL GetArg(1,fname)
  ELSE
     print *, " "
     WRITE(*,'(2a)')"        SBProf (part of CubEx package)   "
     WRITE(*,*) " "
     WRITE(*,'(a)') " produce circular averaged SB Profiles from images"
     WRITE(*,'(a)') " by Sebastiano Cantalupo (cantalupo@phys.ethz.ch)"
     print *, " "
     !WRITE(*,'(a)') "usage: Cube2Spc.x -cube <name> -idcube <name> -id <val> [-option <val>]"
     WRITE(*,'(a)') " "
     WRITE(*,'(a)') "  options:"
     WRITE(*,'(a)') "  -cube          <name>             : datacube name (NO DEFAULT)"
     WRITE(*,'(a)') "  -outname        <string>          : if a name is provided, save the 1d spectrum for the selected object in a ASCII file"
     WRITE(*,'(a)') "  -show           <bol>             : if .true. (default), shows the 1d spectrum using python"
     WRITE(*,'(a)') "  -cen            <int int>         : central pixel coordinates, if not provided uses brightest pixel"
     WRITE(*,'(a)') "  -binsize        <int>             : binsize in unit of image pixels, in dex if logbin=.true. (default=0.3)"
     WRITE(*,'(a)') "  -logbin         <bol>             : if .true. (default) use log bins in x direction"
     WRITE(*,'(a)') "  -logplot        <bol>             : if .true. (default) use log scale in output file and plot"
     WRITE(*,'(a)') "  -logmin         <bol>             : minimum value of y-axis for log plots (default=-22) "
     WRITE(*,'(a)') "  -rmbkg         <string>           : if set as 'mean' ('median') removes avg-sigma-clip (median), default='none'"
     STOP
  END IF

  !..read from command line
  DO i=1,narg,2
     CALL getarg(i,opt)
     CALL getarg(i+1,arg)
     SELECT CASE(TRIM(opt))
     CASE('-cube')       ; READ(arg,'(a)') InpFile
     !CASE('-masktype')   ; READ(arg,'(a)') masktype
     CASE('-outname')       ; READ(arg,'(a)') Output1d
     CASE('-idcube')        ; READ(arg,'(a)') IdCube
     CASE('-varcube')       ; READ(arg,'(a)') VarCube
     CASE('-id')            ; READ(arg,*) SelId
     !CASE('-interactive')   ; READ(arg,*) interactive
     CASE('-show')          ; READ(arg,*) plot
     CASE('-smoothr')       ; READ(arg,*) smooth
     !CASE('-outEPS')        ; READ(arg,'(a)') outEPS
     !CASE('-pixavg')        ; READ(arg,'(a)') pixavg
     !CASE('-bkgfile')       ; READ(arg,'(a)') bkgfile
     !CASE('-bkgthr')        ; READ(arg,*) masklayer_thr
     CASE('-binsize')        ; READ(arg,*) binsize
     CASE('-cen')            ; READ(arg,*) cen
     CASE('-logbin')         ; READ(arg,*) logbin
     CASE('-logmin')         ; READ(arg,*) logmin
     CASE('-logplot')        ; READ(arg,*) logplot
     CASE default
        print *, "command line argument ",TRIM(opt), " not recognized!"
        STOP
     END SELECT
  END DO

  !..convert logmin to linear scale
  logmin=10**(logmin)


!..perform few checks
  IF(TRIM(InpFile)=="??") THEN
     print *, "please provide the input datacube with the -cube option!"
     STOP
  !ELSEIF(TRIM(IdCube)=="??") THEN
  !   print *, "pleased provide the id cube with the -idcube option!"
  !   STOP
  !ELSEIF(SelId==-1) THEN
  !   print *, "please provide the id number with the -id option!"
  !   STOP
  END IF

  !IF(TRIM(masktype)=="3d".and.pixavg) STOP "-pixavg .true. is not allowed with -masktype 3d!"
  

END SUBROUTINE ReadCommandLine


SUBROUTINE showplot(fname)

  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN) :: fname

  
  !..write a python script:

  OPEN(99,file="plot.py",action="write")
  WRITE(99,'(a)') "import matplotlib.pyplot as plt"
  WRITE(99,'(a)') "import numpy as np"
  WRITE(99,'(a)') "from scipy import signal"
  !WRITE(99,'(3a)') "f2 = open('",TRIM(fname),"', 'r')"
  !WRITE(99,'(3a)') "lines = f2.readlines()"
  !WRITE(99,'(3a)') "x1 = []"
  !WRITE(99,'(3a)') "y1 = []"  
  !WRITE(99,'(a)') "for line in lines:"
  !WRITE(99,'(a)') "   p = line.split()"
  !WRITE(99,'(a)') "   x1.append(float(p[1]))"
  !WRITE(99,'(a)') "   y1.append(float(p[2]))"
  !WRITE(99,'(a)') " "  
  !WRITE(99,'(a)') "xv = np.array(x1)"
  !WRITE(99,'(a)') "yv = np.array(y1)"
  WRITE(99, '(3a)') "xv , yv = np.loadtxt('", TRIM(fname), "', usecols=(1,2), unpack=True)"
  IF(smooth>0) THEN
     WRITE(99,'(a,i2,a)') "w = signal.gaussian(20,std=",smooth,")"
     WRITE(99,'(a)') "yv = signal.convolve(yv, w, mode='same')"
  END IF
  WRITE(99,'(a)') "fig = plt.figure() "
  WRITE(99,'(a)') "splot = fig.add_subplot(111) "
  !WRITE(99,'(a)') "splot.set_xlabel('Observed wavelength ($\mathrm{\AA}$)') "
  !WRITE(99,'(a)') "splot.set_ylabel('Flux (10$^{-20}$ erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}$)') "
  WRITE(99,'(a)') "splot.plot(xv, yv) "
  WRITE(99,'(a)') "plt.show()"
  CLOSE(99)

  CALL SYSTEM("python plot.py")

END SUBROUTINE showplot

END PROGRAM SBProf

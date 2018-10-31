PROGRAM CubeStat

  USE StatLib
  USE CubeLib
  IMPLICIT NONE
  CHARACTER(len=500) :: inpcube, thisfile, arg
  INTEGER(kind=4)    :: n_ext, nfiles, ierr, niter, MaxIterations
  REAL(kind=4)       :: this_mean, this_median, this_sigma, clipval

  Verbosity=0

  IF(iargc()<1.or.iargc()>3) THEN
     print *, "Please provide cube/image or list name, e.g. CubeStat.x DATACUBE_FINAL.fits"
     print *, "and, optionally, the number of iterations for sigma clipping (default=5) and abs clipthreshold in sigma (default=3.0)"
     print *, "e.g., CubeStat.x DATACUBE_FINAL.fits 5 3.0"
     print *, "use 0 for no clipping, e.g., CubeStat.x DATACUBE_FINAL.fits 0"
     STOP
  END IF

  niter=5
  clipval=3.0

  CALL GetArg(1,inpcube)
  IF(iargc()==2) THEN
     CALL GetArg(2,arg)
     READ(arg,*) niter
  ELSEIF(iargc()==3) THEN
     CALL GetArg(3,arg)
     READ(arg,*) clipval
  ELSEIF(iargc()>3) THEN
     print *, "extra arguments on command line IGNORED"
     print *, "using: niter=",niter," clipval=", clipval
  END IF
     

  print *, "  filename,        mean,           sigma,            median    "


  !..check if input file is fits or a list
  IF(INDEX(TRIM(inpcube),".fits")>0) THEN

     CALL ReadCube(inpcube)
     IF(niter==0) THEN
        this_mean=Mean(PACK(Cube,MASK=Cube/=UNDEF))
        this_sigma=StdDev(PACK(Cube,MASK=Cube/=UNDEF))
        this_median=Median(PACK(Cube,MASK=Cube/=UNDEF))
     ELSE
        CALL SigmaClip(PACK(Cube,MASK=Cube/=UNDEF), this_mean, this_median, this_sigma, MaxIterations=niter, ClipVal=[-clipval,clipval])     
     END IF
     WRITE(*,*)  TRIM(inpcube), this_mean, this_sigma, this_median

  ELSE

     nfiles=0
     ierr=0
     OPEN(1,file=inpcube,action="read")
     DO
        READ(1,'(a)',iostat=ierr) thisfile
        IF(ierr/=0) EXIT
        nfiles=nfiles+1
        CALL ReadCube(thisfile)
        IF(niter==0) THEN
          this_mean=Mean(PACK(Cube,MASK=Cube/=UNDEF))
          this_sigma=StdDev(PACK(Cube,MASK=Cube/=UNDEF))
          this_median=Median(PACK(Cube,MASK=Cube/=UNDEF))
       ELSE
          CALL SigmaClip(PACK(Cube,MASK=Cube/=UNDEF), this_mean, this_median, this_sigma, MaxIterations=niter, ClipVal=[-clipval,clipval])
       END IF
       WRITE(*,*)  TRIM(thisfile), this_mean, this_sigma, this_median  
    END DO
    CLOSE(1)

  END IF


END PROGRAM CubeStat

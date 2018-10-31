PROGRAM main

  USE StatLib
  USE CubeLib
  USE OMP_LIB
  IMPLICIT NONE
  CHARACTER(len=350) :: InpFile, OutputCube, OutputVarCube, CubeType, ftype, OutputImage
  INTEGER(kind=4)    :: k, i, j, FilterSize(3), xmin(3), xmax(3), n_ext, pos(3), csize(3), ncpu, ngood, bpsize(3), orsize(3), maskpix(2)
  REAL(kind=4), ALLOCATABLE :: BKG(:,:,:), filter(:), Cube_(:,:,:), Var_(:,:,:)
  REAL(kind=4) :: this_med, meanval, sigma, original_median
  REAL(kind=4) :: t1, t2, InitCPUTime
  REAL, PARAMETER  :: st_conv(8)=[0.,0.,86400.,0.,3600.,60.,1.,0.001]
  INTEGER(kind=4)    ::  st_in(8), st1(8), st2(8)
  CHARACTER(len=300) :: date, time, zone
  LOGICAL, ALLOCATABLE :: zmask(:,:,:)

  Verbosity=2
  ncpu=0

  CALL report_time('init')

  !.. collect input parameters
  CALL ReadCommandLine

  !..read datacube
  CALL ReadCube(InpFile)

  ALLOCATE(zmask(SIZE(Cube,DIM=1),SIZE(Cube,Dim=2),SIZE(Cube,DIM=3)))
  zmask=.true.
  IF(ALL(maskpix/=0)) THEN
     zmask(:,:,maskpix(1):maskpix(2))=.false.
  END IF

  !print *, "Original Cube median=",Median(PACK(Cube,MASK=Cube/=UNDEF))

  !..compute original median
  IF(TRIM(OutputImage)/="??") THEN
     IF(ALLOCATED(Image)) DEALLOCATE(Image)
     ALLOCATE(Image(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2)))
     Image(:,:)=SUM(Cube(:,:,:), DIM=3, MASK=(Cube(:,:,:)/=UNDEF))
     original_median=Median(PACK(Image,MASK=Image/=UNDEF))
  END IF

  !..resample cube on a larger 3D mesh, if requested
  IF(ANY(bpsize/=1)) THEN

     FORALL(i=1:3) orsize(i)=SIZE(Cube,DIM=i)
     FORALL(i=1:3) csize(i)=CEILING(REAL(SIZE(Cube,DIM=i))/bpsize(i))
     !FORALL(i=1:3) csize(i)=SIZE(Cube,DIM=i)/bpsize(i)+1
     print *, "BKG cube size=",csize
     ALLOCATE(Cube_(csize(1),csize(2),csize(3)))

     IF(TRIM(OutputVarCube)/="??") THEN
        ALLOCATE(Var_(csize(1),csize(2),csize(3)))
        Var_=0
     END IF

     print *, "producing BKG cube..."

     CALL report_time('in')
     CALL Resample
     CALL report_time('end')

  ELSE
     FORALL(i=1:3) csize(i)=SIZE(Cube,DIM=i)
     ALLOCATE(Cube_(csize(1),csize(2),csize(3)))
     Cube_=Cube
  END IF

  FORALL(i=1:3) csize(i)=SIZE(Cube_,DIM=i)
  ALLOCATE(BKG(csize(1),csize(2),csize(3)))
  IF(ALLOCATED(Var_)) ALLOCATE(Var(orsize(1),orsize(2),orsize(3)))

  print *, "Median filtering with filter *size*=", 2*FilterSize(:)+1
  WRITE(*,'(a,$)') "filtering"

  CALL report_time('in')
  CALL MedianFiltering
  CALL report_time('end')


  !..perform subtraction, if requested
  IF(ANY(bpsize/=1)) THEN
        
     DO k=1,csize(3)
        DO j=1,csize(2)
           DO i=1,csize(1)

              pos(:)=[i,j,k]  !..resampled units

              !..original units
              xmin(:)=(pos(:)-1)*bpsize(:)+1
              xmax(:)=MIN(pos(:)*bpsize(:),orsize(:))

              IF(TRIM(cubetype)=="BKGsub") THEN
                 WHERE(Cube(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))/=UNDEF)
                    Cube(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))=&
                         Cube(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))-BKG(i,j,k)
                 END WHERE
              ELSE
                 WHERE(Cube(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))/=UNDEF)
                    Cube(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))=BKG(i,j,k)
                 END WHERE                
              END IF

              IF(ALLOCATED(Var)) Var(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))=Var_(i,j,k)

           END DO
        END DO
     END DO


  ELSE
     IF(TRIM(cubetype)=="BKGsub") THEN
        Cube=Cube-BKG
     ELSE
        Cube=BKG
     END IF
  END IF

  CALL WriteCube(OutputCube)

  IF(TRIM(OutputImage)/="??") THEN
     IF(ALLOCATED(Image)) DEALLOCATE(Image)
     ALLOCATE(Image(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2)))
     Image(:,:)=SUM(Cube(:,:,:), DIM=3, MASK=(Cube(:,:,:)/=UNDEF))
     CALL WriteImage(OutputImage)
     print *, "Image median=",Median(PACK(Image,MASK=Image/=UNDEF)), "original median=",original_median
  END IF

  IF(ALLOCATED(Var)) THEN
     Cube=Var
     CALL WriteCube(OutputVarCube)
     print *, "Mean Variance=", Mean(PACK(Var,MASK=VAR/=UNDEF))
  END IF

  print *, "end of Job! :-)"

CONTAINS

SUBROUTINE ReadCommandLine

  IMPLICIT NONE
  CHARACTER(len=300) :: ExeName, fname, string, ds, opt, arg, cmdfile
  CHARACTER(len=700) :: cmdline
  INTEGER :: narg, iarg, i, is
  LOGICAL :: savecmd


!..default
  InpFile="??"
  OutputCube="??"
  OutputVarCube="??"
  CubeType="BKGsub"
  Ftype="median"
  FilterSize=[0,0,2]
  bpsize=[1,1,10]
  OutputImage="??"
  savecmd=.true.
  maskpix=[0,0]

  CALL GetArg(0,ExeName)
  narg=iargc()
  IF(narg>0) THEN
     CALL GetArg(1,fname)
  ELSE
     print *, " "
     WRITE(*,'(2a)')"     CubeBKGSub (part of CubEx package) ",TRIM(version)
     WRITE(*,'(a)') "   BacKGround and continuum Subtraction software "
     WRITE(*,'(a)') " "
     WRITE(*,'(a)') " by Sebastiano Cantalupo (cantalupo@phys.ethz.ch)"
     print *, " "
     WRITE(*,'(a)') "usage: CubeBKGSub -cube <name> -out <name> "
     WRITE(*,'(a)') " "
     WRITE(*,'(a)') "example (continuum subtraction): CubeBKGSub -cube DATACUBE.fits -out DATACUBE.CSub.fits -bpsize '1 1 20' -bfrad '0 0 2' "
     WRITE(*,'(a)') "NB: for continuum subtraction the first two parameters in bpsize should be '1 1' and in bfrad '0 0' "
     WRITE(*,'(a)') "  "
     WRITE(*,'(a)') "  options:"
     WRITE(*,'(a)') "  -ftype              <string>           : filter type, 'sigclip' (median sigma clipping; default) or 'median' (faster)"
     WRITE(*,'(a)') "  -bpsize             <3 integer values> : size of one cell of the background/continuum cube in units of the original cube (default='1 1 10')."  
     WRITE(*,'(a)') "  -bfrad              <3 integer values> : background/continuum cube filter *radius* in x, y and z in units of the bkg cube spaxels size (default='0 0 2')" 
     WRITE(*,'(a)') "  -otype              <string>           : 'BKGsub' (default) outputs the bkg/continuum subtracted cube, 'BKG' outputs the estimated background/continuum" 
     WRITE(*,'(a)') "  -image              <string>           : if defined, produces a 2D image file of the output cube projecting along the z-direction" 
     WRITE(*,'(a)') "  -outvarcube         <string>           : if defined, produces a variance cube file estimated from sigmaclipping (NB: it sets ftype=sigclip!)"
     WRITE(*,'(a)') "                                            this is useful in case there is no variance associated with original cube (DON'T use that for MUSE cubes!)"
     WRITE(*,'(a)') "  -savecmd            <bol>              : if .true. (default) saves the command line to the file <outputcube>.cmd"
     WRITE(*,'(a)') "  -maskpix            <2 int values>     : initial and final pixel in z direction for line masking (default=0, no masking)"
     STOP
  END IF

  IF(savecmd) WRITE(11214,'(2a,$)') TRIM(Exename), " "

  !..read from command line
  DO i=1,narg,2
     CALL getarg(i,opt)
     CALL getarg(i+1,arg)
     SELECT CASE(TRIM(opt))
     CASE('-cube')          ; READ(arg,'(a)') InpFile
     CASE('-out')           ; READ(arg,'(a)') OutputCube
     CASE('-otype')         ; READ(arg,'(a)') CubeType
     CASE('-ftype')         ; READ(arg,'(a)') fType
     CASE('-bfrad')         ; READ(arg,*) FilterSize(:)
     CASE('-bpsize')        ; READ(arg,*) bpsize
     CASE('-maskpix')       ; READ(arg,*) maskpix
     CASE('-image')         ; READ(arg,*) OutputImage
     CASE('-outvarcube')    ; READ(arg,*) OutputVarCube
     CASE default
        print *, "command line argument ",TRIM(opt), " not recognized!"
        STOP
     END SELECT
     IF(savecmd) WRITE(11214,'(4a,$)') TRIM(opt),' "',TRIM(arg),'" '
  END DO

!..perform few checks
  IF(TRIM(InpFile)=="??") THEN
     print *, "please provide the input datacube with the -cube option!"
     STOP
  ELSEIF(TRIM(OutputCube)=="??") THEN
     STOP "please provide OutputCube filename with the -out option!" 
  END IF
  IF(TRIM(Cubetype)/="BKGsub".and.TRIM(Cubetype)/="BKG") STOP "Valid option for cubetype are 'BKGsub' or 'BKG'!"
 
  IF(TRIM(OutputVarCube)/="??") THEN
     IF(TRIM(ftype)=="median") THEN
        STOP " ftype=median not allowed because a variance file is requested!"
     END IF
  END IF

  IF(savecmd) THEN !..save the command line to a file
     is=INDEX(TRIM(OutputCube),".fits")
     IF(is==0) THEN
        cmdfile=TRIM(OutputCube)//".cmd"
     ELSE
        cmdfile=TRIM(OutputCube(1:is-1))//".cmd"
     END IF
     print *, "Cmd line stored in this file:", TRIM(cmdfile)
     REWIND(11214)
     READ(11214,'(a)') cmdline
     CALL SYSTEM('rm -f fort.11214')
     OPEN(1,file=cmdfile,action="write")
     WRITE(1,'(a)') TRIM(cmdline)
     CLOSE(1)
  END IF

END SUBROUTINE ReadCommandLine

!--------------------------------------------------------------

SUBROUTINE Resample

  IMPLICIT NONE

!$OMP PARALLEL shared(bpsize,orsize,Cube_,ftype,zmask) private(i,j,k,pos,xmin,xmax,ngood,filter,this_med,meanval,sigma)

  IF(ALLOCATED(filter)) DEALLOCATE(filter)
  ALLOCATE(filter(PRODUCT(bpsize(:))))
  filter=0

!$OMP DO SCHEDULE(DYNAMIC)
  DO k=1,csize(3)
     !WRITE(*,'(a,$)') "."
     DO j=1,csize(2)
        DO i=1,csize(1)
  
           pos(:)=[i,j,k]  !..resampled units

           !..convert to original units
           xmin(:)=(pos(:)-1)*bpsize(:)+1
           xmax(:)=MIN(pos(:)*bpsize(:),orsize(:))

           ngood=COUNT(Cube(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))/=UNDEF.and.zmask(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3)))

           IF(ngood>0) THEN

              filter(1:ngood)=PACK(Cube(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3)),&
                   MASK=Cube(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))/=UNDEF.and.zmask(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3)))
              
              IF(TRIM(ftype)=="median") THEN
                 this_med=Median(filter(1:ngood))
              ELSEIF(TRIM(ftype)=="sigclip") THEN
                 CALL SigmaClip(filter(1:ngood), meanval, this_med, sigma, Verbosity=1)
              END IF

!$OMP CRITICAL              
              Cube_(i,j,k)=this_med
              IF(ALLOCATED(Var_)) Var_(i,j,k)=sigma*sigma
!$OMP END CRITICAL

           ELSE

!$OMP CRITICAL
              Cube_(i,j,k)=UNDEF
!$OMP END CRITICAL

           END IF

        END DO
     END DO
  END DO
!$OMP END DO
  DEALLOCATE(filter)
!$OMP END PARALLEL

END SUBROUTINE RESAMPLE

!-------------------------------------------------------------

SUBROUTINE MedianFiltering

  IMPLICIT NONE

!$OMP PARALLEL shared(csize,FilterSize,BKG,ftype) private(i,j,k,pos,xmin,xmax,ngood,filter,this_med,meanval,sigma)

  IF(ALLOCATED(filter)) DEALLOCATE(filter)
  ALLOCATE(filter(PRODUCT(2*Filtersize(:)+1)))
  filter=0

!$OMP DO SCHEDULE(static)
  DO k=1,csize(3)
     !print *, "zpixel=",k,"over",csize(3)
     !WRITE(*,'(a,$)') "."
     DO j=1,csize(2)
        DO i=1,csize(1)
           
           IF(Cube_(i,j,k)/=UNDEF.or.TRIM(cubetype)=="BKG") THEN

              pos(:)=[i,j,k]
              xmin(:)=MAX(pos(:)-filtersize(:),1)
              xmax(:)=MIN(pos(:)+filtersize(:),csize(:))

              ngood=COUNT(Cube_(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))/=UNDEF)

              filter(1:ngood)=PACK(Cube_(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3)),&
                   MASK=Cube_(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))/=UNDEF)
              

              IF(TRIM(ftype)=="median") THEN
                 this_med=Median(filter(1:ngood))
              ELSE 
                 CALL SigmaClip(filter(1:ngood), meanval, this_med, sigma, Verbosity=1)
              END IF
             
!$OMP CRITICAL
              IF(TRIM(ftype)=="avgsigclip") THEN
                 BKG(i,j,k)=meanval
              ELSE
                 BKG(i,j,k)=this_med
              END IF
!$OMP END CRITICAL

           END IF

        END DO

     END DO
  END DO
!$OMP END DO

  DEALLOCATE(filter)

!$OMP END PARALLEL
           
  print *, "done"

END SUBROUTINE MedianFiltering

!--------------------------------------------------------------------------

SUBROUTINE report_time(what)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: what

    SELECT CASE(TRIM(what))
    CASE("init")
       CALL CPU_TIME(InitCPUTime)
       CALL DATE_AND_TIME(date,time,zone,st_in)
    CASE("in")
       CALL CPU_TIME(t1)
       CALL DATE_AND_TIME(date,time,zone,st1)
    CASE("end")
       CALL CPU_TIME(t2)          
       CALL DATE_AND_TIME(date,time,zone,st2)
       IF(Verbosity>=2) THEN
          print *, " "
          print *, "CPU_TIME        = ",t2-t1
          print *, "CPU_TIME so far = ",t2-InitCPUTime
          print *, "SYSTEM_TIME so far =", SUM((st2-st_in)*st_conv)
          print *, " "
       ENDIF
    END SELECT

END SUBROUTINE report_time


END PROGRAM main

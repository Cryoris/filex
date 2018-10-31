PROGRAM CubeCombine

  USE CubeLib
  USE StatLib
  USE OMP_LIB
  IMPLICIT NONE
  INTEGER, PARAMETER   :: MaxDimX=3000, MaxNCubes=1000
  CHARACTER(len=250) :: listname, combcube, fname(MaxNCubes), wname(MaxNCubes), NExpCube, comb, z1s, z2s, cutname, weightlist, masklist, maskname, &
       offsetlist
  INTEGER :: ncubes, ierr, i, n_ext, dx, dy, dz, x, y, z, MaxIterations, iargc, maxiter, MinClip, ncpus, CubeDims(3), status, &
       zbinsize, last_zbinsize, myid, zmin, zmax, nreads, ii, thisfile, lun, undef_var_counts, start_pixel(3), end_pixel(3), incs(3), group, &
       xmin, xmax, ymin, ymax, naxes, x1, x2, y1, y2, nHDUs
  INTEGER, ALLOCATABLE :: data_unit(:), var_unit(:), weight_unit(:), TotGoodVals(:,:), xoff(:), yoff(:), ThisCubeDims(:,:)
  REAL, ALLOCATABLE    :: meanval(:,:), this_med(:,:), sigma(:,:), this_var(:,:) 
  REAL :: NValClipThreshold, ClipVal_
  REAL(kind=4), ALLOCATABLE :: NExp(:,:,:), DataLayer(:,:,:), VarLayer(:,:,:), imMask(:,:,:), Weight(:,:,:), NExpLayer(:,:,:)
  LOGICAL :: reallocate, usePropVar, var_weight, anyf, writeNaN, updatec
  LOGICAL :: GoodMask(MaxNCubes)
  INTEGER(kind=4)    ::  st_in(8), st1(8), st2(8)
  CHARACTER(len=300) :: date, time, zone
  REAL(kind=4) :: t1, t2, InitCPUTime
  REAL, PARAMETER  :: st_conv(8)=[0.,0.,86400.,0.,3600.,60.,1.,0.001]

!--- CPU_Time and System Time initialization
  CALL report_time("init")
 
  CALL ReadCommandLine

  VERBOSITY=1

  !..count cubes
  ncubes=0
  OPEN(1,file=listname)
  ierr=0
  print *, "reading: ", TRIM(listname)
  DO
     READ(1,'(a)',iostat=ierr) fname(ncubes+1)
     IF(ierr/=0) EXIT
     ncubes=ncubes+1
  END DO
  IF(TRIM(weightlist)/="??") THEN
     print *, "apply weights from files in this list: ", TRIM(weightlist)
     OPEN(11,file=weightlist,action="read")
     DO i=1,ncubes
        READ(11,'(a)') wname(i)
     END DO
     CLOSE(11)
  END IF

  print *, "ncubes=", ncubes

  ALLOCATE(xoff(ncubes),yoff(ncubes))
  xoff=0 ; yoff=0
  IF(TRIM(offsetlist)/="??") THEN
     OPEN(11,file=offsetlist,action="read")
     DO i=1,ncubes
        READ(11,*) xoff(i),yoff(i)
     END DO
     CLOSE(11)
     print *, "applying offsets from file: ", TRIM(offsetlist)
  END IF


 !..get number of threads
!$OMP PARALLEL
  myid=omp_get_thread_num()+1
  IF(myid==1) THEN
     ncpus=omp_get_num_threads()
     print *, "ncpus=",ncpus
  END IF
!$OMP END PARALLEL

 
  !..get cube size 
  ALLOCATE(ThisCubeDims(ncubes,3))
  REWIND(1); READ(1,'(a)') fname(1)
  IF(TRIM(offsetlist)=="??") THEN !..all cubes should have the same dimension
     !..get the size from the first one
     CALL GetCubeSize(fname(1), CubeDims)
     print *, "CubeDims=", CubeDims
     CLOSE(1)
     xmin=1 ; ymin=1 ; xmax=CubeDims(1) ; ymax=CubeDims(2)
     DO i=1,ncubes ; ThisCubeDims(i,:)=CubeDims(:) ; END DO
  ELSE   !..get cube dimension with offsets
     xmin=9999 ; xmax=1
     ymin=9999 ; ymax=1
     DO i=1,ncubes
        CALL GetCubeSize(fname(i), CubeDims)
        ThisCubeDims(i,:)=CubeDims(:)
        xmin=MIN(1+xoff(i),xmin) ; ymin=MIN(1+yoff(i),ymin)
        xmax=MAX(CubeDims(1)+xoff(i),xmax) ; ymax=MAX(CubeDims(2)+yoff(i),ymax)
     END DO
     print *, "xmin=",xmin, "xmax=",xmax, "ymin=",ymin, "ymax=",ymax
     CubeDims(1)=xmax-xmin+1 ; CubeDims(2)=ymax-ymin+1
     print *, "Final CubeDims with offsets=", CubeDims
  END IF

  IF(CubeDims(3)==1) THEN
     naxes=2
  ELSE
     naxes=3
  END IF
  print *, "naxes=",naxes

  !..allocate and read mask files if needed
  IF(TRIM(masklist)/="??") THEN
     ALLOCATE(imMask(ncubes,CubeDims(1),CubeDims(2)))
     print *, "reading image masks from file:", TRIM(masklist)
     OPEN(1,file=masklist,action="read")
     DO i=1,ncubes
        READ(1,'(a)') maskname
        CALL ReadCube(maskname)
        imMask(i,:,:)=Cube(:,:,1)
     END DO
     DEALLOCATE(Cube)
     CLOSE(1)
  END IF

  CALL report_time("in")

  !..allocate final cubes on the stack                                                                                                                      
  ALLOCATE(Cube(CubeDims(1),CubeDims(2),CubeDims(3)),&
       Var(CubeDims(1),CubeDims(2),CubeDims(3)))
  IF(TRIM(NExpCube)/="??") THEN
     ALLOCATE(NExp(CubeDims(1),CubeDims(2),CubeDims(3)))
  ELSE !..allocate dummy array to avoid problems in the parallel section
     ALLOCATE(NExp(1,1,1))
  END IF

  !..start parallel section

  !$OMP PARALLEL DEFAULT(NONE) SHARED(CubeDims,ncubes,fname,wname, minclip, maxiter, comb, ClipVal_,&
  !$OMP& NValClipThreshold, Verbosity, weightlist, var_weight, usePropVar, NExpCube, NExp, masklist, imMask, &
  !$OMP& Cube, Var, ThisCubeDims, xoff, yoff, naxes, xmin, ymin) &
  !$OMP& PRIVATE(DataLayer, VarLayer, data_unit, var_unit, weight_unit, Weight, NExpLayer, status, &
  !$OMP& i, incs, group, x, y, z, start_pixel, end_pixel, anyf, TotGoodVals, this_med, meanval, &
  !$OMP& sigma, reallocate, this_var, GoodMask, myid, x1, x2, y1, y2)

  myid=omp_get_thread_num()+1

  !..allocate layers and unit array for each process
  ALLOCATE(DataLayer(ncubes,CubeDims(1),CubeDims(2)),VarLayer(ncubes,CubeDims(1),CubeDims(2)), data_unit(ncubes), var_unit(ncubes), &
       weight_unit(ncubes))
  ALLOCATE(Weight(ncubes,CubeDims(1),CubeDims(2)),NExpLayer(ncubes,CubeDims(1),CubeDims(2)))
  ALLOCATE(TotGoodVals(CubeDims(1),CubeDims(2)),this_med(CubeDims(1),CubeDims(2)),meanval(CubeDims(1),CubeDims(2)),&
       sigma(CubeDims(1),CubeDims(2)),this_var(CubeDims(1),CubeDims(2)))

  !..initialize DataLayer
  DataLayer=UNDEF

  !..open data and var extensions for all files
  status=0

  !..get layer data and perform combination
  incs(:)=1
  group=1
  IF(myid==1) print *, "working on layers (%): "
 

  !$OMP DO SCHEDULE(DYNAMIC)
  layerloop: DO z=1,CubeDims(3)

     !print *, myid, z

     IF(MOD(z,50)==0) THEN

        WRITE(*,'(f5.1,$)') REAL(z)/CubeDims(3)*100. 

     END IF

     start_pixel(:)=[1,1,z]

     DO i=1,ncubes

        end_pixel(:)=[ThisCubeDims(i,1),ThisCubeDims(i,2),z]

        x1=1+xoff(i)-(xmin-1) 
        x2=ThisCubeDims(i,1)+x1-1
        y1=1+yoff(i)-(ymin-1)
        y2=ThisCubeDims(i,2)+y1-1

        !print *, x1,x2,y1,y2

        data_unit(i)=myid*1000+500+i
        CALL ftdopn(data_unit(i),TRIM(fname(i))//"[1]",0,status)
        IF(status/=0) THEN
           !..try to open a single extension file
           status=0
           CALL ftdopn(data_unit(i),TRIM(fname(i)),0,status)
           IF(status/=0) THEN
              print *, "problem reading file:",TRIM(fname(i)),"!"
              print *, "data_unit=",data_unit(i)
              STOP
           END IF
        END IF
        CALL ftgsve(data_unit(i),group,naxes,ThisCubeDims(i,1:naxes),start_pixel(1:naxes),end_pixel(1:naxes),incs(1:naxes),UNDEF,&
             DataLayer(i,x1:x2,y1:y2),anyf,status)
        IF(status/=0) STOP " Problem reading data subset!"
        CALL ftclos(data_unit(i),status)

        IF(UsePropVar.or.var_weight) THEN

           var_unit(i)=myid*1000+700+i
           CALL ftdopn(var_unit(i),TRIM(fname(i))//"[2]",0,status)
           IF(status/=0) THEN
              print *, "problem reading file:",TRIM(fname(i)),"[2] !"
              print *, "var_unit=",var_unit(i)
              STOP
           END IF
           CALL ftgsve(var_unit(i),group,naxes,ThisCubeDims(i,1:naxes),start_pixel(1:naxes),end_pixel(1:naxes),incs(1:naxes),UNDEF,&
                VarLayer(i,x1:x2,y1:y2),anyf,status)
           IF(status/=0) STOP " Problem reading var subset!"
           CALL ftclos(var_unit(i),status)

        ELSE

           VarLayer(i,:,:)=1.

        END IF

        IF(TRIM(masklist)/="??") THEN
           WHERE(imMask(i,:,:)==0.) &
             DataLayer(i,x1:x2,y1:y2)=UNDEF
        END IF

        IF(TRIM(weightlist)/="??") THEN
           weight_unit(i)=myid*1000+900+i
           CALL ftdopn(weight_unit(i),TRIM(wname(i))//"[1]",0,status)
           IF(status/=0) THEN
              print *, "problem reading file:",TRIM(wname(i)),"[2] !"
              print *, "weight_unit=",weight_unit(i)
              STOP
           END IF
        !END IF
        !IF(TRIM(weightlist)/="??") THEN
           CALL ftgsve(weight_unit(i),group,naxes,ThisCubeDims(i,1:naxes),start_pixel(1:naxes),end_pixel(1:naxes),incs(1:naxes),UNDEF,&
             Weight(i,x1:x2,y1:y2),anyf,status)
           IF(status/=0) STOP " Problem reading weight subset!"           
           CALL ftclos(weight_unit(i),status)
        ELSEIF(var_weight) THEN
           WHERE(VarLayer(i,:,:)>0) 
              Weight(i,:,:)=1./VarLayer(i,:,:)
           ELSEWHERE
              Weight(i,:,:)=0.
           END WHERE
        END IF

     END DO

     !print *, "combining..."

     yloop:DO y=1,CubeDims(2)
        xloop:DO x=1,CubeDims(1)  !..buffer x dimension

           IF(ALL(DataLayer(:,x,y)==UNDEF)) THEN
              
              meanval(x,y)=UNDEF
              this_med(x,y)=UNDEF
              sigma(x,y)=UNDEF
              this_var(x,y)=UNDEF
              TotGoodVals(x,y)=0

           ELSE

              IF(TRIM(weightlist)/="??".or.var_weight) THEN
                 CALL SigmaClip(array=PACK(DataLayer(:,x,y),MASK=DataLayer(:,x,y)/=UNDEF), meanclip=meanval(x,y), medianclip=this_med(x,y), finalsigma=sigma(x,y), &
                      ClipVal=[-ClipVal_,ClipVal_], NValClipThreshold=MinClip, MaxIterations=maxiter, Verbosity=0, TotGoodVals=TotGoodVals(x,y), &
                      weight=PACK(Weight(:,x,y),MASK=DataLayer(:,x,y)/=UNDEF), GoodMask=GoodMask(1:COUNT(DataLayer(:,x,y)/=UNDEF)), &
                      VarArray=PACK(VarLayer(:,x,y),MASK=DataLayer(:,x,y)/=UNDEF), PropVar=this_var(x,y))
              ELSE
                 CALL SigmaClip(array=PACK(DataLayer(:,x,y),MASK=DataLayer(:,x,y)/=UNDEF), meanclip=meanval(x,y), medianclip=this_med(x,y), finalsigma=sigma(x,y), &
                      ClipVal=[-ClipVal_,ClipVal_], NValClipThreshold=MinClip, MaxIterations=maxiter, Verbosity=0, TotGoodVals=TotGoodVals(x,y), &
                      GoodMask=GoodMask(1:COUNT(DataLayer(:,x,y)/=UNDEF)), VarArray=PACK(VarLayer(:,x,y),MASK=DataLayer(:,x,y)/=UNDEF), PropVar=this_var(x,y))
              END IF

           END IF

        END DO xloop

     END DO yloop

     !print *, "assigning values..."

     !$OMP CRITICAL
     IF(TRIM(comb)=="mean") THEN
        Cube(:,:,z)=meanval(1:CubeDims(1),1:CubeDims(2))
     ELSE
        Cube(:,:,z)=this_med(1:CubeDims(1),1:CubeDims(2))
     END IF
     
     !print *, "assigning variance..."

     IF(usePropVar) THEN
        Var(:,:,z)=this_var(1:CubeDims(1),1:CubeDims(2))
     ELSE  !..approximate variance of the sample mean with the square of the standard error of the mean  
        WHERE(TotGoodVals(1:CubeDims(1),1:CubeDims(2))>1.and.sigma(1:CubeDims(1),1:CubeDims(2))/=UNDEF)
           Var(:,:,z)=sigma(1:CubeDims(1),1:CubeDims(2))**2/TotGoodVals(1:CubeDims(1),1:CubeDims(2))
        ELSEWHERE
           Var(:,:,z)=UNDEF
        END WHERE
     END IF
     

     IF(TRIM(NExpCube)/="??") THEN
        NExp(:,:,z)=TotGoodVals(1:CubeDims(1),1:CubeDims(2))
     END IF
     !$OMP END CRITICAL

     !print *, "done"

  END DO layerloop
  !$OMP END DO

  !..close files
  !DO i=1,ncubes
  !   CALL ftclos(data_unit(i),status)
  !   CALL ftclos(var_unit(i),status)
  !END DO

  !print *, "done"

  !$OMP END PARALLEL

  !print *, "done"

  print *, " "

  Verbosity=2
  CALL report_time("end")

  print *, "writing combined cube..."
  print *, "minmax=",MINVAL(Cube,MASK=Cube/=UNDEF),MAXVAL(Cube,MASK=Cube/=UNDEF)
  IF(ALLOCATED(Var)) THEN
     print *, "minmax variance=",MINVAL(Var,MASK=Cube/=UNDEF.and.Var/=UNDEF), MAXVAL(Var,MASK=Cube/=UNDEF.and.Var/=UNDEF)
  END IF

  IF(updatec) THEN
     CALL UpdateCube(InpFile=fname(1),OutFile=combcube,writeNaN=writeNaN,author="CubEx_CubeCombine_"//TRIM(version),nHDUs=nHDUs)
     IF(ALLOCATED(Var).and.nHDUs==1) THEN
        Cube=Var
        combcube=combcube(1:INDEX(combcube,".fits")-1)//".VAR.fits"
        CALL UpdateCube(InpFile=fname(1),OutFile=combcube,writeNaN=writeNaN,author="CubEx_CubeCombine_"//TRIM(version))
     END IF
  ELSE
     CALL WriteCube(combcube,writeNaN=writeNaN,author="CubEx_CubeCombine_"//TRIM(version),n_ext=2)
  END IF

  IF(TRIM(NExpCube)/="??") THEN
     Cube=NExp
     CALL WriteCube(NExpCube, n_ext=1)
  END IF


CONTAINS

 SUBROUTINE ReadCommandLine

  IMPLICIT NONE
  CHARACTER(len=300) :: ExeName, fname, string, ds, opt, arg 
  INTEGER :: narg, iargc, i, is
  LOGICAL :: ex


!..default
  listname="??"
  combcube="??"
  NExpCube="??"
  weightlist="??"
  masklist="??"
  ClipVal_=3.
  MinClip=4
  maxiter=5
  comb="mean"
  usePropVar=.true.
  var_weight=.false.
  writeNaN=.true.
  updatec=.true.
  offsetlist="??"

  CALL GetArg(0,ExeName)
  narg=iargc()
  IF(narg>0) THEN
     CALL GetArg(1,fname)
  ELSE
     print *, " "
     WRITE(*,'(2a)')"        CubeCombine (part of CubEx package) ",TRIM(version)
     WRITE(*,'(a)') " "
     WRITE(*,'(a)') " by Sebastiano Cantalupo (cantalupo@phys.ethz.ch)"
     print *, " "
     WRITE(*,'(a)') "usage: CubeComb.x -list <name> -out <name> [-option <val>]"
     WRITE(*,'(a)') " "
     WRITE(*,'(a)') "  options:"
     WRITE(*,'(a)') "  -list               <string>          : list of files to combine (NO default)"
     WRITE(*,'(a)') "  -out                <string>          : name of output cube      (NO default)"
     WRITE(*,'(a)') "  -offsetlist         <string>          : name of the file with spatial offsets in pixels, if not provided no offsets are applied"
     WRITE(*,'(a)') "  -weightlist         <string>          : if defined, applies weight using values from files in this list (inverse squared)"
     WRITE(*,'(a)') "  -varweight          <bol>             : if .true.,  applies weight using variance cube (default=.false.)"
     WRITE(*,'(a)') "  -masklist           <string>          : if defined, applies a 2D mask during combination using values from files in this list (bad pixels=0)"
     WRITE(*,'(a)') "  -outexp             <string>          : if a filename is provided, it produces a cube with nexp"
     WRITE(*,'(a)') "  -comb               <string>          : options: mean/median (default=mean)"
     WRITE(*,'(a)') "  -clipval            <real>            : abs of clipping value in sigma (default=3.) "
     WRITE(*,'(a)') "  -minclip            <int>             : minimum value of pixels to stop iterative clipping if convergence is not reachead (default=4)"
     WRITE(*,'(a)') "  -maxiter            <int>             : maximum number of iterations (default=5)"
     WRITE(*,'(a)') "  -writeNaN           <bol>             : if .true. writes NaN instead of UNDEF (default=.true.)"
     WRITE(*,'(a)') "  -propvar            <bol>             : if .true., the final variance will be the propagated variance, otherwise, the statistical variance from "
     WRITE(*,'(a)') "                                           cube combination will be provided (default=.true.) and original variance will not be read"
     STOP
  END IF

  !..read from command line
  DO i=1,narg,2
     CALL getarg(i,opt)
     CALL getarg(i+1,arg)
     SELECT CASE(TRIM(opt))
     CASE('-list')          ; READ(arg,'(a)') listname
     CASE('-out')           ; READ(arg,'(a)') combcube
     CASE('-weightlist')    ; READ(arg,'(a)') weightlist
     CASE('-offsetlist')    ; READ(arg,'(a)') offsetlist
     CASE('-outexp')        ; READ(arg,'(a)') NExpCube
     CASE('-clipval')       ; READ(arg,*) ClipVal_
     CASE('-minclip')       ; READ(arg,*) MinClip
     CASE('-maxiter')       ; READ(arg,*) maxiter
     CASE('-comb')          ; READ(arg,'(a)') comb
     CASE('-propvar')       ; READ(arg,*) usePropVar
     CASE('-masklist')      ; READ(arg,'(a)') masklist
     CASE('-varweight')     ; READ(arg,*) var_weight
     CASE('-writeNaN')      ; READ(arg,*) writeNaN
     CASE('-update')        ; READ(arg,*) updatec
     CASE default
        print *, "command line argument ",TRIM(opt), " not recognized!"
        STOP
     END SELECT
  END DO

!..perform few checks
  IF(TRIM(comb)/="mean".and.TRIM(comb)/="median") STOP "-comb must be 'mean' OR 'median'"


END SUBROUTINE ReadCommandLine

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


!--------------------------------------------

  


END PROGRAM CubeCombine

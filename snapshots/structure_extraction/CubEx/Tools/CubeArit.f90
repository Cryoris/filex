PROGRAM CubeArit

  USE StatLib
  USE CubeLib
  IMPLICIT NONE
  CHARACTER(len=500) :: inpcube1, inpcube2, op, resultcube, thisfile1, thisfile2, resultcube_, string
  INTEGER(kind=4)    :: n_ext, nfiles, ierr, niter, MaxIterations, i
  REAL(kind=4)       :: this_mean, this_median, this_sigma, const, rdum
  REAL(kind=4), ALLOCATABLE :: Cube1(:,:,:), Cube2(:,:,:), Var1(:,:,:), Var2(:,:,:)
  LOGICAL :: switch, ex
  
  IF(iargc()<1.or.iargc()>4) THEN
     print *, " usage:"
     print *, " "
     print *, " CubeArit.x Cube1.fits <operator> Cube2.fits ResultCube.fits "
     print *, " operators= '+', '-', '*', '/' 'SNR' "
     print *, " e.g., CubeArit Cube1.fits - Cube2.fits Result.fits "
     print *, " "
     print *, " input files may be lists (or just the first one), in this case the last argument"
     print *, " is the output base (output fits will have a running index),"
     print *, " e.g., CubeArit Cube1.list - Cube2.list ResultCube "
     print *, " "
     print *, " third argument may also be a REAL constant, "
     print *, " e.g., CubeArit Cube1.fits * 1.e-20 ResultCube.fits "
     print *, " "
     print *, " or a ASCII file with a multiplicative factor for each layer (two columns format), "
     print *, " e.g., CubeArit Cube1.fits * LayerCorr.txt ResultCube.fits "
     print *, " in this case: append + to the file names if you want to apply the operators to"
     print *, " multiextension files and in order to produce a multiext cube"
     STOP
  END IF

  switch=.false.

  CALL GetArg(1,inpcube1)
  CALL GetArg(2,op)
  CALL GetArg(3,inpcube2)
  CALL GetArg(4,resultcube)

  !..check if input file is fits or a list
  IF(INDEX(TRIM(inpcube1),".fits")>0) THEN !..fits file

     IF(INDEX(TRIM(inpcube1),"+")==0) THEN !..read only one extension
        CALL ReadCube(inpcube1)
        ALLOCATE(Cube1(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2),SIZE(Cube,DIM=3)))
        Cube1=Cube
     ELSE
        inpcube1=inpcube1(1:INDEX(TRIM(inpcube1),"+")-1) !..cut the + at the end
        CALL ReadCube(inpcube1, n_ext=2)
        ALLOCATE(Cube1(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2),SIZE(Cube,DIM=3)))
        Cube1=Cube       
        ALLOCATE(Var1(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2),SIZE(Cube,DIM=3)))
        Var1=Var
     END IF

     !..check if 3rd argument is a fits, a constant or a ascii file
     IF(INDEX(TRIM(inpcube2),".fits")==0) THEN !..not a fits file 
        !..check if the file is ascii
        INQUIRE(file=TRIM(inpcube2),exist=ex)
        IF(.not.ex) THEN
           !..try to read a number
           READ(inpcube2,*,iostat=ierr) const
           IF(ierr/=0) STOP "third argument must be a file or a REAL constant!"
           Cube=const
        ELSE
           OPEN(99,file=inpcube2,action='read')
           DO i=1,SIZE(Cube,DIM=3)
              READ(99,*) rdum,const
              Cube(:,:,i)=const
              Var(:,:,i)=const*const
           END DO
           CLOSE(99)
        END IF
     ELSE
        CALL ReadCube(inpcube2)
     END IF

     CALL ApplyOp

     IF(INDEX(TRIM(resultcube),'+')==0) THEN
        !CALL WriteCube(resultcube)
        CALL UpdateCube(InpFile=inpcube1,OutFile=resultcube,writeNaN=.true.)
     ELSE
        resultcube=resultcube(1:INDEX(TRIM(resultcube),"+")-1) !..cut the + at the end
        CALL UpdateCube(InpFile=inpcube1,OutFile=resultcube,writeNaN=.true.)
        !CALL WriteCube(resultcube,multiext=.true.,n_ext=2)
     END IF

  ELSEIF(INDEX(TRIM(inpcube2),".fits")>0) THEN !..second file is fits and first file is a list

     !..switch file order
     switch=.true.
     CALL ReadCube(inpcube2)
     ALLOCATE(Cube1(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2),SIZE(Cube,DIM=3)))
     Cube1=Cube    
     
     nfiles=0
     ierr=0
     OPEN(1,file=inpcube1,action="read")
     DO

        READ(1,'(a)',iostat=ierr) thisfile1
        IF(ierr/=0) EXIT

        nfiles=nfiles+1

        CALL ReadCube(thisfile1)

        CALL ApplyOp

        WRITE(string,'(i4.4)') nfiles
        resultcube_=TRIM(resultcube)//"."//TRIM(string)//".fits"

        WRITE(*,'(7a)') TRIM(thisfile1), " ", TRIM(op), " ", TRIM(inpcube2), " --> ", TRIM(resultcube_)

        !CALL WriteCube(resultcube_)
        CALL UpdateCube(InpFile=thisfile1,OutFile=resultcube_,writeNaN=.true.)
     END DO

  ELSE  !..both files are lists, or first file is a list and 3rd argument is a constant

     !..check if 3rd argument is a constant
     READ(inpcube2,*,iostat=ierr) const

     IF(ierr/=0) THEN !..both files are lists

        nfiles=0
        ierr=0
        OPEN(1,file=inpcube1,action="read")
        OPEN(2,file=inpcube2,action="read")
        DO

           READ(1,'(a)',iostat=ierr) thisfile1
           IF(ierr/=0) EXIT
           READ(2,'(a)',iostat=ierr) thisfile2
           IF(ierr/=0) EXIT

           nfiles=nfiles+1
        
           CALL ReadCube(thisfile1)

           IF(ALLOCATED(Cube1)) DEALLOCATE(Cube1)
           ALLOCATE(Cube1(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2),SIZE(Cube,DIM=3)))
           Cube1=Cube

           CALL ReadCube(thisfile2)

           CALL ApplyOp

           WRITE(string,'(i4.4)') nfiles
           resultcube_=TRIM(resultcube)//"."//TRIM(string)//".fits"
           WRITE(*,'(7a)') TRIM(thisfile1), " ", TRIM(op), " ", TRIM(thisfile2), " --> ", TRIM(resultcube_)

           !CALL WriteCube(resultcube_)
           CALL UpdateCube(InpFile=thisfile1,OutFile=resultcube_,writeNaN=.true.)

        END DO

     ELSE !..first arg is a list, 3rd argument is a constant

        nfiles=0
        ierr=0
        OPEN(1,file=inpcube1,action="read")
        DO

           READ(1,'(a)',iostat=ierr) thisfile1
           IF(ierr/=0) EXIT

           nfiles=nfiles+1
        
           CALL ReadCube(thisfile1)

           IF(ALLOCATED(Cube1)) DEALLOCATE(Cube1)
           ALLOCATE(Cube1(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2),SIZE(Cube,DIM=3)))
           Cube1=Cube

           Cube=const

           CALL ApplyOp

           WRITE(string,'(i4.4)') nfiles
           resultcube_=TRIM(resultcube)//"."//TRIM(string)//".fits"
           WRITE(*,'(7a)') TRIM(thisfile1), " ", TRIM(op), " ", const, " --> ", TRIM(resultcube_)

           !CALL WriteCube(resultcube_)
           CALL UpdateCube(InpFile=thisfile1,OutFile=resultcube_,writeNaN=.true.)

        END DO    

     END IF

  END IF


CONTAINS

  SUBROUTINE ApplyOp

    IMPLICIT NONE


     IF(SIZE(Cube,DIM=1)/=SIZE(Cube1,DIM=1).or.SIZE(Cube,DIM=2)/=SIZE(Cube1,DIM=2).or.&
          SIZE(Cube,DIM=3)/=SIZE(Cube1,DIM=3)) STOP "Input Cube have different dimensions!"


     SELECT CASE(TRIM(op))
     CASE("+")
        WHERE(Cube/=UNDEF.and.Cube1/=UNDEF) 
           Cube=Cube1+Cube
        ELSEWHERE
           Cube=UNDEF
        END WHERE
     CASE("-")
        IF(.not.switch) THEN
           WHERE(Cube/=UNDEF.and.Cube1/=UNDEF) 
              Cube=Cube1-Cube
           ELSEWHERE
              Cube=UNDEF
           END WHERE
        ELSE
           WHERE(Cube/=UNDEF.and.Cube1/=UNDEF) 
              Cube=Cube-Cube1
           ELSEWHERE
              Cube=UNDEF
           END WHERE
        END IF
     CASE("*")
        WHERE(Cube/=UNDEF.and.Cube1/=UNDEF) 
           Cube=Cube1*Cube
        ELSEWHERE
           Cube=UNDEF
        END WHERE
        IF(ALLOCATED(Var).and.ALLOCATED(Var1)) THEN
           WHERE(Cube/=UNDEF.and.Cube1/=UNDEF) 
              Var=Var1*Var
           ELSEWHERE
              Var=UNDEF
           END WHERE
        END IF
     CASE("/")
        IF(.not.switch) THEN
           WHERE(Cube/=UNDEF.and.Cube1/=UNDEF.and.Cube/=0.) 
              Cube=Cube1/Cube
           ELSEWHERE
              Cube=UNDEF
           END WHERE
           IF(ALLOCATED(Var).and.ALLOCATED(Var1)) THEN
              WHERE(Cube/=UNDEF.and.Cube1/=UNDEF.and.Var/=UNDEF.and.Var1/=UNDEF.and.Var/=0) 
                 Var=Var1/Var
              ELSEWHERE
                 Var=UNDEF
              END WHERE
           END IF
        ELSE
           WHERE(Cube/=UNDEF.and.Cube1/=UNDEF.and.Cube/=0.) 
              Cube=Cube/Cube1
           ELSEWHERE
              Cube=UNDEF
           END WHERE      
           IF(ALLOCATED(Var).and.ALLOCATED(Var1)) THEN
              WHERE(Cube/=UNDEF.and.Cube1/=UNDEF.and.Var/=UNDEF.and.Var1/=UNDEF.and.Var1/=0) 
                 Var=Var/Var1
              ELSEWHERE
                 Var=UNDEF
              END WHERE
           END IF
        END IF
     CASE("SNR")
        WHERE(Cube/=UNDEF.and.Cube1/=UNDEF.and.Cube>0.) 
           Cube=Cube1/sqrt(Cube)
        ELSEWHERE
           Cube=UNDEF
        END WHERE
     CASE default
        STOP "selected operator not available!"
     END SELECT

END SUBROUTINE ApplyOp

END PROGRAM CubeArit

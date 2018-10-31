PROGRAM CubeReplace

  USE CubeLib
  IMPLICIT NONE
  CHARACTER(len=500) :: inpcube, outcube, val_old_, val_new_, NaN
  REAL(kind=4) :: val_old, val_new, neg_number
  INTEGER(kind=4) :: ierr, xmin, xmax, ymin, ymax, zmin, zmax
  LOGICAL :: put_NaN, replace_NaN, replace_all

  
  ierr=0
  neg_number=-1.

  CALL ReadCommandLine

  CALL ReadCube(inpcube)

  !..adjust region if necessary
  IF(xmax==-1) xmax=SIZE(Cube,DIM=1)
  IF(ymax==-1) ymax=SIZE(Cube,DIM=2)
  IF(zmax==-1) zmax=SIZE(Cube,DIM=3)
  
  IF(replace_NaN) THEN
     IF(replace_all) Cube(xmin:xmax,ymin:ymax,zmin:zmax)=UNDEF
     WHERE(Cube(xmin:xmax,ymin:ymax,zmin:zmax)==UNDEF) Cube(xmin:xmax,ymin:ymax,zmin:zmax)=val_new
  ELSE
     IF(replace_all) Cube(xmin:xmax,ymin:ymax,zmin:zmax)=val_old
     IF(put_NaN) THEN
        WHERE(Cube(xmin:xmax,ymin:ymax,zmin:zmax)==val_old) Cube(xmin:xmax,ymin:ymax,zmin:zmax)=sqrt(neg_number)
     ELSE
        WHERE(Cube(xmin:xmax,ymin:ymax,zmin:zmax)==val_old) Cube(xmin:xmax,ymin:ymax,zmin:zmax)=val_new
     END IF
  END IF
     
  CALL WriteCube(outcube)


CONTAINS

  SUBROUTINE ReadCommandLine

      IMPLICIT NONE
  CHARACTER(len=300) :: ExeName, fname, string, ds, opt, arg 
  INTEGER :: narg, iarg, i, is
  LOGICAL :: ex, set_imtype


!..default
  inpcube="??"
  outcube="??"
  val_old_="all"
  val_new_="none"
  xmin=1
  ymin=1
  zmin=1
  xmax=-1
  ymax=-1
  zmax=-1

  put_NaN=.false.
  replace_NaN=.false.
  replace_all=.false.

  CALL GetArg(0,ExeName)
  narg=iargc()
  IF(narg<1) THEN
     print *, " "
     WRITE(*,'(2a)')"        CubeReplace (part of CubEx package)   "
     WRITE(*,'(a)') " "
     WRITE(*,'(a)') " by Sebastiano Cantalupo (cantalupo@phys.ethz.ch)"
     print *, " "
     WRITE(*,'(a)') "usage: CubeReplace -cube <name> -out <name> [-option <val>]"
     WRITE(*,'(a)') " "
     WRITE(*,'(a)') "  options:"
     WRITE(*,'(a)') "  -cube               <string>          : input datacube file name (NO DEFAULT)"
     WRITE(*,'(a)') "  -out                <string>          : output file name (NO DEFAULT)"
     WRITE(*,'(a)') "  -oldval             <real/string>     : values to be replaced, use 'N' for NaNs, if not selected, all values will be replaced (DEFAULT)"
     WRITE(*,'(a)') "  -newval             <real/string>     : new values, use 'N' for NaNs (NO DEFAULT)"
     WRITE(*,'(a)') "  -xmin               <int>             : xmin of the region for the replacement (default=1)"
     WRITE(*,'(a)') "  -xmax               <int>             : xmax of the region for the replacement (default=maximum x dimension)"
     WRITE(*,'(a)') "  -ymin               <int>             : ymin of the region for the replacement (default=1)"
     WRITE(*,'(a)') "  -ymax               <int>             : ymax of the region for the replacement (default=maximum y dimension)"
     WRITE(*,'(a)') "  -zmin               <int>             : ymin of the region for the replacement (default=1)"
     WRITE(*,'(a)') "  -zmax               <int>             : ymax of the region for the replacement (default=maximum y dimension)"
  ELSE
     !..read options from command line
     DO i=1,narg,2
        CALL getarg(i,opt)
        CALL getarg(i+1,arg)
        SELECT CASE(TRIM(opt))
        CASE('-cube')          ; READ(arg,'(a)') inpcube
        CASE('-out')           ; READ(arg,'(a)') outcube
        CASE('-oldval')        ; READ(arg,'(a)') val_old_
        CASE('-newval')        ; READ(arg,'(a)') val_new_
        CASE('-xmin')          ; READ(arg,*) xmin
        CASE('-xmax')          ; READ(arg,*) xmax
        CASE('-ymin')          ; READ(arg,*) ymin
        CASE('-ymax')          ; READ(arg,*) ymax
        CASE('-zmin')          ; READ(arg,*) zmin
        CASE('-zmax')          ; READ(arg,*) zmax
        CASE default
           print *, "command line argument ",TRIM(opt), " not recognized!"
           STOP
        END SELECT
     END DO
  END IF

  IF(TRIM(val_old_)=="all") THEN
     replace_all=.true.
  ELSE
     replace_all=.false.
     READ(val_old_,*,iostat=ierr) val_old
     IF(ierr/=0) THEN
        !..check if string contains "NaN"
        IF(INDEX(val_old_,"N")==0) THEN
           print *, "argument:",TRIM(val_old_),"not recognized! Please use a real number or 'N' "
           STOP
        ELSE
           replace_NaN=.true.
        END IF
     END IF
  END IF

  IF(TRIM(val_new_)=="none") STOP "Please provide the new values with the -newval option!"
  
  ierr=0
  READ(val_new_,*,iostat=ierr) val_new
  IF(ierr/=0) THEN
     !..check if string contains "NaN"
     IF(INDEX(val_new_,"N")==0) THEN
        print *, "argument:",TRIM(val_new_),"not recognized! Please use a real number or 'N' "
        STOP
     ELSE
        put_NaN=.true.
     END IF
  END IF

  IF(replace_all) THEN

     IF(.not.put_NaN) THEN
        print *, "replacing all values in selected region with ", val_new
     ELSE
        print *, "replacing all values in selected regions with NaNs"
     END IF

  ELSE
  
     IF(.not.put_NaN) THEN
        IF(.not.replace_NaN) THEN
           print *, "substituing", val_old, "with", val_new, " in file ", TRIM(inpcube), " new file name= ", TRIM(outcube)
        ELSE
           print *, "replacing NaNs with", val_new, "in file ", TRIM(inpcube), " new file name= ", TRIM(outcube)
        END IF
     ELSE
        print *, "substituing", val_old, "with NaN  in file ", TRIM(inpcube), " new file name= ", TRIM(outcube)
     END IF

  END IF


END SUBROUTINE ReadCommandLine


  
END PROGRAM CubeReplace

  

!---------------------------
! Shows and/or edit keywords in
! fits header for Cubes/images
!
!
PROGRAM CubeHeader

  IMPLICIT NONE
  CHARACTER(len=500) :: inpfile, outfile, keyword, keyval, keycomment, show, edit(2), stype
  CHARACTER(len=80)  :: comment
  LOGICAL :: showall
  INTEGER(kind=4) :: i, ierr, itype, decimals
  INTEGER :: status, inunit, or_blocksize, nkeys, keysadd, rwmode
  REAL(kind=4)    :: rtype
  REAL(kind=4)    :: dtype

  CALL ReadCommandLine

!..read header from input file

  status=0
  !..get an unused unit for the input file
  CALL ftgiou(inunit,status)
  IF(status/=0) STOP "problem with ftgiou"

  !..open input file 
  IF(TRIM(edit(1))=="??") THEN
     rwmode=0
  ELSE
     rwmode=1
  END IF       
  CALL ftnopn(inunit, InpFile, rwmode, status)
  IF(status/=0) THEN
     print *, "problem opening inpfile: ", TRIM(inpfile)
     STOP
  END IF

  IF(showall) THEN  !..simply go through the header and print all information
     
     !..get number of existing keywords in current header
     CALL ftghsp(inunit, nkeys, keysadd, status)

     !..go through header
     DO i=1,nkeys
        CALL ftgkyn(inunit, i, keyword, keyval, keycomment, status)
        print *, TRIM(keyword), " = ", TRIM(keyval), " / ", TRIM(keycomment)
     END DO

  ELSE

     IF(TRIM(show)/="??") THEN

        CALL ftgkey(inunit, show, keyval, keycomment, status)
        print *, TRIM(keyval)

     END IF
     

     IF(TRIM(edit(1))/="??") THEN

        !..get keyword type
        ierr=0
        !..try to read as integer
        READ(edit(2),*,iostat=ierr) itype
        IF(ierr==0) THEN !..write integer
           CALL ftukyj(inunit,edit(1), itype, comment, status)
           print *, "changing/adding integer keyword: ", TRIM(edit(1)), itype
        ELSE
           ierr=0
           !..try to read as real
           READ(edit(2),*,iostat=ierr) rtype
           IF(ierr==0) THEN !..write real
              IF(abs(rtype)>1.e5) THEN
                 decimals=8
              ELSE
                 decimals=-8
              END IF
              CALL ftukye(inunit,edit(1), rtype, decimals, comment, status)
              print *, "changing/adding real keyword: ", TRIM(edit(1)), rtype
           ELSE
              !..try to read as double
              ierr=0
              READ(edit(2),*,iostat=ierr) dtype
              IF(ierr==0) THEN !..write double
                 IF(abs(dtype)>1.d5) THEN
                    decimals=8
                 ELSE
                    decimals=-8
                 END IF
                 CALL ftukyd(inunit, edit(1), dtype, decimals, comment, status)
                 print *, "changing/adding double keyword: ", TRIM(edit(1)), dtype
              ELSE
                 !..try to read as string
                 ierr=0
                 READ(edit(2),'(a)',iostat=ierr) stype
                 IF(ierr==0) THEN !..write string
                    CALL ftukys(inunit, edit(1), stype, comment, status)
                    print *, "changing/adding string keyword: ", TRIM(edit(1)), " ", TRIM(stype)
                 ELSE
                    print *, "datatype of keywords to edit not recognized: ", TRIM(edit(2))
                    STOP
                 END IF
              END IF
           END IF
        END IF

     END IF

  END IF

  CALL ftclos(inunit, status)

CONTAINS

  SUBROUTINE ReadCommandLine

    IMPLICIT NONE
    CHARACTER(len=500) :: ExeName, string, ds, opt, arg , fname, longstring
    INTEGER :: narg, iarg, i, is, ierr
    LOGICAL :: ex

!..default
    inpfile="??"
    outfile="??"
    showall=.false.
    show="??"
    edit(1)="??"
    edit(2)="??"
    comment=" "

    CALL GetArg(0,ExeName)
    narg=iargc()
    IF(narg==0) THEN
       print *, " "
       WRITE(*,'(2a)')"        CubeHeader (part of CubEx package)   "
       WRITE(*,'(a)') "   Print and/or modify keywords in fits header  "
       WRITE(*,'(a)') " "
       WRITE(*,'(a)') " by Sebastiano Cantalupo (cantalupo@phys.ethz.ch)"
       print *, " "
       WRITE(*,'(a)')"usage: CubeHeader <fitsfile>   to simply show the header, use [] in filename to select extension"
       WRITE(*,'(a)')"  or "
       WRITE(*,'(a)')"usage: CubeHeader -cube <name> [-option <val>]"
       WRITE(*,'(a)')" "
       WRITE(*,'(a)')" options:"
       WRITE(*,'(a)')"   -cube            <name>           : fits file name "
       WRITE(*,'(a)')"   -show            <string>         : if provided, print the value of the selected keyword"
       WRITE(*,'(a)')"   -edit            <string>         : if provided, change/add the value of the selected keyword "
       WRITE(*,'(a)')"                                       to the value given by the -val option below. Datatype is determined automatically"
       WRITE(*,'(a)')"   -val             <string>         : value of the keyword to change/add "
       WRITE(*,'(a)')"   -comment         <string>         : comment of the keyword to change/add "
       STOP
    ELSEIF(narg==1) THEN
       CALL GetArg(1,inpfile)
       showall=.true.
       RETURN
    END IF

      !..read from command line
  DO i=1,narg,2
     CALL getarg(i,opt)
     CALL getarg(i+1,arg)
     SELECT CASE(TRIM(opt))
     CASE('-InpFile')       ; READ(arg,'(a)') InpFile
     CASE('-cube')          ; READ(arg,'(a)') InpFile
     CASE('-show')          ; READ(arg,'(a)') show
     CASE('-edit')          ; READ(arg,'(a)') edit(1)
     CASE('-val')           ; READ(arg,'(a)') edit(2)
     CASE('-comment')       ; READ(arg,'(a)') comment
      CASE default
        print *, "command line argument ",TRIM(opt), " not recognized!"
        STOP
     END SELECT
  END DO

!..perform few checks
  IF(TRIM(InpFile)=="??") THEN
     print *, "please provide the input datacube with the -InpFile or -cube option!"
     STOP
  END IF

  IF(TRIM(edit(1))/="??".and.TRIM(edit(2))=="??") THEN
     STOP "Please provide the value of the keyword to edit in the -val option!"
  END IF

  IF(TRIM(show)=="??".and.TRIM(edit(1))=="??") showall=.true.


  END SUBROUTINE ReadCommandLine


END PROGRAM CubeHeader

PROGRAM CubEx

  USE GlobalModule
  IMPLICIT NONE
  INTEGER(kind=4)    ::  st_in(8), st1(8), st2(8)
  CHARACTER(len=300) :: date, time, zone
  REAL(kind=4) :: t1, t2, NPositive
  REAL, PARAMETER  :: st_conv(8)=[0.,0.,86400.,0.,3600.,60.,1.,0.001]

!--- CPU_Time and System Time initialization
  CALL report_time("init")

  !.. read input parameter file
  CALL ReadParameterFile

  !.. read/create input datacube
  CALL ReadInputFile
  !CALL CreateCube

  !..if a previous IdCube is provided performs only photometry
  !..and write catalogue
  IF(TRIM(IdCube)/="??") THEN
     CALL UseIdCube
     CALL report_time("end")
     IF(Verbosity>=1) print *, "end of Job :-)"
     STOP
  END IF
     
  !..apply filter if requested
  IF(ApplyFilter) THEN
     CALL report_time("in")
     CALL Filter
     CALL report_time("end")
  END IF

  !..extraction and minimal detection
  CALL report_time("in")
  CALL Extract
  CALL report_time("end")

  !..if requested, perform deblending
  IF(Deblend_) THEN
     CALL report_time("in")
     CALL Deblend
     CALL report_time("end")
  END IF

  !..perform photometry
  CALL report_time("in")
  CALL Photometry
  !IF(Verbosity>=2) print *, "NObj=",COUNT(Obj(:)%Id>0)
  CALL report_time("end")

  !..associate multiple detection in z, if requested
  IF(AssocFrac>0.) THEN
     CALL report_time("in")
     CALL Associate
     CALL report_time("end")
  END IF

  !..write check cube if requested
  IF(NCheckCubes>0) THEN
     CALL report_time("in")
     CALL WriteCheckCube
     CALL report_time("end")
  END IF

  IF(VERBOSITY>=2) print *, "Writing catalogue:", TRIM(Catalogue)
  CALL WriteCatalogue(Catalogue)

  !..check reliability of the catalogue, if requested
  IF(ReliabCheck) THEN
     
     CALL report_time("in")

     !..gather number of objects found in the previous catalogue
     NPositive=COUNT(Obj(:)%Id>0)

     IF(Verbosity>=2) print *, "Performing extraction and photometry checks on 'negative' cube..."

     WHERE(Cube/=UNDEF) Cube=-Cube
     IF(ApplyFilter) THEN
        WHERE(CubeF/=UNDEF) CubeF=-CubeF
     END IF
     CALL Extract
     CALL Photometry

     !..print information
     print *, " "
     print *, "The 'blind' reliability factor of this catalogue is:",(1.-COUNT(Obj(:)%Id>0)/REAL(NPositive))*100.,"%"
     print *, " "

     IF(TRIM(NegCatalogue)/="??") CALL WriteCatalogue(NegCatalogue)
        

     CALL report_time("end")

  END IF

  IF(AssocFrac>0.and.TRIM(AssocCatalogue)/="??") THEN 

     IF(Verbosity>=2) THEN
        print *, " "
        print *, "Producing merged catalogue with associated objects: ", TRIM(AssocCatalogue)
        print *, " "
     END IF

     !..merge associated objects and update mask
     CALL MergeAssoc

     !..performs photometry on merged list
     CALL Photometry
     IF(Verbosity>=2) THEN
        print *, " "
        print *, "FINAL NObj in merged catalogue=",COUNT(Obj(:)%Id>0)
        print *, " "
     END IF

     !..write catalogue
     CALL WriteCatalogue(AssocCatalogue)

  END IF

  CALL report_time("end")

  IF(Verbosity>=1) print *, "end of Job :-)"

CONTAINS

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

END PROGRAM CubEx

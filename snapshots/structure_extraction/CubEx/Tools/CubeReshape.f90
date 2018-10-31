PROGRAM CubeReshape

  USE CubeLib
  IMPLICIT NONE
  CHARACTER(len=500) :: inpcube, outcube, arg, opt
  REAL(kind=4) :: npsize, opsize, npos(2), opos(2), dist, mindist
  INTEGER(kind=4) :: oDimX, oDimY, oDimZ, nDimX, nDimY, nDimZ, i, j, k, ii, jj, thispix(2), rfact(2), narg
  REAL, ALLOCATABLE :: nCube(:,:,:)
  LOGICAL :: sm

  narg=iargc()
  IF(narg<1) THEN
     print *, "usage: CubeReshape -cube <inpcube> -out <outcube> [options]"
     print *, "options:"
     print *, "-rx <int> : rescaling factor in the x direction (default=1)"
     print *, "-ry <int> : rescaling factor in the y direcion (default=1)"
     print *, "-sm <.true./.false.> : if .true. apply gaussian smoothing with radius         "
     print *, "                       (rx-1) and/or (ry-1) after rescaling if rx and/or ry/=1"
     STOP
  END IF

  !...defaults
  rfact=1
  sm=.false.

  !..read options from command line
  DO i=1,narg,2
     CALL getarg(i,opt)
     CALL getarg(i+1,arg)
     SELECT CASE(TRIM(opt))
     CASE('-cube')          ; READ(arg,'(a)') inpcube
     CASE('-out')           ; READ(arg,'(a)') outcube
     CASE('-rx')            ; READ(arg,*) rfact(1)
     CASE('-ry')            ; READ(arg,*) rfact(2)
     CASE('-sm')            ; READ(arg,*) sm
     END SELECT
  END DO

  IF(ALL(rfact/=1)) STOP "rescaling on multiple axes not implemented yet!"
  IF(ALL(rfact==1)) STOP "select rescaling factor with the options -rx or -ry"
  
  print *, "reading: ", TRIM(inpcube)
  CALL ReadCube(inpcube)
  oDimX=SIZE(Cube,Dim=1)
  oDimY=SIZE(Cube,Dim=2)
  oDimZ=SIZE(Cube,Dim=3)

  print *, "old cube/image size: ", oDimX, oDimY, oDimZ

  !..get new sizes
  nDimX=oDimX*rfact(1)
  nDimY=oDimY*rfact(2)
  nDimZ=oDimZ
  print *, "new cube/image size: ", nDimX, nDimY, nDimZ

  ALLOCATE(nCube(nDimX,nDimY,nDimZ))
  
  IF(rfact(1)/=1) THEN
     DO i=1,oDimX
        ii=(i-1)*rfact(1)+1
        DO j=0,rfact(1)-1,1
           nCube(ii+j,:,:)=Cube(i,:,:)
        END DO
     END DO
  ELSEIF(rfact(2)/=1) THEN
     DO i=1,oDimY
        ii=(i-1)*rfact(2)+1
        DO j=0,rfact(2)-1,1
           nCube(:,ii+j,:)=Cube(:,i,:)
        END DO
     END DO
  END IF

  !..write cube
  print *, "writing: ", TRIM(outcube)
  CALL WriteLocalCube(nCube, outcube)

END PROGRAM CubeReshape

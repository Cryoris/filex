SUBROUTINE UseIdCube

  USE Globalmodule
  USE StatLib
  IMPLICIT NONE
  REAL(kind=4), ALLOCATABLE :: IdMask(:,:,:), IdMask_2D(:,:)
  INTEGER(kind=4) :: unit, status, rank, naxes(3), ival, nfound, i, IdMin, IdMax, group, in, end, &
       j, k, nval, val(1:1000), ierr, x, y
  CHARACTER(len=500) :: comment
  LOGICAL :: anyf
  REAL(kind=4) :: thisflux, thiserr, apr_2


  status=0

  ! ---- read IdCube

  IF(VERBOSITY>=1) print *, "Reading IdCube= ", TRIM(IdCube)
  

  !..get an unused unit
  CALL ftgiou(unit,status)
  IF(status/=0) STOP "problem with ftgiou"

  !..open file in read-only mode
  CALL ftdopn(unit,IdCube,0,status)
  IF(status/=0) STOP "problem with ftgiou"

  !..get cube rank and dimensions
  status=0
  CALL ftgkyj(unit,'NAXIS',rank,comment,status)
  IF(status/=0) STOP "problem reading NAXIS keyword"

  status=0
  IF(rank==3) THEN
     !..get cube size
     CALL ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
  ELSEIF(rank==2) THEN
     !..get image size
     CALL ftgkyj(unit,'NAXIS1',ival,comment,status)
     naxes(1)=ival
     CALL ftgkyj(unit,'NAXIS2',ival,comment,status)
     naxes(2)=ival
     naxes(3)=1
  ELSE
     print *, "input file is not a cube or a image according to NAXIS keyword:", rank
     STOP
  END IF

  !..check cube dimension
  IF(naxes(1)/=DimX.or.naxes(2)/=DimY) THEN
     print *, TRIM(IdCube), " and ", TRIM(InpFile), " have different spatial dimensions!"
     print *, naxes(1:2), " vs ", DimX, DimY
     STOP
  END IF
  IF(naxes(3)/=DimZ) THEN
     print *, "WARNING: ", TRIM(IdCube)," and ", TRIM(InpFile), " have different z-dimensions"
     print *, "            IsoFlux measurement will be performed only in the spatial dimension"
     print *, naxes(3), " vs ", DimZ
  ENDIF

  ALLOCATE(IdMask(naxes(1),naxes(2),naxes(3)),IdMask_2D(naxes(1),naxes(2)))

  !..read data
  group=1
  in=1  
  end=PRODUCT(naxes)
  CALL ftgpve(unit,group,in,end,UNDEF,IdMask,anyf,status)
  IF(status/=0) STOP "problem reading data in IdCube"
  
  CALL ftclos(unit,status)
  IF(status/=0) STOP "problem closing fits file"

  !..free all allocated units
  call ftfiou(-1,status)
  IF(status/=0) STOP "problem with ftfiou"


  !--- get min and max Id values
  IdMin=INT(MINVAL(IdMask,MASK=IdMask/=UNDEF.and.IdMask/=0))
  IdMax=INT(MAXVAL(IdMask))
  
  IF(Verbosity>=1) print *, "original minmax IdCube=", IdMin, IdMax

  !..allocate Obj array
  ALLOCATE(Obj(IdMin:IdMax))
  Obj(:)%Id=0

  !..read InpCat
  CALL ReadInpCat

  IF(TRIM(IdCubeOnlyList)/="-1") THEN
     !..finds how many entry there are in the string
     DO i=1,1000
        ierr=0
        READ(IdCubeOnlyList,*,iostat=ierr) val(1:i)
        IF(ierr/=0) THEN
           nval=i-1
           EXIT
        END IF
     END DO
     IF(i>=1000) STOP "problem reading IdCubeOnly option!"
     !..read entries and apply values
     READ(IdCubeOnlyList,*) val(1:nval)
     print *, "performing photometry only on the following objects: ", val(1:nval)
     DO i=IdMin, IdMax
        IF(ALL(val(1:nval)/=Obj(i)%Id)) Obj(i)%Id=0
     END DO

  END IF

  !..get aperRadius square for later
  apr_2=AperRadius**2

  !--- perform photometry and write catalogue
  print *, "performing photometry..."
  OPEN(1,file=Catalogue,action="write")
  IF(printheader) THEN
     WRITE(1,'(a)') "###################################################"
     WRITE(1,'(a)') "# col 1:     Id "
     WRITE(1,'(a)') "# col 2:     IsoFlux"
     WRITE(1,'(a)') "# col 3:     Err_IsoFlux"
     WRITE(1,'(a)') "# col 4:     AperFlux"
     WRITE(1,'(a)') "# col 5:     Err_AperFlux"
     WRITE(1,'(a)') "# col 6:     aperture radius (pixel)"
     WRITE(1,'(a)') "# col 7:     zmin (pixel)"
     WRITE(1,'(a)') "# col 8:     zmax (pixel)"
  END IF
  DO i=IdMin,IdMax

     IF(Obj(i)%Id==0) CYCLE

     !..get IsoFlux and Err
     IF(naxes(3)==DimZ) THEN !..cube and mask have same dimensions
        thisflux=SUM(Cube(1:DimX,1:DimY,1:DimZ),MASK=(Cube(1:DimX,1:DimY,1:DimZ)/=UNDEF.and.IdMask==i))
        thiserr=sqrt(SUM(Var(1:DimX,1:DimY,1:DimZ),MASK=(Var(1:DimX,1:DimY,1:DimZ)/=UNDEF.and.IdMask==i)))
     ELSE !..sum spaxel by spaxel instead
        thisflux=0
        thiserr=0
        DO y=Obj(i)%boxmin(2),Obj(i)%boxmax(2)
           DO x=Obj(i)%boxmin(1),Obj(i)%boxmax(1)
              IF(ANY(IdMask(x,y,:)==i)) THEN
                 thisflux=thisflux+SUM(Cube(x,y,1:DimZ),MASK=Cube(x,y,1:DimZ)/=UNDEF)
                 thiserr=thiserr+SUM(Var(x,y,1:DimZ),MASK=Var(x,y,1:DimZ)/=UNDEF)
              END IF
           END DO
        END DO
        thiserr=sqrt(thiserr)
     END IF

     !..get aperture photometry for this object
     CALL this_AperPhotometry(this_ID=i)

     !..write into catalogue
     WRITE(1,*) i,thisflux,thiserr, Obj(i)%AperFlux, Obj(i)%AperErr, AperRadius, Obj(i)%boxmin(3), Obj(i)%boxmax(3)

  END DO
  CLOSE(1)

  print *, " "
  print *, "output catalogue: ", TRIM(Catalogue)
  
CONTAINS

!---------------------------------------- 

  SUBROUTINE ReadInpCat

    IMPLICIT NONE
    INTEGER :: ierr, id, idum, ii
    REAL :: rdum
    CHARACTER(len=10) :: string

    print *, "reading InpCat: ", TRIM(InpCat)
    OPEN(11,file=InpCat,action="read")
    ierr=0
    DO
       READ(11,*,IOSTAT=ierr) string
       IF(ierr/=0) EXIT
       IF(SCAN(TRIM(string),'!#%$')>0) CYCLE 
       BACKSPACE(11)
       READ(11,*) id
       BACKSPACE(11)
       READ(11,*) Obj(id)%Id, idum, rdum, rdum, rdum, Obj(id)%lcen(1:3), Obj(id)%boxmin(1:3), Obj(id)%boxmax(1:3)
    END DO

    IF(naxes(3)/=DimZ) THEN !..change z-box values to the full wavelenght range of current cube
       Obj(:)%boxmin(3)=1
       Obj(:)%boxmax(3)=DimZ
    END IF


  END SUBROUTINE ReadInpCat

!-----------------------------------------------------

  SUBROUTINE this_AperPhotometry(this_ID)

    IMPLICIT NONE
    INTEGER(kind=4), INTENT(IN) :: this_ID
    REAL(kind=4)    :: lcen, sum_w, dist, pos(2), oneD(1:DimZ), meanclip, medianclip, this_sigma
    INTEGER(kind=4) :: x, y, z, i, xmin(3), xmax(3)
    LOGICAL :: oneD_mask(1:DimZ)

    i=this_ID

    !..recover original bounding box
    xmin=Obj(i)%boxmin
    xmax=Obj(i)%boxmax

    ! -------- perform cylindrical aperture photometry centered on the light-weighted centroid
    !.. NB: correction for fractional area of pixels TO BE IMPLEMENTED

    IF(naxes(3)==DimZ) THEN !..if the cube and IdMask have different dimensions we skip this part
                            !..and use the full wavelength range instead
   
       !..find new bounding box
       IF(AperDz<0) THEN !..get automatic Dz from a 1D spectrum

          !..produce a 2D mask first
          DO y=xmin(2),xmax(2)
             DO x=xmin(1),xmax(1)
                IF(ANY(IdMask(x,y,:)==this_ID)) THEN
                   IdMask_2D(x,y)=this_ID
                ELSE
                   IdMask_2D(x,y)=0
                END IF
             END DO
          END DO

          !..get optimally extracted 1d spectrum using the 2D mask
          DO z=1,DimZ
             oneD(z)=SUM(Cube(xmin(1):xmax(1),xmin(2):xmax(2),z),MASK=Cube(xmin(1):xmax(1),xmin(2):xmax(2),z)/=UNDEF.and.&
                  IdMask_2D(xmin(1):xmax(1),xmin(2):xmax(2))==this_ID)       
          END DO

          !..get noise of the spectrum (excluding line)
          oneD_mask=.true.
          oneD_mask(xmin(3):xmax(3))=.false.
          CALL  SigmaClip(PACK(oneD,MASK=oneD_mask), MeanClip, MedianClip, this_sigma)

          !..find new bounding box
          DO z=xmin(3),1,-1; IF(oneD(z)<abs(AperDz)*this_sigma+meanclip) EXIT; END DO; xmin(3)=MAX(1,z)
          DO z=xmax(3),DimZ; IF(oneD(z)<abs(AperDz)*this_sigma+meanclip) EXIT; END DO; xmax(3)=MIN(z,DimZ)
                
       ELSEIF(AperDz>0) THEN
                
          xmin(3)=MAX(1,NINT(Obj(i)%lcen(3)-AperDz/2.))
          xmax(3)=MIN(DimZ,NINT(Obj(i)%lcen(3)+AperDz/2.))
                
       END IF !..if AperDz==0, use the original boxmin and boxmax values in z
       
    ENDIF  

    !..adjust bounding box in space
    xmin(1:2)=MAX(1,NINT(Obj(i)%lcen(1:2)-AperRadius))
    xmax(1)=MIN(DimX,NINT(Obj(i)%lcen(1)+AperRadius))
    xmax(2)=MIN(DimY,NINT(Obj(i)%lcen(2)+AperRadius))

    !..perform aperture photometry
    Obj(i)%AperFlux=0.
    DO y=xmin(2),xmax(2)
       DO x=xmin(1),xmax(1)
          pos=[x-0.5,y-0.5]
          dist=SUM((pos-Obj(i)%lcen(1:2))**2)
          IF(dist<apr_2) THEN
             Obj(i)%AperFlux=Obj(i)%AperFlux+SUM(Cube(x,y,xmin(3):xmax(3)),MASK=(Cube(x,y,xmin(3):xmax(3))/=UNDEF))
             Obj(i)%AperErr=Obj(i)%AperFlux+SUM(Var(x,y,xmin(3):xmax(3)),MASK=(Cube(x,y,xmin(3):xmax(3))/=UNDEF).and.Var(x,y,xmin(3):xmax(3))/=UNDEF)
          END IF
       END DO
    END DO
    IF(Obj(i)%AperErr>0) THEN
       Obj(i)%AperErr=sqrt(Obj(i)%AperErr)
    ELSE
       Obj(i)%AperErr=0.
    END IF

    !..update photometry bounding box in z
    Obj(i)%boxmin(3)=xmin(3)
    Obj(i)%boxmax(3)=xmax(3)


  END SUBROUTINE this_AperPhotometry

END SUBROUTINE UseIdCube

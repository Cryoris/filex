!-----------------------------------------------------
!
! performs simple photometry and additional detection checks
! using the 3D mask defined in the previous steps.
!
! Author: SC
! Last Mod: Jun 2nd, 2016
!
!-----------------------------------------------------

SUBROUTINE Photometry

  USE Globalmodule
  USE StatLib
  IMPLICIT NONE
  INTEGER(kind=4) :: nobj, i, xmax(3), xmin(3), area, dz, x, y, z, removed, kept, this_ID
  REAL(kind=4)    :: sum_w, pos(2), apr_2, dist, err
  INTEGER(kind=4), ALLOCATABLE :: IdMask_2D(:,:)
  REAL(kind=4)    :: oneD(1:DimZ), meanclip, medianclip, this_sigma
  LOGICAL         :: oneD_mask(1:DimZ)


  nobj=SIZE(Obj)

  apr_2=AperRadius**2

  IF(VERBOSITY>=1) print *, "Performing photometry and additional detection checks..."

  removed=0
  kept=0

  IF(SSN_Threshold>0) ALLOCATE(IdMask_2D(0:DimX+1,0:DimY+1)) !..same spatial dimensions of Mask

  DO i=1,nobj
     
     IF(Obj(i)%Id<=0) CYCLE

     xmin=Obj(i)%boxmin
     xmax=Obj(i)%boxmax

     !..check area and dz
     Obj(i)%dz=xmax(3)-xmin(3)+1
     Obj(i)%area=COUNT(SUM(Mask(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3)),&
          MASK=(Mask(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))==IdToLabel(Obj(i)%Id)),DIM=3)/=0)

     Obj(i)%IsoFlux=SUM(Cube(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3)),&
          MASK=(Mask(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))==IdToLabel(Obj(i)%Id)).and.&
          Cube(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))/=UNDEF)

     err=SUM(Var(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3)),&
          MASK=(Mask(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))==IdToLabel(Obj(i)%Id)).and.&
          Var(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))/=UNDEF.and.&
          Cube(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))/=UNDEF)

     IF(err/=0.) THEN
        Obj(i)%IsoErr=sqrt(err)
     ELSE
        Obj(i)%IsoErr=0.
     END IF


     IF(Obj(i)%IsoFlux<Obj(i)%IsoErr*ISN_Threshold.or.Obj(i)%dz<MinDZ.or.Obj(i)%dz>MaxDz.or.Obj(i)%area<MinArea) THEN
        Obj(i)%Id=-Obj(i)%Id   !..flag this object as non-detection
        IF(TRIM(CheckCube(1))/="??") THEN
           WHERE(Mask==IdToLabel(-Obj(i)%Id)) Mask=0   !..remove object from mask for checkcube, if requested
        END IF
        removed=removed+1
        CYCLE
     END IF

     IF(SSN_Threshold>0.) THEN !..perform spectral SN Threshold check

        this_ID=Obj(i)%Id

        !..produce a 2D mask first
        DO y=xmin(2),xmax(2)
           DO x=xmin(1),xmax(1)
              IF(ANY(Mask(x,y,:)==IdToLabel(this_ID))) THEN
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

        !..readjust the noise considering the number of z-layers (does not include co-variance)
        this_sigma=this_sigma*sqrt(xmax(3)-xmin(3)+1.)
          

        !IF(Obj(i)%Id==42) THEN
        !   print *, xmin, xmax, SUM(oneD(xmin(3):xmax(3))), MeanClip*(xmax(3)-xmin(3)+1), this_sigma
        !END IF

        !..compare the flux and the noise
        IF(SUM(oneD(xmin(3):xmax(3)))<this_sigma*SSN_Threshold) THEN             
           Obj(i)%Id=-Obj(i)%Id   !..flag this object as non-detection
           IF(TRIM(CheckCube(1))/="??") THEN
              WHERE(Mask==IdToLabel(-Obj(i)%Id)) Mask=0   !..remove object from mask for checkcube, if requested
           END IF
           removed=removed+1
           CYCLE
        END IF
        
     END IF

     kept=kept+1

     Obj(i)%BoxFlux=SUM(Cube(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3)), &
          MASK=(Cube(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))/=UNDEF))

     err=SUM(Var(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3)), &
          MASK=(Cube(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))/=UNDEF).and.&
          Var(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))/=UNDEF)
                    
     IF(err>0.) THEN
        Obj(i)%BoxErr=sqrt(err)
     ELSE
        Obj(i)%BoxErr=0.
     END IF

     !..compute light-weighted centroid
     Obj(i)%lcen=0.
     sum_w=0.
     DO z=xmin(3),xmax(3)
        DO y=xmin(2),xmax(2)
          DO x=xmin(1),xmax(1)

             IF(Mask(x,y,z)/=IdToLabel(Obj(i)%Id).or.Cube(x,y,z)==UNDEF) CYCLE

             Obj(i)%lcen(:)=Obj(i)%lcen(:)+([x,y,z]-0.5)*Cube(x,y,z)
             sum_w=sum_w+Cube(x,y,z)

          END DO
       END DO
    END DO
    Obj(i)%lcen=Obj(i)%lcen/sum_w

    ! -------- perform cylindrical aperture photometry centered on the light-weighted centroid
    !.. NB: correction for fractional area of pixels TO BE IMPLEMENTED

    !..find bounding box
    xmin(1:2)=MAX(1,NINT(Obj(i)%lcen(1:2)-AperRadius))
    xmax(1)=MIN(DimX,NINT(Obj(i)%lcen(1)+AperRadius))
    xmax(2)=MIN(DimY,NINT(Obj(i)%lcen(2)+AperRadius))
    xmin(3)=MAX(1,NINT(Obj(i)%lcen(3)-AperDz/2.))
    xmax(3)=MIN(DimZ,NINT(Obj(i)%lcen(3)+AperDz/2.))

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

 END DO

 !IF(kept==0) STOP "All objects removed after additional checks! I stop here..."

 IF(VERBOSITY>=2) THEN
    print *, " removed:", removed, " objects after photometry"
    print *, " remaining objects;", kept
 END IF
             

END SUBROUTINE Photometry

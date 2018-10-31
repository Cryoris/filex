SUBROUTINE Filter

  USE Globalmodule
  IMPLICIT NONE
  INTEGER, PARAMETER :: ns=3          !.. this parameter determines the filtering mask size (=2*ns*FilterRad)
  INTEGER :: nmin, nmax, dd, ix, iy, iz, i, dd_z, xmin, x1, x2, y1, y2, z1, z2, j, k
  REAL(kind=4) :: S
  REAL(kind=4), ALLOCATABLE :: w(:), tmp(:,:,:), w_norm(:)

  IF(Verbosity>=2) THEN
     print *, " "
     print *, "Applying Gaussian Filter with spatial size (radius)=", FilterXYRad, " and spectral size (radius)=", FilterZRad, "pixels"
  END IF

  dd=ns*FilterXYRad
  dd_z=ns*FilterZRad

  IF(Verbosity>=2) print *, "Allocating arrays..."

!..allocate cubes taking into account padding for filtering (and ghost zones)
  ALLOCATE(CubeF(-dd:DimX+MAX(dd,1),-dd:DimY+MAX(dd,1),-dd_z:DimZ+MAX(dd_z,1)),&
       VarF(-dd:DimX+MAX(dd,1),-dd:DimY+MAX(dd,1),-dd_z:DimZ+MAX(dd_z,1)),&
       tmp(-dd:DimX+MAX(dd,1),-dd:DimY+MAX(dd,1),-dd_z:DimZ+MAX(dd_z,1)))
  CubeF=UNDEF
  CubeF(1:DimX,1:DimY,1:DimZ)=Cube(1:DimX,1:DimY,1:DimZ)  !..use this as an initial temporal array (with the right size)
  VarF=UNDEF
  VarF(1:DimX,1:DimY,1:DimZ)=Var(1:DimX,1:DimY,1:DimZ)   !..use this as an initial temporal array (with the right size)
  tmp=UNDEF

  IF(Verbosity>=2) print *, "done"

  IF(Verbosity>=2) print *, "Filtering Cube..."

!--- CUBE ------------

  !..first perform spatial filtering, if requested

  IF(FilterXYRad>0) THEN

     nmin=-ns*FilterXYRad
     nmax=ns*FilterXYRad
     
     IF(ALLOCATED(w)) DEALLOCATE(w,w_norm)
     ALLOCATE(w(nmin:nmax),w_norm(nmin:nmax))

     !..generate 1D filter  
     DO i=nmin,nmax,1
        w(i)=exp(-(REAL(i)**2/(2.*FilterXYRad**2)))
     END DO
     !..normalized version
     w_norm(:)=w(:)/SUM(w(:))

     !..smooth first in x
     DO iz=1,DimZ
        DO iy=1,DimY
           DO ix=1,DimX
              
              IF(ALL(CubeF(ix-dd:ix+dd,iy,iz)==UNDEF)) CYCLE !..no data to filter

              IF(ANY(CubeF(ix-dd:ix+dd,iy,iz)==UNDEF)) THEN !..we can't use previously normalized weight here
                 tmp(ix,iy,iz)=SUM(CubeF(ix-dd:ix+dd,iy,iz)*w(:),MASK=CubeF(ix-dd:ix+dd,iy,iz)/=UNDEF)/&
                      SUM(w(:),MASK=CubeF(ix-dd:ix+dd,iy,iz)/=UNDEF)
              ELSE
                 tmp(ix,iy,iz)=SUM(CubeF(ix-dd:ix+dd,iy,iz)*w_norm(:))
              END IF

           END DO
        END DO
     END DO
  
     !..smooth now in y
     DO iz=1,DimZ
        DO iy=1,DimY
           DO ix=1,DimX
                            
              IF(ALL(tmp(ix,iy-dd:iy+dd,iz)==UNDEF)) THEN
                 CubeF(ix,iy,iz)=UNDEF
                 CYCLE
              END IF

              IF(ANY(tmp(ix,iy-dd:iy+dd,iz)==UNDEF)) THEN !..we can't use normalized weight here
                 CubeF(ix,iy,iz)=SUM(tmp(ix,iy-dd:iy+dd,iz)*w(:),MASK=tmp(ix,iy-dd:iy+dd,iz)/=UNDEF)/&
                      SUM(w(:),MASK=tmp(ix,iy-dd:iy+dd,iz)/=UNDEF)
              ELSE
                 CubeF(ix,iy,iz)=SUM(tmp(ix,iy-dd:iy+dd,iz)*w_norm(:))
              END IF


           END DO
        END DO
     END DO

  END IF

  !..now performs spectral filtering, if requested
  IF(FilterZRad>0) THEN
  
     !..reassign temporary array
     tmp=CubeF

     nmin=-ns*FilterZRad
     nmax=ns*FilterZRad

     IF(ALLOCATED(w)) DEALLOCATE(w,w_norm)
     ALLOCATE(w(nmin:nmax),w_norm(nmin:nmax))
     !..generate 1D filter (weighted) 
     DO i=nmin,nmax,1
        w(i)=exp(-(REAL(i)**2/(2.*FilterZRad**2)))
     END DO
     !..normalized version
     w_norm=w/SUM(w)

     !..smooth in z: NB: this simple implementation is computationally expensive for large cubes (it runs on the "wrong" index for Cube array)
     DO iz=1,DimZ
        DO iy=1,DimY
           DO ix=1,DimX
          
              IF(ALL(tmp(ix,iy,iz-dd_z:iz+dd_z)==UNDEF)) THEN
                 CubeF(ix,iy,iz)=UNDEF
                 CYCLE 
              END IF

              IF(ANY(tmp(ix,iy,iz-dd_z:iz+dd_z)==UNDEF)) THEN 
                 CubeF(ix,iy,iz)=SUM(tmp(ix,iy,iz-dd_z:iz+dd_z)*w(:),MASK=tmp(ix,iy,iz-dd_z:iz+dd_z)/=UNDEF)/&
                      SUM(w(:),MASK=tmp(ix,iy,iz-dd_z:iz+dd_z)/=UNDEF)
              ELSE
                 CubeF(ix,iy,iz)=SUM(tmp(ix,iy,iz-dd_z:iz+dd_z)*w_norm(:))
              END IF

           END DO
        END DO
     END DO

  END IF

  IF(Verbosity>=2) print *, "done"

  IF(.not.ApplyFilterVar) RETURN

!------ VAR --------------------

  IF(Verbosity>=2) print *, "Filtering Variance Cube..."

  !..reset tmp array
  tmp=UNDEF

  !..first perform spatial filtering, if requested

  IF(FilterXYRad>0) THEN

     nmin=-ns*FilterXYRad
     nmax=ns*FilterXYRad
     
     IF(ALLOCATED(w)) DEALLOCATE(w,w_norm)
     ALLOCATE(w(nmin:nmax),w_norm(nmin:nmax))

     !..generate 1D filter  
     DO i=nmin,nmax,1
        w(i)=exp(-(REAL(i)**2/(2.*FilterXYRad**2)))
     END DO
     !..normalized version
     w_norm(:)=w(:)/SUM(w(:))

     !..smooth first in x
     DO iz=1,DimZ
        DO iy=1,DimY
           DO ix=1,DimX
            
              IF(ALL(VarF(ix-dd:ix+dd,iy,iz)==UNDEF)) CYCLE !..no data to filter

              IF(ANY(VarF(ix-dd:ix+dd,iy,iz)==UNDEF)) THEN !..we can't use previously normalized weight here
                 tmp(ix,iy,iz)=SUM(VarF(ix-dd:ix+dd,iy,iz)*w(:)**2,MASK=VarF(ix-dd:ix+dd,iy,iz)/=UNDEF)/&
                      (SUM(w(:),MASK=VarF(ix-dd:ix+dd,iy,iz)/=UNDEF)**2)
              ELSE
                 tmp(ix,iy,iz)=SUM(VarF(ix-dd:ix+dd,iy,iz)*w_norm(:)**2)
              END IF

           END DO
        END DO
     END DO
  
     !..smooth now in y
     DO iz=1,DimZ
        DO iy=1,DimY
           DO ix=1,DimX
             
              IF(ALL(tmp(ix,iy-dd:iy+dd,iz)==UNDEF)) THEN
                 VarF(ix,iy,iz)=UNDEF
                 CYCLE
              END IF
              
              IF(ANY(tmp(ix,iy-dd:iy+dd,iz)==UNDEF)) THEN !..we can't use normalized weight here
                 VarF(ix,iy,iz)=SUM(tmp(ix,iy-dd:iy+dd,iz)*w(:)**2,MASK=tmp(ix,iy-dd:iy+dd,iz)/=UNDEF)/&
                      (SUM(w(:),MASK=tmp(ix,iy-dd:iy+dd,iz)/=UNDEF)**2)
              ELSE
                 VarF(ix,iy,iz)=SUM(tmp(ix,iy-dd:iy+dd,iz)*w_norm(:)**2)
              END IF

           END DO
        END DO
     END DO

  END IF

  !..now performs spectral filtering, if requested
  IF(FilterZRad>0) THEN

     !..reassign temporary array
     tmp=VarF
  
     nmin=-ns*FilterZRad
     nmax=ns*FilterZRad


     IF(ALLOCATED(w)) DEALLOCATE(w,w_norm)
     ALLOCATE(w(nmin:nmax),w_norm(nmin:nmax))
     !..generate 1D filter (weighted) 
     DO i=nmin,nmax,1
        w(i)=exp(-(REAL(i)**2/(2.*FilterZRad**2)))
     END DO
     !..normalized version
     w_norm=w/SUM(w)

     !..smooth in z: NB: this simple implementation is computationally expensive for large cubes (it runs on the "wrong" index for Cube array)
     DO iz=1,DimZ
        DO iy=1,DimY
           DO ix=1,DimX
              
              IF(ALL(tmp(ix,iy,iz-dd_z:iz+dd_z)==UNDEF)) THEN
                 VarF(ix,iy,iz)=UNDEF
                 CYCLE 
              END IF

              IF(ANY(tmp(ix,iy,iz-dd_z:iz+dd_z)==UNDEF)) THEN 
                 VarF(ix,iy,iz)=SUM(tmp(ix,iy,iz-dd_z:iz+dd_z)*w(:)**2,MASK=tmp(ix,iy,iz-dd_z:iz+dd_z)/=UNDEF)/&
                      (SUM(w(:),MASK=tmp(ix,iy,iz-dd_z:iz+dd_z)/=UNDEF)**2)
              ELSE
                 VarF(ix,iy,iz)=SUM(tmp(ix,iy,iz-dd_z:iz+dd_z)*(w_norm(:)**2))
              END IF

           END DO
        END DO
     END DO

  END IF

  IF(Verbosity>=2) print *, "done."

  DEALLOCATE(tmp)

END SUBROUTINE Filter

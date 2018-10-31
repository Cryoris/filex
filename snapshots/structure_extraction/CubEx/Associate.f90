!----------------------------------------------
!
! associates object in the z-direction based on 
! the fraction of the overlapping area. Objects
! are sorted by area before associations are done.
!
! Author:SC
! Last Mod: Oct 16th, 2014
!
!---------------------------------------------

SUBROUTINE Associate

  USE Globalmodule
  USE StatLib
  IMPLICIT NONE
  INTEGER :: nobj, i, j, k, area, xmin(3), xmax(3), il, jl, overlap
  INTEGER(kind=4), ALLOCATABLE :: AssocIm_main(:,:), AssocIm_sec(:,:), order_list(:)
  REAL(kind=4), ALLOCATABLE :: obj_size(:)
  REAL(kind=4) :: frac

  IF(VERBOSITY>=2) print *, "Associating objects..."

!..create auxiliary lists
  nobj=SIZE(Obj)
  ALLOCATE(order_list(nobj),obj_size(nobj))
  DO i=1,nobj
!     IF(Obj(i)%Id<=0) CYCLE
     order_list(i)=i 
     obj_size(i)=Obj(i)%area
     Obj(i)%Assoc=0
  END DO

!..sort objects by area 
  CALL quick_sort(list=obj_size, order=order_list)

!..allocate auxiliary arrays
  ALLOCATE(AssocIm_main(DimX,DimY),AssocIm_sec(DimX,DimY))
  AssocIm_main=0
  AssocIm_sec=0

!..loop starting from largest objects
  DO il=nobj,1,-1

     i=order_list(il)

     IF(Obj(i)%Id<=0) CYCLE

     IF(Obj(i)%Assoc/=0) CYCLE !..this object was already associated to a larger one

     xmin=Obj(i)%boxmin
     xmax=Obj(i)%boxmax


     AssocIm_main(xmin(1):xmax(1),xmin(2):xmax(2))=&
          SUM(mask(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3)),&
          MASK=(mask(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))==IdToLabel(Obj(i)%Id)),DIM=3)

     DO jl=il-1,1,-1 !..loop through detections and check overlapping area

        j=order_list(jl)

        IF(Obj(j)%Id<=0) CYCLE

        xmin=Obj(j)%boxmin
        xmax=Obj(j)%boxmax

        AssocIm_sec(xmin(1):xmax(1),xmin(2):xmax(2))=&
          SUM(mask(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3)),&
          MASK=(mask(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))==IdToLabel(Obj(j)%Id)),DIM=3)

        area=COUNT(AssocIm_sec(xmin(1):xmax(1),xmin(2):xmax(2))/=0)

        overlap=COUNT((AssocIm_main(xmin(1):xmax(1),xmin(2):xmax(2))*AssocIm_sec(xmin(1):xmax(1),xmin(2):xmax(2)))/=0)

        IF(REAL(overlap,KIND=4)>=AssocFrac*REAL(area)) THEN
           Obj(j)%Assoc=Obj(i)%Id
           !..flag also main object (with his own Id)
           Obj(i)%Assoc=Obj(i)%Id
        END IF

        !..reset AssocIm_sec for next object
        AssocIm_sec(xmin(1):xmax(1),xmin(2):xmax(2))=0

     END DO

     !..reset AssocIm_main for next object
     xmin(1:2)=Obj(i)%boxmin(1:2)
     xmax(1:2)=Obj(i)%boxmax(1:2)
     AssocIm_main(xmin(1):xmax(1),xmin(2):xmax(2))=0


  END DO

      
END SUBROUTINE Associate

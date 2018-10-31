!-------------------------------------------
!
! merges associates object (see routine Associate) and
! updates mask.
!
! Author: SC
! Last mod: Oct 16, 2014
!
!-------------------------------------------

SUBROUTINE MergeAssoc

  USE GlobalModule
  USE StatLib
  IMPLICIT NONE
  INTEGER :: nobj, i, j, k, il, jl, xmin_main(3), xmax_main(3), xmin_sec(3), xmax_sec(3), ii
  INTEGER(kind=4), ALLOCATABLE :: list(:)
  REAL(kind=4), ALLOCATABLE :: obj_size(:)


  nobj=SIZE(Obj)

!..create auxiliary lists
  nobj=SIZE(Obj)
  ALLOCATE(list(nobj),obj_size(nobj))
  DO i=1,nobj
     list(i)=i 
     obj_size(i)=Obj(i)%area
  END DO

!..sort objects by area 
  CALL quick_sort(obj_size, order=list)

!..loop starting from largest objects
  DO il=nobj,1,-1

     i=list(il)

     IF(Obj(i)%Id<=0) CYCLE

     IF(Obj(i)%Assoc/=Obj(i)%Id) CYCLE !..either without Assoc or Associated with other object


     DO jl=il-1,1,-1 !..loop trough smaller objects to collect assoc

        j=list(jl)

        IF(Obj(j)%Assoc==Obj(i)%Id) THEN

           xmin_sec=Obj(j)%boxmin
           xmax_sec=Obj(j)%boxmax

           !..update boxmin and boxmax
           DO ii=1,3
              Obj(i)%boxmin(ii)=MIN(Obj(i)%boxmin(ii),xmin_sec(ii))
              Obj(i)%boxmax(ii)=MAX(Obj(i)%boxmax(ii),xmax_sec(ii))
           END DO

           !..merge mask
           WHERE(Mask(xmin_sec(1):xmax_sec(1),xmin_sec(2):xmax_sec(2),xmin_sec(3):xmax_sec(3))==IdToLabel(Obj(j)%Id)) &
                Mask(xmin_sec(1):xmax_sec(1),xmin_sec(2):xmax_sec(2),xmin_sec(3):xmax_sec(3))=IdToLabel(Obj(i)%Id)

           !..remove object from catalogue
           Obj(j)%Id=-Obj(j)%Id

        END IF

     END DO

  END DO

END SUBROUTINE MergeAssoc

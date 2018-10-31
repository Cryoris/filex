!---------------------------------------------------------------
!
! performs deblending using a progressive SN threshold with min, 
! max and number of steps defined by the user (Deblend_* parameters).
! Final ISO volume for each deblended object will correspond to the 
! lowest SN threshold (within selected range) at which the object may be split.
! Deblending is performed on a 2D mask defined by collapsing along the 
! z-direction the the minimum volume occupyied by the object in the datacube.
!
! Author:SC
! Last modification: Oct 16th, 2014
!
!----------------------------------------------------------------

SUBROUTINE Deblend

  USE Globalmodule
  IMPLICIT NONE
  REAL(kind=4), ALLOCATABLE :: SN2d(:,:) 
  INTEGER(kind=4), ALLOCATABLE :: mask2d(:,:), DBmask(:,:)
  INTEGER(kind=4) :: i, xmin(3), xmax(3), xp, yp, deb_nsteps, db, o, j, nblend, this_id
  REAL(kind=4) :: deb_min, deb_max, deb_step, this_deb
  INTEGER(kind=4) :: nobj_orig, id_new, nsplit, orig_Id, label_new, ii, k, id
  INTEGER(kind=4), ALLOCATABLE :: parent(:), NSpax_(:), DBtoId(:)

  IF(VERBOSITY>=2) print *, "Deblending..."

!..allocate local arrays and label (these are deallocated on exit)
  ALLOCATE(parent(Deblend_BufferSize), NSpax_(Deblend_BufferSize), DBtoId(Deblend_BufferSize))

  deb_min=Deblend_MinSNR*SN_Threshold
  deb_max=Deblend_MaxSNR*SN_Threshold
  deb_nsteps=Deblend_NSteps
  deb_step=(deb_max-deb_min)/deb_nsteps

  !..allocate auxiliary mask
  ALLOCATE(mask2d(SIZE(mask,DIM=1),SIZE(mask,DIM=2)),&
       SN2d(SIZE(mask,DIM=1),SIZE(mask,DIM=2)),&
       DBmask(SIZE(mask,DIM=1),SIZE(mask,DIM=2)))
  mask2d=0.
  DBmask=0.
  SN2d=0.

  nobj_orig=COUNT(Obj(:)%Id>0)
  label_new=MAXVAL(Mask)
  id_new=nobj_orig

  !..loop over detected objects to check for blending
  DO o=1,nobj_orig
     
     !..original extent of this object
     xmin=Obj(o)%boxmin
     xmax=Obj(o)%boxmax

     !..create SN image of this object
     IF(ApplyFilter) THEN
        SN2d(xmin(1):xmax(1),xmin(2):xmax(2))=&
             SUM(CubeF(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3)),&
             MASK=(mask(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))==IdToLabel(Obj(o)%Id)),DIM=3)/&
             sqrt(SUM(VarF(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3)),&
             MASK=(mask(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))==IdToLabel(Obj(o)%Id)),DIM=3))
     ELSE
        SN2d(xmin(1):xmax(1),xmin(2):xmax(2))=&
             SUM(Cube(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3)),&
             MASK=(mask(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))==IdToLabel(Obj(o)%Id)),DIM=3)/&
             sqrt(SUM(Var(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3)),&
             MASK=(mask(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3))==IdToLabel(Obj(o)%Id)),DIM=3))
     END IF
     
     !..assign a temporary Id of 1 to the DBmask as a starting point
     WHERE(SN2d/=0.) DBmask=1
     this_id=1

     !..loop over deb thresholds
     deb_loop:DO db=1,deb_nsteps+1
 
        this_deb=deb_min+(db-1)*deb_step

        !..get number of possibly blended object and associated mask2d, DBmask
        !..looping over previously deblended objects
        DO this_id=1,MAXVAL(DBmask(xmin(1):xmax(1),xmin(2):xmax(2)))
           CALL GetNBlend(this_deb, this_id)
        END DO

     END DO deb_loop

     nsplit=MAXVAL(DBmask(xmin(1):xmax(1),xmin(2):xmax(2)))

     !---- DEBLENDING -----------------

     IF(nsplit>1) THEN

        !..create new objects and update auxiliary lists
        DO i=1,nsplit

           id_new=id_new+1
           label_new=label_new+1

           IF(id_new>SIZE(Obj).or.label_new>SIZE(LabelToId)) STOP "Increase Deblend Buffer Size (option Deblend_Buffer)!"

           Obj(id_new)%Id=id_new
           Obj(id_new)%Assoc=0
           Obj(id_new)%NSpax=0

           IdToLabel(id_new)=label_new
           LabelToId(label_new)=id_new

           DBtoId(i)=id_new

           !..initialize bounding boxes values
           DO ii=1,3
              Obj(id_new)%boxmin(ii)=100000
              Obj(id_new)%boxmax(ii)=-1
              Obj(id_new)%xcen(ii)=0
           END DO

        END DO

        !..loop over spaxel of the original object and change global Mask
        orig_Id=Obj(o)%Id

        DO k=xmin(3),xmax(3)
           DO j=xmin(2),xmax(2)
              DO i=xmin(1),xmax(1)

                 IF(Mask(i,j,k)/=IdToLabel(orig_Id)) CYCLE

                 IF(DBMask(i,j)==0) THEN !..this spaxel does not belong to any object now
                    
                    !..remove spaxel from mask
                    Mask(i,j,k)=0

                 ELSE

                    !..associate to a split object with this id:
                    id=DBtoId(DBMask(i,j))
                    Obj(id)%NSpax=Obj(id)%NSpax+1

                    !..update mask with a new label associated to this object
                    Mask(i,j,k)=IdToLabel(id)

                    !..calculate new xcen and boxmin/max:
                    Obj(id)%xcen(:)=Obj(id)%xcen(:)+[i,j,k]-0.5
                
                    Obj(id)%boxmin(1)=MIN(Obj(id)%boxmin(1),i)
                    Obj(id)%boxmin(2)=MIN(Obj(id)%boxmin(2),j)
                    Obj(id)%boxmin(3)=MIN(Obj(id)%boxmin(3),k)
                
                    Obj(id)%boxmax(1)=MAX(Obj(id)%boxmax(1),i)
                    Obj(id)%boxmax(2)=MAX(Obj(id)%boxmax(2),j)
                    Obj(id)%boxmax(3)=MAX(Obj(id)%boxmax(3),k)

                 END IF

              END DO
           END DO
        END DO

        !..remove splitted object from original list
        Obj(orig_Id)%Id=-orig_Id
              
     END IF
     
     !..cleanup SN2d and DBMask for next object
     SN2d(xmin(1):xmax(1),xmin(2):xmax(2))=0.
     DBmask(xmin(1):xmax(1),xmin(2):xmax(2))=0

  END DO

  !..finalize centroid calculation
  DO ii=1,3
     Obj(nobj_orig+1:id_new)%xcen(ii)=Obj(nobj_orig+1:id_new)%xcen(ii)/Obj(nobj_orig+1:id_new)%NSpax
  END DO

!..deallocate local arrays
  DEALLOCATE(parent, NSpax_, mask2d, SN2d)

IF(VERBOSITY>=2) THEN
   print *, "NObj after deblending=", COUNT(Obj(:)%Id>0)    
END IF

!-----------------------------
           
CONTAINS


!-----------------------------

  SUBROUTINE GetNBlend(this_deb, this_id)

    IMPLICIT NONE
    REAL(kind=4), INTENT(IN)     :: this_deb
    INTEGER(kind=4), INTENT(IN)  :: this_id
    INTEGER(kind=4) :: label, prior_labels(9), p, this_label, this_NSpax, ii, x, y, z, id, nlabels, ndet, prev_nblend, init_label, nblend                    
    

    !..initialize arrays
    parent=0
    prior_labels=0

    !..set this to the largest label already present in the mask (if any)
    init_label=MAXVAL(DBmask)
    label=init_label

    !..reset mask2d for this object
    mask2d(xmin(1):xmax(1),xmin(2):xmax(2))=0

    DO j=xmin(2),xmax(2)
       DO i=xmin(1),xmax(1)

          IF(DBmask(i,j)/=this_id) CYCLE

          IF(SN2d(i,j)>this_deb) THEN   !..flag pixel

             !..check (prior) neighbors of this pixel, for simplicity we actually check ALL neighbors here  
             prior_labels=RESHAPE(mask2d(i-1:i+1,j-1:j+1),(/9/))
           
             IF(ALL(prior_labels==0)) THEN   !..new component --> new label
                label=label+1
                IF(label>Deblend_BufferSize) STOP "Increase Deblend buffer stack size (option Deblend_Buffer)!"
                mask2d(i,j)=label                
             ELSE !..this spaxel is connected to another one
                this_label=MINVAL(prior_labels,MASK=prior_labels/=0) 
                mask2d(i,j)=this_label
                !..update parent tree 
                DO p=1,SIZE(prior_labels)
                   IF(prior_labels(p)/=0.and.prior_labels(p)/=this_label) THEN
                      CALL union(this_label, prior_labels(p))
                   END IF
                END DO
             END IF

          END IF

       END DO
    END DO
    
    IF(label<=init_label+1) THEN !..only one single connected object, no need to continue...
       nblend=1
       RETURN
    END IF

    nlabels=MAXVAL(mask2d(xmin(1):xmax(1),xmin(2):xmax(2)))

    NSpax_=0
    !..second pass:
    !... replace labels using the parent tree
    !... get NSpax for each individual connected component
    DO j=xmin(2),xmax(2)
       DO i=xmin(1),xmax(1)

          this_label=mask2d(i,j)

          IF(this_label/=0) THEN
 
             !..assign value from parent tree
             p=this_label
             DO WHILE(parent(p)/=0)
                p=parent(p)
             END DO
             mask2d(i,j)=p
             
             !..update NSpax counter associated with this label
             NSpax_(p)=NSpax_(p)+1

          END IF

       END DO
    END DO

    !..save number of previously detected objects
    prev_nblend=MAXVAL(DBmask(xmin(1):xmax(1),xmin(2):xmax(2)))
    nblend=0

    !..cleanup and count objecs above Deblend_MinNPix
    DO i=init_label,nlabels
       
       IF(parent(i)==0) THEN  
       
          this_label=i
          this_NSpax=NSpax_(this_label)
          
          IF(this_NSpax>=Deblend_MinNPix) THEN
             
             nblend=nblend+1               ! update ndet
            
             IF(nblend==1) THEN !..keep the same id as original object
                WHERE(mask2d(xmin(1):xmax(1),xmin(2):xmax(2))==this_label) &
                     mask2d(xmin(1):xmax(1),xmin(2):xmax(2))=this_id
             ELSE !..create a new id
                WHERE(mask2d(xmin(1):xmax(1),xmin(2):xmax(2))==this_label) &
                     mask2d(xmin(1):xmax(1),xmin(2):xmax(2))=prev_nblend+nblend-1                
             END IF

          END IF
       
       END IF
    END DO

    !..if we have more than one object, split them updating DBmask
    IF(nblend>1) THEN

       !..clean previous object
       WHERE(DBmask(xmin(1):xmax(1),xmin(2):xmax(2))==this_id) DBmask(xmin(1):xmax(1),xmin(2):xmax(2))=0

       !..assign new masks 
       WHERE(mask2d(xmin(1):xmax(1),xmin(2):xmax(2))==this_id) DBmask(xmin(1):xmax(1),xmin(2):xmax(2))=this_id
       DO i=2,nblend
          WHERE(mask2d(xmin(1):xmax(1),xmin(2):xmax(2))==prev_nblend+i-1) DBmask(xmin(1):xmax(1),xmin(2):xmax(2))=prev_nblend+i-1
       END DO

    END IF


  END SUBROUTINE GetNBlend


!------------------------------------------------

  SUBROUTINE union(x,y)

    IMPLICIT NONE
    INTEGER(kind=4) :: x, y

    !..find root labels
    DO WHILE(parent(x)/=0)
       x=parent(x)
    END DO
    DO WHILE(parent(y)/=0)
       y=parent(y)
    END DO

    IF(x<y) THEN
       parent(y)=x
    ELSEIF(x>y) THEN
       parent(x)=y
    END IF

  END SUBROUTINE union


END SUBROUTINE Deblend

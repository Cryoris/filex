!-----------------------------------------------
!
! Compute (extract) connected components in the datacube with a 3d extension 
! of the row-by-row and  Union-Find algorithms (e.g. Shapiro & Stockman, Computer Vision, Mar 2000). 
! This routine performs also minimal detection and fills in the non-photometric 
! fields in Obj derived type array.
!
! Author: SC
! Last Mod: Oct 16, 2014
!-----------------------------------------------


SUBROUTINE Extract

 USE Globalmodule
 IMPLICIT NONE
 INTEGER(kind=4) :: i, j, k, label, prior_labels(27), p, this_label, nobj, this_NSpax, ii, x, y, z, id, nlabels, ndet
 INTEGER(kind=4), ALLOCATABLE :: parent(:), NSpax_(:)

!..initialize local arrays and label (these are deallocated on exit)
 ALLOCATE(parent(maxnlabels), NSpax_(maxnlabels))
 parent=0
 label=0
 Mask=0

 IF(Verbosity>=2) print *, "Extracting objects (first pass)..."

 !..first pass:
 !... assign labels using 8-pixel-connectivity, that becomes 27-spaxel-connectivity in 3D)
 !... build "parent" tree with union-find structure


 IF(.not.ApplyFilter) THEN


    !..make sure that undefined variance pixels that are defined in the datacube have a dummy high value
    WHERE(Cube/=UNDEF.and.(Var==UNDEF.or.Var==0.)) Var=1.e30

    DO k=1,DimZ
       DO j=1,DimY
          DO i=1,DimX

             IF(Cube(i,j,k)==UNDEF) CYCLE

             IF(Cube(i,j,k)>SN_Threshold*sqrt(Var(i,j,k))) THEN   !..flagged spaxel

                !..check (prior) neighbors of this spaxel, for simplicity we actually check ALL neighbors here  
                prior_labels=RESHAPE(mask(i-1:i+1,j-1:j+1,k-1:k+1),(/27/))

                IF(ALL(prior_labels==0)) THEN   !..new component --> new label
                   label=label+1
                   IF(label>maxnlabels) STOP "Increase stack size (maxnlabels)!"
                   mask(i,j,k)=label                
                ELSE !..this spaxel is connected to another one
                   this_label=MINVAL(prior_labels,MASK=prior_labels/=0) 
                   mask(i,j,k)=this_label
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
    END DO

 ELSE
    print *, "Help im in FCube!"

    !..make sure that undefined variance pixels that are defined in the datacube have a dummy high value
    WHERE(CubeF/=UNDEF.and.(VarF==UNDEF.or.VarF==0.)) VarF=1.e30

    DO k=1,DimZ
       DO j=1,DimY
          DO i=1,DimX

             IF(CubeF(i,j,k)==UNDEF) CYCLE

             IF(CubeF(i,j,k)>SN_Threshold*sqrt(VarF(i,j,k))) THEN   !..flagged spaxel

                !..check (prior) neighbors of this spaxel, for simplicity we actually check ALL neighbors here  
                prior_labels=RESHAPE(mask(i-1:i+1,j-1:j+1,k-1:k+1),(/27/))

                IF(ALL(prior_labels==0)) THEN   !..new component --> new label
                   label=label+1
                   IF(label>maxnlabels) STOP "Increase stack size (maxnlabels)!"
                   mask(i,j,k)=label                
                ELSE !..this spaxel is connected to another one
                   this_label=MINVAL(prior_labels,MASK=prior_labels/=0) 
                   mask(i,j,k)=this_label
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
    END DO

 END IF

 IF(SN_Threshold_Conn>0.and.SN_Threshold_Conn<SN_Threshold) THEN
    print *, "connecting at lower SN Threhold..."
    CALL Connect_at_lower_SNT(ApplyFilter)
    print *, "done"
 END IF

 !..mark all pixels > Halo thres as haloes (JG)
 HaloLabel=MAXVAL(mask)+1
 DO k=1,DimZ
    DO j=1,DimY
       DO i=1,DimX
        IF(Cube(i,j,k)>Halo_Threshold*sqrt(Var(i,j,k))) THEN
          mask(i,j,k)=HaloLabel
        END IF
      END DO
    END DO
  END DO

 !..one more due to haloes
 nlabels=MAXVAL(mask)
 IF(Verbosity>=3) print *, "numer of labels=",nlabels
 IF(Verbosity>=2) print *, "Extracting objects (second pass)..."

 NSpax_=0
!..second pass:
!... replace labels using the parent tree
!... get NSpax for each individual connected component
 DO k=1,DimZ
    DO j=1,DimY
       DO i=1,DimX

          this_label=mask(i,j,k)

          IF(this_label==HaloLabel) THEN
            NSpax_(this_label)=NSpax_(this_label)+1

          ELSE IF(this_label/=0) THEN
 
             !..assign value from parent tree
             p=this_label
             DO WHILE(parent(p)/=0)
                p=parent(p)
             END DO
             mask(i,j,k)=p

             !..update NSpax counter associated with this label
             NSpax_(p)=NSpax_(p)+1

          END IF

       END DO
    END DO
 END DO

 !..this is the number of individual connected components found in the cube:
 nobj=COUNT(parent(1:nlabels)==0)
 IF(Verbosity>=1) print *, "NObj Extracted=",nobj

 !...create auxiliary arrays
 IF(ALLOCATED(LabelToId)) DEALLOCATE(LabelToId,IdToLabel)

 IF(.not.Deblend_) THEN
    ALLOCATE(LabelToId(nlabels),IdToLabel(nobj))
 ELSE
    ALLOCATE(LabelToId(nlabels+Deblend_BufferSize),IdToLabel(nobj+Deblend_BufferSize))
 END IF
 LabelToId=0
 IdToLabel=0

 !----- DETECTION (using NSpax) -------------

 !..build auxiliary arrays and count detections
 ndet=0               
 DO i=1,nlabels

    IF(parent(i)==0) THEN
       
       this_label=i
       this_NSpax=NSpax_(this_label)
       
       IF(this_NSpax>MinNSpax) THEN

          ndet=ndet+1               ! update ndet
 
          IdToLabel(ndet)=this_label
          LabelToId(this_label)=ndet
         
       END IF
       
    END IF
 END DO

 IF(Verbosity>=1) print *, "NObj Detected=", Ndet

 ! --- Create object list

 IF(ALLOCATED(Obj)) DEALLOCATE(Obj)
 IF(.not.Deblend_) THEN
    ALLOCATE(Obj(ndet))
 ELSE
    ALLOCATE(Obj(ndet+Deblend_BufferSize))
    Obj(:)%Id=0
 END IF
 

 !..initialize Id, Assoc and fill in NSpax 
 DO i=1,ndet
    Obj(i)%Id=i
    Obj(i)%Assoc=0
    Obj(i)%NSpax=NSpax_(IdToLabel(i))
 END DO

 !..initialize bounding boxes values
 DO ii=1,3
    Obj(:)%boxmin(ii)=100000
    Obj(:)%boxmax(ii)=-1
    Obj(:)%xcen(ii)=0
 END DO

 !..find bounding boxes and centroid for each objects
 DO k=1,DimZ
    DO j=1,DimY
       DO i=1,DimX

          this_label=Mask(i,j,k)

          IF(this_label/=0) THEN

             id=LabelToId(this_label) !..get object associated with pixel

             IF(id/=0) THEN
             
                Obj(id)%xcen(:)=Obj(id)%xcen(:)+[i,j,k]-0.5
                
                Obj(id)%boxmin(1)=MIN(Obj(id)%boxmin(1),i)
                Obj(id)%boxmin(2)=MIN(Obj(id)%boxmin(2),j)
                Obj(id)%boxmin(3)=MIN(Obj(id)%boxmin(3),k)
                
                Obj(id)%boxmax(1)=MAX(Obj(id)%boxmax(1),i)
                Obj(id)%boxmax(2)=MAX(Obj(id)%boxmax(2),j)
                Obj(id)%boxmax(3)=MAX(Obj(id)%boxmax(3),k)         
                
             ELSE  !..cleanup Mask

                Mask(i,j,k)=0

             END IF
                

          END IF

       END DO
    END DO
 END DO

 !..finalize geometrical centroid calculation
 DO ii=1,3
    Obj(1:ndet)%xcen(ii)=Obj(1:ndet)%xcen(ii)/Obj(1:ndet)%NSpax
 END DO
    
 DEALLOCATE(parent,NSpax_)

CONTAINS

!--------------------------------------------

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

!----------------------------------------------

  SUBROUTINE Connect_at_lower_SNT(ApplyFilter_)

    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: ApplyFilter_
    REAL, ALLOCATABLE :: SNR(:,:,:)
    INTEGER :: nfound, it

    ALLOCATE(SNR(0:DimX+1,0:DimY+1,0:DimZ+1))
    SNR=UNDEF
    

    print *, "producing SNR array..."
    IF(.not.ApplyFilter_) THEN

       WHERE(Cube(1:DimX,1:DimY,1:DimZ)/=UNDEF.and.&
            Var(1:DimX,1:DimY,1:DimZ)/=UNDEF.and.&
            Var(1:DimX,1:DimY,1:DimZ)>=0.)
          SNR(1:DimX,1:DimY,1:DimZ)=Cube(1:DimX,1:DimY,1:DimZ)/sqrt(Var(1:DimX,1:DimY,1:DimZ))
       END WHERE

    ELSE

       WHERE(CubeF(1:DimX,1:DimY,1:DimZ)/=UNDEF.or.&
            VarF(1:DimX,1:DimY,1:DimZ)/=UNDEF.or.&
            VarF(1:DimX,1:DimY,1:DimZ)>=0.)
          SNR(1:DimX,1:DimY,1:DimZ)=CubeF(1:DimX,1:DimY,1:DimZ)/sqrt(VarF(1:DimX,1:DimY,1:DimZ))
       END WHERE

    END IF
    print *, "done"

    it=0
    iter_loop: DO
    
       nfound=0
       DO k=1,DimZ
          DO j=1,DimY
             DO i=1,DimX

                IF(SNR(i,j,k)==UNDEF.or.mask(i,j,k)/=0) CYCLE
                
                IF(SNR(i,j,k)>SN_Threshold_Conn) THEN   !..check connection

                   !..check (prior) neighbors of this spaxel, for simplicity we actually check ALL neighbors here  
                   prior_labels=RESHAPE(mask(i-1:i+1,j-1:j+1,k-1:k+1),(/27/))
                   
                   IF(ANY(prior_labels/=0)) THEN  
                      
                      !IF(ANY(SNR(i-1:i+1,j-1:j+1,k-1:k+1)>SN_Threshold)) THEN
                         
                         nfound=nfound+1
                         this_label=MINVAL(prior_labels,MASK=prior_labels/=0) 
                         mask(i,j,k)=this_label
                         !..update parent tree 
                         DO p=1,SIZE(prior_labels)
                            IF(prior_labels(p)/=0.and.prior_labels(p)/=this_label) THEN
                               CALL union(this_label, prior_labels(p))
                            END IF
                         END DO

                      !kEND IF
                   
                   END IF

                END IF
                   
             END DO
          END DO
       END DO

       it=it+1
       print *, "iteration", it, "new connected voxels found=",nfound
          
       IF(nfound==0) EXIT iter_loop

    END DO iter_loop

    DEALLOCATE(SNR)

  END SUBROUTINE Connect_at_lower_SNT


END SUBROUTINE EXTRACT



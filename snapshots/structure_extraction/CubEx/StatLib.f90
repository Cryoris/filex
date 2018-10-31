!------------------------------------------------------------------------
!
! A series of statistical functions and subroutines to make life easier.
!
! AUTHOR: SC
! LAST MOD: Mar 6th, 2015
!
! CONTENTS:
!
!   SUBROUTINES 
!       SigmaClip : median sigmaclipping routine
! 
!   FUNCTIONS
!       Mean
!       Median
!       StdDev
!       quick_sort
!       random_gaussian
!       ImSmooth
!
MODULE StatLib


CONTAINS

!------------------------------------------------------------
!
! SUBROUTINE SigmaClip --- Median Sigma Clipping
!
! This routine returns a sigmaclipped mean, median, and stddev of a real array
! clipping is performed around median iteratively until full convergence (desidered convergence level 
! may be adjusted with the optional parameter ConvergenceNum, see below). 
!
! AUTHOR: SC
! LAST MOD: Aug 29th, 2014
!
! INPUT: 
!        Array        (real, input array)
!
! OUTPUT: 
!        MeanClip     (real, clipped mean)
!        MedianClip   (real, clipped median)
!        FinalSigma   (real, clipped stddev)
!
! OPTIONAL INPUT:
!        ClipVal             (real(2), two-sided clipping threshold in unit of stddev, default=[-3.,3.])
!        MaxIterations       (integer, maximum number of clipping iterations, default=1000)
!        ConvergenceNum      (real,    iterate until the fraction of removed data points from previous iteration is larger than this value, default=0.)
!        NValClipThreshold   (integer, iterations will stop when the number of remaining values decreases below this threshold, default=3)
!        CRrejectFact        (real,    if present, the largest data point in the array will be removed if the ratio between 
!                                      the stddev including/excluding this data point is larger than CRrejectFact. It may be used to remove cosmic rays)
!        Threshold           (real(2), if present, all the values smaller than Threshold(1) and larger than Threshold(2) in the original array
!                                      will not be considered)
!        weight              (real array with weights (same size as Array), if present performs a weighted clipped mean) 
!
! OPTIONAL OUTPUT:
!        TotIterations       (integer, returns the total number of iterations performed)
!        TotGoodVals         (integer, returns the number of remaining values in the array after clipping)
!        GoodMask            (logical array, .false. for clipped pixels)
!
! 
!

RECURSIVE SUBROUTINE SigmaClip(array, MeanClip, MedianClip, FinalSigma, ClipVal, MaxIterations, ConvergenceNum, NValClipThreshold, TotIterations, TotGoodVals, CRrejectFact, Threshold, Verbosity, weight, GoodMask, VarArray, PropVar, UNDEF_val)

  IMPLICIT NONE
  REAL, DIMENSION(:), INTENT(IN) :: array
  REAL, INTENT(OUT) :: MeanClip, MedianClip, FinalSigma
  REAL, OPTIONAL, INTENT(IN) :: ClipVal(2), ConvergenceNum, CRrejectFact, Threshold(2), weight(SIZE(array))
  REAL, OPTIONAL, INTENT(IN) :: VarArray(SIZE(array)), UNDEF_val
  REAL, OPTIONAL, INTENT(OUT) :: PropVar
  INTEGER, OPTIONAL, INTENT(IN) :: MaxIterations, Verbosity, NValClipThreshold
  INTEGER, OPTIONAL, INTENT(OUT) :: TotIterations, TotGoodVals
  LOGICAL, OPTIONAL, INTENT(OUT) :: GoodMask(SIZE(array))
  !..local variables
  REAL :: a(SIZE(array)), UNDEF_
  REAL :: sigclip(2), conv, thisMedian, thisMean, thisSigma, ThisConv, nmin, sorted_weight(SIZE(array)), sorted_var(SIZE(array))
  INTEGER :: maxiter, ThisIter, first, last, ngood, ngood_old, Verbosity_, order(SIZE(array)), i
  LOGICAL :: KeepMask(SIZE(array))

!..initialize local arrays
  KeepMask=.true.
  a=array

!..initiliaze default/optional values
  IF(present(ClipVal)) THEN
     sigclip=ClipVal
  ELSE
     sigclip=[-3.,3.]
  END IF

  IF(present(MaxIterations)) THEN
     maxiter=MaxIterations
  ELSE
     maxiter=1000
  END IF
  
  IF(present(ConvergenceNum)) THEN
     conv=ConvergenceNum
  ELSE
     conv=0.   ! 0.=continue until no values are removed between iterations
  END IF

  IF(present(NValClipThreshold)) THEN
     nmin=NValClipThreshold
  ELSE
     nmin=4
  END IF  

  IF(present(Verbosity)) THEN
     Verbosity_=Verbosity
  ELSE
     Verbosity_=1
  END IF

  IF(present(UNDEF_val)) THEN
     UNDEF_=UNDEF_val
  ELSE
     UNDEF_=-999.0
  END IF

!..initial sorting of the array and initial values
  IF(present(weight).or.present(VarArray)) THEN
     DO i=1,SIZE(array) 
        order(i)=i 
     END DO
     CALL quick_sort(a,order=order)
     DO i=1,SIZE(array)
        IF(present(weight)) sorted_weight(i)=weight(order(i))
        IF(present(VarArray)) sorted_var(i)=VarArray(order(i))
     END DO
  ELSE
     CALL quick_sort(a)
  END IF
          
!---- Main Iteration ----------

!..initialize variables
  ThisConv=1
  ThisIter=0
  first=1
  last=SIZE(a)
  ngood=SIZE(a)

!..if requested, remove positive large outliers from array (e.g., cosmic rays)
!..comparing sigma including/excluding largest value. If ratio is larger than CRrejectFact
!..then the last pixel is removed.
  IF(present(CRrejectFact)) THEN
      IF(StdDev(a(first:last),usemedian=.true.)/StdDev(a(first:last-1),usemedian=.true.)>CRrejectFact) THEN
         last=last-1
         ngood=ngood-1
         KeepMask(last)=.false.
      END IF
   END IF

!..if requested, apply data thresholding
   IF(present(Threshold)) THEN
      WHERE(a<Threshold(1)) KeepMask=.false.
      WHERE(a>Threshold(2)) KeepMask=.false.
      !..find new first good value starting from the left
      DO WHILE(.not.KeepMask(first))
         first=first+1
      END DO
      !..find new last good value starting from the right
      DO WHILE(.not.KeepMask(last))
         last=last-1
      END DO
      !..update number of remaining good values:
      ngood=last-first+1
      IF(ngood<1) THEN !..set everything UNDEF_val
         MeanClip=UNDEF_; MedianClip=UNDEF_; FinalSigma=UNDEF_ 
         IF(present(TotIterations)) TotIterations=0
         IF(present(TotGoodVals)) TotGoodVals=0
         IF(present(VarArray)) PropVar=UNDEF_
         RETURN
      END IF
  END IF

!..clipping loop
  iterloop: DO WHILE(ThisConv>conv)
 
     ngood_old=ngood

     ThisIter=ThisIter+1
     IF(ngood<=nmin) THEN
        IF(Verbosity_>=2) THEN
           print *, "WARNING: convergence not achieved in SigClip to the required level:", conv
           print *, "Actual convergenve level is:", ThisConv, "i.e.", ngood, "/", SIZE(a)
           print *, "Remaining values from clipping:", ngood, " ; mininum acceptable number is:", nmin
        END IF
        EXIT iterloop
     END IF

     !..midpoint is:
     IF(MOD(ngood,2)==0) THEN
        ThisMedian=(a(ngood/2+first-1)+a(ngood/2+first))*0.5
     ELSE
        ThisMedian=a(ngood/2+first)
     END IF
     
     !..new stddev
     ThisSigma=StdDev(a(first:last),usemedian=.true.)

     !..update mask
     WHERE((a(first:last)-ThisMedian)<sigclip(1)*ThisSigma) KeepMask(first:last)=.false.
     WHERE((a(first:last)-ThisMedian)>sigclip(2)*ThisSigma) KeepMask(first:last)=.false.   

     !..find new first good value starting from the left
     DO WHILE(.not.KeepMask(first))
        first=first+1
     END DO
     !..find new last good value starting from the right
     DO WHILE(.not.KeepMask(last))
        last=last-1
     END DO

     !..update number of remaining good values:
     ngood=last-first+1
    
     !..update convergence
     ThisConv=(ngood_old-ngood)/REAL(ngood_old)

  END DO iterloop

  !print *, "ngood=",ngood
  !print *, "SIZE(a)=",SIZE(a)
  
  !..final clipped median
  IF(MOD(ngood,2)==0) THEN
     MedianClip=(a(ngood/2+first-1)+a(ngood/2+first))*0.5
  ELSE
     MedianClip=a(ngood/2+first)
  END IF
     
  !..final clipped mean and stddev around mean
  IF(present(weight)) THEN
     MeanClip=Mean(a(first:last),weight=sorted_weight(first:last))
     FinalSigma=StdDev(a(first:last),usemedian=.false.,weight=sorted_weight(first:last))
  ELSE
     MeanClip=Mean(a(first:last))
     FinalSigma=StdDev(a(first:last),usemedian=.false.)
  END IF

  !..propagate variance and return value, if requested
  IF(present(VarArray)) THEN
     IF(COUNT(sorted_var(first:last)/=UNDEF_)>0) THEN
        PropVar=SUM(sorted_var(first:last),MASK=sorted_var(first:last)/=UNDEF_)/(COUNT(sorted_var(first:last)/=UNDEF_)**2)
     ELSE
        PropVar=UNDEF_
     END IF
  END IF

  !..if requested, return number of iterations performed 
  IF(present(TotIterations)) TotIterations=ThisIter

  !..if requested, return number of remaining good values
  IF(present(TotGoodVals)) TotGoodVals=ngood

  !..if requested, return KeepMask
  IF(present(GoodMask)) GoodMask=KeepMask

END SUBROUTINE SigmaClip

!------------------------------

REAL FUNCTION Mean(array, weight)

  IMPLICIT NONE
  REAL, DIMENSION(:), INTENT(IN) :: array
  REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: weight

  IF(present(weight)) THEN
     Mean=SUM(array*weight)/SUM(weight)
  ELSE
     Mean=SUM(array)/SIZE(array)
  END IF

END FUNCTION Mean

!---------------------------------

REAL FUNCTION StdDev(array,usemedian, weight)

  IMPLICIT NONE
  REAL, DIMENSION(:), INTENT(IN) :: array
  REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: weight
  LOGICAL, OPTIONAL :: usemedian
  REAL :: m

  IF(present(usemedian)) THEN
     IF(usemedian) THEN
        m=Median(array)
     ELSE
        IF(present(weight)) THEN
           m=Mean(array,weight=weight)
        ELSE
           m=Mean(array)
        END IF
     END IF
  ELSE
     IF(present(weight)) THEN
        m=Mean(array,weight=weight)
     ELSE
        m=Mean(array)
     END IF
  END IF

  IF(present(weight)) THEN
     StdDev=SQRT(SUM(((array(:)-m)*weight(:))**2))/SUM(weight)
  ELSE
     IF(SIZE(array)>1) THEN
        StdDev=SQRT(SUM((array-m)**2)/SIZE(array-1)) !..corrected sample standard deviation
     ELSE
        StdDev=SQRT(SUM((array-m)**2))
     END IF
  END IF


END FUNCTION StdDev

!----------------------------------------------
!
! Median()
! Input: real array
! Output: median
!

RECURSIVE REAL FUNCTION Median(array)

  IMPLICIT NONE
  REAL, DIMENSION(:), INTENT(IN) :: array
  REAL, DIMENSION(:), ALLOCATABLE :: array_copy
  INTEGER :: n

  n=SIZE(array)

  ALLOCATE(array_copy(n))
  array_copy=array
  CALL quick_sort(array_copy)
  IF(MOD(n,2)==0) THEN
     Median=(array_copy(n/2)+array_copy(n/2+1))*0.5
  ELSE
     Median=array_copy(n/2+1)
  END IF

END FUNCTION Median


!-----------------------------------------------
!
! QUICK SORT routine for real arrays
! from: Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! INPUT: real array (list), associated integer array with position of elements in the original order (optional)
! OUTPUT: sorted real array in ascending order
!

RECURSIVE SUBROUTINE quick_sort(list, order)

IMPLICIT NONE
REAL, DIMENSION (:), INTENT(IN OUT)  :: list
INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: order

CALL quick_sort_1(1, SIZE(list))

CONTAINS

  RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

    INTEGER, INTENT(IN) :: left_end, right_end

    !     Local variables
    INTEGER             :: i, j, itemp
    REAL                :: reference, temp
    INTEGER, PARAMETER  :: max_simple_sort_size = 6

    IF (right_end < left_end + max_simple_sort_size) THEN
       ! Use interchange sort for small lists
       CALL interchange_sort(left_end, right_end)
    ELSE
       ! Use partition ("quick") sort
       reference = list((left_end + right_end)/2)
       i = left_end - 1; j = right_end + 1
       DO
          ! Scan list from left end until element >= reference is found
          DO
             i = i + 1
             IF (list(i) >= reference) EXIT
          END DO
          ! Scan list from right end until element <= reference is found
          DO
             j = j - 1
             IF (list(j) <= reference) EXIT
          END DO


          IF (i < j) THEN
             ! Swap two out-of-order elements
             temp = list(i); list(i) = list(j); list(j) = temp
             IF(present(order)) THEN
                itemp = order(i); order(i) = order(j); order(j) = itemp
             ENDIF
          ELSE IF (i == j) THEN
             i = i + 1
             EXIT
          ELSE
             EXIT
          END IF
       END DO

       IF (left_end < j) CALL quick_sort_1(left_end, j)
       IF (i < right_end) CALL quick_sort_1(i, right_end)
    END IF
    
  END SUBROUTINE quick_sort_1


  RECURSIVE SUBROUTINE interchange_sort(left_end, right_end)

    INTEGER, INTENT(IN) :: left_end, right_end

    !     Local variables
    INTEGER             :: i, j, itemp
    REAL                :: temp
    
    DO i = left_end, right_end - 1
       DO j = i+1, right_end
          IF (list(i) > list(j)) THEN
             temp = list(i); list(i) = list(j); list(j) = temp
             IF(present(order)) THEN
                itemp = order(i); order(i) = order(j); order(j) = itemp
             ENDIF
          END IF
       END DO
    END DO

  END SUBROUTINE interchange_sort

END SUBROUTINE quick_sort

!-----------------------------------------
!
! produce a gaussian random distribution with stddev equal to sigma centered on 0
!
! INOUT: real array
! INPUT: stddev

SUBROUTINE random_gaussian(array, sigma)

  IMPLICIT NONE
  REAL(kind=4), DIMENSION(:), INTENT(INOUT) :: array
  REAL(kind=4), INTENT(IN) :: sigma
  !..local variables
  REAL(kind=8) :: v1, v2, R(2), rsq
  INTEGER :: i

  DO i=1,SIZE(array)
     rsq=1.d0
     DO WHILE(rsq>=1.d0)
        CALL RANDOM_NUMBER(R)
        v1=2.d0*R(1)-1.d0
        v2=2.d0*R(2)-1.d0                        
        rsq=v1*v1+v2*v2
     END DO
     array(i)=REAL(v1*dSQRT(-2.d0*dlog(rsq)/rsq),KIND=4)*sigma
  END DO

END SUBROUTINE random_gaussian

!---------------------------------------------
!
! image smoothing with gaussian filter
! INPUT: image, filter radius
! OUTPUT: smoothed image

FUNCTION ImSmooth(InputImage, FilterXYRad, UNDEF)

  IMPLICIT NONE
  REAL(kind=4), INTENT(IN)    :: InputImage(:,:), UNDEF
  INTEGER(kind=4), INTENT(IN) :: FilterXYRad
  REAL(kind=4) :: ImSmooth(SIZE(InputImage,DIM=1),SIZE(InputImage,DIM=2))
  REAL(kind=4), ALLOCATABLE :: w(:),w_norm(:), tmp(:,:), im_f(:,:) 
  INTEGER :: nmin, nmax, i, ix, iy, iz, DimX, DimY,  dd, xmin, x1, x2, y1, y2, z1, z2, j, k
  INTEGER, PARAMETER :: ns=3          !.. this parameter determines the filtering mask size (=2*ns*FilterRad)

  DimX=SIZE(InputImage,DIM=1); DimY=SIZE(InputImage,DIM=2)

  !..assign default for edges and undefined regions
  ImSmooth=UNDEF

  dd=ns*FilterXYRad
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

  !..allocate temporary images including padding
  ALLOCATE(tmp(-dd:DimX+dd,-dd:DimY+dd),im_f(-dd:DimX+dd,-dd:DimY+dd))
  im_f=UNDEF
  im_f(1:DimX,1:DimY)=InputImage(1:DimX,1:DimY)
  tmp=UNDEF

  !..smooth first in x
  DO iy=1,DimY
     DO ix=1,DimX
              
        IF(ALL(Im_f(ix-dd:ix+dd,iy)==UNDEF)) CYCLE !..no data to filter

        IF(ANY(Im_f(ix-dd:ix+dd,iy)==UNDEF)) THEN !..we can't use previously normalized weight here
           tmp(ix,iy)=SUM(Im_f(ix-dd:ix+dd,iy)*w(:),MASK=Im_f(ix-dd:ix+dd,iy)/=UNDEF)/&
                      SUM(w(:),MASK=Im_f(ix-dd:ix+dd,iy)/=UNDEF)
        ELSE
           tmp(ix,iy)=SUM(Im_f(ix-dd:ix+dd,iy)*w_norm(:))
        END IF

     END DO
  END DO
  

!..smooth now in y
  DO iy=1,DimY
     DO ix=1,DimX
                            
        IF(ALL(tmp(ix,iy-dd:iy+dd)==UNDEF)) THEN
           ImSmooth(ix,iy)=UNDEF
           CYCLE
        END IF

        IF(ANY(tmp(ix,iy-dd:iy+dd)==UNDEF)) THEN !..we can't use normalized weight here
           ImSmooth(ix,iy)=SUM(tmp(ix,iy-dd:iy+dd)*w(:),MASK=tmp(ix,iy-dd:iy+dd)/=UNDEF)/&
                SUM(w(:),MASK=tmp(ix,iy-dd:iy+dd)/=UNDEF)
        ELSE
           ImSmooth(ix,iy)=SUM(tmp(ix,iy-dd:iy+dd)*w_norm(:))
        END IF


     END DO
  END DO

END FUNCTION ImSmooth



END MODULE StatLib

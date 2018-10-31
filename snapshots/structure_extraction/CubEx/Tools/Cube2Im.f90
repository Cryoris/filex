PROGRAM main

  USE StatLib
  USE CubeLib
  IMPLICIT NONE
  CHARACTER(len=350) :: InpFile, IdCube, OutputImage, VarCube, ImType, Output1d, vzerotype, ThisOutputImage, skycubename, &
       imageunits, string, bkgrm_outfile, snrmap, rvarfile
  INTEGER(kind=4)    :: SelId, k, i, j, zstart, zend, nend, nint, max_pixel(3), boxsm, boxsm_z, is, NL, DimX, DimY, DimZ, nlpad, gsm, &
       maskshift(3), zshift, boxmin(3), boxmax(3), xmin, ymin, idpad, cc_zmin, cc_zmax, xmax, ymax, sky_zmin, sky_zmax, OrCubeDims(3), &
       obj_zstart, obj_zend, ierr, masklayer
  REAL(kind=4) :: pixscale, pixAng, lmin, lmax, thismax, MeanClip, MedianClip, FinalSigma, rfact
  REAL(kind=4), ALLOCATABLE :: CheckCube(:,:,:), vel_0(:,:), skycube(:,:,:), sky(:), NoiseLayer(:,:), ImageVar(:,:), SmoothIm(:,:), &
       NoiseLayerVar(:,:), ImageSNR(:,:)
  REAL(kind=8) :: CheckCubeWCS(SIZE(WCS))
  CHARACTER(len=3) :: proj, CrossProj
  REAL(kind=4) :: vzero, this_vel, this_sum, vzero_lambda, this_sum2, s2, skyfthr, skylev, negnum
  REAL(kind=4), PARAMETER :: c_km_s=299792.458
  LOGICAL, ALLOCATABLE :: skyfilter(:)
  LOGICAL :: sbscale, writeNaN, nscale, bkgrm
  
 
  negnum=-1.0

  !.. collect input parameters
  CALL ReadCommandLine

  !..get original cube dimension
  CALL GetCubeSize(InpFile,OrCubeDims)

  Verbosity=2
  
  !..read checkcube if requested
  IF(TRIM(IdCube)/="??") THEN

     CALL ReadCube(IdCube)
     !..save WCS values for later
     CheckCubeWCS=WCS

     IF(SIZE(Cube,DIM=1)/=OrCubeDims(1).or.SIZE(Cube,DIM=2)/=OrCubeDims(2)) &        
          STOP "Original Cube and IdCube must have the same SPATIAL dimensions!"
     
     !..perform mask shift if requested
     zshift=0
     IF(ANY(maskshift/=0)) THEN
        print *, " shifting mask..."
        DO i=1,2
           Cube=CSHIFT(Cube,SHIFT=-maskshift(i),DIM=i)
        END DO
        IF(maskshift(3)<=SIZE(Cube,DIM=3)) THEN !..do a shift in the mask
           Cube=CSHIFT(Cube,SHIFT=-maskshift(3),DIM=3)
        ELSE !..do a shift in the cube
           zshift=maskshift(3)
        END IF
        print *, "done"
     END IF

     !..find obj_zstart and obj_zend
     CALL FindBoundingBox  !..this returns the boxmin and boxmax of the selected object     
     obj_zstart=boxmin(3)
     obj_zend=boxmax(3)

   
     !..if padding is requested, change filename to a cube section in x and y
     !..and allocate a smaller checkcube
     IF(idpad/=9999) THEN
        
        xmin=MAX(1,boxmin(1)-idpad)
        xmax=MIN(SIZE(Cube,DIM=1),boxmax(1)+idpad)
        ymin=MAX(1,boxmin(2)-idpad)
        ymax=MIN(SIZE(Cube,DIM=2),boxmax(2)+idpad)
        !cc_zmin=MAX(1,boxmin(3)-idpad(3))
        !cc_zmax=MIN(SIZE(Cube,DIM=3),boxmax(3)+idpad(3))
        !lmin=(cc_zmin-WCS(3))*WCS(11)+WCS(6)
        !lmax=(cc_zmax-WCS(3))*WCS(11)+WCS(6)

        WRITE(string,'(4(a,i4.4),a)') "[",xmin,":",xmax,",",ymin,":",ymax,",*]"
        InpFile=TRIM(InpFile)//TRIM(string)
        print *, "updating InpFile name to: ",TRIM(InpFile)
        IF(TRIM(VarCube)/="??") THEN
           VarCube=TRIM(VarCube)//TRIM(string)
           print *, "updating VarCube name to: ",TRIM(VarCube)
        END IF
        
        ALLOCATE(CheckCube(xmax-xmin+1,ymax-ymin+1,SIZE(Cube,DIM=3)))
        CheckCube(:,:,:)=Cube(xmin:xmax,ymin:ymax,:)

     ELSE

        ALLOCATE(CheckCube(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2),SIZE(Cube,DIM=3)))
        CheckCube=Cube

     END IF

  ELSE

     CheckCubeWCS=0.

  END IF

  !..read variance cube if requested
  IF(TRIM(VarCube)/="??") THEN
     CALL ReadCube(VarCube)
     ALLOCATE(Var(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2),SIZE(Cube,DIM=3))) !..NB: Var is a global array
     Var=Cube
     print *, "minmax var=", MINVAL(Var,MASK=Var/=UNDEF), MAXVAL(Var)
     IF(TRIM(rvarfile)/="??") THEN
        print *, "rescaling var using: ", TRIM(rvarfile)
        OPEN(11,file=rvarfile,action="read")
        !..skip header
        ierr=1
        DO WHILE(ierr/=0)
           READ(11,*,iostat=ierr) i, rfact
        END DO
        BACKSPACE(11)
        ierr=0
        DO 
           READ(11,*,iostat=ierr) i, rfact
           IF(ierr/=0) EXIT
           WHERE(Var(:,:,i)/=UNDEF) Var(:,:,i)=Var(:,:,i)*rfact
        END DO
        CLOSE(11)
     END IF
  END IF

  !..read sky cube (a cube before sky-sub), if requested
  !..and produce skyfilter
  IF(TRIM(skycubename)/="??") THEN

     CALL ReadCube(skycubename)
     
     !..perform dimnesionality check
     IF(SIZE(Cube,DIM=3)/=OrCubeDims(3)) STOP "SkyCube z-dimension must be equal to the Original Cube one!"     
     sky_zmin=1
     sky_zmax=SIZE(Cube,DIM=3)     
     !print *, "sky_zmin=",sky_zmin,"sky_zmax=",sky_zmax
     ALLOCATE(sky(sky_zmax-sky_zmin+1))

     !..produce sky, layer by layer
     print *, "computing sky..."
     DO i=1,SIZE(sky)
        sky(i)=Median(PACK(Cube(:,:,sky_zmin+i-1),MASK=Cube(:,:,sky_zmin+i-1)/=UNDEF))
     END DO

     !..produce skyfilter
     skylev=Median(sky)
     skyfilter = (sky>=skyfthr*skylev) 
     print *, "done"

  END IF

  !..read datacube
  !..NB: Cube WCS are automatically adjusted in the CFTISIO section
  CALL ReadCube(InpFile)
  
  !..mask layer if requested
  IF(masklayer>0) THEN
     Cube(:,:,masklayer)=UNDEF
  END IF

  !..remove layer by layer avg sigma clip if requested
  IF(bkgrm) THEN
     IF(TRIM(bkgrm_outfile)/="??") OPEN(11,file=bkgrm_outfile,action="write")
     DO i=1,SIZE(Cube,DIM=3)
        IF(ANY(Cube(:,:,i)/=UNDEF)) THEN
           CALL SigmaClip(PACK(Cube(:,:,i),MASK=Cube(:,:,i)/=UNDEF), MeanClip, MedianClip, FinalSigma)
           WHERE(Cube(:,:,i)/=UNDEF) Cube(:,:,i)=Cube(:,:,i)-MeanClip
           IF(TRIM(bkgrm_outfile)/="??") WRITE(11,*) i, MeanClip
        END IF
     END DO
     IF(TRIM(bkgrm_outfile)/="??") CLOSE(11)
  END IF

  !..apply skyfilter if requested
  IF(TRIM(skycubename)/="??") THEN
     print *, "applying sky filter..."
     DO i=1,SIZE(Cube,DIM=3)
        IF(skyfilter(i)) Cube(:,:,i)=UNDEF
     END DO
     print *, "done"
  END IF

  !..if Cube and CheckCube z-dimensions are different and if WCS are present, get initial pixel value of CheckCube in Cube
  IF(TRIM(IdCube)/="??") THEN
     IF(SIZE(Cube,DIM=3)>SIZE(CheckCube,DIM=3)) THEN
        IF(ALL(CheckCubeWCS==0.d0).and.ALL(WCS==0.d0)) &
             STOP "CheckCube and Cube have different z-dimension and WCS were not found to correct for that!"
        zstart=INT((CheckCubeWCS(6)-WCS(6))/WCS(11)+1)+zshift
        zend=zstart+SIZE(CheckCube,DIM=3)-1
        IF(zstart>SIZE(Cube,DIM=3).or.zend>SIZE(Cube,DIM=3)) STOP "mask outside of cube! Check mask or maskshift"
     ELSE
        IF(SIZE(Cube,DIM=3)<SIZE(CheckCube,DIM=3)) print *, "WARNING: CheckCube is larger than datacube!"
        zstart=1
        zend=SIZE(Cube,DIM=3)
     END IF
  END IF

  IF(boxsm/=0) THEN !..perform spatial smoothing
     IF(boxsm_z==0) THEN
        print *, "performing boxcar spatial filtering with size (pixels)=", boxsm
        CALL BoxCar(bcsmooth=boxsm,bcsmooth_z=1)
        IF(ALLOCATED(Var)) CALL BoxCarVar(bcsmooth=boxsm,bcsmooth_z=1)
     ELSE
        print *, "performing boxcar spatial and redshift filtering with size (pixels)=", boxsm, boxsm_z
        CALL BoxCar(bcsmooth=boxsm,bcsmooth_z=boxsm_z)
        IF(ALLOCATED(Var)) CALL BoxCarVar(bcsmooth=boxsm,bcsmooth_z=boxsm_z)
     END IF
  ELSEIF(boxsm_z/=0) THEN
     print *, "performing boxcar redshift filtering with size (pixels)=", boxsm_z
     CALL BoxCar(bcsmooth=1,bcsmooth_z=boxsm_z)
     IF(ALLOCATED(Var)) CALL BoxCarVar(bcsmooth=1,bcsmooth_z=boxsm_z)
  END IF

  IF(SelId/=0) THEN !..cut out everything but this object from final image

     IF(TRIM(CrossProj)=="no") THEN

        !..if requested, save a layer to use where there are no object voxels 
        IF(NL/=0) THEN

           DimX=SIZE(Cube,DIM=1); DimY=SIZE(Cube,DIM=2); DimZ=SIZE(Cube,DIM=3)

           IF(DimZ==1) THEN
              print *, "NB: input Cube is an image, setting 1 as default noise layer"
              NL=1 ; nlpad=0
              !..readjust mask parameters as well
              obj_zstart=1 ; obj_zend=1
           END IF

           IF(NL==-1) THEN !..set default NL as central layer of selected object
              SELECT CASE(TRIM(proj))
              CASE('x')
                 NL=DimX/2+1
              CASE('y')
                 NL=DimY/2+1
              CASE('z')
                 NL=(MIN(obj_zend,DimZ)-MAX(obj_zstart,1)+1)/2+obj_zstart
              END SELECT
              print *, "NB: setting central object layer as default noise layer, nl=",NL
           END IF

           IF(TRIM(proj)=="z") THEN
              IF(zstart/=1) THEN
                 print *, "readjusting noise layer in z-direction"
                 NL=NL+zstart-1  !..readjust NL in case cube is larger than current selection
                 print *, "new nl=",NL
                 !..readjust obj_zstart, obj_zend
                 obj_zstart=obj_zstart+zstart-1
                 obj_zend=obj_zend+zstart-1
              END IF
           END IF

           SELECT CASE(TRIM(proj))
           CASE('x')
              ALLOCATE(NoiseLayer(DimY,DimZ))
              NoiseLayer(:,:)=SUM(Cube(MAX(1,NL-nlpad):MIN(NL+nlpad,DimX),:,:),DIM=1,MASK=Cube(MAX(1,NL-nlpad):MIN(NL+nlpad,DimX),:,:)/=UNDEF)
              IF(ALLOCATED(Var)) THEN
                 ALLOCATE(NoiseLayerVar(DimY,DimZ))
                 NoiseLayerVar(:,:)=SUM(Var(MAX(1,NL-nlpad):MIN(NL+nlpad,DimX),:,:),DIM=1,MASK=Var(MAX(1,NL-nlpad):MIN(NL+nlpad,DimX),:,:)/=UNDEF)
              END IF
           CASE('y')
              ALLOCATE(NoiseLayer(DimX,DimZ))
              NoiseLayer(:,:)=SUM(Cube(:,MAX(1,NL-nlpad):MIN(NL+nlpad,DimY),:),DIM=2,MASK=Cube(:,MAX(1,NL-nlpad):MIN(NL+nlpad,DimY),:)/=UNDEF)
              IF(ALLOCATED(Var)) THEN
                 ALLOCATE(NoiseLayerVar(DimX,DimZ))
                 NoiseLayerVar(:,:)=SUM(Var(:,MAX(1,NL-nlpad):MIN(NL+nlpad,DimY),:),DIM=2,MASK=Var(:,MAX(1,NL-nlpad):MIN(NL+nlpad,DimY),:)/=UNDEF)                
              END IF
           CASE('z')
              ALLOCATE(NoiseLayer(DimX,DimY))
              NoiseLayer(:,:)=SUM(Cube(:,:,MAX(1,obj_zstart,NL-nlpad):MIN(NL+nlpad,DimZ,obj_zend)),DIM=3, &
                   MASK=Cube(:,:,MAX(1,obj_zstart,NL-nlpad):MIN(NL+nlpad,DimZ,obj_zend))/=UNDEF)
              IF(ALLOCATED(VAR)) THEN
                 ALLOCATE(NoiseLayerVar(DimX,DimY))
                 NoiseLayerVar(:,:)=SUM(Var(:,:,MAX(1,obj_zstart,NL-nlpad):MIN(NL+nlpad,DimZ,obj_zend)),DIM=3, &
                      MASK=Var(:,:,MAX(1,obj_zstart,NL-nlpad):MIN(NL+nlpad,DimZ,obj_zend))/=UNDEF)                 
              END IF
              !print *, obj_zstart, NL-nlpad
              print *, "NoiseLayer zmin=",MAX(1,obj_zstart,NL-nlpad)," zmax=",MIN(NL+nlpad,DimZ,obj_zend)
           END SELECT

           IF(TRIM(imtype)=="Mean".or.TRIM(imtype)=="Median") THEN
              NoiseLayer=NoiseLayer/(2*nlpad+1)
              IF(ALLOCATED(NoiseLayerVar)) THEN
                 NoiseLayerVar=NoiseLayerVar/((2*nlpad+1)**2)
              END IF
           END IF


        END IF

        WHERE(CheckCube/=REAL(SelId)) Cube(:,:,zstart:zend)=UNDEF
        IF(zstart>1) Cube(:,:,1:zstart-1)=UNDEF
        IF(zend<SIZE(Cube,DIM=3)) Cube(:,:,zend+1:SIZE(Cube,DIM=3))=UNDEF
        IF(ALLOCATED(Var)) THEN
           WHERE(CheckCube/=REAL(SelId)) Var(:,:,zstart:zend)=UNDEF
           IF(zstart>1) Var(:,:,1:zstart-1)=UNDEF
           IF(zend<SIZE(Cube,DIM=3)) Var(:,:,zend+1:SIZE(Cube,DIM=3))=UNDEF
        END IF

        
        !..if requested, apply the layer saved before 
        IF(NL/=0) THEN
           SELECT CASE(TRIM(proj))
           CASE('x')
              IF(ALLOCATED(Var)) THEN
                 WHERE(Cube(NL,:,:)==UNDEF) Var(NL,:,:)=NoiseLayerVar(:,:)
              END IF
              WHERE(Cube(NL,:,:)==UNDEF) Cube(NL,:,:)=NoiseLayer(:,:)
           CASE('y')
              IF(ALLOCATED(Var)) THEN
                 WHERE(Cube(:,NL,:)==UNDEF) Var(:,NL,:)=NoiseLayerVar(:,:)
              END IF
              WHERE(Cube(:,NL,:)==UNDEF) Cube(:,NL,:)=NoiseLayer(:,:)
           CASE('z')
              IF(ALLOCATED(Var)) THEN
                 WHERE(Cube(:,:,NL)==UNDEF) Var(:,:,NL)=NoiseLayerVar(:,:)
              END IF
              WHERE(Cube(:,:,NL)==UNDEF) Cube(:,:,NL)=NoiseLayer(:,:)
           END SELECT
           DEALLOCATE(NoiseLayer)
           IF(ALLOCATED(NoiseLayerVar)) DEALLOCATE(NoiseLayerVar)
        END IF       
              
        
     ELSE

        SELECT CASE(TRIM(CrossProj))
        CASE('yz')
           DO k=1,SIZE(CheckCube,DIM=3)
              DO j=1,SIZE(CheckCube,DIM=2)
                 IF(ALL(CheckCube(:,j,k)/=REAL(SelId))) THEN
                    Cube(:,j,k+zstart-1)=UNDEF
                    IF(ALLOCATED(Var))  Var(:,j,k+zstart-1)=UNDEF
                 END IF
              END DO
           END DO
        CASE('xz')
           DO k=1,SIZE(CheckCube,DIM=3)
              DO i=1,SIZE(CheckCube,DIM=1)
                 IF(ALL(CheckCube(i,:,k)/=REAL(SelId))) THEN
                    Cube(i,:,k+zstart-1)=UNDEF
                    IF(ALLOCATED(Var)) Var(i,:,k+zstart-1)=UNDEF
                 END IF
              END DO
           END DO
        CASE('xy')
           DO j=1,SIZE(CheckCube,DIM=2)
              DO i=1,SIZE(CheckCube,DIM=1)
                 IF(ALL(CheckCube(i,j,:)/=REAL(SelId))) THEN
                    Cube(i,j,:)=UNDEF
                    IF(ALLOCATED(Var)) Var(i,j,:)=UNDEF
                 END IF
              END DO
           END DO     
        CASE default
           print *, "selected value for -mask2d :",TRIM(CrossProj), "  not recognized! Type Cube2Im.x for help."
           STOP
        END SELECT

     END IF

  END IF


  IF(TRIM(ImType)=="VMap".or.TRIM(ImType)=="vmap") THEN

     !----- produce a velocity map ---------

     IF(WCS(11)==0.) STOP "No WCS found in input file!"

     CALL AssignWCS(type="Image")
     
     !..first gets the vzero values in the right units

     IF(vzero==-1.0) THEN !..find vzero pixel
        max_pixel=MAXLOC(Cube, MASK=Cube/=UNDEF)
        print *, "max pixel location=", max_pixel
        vzero=max_pixel(3)-0.5
        vzero_lambda=vzero*WCS(11)+WCS(6)
     ELSE
        !..convert vzero to pixel units if necessary
        IF(TRIM(vzerotype)=="lambda") THEN
           vzero_lambda=vzero
           vzero=(vzero-WCS(6))/WCS(11)
        ELSEIF(TRIM(vzerotype)=="kms") THEN
           vzero_lambda=c_km_s
           vzero=(vzero-WCS(6))/WCS(11)
           print *, WCS(6), WCS(11), vzero
        ELSE
           vzero_lambda=vzero*WCS(11)+WCS(6)
        END IF
     END IF

     ALLOCATE(Image(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2)))
     ALLOCATE(vel_0(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2)))

     !..loop trhought the cube or selected objet and produce zeroth-moment velocity map
     DO j=1,SIZE(Cube,DIM=2)
        DO i=1,SIZE(Cube,DIM=1)
           this_vel=0.
           this_sum=0.
           DO k=1,SIZE(Cube,DIM=3)
              IF(Cube(i,j,k)/=UNDEF.and.Cube(i,j,k)>0.) THEN
                 this_vel=this_vel+(k-0.5-vzero)*Cube(i,j,k)
                 this_sum=this_sum+Cube(i,j,k)
              END IF
           END DO
           IF(this_sum>0.) THEN
              vel_0(i,j)=this_vel/this_sum
              Image(i,j)=this_vel/this_sum*(WCS(11)/vzero_lambda)*c_km_s
           ELSE
              vel_0(i,j)=UNDEF
              IF(writeNaN) THEN
                 Image(i,j)=sqrt(negnum)
              ELSE
                 Image(i,j)=UNDEF
              END IF
           END IF
        END DO
     END DO

!..write output image (append _0 to the file name)
     is=INDEX(TRIM(OutputImage),".fits")
     IF(is==0) is=LEN_TRIM(OutputImage)+1
     ThisOutputImage=TRIM(OutputImage(1:is-1))//"_0.fits"    
     CALL WriteImage(ThisOutputImage)
  
     !..calculate now the first moment (dispersion)
     DO j=1,SIZE(Cube,DIM=2)
        DO i=1,SIZE(Cube,DIM=1)
           this_vel=0.
           this_sum=0.
           this_sum2=0.
           IF(vel_0(i,j)==UNDEF) THEN
              this_vel=UNDEF
              this_sum=0
           ELSE
              DO k=1,SIZE(Cube,DIM=3)
                 IF(Cube(i,j,k)/=UNDEF.and.Cube(i,j,k)>0.) THEN
                    this_vel=this_vel+(k-0.5-vzero-vel_0(i,j))**2*Cube(i,j,k)
                    this_sum=this_sum+Cube(i,j,k)
                    this_sum2=this_sum2+Cube(i,j,k)**2
                 END IF
              END DO
           END IF
           IF(this_sum>0..and.this_sum*this_sum>this_sum2) THEN
              !..unbiased weighted estimator of the variance:
              s2=this_vel*this_sum/(this_sum*this_sum-this_sum2)
              Image(i,j)=sqrt(s2)*(WCS(11)/vzero_lambda)*c_km_s
           ELSE
              IF(writeNaN) THEN
                 Image(i,j)=sqrt(negnum)
              ELSE
                 Image(i,j)=UNDEF
              END IF
           END IF
        END DO
     END DO     

!..write output image (append _1 to the file name)
     is=INDEX(TRIM(OutputImage),".fits")
     IF(is==0) is=LEN_TRIM(OutputImage)+1
     ThisOutputImage=TRIM(OutputImage(1:is-1))//"_1.fits"    
     CALL WriteImage(ThisOutputImage)
  

  ELSE

     !.. produce image along requested projection excluding UNDEF pixels
     SELECT CASE(TRIM(proj))

     CASE('x')

        CALL AssignWCS(type="2Dspectrum")
                      
        ALLOCATE(Image(SIZE(Cube,DIM=3),SIZE(Cube,DIM=2)))
        Image=UNDEF
        IF(TRIM(ImType)=="Flux".or.TRIM(ImType)=="flux") THEN
           Image(:,:)=TRANSPOSE(SUM(Cube(:,:,:), DIM=1, MASK=(Cube(:,:,:)/=UNDEF)))
        ELSEIF(TRIM(ImType)=="Mean".or.TRIM(ImType)=="mean") THEN
           Image(:,:)=TRANSPOSE(SUM(Cube(:,:,:), DIM=1, MASK=(Cube(:,:,:)/=UNDEF))/MAX(1,COUNT(Cube(:,:,:)/=UNDEF,DIM=1)))
        ELSEIF(TRIM(ImType)=="Median".or.TRIM(ImType)=="median") THEN
           DO j=1,SIZE(Cube,DIM=2)
              DO i=1,SIZE(Cube,DIM=3)
                 IF(ANY(Cube(:,j,i)/=UNDEF)) THEN
                    Image(i,j)=Median(PACK(Cube(:,j,i),MASK=Cube(:,j,i)/=UNDEF))
                 ELSE
                    Image(i,j)=UNDEF
                 END IF
              END DO
           END DO
        END IF

        IF(ALLOCATED(Var)) THEN
           ALLOCATE(ImageVar(SIZE(Cube,DIM=3),SIZE(Cube,DIM=2)))
           ImageVar=UNDEF
           IF(TRIM(ImType)=="Flux".or.TRIM(ImType)=="flux") THEN
              ImageVar(:,:)=TRANSPOSE(SUM(Var(:,:,:), DIM=1, MASK=(Var(:,:,:)/=UNDEF)))
           ELSEIF(TRIM(ImType)=="Mean".or.TRIM(ImType)=="mean") THEN
              ImageVar(:,:)=TRANSPOSE(SUM(Var(:,:,:), DIM=1, MASK=(Var(:,:,:)/=UNDEF))/MAX(1,COUNT(Var(:,:,:)/=UNDEF,DIM=1)))
           ELSEIF(TRIM(ImType)=="Median".or.TRIM(ImType)=="median") THEN
              DO j=1,SIZE(Var,DIM=2)
                 DO i=1,SIZE(Var,DIM=3)
                    IF(ANY(Var(:,j,i)/=UNDEF)) THEN
                       ImageVar(i,j)=Median(PACK(Var(:,j,i),MASK=Cube(:,j,i)/=UNDEF))
                    ELSE
                       ImageVar(i,j)=UNDEF
                    END IF
                 END DO
              END DO
           END IF
        END IF

     CASE('y')

        CALL AssignWCS(type="2Dspectrum")
        
        ALLOCATE(Image(SIZE(Cube,DIM=3),SIZE(Cube,DIM=1)))
        Image=UNDEF
        IF(TRIM(ImType)=="Flux".or.TRIM(ImType)=="flux") THEN
           Image(:,:)=TRANSPOSE(SUM(Cube(:,:,:), DIM=2, MASK=(Cube(:,:,:)/=UNDEF)))
        ELSEIF(TRIM(ImType)=="Mean".or.TRIM(ImType)=="mean") THEN
           Image(:,:)=TRANSPOSE(SUM(Cube(:,:,:), DIM=2, MASK=(Cube(:,:,:)/=UNDEF))/MAX(1,COUNT(Cube(:,:,:)/=UNDEF,DIM=2)))
        ELSEIF(TRIM(ImType)=="Median".or.TRIM(ImType)=="median") THEN
           DO j=1,SIZE(Cube,DIM=1)
              DO i=1,SIZE(Cube,DIM=3)
                 IF(ANY(Cube(j,:,i)/=UNDEF)) THEN
                    Image(i,j)=Median(PACK(Cube(j,:,i),MASK=Cube(j,:,i)/=UNDEF))
                 ELSE
                    Image(i,j)=UNDEF
                 END IF
              END DO
           END DO
        END IF

        IF(ALLOCATED(Var)) THEN
           ALLOCATE(ImageVar(SIZE(Cube,DIM=3),SIZE(Cube,DIM=1)))
           ImageVar=UNDEF
           IF(TRIM(ImType)=="Flux".or.TRIM(ImType)=="flux") THEN
              ImageVar(:,:)=TRANSPOSE(SUM(Var(:,:,:), DIM=2, MASK=(Var(:,:,:)/=UNDEF)))
           ELSEIF(TRIM(ImType)=="Mean".or.TRIM(ImType)=="mean") THEN
              ImageVar(:,:)=TRANSPOSE(SUM(Var(:,:,:), DIM=2, MASK=(Var(:,:,:)/=UNDEF))/MAX(1,COUNT(Var(:,:,:)/=UNDEF,DIM=2)))
           ELSEIF(TRIM(ImType)=="Median".or.TRIM(ImType)=="median") THEN
              DO j=1,SIZE(Var,DIM=1)
                 DO i=1,SIZE(Var,DIM=3)
                    IF(ANY(Var(j,:,i)/=UNDEF)) THEN
                       ImageVar(i,j)=Median(PACK(Var(j,:,i),MASK=Var(j,:,i)/=UNDEF))
                    ELSE
                       ImageVar(i,j)=UNDEF
                    END IF
                 END DO
              END DO
           END IF
        END IF
        

     CASE('z')

        CALL AssignWCS(type="Image")
        
        ALLOCATE(Image(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2)))
        Image=UNDEF
        IF(TRIM(ImType)=="Flux".or.TRIM(ImType)=="flux") THEN
           Image(:,:)=SUM(Cube(:,:,:), DIM=3, MASK=(Cube(:,:,:)/=UNDEF))
        ELSEIF(TRIM(ImType)=="Mean".or.TRIM(ImType)=="mean") THEN
           Image(:,:)=SUM(Cube(:,:,:), DIM=3, MASK=(Cube(:,:,:)/=UNDEF))/MAX(1,COUNT(Cube(:,:,:)/=UNDEF,DIM=3))
        ELSEIF(TRIM(ImType)=="Median".or.TRIM(ImType)=="median") THEN
           DO j=1,SIZE(Cube,DIM=2)
              DO i=1,SIZE(Cube,DIM=1)
                 IF(ANY(Cube(i,j,:)/=UNDEF)) THEN
                    Image(i,j)=Median(PACK(Cube(i,j,:),MASK=Cube(i,j,:)/=UNDEF))
                 ELSE
                    Image(i,j)=UNDEF
                 END IF
              END DO
           END DO
        END IF

        IF(ALLOCATED(Var)) THEN
           ALLOCATE(ImageVar(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2)))
           ImageVar=UNDEF
           IF(TRIM(ImType)=="Flux".or.TRIM(ImType)=="flux") THEN
              ImageVar(:,:)=SUM(Var(:,:,:), DIM=3, MASK=(Var(:,:,:)/=UNDEF))
           ELSEIF(TRIM(ImType)=="Mean".or.TRIM(ImType)=="mean") THEN
              ImageVar(:,:)=SUM(Var(:,:,:), DIM=3, MASK=(Var(:,:,:)/=UNDEF))/MAX(1,COUNT(Var(:,:,:)/=UNDEF,DIM=3))
           ELSEIF(TRIM(ImType)=="Median".or.TRIM(ImType)=="median") THEN
              DO j=1,SIZE(Var,DIM=2)
                 DO i=1,SIZE(Var,DIM=1)
                    IF(ANY(Var(i,j,:)/=UNDEF)) THEN
                       ImageVar(i,j)=Median(PACK(Var(i,j,:),MASK=Var(i,j,:)/=UNDEF))
                    ELSE
                       ImageVar(i,j)=UNDEF
                    END IF
                 END DO
              END DO
           END IF
        END IF

     CASE default
     print *, "selected value for -proj :",TRIM(proj), "  not recognized! Type Cube2Im for help."
     STOP
  END SELECT

!..if requested, produce a SNRmap
  IF(TRIM(snrmap)/="??") THEN
     ALLOCATE(ImageSNR(SIZE(Image,DIM=1),SIZE(Image,DIM=2)))
     WHERE(Image/=UNDEF.and.ImageVar>0.) 
        ImageSNR=Image/sqrt(ImageVar)
     ELSEWHERE
        ImageSNR=UNDEF
     END WHERE
  END IF


  !..rescale image values if requested
  IF(sbscale) THEN
     IF(TRIM(proj)=="z") THEN
        IF(TRIM(ImType)=="Flux".or.TRIM(ImType)=="flux") THEN
           imageunits="10**(-18)*erg/s/cm**2/arcsec**2"
           print *, "rescaling units to: ", TRIM(imageunits)
           print *, "Total (Obj) Flux=",SUM(Image,MASK=Image/=UNDEF)*pixAng*1.e-20," erg/s/cm**2"
           IF(ALLOCATED(ImageVar)) THEN
              print *, "Total Obj std= ",sqrt(SUM(ImageVar,MASK=ImageVar/=UNDEF))*pixAng*1.e-20," erg/s/cm**2"
           END IF
           WHERE(Image/=UNDEF) Image=Image/(pixscale*pixscale)*pixAng
        ELSE
           imageunits="10**(-18)*erg/s/cm**2/arcsec**2/Angstrom"       
           print *, "rescaling units to: ", TRIM(imageunits)
           print *, "Total (Obj) Flux density=",SUM(Image,MASK=Image/=UNDEF)*1.e-20," erg/s/cm**2/Angstrom"
           IF(ALLOCATED(ImageVar)) THEN
              print *, "Total (Obj) Flux density std= ",sqrt(SUM(ImageVar,MASK=ImageVar/=UNDEF))*1.e-20," erg/s/cm**2/Angstrom"
           END IF
           WHERE(Image/=UNDEF) Image=Image/(pixscale*pixscale)
        END IF
     ELSE
        IF(TRIM(ImType)=="Flux".or.TRIM(ImType)=="flux") THEN
           imageunits="10**(-18)*erg/s/cm**2/arcsec/Angstrom"
           print *, "rescaling units to: ", TRIM(imageunits)
           WHERE(Image/=UNDEF) Image=Image/pixscale
        ELSE
           imageunits="10**(-18)*erg/s/cm**2/arcsec**2/Angstrom"       
           print *, "rescaling units to: ", TRIM(imageunits)
           WHERE(Image/=UNDEF) Image=Image/(pixscale*pixscale)
        END IF
     END IF       
     WHERE(Image/=UNDEF) Image=Image*0.01 !..convert from MUSE standard (10^-20) to (10^-18)
  END IF

  !..apply gaussian smoothing if requested
  IF(gsm>0) THEN
     print *, "gaussian smoothing final image..."
     Image=ImSmooth(Image,gsm,UNDEF)
     print *, "done"
  END IF

  !..if requested, normalize 2D spectra to the maximum (after smoothing)
  IF(nscale) THEN
     DO j=1,SIZE(Image,DIM=2)
        IF(ANY(Image(:,j)/=UNDEF)) THEN
           thismax=MAXVAL(Image(:,j))
           IF(thismax>0.) THEN
              Image(:,j)=Image(:,j)/thismax
           ELSE
              print *, "NB: maximum value in row=",j," is <=0! Setting row UNDEF..."
              Image(:,j)=UNDEF
           END IF
        END IF
     END DO
  END IF

  !..write output image
  IF(writeNaN) THEN
     WHERE(Image==UNDEF.or.Image==0.0) Image=sqrt(negnum)
  ELSE
     WHERE(Image==UNDEF) Image=0.
  END IF
  CALL WriteImage(OutputImage,imageunits=imageunits)

  !..write SNR map if requested
  IF(TRIM(snrmap)/="??") THEN
     Image=ImageSNR
     imageunits="SNR"
     IF(writeNaN) THEN
        WHERE(Image==UNDEF.or.Image==0.0) Image=sqrt(negnum)
     ELSE
        WHERE(Image==UNDEF) Image=0.
     END IF
     CALL WriteImage(snrmap,imageunits=imageunits)    
  END IF

END IF





CONTAINS

SUBROUTINE ReadCommandLine

  IMPLICIT NONE
  CHARACTER(len=300) :: ExeName, fname, string, ds, opt, arg 
  INTEGER :: narg, iarg, i, is
  LOGICAL :: ex, set_imtype


!..default
  InpFile="??"
  IdCube="??"
  OutputImage="??"
  VarCube="??"
  SelId=0
  proj='z'
  CrossProj='no'
  ImType="Mean"
  vzero=-1.0
  vzerotype="??"
  boxsm=0
  boxsm_z=0
  skyfthr=1.5
  skycubename="??"
  NL=0
  nlpad=0
  maskshift=0
  sbscale=.false.
  nscale=.false.
  pixscale=0.2
  pixAng=1.25
  writeNaN=.false.
  idpad=0
  gsm=0
  bkgrm=.false.
  bkgrm_outfile="??"
  snrmap="??"
  rvarfile="??"
  masklayer=-1

  set_imtype=.false.

  CALL GetArg(0,ExeName)
  narg=iargc()
  IF(narg<1) THEN
     print *, " "
     WRITE(*,'(2a)')"        Cube2Im (part of CubEx package)   "
     WRITE(*,'(a)') " "
     WRITE(*,'(a)') " by Sebastiano Cantalupo (cantalupo@phys.ethz.ch)"
     print *, " "
     WRITE(*,'(a)') "usage: Cube2Im.x -cube <name> -out <name> [-option <val>]"
     WRITE(*,'(a)') " "
     WRITE(*,'(a)') "  options:"
     WRITE(*,'(a)') "  -cube               <string>          : input datacube file name (NO DEFAULT)"
     WRITE(*,'(a)') "  -out                <string>          : output image file name (default=<cubename>.IM.fits)"
     WRITE(*,'(a)') "  -snrmap             <string>          : if requested, produce a snrmap with this name"
     WRITE(*,'(a)') "  -varcube            <string>          : associated Variance cube, if a SNR map is requested"
     WRITE(*,'(a)') "  -rvarfile           <string>          : if provided, rescale the variance cube layer by layer by the values in this file"
     WRITE(*,'(a)') "  -proj               'x'/'y'/'z'       : image projection direction (default=z)"
     WRITE(*,'(a)') "  -imtype             <string>          : options: flux (sum over not UNDEF pixels, this is the DEFAULT if -id is set)"
     WRITE(*,'(a)') "                                        :          mean (mean over not UNDEF pixels, this is the DEFAULT if -id is not set)"
     WRITE(*,'(a)') "                                        :          median (median over not UNDEF pixels)"
     WRITE(*,'(a)') "                                        :          vmap (velocity map with respect to vzero, proj=z by default)"
     WRITE(*,'(a)') "  -vzero              <real>            : zero velocity value for imtype=vmap. Value is in pixel units if <4000.0, lambda otherwise."
     WRITE(*,'(a)') "                                        :  If not defined, default is given by the pixel with the largest flux in the cube or selected object"
     WRITE(*,'(a)') "  -vzerotype          <string>          : force vzero to be 'pixel' or 'lambda' (Angstrom) "
     WRITE(*,'(a)') "  -idcube             <string>          : checkcube from CubEx with object Id. Default is cubename+'Objects_Id'.fits "
     WRITE(*,'(a)') "                                          or cubename+'Objects_Id_Assoc'.fits if the former does not exist in the current dir."
     WRITE(*,'(a)') "  -id                 <int>             : Id of selected object for output image (default=0, no selection). All the pixels"
     WRITE(*,'(a)') "                                          that are not associated with this object will be removed before projection."
     WRITE(*,'(a)') "  -idpad              <int>             : box padding around selected object in both SPATIAL directions (default=0). "
     WRITE(*,'(a)') "                                           NB: larger than cube dim is allowed. Use -1 for whole cube."
     WRITE(*,'(a)') "  -mask2d          'no'/'xy'/'xz'/'yz'  : if /='no', apply a 2d mask obtained from the selected plane rather than a 3d mask (default='no')"
     WRITE(*,'(a)') "  -boxsm              <int>             : if /=0, performs a spatial  boxcar smoothing with this size to the cube *before* any projection, (def=0)"
     WRITE(*,'(a)') "  -boxsm_z            <int>             : if /=0, performs a redshift boxcar smoothing with this size to the cube *before* any projection, (def=0)"
     WRITE(*,'(a)') "  -gsm                <int>             : if /=0, performs spatial gaussian smoothing with this 'radius' on final image *after* projection (def=0)"
     WRITE(*,'(a)') "                                          NB: a gsm of 1 means that the gaussian smoothing filter size is 3 pixels."  
     WRITE(*,'(a)') "  -skycube            <name>            : if defined, use this cube with unsubtracted sky with the option -skyfthr below to filter out"
     WRITE(*,'(a)') "                                          regions in the spectra with sky-lines before producing images."
     WRITE(*,'(a)') "  -skyfthr            <real>            : threshold for sky-line selection for the sky-line-filter (default=1.5 times the avg sky continuum)."
     WRITE(*,'(a)') "  -nl                 <int>             : for 3d mask, substitute pixels outside of mask with values taken from this layer if /=0 (default=0)."
     WRITE(*,'(a)') "                                           use -1 to automatically select the central layer of the object according to selected projection."
     WRITE(*,'(a)') "  -nlpad              <int>             : if /=0, add to the noise layer the adjacent nlpad layers in both directions (e.g., '-nlpad 1' -> 3 layers)"
     WRITE(*,'(a)') "                                          NB: using a very large number will give, e.g. the pseudo-NB cut around the object." 
     WRITE(*,'(a)') '  -maskshift          <int array>       : shift the checkcube mask in 3D by this value (default="0 0 0"), useful to check "noise" '
     WRITE(*,'(a)') '  -sbscale            <bol>             : if .true. rescales image values for MUSE cubes to 10**-18cgs/arcsec^2[/A] (default=.false.)'
     WRITE(*,'(a)') '  -nscale             <bol>             : if .true. rescales 2D spectral image values to the maximum in spectral direction (default=.false.)'
     WRITE(*,'(a)') '                                            NB: this is useful to understand the variation of FWHM with spatial position'
     WRITE(*,'(a)') "  -writeNaN           <bol>             : if .true., write NaNs instead of 0 for UNDEF pixels in the image (default=.false.)"
     WRITE(*,'(a)') "  -bkgrm              <bol>             : if .true., removes from each layer the layer average sigma clip (default=.false.)"
     WRITE(*,'(a)') "  -bkgrm_oufile       <name>            : if defined, save the layer by layer average sigma clip values in this file"
     STOP
  ELSEIF(narg==1) THEN !..use defaults
     CALL GetArg(1,InpFile)
     is=INDEX(TRIM(InpFile),".fits")
     IF(is==0) is=LEN_TRIM(InpFile)+1
     OutputImage=TRIM(InpFile(1:is-1))//".IM.fits"
     RETURN
  ELSE
     !..read options from command line
     DO i=1,narg,2
        CALL getarg(i,opt)
        CALL getarg(i+1,arg)
        SELECT CASE(TRIM(opt))
        CASE('-cube')          ; READ(arg,'(a)') InpFile
        CASE('-out')           ; READ(arg,'(a)') OutputImage
        CASE('-snrmap')        ; READ(arg,'(a)') snrmap
        CASE('-idcube')        ; READ(arg,'(a)') IdCube
        CASE('-varcube')       ; READ(arg,'(a)') VarCube
        CASE('-imtype')        ; READ(arg,'(a)') ImType ; set_imtype=.true.
        CASE('-vzerotype')     ; READ(arg,'(a)') vzerotype
        CASE('-id')            ; READ(arg,*) SelId
        CASE('-proj')          ; READ(arg,*) proj
        CASE('-mask2d')        ; READ(arg,*) CrossProj
        CASE('-vzero')         ; READ(arg,*) vzero
        CASE('-boxsm')         ; READ(arg,*) boxsm
        CASE('-boxsm_z')       ; READ(arg,*) boxsm_z
        CASE('-skycube')       ; READ(arg,'(a)') skycubename
        CASE('-skyfthr')       ; READ(arg,*) skyfthr
        CASE('-nl')            ; READ(arg,*) NL
        CASE('-nlpad')         ; READ(arg,*) nlpad
        CASE('-maskshift')     ; READ(arg,*) maskshift(:)
        CASE('-sbscale')       ; READ(arg,*) sbscale
        CASE('-nscale')        ; READ(arg,*) nscale
        CASE('-writeNaN')      ; READ(arg,*) writeNaN   
        CASE('-idpad')         ; READ(arg,*) idpad
        CASE('-gsm')           ; READ(arg,*) gsm
        CASE('-bkgrm')         ; READ(arg,*) bkgrm
        CASE('-bkgrm_outfile') ; READ(arg,'(a)') bkgrm_outfile
        CASE('-rvarfile')      ; READ(arg,'(a)') rvarfile
        CASE('-masklayer')     ; READ(arg,*) masklayer
        CASE default
           print *, "command line argument ",TRIM(opt), " not recognized!"
           STOP
        END SELECT
     END DO
  END IF

!..perform few checks
  IF(TRIM(InpFile)=="??") THEN
     print *, "please provide the input datacube with the -InpFile option!"
     STOP
  END IF
  IF(TRIM(OutputImage)=="??") THEN
     is=INDEX(TRIM(InpFile),".fits")
     IF(is==0) is=LEN_TRIM(InpFile)+1
     OutputImage=TRIM(InpFile(1:is-1))//".IM.fits" 
  END IF
  IF(CrossProj==proj) STOP "CrossProj and Proj axis cannot be the same! Please adjust the -Proj and/or -CrossProj options."
!  IF(TRIM(ImType)/="Flux".and.TRIM(ImType)/="Mean".and.TRIM(ImType)/="Median") STOP "Requested ImType is not a valid option (use: Flux, Mean or Median)"

  IF(SelId/=0.and.TRIM(IdCube)=="??") THEN !..associate default name
     is=INDEX(TRIM(InpFile),".fits")
     IF(is==0) is=LEN_TRIM(InpFile)+1
     IdCube=TRIM(InpFile(1:is-1))//".Objects_Id.fits"
     INQUIRE(FILE=IdCube,EXIST=ex)
     IF(.not.ex) IdCube=TRIM(InpFile(1:is-1))//".Objects_Id_Assoc.fits"
     INQUIRE(FILE=IdCube,EXIST=ex)
     IF(.not.ex) STOP "Please provide the name of the associated checkcube with the objects Id with the option -idcube!"
  END IF

  IF(TRIM(ImType)=="VMap".or.TRIM(ImType)=="vmap") THEN
     !..check if vzerotype has a reasonable value
     IF(TRIM(vzerotype)/="??".and.TRIM(vzerotype)/="pixel".and.TRIM(vzerotype)/="lambda".and.TRIM(vzerotype)/="kms") THEN
        print *, "Selected value for vzerotype: ", TRIM(vzerotype), " not recognized!"
        STOP
     END IF
     !..force z-projection
     proj="z"
     !..check vzero units, if defined
     IF(vzero/=-1.0) THEN
        IF(TRIM(vzerotype)=="??") THEN
           IF(vzero<4000.0) THEN
              vzerotype="pixel"
           ELSE
              vzerotype="lambda"
           END IF
        END IF
     ELSE
        !..use pixels by default
        vzerotype="pixel"
     END IF

     IF(TRIM(OutputImage)=="??") THEN
        is=INDEX(TRIM(InpFile),".fits")
        IF(is==0) is=LEN_TRIM(InpFile)+1
        OutputImage=TRIM(InpFile(1:is-1))//".VMAP.fits" 
     END IF
   
  END IF

  !..check that nscale is not requested with proj=z
  IF(TRIM(proj)=="z".and.nscale) STOP "-nscale .true. is not allowed for images (-proj z), only for 2D spectra (-proj x or y)!"

  !..adjust idpad to a very large number if option -1 was selected
  IF(idpad==-1) idpad=9999
  

!..select imtype default if necessary
  IF(.not.set_imtype) THEN
     IF(SelId/=0) THEN
        imtype="flux"
        print *, "NB: Setting ImType='flux' "
     ELSE
        imtype="mean"
        print *, "NB: Setting ImType='mean' "
     END IF
  END IF

!..check that a var cube is provided if a snrmap is requested
  IF(TRIM(snrmap)/="??") THEN
     IF(TRIM(VarCube)=="??") STOP "Please provide a variance cube with the option -varcube!"
  END IF


END SUBROUTINE ReadCommandLine


SUBROUTINE showplot(fname)

  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN) :: fname


  !..write a python script:

  OPEN(99,file="plot.py",action="write")
  WRITE(99,'(a)') "import matplotlib.pyplot as plt"
  WRITE(99,'(a)') "import numpy as np"
  WRITE(99,'(3a)') "f2 = open('",TRIM(fname),"', 'r')"
  WRITE(99,'(3a)') "lines = f2.readlines()"
  WRITE(99,'(3a)') "x1 = []"
  WRITE(99,'(3a)') "y1 = []"  
  WRITE(99,'(a)') "for line in lines:"
  WRITE(99,'(a)') "\t p = line.split()"
  WRITE(99,'(a)') "\t x1.append(float(p[1]))"
  WRITE(99,'(a)') "\t y1.append(float(p[2]))"
  WRITE(99,'(a)') " "  
  WRITE(99,'(a)') "xv = np.array(x1)"
  WRITE(99,'(a)') "yv = np.array(y1)"
  WRITE(99,'(a)') "plt.plot(xv, yv)"
  WRITE(99,'(a)') "plt.show()"
  CLOSE(99)

  CALL SYSTEM("python plot.py")

END SUBROUTINE showplot

SUBROUTINE AssignWCS(type)

  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN) :: type

  SELECT CASE(TRIM(type))
  CASE("2Dspectrum")
     !..assign WCS for 2D spectrum
     ALLOCATE(imWCS(6),imWCSstrings(4), imWCSlabels(6), imWCSlabels_strings(4))
     imWCSlabels(:)=["CRPIX1","CRPIX2","CRVAL1","CRVAL2","CD1_1 ","CD2_2 "]
     imWCSlabels_strings(:)=["CTYPE1 ","CTYPE2 ","CUNIT1 ","CUNIT2 "]
     imWCSstrings(:)=["AWAV    ","PIXEL   ","Angstrom","Pixels  "]
     imWCS(1)=WCS(3); imWCS(2)=1; imWCS(3)=WCS(6); imWCS(4)=1; imWCS(5)=WCS(11); imWCS(6)=1
  CASE("Image")
     !..assign WCS for image
     ALLOCATE(imWCS(8),imWCSstrings(4), imWCSlabels(8), imWCSlabels_strings(4))
     imWCSlabels(:)=["CRPIX1","CRPIX2","CRVAL1","CRVAL2","CD1_1 ","CD1_2 ","CD2_1 ","CD2_2 "]
     imWCSlabels_strings(:)=["CTYPE1 ","CTYPE2 ","CUNIT1 ","CUNIT2 "]
     imWCS(1)=WCS(1); imWCS(2)=WCS(2); imWCS(3)=WCS(4); imWCS(4)=WCS(5); imWCS(5:8)=WCS(7:10)
     imWCSstrings(1)=WCSstrings(1); imWCSstrings(2)=WCSstrings(2); imWCSstrings(3)=WCSstrings(4); imWCSstrings(4)=WCSstrings(5)
  CASE default
     print *, "type= ",TRIM(type)," in AssignWCS not recognized!"
  END SELECT

END SUBROUTINE AssignWCS

SUBROUTINE FindBoundingBox

  IMPLICIT NONE
  INTEGER(kind=4) :: i,j,k

  boxmin=100000
  boxmax=-100000
  
  !print *, "SelId=", SelId
  !print *, MINVAL(Cube), MAXVAL(Cube)

  DO k=1,SIZE(Cube,DIM=3)
     DO j=1,SIZE(Cube,DIM=2)
        DO i=1,SIZE(Cube,DIM=1)

           IF(Cube(i,j,k)==SelId) THEN
              boxmin(1)=MIN(boxmin(1),i)
              boxmin(2)=MIN(boxmin(2),j)
              boxmin(3)=MIN(boxmin(3),k)
              boxmax(1)=MAX(boxmax(1),i)
              boxmax(2)=MAX(boxmax(2),j)
              boxmax(3)=MAX(boxmax(3),k)
           END IF

        END DO
     END DO
  END DO

  IF(ALL(boxmin==100000).or.ALL(boxmax==100000)) STOP "Object with request id not found!"

END SUBROUTINE FindBoundingBox


END PROGRAM main

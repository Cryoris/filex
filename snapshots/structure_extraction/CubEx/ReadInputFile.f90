RECURSIVE SUBROUTINE Readinputfile

  USE Globalmodule
  USE StatLib
  IMPLICIT NONE
  INTEGER(kind=4) :: status, unit, readwrite, blocksize, naxes(3), nfound, group, in, end, ext, &
       cnt, dist2, r, ix, iy, ierr, x, y, is, i, n_ext_to_read, rank, ival, z, j, naxes4(4), fmin, fmax, &
       zval, x1, x2, y1, y2
  CHARACTER(len=350) :: fname, comment
  REAL(kind=4), ALLOCATABLE :: DATA(:,:,:), buffer(:), var_rescale(:), var_rescale_f(:)
  LOGICAL :: anyf, ex
  REAL(Kind=4) :: MeanClip, MedianClip, this_sigma, var_meanclip, var_medianclip, var_this_sigma, rescalef
  REAL(kind=8) :: dum

  IF(.not.MultiExt.and.TRIM(VarFile)=="??") THEN
     n_ext_to_read=1
  ELSE
     n_ext_to_read=2
  END IF

  DO ext=1,n_ext_to_read

     status=0

     !..get an unused unit
     CALL ftgiou(unit,status)

     IF(status/=0) STOP "problem with ftgiou"

     IF(ext==1) THEN
        IF(MultiExt) THEN
           !..open first extension: DATA
           fname=TRIM(InpFile)//"[1]"
        ELSE
           fname=TRIM(InpFile)
        END IF
     ELSE
        IF(MultiExt) THEN
           !..open second extension: STAT
           fname=TRIM(InpFile)//"[2]"
        ELSE
           fname=TRIM(VarFile)
        END IF
     END IF
     readwrite=0
     CALL ftdopn(unit,fname,readwrite,status)

     !..perform few checks
     IF(status/=0) THEN
        !..try to open a single fits file with no extension
        fname=TRIM(InpFile)
        status=0
        CALL ftgiou(unit,status)
        CALL ftdopn(unit,fname,readwrite,status)
        IF(status/=0) THEN
           print *, "Problem reading file:",TRIM(fname)
           STOP
        ELSEIF(ext==2) THEN !..maybe VarFile and/or MultiExt are not set in the options
           IF(TRIM(VarFile)=="??") THEN
              !..check if there is a variance file with a default name:
              is=INDEX(TRIM(InpFile),".fits")
              VarFile=TRIM(InpFile(1:is-1))//".VAR.fits"
              INQUIRE(File=VarFile, Exist=ex)
              IF(.not.ex) STOP "Please provide the variance cube with the option VarFile!"
           END IF
           MultiExt=.false.
        END IF
     END IF


     IF(Verbosity>=2) print *, "Reading:",TRIM(fname)

     !..read and store WCS information, if present
     IF(ext==1) THEN !..only need to do this once

        !..check that keywords are there...
        status=0
        CALL ftgkyd(unit,"CRPIX1",dum,comment,status)
        
        IF(status==0) THEN !read keywords

           DO i=1,SIZE(WCS)
              CALL ftgkyd(unit,TRIM(WCSlabels(i)),WCS(i),comment,status)
              IF(status/=0) THEN
                 print *, "WARNING:: WCS keyword ", TRIM(WCSlabels(i)), " not found!"
                 IF(i==11) THEN
                    !..sometimes it is because CD3_3 is called in a different way... in this case use
                    !..the default value for MUSE:
                    WCS(11)=1.25
                    print *, "          Using default MUSE value ", WCS(11)                   
                 END IF
                 status=0
              END IF
           END DO

           DO i=1,SIZE(WCSstrings)
              CALL ftgkys(unit,TRIM(WCSlabels_strings(i)),WCSstrings(i),comment,status)
              IF(status/=0) THEN
                 print *, "WARNING:: WCS keyword ", TRIM(WCSlabels_strings(i)), " not found!"
                 print *, "          Using default value ", TRIM(WCSlabels_strings(i)), " = ", TRIM(WCSstrings(i))
                 status=0 
              END IF
           END DO

           !..adjust zmin and zmax in CUNIT3 if requested
           IF(lmin/=0) THEN
              zmin=(lmin-WCS(6))/WCS(11)+WCS(3)
           END IF
           IF(lmax/=0) THEN
              zmax=(lmax-WCS(6))/WCS(11)+WCS(3)
           END IF

           !..update WCS(6) info:
           WCS(6)=WCS(6)+(zmin-WCS(3))*WCS(11)

        ELSE
           IF(lmin/=0.or.lmax/=0) STOP "Can't use lmin and lmax command without WCS info!"
           print *, " "
           print *, "NB: no WCS info found in input file"
           print *, "--> no ra & dec will be produced in the output catalogue"
           print *, " "
           WCS(:)=0.d0
        END IF

     END IF

     !..get cube rank
     status=0
     CALL ftgkyj(unit,'NAXIS',rank,comment,status)
     IF(status/=0) STOP "problem reading NAXIS keyword"

     status=0
     IF(rank==4) THEN
        CALL ftgknj(unit,'NAXIS',1,4,naxes4,nfound,status)
        naxes(1:3)=naxes4(1:3)
     ELSEIF(rank==3) THEN
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

     !..get cube size
     !CALL ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)

     IF(ext==1) THEN !..we only need to do this once:
        zmax=MIN(zmax,naxes(3))
        IF(Verbosity>=2) THEN
           print *, "Original Cube dimensions =", naxes
           print *, "Selected subset, zmin=", zmin,"zmax=",zmax,"pixels"
        END IF
        !..dimension without ghost zones
        DimX=naxes(1)
        DimY=naxes(2)
        DimZ=zmax-zmin+1
        ALLOCATE(Cube(0:DimX+1,0:DimY+1,0:DimZ+1),Var(0:DimX+1,0:DimY+1,0:DimZ+1),&
             Mask(0:DimX+1,0:DimY+1,0:DimZ+1))
        IF(.not.CFITSIOLibFix) THEN !..see note below
           ALLOCATE(DATA(DimX,DimY,DimZ))
        ELSE
           ALLOCATE(DATA(DimX,DimY,1:zmax))
        END IF
        Cube=0.
        Var=0.
        Mask=0.
     END IF

     group=1

     !..select subset to read
     IF(CFITSIOLibFix) THEN
        in=1   ! some version of CFITSIO library have problems starting from in/=1, so we have to read the array to the end and cut it
     ELSE
        in=(zmin-1)*naxes(1)*naxes(2)+1
     END IF
     end=zmax*naxes(1)*naxes(2)

     !..read cube subset
     CALL ftgpve(unit,group,in,end,UNDEF,DATA,anyf,status)
     IF(status/=0) THEN
        IF(.not.CFITSIOLibFix) THEN
           !..let's try to fix the lib and start again...
           DEALLOCATE(Cube,Var,Mask,DATA)
           CFITSIOLibFix=.true.
           print *, "Problem reading file!"
           print *, "NB: automatically setting CFITSIOLibFix=.true. and reading the file again..."
           CALL ReadInputFile
           EXIT
        ELSE
           STOP "Problem reading datacube (and CFISTIOLibFix doesn't solve the problem...)"
        END IF
     END IF

     !..copy within ghost zones
     IF(ext==1) THEN
        Cube(1:DimX,1:DimY,1:DimZ)=DATA(1:DimX,1:DimY,zmin:zmax)
        IF(Verbosity>=2) THEN
           print *, "Min and Max=",MINVAL(Cube,MASK=(Cube/=UNDEF)), MAXVAL(Cube)
        END IF
     ELSE
        Var(1:DimX,1:DimY,1:DimZ)=DATA(1:DimX,1:DimY,zmin:zmax)
        IF(Verbosity>=2) print *, "Min and Max=",MINVAL(Var,MASK=(Var/=UNDEF)), MAXVAL(Var)
     END IF

     IF(Verbosity>=2.and.anyf) print *, "NB: undefined pixels are present"

     CALL ftclos(unit,status)

     !..free all allocated units
     call ftfiou(-1,status)
     IF(status/=0) STOP "problem with ftfiou"


  END DO

!..apply data thresholding if requested
  IF(ApplyDataThreshold) THEN
     WHERE(Cube/=UNDEF.and.(Cube<Threshold(1).or.Cube>Threshold(2))) Cube=UNDEF
  END IF

  
!..if variance was not read, estimate it from the Cube layer by layer with a simple
!..sigma clipping 
  IF(n_ext_to_read==1) THEN
     print *, "Estimating variance layer by layer from datacube..."
     IF(TRIM(EstimatedVarOutFile)/="??")  THEN
        OPEN(22,file=EstimatedVarOutFile)
        WRITE(22,*) "         z     meanclip           sigma             sigma^2"
     END IF
     DO z=1,DimZ
        IF(ANY(Cube(1:DimX,1:DimY,z)/=UNDEF)) THEN
           CALL SigmaClip(PACK(Cube(1:DimX,1:DimY,z),MASK=Cube(1:DimX,1:DimY,z)/=UNDEF),meanclip,medianclip,this_sigma)
           Var(:,:,z)=this_sigma*this_sigma
           IF(TRIM(EstimatedVarOutFile)/="??") WRITE(22,*) z, meanclip, this_sigma, this_sigma*this_sigma
        ELSE
           Var(:,:,z)=UNDEF
        END IF
        WHERE(Cube(:,:,z)==UNDEF) Var(:,:,z)=UNDEF
     END DO
     IF(TRIM(EstimatedVarOutFile)/="??") CLOSE(22)
  END IF

!..if requested, rescale original variance by the one estimated layer by layer
!  IF(RescaleVar) THEN
!     print *, "Rescaling Variance layer by layer..."
!    IF(TRIM(RescalingVarOutFile)/="??")  THEN
!       OPEN(22,file=RescalingVarOutFile)
!       WRITE(22,*) "         z     rescaling_fact"
!    END IF
!     DO z=1,DimZ
!        IF(ANY(Cube(1:DimX,1:DimY,z)/=UNDEF)) THEN
!           CALL SigmaClip(PACK(Cube(1:DimX,1:DimY,z),MASK=Cube(1:DimX,1:DimY,z)/=UNDEF),meanclip,medianclip,this_sigma) 
!           CALL SigmaClip(PACK(Var(1:DimX,1:DimY,z),MASK=(Cube(1:DimX,1:DimY,z)/=UNDEF.and.Var(1:DimX,1:DimY,z)/=UNDEF)),&
!                var_meanclip,var_medianclip,var_this_sigma) 
!           IF(var_medianclip>0.) THEN
!              var_rescale=(this_sigma*this_sigma)/var_medianclip
!              var_rescale=MAX(MIN(var_rescale,RescaleVarMax),RescaleVarMin)
!              IF(TRIM(RescalingVarOutFile)/="??")  WRITE(22,*) z, var_rescale
!              WHERE(Var(:,:,z)/=UNDEF) Var(:,:,z)=Var(:,:,z)*var_rescale
!           ENDIF
!        END IF
!     END DO
!     IF(TRIM(RescalingVarOutFile)/="??") CLOSE(22)
!  END IF

!..if requested, rescale original variance by the one estimated layer by layer
  IF(RescaleVar) THEN
     IF(TRIM(RescalingVarInpFile)/="??") THEN
        print *, "Using Rescaling Variance factors from: ",TRIM(RescalingVarInpFile) 
        OPEN(22,file=RescalingVarInpFile,action="read")
        !..skip header
        ierr=1
        DO
           READ(22,*,iostat=ierr) zval, rescalef
           IF(ierr==0) THEN
              BACKSPACE(22)
              EXIT
           END IF
        END DO
        !..read file and rescale, layer by layer
        DO z=1,DimZ
           READ(22,*) zval, rescalef
           IF(zval/=z) STOP "RescalingVarInpFile does not match Variance cube!"
           WHERE(Var(:,:,z)/=UNDEF) Var(:,:,z)=Var(:,:,z)*rescalef
        END DO
        CLOSE(22)
     ELSE          
        print *, "Calculating Rescaling Variance factors layer by layer..."
        IF(ALL(RescaleVarArea/=-1)) THEN
           print *, "Using values within this area (xmin,xmax,ymin,ymax)=", RescaleVarArea
           x1=RescaleVarArea(1); x2=RescaleVarArea(2); y1=RescaleVarArea(3); y2=RescaleVarArea(4) 
        ELSE
           x1=1; x2=DimX; y1=1 ; y2=DimY
        END IF           
        ALLOCATE(var_rescale(DimZ))
        !..get rescaling factors
        DO z=1,DimZ
           IF(ANY(Cube(x1:x2,y1:y2,z)/=UNDEF)) THEN
              CALL SigmaClip(PACK(Cube(x1:x2,y1:y2,z),MASK=Cube(x1:x2,y1:y2,z)/=UNDEF),meanclip,medianclip,this_sigma) 
              CALL SigmaClip(PACK(Var(x1:x2,y1:y2,z),MASK=(Cube(x1:x2,y1:y2,z)/=UNDEF.and.Var(x1:x2,y1:y2,z)/=UNDEF)),&
                   var_meanclip,var_medianclip,var_this_sigma) 
              IF(var_medianclip>0.) THEN
                 var_rescale(z)=(this_sigma*this_sigma)/var_medianclip    
              ELSE
                 var_rescale(z)=UNDEF
              END IF
           ELSE
              var_rescale(z)=UNDEF
           END IF
        END DO
     !..if requested perform a median filtering of the rescaling factors
        IF(RescaleVarFR>0) THEN
           print *, "Applying requested median filtering of the rescaling factor..."
           ALLOCATE(var_rescale_f(DimZ))
           DO z=1,DimZ
              fmin=MAX(1,z-RescaleVarFR) ; fmax=MIN(DimZ,z+RescaleVarFR)
              IF(ANY(var_rescale(fmin:fmax)/=UNDEF.and.var_rescale(fmin:fmax)>RescaleVarMin&
                   .and.var_rescale(fmin:fmax)<RescaleVarMax)) THEN
                 var_rescale_f(z)=Median(PACK(var_rescale(fmin:fmax),MASK=var_rescale(fmin:fmax)/=UNDEF.and.&
                      var_rescale(fmin:fmax)>RescaleVarMin.and.var_rescale(fmin:fmax)<RescaleVarMax))
              ELSE
                 var_rescale_f(z)=UNDEF
              END IF
           END DO
           !..replace arrays and deallocate temporary one
           var_rescale=var_rescale_f
           DEALLOCATE(var_rescale_f)
        END IF
        !..replace values out of boundary
        WHERE(var_rescale==UNDEF) var_rescale=1.0
        WHERE(var_rescale<RescaleVarMin) var_rescale=RescaleVarMin
        WHERE(var_rescale>RescaleVarMax) var_rescale=RescaleVarMax
        !..............................................................
        !..finally, rescale variance cube layer by layer and write output file if requested
        IF(TRIM(RescalingVarOutFile)/="??")  THEN
           print *, "Saving Rescaling Variance factors here: ", TRIM(RescalingVarOutFile)
           OPEN(22,file=RescalingVarOutFile)
           WRITE(22,*) "#         z     rescaling_fact"
        END IF
        DO z=1,DimZ
           IF(TRIM(RescalingVarOutFile)/="??")  WRITE(22,*) z, var_rescale(z)
           WHERE(Var(:,:,z)/=UNDEF) Var(:,:,z)=Var(:,:,z)*var_rescale(z)
        END DO
        IF(TRIM(RescalingVarOutFile)/="??") CLOSE(22)
        DEALLOCATE(var_rescale)
     END IF
  END IF


!---------- MASKING ----------------------------

  IF(TRIM(ObjMaskList)/="??") THEN

     IF(Verbosity>=2) print *, "Applying an object mask list from this file:", TRIM(ObjMaskList),"..."
     
     cnt=0
     OPEN(1,file=ObjMaskList,action="read")
     DO 
        READ(1,*,iostat=ierr) x, y, r
        IF(ierr/=0) EXIT
        cnt=cnt+1
        DO iy=y-r,y+r
           DO ix=x-r,x+r
              dist2=(ix-x)**2+(iy-y)**2
              IF(dist2<=r*r) THEN
                 Cube(ix,iy,:)=UNDEF
                 Var(ix,iy,:)=UNDEF
              END IF
           END DO
        END DO
     END DO
     CLOSE(1)
     
     IF(Verbosity>=2) print *, "done. Found n=",cnt," objects."
     
  END IF

  IF(TRIM(LayerMaskList)/="??") THEN
     
     IF(Verbosity>=2) print *, "Applying a layer mask list from this file:", TRIM(LayerMaskList),"..."
     cnt=0
     OPEN(1,file=LayerMaskList,action="read")
     DO 
        READ(1,*,iostat=ierr) z
        IF(ierr/=0) EXIT
        cnt=cnt+1
        Cube(:,:,z)=UNDEF
        Var(:,:,z)=UNDEF
     END DO
     CLOSE(1)

     IF(Verbosity>=2) print *, "done. Found n=",cnt," layers to mask."

  END IF

  IF(ALLOCATED(SourceMask)) THEN
     
     IF(Verbosity>=2) print *, "Applying SourceMask..."
     IF(SIZE(SourceMask,DIM=1)/=SIZE(Cube,DIM=1)-2.or.&
         SIZE(SourceMask,DIM=2)/=SIZE(Cube,DIM=2)-2) STOP "SourceMask and Cube have different spatial dimensions!" !..NB: cubes have ghost zones
     DO j=1,SIZE(SourceMask,DIM=2)
        DO i=1,SIZE(SourceMask,DIM=1)
           IF(SourceMask(i,j)/=0) THEN
              Cube(i,j,:)=UNDEF
              Var(i,j,:)=UNDEF
           END IF
        END DO
     END DO
     
    IF(Verbosity>=2) print *, "done. Number of masked spaxels=",COUNT(SourceMask/=0)

 END IF

!..if requested, trim the edges of the cube by XYedge
  IF(XYedge/=0) THEN
     Cube(1:XYedge,:,:)=UNDEF
     Cube(DimX-XYedge:DimX,:,:)=UNDEF
     Cube(:,1:XYedge,:)=UNDEF
     Cube(:,DimY-XYedge:DimY,:)=UNDEF
     Cube(1:XYedge,1:XYedge,:)=UNDEF
     Cube(DimX-XYedge:DimX,DimY-XYedge:DimY,:)=UNDEF
     !....
     Var(1:XYedge,:,:)=UNDEF
     Var(DimX-XYedge:DimX,:,:)=UNDEF
     Var(:,1:XYedge,:)=UNDEF
     Var(:,DimY-XYedge:DimY,:)=UNDEF
     Var(1:XYedge,1:XYedge,:)=UNDEF
     Var(DimX-XYedge:DimX,DimY-XYedge:DimY,:)=UNDEF
  END IF


END SUBROUTINE Readinputfile

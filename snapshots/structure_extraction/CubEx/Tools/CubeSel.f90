PROGRAM CubeSel

  USE StatLib
  USE CubeLib
  IMPLICIT NONE
  CHARACTER(len=350) :: InpFile, OutputCube, OutputImage, idcube, string
  INTEGER(kind=4)    :: zmin, zmax, n_ext, is, i, j, k, xmax, ymax, zmin_large, zmax_large, bcsmooth, bcsmooth_z,id, zstart, zend, &
       boxmin(3), boxmax(3), xmin, ymin, idpad(3), cc_zmin, cc_zmax, DimX, DimY, DimZ
  LOGICAL :: OutVarCube, CFITSIOLibFix, BKGcorr, CONTsub, multiext, writeNaN, idclean, writeFloor, transpose
  REAL(kind=4)       :: bkgc, medtot, sigma, meanval, zbuf_blue, zbuf_red, this_med, ContClipVal, lmin, lmax, lmin_large, lmax_large, negnum, floor
  REAL(kind=4), ALLOCATABLE :: Smooth(:,:,:), CheckCube(:,:,:), CubeCut(:,:,:)
  REAL(kind=8) :: CheckCubeWCS(SIZE(WCS))

  Verbosity=2
  negnum=-1.0
  
  !.. collect input parameters
  CALL ReadCommandLine

  !.get number of extensions in the original file
  CALL GetNExt(InpFile,n_ext)
  
  !..issue a warning
  IF(n_ext==1) THEN
     IF(OutVarCube.or.multiext) THEN
        print *, " "
        print *, " WARNING: only one extension found in original file. Variance cube will NOT be produced."
        print *, " "
     END IF
     !..turn OutVarCube and multiext off
     OutVarCube=.false.
     multiext=.false.
  ELSE     
     !..force n_ext to 1 if no OutVarCube and no multiext is requested
     IF(.not.OutVarCube.and..not.multiext) n_ext=1
  END IF
     
  !..read idcube, if requested
  IF(id/=-1) THEN
     
     CALL ReadCube(IdCube)
     !..save WCS values for later
     CheckCubeWCS=WCS
     
     !..if padding is requested, change filename to a cube section in x and y
     !..override lmin, lmax and allocate a smaller checkcube
     IF(ANY(idpad/=9999)) THEN
        
        CALL FindBoundingBox

        xmin=MAX(1,boxmin(1)-idpad(1))
        xmax=MIN(SIZE(Cube,DIM=1),boxmax(1)+idpad(1))
        ymin=MAX(1,boxmin(2)-idpad(2))
        ymax=MIN(SIZE(Cube,DIM=2),boxmax(2)+idpad(2))
        cc_zmin=MAX(1,boxmin(3)-idpad(3))
        cc_zmax=MIN(SIZE(Cube,DIM=3),boxmax(3)+idpad(3))
        lmin=(cc_zmin-WCS(3))*WCS(11)+WCS(6)
        lmax=(cc_zmax-WCS(3))*WCS(11)+WCS(6)

        WRITE(string,'(4(a,i4.4),a)') "[",xmin,":",xmax,",",ymin,":",ymax,",*]"
        InpFile=TRIM(InpFile)//TRIM(string)
        print *, "updating InpFile name to: ",TRIM(InpFile)

        ALLOCATE(CheckCube(xmax-xmin+1,ymax-ymin+1,cc_zmax-cc_zmin+1))
        CheckCube(:,:,:)=Cube(xmin:xmax,ymin:ymax,cc_zmin:cc_zmax)

     ELSE

        ALLOCATE(CheckCube(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2),SIZE(Cube,DIM=3)))
        CheckCube=Cube

     END IF

  ELSE

     CheckCubeWCS=0.

  END IF
  
  !.. read and cut input datacube including zbuffer if requested
  IF(lmin==-1.and.lmax==-1) THEN !..cut using pixel units

     zmin_large=zmin-zbuf_blue
     IF(zmax/=-1) THEN
        zmax_large=zmax+zbuf_red
     ELSE 
        !..set a dummy large value (in order to read the whole cube)
        zmax_large=9999
     END IF

     CALL ReadCube(InpFile=InpFile, n_ext=n_ext, CFITSIOLibFix=CFITSIOLibFix, zmin=zmin_large, zmax=zmax_large)
        
  ELSE !..cut using lmin and lmax

     lmin_large=lmin-zbuf_blue
     IF(lmax/=-1) THEN
        lmax_large=lmax+zbuf_red
     ELSE 
        !..set a dummy large value (in order to read the whole cube)
        lmax_large=9999
     END IF

     CALL ReadCube(InpFile=InpFile, n_ext=n_ext, CFITSIOLibFix=CFITSIOLibFix, lmin=lmin_large, lmax=lmax_large)

     !..readjust zmin and zmax
     IF(lmin==-1) THEN
        zmin=1
     ELSE
        zmin=NINT((lmin-WCS(6))/WCS(11)+WCS(3))
     END IF
     IF(lmax==-1) THEN
        zmax=SIZE(Cube,DIM=3)+zmin-1
     ELSE
        zmax=MIN(SIZE(Cube,DIM=3),NINT((lmax-WCS(6))/WCS(11)+WCS(3)))+zmin-1
     END IF

     !..readjust zbuf to pixel units
     zbuf_blue=zbuf_blue/WCS(11)
     zbuf_red = zbuf_red/WCS(11)

  END IF

  !..remove continuum and/or background from Cube using median over extended wvl range
  IF(CONTsub.or.BKGcorr) CALL BKGSubtraction(ContSub=ContSub, ContClipVal=ContClipVal, zmin=zmin, zmax=zmax)

  !..performs a boxcar smoothing over a box with side 2*bcsmooth+1 in the spatial direction, and 2*bcsmooth_z+1 in the redshift direction
  IF(bcsmooth>0) CALL Boxcar(bcsmooth, bcsmooth_z=bcsmooth_z, zmin=zmin, zmax=zmax)

  
  !..set voxels not associated with selected id to UNDEF, if requested
  IF(id/=-1.and.idclean) THEN
     !..if Cube and CheckCube dimensions are different and if WCS are present, get initial pixel value of CheckCube in Cube
     IF(SIZE(Cube,DIM=3)>SIZE(CheckCube,DIM=3)) THEN
        IF(ALL(CheckCubeWCS==0.d0).and.ALL(WCS==0.d0)) &
             STOP "CheckCube and Cube have different z-dimension and WCS were not found to correct for that!"
        zstart=INT((CheckCubeWCS(6)-WCS(6))/WCS(11)+1)
        zend=zstart+SIZE(CheckCube,DIM=3)-1
        IF(zstart>SIZE(Cube,DIM=3).or.zend>SIZE(Cube,DIM=3)) STOP "mask outside of cube! Check mask or maskshift"
     ELSE
        IF(SIZE(Cube,DIM=3)<SIZE(CheckCube,DIM=3)) THEN
           IF(ANY(idpad/=9999)) print *, "WARNING: CheckCube is larger than datacube!"
        END IF
        zstart=1
        zend=SIZE(Cube,DIM=3)
     END IF
     !..if writeFloor is requested, get avg value outside of objects
     floor=Mean(PACK(Cube(:,:,zstart:zend),MASK=CheckCube==0.and.Cube(:,:,zstart:zend)/=UNDEF))
     !..remove/update voxels
     IF(writeFloor) THEN
        WHERE(CheckCube/=REAL(id)) Cube(:,:,zstart:zend)=floor
        IF(zstart>1) Cube(:,:,1:zstart-1)=floor
        IF(zend<SIZE(Cube,DIM=3)) Cube(:,:,zend+1:SIZE(Cube,DIM=3))=floor
     ELSE
        WHERE(CheckCube/=REAL(id)) Cube(:,:,zstart:zend)=UNDEF
        IF(zstart>1) Cube(:,:,1:zstart-1)=UNDEF
        IF(zend<SIZE(Cube,DIM=3)) Cube(:,:,zend+1:SIZE(Cube,DIM=3))=UNDEF
     END IF
  END IF
  
  !..if WCS are present, adjust WCS CRVAL3 value:
  !..NB: CRPIX1 anx CRPIX2 are automatically adjusted in the CFTSIO section
  IF(WCS(6)/=0.d0) THEN
     WCS(6)=WCS(6)+(zmin-WCS(3))*WCS(11)
  END IF

 
  !..convert zmax and zmin to the CutCube units, removing buffer zones if necessary (zmax-zmin+1 is the actual zdimension of the Cube we want to write)
  IF(zmax/=-1) THEN
     zmax=zbuf_blue+zmax-zmin+1
  ELSE
     zmax=SIZE(Cube,DIM=3)
  END IF
  zmin=zbuf_blue+1
     
  !.. produce whitelight image if requested, excluding buffer
  IF(TRIM(OutputImage)/="??") THEN
     IF(ALLOCATED(Image)) DEALLOCATE(Image)
     ALLOCATE(Image(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2)))
     Image=UNDEF
     Image(:,:)=SUM(Cube(:,:,zmin:zmax),MASK=Cube(:,:,zmin:zmax)/=UNDEF,DIM=3)
     !..write output image
     IF(writeNaN) THEN
        WHERE(Image==UNDEF.or.Image==0.0) Image=sqrt(negnum)
     ELSE
        WHERE(Image==UNDEF) Image=0.
     END IF
     CALL AssignWCS(type="Image")
     CALL WriteImage(OutputImage)
  END IF

  IF(ALLOCATED(CheckCube)) DEALLOCATE(CheckCube)

  !..if requested transpose the cube inverting X and Y
  IF(Transpose) THEN
     DimX=Size(CUBE,DIM=1); DimY=Size(CUBE,DIM=2);  DimZ=Size(CUBE,DIM=3)
     !..use CheckCube as a temporary array
     ALLOCATE(CheckCube(DimY,DimX,DimZ))
     DO j=1,DimY
        DO i=1,DimX
           CheckCube(j,i,:)=Cube(i,j,:)
        END DO
     END DO
     DEALLOCATE(Cube) ; ALLOCATE(Cube(DimY,DimX,DimZ))
     Cube=CheckCube
     DEALLOCATE(CheckCube)
  END IF
     

  !.. write output cube(s) if requested, excluding buffer
  IF(TRIM(OutputCube)/="??") THEN
     CALL WriteCube(OutputCube, n_ext=n_ext, zmin=zmin, zmax=zmax, multiext=multiext, writeNaN=writeNaN)
  END IF


  
CONTAINS

SUBROUTINE ReadCommandLine

  IMPLICIT NONE
  CHARACTER(len=300) :: ExeName, fname, string, ds, opt, arg 
  INTEGER :: narg, iarg, i, is
  LOGICAL :: ex

!..default
  InpFile="??"
  OutputCube="??"
  zmin=1
  zmax=-1    !..default will be adjusted to the size of the original datacube later, if zmax not selected
  lmin=-1
  lmax=-1
  CFITSIOLibFix=.true.
  OutVarCube=.true.
  OutputImage="??"
  BKGcorr=.false.
  CONTsub=.false.
  zbuf_blue=0
  zbuf_red=0
  ContClipVal=3.
  bcsmooth=0
  bcsmooth_z=-1 !..default value will be adjusted to the value of bcsmooth later, if bcsmooth/=0
  multiext=.true.
  writeNaN=.false.
  id=-1
  idcube="??"
  idpad=0
  idclean=.false.
  Transpose=.false.

  CALL GetArg(0,ExeName)
  narg=iargc()
  IF(narg>0) THEN
     CALL GetArg(1,fname)
  ELSE
     print *, " "
     WRITE(*,'(2a)')"        CubeSel (part of CubEx package)   "
     WRITE(*,'(a)') " SubCubes Selection software (and more! see options below) "
     WRITE(*,'(a)') " "
     WRITE(*,'(a)') " by Sebastiano Cantalupo (cantalupo@phys.ethz.ch)"
     print *, " "
     WRITE(*,'(a)')"usage: CubeSel.x -InpFile <name> [-option <val>]"
     WRITE(*,'(a)')" "
     WRITE(*,'(a)')"  ----- input/output options:"
     WRITE(*,'(a)')"  -InpFile          <name>           : datacube file name "
     WRITE(*,'(a)')"   -cube            <name>           : identical to -InpFile "
     WRITE(*,'(a)')"  -OutputCube       <name>           : name of the output cube file. If not requested, an OutputImage file name is required (see below)"
     WRITE(*,'(a)')"   -out             <name>           : identical to -OutputCube "
     WRITE(*,'(a)')"  -OutputImage      <name>           : if requested, produces a 'whitelight' image with the selected datacube"
     WRITE(*,'(a)')"   -outim           <name>           : identical to -OutputImage "
     WRITE(*,'(a)')"  -multiext         <.true./.false.> : if true (default) and original file is multiext writes a multiext file instead of two separate files" 
     WRITE(*,'(a)')"  -OutVarCube       <.true./.false.> : if true (default), variance is present and multiext=.false. produces a separate file with variance cube. "
     WRITE(*,'(a)')"                                       Filename is <outputcube>.VAR.fits"
     WRITE(*,'(a)')"  ----- selection by layers:"
     WRITE(*,'(a)')"  -zmin             <int>            : initial pixel in the z direction for datacube selection (default=1)"
     WRITE(*,'(a)')"  -zmax             <int>            : final pixel in the z direction for datacube selection (default=original cube size)"
     WRITE(*,'(a)')"  -lmin             <real>           : initial value in the z direction for datacube selection in CUNIT3 (e.g., Angstrom)"
     WRITE(*,'(a)')"  -lmax             <real>           : final value in the z direction for datacube selection in CUNIT3 (e.g. Angstrom)"
     WRITE(*,'(a)')"  "
     WRITE(*,'(a)')"  ----- selection by object id:"
     WRITE(*,'(a)')"  -idcube           <name>           : name of the 3d mask with id values (e.g., Objet_Id Checkcube produced by CubEx)."
     WRITE(*,'(a)')"                                       If not provided and -id is selected, the standard associated CubEx name will be used."
     WRITE(*,'(a)')"  -id               <int>            : selected obj id as in the 3d mask provided with the -idcube option or with standard CubEx name."
     WRITE(*,'(a)')"  -idclean          <bol>            : if .true., set to UNDEF all voxels not associated with this object id (default=.false.)"
     WRITE(*,'(a)')"  -idpad            <int int int>    : box padding around selected object in each direction (default=0). NB: larger than cube dim is allowed."
     WRITE(*,'(a)')"                                       Use -1 if you want to get everything in a specific direction or '-1 -1 -1' for whole cube."
     WRITE(*,'(a)')" "
     WRITE(*,'(a)')" ------ other options:"
!     WRITE(*,'(a)')"  -bkgcorr          <.true./.false.> : if true, subtracts the estimated residual bkg from datacube pixel by pixel (default=.false.)"
!     WRITE(*,'(a)')"                                       NB: this option is deprecated, use CubeBKGSub instead. "
!     WRITE(*,'(a)')"  -contsub          <.true./.false.> : if true, subtracts continuum objects from datacube (default=.false.) AND residual bkg"
!     WRITE(*,'(a)')"                                       NB: this option is deprecated, use CubeBKGSub instead."
!     WRITE(*,'(a)')"  -contclipval      <val>            : value in sigma to consider a pixel a continuum source rather than bkg residual (default=3.), only applies if contsub=.false."
!     WRITE(*,'(a)')"  -zbuf_blue        <val>            : extra buffer (pixels) for continuum and bkg subtraction on the blue side of selected range (default=0)"
!     WRITE(*,'(a)')"  -zbuf_red         <val>            : extra buffer (pixels) for continuum and bkg subtraction on the red side of selected range (default=0)"
     WRITE(*,'(a)')"  -bcsmooth         <int>            : if /=0, it does a boxcar average over a box with side 2*bcsmooth+1 pixels in the spatial "
     WRITE(*,'(a)')"                                        and z direction (default=0, i.e., no boxcar smoothing)"
     WRITE(*,'(a)')"  -bcsmooth_z       <int>            : boxcar smoothing lenght in the z direction if bcsmooth is requested (default=bcsmooth)"
     WRITE(*,'(a)')"  -writeNaN         <.true./.false.> : if .true. writes NaN in the cube instead of UNDEF (-999.0). NB: NaN handling may be compiler dependent! (default=.false.)"
     WRITE(*,'(a)')"  -writefloor       <.true./.false.> : if .true. uses the avg value of voxels outside of objects instead of UNDEF. Use for Objects_SNR_F (def=.false.)."
     WRITE(*,'(a)')"  -CFITSIOLibFix    <.true./.false.> : set this to .true. (default) if CFITSIO library crashes while reading the file " 
     stop
  END IF


  !..read from command line
  DO i=1,narg,2
     CALL getarg(i,opt)
     CALL getarg(i+1,arg)
     SELECT CASE(TRIM(opt))
     CASE('-InpFile')       ; READ(arg,'(a)') InpFile
     CASE('-cube')          ; READ(arg,'(a)') InpFile
     CASE('-zmin')          ; READ(arg,*) zmin
     CASE('-zmax')          ; READ(arg,*) zmax
     CASE('-lmin')          ; READ(arg,*) lmin
     CASE('-lmax')          ; READ(arg,*) lmax
     CASE('-zbuf_blue')     ; READ(arg,*) zbuf_blue
     CASE('-zbuf_red')      ; READ(arg,*) zbuf_red
     CASE('-CFITSIOLibFix') ; READ(arg,*) CFITSIOLibFix
     CASE('-OutputCube')    ; READ(arg,'(a)') OutputCube
     CASE('-out')           ; READ(arg,'(a)') OutputCube
     CASE('-OutVarCube')    ; READ(arg,*) OutVarCube
     CASE('-OutputImage')   ; READ(arg,'(a)') OutputImage
     CASE('-outim')         ; READ(arg,'(a)') OutputImage
     CASE('-bkgcorr')       ; READ(arg,*) BKGcorr
     CASE('-contsub')       ; READ(arg,*) CONTsub
     CASE('-contclipval')   ; READ(arg,*) ContClipVal
     CASE('-bcsmooth')      ; READ(arg,*) bcsmooth
     CASE('-bcsmooth_z')    ; READ(arg,*) bcsmooth_z
     CASE('-writeNaN')      ; READ(arg,*) writeNaN
     CASE('-multiext')      ; READ(arg,*) multiext
     CASE('-idcube')        ; READ(arg,'(a)') idcube
     CASE('-idpad')         ; READ(arg,*) idpad
     CASE('-id')            ; READ(arg,*) id
     CASE('-idclean')       ; READ(arg,*) idclean
     CASE('-writefloor')    ; READ(arg,*) writeFloor
     CASE('-transpose')     ; READ(arg,*) Transpose
     CASE default
        print *, "command line argument ",TRIM(opt), " not recognized!"
        STOP
     END SELECT
  END DO

!..perform few checks
  IF(TRIM(InpFile)=="??") THEN
     print *, "please provide the input datacube with the -InpFile or -cube option!"
     STOP
  ELSEIF(TRIM(OutputCube)=="??".and.TRIM(OutputImage)=="??") THEN
     STOP "please provide OutputCube and/or OutputImage filenames with the -OutputCube/-out and/or -OutputImage/-outim options!" 
  END IF

  IF(bcsmooth/=0.and.bcsmooth_z==-1) bcsmooth_z=bcsmooth

  IF(Id/=-1.and.TRIM(IdCube)=="??") THEN !..associate default name
     is=INDEX(TRIM(InpFile),".fits")
     IF(is==0) is=LEN_TRIM(InpFile)+1
     IdCube=TRIM(InpFile(1:is-1))//".Objects_Id.fits"
     INQUIRE(FILE=IdCube,EXIST=ex)
     IF(.not.ex) IdCube=TRIM(InpFile(1:is-1))//".Objects_Id_Assoc.fits"
     INQUIRE(FILE=IdCube,EXIST=ex)
     IF(.not.ex) STOP "Please provide the name of the associated checkcube with the objects Id with the option -idcube!"
  END IF


  !..adjust idpad to a very large number if option -1 was selected
  WHERE(idpad==-1) idpad=9999
 
  IF(BKGCorr) STOP "-bkgcorr is a depracated option, use CubeBKGSub instead!"
  IF(CONTsub) STOP "-contsub is a deprecated option, use CubeBKGSub instead!"
  
  
END SUBROUTINE ReadCommandLine

SUBROUTINE FindBoundingBox

  IMPLICIT NONE
  INTEGER(kind=4) :: i,j,k

  boxmin=100000
  boxmax=-100000
  
  DO k=1,SIZE(Cube,DIM=3)
     DO j=1,SIZE(Cube,DIM=2)
        DO i=1,SIZE(Cube,DIM=1)

           IF(Cube(i,j,k)==id) THEN
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

END SUBROUTINE FindBoundingBox

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


END PROGRAM CubeSel

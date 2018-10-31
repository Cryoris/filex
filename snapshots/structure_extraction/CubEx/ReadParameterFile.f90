SUBROUTINE ReadParameterFile

  USE Globalmodule
  IMPLICIT NONE
  CHARACTER(len=500) :: ExeName, fname, string, ds, opt, arg, comment, RescaleVarArea_string
  INTEGER(kind=4)    :: narg, myid, nthreads, ierr, i, ii, n, is, dims(3), j, unit, status, group, &
       in, end, rank, ival, naxes(3), nfound
  INTEGER            :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM, first_arg, val(1000), nval
  LOGICAL            :: ex, readparfile, anyf
  INTEGER(kind=4), ALLOCATABLE :: MaskOnlyArray(:), UnMaskArray(:)


  CALL GetArg(0,ExeName)
  narg=iargc()
  IF(narg>0) THEN
     CALL GetArg(1,fname)
  ELSE
     CALL Usage
     stop
  END IF

!..set default value for the mandatory fields
  InpFile="??"            !.. input file name of the datacube

!..check if parameter file exists and that syntax is correct
  INQUIRE(File=fname, Exist=ex)
  IF(ex) THEN
     !..check if this is the name of the input file instead...
     IF(INDEX(fname,".fits")>0) THEN
        !..interpret this as an input file and continue without a parameter file
        InpFile=TRIM(fname)
        readparfile=.false.
        first_arg=2
     ELSE
        readparfile=.true.
        first_arg=2
     END IF
  ELSE
     IF(SCAN(fname,"-")==0) THEN
        !..print a short Usage screen
        WRITE(*,'(a)') " "
        WRITE(*,'(a)') "please either provide a valid parameter file with the syntax:"
        print *, " "
        WRITE(*,'(a,1x,a)')TRIM(ExeName), " <parameter file> [-option [arg]] "
        print *, " "
        print *, "where command-line options have precedence over the parameter file value"
        print *, " "
        print *, "or provide at least the name of the imput file on the command line:"
        print *, " "
        WRITE(*,'(a,1x,a)')TRIM(ExeName), " [-InpFile] <name.fits> [-option [arg]] "
        print *, " "
        print *, "for the full list of options, type ", TRIM(ExeName)
        print *, " "
        STOP
     ELSE
        IF(TRIM(fname)=="-help".or.TRIM(fname)=="--help".or.TRIM(fname)=="-h".or.TRIM(fname)=="--h") THEN
           CALL Usage
           STOP
        END IF
        readparfile=.false.
        first_arg=1
     END IF
  END IF



  !..default values 
  MultiExt=.true.         !.. if .true. interprets InpFile as a multiextension datacube with data on ext 1 and variance on ext 2
  VarFile="??"            !.. name of the variance file if MultiExt=.false., if a name is provided for this file, then MultiExt is turned to .false.
                          !.. if MultiExt=.false. and VarFile is not provided, the variance will be estimated from the Cube (layer by layer).
  RescaleVar=.false.      !.. if .true., rescale the original variance in datacube by the estimated variance layer by layer (default=.false.)
  RescaleVarArea=-1       !.. array of 4 values (xmin,xmax,ymin,ymax), if all /=-1 rescales variance using only the area defined here
  RescaleVarArea_string="all"  !.. this is used to read the RescaleVarArea values, all=-1
  RescaleVarMin=0.9       !.. minimum rescaling factor for variance rescaling if RescaleVar=.true. (default=0.9)
  RescaleVarMax=10.       !.. maximum rescaling factor for variance rescaling if RescaleVar=.true. (default=10.)
  RescaleVarFR=150        !.. if >0, it will median-filter the variance rescaling factors using a filter with this radius (default=150)
                          !.. only values between RescaleVarMin and RescaleVarMax will be considered in the median filtering process.
  MinArea=1               !.. minimum number of spatial pixels for detection (projection over z)
  MinDz=1                 !.. minimum number of spectral pixels for detection (projection over x,y)
  MaxDz=10000             !.. max number of spectral pixels for detection (projection over x,y). Useful for removing continuum objs
  MinNSpax=30             !.. minimum number of voxels for detection. ALL three conditions above need to be met for detection
  SN_Threshold=3          !.. voxel individual threshold in SNR for detection
  Halo_Threshold=100      !.. voxel individual threshold in SNR for detection
  ISN_Threshold=0         !.. integrated threshold in SNR for detection
  SSN_Threshold=0         !.. integrated spectral threshold (using the optimally extracted 1d spectrum) in SNR for detection
  SN_Threshold_Conn=0     !.. voxel individual threshold in SNR for detection if connected to a previously detected voxel (default=0, not used)
  Catalogue="??"          !.. name of the final catalogue with objects, if not defined it will be InpFile with .cat extension instead of .fits
  ReliabCheck=.false.     !.. if .true. checks the reliability of the produced catalogue performing extraction and photometry on the "negative" Cube
  NegCatalogue="??"       !.. if ReliabCheck=.true., produces a file with the object detected in the 'negative catalogue' (default=no output file)
  CheckCube="??"          !.. name of CheckCubes requested. Default is "InpFile" + "CheckCubeType" + ".fits"
  NCheckCubes=1           !.. number of checkcubes requested (default=1), use 0 for no checkcubes
  Verbosity=2             !.. verbosity level (0=no screen output, 3=maximum verbosity)
  zmin=1                  !.. minimum spectral pixel for datacube selection (units: pixels)
  zmax=9999               !.. maximum spectral pixel for datacube selection      "
  lmin=0                  !.. minimum z-coordinate for datacube selection in CUNITS3, e.g. Angstrom, if =/0 overrides zmin (requires WCS in the fits!)
  lmax=0                  !.. maximum z-coordinate for datacube selection in CUNITS3, e.g. Angstrom, if =/0 overrides zmax (requires WCS in the fits!)
  XYedge=0                !.. trim XY edges of the cube by this amount (pixels)
  IdCube="??"             !.. if a filename is provided, performs only photometry using this cube with objects id skipping detection and extraction
  IdCubeOnlyList="-1"     !.. if/=-1, it will only use the selected objects from IdCube, ids must be all in the same line/string (e.g., IdCubeOnly = "35 12 47")
  InpCat="??"             !.. if a filename is provided, performs only photometry using the object from this CubEx catalogue skipping detection and extraction
  InpCatOnlyList="-1"     !.. if/=-1, it will only use the selected objects from InpCat, ids must be all in the same line/string (e.g., IdCubeOnly = "35 12 47")
  CheckCubeFMT="fits"     !.. format for CheckCube: fits/bov
  CheckCubeType="Objects_Id" !.. CheckCube variable; options: 
                          !..    Objects     = original cube pixels above detection threshold
                          !..    SNR         = Signal to Noise cube (simply Cube/sqrt(Var)) 
                          !..    Objects_SNR = SNR cube including only pixels of detected objects
                          !..    Objects_Id  = cube with the pixel labelled with the Id of detected objects
  Objects_SNR_Floor=0.99  !.. Floor value in units of SN_Threshold for CheckCubes Objects_SNR and Objects_SNR_F
  maxnlabels=1000000      !.. stack size for the maximum number of labels in the connected components labeling algorithm
  AperRadius=3            !.. Cylindrical aperture area. units: pixels
  AperDz=10               !.. Cylindrical aperture depth. units: spectral pixels 
  CFITSIOLibFix=.false.   !.. set this to .true. if CFITSIO library crashes while reading the file       
  ApplyFilter=.false.     !.. if true, applies a gaussian filter to the cube for extraction/detection (photometry is performed on original cube)
  ApplyFilterVar=.false.  !.. if true, applies a gaussian filter to the variance cube as well, assuming that covariance is 0
  FilterXYRad=0           !.. gaussian spatial filter size (radius) in pixels, 0=no filtering
  FilterZRad=0            !.. gaussian wavelenght filter size (radius) in spectral pixels, 0=no filtering
  ObjMaskList="??"        !.. ASCII file with list of (continuum) objects to be masked before extraction, detection and photometry. 
                          !....entries are "x y radius", one object per row. Mask is cylindrical over all datacube. Units: integer pixels.
  LayerMaskList="??"      !.. ASCII file with list of layers (one per row) to be masked before extraction, detection and photometry. Units: integer pixels.
  PrintHeader=.true.      !..print header in the catalogue file
  AssocFrac=-1.0          !..fraction of overlapping projected XY area to "associate" individual detections separated in wavelength, default means no association 
  AssocCatalogue="??"     !..name of assoc catalogue with merged objects
  Deblend_=.false.        !..flag to activate deblending
  Deblend_MinSNR=1.       !..minimum SNR threshold for deblending in units of SN_Threshold
  Deblend_MaxSNR=10.      !..maximum SNR threshold for deblending in units of SN_Threshold
  Deblend_NSteps=30       !..number of deblending thresholds between Deblend_Min_SNR and Deblend_Max_SNR
  Deblend_MinNPix=4       !..minimum number of spatial pixels above deblending threshold to trigger object splitting
  Deblend_BufferSize=100000 !..buffer size for the number of deblended object and deblending labels  
  SourceMaskName="??"     !..IMAGE with objects id of source pixels to mask before detection and extraction (e.g., the result of a CubEx run on the white-light image)
  MaskOnlyList="-1"       !..if/=-1, it will only mask the selected objects from SourceMask, ids must be all in the same line/string (e.g., MaskOnly = "35 12 47")
  UnMaskList="-1"         !..if/=-1, it will not consider the selected objects from SourceMask, ids must be all in the same line/string (e.g., UnMask = "23 145")
  RescalingVarOutFile="??" !..name of the output file with rescaling factors applied to the variance cube (see RescaleVar option).
  EstimatedVarOutFile="??"
  RescalingVarInpFile="??" !..if a file if provided, it will use the values found here (see RescaleVar option). A previous RescalingVarOutFile can be used.
  ApplyDataThreshold=.false.
  Threshold=[-1.e19,1.e19]
  
  !..read from parameter file, if provided
  IF(readparfile) THEN
     ierr=0
     OPEN(1,file=fname,action='read',form='formatted')
     DO
        READ(1,*,IOSTAT=ierr) string
        IF(ierr/=0) EXIT
        IF(SCAN(TRIM(string),'!#%$')>0) CYCLE 
        BACKSPACE(1)
        IF(Verbosity>=3) print *, "reading parameter:",TRIM(string)
        SELECT CASE(TRIM(string))  
        CASE('InpFile')       ; READ(1,*) ds,ds,InpFile
        CASE('VarFile')       ; READ(1,*) ds,ds,VarFile
        CASE('IdCube')        ; READ(1,*) ds,ds,IdCube
        CASE('idcube')        ; READ(1,*) ds,ds,IdCube
        CASE('IdCubeOnly')    ; READ(1,*) ds,ds,IdCubeOnlyList
        CASE('idcubeonly')    ; READ(1,*) ds,ds,IdCubeOnlyList
        CASE('InpCat')        ; READ(1,*) ds,ds,InpCat
        CASE('inpcat')        ; READ(1,*) ds,ds,inpcat
        CASE('InpCatOnly')    ; READ(1,*) ds,ds,InpCatOnlyList
        CASE('inpcatonly')    ; READ(1,*) ds,ds,InpCatOnlyList
        CASE('RescaleVar')    ; READ(1,*) ds,ds,RescaleVar
        CASE('RescaleVarArea'); READ(1,*) ds,ds,RescaleVarArea_string
        CASE('RescaleVarMin') ; READ(1,*) ds,ds,RescaleVarMin
        CASE('RescaleVarMax') ; READ(1,*) ds,ds,RescaleVarMax
        CASE('RescaleVarFR' ) ; READ(1,*) ds,ds,RescaleVarFR
        CASE('MultiExt')      ; READ(1,*) ds,ds,MultiExt
        CASE('ReliabCheck')   ; READ(1,*) ds,ds,ReliabCheck
        CASE('MinArea')       ; READ(1,*) ds,ds,MinArea
        CASE('MinDz')         ; READ(1,*) ds,ds,MinDz
        CASE('MaxDz')         ; READ(1,*) ds,ds,MaxDz
        CASE('MinNSpax')      ; READ(1,*) ds,ds,MinNSpax
        CASE('MinNVox')       ; READ(1,*) ds,ds,MinNSpax
        CASE('SN_Threshold')  ; READ(1,*) ds,ds,SN_Threshold
        CASE('Halo_Threshold')  ; READ(1,*) ds,ds,Halo_Threshold
        CASE('ISN_Threshold') ; READ(1,*) ds,ds,ISN_Threshold
        CASE('SSN_Threshold') ; READ(1,*) ds,ds,SSN_Threshold
        CASE('SN_Threshold_Conn') ; READ(1,*) ds,ds,SN_Threshold_Conn
        CASE('Catalogue')     ; READ(1,*) ds,ds,Catalogue
        CASE('AssocCatalogue'); READ(1,*) ds,ds,AssocCatalogue
        CASE('NegCatalogue')  ; READ(1,*) ds,ds,NegCatalogue
        CASE('NCheckCubes')   ; READ(1,*) ds,ds,NCheckCubes
        CASE('CheckCube') 
            READ(1,*) ds,ds,CheckCube(1)
            IF(TRIM(CheckCube(1))/="??") THEN
               BACKSPACE(1)
               READ(1,*) ds,ds,CheckCube(1:NCheckCubes)
            END IF
        CASE('Verbosity')     ; READ(1,*) ds,ds,Verbosity
        CASE('zmin')          ; READ(1,*) ds,ds,zmin
        CASE('zmax')          ; READ(1,*) ds,ds,zmax
        CASE('lmin')          ; READ(1,*) ds,ds,lmin
        CASE('lmax')          ; READ(1,*) ds,ds,lmax
        CASE('XYedge')        ; READ(1,*) ds,ds,XYedge
        CASE('CheckCubeFMT')  ; READ(1,*) ds,ds,CheckCubeFMT
        CASE('maxnlabels')    ; READ(1,*) ds,ds,maxnlabels
        CASE('AperRadius')    ; READ(1,*) ds,ds,AperRadius
        CASE('AperDz')        ; READ(1,*) ds,ds,AperDz
        CASE('CheckCubeType') ; READ(1,*) ds,ds,CheckCubeType(1:NCheckCubes)
        CASE('CFITSIOLibFix') ; READ(1,*) ds,ds,CFITSIOLibFix
        CASE('ApplyFilter')   ; READ(1,*) ds,ds,ApplyFilter
        CASE('ApplyFilterVar'); READ(1,*) ds,ds,ApplyFilterVar
        CASE('FilterXYRad')   ; READ(1,*) ds,ds,FilterXYRad
        CASE('FilterZRad')    ; READ(1,*) ds,ds,FilterZRad
        CASE('ObjMaskList')   ; READ(1,*) ds,ds,ObjMaskList
        CASE('LayerMaskList') ; READ(1,*) ds,ds,LayerMaskList
        CASE('PrintHeader')   ; READ(1,*) ds,ds,PrintHeader
        CASE('AssocFrac')     ; READ(1,*) ds,ds,AssocFrac
        CASE('Deblend')       ; READ(1,*) ds,ds,Deblend_
        CASE('Deblend_MinSNR'); READ(1,*) ds,ds,Deblend_MinSNR
        CASE('Deblend_MaxSNR'); READ(1,*) ds,ds,Deblend_MaxSNR
        CASE('Deblend_NSteps'); READ(1,*) ds,ds,Deblend_NSteps
        CASE('Deblend_MinNPix');READ(1,*) ds,ds,Deblend_MinNPix
        CASE('Deblend_Buffer'); READ(1,*) ds,ds,Deblend_BufferSize
        CASE('Objects_SNR_Floor'); READ(1,*) ds,ds,Objects_SNR_Floor
        CASE('SourceMask')    ; READ(1,*) ds,ds,SourceMaskName
        CASE('MaskOnly')      ; READ(1,*) ds,ds,MaskOnlyList
        CASE('UnMask')        ; READ(1,*) ds,ds,UnMaskList
        CASE('ApplyDataThreshold') ; READ(1,*) ds,ds,ApplyDataThreshold
        CASE('Threshold')     ; READ(1,*) Threshold(:)
        CASE('RescalingVarOutFile') ; READ(1,*) ds,ds,RescalingVarOutFile
        CASE('RescalingVarInpFile') ; READ(1,*) ds,ds,RescalingVarInpFile
        CASE('EstVarOutFile') ; READ(1,*) ds,ds,EstimatedVarOutFile
        CASE default
           print *, "parameter ",TRIM(string), " in parameter file not recognized!"
           IF(SCAN(string,"=")>0) THEN
              WRITE(*,'(a)') "NB: you must leave at least one blank space between the parameter name, the '=' and the value"
              print *, "  e.g.:  ", string(1:INDEX(string,"=")-1)," = ",string(INDEX(string,"=")+1:LEN_TRIM(string))
           END IF
           print *, "Please type", TRIM(EXEname), " to get the full list of valid options"
           STOP
        END SELECT
     END DO
     
  END IF

  !..read from command line
  ii=first_arg-1
  DO WHILE(ii<narg)
     ii=ii+1; CALL getarg(ii,opt)
     ii=ii+1; CALL getarg(ii,arg)
     IF(Verbosity>=3) print *, TRIM(opt), " ", TRIM(arg)
     SELECT CASE(TRIM(opt))
     CASE('-InpFile')       ; READ(arg,'(a)') InpFile
     CASE('-cube')          ; READ(arg,'(a)') InpFile
     CASE('-VarFile')       ; READ(arg,'(a)') VarFile
     CASE('-var')           ; READ(arg,'(a)') VarFile
     CASE('-IdCube')        ; READ(arg,'(a)') IdCube
     CASE('-idcube')        ; READ(arg,'(a)') IdCube
     CASE('-IdCubeOnly')    ; READ(arg,'(a)') IdCubeOnlyList
     CASE('-idcubeonly')    ; READ(arg,'(a)') IdCubeOnlyList
     CASE('-InpCat')        ; READ(arg,'(a)') InpCat
     CASE('-inpcat')        ; READ(arg,'(a)') InpCat
     CASE('-InpCatOnly')    ; READ(arg,'(a)') InpCatOnlyList
     CASE('-inpcatonly')    ; READ(arg,'(a)') InpCatOnlyList
     CASE('-RescaleVar')    ; READ(arg,*) RescaleVar
     CASE('-RescaleVarArea'); READ(arg,'(a)') RescaleVarArea_string    
     CASE('-RescaleVarMin') ; READ(arg,*) RescaleVarMin
     CASE('-RescaleVarMax') ; READ(arg,*) RescaleVarMax
     CASE('-RescaleVarFR')  ; READ(arg,*) RescaleVarFR
     CASE('-ApplyDataThreshold') ; READ(arg,*) ApplyDataThreshold
     CASE('-Threshold')    ; READ(arg,*) Threshold
     CASE('-MinArea')       ; READ(arg,*) MinArea
     CASE('-ReliabCheck')   ; READ(arg,*) ReliabCheck
     CASE('-MultiExt')      ; READ(arg,*) MultiExt
     CASE('-m')             ; READ(arg,*) MultiExt
     CASE('-MinDz')         ; READ(arg,*) MinDz
     CASE('-MaxDz')         ; READ(arg,*) MaxDz
     CASE('-MinNSpax')      ; READ(arg,*) MinNSpax
     CASE('-MinNVox')       ; READ(arg,*) MinNSpax
     CASE('-n')             ; READ(arg,*) MinNSpax
     CASE('-SN_Threshold')  ; READ(arg,*) SN_Threshold
     CASE('-Halo_Threshold'); READ(arg,*) Halo_Threshold
     CASE('-ISN_Threshold') ; READ(arg,*) ISN_Threshold
     CASE('-SSN_Threshold') ; READ(arg,*) SSN_Threshold
     CASE('-sn')            ; READ(arg,*) SN_Threshold
     CASE('-isn')           ; READ(arg,*) ISN_Threshold
     CASE('-ssn')           ; READ(arg,*) SSN_Threshold
     CASE('-SN_Threshold_Conn')  ; READ(arg,*) SN_Threshold_Conn
     CASE('-Catalogue')     ; READ(arg,'(a)') Catalogue
     CASE('-AssocCatalogue'); READ(arg,'(a)') AssocCatalogue
     CASE('-NegCatalogue')  ; READ(arg,'(a)') NegCatalogue
     CASE('-NCheckCubes')   ; READ(arg,*) NCheckCubes    
     CASE('-Verbosity')     ; READ(arg,*) Verbosity
     CASE('-zmin')          ; READ(arg,*) zmin
     CASE('-zmax')          ; READ(arg,*) zmax
     CASE('-lmin')          ; READ(arg,*) lmin
     CASE('-lmax')          ; READ(arg,*) lmax
     CASE('-XYedge')        ; READ(arg,*) XYedge
     CASE('-CheckCubeFMT')  ; READ(arg,*) CheckCubeFMT
     CASE('-maxnlabels')    ; READ(arg,*) maxnlabels
     CASE('-AperRadius')    ; READ(arg,*) AperRadius
     CASE('-AperDz')        ; READ(arg,*) AperDz
     CASE('-CFITSIOLibFix') ; READ(arg,*) CFITSIOLibFix
     CASE('-ApplyFilter')   ; READ(arg,*) ApplyFilter
     CASE('-f')             ; READ(arg,*) ApplyFilter
     CASE('-ApplyFilterVar'); READ(arg,*) ApplyFilterVar
     CASE('-fv')            ; READ(arg,*) ApplyFilterVar
     CASE('-FilterXYRad')   ; READ(arg,*) FilterXYRad
     CASE('-fsr')           ; READ(arg,*) FilterXYRad
     CASE('-FilterZRad')    ; READ(arg,*) FilterZRad
     CASE('-fzr')           ; READ(arg,*) FilterZRad
     CASE('-ObjMaskList')   ; READ(arg,*) ObjMaskList
     CASE('-LayerMaskList') ; READ(arg,*) LayerMaskList
     CASE('-PrintHeader')   ; READ(arg,*) PrintHeader
     CASE('-AssocFrac')     ; READ(arg,*) AssocFrac
     CASE('-Deblend')       ; READ(arg,*) Deblend_
     CASE('-Deblend_MinSNR'); READ(arg,*) Deblend_MinSNR
     CASE('-Deblend_MaxSNR'); READ(arg,*) Deblend_MaxSNR
     CASE('-Deblend_NSteps'); READ(arg,*) Deblend_NSteps
     CASE('-Deblend_MinNPix');READ(arg,*) Deblend_MinNPix
     CASE('-Deblend_Buffer'); READ(arg,*) Deblend_BufferSize
     CASE('-Objects_SNR_Floor'); READ(arg,*) Objects_SNR_Floor
     CASE('-SourceMask')    ; READ(arg,'(a)') SourceMaskName
     CASE('-sourcemask')    ; READ(arg,'(a)') SourceMaskName
     CASE('-MaskOnly')      ; READ(arg,'(a)') MaskOnlyList
     CASE('-maskonly')      ; READ(arg,'(a)') MaskOnlyList
     CASE('-UnMask')        ; READ(arg,'(a)') UnMaskList
     CASE('-unmask')        ; READ(arg,'(a)') UnMaskList
     CASE('-RescalingVarOutFile'); READ(arg,'(a)') RescalingVarOutFile
     CASE('-RescalingVarInpFile'); READ(arg,'(a)') RescalingVarInpFile
     CASE('-EstVarOutFile'); READ(arg,'(a)') EstimatedVarOutFile
     CASE('-CheckCube')     
        READ(arg,'(a)') CheckCube(1)
        DO n=2,NCheckCubes
           ii=ii+1
           CALL getarg(ii,arg)
           READ(arg,'(a)') CheckCube(n)
        END DO  
     CASE('-CheckCubeType')
         READ(arg,'(a)') CheckCubeType(1)       
        DO n=2,NCheckCubes
           ii=ii+1
           CALL getarg(ii,arg)
           READ(arg,'(a)') CheckCubeType(n)
        END DO  
     CASE('-cct')
         READ(arg,'(a)') CheckCubeType(1)       
        DO n=2,NCheckCubes
           ii=ii+1
           CALL getarg(ii,arg)
           READ(arg,'(a)') CheckCubeType(n)
        END DO  
     CASE default
        print *, "command line argument ",TRIM(opt), " not recognized!"
        print *, "Please type ", TRIM(EXEname), " to get the full list of valid options"
        STOP
     END SELECT
  END DO

!.. initial printout
  IF(Verbosity>=2) THEN
     WRITE(*,'(a)') "------------------------------------------------------------ "
     WRITE(*,'(2a)')"       CubExtractor       version: ", TRIM(version)
     WRITE(*,'(a)') " by Sebastiano Cantalupo (cantalupo@phys.ethz.ch)"
     WRITE(*,'(a)') " "
  END IF

! ---- performs few checks

  !..check for non-default parameters
  IF(TRIM(InpFile)=="??") STOP "Please provide InpFile name in the parameter file or with the command line option -InpFile <name> "

  !..associate default file names if necessary
  is=INDEX(TRIM(InpFile),".fits")
  IF(is==0) is=LEN_TRIM(InpFile)+1
  IF(TRIM(Catalogue)=="??") Catalogue=TRIM(InpFile(1:is-1))//".cat"
  IF(TRIM(AssocCatalogue)=="??".and.AssocFrac>0.) AssocCatalogue=TRIM(InpFile(1:is-1))//".assoc.cat"
  IF(NCheckCubes>0.and.TRIM(CheckCube(1))=="??") THEN
     DO i=1,NCheckCubes
        CheckCube(i)=TRIM(InpFile(1:is-1))//"."//TRIM(CheckCubeType(i))//".fits"
     END DO
  END IF

  !..general checks
  IF(ApplyFilter.and.FilterXYRad==0.and.FilterZRad==0) STOP "You must define a spatial (FilterXYRad) and/or spectral (FilterZRad) filter radius for filtering!"
  IF(ANY(CheckCubeType=="Objects_SNR_F").or.ANY(CheckCubeType=="SNR_F")) THEN
     IF(.not.ApplyFilter.or..not.ApplyFilterVar) THEN
        print *, "You cannot request smoothed checkcubes without performing smoothing on both data and variance!" 
        STOP "Please change the parameter file"
     END IF
  END IF
  IF(.not.MultiExt.and.TRIM(VarFile)=="??") THEN
     !STOP "You must provide the name of the VarFile if not contained in the InpFile (in this case, set MultiExt=.true.)"
     print *, "NB: VarFile not provided and MultiExt=.false. "
     print *, "---> Variance will be estimated layer by layer from the DATACUBE"
  END IF
  IF(TRIM(VarFile)/="??") MultiExt=.false.

  !..adjust Objects_SNR_Floor
  Objects_SNR_Floor=Objects_SNR_Floor*SN_Threshold

  !..read sourcemask if provided, check and split associated parameter strings if necessary
  !..and read/update SourceMask
  IF(TRIM(SourceMaskName)/="??") THEN

     IF(Verbosity>=2) print *, "reading SourceMask=",TRIM(SourceMaskName)

     !..open file
     status=0
     CALL ftgiou(unit,status)
     IF(status/=0) STOP "problem with ftgiou"
     CALL ftdopn(unit,SourceMaskName,0,status)
     IF(status/=0) STOP "problem with ftdopn"

     !..get cube rank
     status=0
     CALL ftgkyj(unit,'NAXIS',rank,comment,status)
     IF(status/=0) STOP "problem reading NAXIS keyword"

     !..read data
     status=0
     IF(rank==3) THEN
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
     ALLOCATE(Cube(naxes(1),naxes(2),naxes(3)),SourceMask(naxes(1),naxes(2)))
     group=1; in=1; end=PRODUCT(naxes)
     CALL ftgpve(unit,group,in,end,UNDEF,Cube,anyf,status)
     IF(status/=0) STOP "problem reading data"
     SourceMask(:,:)=INT(Cube(:,:,1))
     DEALLOCATE(Cube)

     !..close file and units
     CALL ftclos(unit,status)
     CALL ftfiou(unit,status)

     

     !..check and split associated parameter strings and update SourceMask if necessary
     IF(TRIM(MaskOnlyList)/="-1") THEN
        !..finds how many entry there are in the string
        DO i=1,1000
           ierr=0
           READ(MaskOnlyList,*,iostat=ierr) val(1:i)
           IF(ierr/=0) THEN
              nval=i-1
              EXIT
           END IF
        END DO
        IF(i>=1000) STOP "problem reading MaskOnly option!"
        !..read entries and apply values
        READ(MaskOnlyList,*) val(1:nval)
        DO j=1,SIZE(SourceMask,DIM=2)
           DO i=1,SIZE(SourceMask,DIM=1)
              IF(ALL(val(1:nval)/=SourceMask(i,j))) SourceMask(i,j)=0
           END DO
        END DO
     END IF
     IF(TRIM(UnMaskList)/="-1") THEN
        !..finds how many entry there are in the string
        DO i=1,1000
           ierr=0
           READ(UnMaskList,*,iostat=ierr) val(1:i)
           IF(ierr/=0) THEN
              nval=i-1
              EXIT
           END IF
        END DO
        IF(i>=1000) STOP "problem reading UnMask option!"
        !..read entries and apply values
        READ(UnMaskList,*) val(1:nval)
        DO j=1,SIZE(SourceMask,DIM=2)
           DO i=1,SIZE(SourceMask,DIM=1)
              IF(ANY(val(1:nval)==SourceMask(i,j))) SourceMask(i,j)=0
           END DO
        END DO
     END IF
 
  END IF

  !..perform a check
  IF(TRIM(InpCat)/="??".and.TRIM(IdCube)=="??") STOP "Please provide the IdCube associated with InpCat"
  IF(TRIM(InpCat)=="??".and.TRIM(IdCube)/="??") STOP "Please provide the InpCat associated with IdCube"
  IF(TRIM(InpCat)==TRIM(Catalogue)) STOP "Output Catalogue and InpCat have the same name! Please change output catalogue name."

  !..assign alias
  IF(IdCubeOnlyList=="-1".and.InpCatOnlyList/="-1") THEN
     IdCubeOnlyList=InpCatOnlyList
  END IF

  !..read RescaleVarArea values if necessary
  IF(TRIM(RescaleVarArea_string)/="all") THEN
     READ(RescaleVarArea_string,*) RescaleVarArea(1:4)
  END IF

CONTAINS

SUBROUTINE Usage

  IMPLICIT NONE
  LOGICAL :: ex


  !---------------------------------------------------------------------
  !..print some info
  WRITE(*,'(a)') "#------------------------------------------------------------ "
  WRITE(*,'(2a)')"#       CubExtractor       version: ", TRIM(version)
  WRITE(*,'(a)') "# "
  WRITE(*,'(a)') "# by Sebastiano Cantalupo (cantalupo@phys.ethz.ch)"
  WRITE(*,'(a)') "# "
  WRITE(*,'(a)') "# usage:"
  WRITE(*,'(2(a,1x),a)') "#", TRIM(ExeName), " <parameter file> [-option [arg]] "
  WRITE(*,'(a)') "# or (if you don't want to provide a parameter file):"
  WRITE(*,'(2(a,1x),a)') "#", TRIM(ExeName), " [-InpFile] <name> [-option [arg]] "
  WRITE(*,'(a)') "# "
  WRITE(*,'(a)') "# command line options:"
  WRITE(*,'(a)') "#  -<ParameterName> <ParameterValue> : overrides parameter values found in parameter file"
  WRITE(*,'(3a)')"#   e.g.:", TRIM(ExeName), " CubEx.par -Catalogue test.cat -MinNSPax 20 -SN_Threshold 2"
  WRITE(*,'(a)') "# "
  WRITE(*,'(a)') "# Here is the full list of available parameters with some explanation and default values:"
  WRITE(*,'(a)') "# "
  WRITE(*,'(a)') "# NB: you can dump this screen output to a file and use it as a parameter file!"
  WRITE(*,'(a)') "#     Parameter file rules: lines beginning with #%$! are comments. TAB not allowed. No specific order."
  WRITE(*,'(a)') "#                           At least one blank space must be present between the parameter name, '=' and the value!"
  WRITE(*,'(a)') "#     To restore default values, just comment with, e.g. #" 
  WRITE(*,'(a)') "# "
  WRITE(*,'(a)') " "
  WRITE(*,'(a)') "#--- Cube-related parameters:"
  WRITE(*,'(a)')'InpFile = "??"             !.. input file name of the datacube (NO DEFAULT, you need to provide one!) [short alias parameter "-cube"]'
  WRITE(*,'(a)')'Catalogue = "??"           !.. name of the final catalogue with objects (default= same as InpFile with .cat extension)'
  WRITE(*,'(a)')'zmin = 1                   !.. minimum spectral pixel for datacube selection (units: pixels)'
  WRITE(*,'(a)')'zmax = 9999                !.. maximum spectral pixel for datacube selection      "'
  WRITE(*,'(a)')'lmin = 0                   !.. minimum z-coordinate for datacube selection in CUNITS3, e.g. Angstrom, if =/0 overrides zmin (requires WCS in the fits!)'
  WRITE(*,'(a)')'lmax = 0                   !.. maximum z-coordinate for datacube selection in CUNITS3, e.g. Angstrom, if =/0 overrides zmax (requires WCS in the fits!)'
  print *, " "
  WRITE(*,'(a)') "#--- Variance-related parameters:"
  WRITE(*,'(a)')'VarFile = "??"             !.. name of the variance file if not contained in InpFile. Needed only if the file name is non standard (i.e., InpFile+".VAR.fits"). [alias "-var"]'
  WRITE(*,'(a)')'MultiExt = .true.          !.. if .true. interprets InpFile as a multiextension datacube with data on ext 1 and variance on ext 2 '
  WRITE(*,'(a)')'                           !... unless VarFile is provided or default VarFile name exists in the current directory. [alias "-m"]'
  WRITE(*,'(a)')'                           !... if .false. and VarFile=??, then a simple variance will be estimated from the cube itself, layer by layer.'
  WRITE(*,'(a)')'RescaleVar = .false.       !.. if .true., rescale the original variance in datacube using the estimated variance layer by layer'
  WRITE(*,'(a)')'RescaleVarArea = "all"     !.. area for variance rescaling, four integeer values="xmin xmax ymin ymax", use "all" for all FoV '
  WRITE(*,'(a)')'RescaleVarMin = 0.9        !.. minimum rescaling factor for variance rescaling if RescaleVar=.true.'
  WRITE(*,'(a)')'RescaleVarMax = 10.        !.. maximum rescaling factor for variance rescaling if RescaleVar=.true.'
  WRITE(*,'(a)')'RescaleVarFR = 150         !.. if >0, it will median-filter the variance rescaling factors using a filter with this radius'
  WRITE(*,'(a)')'                           !... only values between RescaleVarMin and RescaleVarMax will be considered in the median filtering.'
  WRITE(*,'(a)')'RescalingVarOutFile = "??" !.. if provided, writes down in this file the rescaling factor of the variance (see RescaleVar option)'
  WRITE(*,'(a)')'RescalingVarInpFile = "??" !.. if provided, uses the rescaling factor of the variance from this file (see RescaleVar option)'
  WRITE(*,'(a)')'EstVarOutFile = "??"       !.. if provided, writes down in this file the estimated variance (see MulitExt option)'   
  WRITE(*,'(a)') " "
  WRITE(*,'(a)') "#--- Filtering parameters:"
  WRITE(*,'(a)')'ApplyFilter = .false.      !.. if true, applies a gaussian filter to the cube for detection (photometry is performed on original cube) [alias "-f"] '
  WRITE(*,'(a)')'ApplyFilterVar = .false.   !.. if true, applies the same gaussian filter as above aslo to the variance cube for detection (assuming null covariance) [alias "-fv"]'
  WRITE(*,'(a)')'FilterXYRad = 0            !.. gaussian spatial filter size (radius) in pixels, 0=no filtering [alias "-fsr"]'
  WRITE(*,'(a)')'FilterZRad = 0             !.. gaussian wavelenght filter size (radius) in spectral pixels, 0=no filtering [alias "-fzr"]'
  print *, " "
  WRITE(*,'(a)')"#--- Detection parameters:"
  WRITE(*,'(a)')'SN_Threshold = 3           !.. voxel individual threshold in SNR for detection [alias "-sn"] '
  WRITE(*,'(a)')'Halo_Threshold = 100       !.. voxel individual threshold in SNR for detection [alias "-sn"] '
  WRITE(*,'(a)')'ISN_Threshold = 0          !.. integrated SNR over the 3d-mask for detection (0=not used) [alias "-isn"]'
  WRITE(*,'(a)')'SSN_Threshold = 0          !.. integrated spectral SNR (using 1d optimally extracted spectra) for detection (0=not used) [alias "-ssn"]'
  WRITE(*,'(a)')'SN_Threshold_Conn = 0      !.. voxel individual threshold in SNR for detection if connected to a previously detected voxel (0=not used)'
  WRITE(*,'(a)')'MinNVox = 30               !.. minimum number of voxels for detection. ALL three conditions above need to be met for detection [alias "-n" and "MinNSpax"] '
  WRITE(*,'(a)')'MinArea = 1                !.. minimum number of spatial pixels for detection (projection over z)'
  WRITE(*,'(a)')'MinDz = 1                  !.. minimum number of spectral pixels for detection (projection over x,y)'
  WRITE(*,'(a)')'MaxDz = 10000              !.. max number of spectral pixels for detection (projection over x,y). Useful for removing continuum objs'
  WRITE(*,'(a)')'ReliabCheck = .false.      !.. if .true. checks the reliability of the catalogue performing extraction and photometry on the "negative" Cube'
  WRITE(*,'(a)')'NegCatalogue = "??"        !.. if set and if ReliabCheck=.true., it produces a file with the objects detected in the "negative" Cube'
  print *, " "
  WRITE(*,'(a)')"#--- Deblending and Merging parameters:"
  WRITE(*,'(a)')'Deblend = .false.          !.. flag to activate deblending '
  WRITE(*,'(a)')'Deblend_MinSNR = 1.        !.. minimum SNR threshold for deblending in units of SN_Threshold '
  WRITE(*,'(a)')'Deblend_MaxSNR = 10.       !.. maximum SNR threshold for deblending in units of SN_Threshold '
  WRITE(*,'(a)')'Deblend_NSteps = 30        !.. number of deblending thresholds between Deblend_Min_SNR and Deblend_Max_SNR '
  WRITE(*,'(a)')'Deblend_MinNPix = 4        !.. minimum number of spatial pixels above deblending threshold to trigger object splitting '
  WRITE(*,'(a)')'AssocFrac = -1             !.. fraction of overlapping projected XY area to "associate/merge" individual detections separated' 
  WRITE(*,'(a)')'                           !... in wavelength, default means no association (and no assoc catalogue produced)'
  WRITE(*,'(a)')'AssocCatalogue = "??"      !.. name of the merged catalogue with associated objects (default= same as InpFile with assoc.cat extension)'
  print *, " "
  WRITE(*,'(a)')'#--- Photometry parameters:'
  WRITE(*,'(a)')'AperRadius = 3             !.. Cylindrical aperture area. units: pixels'
  WRITE(*,'(a)')'AperDz = 0                 !.. If defined > 0 : Cylindrical aperture spectral depth around light-weighted centroid. units: spectral pixels.'
  WRITE(*,'(a)')'                           !.. If defined < 0 : automatic spectral extraction is performed above abs(AperDz) sigma using 1d spectra.'
  WRITE(*,'(a)')'                           !.. If defined = 0 : simply use the spectral width from the segmentation map.'
  WRITE(*,'(a)')'InpCat = "??"              !.. if a filename is provided, performs only photometry using the positions of the objects'
  WRITE(*,'(a)')'                           !... in this CubEx catalogue skipping detection and extraction, a IdCube (see below) is also needed for the IsoFlux measurement'
  WRITE(*,'(a)')'InpCatOnly = "-1"          !.. if defined and /=-1, it will only perform photometry for the selected objects in InpCat (e.g., InpCatOnly = "35 12 47")'
  WRITE(*,'(a)')'IdCube = "??"              !.. associated IdCube to the InpCat for the IsoFlux measurement'
!  WRITE(*,'(a)')'IdCubeOnly = "-1"          !.. if defined and /=-1, it will only perform photometry for the selected objects in IdCube (e.g., IdCubeOnly = "35 12 47")'
  print *, " "
  WRITE(*,'(a)')"#--- CheckCube parameters (NB: these are arrays if NCheckCubes>1):"
  WRITE(*,'(a)')'NCheckCubes = 1              !.. number of CheckCubes produced by CubEx (1 by default), use 0 for no checkcubes'
  WRITE(*,'(a)')'                             !.. NB: if /=1, NCheckCubes parameter MUST preced any other CheckCube options in the parameter file!'
  WRITE(*,'(a)')'CheckCube = "??"             !.. name of CheckCube(s) if requested. Default is "InpFile"+"CheckCubeType"+".fits". NB: all names on the same line!'
  WRITE(*,'(a)')'CheckCubeFMT = "fits"        !.. format for CheckCube(s): either fits/bov'
  WRITE(*,'(a)')'CheckCubeType = "Objects_Id" !.. CheckCube(s) variable; options (one per checkcube requested, all on the same line): [alias "-cct"]'
  WRITE(*,'(a)')'                             !... Objects          = original cube pixels associated with detected objects'
  WRITE(*,'(a)')'                             !... Objects_F        = filtered cube pixels associated with detected objects'
  WRITE(*,'(a)')'                             !... Objects_SNR      = SNR cube including only pixels of detected objects'
  WRITE(*,'(a)')'                             !... Objects_SNR_F    = filtered SNR cube including only pixels of detected objects'
  WRITE(*,'(a)')'                             !... Objects_Id       = cube with the pixel labeled with the Id of detected objects (3d segmentation map)'
  WRITE(*,'(a)')'                             !... Objects_Id_Assoc = cube with the pixel labeled with the Id or Assoc Id in case of associated objects'
  WRITE(*,'(a)')'                             !... SNR              = Signal to Noise cube (simply Cube/sqrt(Var)) '
  WRITE(*,'(a)')'                             !... SNR_F            = filtered Signal to Noise cube (simply CubeF/sqrt(VarF)) '
  WRITE(*,'(a)')'                             !... Residuals        = original cube pixels that are not in detected objects'
  WRITE(*,'(a)')'                             !... EXAMPLE: set NCheckCubes = 3 and CheckCubeType = "Objects_Id" "Objects_SNR_F" "SNR_F" '
  WRITE(*,'(a)')'Objects_SNR_Floor = 0.99     !.. Floor value **in units of SN_Threshold** for CheckCubes Objects_SNR and Objects_SNR_F, '
  WRITE(*,'(a)')'                             !.. i.e., all voxels not associated to objects will have this value in the SNR CheckCubes.'
  WRITE(*,'(a)')'                             !.. Parameter for visualization purposes if Contours on VisIt are used (because of VisIt interpolation).'   
  print *, " "
  WRITE(*,'(a)')"#--- Masking parameters:"
  WRITE(*,'(a)')'SourceMask = "??"          !.. IMAGE with objects id of sources to mask before detection (e.g., the result of a CubEx run on the white-light image)'
  WRITE(*,'(a)')'MaskOnly   = "-1"          !.. if defined and /=-1, it will only mask the selected objects in sourcemask (e.g., MaskOnly = "35 12 47")'
  WRITE(*,'(a)')'UnMask     = "-1"          !.. if defined and /=-1, it will not  mask the selected objects in sourcemask (e.g., UnMask = "23 145")'
  WRITE(*,'(a)')'XYedge = 0                 !.. mask FOV edges of the cube by this amount (pixels), useful for cubes with noisy edges'
  WRITE(*,'(a)')'LayerMaskList = "??"       !.. ASCII file with list of layers (one per row) to be masked before extraction, detection and photometry. '
  WRITE(*,'(a)')'ObjMaskList = "??"         !.. ASCII file with list of pixels to be masked before extraction, detection and photometry. '
  WRITE(*,'(a)')'                           !... entries are "x y radius", one object per row. Mask is cylindrical over all datacube. Units: integer pixels.'
  print *, " "
  WRITE(*,'(a)')"#--- OTHER:" 
  WRITE(*,'(a)')'Verbosity = 2              !.. verbosity level (0=no screen output, 3=maximum verbosity)'
  WRITE(*,'(a)')'PrintHeader = .true.       !.. print header in the catalogue file'
  WRITE(*,'(a)')'maxnlabels = 1000000       !.. stack size for the maximum number of labels in the connected components labeling algorithm'
  WRITE(*,'(a)')'Deblend_Buffer = 100000    !.. stack size for the number of deblended object and deblending labels'  
  WRITE(*,'(a)')'CFITSIOLibFix = .false.    !.. set this to .true. if CFITSIO library crashes while reading a subsection of the datacube'   
  print *, " "

END SUBROUTINE Usage

END SUBROUTINE ReadParameterFile

MODULE Globalmodule

!..flag parameter
  INTEGER(kind=4), PARAMETER :: HaloId=-1
  INTEGER(kind=4) :: HaloLabel

!..global arrays
  TYPE Object_type
     INTEGER(kind=4) :: Id, area, dz, boxmin(3), boxmax(3), NSpax, Assoc
     REAL(kind=4)    :: IsoFlux, IsoErr, xcen(3), lcen(3), BoxFlux, BoxErr, AperFlux, AperErr
     ! NB: xcen= geometrical centroid, lcen= centroid weighted by flux
  END type Object_type
  TYPE(Object_type), ALLOCATABLE :: Obj(:)
  REAL(kind=4), ALLOCATABLE :: Cube(:,:,:), Var(:,:,:), CubeF(:,:,:), VarF(:,:,:)
  INTEGER(kind=4), ALLOCATABLE :: Mask(:,:,:), SourceMask(:,:)

  !..global variables/parameters
  CHARACTER(len=10), PARAMETER  :: version="1.7"
  CHARACTER(len=500) :: InpFile, Catalogue, CheckCube(10), CheckCubeFMT, CheckCubeType(10), ObjMaskList, VarFile, NegCatalogue, AssocCatalogue, LayerMaskList, &
       SourceMaskName, MaskOnlyList, UnMaskList, EstimatedVarOutFile, RescalingVarOutFile, IdCube, RescalingVarInpFile, IdCubeOnlyList, InpCat, InpCatOnlyList
  REAL(kind=4)       :: InitCPUTime, AssocFrac, lmin, lmax, Objects_SNR_Floor, Threshold(2), RescaleVarMin, RescaleVarMax
  INTEGER(kind=4)    :: Verbosity, zmin, zmax, DimX, DimY, DimZ, maxnlabels, AperRadius, AperDz, FilterXYRad, FilterZRad, NCheckCubes, XYedge, RescaleVarFR, RescaleVarArea(4)
  INTEGER(kind=4), ALLOCATABLE :: LabelToId(:)              ! auxiliary array that gives the Id in Obj array that correspond to a given label in Mask
  INTEGER(kind=4), ALLOCATABLE :: IdToLabel(:)              ! auxiliary array that gives the label in Mask from the corresponding Id in Obj array
  REAL(kind=4), PARAMETER :: UNDEF=-999.0                   ! value for undefined pixels
  LOGICAL :: CFITSIOLibFix                                  ! flag to correct for a bug in some version of CFITSIO Libraries
  LOGICAL :: ApplyFilter, ApplyFilterVar, PrintHeader, MultiExt, ReliabCheck, Deblend_, RescaleVar, ApplyDataThreshold                                    
  CHARACTER(len=200) :: WCSlabels(11)=["CRPIX1","CRPIX2","CRPIX3","CRVAL1","CRVAL2","CRVAL3","CD1_1 ","CD1_2 ","CD2_1 ","CD2_2 ","CD3_3 "]
  REAL(kind=8) :: WCS(11)
  CHARACTER(len=200) :: WCSlabels_strings(6)=["CTYPE1","CTYPE2","CTYPE3","CUNIT1","CUNIT2","CUNIT3"]
  CHARACTER(len=8) :: WCSstrings(6)=["RA---TAN","DEC--TAN","LINEAR  ","deg     ","deg     ","Angstrom"]  !..default value for MUSE datacubes

!..Detection parameters
  INTEGER(kind=4) :: MinArea, MinDz, MinNSpax, MaxDz, Deblend_NSteps, Deblend_MinNPix, Deblend_BufferSize
  REAL(kind=4)    :: SN_Threshold, Halo_Threshold, Deblend_MinSNR, Deblend_MaxSNR, SN_Threshold_Conn, ISN_Threshold, SSN_Threshold


END MODULE Globalmodule
  

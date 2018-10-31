MODULE CubeLib

  IMPLICIT NONE
  REAL(kind=4), ALLOCATABLE :: Cube(:,:,:), Var(:,:,:), Image(:,:) 
  REAL(kind=4), PARAMETER :: UNDEF=-999.0                   ! value for undefined pixels
  INTEGER(kind=4)    :: Verbosity
  CHARACTER(len=200) :: WCSlabels(11)=["CRPIX1","CRPIX2","CRPIX3","CRVAL1","CRVAL2","CRVAL3","CD1_1 ","CD1_2 ","CD2_1 ","CD2_2 ","CD3_3 "]!, imWCSlabels(8)=["CRPIX1","CRPIX2","CRVAL1","CRVAL2","CD1_1 ","CD1_2 ","CD2_1 ","CD2_2 "]
  CHARACTER(len=200) :: WCSlabels_strings(6)=["CTYPE1 ","CTYPE2 ","CTYPE3 ","CUNIT1 ","CUNIT2 ","CUNIT3 "]!, &
      !imWCSlabels_strings(4)=["CTYPE1 ","CTYPE2 ","CUNIT1 ","CUNIT2 "] 
  CHARACTER(len=8) :: fits_varname(2) 
  REAL(kind=8) :: WCS(11)
  REAL(kind=8), ALLOCATABLE :: imWCS(:)
  CHARACTER(len=8) :: WCSstrings(6)
  CHARACTER(len=8), ALLOCATABLE :: imWCSstrings(:), imWCSlabels(:), imWCSlabels_strings(:)
  CHARACTER(len=4) :: version="v1.7"

  
CONTAINS

!------------------------------------

RECURSIVE SUBROUTINE ReadCube(InpFile, n_ext, CFITSIOLibFix, zmin, zmax, lmin, lmax, reallocate, LocalCube, lun, read_var_only)

  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN):: InpFile
  INTEGER(kind=4), OPTIONAL, INTENT(IN) :: n_ext
  LOGICAL, OPTIONAL, INTENT(IN):: CFITSIOLibFix, reallocate, read_var_only
  INTEGER(kind=4), OPTIONAL, INTENT(IN) :: zmin, zmax, lun
  REAL(kind=4), OPTIONAL, INTENT(IN) :: lmin, lmax
  REAL(kind=4), OPTIONAL, INTENT(OUT) :: LocalCube(:,:,:)
  INTEGER(kind=4) :: status, unit, rwstatus, blocksize, naxes(3), group, in, end, ext, zmin_, zmax_, DimX, DimY, DimZ, nfound, n_ext_, i, rank, ival, start_ext, is
  CHARACTER(len=350) :: fname, comment, section
  LOGICAL :: anyf, reallocate_
  REAL(kind=4), ALLOCATABLE :: DATA(:,:,:)
  REAL(kind=8) :: dum

  !..checks for optional arguments:

  IF(present(zmin)) THEN
     zmin_=zmin
  ELSE
     zmin_=1
  END IF

  IF(present(n_ext)) THEN
     n_ext_=n_ext
  ELSE
     n_ext_=1 !..read only first extension ("Cube") by default
  END IF

  IF(present(reallocate)) THEN
     reallocate_=reallocate
  ELSE
     reallocate_=.true.
  END IF

  IF(present(lun)) THEN
     unit=lun
  ELSE
     unit=-1
  END IF

  IF(present(read_var_only)) THEN
     IF(read_var_only) THEN
        n_ext_=2
        start_ext=2
     ELSE
        start_ext=1
     END IF
  ELSE
     start_ext=1
  END IF

  !..check if filename contains cube sections and save string separately
  is=INDEX(TRIM(InpFile),".fits")
  section=TRIM(InpFile(is+5:))
  !print *, "section=",TRIM(section)


  DO ext=start_ext,n_ext_

     status=0

     !..get an unused unit
     IF(unit<0) CALL ftgiou(unit,status)

     IF(status/=0) STOP "problem with ftgiou"

     IF(ext==1) THEN
        !..open first extension: DATA
        IF(TRIM(section)/="[1]".and.TRIM(section)/="[2]") THEN        
           fname=TRIM(InpFile(1:is+4))//"[1]"//TRIM(section)
        ELSE
           fname=TRIM(InpFile)
        END IF
     ELSE
        !..open second extension: STAT
        IF(TRIM(section)/="[1]".and.TRIM(section)/="[2]") THEN
           fname=TRIM(InpFile(1:is+4))//"[2]"//TRIM(section)
        ELSE
           fname=TRIM(InpFile)
        END IF
     END IF
     rwstatus=0   ! 0=readonly
     CALL ftdopn(unit,fname,rwstatus,status)
     IF(status/=0) THEN
        IF(present(read_var_only)) THEN
           IF(read_var_only) print *, "problem reading variance in file: ", TRIM(fname)
           STOP
        END IF
        !..try to open a single fits file with no extension, 
        !..unless we were supposed to read the variance in ext=2 here.
        IF(ext==2.and.start_ext==1) THEN
           print *, "problem reading variance in file:",TRIM(fname)
           STOP
        END IF        
        fname=TRIM(InpFile)
        status=0
        CALL ftclos(unit,status)
        IF(unit<0) CALL ftgiou(unit,status)
        CALL ftdopn(unit,fname,rwstatus,status)
        IF(status/=0) THEN
           print *, "problem reading file: ", TRIM(fname)
           STOP 
        END IF
     END IF

     IF(Verbosity>=2) print *, "Reading:",TRIM(fname)

     !..get cube rank
     CALL ftgkyj(unit,'NAXIS',rank,comment,status)

     IF(rank==3) THEN
        !..get cube size
        CALL ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
     ELSEIF(rank==2) THEN
        !print *, "reading as image..."
        !..get image size
        CALL ftgkyj(unit,'NAXIS1',ival,comment,status)
        naxes(1)=ival
        CALL ftgkyj(unit,'NAXIS2',ival,comment,status)
        naxes(2)=ival
        naxes(3)=1
        !print *, "naxes=",naxes
     ELSE
        print *, "input file is not a cube or a image according to NAXIS keyword:", rank
        STOP
     END IF

     !..read variable name
     CALL ftgkys(unit,"OBJECT",fits_varname(ext),comment,status)
     IF(status/=0) THEN
        IF(Verbosity>=1) print *, "Warning: fits variable name not found, using standard name"
        fits_varname(ext)="hdu1"
     END IF

     IF(ext==start_ext) THEN !..we only need to do this once:

        !..read and store WCS information, if present

        !..check that keywords are there...
        status=0
        CALL ftgkyd(unit,"CRPIX1",dum,comment,status)
        
        IF(status==0) THEN !read keywords
           DO i=1,SIZE(WCS)
              CALL ftgkyd(unit,TRIM(WCSlabels(i)),WCS(i),comment,status)
              IF(status/=0) THEN
                 IF(Verbosity>=1) print *, "WARNING:: WCS keyword ", TRIM(WCSlabels(i)), " not found!"
                 IF(TRIM(WCSlabels(i))=="CD3_3") THEN !..try to read "CDELT3" instead
                    status=0
                    CALL ftgkyd(unit,"CDELT3",WCS(i),comment,status)
                    IF(status==0) THEN
                      IF(Verbosity>=1) print *, "Using CDELT3 instead"
                    ELSE
                       IF(Verbosity>=1) print *, "CDELT3 not found, using MUSE default value of CD3_3=1.25"
                       status=0 
                       WCS(i)=1.25
                       !WCS(11)=1. !.. that's the scale factor for z-direction
                    END IF
                 ELSE
                    status=0 
                    WCS=0
                    WCS(11)=1. !.. that's the scale factor for z-direction                   
                 END IF
              END IF
           END DO
           status=0
           DO i=1,SIZE(WCSstrings)
              CALL ftgkys(unit,TRIM(WCSlabels_strings(i)),WCSstrings(i),comment,status)
              IF(status/=0) THEN
                 IF(Verbosity>=1) print *, "WARNING:: WCS string keyword ", TRIM(WCSlabels_strings(i)), " not found!"
                 WCSstrings(:)="??"
                 status=0 
              END IF
           END DO

           !..adjust zmin and zmax in CUNIT3 if requested
           IF(present(lmin)) THEN
                 zmin_=NINT((lmin-WCS(6))/WCS(11)+WCS(3))
           END IF
           IF(present(lmax)) THEN
                 zmax_=MIN(naxes(3),NINT((lmax-WCS(6))/WCS(11)+WCS(3)))
           END IF
              
        ELSE

           IF(Verbosity>=2) THEN
              print *, " "
              print *, "NB: no WCS info found in input file"
              print *, " "
           END IF
           !..reset WCS if necessary
           IF(present(read_var_only)) THEN
              IF(.not.read_var_only) THEN
                 WCS(:)=0.d0
              ENDIF
           ELSE
              WCS(:)=0.d0
           END IF
           status=0

        END IF

        IF(present(zmax).and..not.present(lmax)) THEN
           zmax_=MIN(zmax,naxes(3))
        ELSEIF(.not.present(lmax)) THEN
           zmax_=naxes(3)
        ELSE
           IF(ALL(WCS==0.d0)) zmax_=naxes(3)
        END IF
     
        IF(Verbosity>=2) THEN
           print *, "Cube dimensions =", naxes
           print *, "Selected subset, zmin=", zmin_,"zmax=",zmax_,"pixels"
        END IF

        !..dimension without ghost zones
        DimX=naxes(1)
        DimY=naxes(2)
        DimZ=zmax_-zmin_+1

        IF(zmin_/=1.or.zmax_/=naxes(3)) THEN  !..read a subset with a temporary array

           IF(present(CFITSIOLibFix)) THEN
              IF(.not.CFITSIOLibFix) THEN !..see note below
                 IF(ALLOCATED(DATA)) DEALLOCATE(DATA)
                 ALLOCATE(DATA(DimX,DimY,DimZ))
                 in=(zmin_-1)*naxes(1)*naxes(2)+1
              ELSE
                 IF(ALLOCATED(DATA)) DEALLOCATE(DATA)
                 ALLOCATE(DATA(DimX,DimY,1:zmax_))
                 in=1   ! some version of CFITSIO library have problems starting from in/=1, so we have to read the array to the end and cut it
              END IF
           ELSE
              IF(ALLOCATED(DATA)) DEALLOCATE(DATA)
              ALLOCATE(DATA(DimX,DimY,1:zmax_)) !..use Fix by default
              in=1   ! some version of CFITSIO library have problems starting from in/=1, so we have to read the array to the end and cut it
           END IF

        ELSE

           in=1

        END IF


        IF(.not.present(LocalCube)) THEN
           IF(reallocate_) THEN
              IF(ALLOCATED(Cube)) DEALLOCATE(Cube)
              ALLOCATE(Cube(DimX,DimY,DimZ))
              IF(n_ext_==2) THEN
                 IF(ALLOCATED(Var)) DEALLOCATE(Var)
                 ALLOCATE (Var(DimX,DimY,DimZ))
                 Var=0.
              END IF
              Cube=0.
           END IF
        END IF
          

        end=zmax_*naxes(1)*naxes(2)
        
     END IF

     group=1

     !..read cube subset
     IF(zmin_/=1.or.zmax_/=naxes(3)) THEN     

        CALL ftgpve(unit,group,in,end,UNDEF,DATA,anyf,status)

        IF(status/=0) THEN
           CALL ftclos(unit,status)
           IF(present(CFITSIOLibFix)) THEN
              IF(CFITSIOLibFix) STOP "Problem reading datacube, and CFITSIOLibFix doesn't fix the problem...!"
           END IF
           !..let's try to fix the lib and start again...
           !DEALLOCATE(Cube,DATA)
           IF(present(LocalCube)) THEN
              IF(present(read_var_only)) THEN
                 CALL ReadCube(InpFile, n_ext=n_ext_, CFITSIOLibFix=.true., zmin=zmin_, zmax=zmax_, LocalCube=LocalCube, lun=unit+200, read_var_only=read_var_only)
              ELSE
                 CALL ReadCube(InpFile, n_ext=n_ext_, CFITSIOLibFix=.true., zmin=zmin_, zmax=zmax_, LocalCube=LocalCube, lun=unit+200)
              END IF
           ELSE
              CALL ReadCube(InpFile, n_ext=n_ext_, CFITSIOLibFix=.true., zmin=zmin_, zmax=zmax_, lun=MAX(111,unit+200))
           END IF
           !..if succeeded we can just exit here
           EXIT
        END IF
 
        !..copy over
        IF(.not.present(LocalCube)) THEN
           IF(ext==1) THEN
              Cube(1:DimX,1:DimY,1:DimZ)=DATA(1:DimX,1:DimY,zmin_:zmax_)
           ELSE
              Var(1:DimX,1:DimY,1:DimZ)=DATA(1:DimX,1:DimY,zmin_:zmax_)
           END IF
        ELSE
           LocalCube(1:DimX,1:DimY,1:DimZ)=DATA(1:DimX,1:DimY,zmin_:zmax_)
        END IF

        
    ELSE  !..read the whole cube or variance directly
        
       IF(.not.present(LocalCube)) THEN
          IF(ext==1) THEN
             CALL ftgpve(unit,group,in,end,UNDEF,Cube,anyf,status)          
             IF(status/=0) STOP "Problem reading Cube"
          ELSE
             CALL ftgpve(unit,group,in,end,UNDEF,Var,anyf,status)
             IF(status/=0) STOP "Problem reading Var"
          END IF
       ELSE
          CALL ftgpve(unit,group,in,end,UNDEF,LocalCube,anyf,status)          
          IF(status/=0) STOP "Problem reading LocalCube"
       END IF

    END IF

     IF(Verbosity>=2.and.anyf) print *, "NB: undefined pixels are present"

     CALL ftclos(unit,status)
     IF(status/=0) STOP "problem with ftclos!"
     
  END DO

  IF(ALLOCATED(DATA)) DEALLOCATE(DATA)

  !..free all allocated units
  call ftfiou(-1,status)
  IF(status/=0) STOP "problem with ftfiou"

END SUBROUTINE ReadCube

!-----------------------------------------------

SUBROUTINE WriteCube(OutputCube, n_ext, zmin, zmax, multiext, writeNaN, author)

  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN) :: OutputCube
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: author
  INTEGER(kind=4), INTENT(IN), OPTIONAL :: n_ext, zmin, zmax
  LOGICAL, INTENT(IN), OPTIONAL :: multiext
  CHARACTER(len=250) :: varfile, datafile
  INTEGER :: status, naxis, naxes(3), group, blocksize, bitpix, unit, nelements, fpixel, ext, rwstatus, n_ext_, is, zmin_, zmax_, i
  LOGICAL :: ex, extend, simple, multiext_, writeNaN_
  LOGICAL, INTENT(IN), OPTIONAL :: writeNaN
  REAL(kind=4) :: negnum

  IF(present(multiext)) THEN
     multiext_=multiext
  ELSE
     multiext_=.false.
  END IF
  IF(present(n_ext)) THEN
     n_ext_=n_ext
     IF(multiext_) n_ext_=1
  ELSE
     n_ext_=1 !..write only first extension ("Cube") by default
  END IF

  !..set UNDEF values to NaN, unless requested otherwise
  negnum=-1.0
  IF(present(writeNaN)) THEN
     writeNaN_=writeNaN
  ELSE
     writeNaN_=.true. !..DEFAULT
  END IF
  IF(writeNaN_) THEN
     WHERE(Cube==UNDEF) Cube=sqrt(negnum)
     IF(ALLOCATED(Var)) THEN
        WHERE(Var==UNDEF) Var=sqrt(negnum)
     END IF
  END IF


  zmin_=1
  IF(present(zmin)) zmin_=zmin
  zmax_=SIZE(Cube,DIM=3)
  IF(present(zmax)) zmax_=zmax

  is=INDEX(TRIM(OutputCube),".fits")

  
  IF(n_ext_==2) THEN
     IF(is==0) THEN
        varfile=TRIM(OutputCube)//".VAR.fits"
        datafile=TRIM(OutputCube)//".fits"
     ELSE
        varfile=TRIM(OutputCube(1:is-1))//".VAR.fits"
        datafile=TRIM(OutputCube)
     END IF
     print *, "Writing Output Cubes: ", TRIM(datafile), " ", TRIM(varfile)
  ELSE
     IF(is==0) THEN
        datafile=TRIM(OutputCube)//".fits"
     ELSE
        datafile=TRIM(OutputCube)
     END IF
     print *, "Writing Output Cube: ", TRIM(datafile)
  END IF

  DO ext=1,n_ext_

     status=0

     !..get an unused unit
     CALL ftgiou(unit,status)
     
     IF(status/=0) STOP "problem with ftgiou"
     
     rwstatus=1  ! 1=readwrite
     INQUIRE(file=datafile,EXIST=ex)
     IF(ex) THEN !..delete file
        CALL ftopen(unit,datafile,rwstatus,blocksize,status)
        IF(status/=0) THEN
           print *, "problem with ftopn! Try to remove the file ",TRIM(datafile), " manually..."
           STOP
        END IF
        CALL ftdelt(unit, status)
        IF(status/=0) STOP "problem with ftdelt"
     ENDIF
     blocksize=1
     CALL ftinit(unit,datafile,blocksize,status)
     IF(status/=0) STOP "problem with ftinit"

     
     !  Initialize parameters about the FITS image.
     IF(multiext_) THEN
        simple=.true.
        extend=.true.
        bitpix=8
        naxis=0
        call ftphpr(unit,simple,bitpix,0,1,0,1,extend,status)
        if(status/=0) STOP "problem with keywords writing for extension 0!"
        !..update/create author keyword if requested
        IF(present(author)) CALL ftukys(unit,"AUTHOR  ",TRIM(author), "origin of the file", status)
     END IF


     simple=.true.
     extend=.true.
     bitpix=-32
     naxis=3
     naxes(:)=SHAPE(Cube(:,:,zmin_:zmax_))
     group=1
     fpixel=1
     nelements=PRODUCT(naxes)     

     IF(multiext_) THEN
        call ftiimg(unit,bitpix,naxis,naxes,status)
        if(status/=0) STOP "problem with ftiimg for extension 1"
     END IF

     !  Write the required header keywords to the file
     IF(.not.multiext_) THEN
        call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
        IF(status/=0) STOP "problem with ftphpr"
        !..update/create author keyword if requested
        IF(present(author)) CALL ftukys(unit,"AUTHOR  ",TRIM(author), "origin of the file", status)
     END IF

     !.. if WCS were present in the original file, write them here:
     IF(ANY(WCS/=0.d0)) THEN
        DO i=1,SIZE(WCS)
           CALL ftpkyd(unit,TRIM(WCSlabels(i)),WCS(i),12,' ',status)
        END DO
        DO i=1,SIZE(WCSstrings)
           CALL ftpkys(unit,TRIM(WCSlabels_strings(i)),WCSstrings(i),' ',status)
        END DO
     END IF

     !..write variable names
     CALL ftpkys(unit,"OBJECT",TRIM(fits_varname(ext)),' ',status)
     IF(status/=0) THEN
        print *, "WARNING: problem writing fits_varname"
        status=0
     END IF

     !  Write the array to the FITS file.
     IF(ext==1) THEN
        call ftppre(unit,group,fpixel,nelements,Cube(:,:,zmin_:zmax_),status)
        IF(multiext_) THEN
           call ftiimg(unit,bitpix,naxis,naxes,status) !..insert a new HDU
           call ftppre(unit,group,fpixel,nelements,Var(:,:,zmin_:zmax_),status)
        END IF
     ELSEIF(ext==2) THEN
        call ftppre(unit,group,fpixel,nelements,Var(:,:,zmin_:zmax_),status)
     END IF
     IF(status/=0) STOP "problem with ftppre"
     
     call ftclos(unit, status)
     IF(status/=0) STOP "problem with ftclos!"
     
     !..change name to Var cube if necessary for next iteration
     datafile=TRIM(varfile)
     
  END DO

  call ftfiou(-1,status)
  IF(status/=0) STOP "problem with ftfiou"
  
END SUBROUTINE WriteCube

!-------------------------------------------

!-----------------------------------------------

SUBROUTINE WriteLocalCube(LocalCube, OutputCube, zmin, zmax, author)

  IMPLICIT NONE
  REAL(kind=4), INTENT(IN) :: LocalCube(:,:,:)
  CHARACTER(len=*), INTENT(IN) :: OutputCube
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: author
  INTEGER(kind=4), INTENT(IN), OPTIONAL :: zmin, zmax
  CHARACTER(len=250) :: varfile, datafile
  INTEGER :: status, naxis, naxes(3), group, blocksize, bitpix, unit, nelements, fpixel, ext, rwstatus, n_ext_, is, zmin_, zmax_, i
  LOGICAL :: ex, extend, simple
  REAL(kind=4) :: negnum

  zmin_=1
  IF(present(zmin)) zmin_=zmin
  zmax_=SIZE(LocalCube,DIM=3)
  IF(present(zmax)) zmax_=zmax

  datafile=OutputCube

  status=0

  !..get an unused unit
  CALL ftgiou(unit,status)
     
  IF(status/=0) STOP "problem with ftgiou"
     
  rwstatus=1  ! 1=readwrite
  INQUIRE(file=datafile,EXIST=ex)
  IF(ex) THEN !..delete file
     CALL ftopen(unit,datafile,rwstatus,blocksize,status)
     IF(status/=0) THEN
        print *, "problem with ftopn! Try to remove the file ",TRIM(datafile), " manually..."
        STOP
     END IF
     CALL ftdelt(unit, status)
     IF(status/=0) STOP "problem with ftdelt"
  ENDIF
  blocksize=1
  CALL ftinit(unit,datafile,blocksize,status)
  IF(status/=0) STOP "problem with ftinit"

  simple=.true.
  extend=.true.
  bitpix=-32
  naxis=3
  naxes(:)=SHAPE(LocalCube(:,:,zmin_:zmax_))
  group=1
  fpixel=1
  nelements=PRODUCT(naxes)     

  !  Write the required header keywords to the file
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  IF(status/=0) STOP "problem with ftphpr"

  !..update/create author keyword if requested
  IF(present(author)) CALL ftucrd(unit,"AUTHOR  ",TRIM(author), status)

  !.. if WCS were present in the original file, write them here:
  IF(ANY(WCS/=0.d0)) THEN
     DO i=1,SIZE(WCS)
        CALL ftpkyd(unit,TRIM(WCSlabels(i)),WCS(i),12,' ',status)
     END DO
     DO i=1,SIZE(WCSstrings)
        CALL ftpkys(unit,TRIM(WCSlabels_strings(i)),WCSstrings(i),' ',status)
     END DO
  END IF
  
  !  Write the array to the FITS file.
  call ftppre(unit,group,fpixel,nelements,LocalCube(:,:,zmin_:zmax_),status)
  IF(status/=0) STOP "problem with ftppre"
     
  call ftclos(unit, status)
  IF(status/=0) STOP "problem with ftclos!"

  call ftfiou(-1,status)
  IF(status/=0) STOP "problem with ftfiou"
  
END SUBROUTINE WriteLocalCube

!----------------------------------


SUBROUTINE UpdateCube(InpFile, OutFile, writeNaN, author, nHDUs, checknHDUs)

  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN) :: InpFile, OutFile
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: author
  INTEGER(kind=4), INTENT(OUT), OPTIONAL :: nHDUs
  LOGICAL, INTENT(IN), OPTIONAL :: checknHDUs
  INTEGER :: status, unit, inunit, hdutype, blocksize, or_blocksize, rwstatus, nelements, naxis, nh, nHDUs_
  LOGICAL, INTENT(IN), OPTIONAL :: writeNaN
  LOGICAL :: ex, writeNaN_, checknHDUs_
  REAL(kind=4) :: negnum
  CHARACTER(len=80) :: comment

  IF(VERBOSITY>=2) print *, "Writing output file: ",TRIM(OutFile), " using header information from: ",TRIM(InpFile)

  !..set UNDEF values to NaN, unless requested otherwise
  negnum=-1.0
  IF(present(writeNaN)) THEN
     writeNaN_=writeNaN
  ELSE
     writeNaN_=.true. !..DEFAULT
  END IF
  IF(writeNaN_) THEN
     WHERE(Cube==UNDEF) Cube=sqrt(negnum)
     IF(ALLOCATED(Var)) THEN
        WHERE(Var==UNDEF) Var=sqrt(negnum)
     END IF
  END IF

  status=0

  !..get an unused unit for the input file
  CALL ftgiou(inunit,status)
  IF(status/=0) STOP "problem with ftgiou"

  !..open input file 
  CALL ftopen(inunit, InpFile, 0, or_blocksize, status)
  IF(status/=0) STOP "problem opening InpFile"

  !..check how many HDU there are in the inpfile, if requested
  CALL ftthdu(inunit, nHDUs_, status)
  IF(status/=0) STOP "problem retrieving nHDU"

  IF(present(nHDUs)) nHDUs=nHDUS_

  !..get an unused unit for the output 
  CALL ftgiou(unit,status)   
  IF(status/=0) STOP "problem with ftgiou"

  !..check if the output file exists, in that case, delete the file
  rwstatus=1  ! 1=readwrite
  INQUIRE(file=OutFile,EXIST=ex)
  IF(ex) THEN !..delete file
     CALL ftopen(unit,OutFile,rwstatus,blocksize,status)
     IF(status/=0) THEN
        print *, "problem with ftopn! Try to remove the file ",TRIM(OutFile), " manually..."
        STOP
     END IF
     CALL ftdelt(unit, status)
     IF(status/=0) STOP "problem with ftdelt"
  ENDIF
    
  !..open new file 
  CALL ftinit(unit,OutFile,or_blocksize,status)
  IF(status/=0) STOP "problem with ftinit"


  IF(nHDUs_>1) THEN !..multiextension file (MUSE kind)

     !..move to the extension 1 in inpfile
     CALL ftmahd(inunit, 1, hdutype, status)
     IF(status/=0) STOP "problem moving to 1st extension in InpFile"

     !..move to extension 1 of OutFile
     CALL ftmahd(unit, 1, hdutype, status)
     IF(status/=0) STOP "problem moving to 1st extension in OutFile"
     
     !print *, "copying main header..."

     !..copy over the header from first extension
     CALL ftcphd(inunit,unit,status)
     IF(status/=0) STOP "problem copying header"

     !..update/create author keyword if requested
     IF(present(author)) CALL ftukys(unit,"AUTHOR  ",TRIM(author), "origin of the file", status)

     !..move to extension 2 in InpFile
     CALL ftmahd(inunit, 2, hdutype, status)
     IF(status/=0) STOP "problem moving to 2nd extension in InpFile"

     !..create a new extension in OutFile
     CALL ftcrhd(unit, status)
     IF(status/=0) STOP "problem creating 2nd extension in OutFile"

     !print *, "copying header of hdu=2..."

     !..copy the header
     CALL ftcphd(inunit,unit,status)
     IF(status/=0) STOP "problem copying header of ext 2"  

     !..update NAXIS keywords
     CALL ftukyj(unit, "NAXIS1", SIZE(Cube,DIM=1), "Axis lenght", status)
     IF(status/=0) STOP "problem updating NAXIS1 keyword"
     CALL ftukyj(unit, "NAXIS2", SIZE(Cube,DIM=2), "Axis lenght", status)
     IF(status/=0) STOP "problem updating NAXIS2 keyword"
     CALL ftgkyj(unit, "NAXIS", naxis, comment, status)
     IF(status/=0) STOP "problem reading NAXIS keyword"
     IF(naxis==3) THEN
        CALL ftukyj(unit, "NAXIS3", SIZE(Cube,DIM=3), "Axis lenght", status)
        IF(status/=0) STOP "problem updating NAXIS3 keyword"
     END IF
     
     !print *, "copying data of hdu=2..."
     
     !..write the data
     nelements=SIZE(Cube,DIM=1)*SIZE(Cube,DIM=2)*SIZE(Cube,Dim=3)
     call ftppre(unit,1,1,nelements,Cube,status)
     IF(status/=0) STOP "problem writing Cube"
     
     IF(ALLOCATED(Var)) THEN

        !..move to extension 3
        CALL ftmahd(inunit, 3, hdutype, status)
        IF(status/=0) STOP "problem moving to 3rd extension in InpFile"
        
        !..create a new extension in OutFile
        CALL ftcrhd(unit, status)
        IF(status/=0) STOP "problem creating 3rd extension in OutFile" 
        
        !print *, "copying header of hdu=3..."
        
        !..copy the header
        CALL ftcphd(inunit,unit,status)
        IF(status/=0) STOP "problem copying header of ext 3"  
        
        !print *, "copying data of hdu=3..."

        !..update NAXIS keywords
        CALL ftukyj(unit, "NAXIS1", SIZE(Cube,DIM=1), "Axis lenght", status)
        IF(status/=0) STOP "problem updating NAXIS1 keyword"
        CALL ftukyj(unit, "NAXIS2", SIZE(Cube,DIM=2), "Axis lenght", status)
        IF(status/=0) STOP "problem updating NAXIS2 keyword"
        CALL ftgkyj(unit, "NAXIS", naxis, comment, status)
        IF(status/=0) STOP "problem reading NAXIS keyword"
        IF(naxis==3) THEN
           CALL ftukyj(unit, "NAXIS3", SIZE(Cube,DIM=3), "Axis lenght", status)
           IF(status/=0) STOP "problem updating NAXIS3 keyword"
        END IF

        !..write the data
        call ftppre(unit,1,1,nelements,Var,status)
        IF(status/=0) STOP "problem writing Var"

     END IF

  ELSE  !..single extension file with header in the primary HDU

     !..copy over the header from first extension
     CALL ftcphd(inunit,unit,status)
     IF(status/=0) STOP "problem copying header"

     !..update NAXIS keywords
     CALL ftukyj(unit, "NAXIS1", SIZE(Cube,DIM=1), "Axis lenght", status)
     IF(status/=0) STOP "problem updating NAXIS1 keyword"
     CALL ftukyj(unit, "NAXIS2", SIZE(Cube,DIM=2), "Axis lenght", status)
     IF(status/=0) STOP "problem updating NAXIS2 keyword"
     CALL ftgkyj(unit, "NAXIS", naxis, comment, status)
     IF(status/=0) STOP "problem reading NAXIS keyword"
     IF(naxis==3) THEN
        CALL ftukyj(unit, "NAXIS3", SIZE(Cube,DIM=3), "Axis lenght", status)
        IF(status/=0) STOP "problem updating NAXIS3 keyword"
     END IF

     !..update/create author keyword if requested
     IF(present(author)) CALL ftukys(unit,"AUTHOR  ",TRIM(author), "origin of the file", status)     

     !..write the data
     nelements=SIZE(Cube,DIM=1)*SIZE(Cube,DIM=2)*SIZE(Cube,Dim=3)
     call ftppre(unit,1,1,nelements,Cube,status)
     IF(status/=0) STOP "problem writing Cube"
    
  END IF

  !..close files
  CALL ftclos(unit,status)
  IF(status/=0) STOP "problem closing files"
  CALL ftclos(inunit,status)
  IF(status/=0) STOP "problem closing files"

  call ftfiou(-1,status)
  IF(status/=0) STOP "problem with ftfiou"

END SUBROUTINE UpdateCube

  


!------------------------------------------

SUBROUTINE WriteImage(OutputImage,imageunits)

  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN) :: OutputImage
  INTEGER :: status, naxis, naxes(2), group, blocksize, bitpix, unit, nelements, fpixel, ext, rwstatus, i, j
  LOGICAL :: ex, extend, simple
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: imageunits

  IF(Verbosity>=2) print *, "Writing 2D image file ", TRIM(OutputImage)

  IF(.not.ALLOCATED(Image)) STOP "Image array has not been allocated!"

  status=0

  !..get an unused unit
  CALL ftgiou(unit,status)
     
  IF(status/=0) STOP "problem with ftgiou"
     
  rwstatus=1  ! 1=readwrite
  INQUIRE(file=OutputImage,EXIST=ex)
  IF(ex) THEN !..delete file
     CALL ftdopn(unit,OutputImage,rwstatus,status)
     IF(status/=0) STOP "problem with ftdopn"
     CALL ftdelt(unit, status)
     IF(status/=0) STOP "problem with ftdelt"
  ENDIF
  blocksize=1
  CALL ftinit(unit,OutputImage,blocksize,status)
  IF(status/=0) STOP "problem with ftinit"
          
  !  Initialize parameters about the FITS image.
  simple=.true.
  extend=.true.
  bitpix=-32
  naxis=2
  naxes(:)=SHAPE(Image)
     
  !  Write the required header keywords to the file
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  IF(status/=0) STOP "problem with ftphpr"

  !..update/create units keyword if requested
  IF(present(imageunits)) CALL ftukys(unit,"BUNIT   ",TRIM(imageunits), "units", status)

  !.. if image WCS are passed through the global module, write them here:
  IF(ALLOCATED(imWCS)) THEN
     DO i=1,SIZE(imWCS)
        CALL ftpkyd(unit,TRIM(imWCSlabels(i)),imWCS(i),12,' ',status)
     END DO
     DO i=1,SIZE(imWCSstrings)
        CALL ftpkys(unit,TRIM(imWCSlabels_strings(i)),imWCSstrings(i),' ',status)
     END DO
  END IF
  
  
  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=PRODUCT(naxes)
  call ftppre(unit,group,fpixel,nelements,Image,status)
  IF(status/=0) STOP "problem with ftppre"
  
  call ftclos(unit, status)
  IF(status/=0) STOP "problem with ftclos!"

  CALL ftfiou(-1, status)
  IF(status/=0) STOP "problem with ftfiou!"


END SUBROUTINE WriteImage

!------------------------------------------------------------------------

SUBROUTINE BKGSubtraction(ContSub, ContClipVal, zmin, zmax)

  USE StatLib
  IMPLICIT NONE
  LOGICAL, INTENT(IN), OPTIONAL :: ContSub
  REAL, INTENT(IN), OPTIONAL :: ContClipVal
  INTEGER, INTENT(IN), OPTIONAL :: zmin, zmax
  INTEGER(kind=4) :: zmin_, zmax_, i, j, k
  REAL(kind=4) :: meanval, medtot, sigma, ContClipVal_, this_med
  REAL(kind=4), ALLOCATABLE :: BKG(:,:)
  LOGICAL :: CONTsub_

  !..default values for optional parameter
  ContSub_=.true.
  ContClipVal_=3.
  zmin_=1
  zmax_=SIZE(Cube,DIM=3)

  !..adjust values if necessary
  IF(present(ContSub)) ContSub_=ContSub
  IF(present(ContClipVal)) ContClipVal_=ContClipVal
  IF(present(zmin)) zmin_=zmin
  IF(present(zmax)) zmax_=zmax

  IF(CONTsub_) THEN !..that's the easiest case, we consider continuum and bkg the same

     print *, "removing bkg residual and continuum sources..."
     DO j=1,SIZE(Cube,DIM=2)
        DO i=1,SIZE(Cube,DIM=1)
           IF(ALL(Cube(i,j,:)/=UNDEF)) THEN
              CALL SigmaClip(PACK(Cube(i,j,:),.true.), meanval, this_med, sigma)
              Cube(i,j,zmin_:zmax_)=Cube(i,j,zmin_:zmax_)-this_med
           ELSE
              Cube(i,j,zmin_:zmax_)=UNDEF
           END IF
        END DO
     END DO

  ELSE !..we will separate background from continuum pixels using sigclip

     print *, "removing bkg residual..."
     !..first, estimate the stddev 
     ALLOCATE(BKG(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2)))
     DO j=1,SIZE(Cube,DIM=2)
        DO i=1,SIZE(Cube,DIM=1)
           IF(ALL(Cube(i,j,:)/=UNDEF)) THEN
              CALL SigmaClip(PACK(Cube(i,j,:),.true.), meanval, this_med, sigma)
              BKG(i,j)=this_med
           ELSE
              BKG(i,j)=UNDEF
           END IF
        END DO
     END DO
     CALL SigmaClip(PACK(BKG,MASK=(BKG/=UNDEF)), meanval, medtot, sigma)
     print *, "meanval=",meanval,"medtot=",medtot,"sigma=",sigma
     !..put undef bkg values to 0
     WHERE(BKG==UNDEF) BKG=0.
     !..remove continuum objects from estimated bkg
     WHERE(BKG-meanval>ContClipVal_*sigma) BKG=0.
     !..subtract bkg
     FORALL(i=zmin_:zmax_) Cube(:,:,i)=Cube(:,:,i)-BKG(:,:)

     !..free memory
     DEALLOCATE(BKG)

  END IF

END SUBROUTINE BKGSubtraction

!-------------------------------

SUBROUTINE BoxCar(bcsmooth, bcsmooth_z, zmin, zmax)

  USE StatLib
  IMPLICIT NONE
  INTEGER(kind=4),INTENT(IN) :: bcsmooth
  INTEGER(kind=4),INTENT(IN), OPTIONAL :: bcsmooth_z
  INTEGER, INTENT(IN), OPTIONAL :: zmin, zmax
  INTEGER(kind=4) :: zmin_, zmax_, i, j, k, bcsmoothz, xmax, ymax
  REAL(kind=4), ALLOCATABLE :: Smooth(:,:,:)


  !..default values for optional parameter
  bcsmoothz=bcsmooth
  zmin_=1
  zmax_=SIZE(Cube,DIM=3)

  !..adjust values if necessary
  IF(present(bcsmooth_z)) bcsmoothz=bcsmooth_z
  IF(present(zmin)) zmin_=zmin
  IF(present(zmax)) zmax_=zmax

  print *, "boxcar smoothing..."

  ALLOCATE(Smooth(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2),zmin_:zmax_))
  ymax=SIZE(Cube,DIM=2)
  xmax=SIZE(Cube,DIM=1)

  !print *, xmax, ymax, zmin_, zmax_, bcsmooth, bcsmoothz

  DO k=zmin_,zmax_
     DO j=1,ymax
        DO i=1,xmax

           IF(ANY(Cube(MAX(1,i-bcsmooth):MIN(xmax,i+bcsmooth),&
                MAX(1,j-bcsmooth):MIN(xmax,j+bcsmooth),&
                MAX(1,k-bcsmoothz):MIN(SIZE(Cube,DIM=3),k+bcsmoothz))/=UNDEF)) THEN

              Smooth(i,j,k)=Mean(PACK(Cube(MAX(1,i-bcsmooth):MIN(xmax,i+bcsmooth),&
                                              MAX(1,j-bcsmooth):MIN(xmax,j+bcsmooth),&
                                              MAX(1,k-bcsmoothz):MIN(SIZE(Cube,DIM=3),k+bcsmoothz)), &
                                               MASK=&
                                              Cube(MAX(1,i-bcsmooth):MIN(xmax,i+bcsmooth),&
                                                   MAX(1,j-bcsmooth):MIN(xmax,j+bcsmooth),&
                                                   MAX(1,k-bcsmoothz):MIN(SIZE(Cube,DIM=3),k+bcsmoothz))/=UNDEF))

           ELSE

              Smooth(i,j,k)=UNDEF
              
           END IF

        END DO
     END DO
  END DO

  !..copy over
  Cube(:,:,zmin_:zmax_)=Smooth(:,:,zmin_:zmax_)

  !..free memory
  DEALLOCATE(Smooth)

END SUBROUTINE BoxCar

!----------------------------------------------------------------

SUBROUTINE BoxCarVar(bcsmooth, bcsmooth_z, zmin, zmax)

  USE StatLib
  IMPLICIT NONE
  INTEGER(kind=4),INTENT(IN) :: bcsmooth
  INTEGER(kind=4),INTENT(IN), OPTIONAL :: bcsmooth_z
  INTEGER, INTENT(IN), OPTIONAL :: zmin, zmax
  INTEGER(kind=4) :: zmin_, zmax_, i, j, k, bcsmoothz, xmax, ymax
  REAL(kind=4), ALLOCATABLE :: Smooth(:,:,:)


  !..default values for optional parameter
  bcsmoothz=bcsmooth
  zmin_=1
  zmax_=SIZE(Cube,DIM=3)

  !..adjust values if necessary
  IF(present(bcsmooth_z)) bcsmoothz=bcsmooth_z
  IF(present(zmin)) zmin_=zmin
  IF(present(zmax)) zmax_=zmax

  print *, "boxcar smoothing..."

  ALLOCATE(Smooth(SIZE(Var,DIM=1),SIZE(Var,DIM=2),zmin_:zmax_))
  ymax=SIZE(Var,DIM=2)
  xmax=SIZE(Var,DIM=1)

  !print *, xmax, ymax, zmin_, zmax_, bcsmooth, bcsmoothz

  DO k=zmin_,zmax_
     DO j=1,ymax
        DO i=1,xmax

           IF(ANY(Var(MAX(1,i-bcsmooth):MIN(xmax,i+bcsmooth),&
                MAX(1,j-bcsmooth):MIN(xmax,j+bcsmooth),&
                MAX(1,k-bcsmoothz):MIN(SIZE(Cube,DIM=3),k+bcsmoothz))/=UNDEF)) THEN

              Smooth(i,j,k)=Mean(PACK(Var(MAX(1,i-bcsmooth):MIN(xmax,i+bcsmooth),&
                                              MAX(1,j-bcsmooth):MIN(xmax,j+bcsmooth),&
                                              MAX(1,k-bcsmoothz):MIN(SIZE(Cube,DIM=3),k+bcsmoothz)), &
                                               MASK=&
                                              Var(MAX(1,i-bcsmooth):MIN(xmax,i+bcsmooth),&
                                                   MAX(1,j-bcsmooth):MIN(xmax,j+bcsmooth),&
                                                   MAX(1,k-bcsmoothz):MIN(SIZE(Cube,DIM=3),k+bcsmoothz))/=UNDEF))

              Smooth(i,j,k)=Smooth(i,j,k)/COUNT(Var(MAX(1,i-bcsmooth):MIN(xmax,i+bcsmooth),&
                                                   MAX(1,j-bcsmooth):MIN(xmax,j+bcsmooth),&
                                                   MAX(1,k-bcsmoothz):MIN(SIZE(Cube,DIM=3),k+bcsmoothz))/=UNDEF)

           ELSE

              Smooth(i,j,k)=UNDEF
              
           END IF

        END DO
     END DO
  END DO

  !..copy over
  Var(:,:,zmin_:zmax_)=Smooth(:,:,zmin_:zmax_)

  !..free memory
  DEALLOCATE(Smooth)

END SUBROUTINE BoxCarVar


!------------------------------------------------------------------

RECURSIVE SUBROUTINE GetCubeSize(InpFile, CubeDims)

  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN):: InpFile
  INTEGER, INTENT(OUT) :: CubeDims(3)
  INTEGER(kind=4) :: status, unit, rwstatus, blocksize, naxes(3), group, in, end, ext, nfound, i, rank, ival
  CHARACTER(len=350) :: fname, comment
  LOGICAL :: anyf
  REAL(kind=8) :: dum

  status=0

  !..get an unused unit
  CALL ftgiou(unit,status)

  IF(status/=0) STOP "problem with ftgiou"

  !..open first extension: DATA
  fname=TRIM(InpFile)//"[1]"
  rwstatus=0   ! 0=readonly
  CALL ftdopn(unit,fname,rwstatus,status)
  IF(status/=0) THEN
     !..try to open a single fits file with no extension
     fname=TRIM(InpFile)
     status=0
     CALL ftgiou(unit,status)
     CALL ftdopn(unit,fname,rwstatus,status)
     IF(status/=0) THEN
        print *, "problem reading file: ", TRIM(fname)
        STOP 
     END IF
  END IF

  IF(Verbosity>=2) print *, "Reading:",TRIM(fname)

  !..get cube rank
  CALL ftgkyj(unit,'NAXIS',rank,comment,status)

  IF(rank==3) THEN
     !..get cube size
     CALL ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
  ELSEIF(rank==2) THEN
     !print *, "reading as image..."
     !..get image size
     CALL ftgkyj(unit,'NAXIS1',ival,comment,status)
     naxes(1)=ival
     CALL ftgkyj(unit,'NAXIS2',ival,comment,status)
     naxes(2)=ival
     naxes(3)=1
     !print *, "naxes=",naxes
  ELSE
     print *, "input file is not a cube or a image according to NAXIS keyword:", rank
     STOP
  END IF

  !..read and store WCS information, if present 

  !..check that keywords are there...
  status=0
  CALL ftgkyd(unit,"CRPIX1",dum,comment,status)
  
  IF(status==0) THEN !read keywords
     DO i=1,SIZE(WCS)
        CALL ftgkyd(unit,TRIM(WCSlabels(i)),WCS(i),comment,status)
        IF(status/=0) THEN
           IF(VERBOSITY>=1) print *, "WARNING:: WCS keyword ", TRIM(WCSlabels(i)), " not found!"
           status=0 
           WCS=0
           WCS(11)=1. !.. that's the scale factor for z-direction
        END IF
     END DO
     DO i=1,SIZE(WCSstrings)
        CALL ftgkys(unit,TRIM(WCSlabels_strings(i)),WCSstrings(i),comment,status)
        IF(status/=0) THEN
           IF(VERBOSITY>=1) print *, "WARNING:: WCS keyword ", TRIM(WCSlabels_strings(i)), " not found!"
           WCSstrings(:)="??"
           status=0 
        END IF
     END DO
     
  ELSE

     IF(VERBOSITY>=2) THEN
        print *, " "
        print *, "NB: no WCS info found in input file"
        print *, " "
     END IF
     WCS(:)=0.d0
     status=0
  END IF


  IF(Verbosity>=2) print *, "Cube dimensions =", naxes
       
  CubeDims=naxes
        
  CALL ftclos(unit,status)
  CALL ftfiou(unit,status)


END SUBROUTINE GetCubeSize

!-----------------------------------------------

SUBROUTINE GetNExt(Inpfile, n_ext)

  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN):: InpFile
  CHARACTER(len=500) :: fname
  INTEGER(kind=4) :: n_ext
  INTEGER(kind=4) :: status, unit, rwstatus

  status=0

  !..get an unused unit
  CALL ftgiou(unit, status)
  IF(status/=0) STOP "problem with ftgiou"

  !..try to open second extension
  rwstatus=0   ! 0=readonly
  fname=TRIM(InpFile)//"[2]"
  CALL ftdopn(unit,fname,rwstatus,status)
  IF(status/=0) THEN
     n_ext=1
  ELSE
     n_ext=2
  END IF

  status=0
  CALL ftclos(unit,status)
  CALL ftfiou(unit,status)
  
END SUBROUTINE GetNExt

  
!------------------------------------

RECURSIVE SUBROUTINE ReadLocalCube(InpFile, LocalCube, unit)

  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN):: InpFile
  REAL(kind=4), INTENT(OUT) :: LocalCube(:,:,:)
  INTEGER(kind=4), OPTIONAL, INTENT(IN) :: unit
  INTEGER(kind=4) :: status, unit_, rwstatus, blocksize, naxes(3), group, in, end, ext, nfound, n_ext_, i, rank, ival
  CHARACTER(len=350) :: comment
  LOGICAL :: anyf


  status=0

  !..get an unused unit if not present
  IF(present(unit)) THEN
     unit_=unit
  ELSE
     CALL ftgiou(unit_,status)
  END IF

  IF(status/=0) STOP "problem with ftgiou"

  rwstatus=0   ! 0=readonly
  CALL ftdopn(unit_,InpFile,rwstatus,status)
  IF(status/=0) THEN
     print *, "problem reading file: ", TRIM(InpFile)
     STOP 
  END IF

  IF(Verbosity>=2) print *, "Reading:",TRIM(InpFile)

  !..get original cube rank
  CALL ftgkyj(unit_,'NAXIS',rank,comment,status)
  
  IF(rank==3) THEN
     !..get original cube size
     CALL ftgknj(unit_,'NAXIS',1,3,naxes,nfound,status)
  ELSEIF(rank==2) THEN
     !print *, "reading as image..."
     !..get image size
     CALL ftgkyj(unit_,'NAXIS1',ival,comment,status)
     naxes(1)=ival
     CALL ftgkyj(unit_,'NAXIS2',ival,comment,status)
     naxes(2)=ival
     naxes(3)=1
     !print *, "naxes=",naxes
  ELSE
     print *, "input file is not a cube or a image according to NAXIS keyword:", rank
     STOP
  END IF

  IF(Verbosity>=2) THEN
     print *, "Cube dimensions =", naxes
     print *, "LocalCube dimensions=", SHAPE(LocalCube)
  END IF
  
  in=1
  end=naxes(1)*naxes(2)*naxes(3)
        
  group=1

  !ALLOCATE(DATA(naxes(1),naxes(2),naxes(3)))

  !..read local cube
  CALL ftgpve(unit_,group,in,end,UNDEF,LocalCube,anyf,status)
  IF(status/=0) STOP "Problem reading data array to local cube"
  
  CALL ftclos(unit_,status)

END SUBROUTINE ReadLocalCube

!-----------------------------------------------

END MODULE CubeLib

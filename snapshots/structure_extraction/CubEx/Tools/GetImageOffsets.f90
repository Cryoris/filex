
PROGRAM main

 IMPLICIT NONE
 CHARACTER(len=500) :: string, catlist(2), fname, catchip2, Offsets, OffsetsAll, WCSCube, WCSOut, thisWCS, comment
 INTEGER :: ierr, nfiles, maxsources, xrange(2), yrange(2), i, ii, idum, init, dx, dy, maxoffset, max_chipoffset_y, dx_(1), dy_(1), lun_chip, lun_file, nchips, n, chip, min_edge, max_edge, smoothr, nmin, nmax, pad
 REAL(kind=4), ALLOCATABLE, DIMENSION(:,:) :: x, y, flux
 REAL(kind=4), ALLOCATABLE, DIMENSION(:) :: dx_arr(:), dy_arr(:), dx_smooth(:), dy_smooth(:), w(:), w_norm(:)
 INTEGER(kind=4), ALLOCATABLE :: ns(:)
 REAL(kind=4)     :: df, this_df, pixsize
 LOGICAL :: checksolution, ex
 INTEGER, PARAMETER :: nn=3          !.. this parameter determines the filtering mask size (=2*ns*smoothr)
 !...WCS parameters:
 CHARACTER(len=200) :: WCSlabels(11)=["CRPIX1","CRPIX2","CRPIX3","CRVAL1","CRVAL2","CRVAL3","CD1_1 ","CD1_2 ","CD2_1 ","CD2_2 ","CD3_3 "]
 CHARACTER(len=200) :: WCSlabels_strings(6)=["CTYPE1","CTYPE2","CTYPE3","CUNIT1","CUNIT2","CUNIT3"]
 REAL(kind=8) :: WCS(11), WCSnew(11), dum
 CHARACTER(len=8) :: WCSstrings(6)
 !..CFITSIO variables:
 INTEGER(kind=4) :: status, unit, rwstatus, blocksize, naxes(3), group, nfound


 CALL ReadParameters

 nmin=-nn*smoothr
 nmax=nn*smoothr
     
 IF(ALLOCATED(w)) DEALLOCATE(w,w_norm)
 ALLOCATE(w(nmin:nmax),w_norm(nmin:nmax))

 IF(smoothr>0) THEN
    !..generate 1D filter  
    DO i=nmin,nmax,1
       w(i)=exp(-(REAL(i)**2/(2.*smoothr**2)))
    END DO
    !..normalized version
    w_norm(:)=w(:)/SUM(w(:))
 ELSE
    w=1.
    w_norm=1.
 END IF

!..allocate histogram array with padding
 ALLOCATE(dx_arr(-maxoffset+nmin:maxoffset+nmax), dy_arr(-maxoffset+nmin:maxoffset+nmax), dx_smooth(-maxoffset+nmin:maxoffset+nmax), dy_smooth(-maxoffset+nmin:maxoffset+nmax))


!... count number of files per chip
 nfiles=0
 ierr=0
 OPEN(91,file=catlist(1),action="read")
 print *, "checking number of files in ", TRIM(catlist(1))
 DO
    READ(91,'(a)',iostat=ierr) string
    IF(ierr/=0) EXIT
    nfiles=nfiles+1
 END DO

 print *, "nfiles=",nfiles

 ALLOCATE(x(nfiles,maxsources),y(nfiles,maxsources),flux(nfiles,maxsources),ns(nfiles))
 x=0 ; y=0 ; flux=0 ; ns=0

 ! ---------- assign values and count sources in catalogues

 !..open all chip files
 REWIND(91)
 IF(nchips>1) THEN
    DO i=2,nchips
       lun_chip=90+i
       OPEN(lun_chip, file=catlist(i), action="read")
    END DO
 END IF

 DO i=1,nfiles

    ii=1

    DO chip=1,nchips

       lun_chip=90+chip
       lun_file=10+i

       READ(lun_chip,*) fname ! read chip1
       OPEN(lun_file,file=fname,action="read")
       DO 
          READ(lun_file,*,iostat=ierr) idum, x(i,ii), y(i,ii), flux(i,ii)
          IF(ierr/=0) THEN
             x(i,ii)=0; y(i,ii)=0; flux(i,ii)=0
             EXIT
          END IF
          IF(x(i,ii)<min_edge.or.y(i,ii)<min_edge.or.x(i,ii)>max_edge.or.y(i,ii)>max_edge) THEN
             ii=ii-1 !..remove this source
          END IF
          IF(ii==maxsources) EXIT
          ii=ii+1
       END DO
       CLOSE(lun_file)
       
    END DO

    ns(i)=COUNT(x(i,:)/=0)
    print *, "file=", i, " ns=",ns(i)

 END DO

 !..close all chip files
 CLOSE(91)
 IF(nchips>1) THEN
    DO i=2,nchips
       lun_chip=90+i
       CLOSE(lun_chip)
    END DO
 END IF


!..for each file compare the sources in the first catalogue with all the rest
!..build histograms of offsets between sources with similar fluxes
!..the offsets are the location of the maxima in the histograms
!..sign of offsets are standard iraf offsets

 OPEN(1,file=Offsets,action="write")

 print *, "file #              xoff     yoff     P(x)     P(y)"

 DO i=1,nfiles

    dx_arr=0
    dy_arr=0

    firstcat:DO init=1,ns(1)
       allothers:DO ii=1,ns(i)

          this_df=abs((flux(i,ii)-flux(1,init))/flux(i,ii))

          IF(this_df<df) THEN

             dx=NINT((x(i,ii)-x(1,init))*10)
             dy=NINT((y(i,ii)-y(1,init))*10)

             IF(abs(dx)>maxoffset.or.abs(dy)>maxoffset) CYCLE
             
             !..gaussian distribute over bins centered on dx & dy 
             dx_arr(dx+nmin:dx+nmax)=dx_arr(dx+nmin:dx+nmax)+w_norm(nmin:nmax)
             dy_arr(dy+nmin:dy+nmax)=dy_arr(dy+nmin:dy+nmax)+w_norm(nmin:nmax)

          END IF

       END DO allothers
    END DO firstcat

    dx_smooth=0
    dy_smooth=0
    DO ii=-maxoffset,maxoffset
       dx_smooth(ii)=SUM(dx_arr(ii+nmin:ii+nmax))!*w_norm(:))
       dy_smooth(ii)=SUM(dy_arr(ii+nmin:ii+nmax))!*w_norm(:))
    END DO
    dx_arr=dx_smooth
    dy_arr=dy_smooth

    dx_=-(MAXLOC(dx_arr)+LBOUND(dx_arr)-1)
    dy_=-(MAXLOC(dy_arr)+LBOUND(dy_arr)-1)


    print *, i, dx_*0.1, dy_*0.1, MAXVAL(dx_arr)/REAL(ns(i)), MAXVAL(dy_arr)/ns(i)

    IF(checksolution) THEN
       print *, "offset  dx_arr  dy_arr"
       DO ii=-maxoffset,maxoffset
          print *, ii, dx_arr(ii), dy_arr(ii)
       END DO
    END IF

    WRITE(1,*)  dx_*0.1*pixsize, dy_*0.1*pixsize

    IF(TRIM(WCSCube)/="??") THEN !..apply correction to WCS headers

       IF(i==1) THEN

          !..read reference frame WCS
          CALL ReadWCS(WCSCube)

          !..add padding if requested
          naxes(1:2)=naxes(1:2)+2*pad

       ELSE

          !..write updated WCS

          WRITE(string,'(i4.4)') i
          thisWCS=TRIM(WCSOut)//TRIM(string)//".fits"

          print *, "writing ", TRIM(thisWCS)

          !..update CRPIX1 and CRPIX2 in WCS
          WCSnew=WCS
          WCSnew(1)=WCS(1)+dx_(1)*0.1+pad
          WCSnew(2)=WCS(2)+dy_(1)*0.1+pad

          CALL WriteWCS(thisWCS)

 
       END IF

    END IF


 END DO
 CLOSE(1)

 IF(nchips>1) CALL FindChipOffset

 print *, "end of Job! :-) "

CONTAINS 

SUBROUTINE ReadParameters


  IMPLICIT NONE
  CHARACTER(len=500) :: ExeName, fname, string, ds, opt, arg 
  INTEGER :: narg, iarg, i, is
  LOGICAL :: ex

  nchips=1
  catlist(1)="??"
  catlist(2)="??"
  Offsets="??"
  OffsetsAll="??"
  maxsources=100   ! maximum number of sources to consider in sextractor catalogue 
  df=0.5         ! difference in relative flux for matching sources
  maxoffset=200  ! pixels, this may be an arbitrary large number (limited by computational speed)
  max_chipoffset_y=3  ! pixels
  min_edge=0
  max_edge=10000000
  smoothr=0.5
  checksolution=.false.
  pixsize=5.555556e-5  !.. units: deg for MUSE
  WCSCube="??"
  WCSOut="output_wcs_"
  pad=0

  CALL GetArg(0,ExeName)
  narg=iargc()
  IF(narg>0) THEN
     CALL GetArg(1,fname)
  ELSE
     WRITE(*,'(a)') " "
     WRITE(*,'(a)') " options: "
     WRITE(*,'(a)') " -nchips  <int>             : number of chips per image frame (default=1)"
     WRITE(*,'(a)') " -catlist <name(s)>         : filename (one per chip) containing the list of catalogue names (NO default)"
     WRITE(*,'(a)') " -out     <name>            : name of output offset file per chip (NO default)"
     WRITE(*,'(a)') " -outall  <name>            : name of output offset file combining chips, only active if nchips>1 (NO default)"
     WRITE(*,'(a)') " -maxsources <int>          : maximum number of sources to consider in each catalogue (default=100)"
     WRITE(*,'(a)') " -df <float>                : difference in relative flux for matching sources (default=0.5)"
     WRITE(*,'(a)') " -maxoffset <int>           : maxoffset in pixels (default=200)"
     WRITE(*,'(a)') " -max_chipoffset <int>      : max offset between chips (default=3)"
     WRITE(*,'(a)') " -min_edge  <float>         : exclude sources closer than this to the min image edges (default=0)"
     WRITE(*,'(a)') " -max_edge  <float>         : as above for max edges (default=100000)"
     WRITE(*,'(a)') " -smoothr <int>             : smoothing radius for the solution histogram in pixels (default=minimum=0.5)"
     WRITE(*,'(a)') " -pixsize <float>           : pixel size for output (default=5.5556e-5 deg, MUSE pixel)"
     WRITE(*,'(a)') " -wcscube <name>            : if defined, uses the header in WCSCube to produce a series of header-files,"
     WRITE(*,'(a)') "                              one per catalogue with updated CRPIX[12] according to the offsets"
     WRITE(*,'(a)') "                              NB: reference frame MUST be the first frame in the catalogue!"
     WRITE(*,'(a)') " -wcsout  <name>            : if defined, use this as a base for the updated wcs headers "
     WRITE(*,'(a)') "                              (+'number'.fits), default='output_wcs_<number>.fits"
     WRITE(*,'(a)') " -pad     <int>             : padding size from each edge of the original cube for WCSCube (default=0, no padding)"
     STOP
  END IF

  !..read from command line
  ii=0
  DO WHILE(ii<narg)
     ii=ii+1; CALL getarg(ii,opt)
     ii=ii+1; CALL getarg(ii,arg)
     SELECT CASE(TRIM(opt))
     CASE('-nchips')        
        READ(arg,*) nchips
        IF(nchips>2) STOP "nchips should be less or equal than 2!"
     CASE('-catlist')      
       READ(arg,'(a)') catlist(1)
        DO n=2,nchips
           ii=ii+1
           CALL getarg(ii,arg)
           READ(arg,'(a)') catlist(n)
        END DO  
     CASE('-out')           ; READ(arg,'(a)') Offsets
     CASE('-outall')        ; READ(arg,'(a)') OffsetsAll
     CASE('-maxsources')    ; READ(arg,*) maxsources
     CASE('-df')            ; READ(arg,*) df
     CASE('-maxoffset')     ; READ(arg,*) maxoffset
     CASE('-max_chipoffset'); READ(arg,*) max_chipoffset_y
     CASE('-min_edge')      ; READ(arg,*) min_edge
     CASE('-max_edge')      ; READ(arg,*) max_edge
     CASE('-smoothr')       ; READ(arg,*) smoothr
     CASE('-checksolution') ; READ(arg,*) checksolution
     CASE('-pixsize')       ; READ(arg,*) pixsize
     CASE('-pad')           ; READ(arg,*) pad
     CASE('-wcscube')       ; READ(arg,'(a)') WCSCube
     CASE('-wcsout')        ; READ(arg,'(a)') WCSOut
     CASE default
        print *, "parameter ", TRIM(opt), " not recognized!"
     END SELECT
  END DO

!..perform few checks
  IF(TRIM(catlist(1))=="??") STOP "please provide catalogue list with the option -catlist "
  IF(TRIM(Offsets)=="??") STOP "please provide the output offset file with the option -out "
  IF(nchips>1.and.TRIM(OffsetsAll)=="??") STOP "please provide the offset file combining chips with the option -outall "

!..multiply maxoffset and smoothr by 10 to go to then 1/10nth of pixel size
 maxoffset=maxoffset*10
 smoothr=MAX(5,smoothr*10)


END SUBROUTINE ReadParameters

SUBROUTINE ReadWCS(InpWCS)

  IMPLICIT NONE
  CHARACTER(len=*) :: InpWCS
  INTEGER :: ii, bitpix

  status=0
  !..get an unused unit
  CALL ftgiou(unit,status)
  IF(status/=0) STOP "problem with ftgiou"

  ! 0=readonly
  CALL ftdopn(unit,InpWCS,0,status)

  print *, "reading WCS information from reference frame:", TRIM(InpWCS)

  !..get cube size
  CALL ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)

  !..read and store WCS information

  !..check that keywords are there...
  status=0
  CALL ftgkyd(unit,"CRPIX1",dum,comment,status)
  IF(status/=0) STOP "NO WCS info found in WCSCube!"
  DO ii=1,SIZE(WCS)
     CALL ftgkyd(unit,TRIM(WCSlabels(ii)),WCS(ii),comment,status)
     IF(status/=0) THEN
        print *, "WARNING:: WCS keyword ", TRIM(WCSlabels(ii)), " not found!"
        STOP
     ENDIF
  END DO
  DO ii=1,SIZE(WCSstrings)
     CALL ftgkys(unit,TRIM(WCSlabels_strings(ii)),WCSstrings(ii),comment,status)
     IF(status/=0) THEN
        print *, "WARNING:: WCS keyword ", TRIM(WCSlabels_strings(ii)), " not found!"
        STOP 
     END IF
  END DO

  print *, WCS

  CALL ftclos(unit,status)
  CALL ftfiou(unit,status)

END SUBROUTINE ReadWCS

!----------------------------------------------------------

SUBROUTINE WriteWCS(datafile)

  IMPLICIT NONE
  CHARACTER(len=*) :: datafile
  INTEGER :: ii, naxis, bitpix
  LOGICAL :: simple, extend


  status=0
  !..get an unused unit
  CALL ftgiou(unit,status)
     
  IF(status/=0) STOP "problem with ftgiou"
     
  rwstatus=1  ! 1=readwrite
  INQUIRE(file=datafile,EXIST=ex)
  IF(ex) THEN !..delete file
     CALL ftdopn(unit,datafile,rwstatus,status)
     IF(status/=0) THEN
        print *, "problem with ftdopn! Try to remove the file ",TRIM(datafile), " manually..."
        STOP
     END IF
     CALL ftdelt(unit, status)
     IF(status/=0) STOP "problem with ftdelt"
  ENDIF
  blocksize=1
  CALL ftinit(unit,datafile,blocksize,status)
  IF(status/=0) STOP "problem with ftinit"
          
  !  Initialize parameters about the FITS image.
  simple=.true.
  extend=.true.
  bitpix=-32
  naxis=3
     
  !  Write the required header keywords to the file
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  IF(status/=0) STOP "problem with ftphpr"
  
  !..write updated WCS
  DO ii=1,SIZE(WCS)
     CALL ftpkyd(unit,TRIM(WCSlabels(ii)),WCSnew(ii),12,' ',status)
  END DO
  DO ii=1,SIZE(WCSstrings)
     CALL ftpkys(unit,TRIM(WCSlabels_strings(ii)),WCSstrings(ii),' ',status)
  END DO
  
  CALL ftclos(unit,status)
  CALL ftfiou(unit,status)
  


END SUBROUTINE WriteWCS

!--------------------------------------------------------


SUBROUTINE FindChipOffset

  IMPLICIT NONE

!..find offsets between chips, using same sources landing on different chips (same flux but different sign) after offsets

  OPEN(1,file=Offsets,action="read")
  OPEN(11,file=OffsetsAll,action="write")
  DO i=1,nfiles

     READ(1,*) dx, dy

     !...apply offsets between images for next step, x & y are now registered to the same frame
     x(i,:)=x(i,:)+dx
     y(i,:)=y(i,:)+dy

     !..copy original values for chip1
     WRITE(11,*) dx, dy  

  END DO

  print *, " "

  dx_arr=0
  dy_arr=0

  DO i=2,nfiles

     firstcat2:DO init=1,ns(1)
        allothers2:DO ii=1,ns(i)

           this_df=abs((flux(i,ii)+flux(1,init))/flux(i,ii)) ! this selects the same source on different chip

           IF(this_df<df) THEN

              dx=NINT(x(i,ii)-x(1,init))
              dy=NINT(y(i,ii)-y(1,init))
             
              IF(abs(dy)<max_chipoffset_y) THEN
                 print *, i,x(i,ii),x(1,init),dy,this_df
                 dx_arr(dx)=dx_arr(dx)+1./this_df
                 dy_arr(dy)=dy_arr(dy)+1./this_df
              END IF

 
           END IF

        END DO allothers2
     END DO firstcat2

  END DO

  dx_=(MAXLOC(dx_arr)+LBOUND(dx_arr)-1)
  dy_=(MAXLOC(dy_arr)+LBOUND(dy_arr)-1)


  print *, "chipoffset(x)=", dx_ , "N=",MAXVAL(dx_arr)
  print *, "chipoffset(y)=", dy_ , "N=",MAXVAL(dy_arr) 

  dy_=0
  print *, "NB:chipoffset(y) set to 0"
 
  REWIND(1)
  DO i=1,nfiles
     READ(1,*) dx, dy
     WRITE(11,*) dx+dx_,dy+dy_
  END DO
  CLOSE(1)
  CLOSE(11)

END SUBROUTINE FindChipOffset


END PROGRAM main
 

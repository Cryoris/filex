PROGRAM CubePSFSub

  USE CubeLib
  USE StatLib
  IMPLICIT NONE
  CHARACTER(len=250) :: inpfile, outfile, psfimage, psfim_shape
  INTEGER :: i, n_ext, DimX, DimY, DimZ, xmin, xmax, ymin, ymax,  zPSF_min, zPSF_max, zSUB_min, zSUB_max, &
       x1, x2, y1, y2, j, zbin, nbins, x, y, zPSFsize, zstepsize, &
       half_zstepsize, zcen, half_zPSFsize, percent, old_percent, maskpix(2)
  REAL(kind=4) :: xs, ys, rmin, rmax, normfact, rescalefact, new_xy(2), this_mean, this_median, sigma, noisethr
  REAL(kind=4), ALLOCATABLE :: tmp(:,:)
  LOGICAL :: recenter, withvar, rmPSFnoise
  LOGICAL, ALLOCATABLE :: zmask(:)


  CALL ReadCommandLine

  !..read input cube
  Verbosity=2
  IF(withvar) THEN
     n_ext=2
  ELSE
     n_ext=1
  END IF
  CALL ReadCube(Inpfile=inpfile, n_ext=n_ext)
  DimX=SIZE(Cube,DIM=1)
  DimY=SIZE(Cube,DIM=2)
  DimZ=SIZE(Cube,DIM=3)
  ALLOCATE(Image(DimX,DimY))


  ALLOCATE(zmask(DimZ))
  zmask=.true.
  IF(ALL(maskpix/=0)) THEN
     print *, "masking pixels in range: ", maskpix(1), maskpix(2)
     zmask(maskpix(1):maskpix(2))=.false.
  END IF

  !..set default for nbins=-1 and max number of bins
  IF(nbins==-1) THEN
     nbins=DimZ
  ELSE
     nbins=MIN(nbins,DimZ)
  END IF

  !..check if zPSFsize was set
  IF(nbins==DimZ.and.zPSFsize<=1) STOP "You must set zPSFsize>1 given the large requested nbins!"
     
  !..get zstep size
  zstepsize=DimZ/nbins
  half_zstepsize=zstepsize/2
  IF(zPSFsize>1) THEN
     half_zPSFsize=zPSFsize/2
  ELSE
     half_zPSFsize=0
  END IF

  old_percent=0
  !print *, "Removing PSF "
  DO zbin=1,nbins

     IF(nbins==1) THEN
        !..assign missing default values if necessary and perform a boundary check
        IF(zPSF_max==-1.0) zPSF_max=DimZ
        IF(zSUB_max==-1.0) zSUB_max=DimZ
        zPSF_max=MIN(zPSF_max,DimZ)
        zSUB_max=MIN(zSUB_max,DimZ)
     ELSE
        IF(zPSFsize==-1) THEN
           zPSF_min=(zbin-1)*zstepsize+1
           zPSF_max=zbin*zstepsize
           zSUB_min=zPSF_min
           zSUB_max=zPSF_max
           print *, "zbin=",zbin,"zPSF_min / max=", zPSF_min, zPSF_max
        ELSE
           zcen=zbin*zstepsize-half_zstepsize
           zPSF_min=MAX(1,zcen-half_zPSFsize)
           zPSF_max=MIN(DimZ,zcen+half_zPSFsize)
           zSUB_min=zcen-half_zstepsize
           zSUB_max=zcen+half_zstepsize
           !..write alive signal
           percent=INT(REAL(zbin)/nbins*20)
           IF(percent/=old_percent) WRITE(*,'(i4,1a,$)') percent*5,"%" 
           old_percent=percent
        END IF
        !print *, "zbin=",zbin,"zPSF_min / max=", zPSF_min, zPSF_max
     END IF

     !..perform a check
     IF(ALL(zmask(zPSF_min:zPSF_max).eqv..false.)) STOP "too many masked z-pixels given the selected nbins, and PSF parameters!"

     !..get a first guess on PSF image location (+/-3 pixels) to save computational time
     !..if rmPSFNoise is requested, all the available FOV is used because is needed for the statistic
     IF(rmPSFnoise) THEN
        xmin=1; xmax=DimX; ymin=1; ymax=DimY
     ELSE
        xmin=MAX(1,INT(xs-rmax-3))
        xmax=MIN(DimX,INT(xs+rmax+3))
        ymin=MAX(1,INT(ys-rmax-3))
        ymax=MIN(DimY,INT(ys+rmax+3))
     END IF
        
     !print *, "xmin, xmax, ymin, ymax, DimX, DimY=", xmin,xmax,ymin,ymax, DimX, DimY
     
     !..produce a white-light image for PSF extraction 
     Image=UNDEF
     DO j=ymin,ymax
        DO i=xmin,xmax
           IF(ANY(Cube(i,j,zPSF_min:zPSF_max)/=UNDEF.and.zmask(zPSF_min:zPSF_max))) THEN
              !CALL SigmaClip(PACK(Cube(i,j,zPSF_min:zPSF_max),MASK=Cube(i,j,zPSF_min:zPSF_max)/=UNDEF),this_mean,this_median,sigma)
              this_mean=Mean(PACK(Cube(i,j,zPSF_min:zPSF_max),MASK=Cube(i,j,zPSF_min:zPSF_max)/=UNDEF.and.zmask(zPSF_min:zPSF_max)))
              Image(i,j)=this_mean
           ELSE
              Image(i,j)=UNDEF
           END IF
        END DO
     END DO
     !..mask negative regions in the image, if requested
     IF(rmPSFnoise) THEN
        IF(ANY(Image/=UNDEF)) THEN
           CALL SigmaClip(PACK(Image,MASK=Image/=UNDEF),this_mean,this_median,sigma,VERBOSITY=2)
           !print *, "image mean, median, sigma=",this_mean, this_median, sigma
           WHERE(Image<noisethr*sigma+this_mean.and.Image/=UNDEF) Image=0.
        END IF
     END IF


     !..rescale image to the selected central values
     !print *, "rescaling image..."
     xmin=MAX(1,INT(xs-rmax))
     xmax=MIN(DimX,INT(xs+rmax))
     ymin=MAX(1,INT(ys-rmax))
     ymax=MIN(DimY,INT(ys+rmax))
     IF(recenter) THEN
        new_xy(1:2)=MAXLOC(image(xmin:xmax,ymin:ymax))+[xmin,ymin]-1.0
        xs=new_xy(1)
        ys=new_xy(2)
        xmin=MAX(1,INT(xs-rmax))
        xmax=MIN(DimX,INT(xs+rmax))
        ymin=MAX(1,INT(ys-rmax))
        ymax=MIN(DimY,INT(ys+rmax))
     END IF
     x1=MAX(1,INT(xs-rmin))
     x2=MIN(DimX,INT(xs+rmin))
     y1=MAX(1,INT(ys-rmin))
     y2=MIN(DimY,INT(ys+rmin))

     !print *, "x1,x2,y1,y2", x1, x2, y1, y2
     
     !..cut out a circular region in the PSF image according to rmax
     IF(PSFim_shape=="circular") THEN
        DO j=ymin,ymax
           DO i=xmin,xmax
              IF((i-xs)*(i-xs)+(j-ys)*(j-ys)>rmax*rmax) image(i,j)=0.
           END DO
        END DO
     END IF
        

     !print *, Mean(PACK(im(x1:x2,y1:y2),MASK=im(x1:x2,y1:y2)/=UNDEF)), im(x1,y1)

     !..subtract zero level background, if present
     !CALL SigmaClip(PACK(Image,MASK=Image/=UNDEF),this_mean, this_median, sigma)
     !print *, "zero level background=", this_mean
     !Image(:,:)=Image(:,:)-this_mean

     !WHERE(Image/=UNDEF) Image(:,:)=Image(:,:)/Mean(PACK(Image(x1:x2,y1:y2),MASK=Image(x1:x2,y1:y2)/=UNDEF))
     !print *, "maxvalue at source center (in unit of the central area mean)=", MAXVAL(Image(xmin:xmax,ymin:ymax)), "source position=",MAXLOC(Image(xmin:xmax,ymin:ymax))+[xmin,ymin]-1
     !print *, x1, x2, y1, y2

     IF(TRIM(PSFimage)/="??") THEN
        Image(1:xmin,:)=0. ; Image(xmax:,:)=0.
        Image(:,1:ymin)=0. ; Image(:,ymax:)=0.
        CALL WriteImage(PSFimage)
     END IF

     !..allocate intermediate array
     IF(ALLOCATED(tmp)) DEALLOCATE(tmp)
     ALLOCATE(tmp(x1:x2,y1:y2))

     !..subtract PSF from each wavelenght layer
     DO i=zSUB_min, zSUB_max

        !..calculate ratio image
        WHERE(Cube(x1:x2,y1:y2,i)>0..and.Image(x1:x2,y1:y2)>0.) 
           tmp(:,:)=Cube(x1:x2,y1:y2,i)/Image(x1:x2,y1:y2)
        ELSEWHERE
           tmp(:,:)=UNDEF
        END WHERE

        IF(ANY(tmp/=UNDEF)) THEN
           CALL SigmaClip(PACK(tmp,MASK=tmp/=UNDEF),this_mean, this_median, sigma, weight=PACK(Cube(x1:x2,y1:y2,i),MASK=tmp/=UNDEF))
           normfact=this_mean
           !normfact=Mean(PACK(Cube(x1:x2,y1:y2,i),MASK=Cube(x1:x2,y1:y2,i)/=UNDEF))/Mean(PACK(Image(x1:x2,y1:y2),MASK=Image(x1:x2,y1:y2)/=UNDEF))
           !DO y=y1,y2
           !   DO x=x1,x2
           !      IF(Cube(x,y,i)/=UNDEF) THEN
           !         tmp(x,y)=image(x,y)*normfact-Cube(x,y,i)
           !      ELSE
           !         tmp(x,y)=UNDEF
           !      END IF
           !   END DO
           !END DO
           !..PSF subtract     
           WHERE(Cube(xmin:xmax,ymin:ymax,i)/=UNDEF) Cube(xmin:xmax,ymin:ymax,i)=Cube(xmin:xmax,ymin:ymax,i)-image(xmin:xmax,ymin:ymax)*normfact
        ELSE
           print *, "no normalization factor found for this layer (no PSF subtraction for this layer)"
        END IF
     END DO

  END DO

  print *, " "

  !..write output cube
  CALL WriteCube(outfile, multiext=withvar)

  print *, "end of Job! :-) "

CONTAINS

 SUBROUTINE ReadCommandLine

  IMPLICIT NONE
  CHARACTER(len=300) :: ExeName, fname, string, ds, opt, arg 
  INTEGER :: narg, iargc, i, is
  LOGICAL :: ex


!..default
  inpfile="??"
  outfile="??"
  psfimage="??"
  xs=-1.0
  ys=-1.0
  rmin=2.0       ! pixels
  rmax=25.0      ! pixels
  zPSF_min=1.0
  zPSF_max=-1.0
  zSUB_min=1.0
  zSUB_max=-1.0
  recenter=.true.
  nbins=1
  withvar=.true.
  psfim_shape='circular'
  zPSFsize=-1
  maskpix=[0,0]
  rmPSFnoise=.true.
  noisethr=1.0
  
  CALL GetArg(0,ExeName)
  narg=iargc()
  IF(narg>0) THEN
     CALL GetArg(1,fname)
  ELSE
     print *, " "
     WRITE(*,'(2a)')"        CubePSFSub (part of CubEx package) ",TRIM(version)
     WRITE(*,'(a)') " "
     WRITE(*,'(a)') " by Sebastiano Cantalupo (cantalupo@phys.ethz.ch)"
     print *, " "
     WRITE(*,'(a)') "usage example: CubePSFSub -cube InpCube.fits -x 153.1 -y 155.7 [-option <val>]"
     WRITE(*,'(a)') " "
     WRITE(*,'(a)') "  options:"
     WRITE(*,'(a)') "  -cube               <string>          : input cube (NO default)"
     WRITE(*,'(a)') "  -x                  <real>            : x position of source to remove (NO default)"
     WRITE(*,'(a)') "  -y                  <real>            : y position of source to remove (NO default)"
     WRITE(*,'(a)') "  -out                <string>          : name of output cube (default=<inputcube>.PSFSub.fits)"
     WRITE(*,'(a)') "  -zPSF_min           <real>            : minimum z position (pixels) for PSF calculation (default=1)"
     WRITE(*,'(a)') "  -zPSF_max           <real>            : maximum z position (pixels) for PSF calculation (default= cube z-size)"
     WRITE(*,'(a)') "  -nbins              <int>             : if>1, performs PSF-sub in equal nbins in z-dir instead of using zPSF_min and zPSF_max (default=1)"
     WRITE(*,'(a)') "                                          NB: set this to a large number or -1 to perform optimized PSF image creation/sub layer-by-layer"
     WRITE(*,'(a)') "                                              in this case, you MUST set zPSFsize as well."
     WRITE(*,'(a)') "  -zPSFsize           <int>             : if nbins>1, use this size (pixels) in z-direction for the PSF image instead of using equal bins"
     WRITE(*,'(a)') "  -zSUB_min           <real>            : minimum z position (pixels) for PSF removal (default=1)"
     WRITE(*,'(a)') "  -zSUB_max           <real>            : maximum z position (pixels) for PSF removal (default= cube z-size)"
     WRITE(*,'(a)') "  -withvar            <bol>             : if .true. (default) reads and propagate the variance extension to the output cube"
     WRITE(*,'(a)') "  -rmin               <int>             : radius in pixels of the (square) region used for PSF image normalization (default=2)"
     WRITE(*,'(a)') "  -rmax               <int>             : radius in pixels of the (square/circular) region used for PSF image subtraction (default=25)"
     WRITE(*,'(a)') "  -psfimshape         <string>          : shape of the PSF image to be removed: 'circular'/'square' (default='circular') "
     WRITE(*,'(a)') "  -psfimout           <string>          : if provided, writes the PSF image in a file with this name (only available for nbins=1)"
     WRITE(*,'(a)') "  -recenter           <bol>             : if .true. (default), recenters the source location to image max around the given x and y pos"
     WRITE(*,'(a)') "  -maskpix            <2 int values>    : initial and final pixel in z direction for line masking (default=0, no masking)"
     WRITE(*,'(a)') "  -rmPSFnoise         <bol>             : if .true. (default), sets to zero all pixels below a given sigma threshold (set by noisethr below)"
     WRITE(*,'(a)') "                                          in the PSF image. This parameter avoids producing artificial flux during PSF subtraction"
     WRITE(*,'(a)') "  -noisethr           <real>            : value in sigmas where rmPSFNoise is applied (see -rmPSFnoise) (default=1)"
     STOP
  END IF

  !..read from command line
  DO i=1,narg,2
     CALL getarg(i,opt)
     CALL getarg(i+1,arg)
     SELECT CASE(TRIM(opt))
     CASE('-cube')          ; READ(arg,'(a)') inpfile
     CASE('-out')           ; READ(arg,'(a)') outfile
     CASE('-x')             ; READ(arg,*) xs
     CASE('-y')             ; READ(arg,*) ys
     CASE('-zPSF_min')      ; READ(arg,*) zPSF_min
     CASE('-zPSF_max')      ; READ(arg,*) zPSF_max
     CASE('-zPSFsize')      ; READ(arg,*) zPSFsize
     CASE('-zSUB_min')      ; READ(arg,*) zSUB_min
     CASE('-zSUB_max')      ; READ(arg,*) zSUB_max
     CASE('-rmin')          ; READ(arg,*) rmin
     CASE('-rmax')          ; READ(arg,*) rmax
     CASE('-nbins')         ; READ(arg,*) nbins
     CASE('-psfimout')      ; READ(arg,'(a)') psfimage
     CASE('-psfimshape')    ; READ(arg,'(a)') psfim_shape
     CASE('-withvar')       ; READ(arg,*) withvar
     CASE('-recenter')      ; READ(arg,*) recenter
     CASE('-maskpix')       ; READ(arg,*) maskpix(1:2)
     CASE('-rmPSFnoise')    ; READ(arg,*) rmPSFnoise
     CASE('-noisethr')      ; READ(arg,*) noisethr
     CASE default
        print *, "command line argument ",TRIM(opt), " not recognized!"
        STOP
     END SELECT
  END DO

!..perform few checks
  IF(TRIM(inpfile)=="??") STOP "Please provide input file name with the option -cube "
  IF(xs==-1.0.or.ys==-1.0)  STOP "Please provide source position with the option -x and -y "
  IF(TRIM(outfile)=="??") THEN !..assign default name
     i=INDEX(TRIM(inpfile),".fits")
     outfile=inpfile(1:i-1)//".PSFSub.fits"
  END IF

  IF(nbins>1.and.TRIM(psfimage)/="??") THEN
     print *, "NB: requested PSFimage file will not be produced, use nbins=1 and select zPSF_min/max zSUB_min/max if you want to check the PSFimage!"
     STOP
  END IF

  IF(TRIM(psfim_shape)/="circular".and.TRIM(psfim_shape)/="square") THEN
     print *, "Selected psfimshape value: ",TRIM(psfim_shape), " not recognized!"
     STOP
  END IF


END SUBROUTINE ReadCommandLine

END PROGRAM CubePSFSub

!----------------------------------------------------------
!
! CubeFix - Flat-field Instrument-space Correction Software
! 
! Author: SC
! Last update: Nov 29 2017
!
PROGRAM CubeFix

 USE StatLib
 USE CubeLib
 IMPLICIT NONE
 CHARACTER(len=500) :: fname, pixtable, comment, datacube, arg, outputbase, CorrCubeName, EDGEMaskName, MaskFile, MapName
 CHARACTER(len=8) :: extname, string, extname_(7)=["xpos    ", "ypos    ", "lambda  ","data    ","dq      ","stat    ","origin  "]
 INTEGER :: MUSE_ORIGIN_SHIFT_IFU=6, MUSE_ORIGIN_SHIFT_XSLICE=24, j, UNDEF_INT=-999, nbins=300
 INTEGER(kind=4) :: unit, status, rwstatus, nfound, naxes(2), in=1, end, group, ii, thispix, n_ext, DimX, DimY, i, maxiter, z, DimZ, XSliceMap_MIN, XSliceMap_MAX, &
      xbin, ybin, blocksize, j_slice, xbin_, ybin_, x, y, MaxIterations, niter, min_xslicemap, max_xslicemap, MinPix, min_vertical_cfact, wl_size, line_wl_size, &
      cube_zmin, cube_zmax, spsize, skybin, binsize(1), this_bin, min_xslice, max_xslice, jj, os_fact, os_min, os_max, i_ifu, e1, e2, k, i_slice, i_xslice, step, &
      EdgeMask_Pix, FOV_EdgeNPix, step_to_run, step_in, step_end, MinStatPixels, FOV_EdgeMask_Pix
 INTEGER(kind=4), ALLOCATABLE :: origin(:,:), IFUmap(:,:), dq(:,:), SliceMap(:,:), slice(:,:), ifu(:,:), xslice(:,:), EDGEMask(:,:), IFUSliceMap(:,:)
 REAL(kind=4), ALLOCATABLE :: data(:,:), np(:,:), im_sigma(:,:), corr_image(:,:), XSliceMap(:,:), pos(:,:), lambda(:,:), sky(:), sky_or(:), sky_sigma(:),&
       IFUmap_(:,:,:), SliceMap_(:,:,:), XSliceMap_(:,:,:), full_sourcemask(:,:), CorrCube(:,:,:), image_save(:,:)
 REAL(kind=8), ALLOCATABLE :: xpos(:,:), ypos(:,:), dcos_ypos(:,:), dsin_ypos(:,:), dcos_xpos(:,:)
 REAL :: this_clipmean, this_clipmedian, this_sigma, this_unclipmean, total_clipmean, total_clipmedian, total_sigma, this_zmin, this_zmax, zmin, zmax, &
      xbinsize, ybinsize, aDEC, cfact(24,48), xfrac, yfrac, pixfrac, this_cfact, xoff, yoff, aRA, source_clip_sigma, CR_clip_sigma, zstep, skystep, avg_sky, &
      skybin_min, skybin_max, this_sky, skysel, minEW, this_corr, this_corr_var, CRclip
 LOGICAL :: anyf, savechecks, skysub, multiext, writeNaN, IFUonly, writeCorrCube, EDGEcorr, native_spherical, writemap, writeoutcube
 REAL(kind=8) ::  xmin, xmax, ymin, ymax, dp, xcen, ycen, xd, yd
 REAL(kind=8), PARAMETER :: pi=3.1415926535d0
 REAL(kind=8) :: deg_to_rad=2.d0*pi/360.d0, rad_to_deg=360.d0/2.d0/pi
 TYPE wlarray
    REAL(kind=4)    :: lmin, lmax, binsize, EW
    INTEGER(kind=4) :: zmin, zmax
    LOGICAL         :: skip
 END type wlarray
 TYPE(wlarray) :: wl(0:1000), wl_old(0:1000), line_wl(0:1000)

 anyf=.false.
 status=0
 spsize=4

 CALL ReadCommandLine

 !..read data and WCS info from datacube
 VERBOSITY=2
 CALL ReadCube(datacube,n_ext=2)
 DimX=SIZE(Cube, DIM=1)
 DimY=SIZE(Cube, DIM=2)
 DimZ=SIZE(Cube, DIM=3)
 IF(writeCorrCube) THEN
    ALLOCATE(CorrCube(DimX,DimY,DimZ))
    CorrCube=1.
 END IF
 VERBOSITY=1
 print *, "MinMax Cube=", MINVAL(Cube,MASK=Cube/=UNDEF), MAXVAL(Cube,MASK=Cube/=UNDEF)
 print *, "MinMax Var=", MINVAL(Var,MASK=Var/=UNDEF), MAXVAL(Var,MASK=Var/=UNDEF)

 !..divide the cube in z-direction
 CALL SplitCube

 !..read PIXTABLE 
 CALL ReadPixTable

 !..allocate map arrays
 ALLOCATE(image(DimX,DimY),IFUmap(DimX,DimY),SliceMap(DimX,DimY),im_sigma(DimX,DimY),&
      corr_image(DimX,DimY),XSliceMap(DimX,DimY), image_save(DimX,DimY))
 ALLOCATE(IFUMap_(DimX,DimY,24),SliceMap_(DimX,DimY,48),XSliceMap_(DimX,DimY,min_xslice:max_xslice))

 ! --- perform 3-step corrections or individual corrections, if requested ---

 IF(step_to_run==0) THEN
    step_in=1
    step_end=3
 ELSE
    step_in=step_to_run
    step_end=step_to_run
 END IF

 step_loop:DO step=step_in,step_end

    print *, " "
    print *, "--- Correction STEP=", step
    print *, " "

    SELECT CASE(step)
    CASE(1) !..slice by slice correction using NBs
       IFUonly=.false.
       EDGEcorr=.false.
    CASE(2) !..ifu correction using sky-lines 
       wl=line_wl
       wl_size=line_wl_size
       IFUonly=.true.
       EDGEcorr=.false.
    CASE(3) !..vertical correction per group and position along the slice using total white-light image.
            !..and, if requested, produces also the edge mask.
       wl_size=1
       wl(1)%zmin=MIN(11,DimZ)
       wl(1)%zmax=MAX(DimZ-10,1)
       wl(1)%lmin=wl(1)%zmin*WCS(11)+WCS(6)
       wl(1)%lmax=wl(1)%zmax*WCS(11)+WCS(6)
       IFUonly=.false.
       EDGEcorr=.true.
    END SELECT
       

    !-- Wavelenght bin Loop ---

    wl_loop:DO z=1,wl_size

       !..find minimum and max wavelenght for this bin
       this_zmin=wl(z)%lmin
       this_zmax=wl(z)%lmax

       print *, " "
       print *, "producing image for zbin=", z , " over ", wl_size

       !..reset map arrays
       image=0.
       IFUmap=0
       SliceMap=0
       XSliceMap=0
       corr_image=1.0
       IFUmap_=0.
       SliceMap_=0.
       XSliceMap_=0.  

       !..produce temporary IFU, Slice and XSlice 3D maps using ra (xpos), dec (ypos) and lambda from pixeltable
       xmin=0
       DO j=1,SIZE(xpos)

          !..remove bad pixels
          IF(dq(1,j)/=0) CYCLE

          !..select pixels in the right wavelength range
          IF(lambda(1,j)<this_zmin.or.lambda(1,j)>this_zmax) CYCLE

          !.."bottom-left pixel"
          xbin=INT(xpos(1,j)-xmin)   !..int part
          xfrac=(xpos(1,j)-xmin)-xbin  !..real part
          ybin=INT(ypos(1,j)-ymin)
          yfrac=(ypos(1,j)-ymin)-ybin
       
          !..loop over the 4 destination pixels
          DO x=-1,1,2
             DO y=-1,1,2
                
                pixfrac=(-0.5*(x-1)+x*xfrac)*(-0.5*(y-1)+y*yfrac)

                xbin_=xbin+0.5*(x+1)
                ybin_=ybin+0.5*(y+1)
                IF(xbin_>DimX.or.xbin_<1) CYCLE
                IF(ybin_>DimY.or.ybin_<1) CYCLE

                IF(pixfrac>=0.25) THEN !..don't consider pixels that contribute less than 25%
                   IFUmap_(xbin_,ybin_,IFU(1,j))=IFUmap_(xbin_,ybin_,IFU(1,j))+pixfrac
                   SliceMap_(xbin_,ybin_,slice(1,j))=SliceMap_(xbin_,ybin_,slice(1,j))+pixfrac
                   XSliceMap_(xbin_,ybin_,xslice(1,j))=XSliceMap_(xbin_,ybin_,xslice(1,j))+pixfrac
                END IF

             END DO
          END DO
       END DO

       !..produce image from datacube with oversampling, if requested
       DO j=1,DimY
          DO i=1,DimX
             IF(os_fact>1) THEN
                ii=(i-1)/os_fact+1
                jj=(j-1)/os_fact+1
             ELSE
                ii=i
                jj=j
             END IF
             IF(ANY(Cube(ii,jj,wl(z)%zmin:wl(z)%zmax)/=UNDEF)) THEN
                image(i,j)=Mean(PACK(Cube(ii,jj,wl(z)%zmin:wl(z)%zmax),MASK=Cube(ii,jj,wl(z)%zmin:wl(z)%zmax)/=UNDEF))
             ELSE
                image(i,j)=UNDEF
             END IF
          END DO
       END DO

       print *, "minmax image=", MINVAL(image,MASK=(image/=UNDEF)), MAXVAL(image)

       IF(savechecks) THEN
          WRITE(string, '(i4.4)') z+1000*step
          fname=TRIM(outputbase)//"_CheckImage_ORIG.z"//TRIM(string)//".fits"
          CALL WriteImage(fname)
       END IF

       CALL SigmaClip(PACK(image,MASK=image/=UNDEF.and.full_sourcemask/=UNDEF), total_clipmean, total_clipmedian, total_sigma, MaxIterations=10)
       print *, "image clipmean, sigma:", total_clipmean, total_sigma

       !..reconstruct IFU, Slice and XSlice 2D maps
       print *, "producing IFU, Slice and XSlice maps.."
       DO j=1,DimY
          DO i=1,DimX
             IF(ANY(IFUMap_(i,j,:)/=0.)) IFUMap(i,j)=MAXLOC(IFUMap_(i,j,:),DIM=1)
             IF(ANY(SliceMap_(i,j,:)/=0.)) SliceMap(i,j)=MAXLOC(SliceMap_(i,j,:),DIM=1)
             IF(ANY(XSliceMap_(i,j,:)/=0.)) XSliceMap(i,j)=MAXLOC(XSliceMap_(i,j,:),DIM=1)+min_xslice-1
          END DO
       END DO
       WHERE(image==UNDEF) 
          IFUMap=0
          SliceMap=0
          XSliceMap=0
       END WHERE

       !..mask image pixels associated with bright sources
       !WHERE(image>source_clip_sigma*total_sigma+total_clipmean) IFUmap=-IFUmap
       WHERE(full_sourcemask==UNDEF.and.IFUmap/=UNDEF) IFUmap=-IFUmap

       !..save XSliceMap, IFUMap and SliceMap images if requested
       IF(savechecks) THEN
          image_save=image !..save image to a temporary array

          image=XSliceMap
          fname=TRIM(outputbase)//"_CheckImage_XSliceMap.z"//TRIM(string)//".fits"
          CALL WriteImage(fname)

          image=IFUMap
          fname=TRIM(outputbase)//"_CheckImage_IFUMap.z"//TRIM(string)//".fits"
          CALL WriteImage(fname)      

          image=SliceMap
          fname=TRIM(outputbase)//"_CheckImage_SliceMap.z"//TRIM(string)//".fits"
          CALL WriteImage(fname)      


          image=image_save !..copy back
       END IF
       !STOP

       print *, "calculating correction map..."

       IF(IFUonly) THEN !..apply only average IFU correction (STEP=2)

          DO j=1,24
             IF(COUNT(image/=UNDEF.and.IFUMap==j)>=MinStatPixels) THEN
                CALL SigmaClip(PACK(image,MASK=(image/=UNDEF.and.IFUmap==j)), this_clipmean, this_clipmedian, this_sigma)
                cfact(j,:)=total_clipmean/this_clipmean
                WHERE(IFUmap==j) corr_image=total_clipmean/this_clipmean
             ELSE
                cfact(j,:)=1.
                WHERE(IFUmap==j) corr_image=1.0
             END IF
          END DO

          IF(savechecks) THEN
             WRITE(string,'(i4.4)') z+1000*step
             fname=TRIM(outputbase)//"_CheckImage_CORR.z"//TRIM(string)//".fits"     
             CALL WriteImage(fname)
          END IF

       ELSEIF(EDGEcorr) THEN !..produces edge mask, ifustackmap and perform vertical correction (STEP=3) 

          IF(EDGEMask_Pix>0.or.FOV_EDGEMask_Pix>0) THEN !..produce a mask for slice edges with value=0 for bad pixels

             ALLOCATE(EDGEMask(DimX,DimY))
             EdgeMask=0

             IF(EdgeMask_Pix>0) THEN
                DO j=1,24
                   DO j_slice=1,48
                      XSliceMap_MIN=MINVAL(XSliceMap,MASK=abs(IFUMap)==j.and.SliceMap==j_slice)
                      XSliceMap_MAX=MAXVAL(XSliceMap,MASK=abs(IFUMap)==j.and.SliceMap==j_slice)
                      WHERE(abs(IFUMap)==j.and.SliceMap==j_slice.and.(XSliceMap>=XSliceMap_MIN+EDGEmask_PIX.and.XSliceMap<=XSliceMap_MAX-EdgeMask_Pix)) EdgeMask=1
                   END DO
                END DO
             END IF

             !..increase the mask at FOV edges if requested (slices 1-12, large XSliceMap, slices 37-48, small XSliceMap)
             IF(FOV_EdgeMask_Pix>EDGEmask_Pix) THEN
                DO j=1,24
                   DO j_slice=1,12
                      XSliceMap_MIN=MINVAL(XSliceMap,MASK=abs(IFUMap)==j.and.SliceMap==j_slice)
                      XSliceMap_MAX=MAXVAL(XSliceMap,MASK=abs(IFUMap)==j.and.SliceMap==j_slice)                  
                      WHERE(abs(IFUMap)==j.and.SliceMap==j_slice.and.XSliceMap>XSliceMap_MAX-FOV_EdgeMask_Pix) EdgeMask=0
                   END DO
                   DO j_slice=37,48
                      XSliceMap_MIN=MINVAL(XSliceMap,MASK=abs(IFUMap)==j.and.SliceMap==j_slice)
                      XSliceMap_MAX=MAXVAL(XSliceMap,MASK=abs(IFUMap)==j.and.SliceMap==j_slice)                  
                      WHERE(abs(IFUMap)==j.and.SliceMap==j_slice.and.XSliceMap<XSliceMap_MIN+FOV_EdgeMask_Pix) EdgeMask=0
                   END DO
                END DO
             END IF
             
          END IF

          IF(writemap) THEN !..produce a map with ifu and stacks
             ALLOCATE(IFUSliceMap(DimX,DimY))
             WHERE(IFUMap/=UNDEF.and.IFUMap/=0.and.SliceMap/=UNDEF.and.SliceMap/=0) 
                IFUSliceMap=abs(IFUMap)*100+SliceMap
             ELSEWHERE
                IFUSliceMap=0
             END WHERE                
          END IF

          CALL VerticalCorr

       ELSE                  !..slice-by-slice and IFU vertical correction

          DO j=1,24
             DO j_slice=1,48
                IF(count(image/=UNDEF.and.IFUmap==j.and.SliceMap==j_slice)>MinStatPixels) THEN
                   CALL SigmaClip(PACK(image,MASK=(image/=UNDEF.and.IFUmap==j.and.SliceMap==j_slice)), this_clipmean, this_clipmedian, this_sigma)
                   cfact(j,j_slice)=total_clipmean/this_clipmean
                ELSE
                   cfact(j,j_slice)=UNDEF  !..it will be fixed below
                END IF
             END DO
          END DO
          
          !..check slices for residual source presence (e.g. cosmic rays)
          !..or lack of statistics. In both cases, we use the average IFU correction.
          DO j=1,24
             !..IFU average correction
             CALL SigmaClip(PACK(cfact(j,:),MASK=cfact(j,:)/=UNDEF), total_clipmean, total_clipmedian, total_sigma)
             IF(wl_size==1) print *, "IFU=",j,"mean=",total_clipmean,"sigma=",total_sigma
             WHERE(cfact(j,:)-total_clipmean>CRclip*total_sigma.or.cfact(j,:)==UNDEF) cfact(j,:)=total_clipmean
          END DO
    

          !..apply slice-by-slice correction
          DO j=1,24
             DO j_slice=1,48
                WHERE(abs(IFUmap)==j.and.SliceMap==j_slice) 
                   image=cfact(j,j_slice)*image
                   corr_image=cfact(j,j_slice)
                END WHERE
             END DO
          END DO

          !..apply additional correction for the initial MinPix of each slice for each of the four "sectors" of the IFUs
          !..the correction is only applied if "negative" (to avoid problems with sources and CRs)
          !..and if we have enough pixels (>MinStatPixels)
          print *, "applying IFU 'vertical' correction..."
          CALL IFU_VerticalCorr

          IF(savechecks) THEN
             WRITE(string,'(i4.4)') z+1000*step
             fname=TRIM(outputbase)//"_CheckImage_CORR.z"//TRIM(string)//".fits"     
             CALL WriteImage(fname)
          END IF

       END IF

       print *, "applying correction to the datacube..."
       print *, "minmax of correction factors=", MINVAL(corr_image), MAXVAL(corr_image)

       IF(wl_size>1) THEN
          cube_zmin=wl(z)%zmin
          cube_zmax=wl(z)%zmax
       ELSE
          cube_zmin=1
          cube_zmax=SIZE(Cube,DIM=3)
       END IF

       os_min=-INT(os_fact/2)
       os_max=INT((os_fact-1)/2)

       IF(writeCorrCube) THEN !..calculate min max of XSliceMap
          XSliceMap_MIN=MINVAL(XSliceMap,MASK=XSliceMap>0)
          XSliceMap_MAX=MAXVAL(XSliceMap)
       END IF

       !..apply correction factors
       DO j=1,SIZE(Cube,DIM=2)
          DO i=1,SIZE(Cube,DIM=1)

             ii=i*os_fact
             jj=j*os_fact

             IF(ANY(corr_image(ii+os_min:ii+os_max,jj+os_min:jj+os_max)/=UNDEF)) THEN

                this_corr=Mean(PACK(corr_image(ii+os_min:ii+os_max,jj+os_min:jj+os_max),MASK=corr_image(ii+os_min:ii+os_max,jj+os_min:jj+os_max)/=UNDEF))
                IF(this_corr<1.) THEN !..increase variance also for "negative" correction 
                   this_corr_var=(2.-this_corr)**2
                ELSE
                   this_corr_var=this_corr**2
                END IF
                IF(this_corr_var<0.) print *, "this_corr_var=",this_corr_var

                WHERE(Cube(i,j,cube_zmin:cube_zmax)/=UNDEF) 
                   Cube(i,j,cube_zmin:cube_zmax)=Cube(i,j,cube_zmin:cube_zmax)*this_corr                       
                END WHERE
                 
                WHERE(Var(i,j,cube_zmin:cube_zmax)/=UNDEF)
                   Var(i,j,cube_zmin:cube_zmax)=Var(i,j,cube_zmin:cube_zmax)*this_corr_var
                END WHERE

                IF(writeCorrCube) THEN
                   IF(XSliceMap(i,j)<XSLiceMap_MIN+10.or.XSliceMap(i,j)>XSliceMap_MAX-10) THEN
                      WHERE(Cube(i,j,cube_zmin:cube_zmax)/=UNDEF)
                         CorrCube(i,j,cube_zmin:cube_zmax)=CorrCube(i,j,cube_zmin:cube_zmax)*this_corr_var
                      END WHERE
                   END IF
                END IF

             END IF

          END DO
       END DO

    END DO wl_loop

 END DO step_loop

 IF(writeoutcube) THEN
    fname=TRIM(outputbase)//".fits"
    !..write updated cube using original header information
    print *, "writing output cube: ", TRIM(fname)
    CALL UpdateCube(InpFile=datacube,OutFile=fname,writeNaN=writeNaN)
 END IF

 IF(EDGEMask_Pix>0.and.ALLOCATED(EDGEMask)) THEN
    image=EDGEMask
    print *, "writing mask file: ", TRIM(EDGEMaskName)
    CALL WriteImage(EDGEMaskName)
 END IF

 IF(ALLOCATED(IFUSliceMap)) THEN
    image=IFUSliceMap
    print *, "writing IFU and Stack Map file: ", TRIM(MapName)
    CALL WriteImage(MapName)
 END IF


 IF(writeCorrCube) THEN
    Cube=CorrCube
    print *, "writing corr cube: ", TRIM(CorrCubeName)
    CALL WriteCube(CorrCubeName,n_ext=1,multiext=.false.)
 END IF
 

CONTAINS 

!-------------------------------------

  SUBROUTINE ReadCommandLine


    IMPLICIT NONE
    CHARACTER(len=500) :: ExeName, fname, string, ds, opt, arg 
    INTEGER :: narg, iarg, i, is
    LOGICAL :: ex

    MinPix=1000
    min_vertical_cfact=1.0
    source_clip_sigma=2.0
    CR_clip_sigma=10.0
    datacube="??"
    pixtable="??"
    outputbase="??"
    savechecks=.false.
    skybin_min=5.
    skybin_max=30.
    os_fact=1
    writeNaN=.true.
    skysel=1.5
    minEW=200.
    writeCorrCube=.false.
    CorrCubeName="??"
    EDGEMask_Pix=3
    FOV_EDGEMask_Pix=7
    FOV_EDGENPix=10
    step_to_run=0
    MaskFile="??"
    writeMap=.true.
    MinStatPixels=8
    CRclip=10.0
    writeoutcube=.true.

    CALL GetArg(0,ExeName)
    narg=iargc()
    IF(narg>0) THEN
       CALL GetArg(1,fname)
    ELSE
       print *, " "
       WRITE(*,'(2a)')"       CubeFix (part of CubEx package) ",TRIM(version)
       WRITE(*,'(a)') "  Flat-fielding Instrument-space Correction Software "
       WRITE(*,'(a)') " "
       WRITE(*,'(a)') " by Sebastiano Cantalupo (cantalupo@phys.ethz.ch)"
       print *, " "
       WRITE(*,'(a)') "usage: CubeFix -cube <name> -pixtable <name> -out <name> [-option <val>]"
       WRITE(*,'(a)') " "
       WRITE(*,'(a)') "  ----- main options:"
       WRITE(*,'(a)') "  -cube               <string>          : datacube to correct (NO DEFAULT)"
       WRITE(*,'(a)') "  -pixtable           <string>          : associated PIXTABLE_REDUCED (NO DEFAULT), it can also be a 'cropped' PIXTABLE"
       WRITE(*,'(a)') "  -out                <string>          : name of the corrected datacube (default=<input_cube_name>_FIX.fits)"
       WRITE(*,'(a)') "                                          use 'none' and '-step 3' to produce ONLY auxiliary maps as described below."
       WRITE(*,'(a)') "  -sourcemask         <string>          : if provided, it looks for a sourcemask image with the name given here (sources are pixels with any value /=0)"
       WRITE(*,'(a)') "  -edgemaskpix        <int>             : if /=0, produces a 2d mask where edgemaskpix pixels from the edges of each slices "
       WRITE(*,'(a)') "                                          have value=0 (default=3). Filename is '<output_cube_name>_SliceEdgeMask.fits'."
       WRITE(*,'(a)') "  -FOVmaskpix         <int>             : if /=0, produces a 2d mask where FOVmaskpix pixels from the edges of the FOV "
       WRITE(*,'(a)') "                                          have value=0 (default=7). This is applied on top of edgemaskpix."
       WRITE(*,'(a)') "  -writemap           <bol>             : if .true. (default) writes an image map where pixel values identify ifu and slice, in the form (ifu*100+slice)"
       WRITE(*,'(a)') "                                          e.g., 1402 is ifu 14 slice 2. This map may be used directly in CubeAdd2Mask."
       WRITE(*,'(a)') "                                          Filename is '<output_cube_name>_IFUSliceMap.fits'."
       WRITE(*,'(a)') "  -step                                 : selects which one of the correction steps will be run:" 
       WRITE(*,'(a)') "                                               1 = slice-by-slice and IFU vertical correction using medium-bands "
       WRITE(*,'(a)') "                                               2 = IFU correction on sky-lines using narrow-bands "
       WRITE(*,'(a)') "                                               3 = 'stack' vertical correction using white-light and edgemask+IFUSliceMap production (if requested)"
       WRITE(*,'(a)') "                                                    use -out 'none' if you want only to get the edgemask and IFUSliceMap without correcting the cube."
       WRITE(*,'(a)') "                                               0 = apply all correction steps (default) "
       WRITE(*,'(a)') " "
       WRITE(*,'(a)') "  ----- advanced options: "
       WRITE(*,'(a)') "  -writecorr          <bol>             : if .true. saves a cube with the applied correction factor (default=.false.)"
       WRITE(*,'(a)') "  -corrname           <string>          : name of the file with PARTIAL correction cube excluding edges (default=<cubename>_CorrFact.fits). To get the "
       WRITE(*,'(a)') "                                           TOTAL correction cube you can simply divide the final Fixed cube by the original unFixed one with CubeArit."
       WRITE(*,'(a)') "  -savecheck          <bol>             : if .true., saves various check images with default names taken from the output cube name (default=.false.)"
       WRITE(*,'(a)') "  -bmin               <real>            : minimum size of wvl bin for the flat-field correction in A (default=5.) before joining the bins below minEW"
       WRITE(*,'(a)') "  -bmax               <real>            : maximum size of wvl bin for the flat-field correction in A (default=30.) before joining the bins below minEW"
       WRITE(*,'(a)') "  -minEW              <real>            : minimum EW size of the joined wvl bin for the flat-field correction (default=200. Angstroms)"
       WRITE(*,'(a)') "  -skythr             <real>            : skylines selection thresholds in units of the median sky (default=1.5)"
       WRITE(*,'(a)') "  -sourceclip         <real>            : image clipping threshold for sources in sigma units (default=2.) if sourcemask is not provided"
       WRITE(*,'(a)') "  -CRclip             <real>            : slice-by-slice clipping threshold for CR and bright pixels in sigma units (default=10.)"
       WRITE(*,'(a)') "                                          slice correction will be replaced by average IFU correction if slice flux is above this threshold."
       WRITE(*,'(a)') "  -edgepix            <int>             : number of pixels from the edge of the slices where the 'IFU vertical' correction is applied,"
       WRITE(*,'(a)') "                                          if this number is larger than the slice width, the correction is applied everywhere (default)"
       WRITE(*,'(a)') "  -edgeFOVpix         <int>             : number of pixels from the edge of the instrument FOV where the 'stack vertical' correction "
       WRITE(*,'(a)') "                                          (STEP 3) is applied (default=10). " 
       WRITE(*,'(a)') "  -writeNaN            <bol>            : if .true. (default) produce NaN instead of UNDEF values for output Cube "
       WRITE(*,'(a)') "                                          (turn this .false. if your compiler complains!)"
       WRITE(*,'(a)') "  -minstatpix         <real>            : minimum number of not UNDEF and off-source voxels to calculate and apply correction factors (default=8)"
       STOP
    END IF
    
    
    !..read from command line
    DO i=1,narg,2
       CALL getarg(i,opt)
       CALL getarg(i+1,arg)
       SELECT CASE(TRIM(opt))
       CASE('-cube')          ; READ(arg,'(a)') datacube
       CASE('-pixtable')      ; READ(arg,'(a)') pixtable
       CASE('-out')           ; READ(arg,'(a)') outputbase
       CASE('-savecheck')     ; READ(arg,*) savechecks
       CASE('-sourcemask')    ; READ(arg,'(a)') MaskFile
       CASE('-bmin')          ; READ(arg,*) skybin_min
       CASE('-bmax')          ; READ(arg,*) skybin_max
       CASE('-minEW')         ; READ(arg,*) minEW
       CASE('-skythr')        ; READ(arg,*) skysel
       CASE('-sourceclip')    ; READ(arg,*) source_clip_sigma
       CASE('-CRclip')        ; READ(arg,*) CRclip
       CASE('-edgepix')       ; READ(arg,*) MinPix
       CASE('-os_fact')       ; READ(arg,*) os_fact
       CASE('-writeNaN')      ; READ(arg,*) writeNaN
       CASE('-writecorr')     ; READ(arg,*) writeCorrCube
       CASE('-writemap')      ; READ(arg,*) writemap
       CASE('-corrname')      ; READ(arg,*) CorrCubeName
       CASE('-edgemaskpix')   ; READ(arg,*) EdgeMask_pix
       CASE('-edgeFOVpix')    ; READ(arg,*) FOV_EdgeNPix
       CASE('-step')          ; READ(arg,*) step_to_run
       CASE('-minstatpix')    ; READ(arg,*) MinStatPixels
       CASE('-writecube')     ; READ(arg,*) writeoutcube
       CASE default
          print *, "command line argument ",TRIM(opt), " not recognized!"
          STOP
       END SELECT
    END DO
    
    !..performs few checks
    IF(TRIM(datacube)=="??") STOP "please provide datacube name with the option -cube!"
    IF(TRIM(pixtable)=="??") STOP "please provide pixtable name with the option -pixtable!"
    IF(TRIM(outputbase)=="none") writeoutcube=.false.
    IF(TRIM(outputbase)=="??".or.TRIM(outputbase)=="none") THEN
       is=INDEX(TRIM(datacube),".fits")
       IF(is==0) is=LEN_TRIM(datacube)+1
       outputbase=TRIM(datacube(1:is-1))//"_FIX"
    ELSE
       is=INDEX(TRIM(outputbase),".fits")
       IF(is==0) is=LEN_TRIM(outputbase)+1
       outputbase=TRIM(outputbase(1:is-1))
    END IF
    IF(TRIM(CorrCubeName)=="??") THEN
       is=INDEX(TRIM(datacube),".fits")
       IF(is==0) is=LEN_TRIM(datacube)+1
       CorrCubeName=TRIM(datacube(1:is-1))//"_CorrFact.fits"       
    END IF

    EDGEMaskName=TRIM(outputbase)//"_SliceEdgeMask.fits"
    MapName=TRIM(outputbase)//"_IFUSliceMap.fits"

  END SUBROUTINE ReadCommandLine

!-----------------------------------------------

  SUBROUTINE ReadPixTable

    IMPLICIT NONE
    REAL(kind=4) :: crpix1, crpix2, xc_cube, yc_cube
    REAL(kind=8) :: xoff_d, yoff_d, xmin1, ymin1, xmin2, ymin2
    REAL(KIND=8), PARAMETER :: pi=3.14159265358979d0
    REAL(KIND=8), PARAMETER :: Rad2Deg=180.0d0/pi, Deg2Rad=pi/180.d0
    REAL(kind=8) :: rx, ry, ra_p, dec_p, phi_p, dx, dy, x, y, r_theta, phi, theta, det, CD_inv(4)


    !..get an unused unit
    CALL ftgiou(unit,status)

    !..get declination from primary header
    fname=TRIM(pixtable)//"[0]"
    CALL ftopen(unit,fname,0,blocksize,status)
    IF(status/=0) THEN
       print *, "unable to open primary header of pixeltable: ", TRIM(fname)
       STOP
    END IF


    !..get ra and declination
    CALL ftgkye(unit,'RA      ',aRA,comment,status)
    IF(status/=0) STOP "RA keyword not found in pixeltable!"
    CALL ftgkye(unit,'DEC     ',aDEC,comment,status)
    IF(status/=0) STOP "DEC keyword not found in pixeltable!"

    CALL ftclos(unit,status)

    !..reopen and move to the first extension
    rwstatus=0   ! 0=readonly
    CALL ftdopn(unit,pixtable,rwstatus,status)
    IF(status/=0) STOP "problem reading file"

    !..get cube size
    CALL ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
    print *, "naxes=",naxes


    !..allocate pos variables
    ALLOCATE(xpos(naxes(1),naxes(2)),&
         ypos(naxes(1),naxes(2)),&
         pos(naxes(1),naxes(2)))

    !..read pixtable variables xpos and ypos in native spherical coordinates
    CALL ReadPixTableVar("xpos")
    CALL ReadPixTableVar("ypos")

    !..deallocate dummy single array pos
    DEALLOCATE(pos)

! ------------- convert pos from celestial spherical to pixel units

    print *, "converting coordinates to pixel units..."

    !..firt, compute the inverse of the CD matrix
    det=WCS(7)*WCS(10)-WCS(8)*WCS(9)
    IF(det==0) STOP "CD matrix has 0 determinant!"
    CD_inv(1)=WCS(10)/det
    CD_inv(2)=WCS(8)/det
    CD_inv(3)=WCS(9)/det
    CD_inv(4)=WCS(7)/det

    !..loop through the pixel table
    DO i=1,SIZE(xpos)

       phi=xpos(1,i)              
       theta=ypos(1,i)+0.5d0*pi  !..the pi/2 is due to MUSE pixtable convention

       !..compute the intermediate world coordinates
       rx=(180.d0/pi)*dsin(phi)/dtan(theta)
       ry=-(180.d0/pi)*dcos(phi)/dtan(theta)

       !..update pixel position
       xpos(1,i)=CD_inv(1)*rx-CD_inv(2)*ry+WCS(1)
       ypos(1,i)=-CD_inv(3)*rx+CD_inv(4)*ry+WCS(2)

    END DO

    print *, "done"

    !.. find center of not-UNDEF cube data using white-light image
    ALLOCATE(image(DimX,DimY))
    image=1.
    DO j=1,DimY
       DO i=1,DimX
          IF(ALL(Cube(i,j,:)==UNDEF)) THEN
             image(i,j)=UNDEF
          ENDIF
       END DO
    END DO
    DO i=1,DimX
       IF(ANY(image(i,:)/=UNDEF)) THEN
          xmin=i-0.5
          EXIT
       END IF
    END DO
    DO i=DimX,1,-1
       IF(ANY(image(i,:)/=UNDEF)) THEN
          xmax=i-0.5
          EXIT
       END IF
    END DO
    DO i=1,DimY
       IF(ANY(image(:,i)/=UNDEF)) THEN
          ymin=i-0.5
          EXIT
       END IF
    END DO
    DO i=DimY,1,-1
       IF(ANY(image(:,i)/=UNDEF)) THEN
          ymax=i-0.5
          EXIT
       END IF
    END DO
    DEALLOCATE(image)
    xc_cube=0.5*(xmin+xmax)
    yc_cube=0.5*(ymin+ymax)

    print *, "--- FOV location within datacube (pixels):"
    print *, "minmax x=", xmin, xmax
    print *, "minmax y=", ymin, ymax
    print *, "center coordinates=", xc_cube, yc_cube


    print *, "--- original FOV pixeltable location (pixels):"
    xmin=MINVAL(xpos); xmax=MAXVAL(xpos)
    ymin=MINVAL(ypos); ymax=MAXVAL(ypos)
    print *, "minmax xpos=", xmin, xmax
    print *, "minmax ypos=", ymin, ymax
    print *, "center coordinates=", 0.5*(xmin+xmax),0.5*(ymin+ymax)

    !..apply offset to have center pixel of pixtable at the same location of
    !..MUSE FOV center in the datacube 
    xoff=xc_cube-0.5*(xmin+xmax)
    yoff=yc_cube-0.5*(ymin+ymax)
    xpos(1,:)=xpos(1,:)+xoff
    ypos(1,:)=ypos(1,:)+yoff
    print *, "--- updated FOV pixeltable location (pixels):"
    xmin=MINVAL(xpos); xmax=MAXVAL(xpos)
    ymin=MINVAL(ypos); ymax=MAXVAL(ypos)
    print *, "minmax xpos=", xmin, xmax
    print *, "minmax ypos=", ymin, ymax    
    print *, "center coordinates=", 0.5*(xmin+xmax),0.5*(ymin+ymax)
    !..set xmin and ymin to the minimum coordinate of the cube, i.e. 0
    xmin=0. ; ymin=0.
    xoff=0. ; yoff=0.


!    !..allocate dummy arrays for coordinate transformation
!    ALLOCATE(dcos_ypos(naxes(1),naxes(2)),dsin_ypos(naxes(1),naxes(2)),dcos_xpos(naxes(1),naxes(2)))

!    print *, "transforming coordinates to celestial spherical..."
!    !..convert from native spherical to celestial spherical (see muse_wcs_position_celestial in muse_wcs.c in muse pipeline)
!    dp = aDEC*deg_to_rad ! delta_p in Paper II (in radians) 
!    ypos = ypos + pi/2.d0 ! add pi/2 again 
!    dcos_ypos=dcos(ypos)
!    dsin_ypos=dsin(ypos)
!    dcos_xpos=dcos(xpos)
!    xpos = -datan2(dcos_ypos * dsin(xpos), dsin_ypos * dcos(dp) + dcos_ypos * dsin(dp) * dcos_xpos)*rad_to_deg
!    ypos = dasin(dsin_ypos * dsin(dp) - dcos_ypos * dcos(dp) * dcos_xpos)*rad_to_deg-aDEC

!    DEALLOCATE(dcos_ypos,dsin_ypos,dcos_xpos)

    ALLOCATE(data(naxes(1),naxes(2)),&
         origin(naxes(1),naxes(2)),&
         dq(naxes(1),naxes(2)),&
         IFU(naxes(1),naxes(2)),&
         slice(naxes(1),naxes(2)),&
         xslice(naxes(1),naxes(2)))

    CALL ReadPixTableVar("data")
    CALL ReadPixTableVar("origin")
    CALL ReadPixTableVar("dq")
    ALLOCATE(lambda(naxes(1),naxes(2)))
    CALL ReadPixTableVar("lambda")
    
    CALL ftclos(unit, status)

    !..convert origin to IFU
    print *, "obtaining IFU information..."
    DO j=1,SIZE(origin)
       IFU(1,j)=IAND(ISHFT(origin(1,j),-MUSE_ORIGIN_SHIFT_IFU),X'1f')
       slice(1,j)=IAND(origin(1,j),X'3f')
       xslice(1,j)=IAND(ISHFT(origin(1,j),-MUSE_ORIGIN_SHIFT_XSLICE),X'7f')
    END DO
    DEALLOCATE(origin)
    min_xslice=MINVAL(xslice)
    max_xslice=MAXVAL(xslice)
    print *, "original min max xslice=",min_xslice,max_xslice

    !..increase image size and resolution by a factor os_fact, if requested
    DimX=DimX*os_fact
    DimY=DimY*os_fact
    xbinsize=xbinsize*os_fact
    ybinsize=ybinsize*os_fact
    xoff=xoff*os_fact
    yoff=yoff*os_fact

  END SUBROUTINE ReadPixTable

!-----------------------------------------

  SUBROUTINE ReadPixTableVar(extname)

    IMPLICIT NONE
    CHARACTER(len=*) :: extname
    INTEGER :: extver, ANY_HDU

    !..move to the selected extension
    extver=0
    status=0
    ANY_HDU=-1
    CALL ftmnhd(unit, ANY_HDU, extname, extver, status)
    IF(status/=0) THEN
       print *, "extension: ", TRIM(extname), "not found!"
       STOP
    END IF
    
    group=1

    print *, "reading extension: ", TRIM(extname)

    SELECT CASE(TRIM(extname))
    CASE("xpos")
       CALL ftgpve(unit,group,1,naxes(2),UNDEF,pos,anyf,status)
       xpos=REAL(pos,KIND=8)
    CASE("ypos")
       CALL ftgpve(unit,group,1,naxes(2),UNDEF,pos,anyf,status)
       ypos=REAL(pos,KIND=8)
    CASE("data")
       CALL ftgpve(unit,group,1,naxes(2),UNDEF,data,anyf,status)
    CASE("lambda")
       CALL ftgpve(unit,group,1,naxes(2),UNDEF,lambda,anyf,status)
    CASE("origin")
       CALL ftgpvj(unit,group,1,naxes(2),UNDEF,origin,anyf,status)
    CASE("dq")
       CALL ftgpvj(unit,group,1,naxes(2),UNDEF,dq,anyf,status)
    END SELECT
    IF(status/=0) STOP "problem reading current extension"

  END SUBROUTINE ReadPixTableVar

!---------------------------------------

 SUBROUTINE SplitCube
    
    IMPLICIT NONE
    REAL(kind=4), ALLOCATABLE :: sky(:),sky_or(:),sky_sigma(:), sourcemask(:,:)
    REAL(kind=4) :: skystep, this_sky_val, clipmean, clipmedian, clipsigma, &
         old_percent, percent, this_EW, avg_sky, right_EW, left_EW
    INTEGER :: this_bin, j, zbins, x, y, jj, right_jj
    INTEGER, PARAMETER :: spsize=2
    REAL(kind=4), ALLOCATABLE :: MaskCube(:,:,:)

    print *, "extracting average sky spectrum..."
    ALLOCATE(sky(SIZE(Cube,DIM=3)),sky_or(SIZE(Cube,DIM=3)), sky_sigma(SIZE(Cube,DIM=3)))

    !..produce a object mask, if not provided as input
    IF(TRIM(MaskFile)=="??") THEN
       print *, "NB: producing a simple source mask from the cube itself"
       ALLOCATE(sourcemask(DimX,DimY),full_sourcemask(DimX,DimY))
       sourcemask=UNDEF
       sourcemask=SUM(Cube,DIM=3,MASK=Cube/=UNDEF)
       DO y=1,DimY
          DO x=1,DimX
             IF(ALL(Cube(x,y,:)==UNDEF)) THEN
                sourcemask(x,y)=UNDEF
             ELSE
                sourcemask(x,y)=Mean(PACK(Cube(x,y,:),MASK=Cube(x,y,:)/=UNDEF))
             END IF
          END DO
       END DO
       CALL SigmaClip(PACK(sourcemask,MASK=sourcemask/=UNDEF),clipmean,clipmedian,clipsigma)
       WHERE(sourcemask>source_clip_sigma*clipsigma+clipmean) sourcemask=UNDEF
    ELSE
       ALLOCATE(MaskCube(DimX,DimY,1),sourcemask(DimX,DimY),full_sourcemask(DimX,DImY))
       CALL ReadLocalCube(MaskFile, MaskCube)
       print *, "Min Max sourcemask=", MINVAL(MaskCube), MAXVAL(MaskCube)
       WHERE(MaskCube(:,:,1)==0) 
          sourcemask=0
       ELSEWHERE
          sourcemask=UNDEF
       END WHERE
    END IF
       
    full_sourcemask=sourcemask

    !ALLOCATE(image(DimX,DimY))
    !image=full_sourcemask
    !CALL WriteImage("testmask.fits")
    !STOP

    !..produce sky spectrum    
    DO j=1,SIZE(sky)
       IF(ANY(Cube(:,:,j)/=UNDEF)) THEN
          sky_or(j)=Mean(PACK(Cube(:,:,j),MASK=Cube(:,:,j)/=UNDEF.and.sourcemask/=UNDEF))
          sky_sigma(j)=StdDev(PACK(Cube(:,:,j),MASK=Cube(:,:,j)/=UNDEF.and.sourcemask/=UNDEF))
       ELSE
          sky_or(j)=UNDEF
       END IF
    END DO


    !..spread sky lines accross several z-pixels to take into account z-variations accross the field
    sky=sky_or
    DO j=spsize+1,SIZE(sky)-spsize
       sky(j)=MAXVAL(sky_or(j-spsize:j+spsize))
    END DO
    WHERE(sky_or==UNDEF) sky=UNDEF
    skystep=1./WCS(11)
    avg_sky=Median(PACK(sky,MASK=sky/=UNDEF))
    print *, "continuum sky level=",avg_sky
    print *, "fraction of the cube with sky lines=",REAL(COUNT(sky>skysel*avg_sky))/SIZE(sky)
    IF(ANY(sky==UNDEF)) print *, "fraction of the cube with NaN layers=", REAL(COUNT(sky==UNDEF))/SIZE(sky)
    
!..build zbins array looping over the sky array 
!..and selecting the minloc of sky between skybin_min and 
!..skybin_max sizes
    print *, "producing wavelength splitting array..." 
!    WRITE(*,'(a,$)') "%"
    this_bin=0
    j=1
    old_percent=0
    this_sky=0.
    this_EW=0.
    !..FIRST passage: splitting using sky minima locations
    print *, "wavelength splitting array:" 
    skybin_min=NINT(skybin_min/WCS(11))
    skybin_max=NINT(skybin_max/WCS(11))
    this_bin=0
    j=1
    zbins=0
    wl(0)%zmax=0
    wl(0)%lmax=WCS(6)
    DO 
       IF(j+skybin_max>SIZE(sky)) EXIT
       binsize=MINLOC(sky(j+INT(skybin_min):j+INT(skybin_max)))+skybin_min
       j=j+binsize(1)
       zbins=zbins+1
       IF(zbins>1000) STOP "more than 1000 wavelength bins requested! Please change skybin parameters!"
       wl(zbins)%zmin=wl(zbins-1)%zmax+1
       wl(zbins)%zmax=j
       wl(zbins)%skip=.false.
    END DO
    IF(zbins==0) THEN !..single bin containing the whole datacube
       zbins=1
       wl(1)%zmin=1 
       wl(1)%zmax=SIZE(Cube,DIM=3)
       wl(1)%skip=.false.
    ELSE  !..build last bin with whatever is left
       zbins=zbins+1
       binsize(1)=SIZE(Cube,DIM=3)-wl(zbins-1)%zmax
       wl(zbins)%zmin=wl(zbins-1)%zmax+1
       wl(zbins)%zmax=SIZE(Cube,DIM=3)
       wl(zbins)%skip=.false.
    END IF
    wl_size=zbins

    !..check which bins need to be joined together according to skysel 
    DO j=1,wl_size
       !..check if edge lands on a sky line
       !..only need to check right edge given the progressive sequence
       IF(sky(wl(j)%zmax)>skysel*avg_sky.and.j<wl_size) THEN !..join with next bin
          !..remove current bin
          wl(j)%skip=.true.
          !..update beginning of the next bin
          wl(j+1)%zmin=wl(j)%zmin
       END IF
    END DO

    !..check which bins need to be joined together according to minEW and skylines
    minEW=minEW/WCS(11) !..convert to pixels
    joinloop:DO j=1,wl_size

       IF(wl(j)%skip) CYCLE
      
       IF(SUM(sky_or(wl(j)%zmin:wl(j)%zmax))<minEW*avg_sky.and.j<wl_size) THEN !..try to join with next bins 

          !..find next "valid" bin on the right
          jj=j+1 
          IF(jj>wl_size) CYCLE joinloop
          IF(wl(jj)%skip) THEN
             DO WHILE(wl(jj)%skip)
                jj=jj+1
                IF(jj>wl_size) CYCLE joinloop
             END DO
          END IF

          IF(MAXVAL(sky_or(wl(jj)%zmin:wl(jj)%zmax))<skysel*avg_sky) THEN !..no sky lines on the right bin

             wl(jj)%zmin=wl(j)%zmin
             !..remove current bin
             wl(j)%skip=.true.

          ELSE !..try to join to the first valid left bin

             !..save values on the right bin
             right_jj=jj
             right_EW=SUM(sky_or(wl(jj)%zmin:wl(jj)%zmax))

             jj=j-1
             IF(jj<1) CYCLE joinloop
             IF(wl(jj)%skip) THEN
                DO WHILE(wl(jj)%skip)
                   jj=jj-1
                   IF(jj<1) CYCLE joinloop
                END DO
             END IF

             IF(MAXVAL(sky_or(wl(jj)%zmin:wl(jj)%zmax))<skysel*avg_sky) THEN !..no sky lines on the left bin

                wl(jj)%zmax=wl(j)%zmax
                !..remove current
                wl(j)%skip=.true.

             ELSE !..sky lines on both sides, join the bin with the smallest EW

                left_EW=SUM(sky_or(wl(jj)%zmin:wl(jj)%zmax))
                IF(left_EW<right_EW) THEN
                   wl(jj)%zmax=wl(j)%zmax
                ELSE
                   wl(right_jj)%zmin=wl(j)%zmin
                END IF
                !..remove current
                wl(j)%skip=.true.                 

             END IF

          END IF

       END IF

    END DO joinloop

    wl_old=wl

    !..cleanup loop
    jj=0
    DO j=1,wl_size
       IF(wl_old(j)%skip) CYCLE
       jj=jj+1
       wl(jj)=wl_old(j)
       wl(jj)%lmin=wl(jj)%zmin*WCS(11)+WCS(6)
       wl(jj)%lmax=wl(jj)%zmax*WCS(11)+WCS(6)
    END DO
    zbins=jj
    wl_size=jj

          
    print *, " "
    print *, "joined wavelength splitting array:"
    print *, " "
    print *, "zbin  zmin     zmax    sky@zmin   sky@zmax   maxsky    EW"
    DO j=1,zbins
       IF(wl(j)%skip) CYCLE
       WRITE(*,'(i3,1x,i6,1x,i6,1x,4(f9.3,1x))') j,wl(j)%zmin, wl(j)%zmax, sky_or(wl(j)%zmin), sky_or(wl(j)%zmax), MAXVAL(sky_or(wl(j)%zmin:wl(j)%zmax)),  SUM(sky_or(wl(j)%zmin:wl(j)%zmax))/avg_sky*WCS(11)  
    END DO
    print *, " "        

!..produce skyline array
    zbins=1
    j=1
    line_wl(1)%zmin=1
    !..cross the sky to find skylines and continuum for wl array for step 2 
    mainloop2: DO WHILE(j<=DimZ)

       IF(sky(j)>skysel*avg_sky.and.j>1.and.j<DimZ) THEN

          !..close previous continuum bin
          line_wl(zbins)%zmax=j-1
          line_wl(zbins)%lmax=line_wl(zbins)%zmax*WCS(11)+WCS(6)

          zbins=zbins+1
          IF(zbins>1000) STOP "increase splitting array size!"

          line_wl(zbins)%zmin=j
          line_wl(zbins)%lmin=line_wl(zbins)%zmin*WCS(11)+WCS(6)

          DO WHILE(sky(j)>skysel*avg_sky)
             j=j+1
             IF(j==DimZ) THEN
                line_wl(zbins)%zmax=j
                line_wl(zbins)%lmax=line_wl(zbins)%zmax*WCS(11)+WCS(6)
                EXIT mainloop2
             END IF
          END DO
          j=j-1

          line_wl(zbins)%zmax=j
          line_wl(zbins)%lmax=line_wl(zbins)%zmax*WCS(11)+WCS(6)

          !..open next continuum bin
          zbins=zbins+1
          line_wl(zbins)%zmin=j+1
          line_wl(zbins)%lmin=line_wl(zbins)%zmin*WCS(11)+WCS(6)

       END IF

       j=j+1

    END DO mainloop2

    !..make sure that last bin closes at the end of the cube
    line_wl(zbins)%zmax=SIZE(Cube,DIM=3)
    line_wl_size=zbins

    print *, " "
    print *, "skyline array:"
    print *, " "
    print *, "zbin  zmin     zmax    sky@zmin   sky@zmax   maxsky    EW"
    DO j=1,zbins
       WRITE(*,'(i3,1x,i6,1x,i6,1x,4(f9.3,1x))') j,line_wl(j)%zmin, line_wl(j)%zmax, sky_or(line_wl(j)%zmin), sky_or(line_wl(j)%zmax), MAXVAL(sky_or(line_wl(j)%zmin:line_wl(j)%zmax)),  SUM(sky_or(line_wl(j)%zmin:line_wl(j)%zmax))/avg_sky*WCS(11)  
    END DO
    print *, " "        


    DEALLOCATE(sky,sky_or,sky_sigma)
    DEALLOCATE(sourcemask)

    
  END SUBROUTINE SplitCube

!------------------------------------------------------

  SUBROUTINE IFU_VerticalCorr

    INTEGER :: this_min, this_max

    IF(count(IFUmap>0.and.image/=UNDEF)<MinStatPixels) STOP "Problem with IFUMap in IFU_VerticalCorr! Check pixtable and offsets."
    CALL SigmaClip(PACK(image,MASK=(IFUmap>0.and.image/=UNDEF)), total_clipmean, total_clipmedian, total_sigma)

    DO j=1,24

       DO group=1,4

          in=12*(group-1)+1
          end=12*group

          min_xslicemap=MINVAL(XSliceMap,MASK=(XSliceMap>0.and.abs(IFUmap)==j.and.SliceMap>=in.and.SliceMap<=end))
          max_xslicemap=MAXVAL(XSliceMap,MASK=(abs(IFUmap)==j.and.SliceMap>=in.and.SliceMap<=end))
      
          DO thispix=1,MIN(MinPix,max_xslicemap-min_xslicemap)
      
             this_min=MAX(min_xslicemap,thispix+min_xslicemap)
             this_max=MIN(max_xslicemap,thispix+min_xslicemap)

             IF(COUNT(IFUmap==j.and.SliceMap>=in.and.SliceMap<=end.and.XSliceMap>=this_min.and.XSliceMap<=this_max.and.image/=UNDEF)<MinStatPixels) CYCLE

             CALL SigmaClip(PACK(image,MASK=(IFUmap==j.and.SliceMap>=in.and.SliceMap<=end.and.XSliceMap>=this_min.and.XSliceMap<=this_max.and.image/=UNDEF)), &
                  this_clipmean, this_clipmedian, this_sigma)
             this_cfact=total_clipmean/this_clipmean

             IF(this_cfact>min_vertical_cfact) THEN
                WHERE(abs(IFUmap)==j.and.SliceMap>=in.and.SliceMap<=end.and.XSliceMap==thispix+min_xslicemap) 
                   image=image*this_cfact
                   corr_image=corr_image*this_cfact
                END WHERE
             END IF
             
          END DO

       END DO
    END DO
      
  END SUBROUTINE IFU_VerticalCorr


!------------------------------------------------------

  SUBROUTINE VerticalCorr

    INTEGER :: this_min, this_max

    CALL SigmaClip(PACK(image,MASK=image/=UNDEF.and.IFUMap>0), total_clipmean, total_clipmedian, total_sigma)
    print *, "image mean=", total_clipmean

    corr_image=1.
    DO group=1,4

       in=12*(group-1)+1
       end=12*group

       min_xslicemap=MINVAL(XSliceMap,MASK=XSliceMap>0.and.SliceMap>=in.and.SliceMap<=end)
       max_xslicemap=MAXVAL(XSliceMap,MASK=SliceMap>=in.and.SliceMap<=end)

       DO j=min_xslicemap,max_xslicemap

          !..skip correction at the edge of the FOV 
          IF(group==1.and.j>max_xslicemap-FOV_EdgeNPix) CYCLE
          IF(group==4.and.j<min_xslicemap+FOV_EdgeNPix) CYCLE

          IF(COUNT(XSLiceMap==j.and.SliceMap>=in.and.SliceMap<=end.and.image/=UNDEF.and.IFUMap>0)<12*MinStatPixels) CYCLE
          CALL SigmaClip(PACK(image,MASK=(XSliceMap==j.and.image/=UNDEF.and.SliceMap>=in.and.SliceMap<=end.and.IFUMap>0)), &
               this_clipmean, this_clipmedian, this_sigma)
          this_cfact=total_clipmean/this_clipmean
          WHERE(XSliceMap==j.and.SliceMap>=in.and.SliceMap<=end) corr_image=this_cfact

       END DO
       
    END DO

  END SUBROUTINE VerticalCorr

!-----------------------------------

  SUBROUTINE transform(pix_x, pix_y, pix_ra, pix_dec)

    IMPLICIT NONE
    REAL(kind=8), INTENT(IN)  :: pix_x, pix_y
    REAL(kind=8), INTENT(OUT) :: pix_ra, pix_dec 
    REAL(KIND=8), PARAMETER :: pi=3.14159265358979d0
    REAL(KIND=8), PARAMETER :: Rad2Deg=180.0d0/pi, Deg2Rad=pi/180.d0
    REAL(kind=8) :: rx, ry, ra_p, dec_p, phi_p, dx, dy, x, y, r_theta, phi, theta

!--- compute transformation

    ! use notation more like WCS paper
    rx = WCS(1)
    ry = WCS(2)
    ra_p =  WCS(4)*Deg2Rad
    dec_p = WCS(5)*Deg2Rad

    ! native longitude of the celestial pole in radians
    phi_p = pi
    IF(dec_p >= pi/2.d0) phi_p = 0.0

    ! compute intermediate world coordinates
    dx = pix_x - rx
    dy = pix_y - ry
    x = WCS(7)*dx + WCS(8)*dy
    y = WCS(9)*dx + WCS(10)*dy

    ! compute the "native spherical coordinates" in radians
    r_theta = dsqrt(x*x + y*y)
    phi = datan2(x, -y)
    theta = datan2(180.d0, pi*r_theta)

    pix_ra=phi
    pix_dec=theta+0.5d0*pi  !..MUSE convention

    !..get ra & dec:
    pix_ra=ra_p + datan2(-dcos(theta)*dsin(phi-phi_p), &
         dsin(theta)*dcos(dec_p)-dcos(theta)*dsin(dec_p)*dcos(phi-phi_p))
    pix_dec = dasin(dsin(theta)*dsin(dec_p)+dcos(theta)*dcos(dec_p)*dcos(phi - phi_p))

    ! convert ra, dec to degrees
    pix_ra=pix_ra*Rad2Deg
    IF(pix_ra<0.d0) pix_ra=pix_ra+360.
    pix_dec=pix_dec*Rad2Deg


  END SUBROUTINE transform

!........................................

  SUBROUTINE RaDec_to_pixel(pix_ra, pix_dec, pix_x, pix_y, native_spherical)

    IMPLICIT NONE
    REAL(kind=8), INTENT(OUT)  :: pix_x, pix_y
    REAL(kind=8), INTENT(IN)   :: pix_ra, pix_dec 
    LOGICAL, INTENT(IN)        :: native_spherical
    REAL(KIND=8), PARAMETER :: pi=3.14159265358979d0
    REAL(KIND=8), PARAMETER :: Rad2Deg=180.0d0/pi, Deg2Rad=pi/180.d0
    REAL(kind=8) :: rx, ry, ra_p, dec_p, phi_p, dx, dy, x, y, r_theta, phi, theta, det, CD_inv(4)


!..compute the inverse of the CD matrix
    det=WCS(7)*WCS(10)-WCS(8)*WCS(9)
    IF(det==0) STOP "CD matrix has 0 determinant!"
    CD_inv(1)=WCS(10)/det
    CD_inv(2)=WCS(8)/det
    CD_inv(3)=WCS(9)/det
    CD_inv(4)=WCS(7)/det

    IF(.not.native_spherical) THEN !..convert to native spherical coordinates
       STOP "conversion from projected tangential not implemented yet! Use native_spherical=.true."
    ELSE
       phi=pix_ra+xoff
       theta=pix_dec+0.5d0*pi+yoff  !..the pi/2 is due to MUSE pixtable convention
    END IF

!..compute the intermediate world coordinates
    rx=(180.d0/pi)*dsin(phi)/dtan(theta)
    ry=-(180.d0/pi)*dcos(phi)/dtan(theta)

!..compute pixel position
    pix_x=CD_inv(1)*rx-CD_inv(2)*ry+WCS(1)
    pix_y=-CD_inv(3)*rx+CD_inv(4)*ry+WCS(2)

  END SUBROUTINE RaDec_to_pixel

END PROGRAM CubeFix

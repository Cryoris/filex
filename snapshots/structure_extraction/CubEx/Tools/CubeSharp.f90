!----------------------------------------------
!
! CubeSharp - SHfit and Add sky Removal Procedure
!
! Author: SC
! Last update: Nov 29 2017
!
PROGRAM CubeSharp

  USE CubeLib
  USE StatLib
  !USE ieee_arithmetic
  IMPLICIT NONE
  CHARACTER(len=250) :: inpname, outname, skytype, outshift, outflux, outsources, MaskFile, hSNCube, sourcecube
  REAL(kind=4), ALLOCATABLE :: skyref(:), sky(:), sourcemask(:,:), ShiftCube(:,:,:), shiftmask(:,:), skyor(:), this_sky(:), avgshift(:), &
       full_sourcemask(:,:), skyref_var(:), CubeOR(:,:,:), Cube_sources(:,:,:), totvar_or(:), totvar(:), skyweight(:), refvar(:), &
       UnitCube(:,:,:), VarOR(:,:,:), Mean_UnitCube(:), Sigma_UnitCube(:), SourceContinuum(:,:,:), SkySub0(:,:,:), SkyRatio(:,:,:), MaskCube(:,:,:)
  LOGICAL, ALLOCATABLE :: skipmask(:,:,:)
  LOGICAL :: linecheck
  INTEGER(kind=4) :: DimX, DimY, DimZ, skyref_x, skyref_y, x, y, z, zbins, mbins, n_ext, lc, it, notconv, ng, k, ncount, wl_size, &
       z1, z2, line, zmin, zmax, percent, old_percent, didconv, zfinebins, old_notconv, conv_maxiterations, &
       conv_miniterations, skipped, x1, x2, y1, y2, nneg, nneg_nc, npixcont, nskipObjId
  INTEGER(kind=4), ALLOCATABLE :: SkipObjId(:)
  REAL(kind=4)    :: zshift, maxshift, zstep, maxmag, magstep, magfact, minmag, clipmean, clipmedian, clipsigma, fshiftsel, this_sky_val, &
       normfact, conv, old_lc_sky, this_shift, maxtotalshift, skyref_frac, conv_threshold, maxfineshift, zfinestep, skysel, avg_sky, skip_fact, &
       clipval, convfrac, refnorm, minvar, var_fact, avg_var_over_sky, new_skyref_var, realtype, skydeblend, MaxFluxCorr, sclipWL, sclipNB, this_sky_val_or, &
       source_line_to_sky_ratio, source_line_to_cont_ratio, source_line_significance
  LOGICAL :: multiext, applyfineshift, writeNaN, apply_avg_shift, shiftcheck, update_refsky, update_refsky_var, subtract_sky_continuum
  LOGICAL, ALLOCATABLE :: Converged(:,:)
  REAL(kind=4), PARAMETER :: skyalpha=-1.0
 TYPE wlarray
    INTEGER(kind=4) :: zmin, zmax
 END type wlarray
 TYPE(wlarray) :: wl(1000)

  !..get parameters
  CALL ReadParam
  !---------------

  !..get number of large and fine bins in z direction
  zbins=2*maxshift/zstep
  zfinebins=2*maxfineshift/zfinestep

  !..read cube and variance (if present)
  VERBOSITY=2
  IF(multiext) THEN
     CALL ReadCube(inpname,n_ext=2)
  ELSE
     CALL ReadCube(inpname,n_ext=1)
     ALLOCATE(Var(SIZE(Cube,DIM=1),SIZE(Cube,DIM=2),SIZE(Cube,DIM=3)))
     Var=1.
  END IF
  DimX=SIZE(Cube,DIM=1); DimY=SIZE(Cube,DIM=2); DimZ=SIZE(Cube,DIM=3)
  ALLOCATE(ShiftCube(DimX,DimY,DimZ),shiftmask(DimX,DimY),&
       totvar_or(DimZ),totvar(DimZ),Converged(DimX,DimY))
  IF(shiftcheck) THEN
     ALLOCATE(UnitCube(DimX,DimY,Dimz),Mean_UnitCube(DimZ),Sigma_UnitCube(DimZ))
     UnitCube=1.
     Mean_UnitCube=1.
     Sigma_UnitCube=0.
  END IF
  ShiftCube=0.
  shiftmask=0
  print *, "MinMax Cube=", MINVAL(Cube,MASK=Cube/=UNDEF), MAXVAL(Cube,MASK=Cube/=UNDEF)
  print *, "MinMax Var=", MINVAL(Var,MASK=Var/=UNDEF), MAXVAL(Var,MASK=Var/=UNDEF)

  !..find sky lines, remove sky continuum and remove 0th order sky lines
  CALL SplitCube

  print *, "removing continuum of sources..."
  CALL RemoveContSources
  print *, "done"

  ALLOCATE(sourcemask(DimX,DimY))

  ! ------- MAIN loop through skylines -----------

  DO line=1,wl_size
    

     !..select initial and final spectral pixel
     z1=wl(line)%zmin
     z2=wl(line)%zmax

     print *, " "
     print *, "working on line n=", line, "over", wl_size
     print *, "z1=",z1,"z2=",z2

     IF(ALLOCATED(skyref)) DEALLOCATE(skyref,sky,skyref_var, Cube_sources, skyweight, skipmask, SkyRatio, CubeOR, VarOR)
     ALLOCATE(skyref(z1:z2),sky(z1:z2),skyref_var(z1:z2), Cube_sources(DimX,DimY,z1:z2), skyweight(z1:z2), &
          skipmask(DimX,DimY,z1:z2),SkyRatio(DimX,DimY,z1:z2), CubeOR(DimX,DimY,z1:z2),VarOR(DimX,DimY,z1:z2))
     CubeOR(:,:,z1:z2)=Cube(:,:,z1:z2)
     VarOR(:,:,z1:z2)=Var(:,:,z1:z2)

     !..put temporarily undefined Var pixel to 0 to avoid problems when shifting
     !..we will restore UNDEF values at the end
     WHERE(Var(:,:,z1:z2)==UNDEF) Var(:,:,z1:z2)=0.

     !..compute the sky contribution ratio
     WHERE(Cube(:,:,z1:z2)/=UNDEF.and.Cube(:,:,z1:z2)/=0.and.SourceContinuum(:,:,z1:z2)/=UNDEF) 
        SkyRatio(:,:,z1:z2)=(Cube(:,:,z1:z2)-SourceContinuum(:,:,z1:z2))/Cube(:,:,z1:z2)
     ELSEWHERE
        SkyRatio(:,:,z1:z2)=1.0
     END WHERE
     WHERE(SkyRatio<0.) SkyRatio=0.

     !..remove source continuum level from cube and save them in a different cube (section)
     sourcemask=full_sourcemask
     DO z=z1,z2
        WHERE(sourcemask==1.and.SourceContinuum(:,:,z)/=UNDEF)
           Cube(:,:,z)=Cube(:,:,z)-SourceContinuum(:,:,z)
        END WHERE
        Cube_sources(:,:,z)=CubeOR(:,:,z)
     END DO

     !..restore undefined voxels
     WHERE(CubeOR(:,:,z1:z2)==UNDEF) Cube(:,:,z1:z2)=UNDEF

     !----- build  ref sky --------
     !..finds line center first
     DO z=z1,z2
        sky(z)=Mean(PACK(Cube(:,:,z),MASK=Cube(:,:,z)/=UNDEF.and.sourcemask==0.))
     END DO
     lc=MAXLOC(sky,DIM=1)+z1-1
     print *, "line center= ",lc, lc*WCS(11)+WCS(6),"sky value= ",sky(lc)

     print *, "number of spaxels for ref sky=", COUNT(abs(Cube(:,:,lc)-sky(lc))/sky(lc)<skyref_frac)
     IF(COUNT(abs(Cube(:,:,lc)-sky(lc))/sky(lc)<skyref_frac)<10) STOP "change ref sky selection parameters!"

     !..ref sky
     DO z=z1,z2
        CALL SigmaClip(array=PACK(Cube(:,:,z),MASK=Cube(:,:,z)/=UNDEF.and.sourcemask==0..and.abs(Cube(:,:,lc)-sky(lc))/sky(lc)<skyref_frac),&
             MeanClip=skyref(z), MedianClip=clipmedian, FinalSigma=skyref_var(z), ClipVal=[-10.,10.])
        skyref_var(z)=skyref_var(z)*skyref_var(z)
     END DO
     refnorm=1./SUM(skyref)
     skyweight(:)=skyref(:)**skyalpha

     !..check if we are going to use the fineshift procedure
     applyfineshift=ANY(skyref>fshiftsel*avg_sky)

     conv=1.e9
     old_lc_sky=lc
     it=0
     notconv=0
     Converged(:,:)=.true.

     !--- iteration loop -----

     iteration_loop: DO WHILE(it<conv_miniterations.or.(conv>conv_threshold.and.it<conv_maxiterations))

        it=it+1
        IF(.not.applyfineshift) WRITE(*,'(a,$)') "."


        !..build mask to check which spectra we can skip
        DO z=z1,z2
           IF(skyref_var(z)/skyref(z)<3*avg_var_over_sky) THEN
              WHERE(Cube(:,:,z)/=UNDEF) 
                 skipmask(:,:,z)=((Cube(:,:,z)-skyref(z))**2<skip_fact*skyref_var(z))
              ELSEWHERE
                 skipmask(:,:,z)=.true.
              END WHERE
           ELSE 
              WHERE(Cube(:,:,z)/=UNDEF)
                 skipmask(:,:,z)=.false.
              ELSEWHERE
                 skipmask(:,:,z)=.true.
              END WHERE
           END IF
       END DO

       

        !..skip also spectra with negative values (either because of source continuum removal or problems)
        !..or undefined SourceContinuum
        !..and set them as unconverged
        DO y=1,DimY
           DO x=1,DimX
              IF(ANY(Cube(x,y,z1:z2)<=0..and.Cube(x,y,z1:z2)/=UNDEF).or.ANY(SourceContinuum(x,y,z1:z2)==UNDEF)) THEN
                 skipmask(x,y,z1:z2)=.true.
                 Converged(x,y)=.false.
              END IF
           END DO
        END DO

        
        skipped=0
        notconv=0
        didconv=0
        !..performs an overall shift
        DO y=1,DimY
           DO x=1,DimX
              IF(ANY(Cube(x,y,z1:z2)==UNDEF).or.ANY(ShiftCube(x,y,z1:z2)==UNDEF)) THEN
                 !Converged(x,y)=.false.
                 CYCLE
              END IF
              IF(ALL(skipmask(x,y,z1:z2))) THEN
                 skipped=skipped+1
                 !Converged(x,y)=.false.
                 CYCLE
              END IF
              CALL SpcShift(x,y,zshift)
           END DO
        END DO

        !..recompute the sky contribution ratio
        WHERE(Cube_sources(:,:,z1:z2)/=UNDEF.and.Cube_sources(:,:,z1:z2)/=0.and.SourceContinuum(:,:,z1:z2)/=UNDEF) 
           SkyRatio(:,:,z1:z2)=(Cube_sources(:,:,z1:z2)-SourceContinuum(:,:,z1:z2))/Cube_sources(:,:,z1:z2)
        ELSEWHERE
           SkyRatio(:,:,z1:z2)=1.0
        END WHERE
        WHERE(SkyRatio<0.) SkyRatio=0.



        !..performs individual pixel shifts if skyline is bright enough
        IF(applyfineshift) THEN

           print *, "applying fine shifts..."
           zmin=z1 ; zmax=z2
           print *, "zmin=",zmin, "zmax=",zmax

           old_notconv=notconv
           notconv=0
           didconv=0
           old_percent=0
           WRITE(*,'(a,$)') "% "
           DO y=1,DimY
              percent=INT(REAL(y)/DimY*10)
              IF(percent/=old_percent) WRITE(*,'(i4,$)') percent*10 
              DO x=1,DimX
                 IF(ANY(Cube(x,y,zmin:zmax)==UNDEF).or.ANY(ShiftCube(x,y,zmin:zmax)==UNDEF).or.ALL(skipmask(x,y,z1:z2))) THEN
                    !Converged(x,y)=.false.
                    CYCLE
                 END IF
                 CALL SpcFineShift(x,y,zmin,zmax)
              END DO
              old_percent=percent
           END DO
           print *, " "
           print *, "convergence not reached for n=",notconv, "spectra over",notconv+didconv
        ELSE
           notconv=0
           didconv=1
        END IF

        !..recompute the sky contribution ratio
        WHERE(Cube_sources(:,:,z1:z2)/=UNDEF.and.Cube_sources(:,:,z1:z2)/=0.and.SourceContinuum(:,:,z1:z2)/=UNDEF) 
           SkyRatio(:,:,z1:z2)=(Cube_sources(:,:,z1:z2)-SourceContinuum(:,:,z1:z2))/Cube_sources(:,:,z1:z2)
        ELSEWHERE
           SkyRatio(:,:,z1:z2)=1.0
        END WHERE
        WHERE(SkyRatio<0.) SkyRatio=0.

        !--UPDATE reference SKY if requested

        !..finds line center
        DO z=z1,z2
           sky(z)=Mean(PACK(Cube(:,:,z),MASK=Cube(:,:,z)/=UNDEF.and.sourcemask==0.))
        END DO
        lc=MAXLOC(sky,DIM=1)+z1-1
        IF(applyfineshift) print *, "line center= ",lc,lc*WCS(11)+WCS(6),"sky value= ",sky(lc)
     
        !..update ref sky
        IF(update_refsky) THEN
            DO z=z1,z2
              CALL SigmaClip(array=PACK(Cube(:,:,z),MASK=Cube(:,:,z)/=UNDEF.and.sourcemask==0..and.abs(Cube(:,:,lc)-sky(lc))/sky(lc)<skyref_frac),&
                   MeanClip=skyref(z), MedianClip=clipmedian, FinalSigma=new_skyref_var, ClipVal=[-10.,10.])
              new_skyref_var=new_skyref_var**2
              IF(update_refsky_var) skyref_var(z)=new_skyref_var
           END DO
        ENDIF
        skyweight(:)=skyref(:)**skyalpha
        
        !..check convergence
        IF(.not.applyfineshift) THEN
           conv=abs((old_lc_sky-sky(lc))/sky(lc))
        ELSE
           IF(REAL(notconv)/MAX(1,(notconv+didconv))>conv_threshold) THEN
              conv=MAX(abs((old_lc_sky-sky(lc))/sky(lc)),abs(old_notconv-notconv)/REAL(notconv))
           ELSE
              conv=abs((old_lc_sky-sky(lc))/sky(lc))
           END IF
        END IF
           
        old_lc_sky=sky(lc)

        !..check convergence on skyref variance
        convfrac=0.
        DO z=z1,z2
           convfrac=MAX(convfrac,COUNT((Cube(:,:,z)-skyref(z))**2<skip_fact*skyref_var(z).and.Cube(:,:,z)/=UNDEF)/REAL(COUNT(Cube(:,:,z)/=UNDEF)))
        END DO
        conv=MIN(conv, 1.-convfrac)
     
     END DO iteration_loop

    !..put back in the shifted spectra with sources
     DO y=1,DimY
        DO x=1,DimX
           IF(sourcemask(x,y)==1) Cube(x,y,z1:z2)=Cube_sources(x,y,z1:z2)
        END DO
     END DO

     IF(shiftcheck) THEN
        DO z=z1,z2
           Mean_UnitCube(z)=Mean(PACK(UnitCube(:,:,z),MASK=Cube(:,:,z)/=UNDEF.and.sourcemask(:,:)==0.))
           Sigma_UnitCube(z)=StdDev(PACK(UnitCube(:,:,z),MASK=Cube(:,:,z)/=UNDEF.and.sourcemask(:,:)==0.))
        END DO
     END IF

     DO y=1,DimY
        DO x=1,DimX

           IF(ANY(ShiftCube(x,y,z1:z2)==UNDEF)) THEN
              Cube(x,y,z1:z2)=CubeOR(x,y,z1:z2)
              Var(x,y,z1:z2)=VarOR(x,y,z1:z2)
              IF(shiftcheck) UnitCube(x,y,z1:z2)=1.
              Converged(x,y)=.false.
           END IF

           !..do the same if UnitCube corrections, if available, on 
           !..sources got too large (indicative of an emission line or other problem)
           IF(shiftcheck.and.sourcemask(x,y)==1) THEN
              IF(ANY(abs(UnitCube(x,y,z1:z2)-Mean_UnitCube(z1:z2))>MaxFluxCorr*Sigma_UnitCube(z1:z2))) THEN
                 Cube(x,y,z1:z2)=CubeOR(x,y,z1:z2)
                 Var(x,y,z1:z2)=VarOR(x,y,z1:z2)
                 UnitCube(x,y,z1:z2)=1.
                 IF(line==17) Image(x,y)=1
                 Converged(x,y)=.false.
              END IF
           END IF

        END DO
     END DO


     IF(ALLOCATED(this_sky)) DEALLOCATE(this_sky,refvar)
     ALLOCATE(this_sky(z1:z2),refvar(z1:z2))

     ! ------- FINALLY remove sky --------------------------
     print *, "z       new_sky    old_sky    "
     DO z=z1,z2
         SELECT CASE(TRIM(skytype))
        CASE("mean")
           this_sky_val=Mean(PACK(Cube(:,:,z),MASK=Cube(:,:,z)/=UNDEF.and.sourcemask==0.))
           this_sky_val_or=Mean(PACK(CubeOR(:,:,z),MASK=CubeOR(:,:,z)/=UNDEF.and.sourcemask==0.))
        CASE("median") 
           this_sky_val=Median(PACK(Cube(:,:,z),MASK=Cube(:,:,z)/=UNDEF.and.sourcemask==0.))
           this_sky_val_or=Median(PACK(CubeOR(:,:,z),MASK=CubeOR(:,:,z)/=UNDEF.and.sourcemask==0.))
        CASE("avgsigclip") 
           CALL SigmaClip(PACK(Cube(:,:,z),MASK=Cube(:,:,z)/=UNDEF.and.sourcemask==0.), this_sky_val, clipmedian, clipsigma,ClipVal=[-clipval,clipval])
           CALL SigmaClip(PACK(CubeOR(:,:,z),MASK=CubeOR(:,:,z)/=UNDEF.and.sourcemask==0.), this_sky_val_or, clipmedian, clipsigma,ClipVal=[-clipval,clipval])
        CASE("medsigclip")
           CALL SigmaClip(PACK(Cube(:,:,z),MASK=Cube(:,:,z)/=UNDEF.and.sourcemask==0.), clipmean, this_sky_val, clipsigma,ClipVal=[-clipval,clipval])
           CALL SigmaClip(PACK(CubeOR(:,:,z),MASK=CubeOR(:,:,z)/=UNDEF.and.sourcemask==0.), clipmean, this_sky_val_or, clipsigma,ClipVal=[-clipval,clipval])
        CASE("nosub")
           this_sky_val=0.
        CASE default
           STOP "selected skytype is not available!"
        END SELECT
        this_sky(z)=this_sky_val
        totvar(z)=clipsigma**2
        WHERE(Cube(:,:,z)/=UNDEF.and.Converged(:,:)) Cube(:,:,z)=Cube(:,:,z)-this_sky_val
        WHERE(Cube(:,:,z)/=UNDEF.and..not.Converged(:,:)) Cube(:,:,z)=Cube(:,:,z)-this_sky_val_or
        print *, z,this_sky_val,this_sky_val_or

      END DO

      !..restore UNDEF values in the original Var cube
      WHERE(VarOR(:,:,z1:z2)==UNDEF) Var(:,:,z1:z2)=UNDEF

      print *, "total number of spectra that have not converged / been shifted =", COUNT(Converged.eqv..false.)

   END DO
   
   print *, " "


   !..write Cube and Variance using original header information
   CALL UpdateCube(Inpfile=inpname,OutFile=outname, writeNaN=writeNaN)

   !..write shift cube if requested
   IF(TRIM(outshift)/="??") THEN
      Cube=ShiftCube
      CALL WriteCube(outshift,n_ext=1)
   END IF

   !..write UnitCube for flux shift checks if requested
   IF(shiftcheck.and.TRIM(outflux)/="??") THEN
      Cube=UnitCube
      CALL WriteCube(outflux,n_ext=1)
   END IF

!  IF(TRIM(outsources)/="??") THEN
!     Cube=Cube_Sources
!     CALL WriteCube(outsources,n_ext=1)
!  END IF


 CONTAINS

  SUBROUTINE ReadParam

    IMPLICIT NONE
    CHARACTER(len=500) :: ExeName, fname, string, ds, opt, arg, skipobjid_longstring
    INTEGER :: narg, iarg, i, is, lastarg, i_linecheck, testval(1000), ierr
    LOGICAL :: ex

!..default parameters
    inpname="??"
    outname="??"
    outshift="??"
    outflux="??"
    outsources="??"
    maxshift=0.99
    maxtotalshift=1.5
    MaxFluxCorr=3.0
    zstep=0.01
    maxfineshift=0.99
    zfinestep=0.01
    fshiftsel=1.5
    skysel=1.5
    skydeblend=5.0
    ng=3
    skyref_frac=1.e9
    conv_threshold=0.005
    conv_maxiterations=30
    conv_miniterations=3
    skytype="avgsigclip"
    clipval=3.0
    sclipWL=3.0
    writeNaN=.true.
    skip_fact=1
    apply_avg_shift=.false.
    shiftcheck=.false.
    update_refsky=.true.
    update_refsky_var=.false.
    subtract_sky_continuum=.true.
    MaskFile="??"
    hSNCube="??"
    sourcecube="??"
    source_line_to_sky_ratio=0.001
    source_line_to_cont_ratio=2.0
    !linecheck=.true.
    i_linecheck=0 !..this flags checks if the option is provided
    npixcont=50
    nSkipObjId=0
    skipobjid_longstring="??"

    CALL GetArg(0,ExeName)
    narg=iargc()
    lastarg=0
    IF(narg==1) THEN
       CALL GetArg(1,inpname)
       lastarg=1
    ELSEIF(narg==0) THEN
       print *, " "
       WRITE(*,'(2a)')"       CubeSharp (part of CubEx package) ",TRIM(version)
       WRITE(*,'(a)') "      SHift and Add sky Removal Procedure   "
       WRITE(*,'(a)') " "
       WRITE(*,'(a)') " by Sebastiano Cantalupo (cantalupo@phys.ethz.ch)"
       print *, " "
       print *, " "
       print *, "usage: CubeSharp.x -cube <inpcube> -out <outcube> [-option ...]" 
       print *, "or     CubeSharp.x <inpcube>"
       print *, " "
       print *, "------- main options: "
       WRITE(*,'(a)') " -cube        <name>    : input cube name (NO DEFAULT)"
       WRITE(*,'(a)') " -out         <name>    : output cube name (default=<inpcubename>_SHARP.fits"
       WRITE(*,'(a)') " -outshift    <name>    : if requested, saves a cube with the total shifts with this name"
       WRITE(*,'(a)') " -outflux     <name>    : if requested, saves a cube with the total, normalized shifted flux per pixel for a continuum-flat source with this name,"
       WRITE(*,'(a)') "                        this correction is applied by default to the final cube"
       WRITE(*,'(a)') " -sourcemask  <name>    : if provided, it looks for a sourcemask image with the name given here (sources are pixels with any value /=0)"
       WRITE(*,'(a)') " -hsncube     <name>    : if provided, it uses this (possibly higher-SN) sky-subtracted cube as initial guess to remove continuum sources"
       WRITE(*,'(a)') "                          with median filtering. NB: the result of a previous CubeSharp and CubeCombine run can be used here."
       WRITE(*,'(a)') "                          If not provided, the original cube after a 0th order sky subtraction is used instead. See also -skipobjid below."
!       WRITE(*,'(a)') " -sourcecube  <name>    : if provided, it uses this cube to remove continuum sources instead of estimating the spectrum."
!       WRITE(*,'(a)') "                          typically, the result of a previous CubeSharp and CubeCombine run can be used here."
       WRITE(*,'(a)') " "
       WRITE(*,'(a)') "------- advanced options:"
       WRITE(*,'(a)') " -writeNaN    <bol>     : if .true. (default) produce NaN instead of UNDEF values for output Cube (turn this .false. if your compiler complains!)"
       WRITE(*,'(a)') " -skytype     <string>  : sky wavelength averaging method after spectral shifting, "
       WRITE(*,'(a)') "                        options: mean, median, avgsigclip (default), medsigclip"
       WRITE(*,'(a)') " -clipval     <real>    : abs sigma clip value for sky averaging method after spectral shifting (default=3.0)"
       WRITE(*,'(a)') " -sclipWL     <real>    : source SNR threshold in White Light image (default=3.0) if sourcemask is not provided."
       WRITE(*,'(a)') " -skythr      <real>    : skyline selection threshold in units of median sky (default=1.5)"
       WRITE(*,'(a)') " -finethr     <real>    : skyline selection threshold in units of median sky for additional fine shift method (default=1.5)"
       WRITE(*,'(a)') " -skydeb      <real>    : skyline deblending threshold in units of median sky (default=5.0), skylines will be deblended if "
       WRITE(*,'(a)') "                        there is a minimum between skylines that falls below this threshold."
       WRITE(*,'(a)') " -fskyref     <real>    : selection parameter for the spectra that contribute to the 'reference sky' "
       WRITE(*,'(a)') "                        spectra are selected if their line center flux differs less than fskyref*skyavg (default=use all cube)"
       WRITE(*,'(a)') " -convthr     <real>    : convergence threshold for the line center sky to stop iterative procedure (default=0.005)"
       WRITE(*,'(a)') " -skipfac     <real>    : threshold in units of **sigma^2** with respect to reference sky to skip shifting each individual pixel (default=1)"
       WRITE(*,'(a)') " -maxiter     <int>     : maximum number of iterations per sky line (default=30)"
       WRITE(*,'(a)') " -miniter     <int>     : minimum number of iterations per sky line (default=3)"
       WRITE(*,'(a)') " -maxtotshift <real>    : maximum total shift allowed in pixel units (default=1.5)"
       WRITE(*,'(a)') " -maxfshift   <real>    : maximum fine shift per iteration in pixel units (default=0.99). NB: cannot be larger than 1 pixel!"
       WRITE(*,'(a)') " -zstep       <real>    : fine shift step size in pixel units (default=0.01)"
       WRITE(*,'(a)') " -ncont       <int>     : number of spectral pixels for median-filtering continuum estimation of sources (default=50; radius)"
       WRITE(*,'(a)') " -skipobjid <int array> : list of ids (on the same line) of objects for which continuum-subtraction and line feature check "
       WRITE(*,'(a)') "                          from -hsncube ARE DISABLED. The id values found in the sourcemaks (-sourcemask) are used. "
       WRITE(*,'(a)') "                          NB: This option is needed, e.g., for variable objects like QSOs. "
       WRITE(*,'(a)') '                          usage example: ... -sourcemaks SourceMask.fits -hsncube DATACUBE.fits -skipobjid "363 256 780" '
       WRITE(*,'(a)') "  "
       WRITE(*,'(a)') " ------- line feature check options. All the conditions below need to be met to define a line feature (i.e., no shifts and 0th-order sky-sub):"
       WRITE(*,'(a)') "          NB: use of linecheck option is discouraged if a higher SN cube is not provided (option -hsncube above) "
       WRITE(*,'(a)') "              if -hsncube is used, be careful with variable sources (use -skipobjid option above). "
       WRITE(*,'(a)') " -lcheck      <bol>     : main switch for line feature check (default=.true. if hsncube is provided, .false. otherwise)"
       WRITE(*,'(a)') " -lsig        <real>    : source line significance (in terms of sigma layer-by-layer) to trigger line feature checks (default=3.0)"
       WRITE(*,'(a)') " -l2sky       <real>    : ratio between source and sky flux to trigger line feature checks (default=0.001)"
       WRITE(*,'(a)') " -l2con       <real>    : ratio between line and continuum flux to trigger line feature checks (default=2.0)"
       WRITE(*,'(a)') " "
       STOP
    ENDIF

    !..read from command line
    DO i=lastarg+1,narg,2
       CALL getarg(i,opt)
       CALL getarg(i+1,arg)
       SELECT CASE(TRIM(opt))
       CASE('-cube')          ; READ(arg,'(a)') inpname
       CASE('-out')           ; READ(arg,'(a)') outname
       CASE('-skytype')       ; READ(arg,'(a)') skytype
       CASE('-skythr')        ; READ(arg,*) skysel
       CASE('-skydeb')        ; READ(arg,*) skydeblend
       CASE('-clipval')       ; READ(arg,*) clipval   
       CASE('-sclipWL')       ; READ(arg,*) sclipWL
       CASE('-finethr')       ; READ(arg,*) fshiftsel
       CASE('-fskyref')       ; READ(arg,*) skyref_frac
       CASE('-convthr')       ; READ(arg,*) conv_threshold
       CASE('-maxiter')       ; READ(arg,*) conv_maxiterations
       CASE('-miniter')       ; READ(arg,*) conv_miniterations
       CASE('-skipfac')       ; READ(arg,*) skip_fact
       CASE('-maxtotshift')   ; READ(arg,*) maxtotalshift
       CASE('-maxfshift')     ; READ(arg,*) maxfineshift
       CASE('-zstep')         ; READ(arg,*) zfinestep
       CASE('-writeNaN')      ; READ(arg,*) writeNaN
       CASE('-outshift')      ; READ(arg,*) outshift
       CASE('-outflux')       ; READ(arg,*) outflux
       CASE('-outsources')    ; READ(arg,*) outsources
       CASE('-shiftcheck')    ; READ(arg,*) shiftcheck
       CASE('-maxfluxcorr')   ; READ(arg,*) MaxFluxCorr
       CASE('-sourcemask')    ; READ(arg,'(a)') MaskFile
       CASE('-hsncube')       ; READ(arg,'(a)') hSNCube
       CASE('-sourcecube')    ; READ(arg,'(a)') sourcecube
       CASE('-lcheck')        ; READ(arg,'(a)') linecheck; i_linecheck=1
       CASE('-l2sky')         ; READ(arg,*) source_line_to_sky_ratio
       CASE('-l2con')         ; READ(arg,*) source_line_to_cont_ratio
       CASE('-lsig')          ; READ(arg,*) source_line_significance
       CASE('-ncont')         ; READ(arg,*) npixcont
       CASE('-skipobjid')     ; READ(arg,'(a)') skipobjid_longstring
       CASE default
          print *, "command line argument ",TRIM(opt), " not recognized!"
          STOP
       END SELECT
    END DO
    
    !..performs few checks
    IF(TRIM(inpname)=="??") STOP "please provide datacube name with the option -cube!"
    IF(TRIM(outname)=="??") THEN
       is=INDEX(TRIM(inpname),".fits")
       IF(is==0) is=LEN_TRIM(inpname)+1
       outname=TRIM(inpname(1:is-1))//"_SHARP.fits"
    ENDIF
    
    IF(i_linecheck==0) THEN !..linecheck option was not provided:
       !..thus decide if we use linecheck or not depending on hSNCube
       IF(TRIM(hSNCube)=="??") THEN
          linecheck=.false.
          print *, "NB: setting linecheck=.false."
       ELSE
          linecheck=.true.
          print *, "NB: setting linecheck=.true."
       END IF
    END IF

    IF(TRIM(skipobjid_longstring)/="??") THEN

       IF(TRIM(MaskFile)=="??") STOP "You must provide a sourcemaks with objects id (-sourcemask)"

       !..split longstring and assign values
       !..finds how many entry there are in the string
       DO i=1,1000
          ierr=0
          READ(skipobjid_longstring,*,iostat=ierr) testval(1:i)
          IF(ierr/=0) THEN
             nSkipObjId=i-1
             EXIT
          END IF
       END DO
       IF(i>=1000) STOP "problem reading -skipobjid option!"
       ALLOCATE(SkipObjId(nSkipObjId))
       !..read entries
       READ(skipobjid_longstring,*) SkipObjId(:)

       print *, "hsncube will not be used for spaxel corresponding to these Object Ids: ", SkipObjId(1:nSkipObjId)
       print *, "in the sourcemask file:", TRIM(MaskFile)
       print *, " "

    END IF
       

  END SUBROUTINE ReadParam

!-------------------------------------------------------------

  SUBROUTINE SplitCube
    
    IMPLICIT NONE
    REAL(kind=4), ALLOCATABLE :: sky(:),sky_or(:),sky_sigma(:), local_sourcemask(:,:)
    REAL(kind=4) :: skystep, binsize, this_sky_val
    INTEGER :: this_bin, j, zbin, cbin_count
    INTEGER, PARAMETER :: spsize=2


    print *, "extracting average sky spectrum..."
    ALLOCATE(sky(SIZE(Cube,DIM=3)),sky_or(SIZE(Cube,DIM=3)), sky_sigma(SIZE(Cube,DIM=3)))

    !..produce a local object mask
    IF(TRIM(MaskFile)=="??") THEN
       ALLOCATE(local_sourcemask(DimX,DimY))
       DO y=1,DimY
          DO x=1,DimX
             IF(ALL(Cube(x,y,:)==UNDEF)) THEN
                local_sourcemask(x,y)=UNDEF
             ELSE
                local_sourcemask(x,y)=Mean(PACK(Cube(x,y,:),MASK=Cube(x,y,:)/=UNDEF))
             END IF
          END DO
       END DO
       CALL SigmaClip(PACK(local_sourcemask,MASK=local_sourcemask/=UNDEF),clipmean,clipmedian,clipsigma)
       WHERE(local_sourcemask>sclipWL*clipsigma+clipmean.and.local_sourcemask/=UNDEF) 
          local_sourcemask=1.0
       ELSEWHERE(local_sourcemask/=UNDEF)
          local_sourcemask=0.0
       END WHERE
    ELSE
       ALLOCATE(MaskCube(DimX,DimY,1),local_sourcemask(DimX,DImY))
       ALLOCATE(full_sourcemask(DimX,DimY))
       CALL ReadLocalCube(MaskFile, MaskCube)
       print *, "Min Max sourcemask=", MINVAL(MaskCube), MAXVAL(MaskCube)
       WHERE(MaskCube(:,:,1)/=0) 
          local_sourcemask=1.0
          full_sourcemask=1.0
       ELSEWHERE
          local_sourcemask=0.0
          full_sourcemask=0.0
       END WHERE
       !..deallocate MaskCube here if not needed later
       IF(nSkipObjId==0) DEALLOCATE(MaskCube)
    END IF

   !..produce sky spectrum    
    DO j=1,SIZE(sky)
       IF(ANY(Cube(:,:,j)/=UNDEF)) THEN
          sky_or(j)=Mean(PACK(Cube(:,:,j),MASK=Cube(:,:,j)/=UNDEF.and.local_sourcemask/=UNDEF))
          sky_sigma(j)=StdDev(PACK(Cube(:,:,j),MASK=Cube(:,:,j)/=UNDEF.and.local_sourcemask/=UNDEF))
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

!..build zbins array looping over the sky array 
!..and selecting the minloc of sky between skybin_min and 
!..skybin_max sizes
    print *, "removing sky continuum and producing wavelength splitting array..." 
    WRITE(*,'(a,$)') "%"
    this_bin=0
    j=1
    zbin=0
    old_percent=0
    avg_var_over_sky=0.
    cbin_count=0
    !..cross the sky to find skylines for wl array and subtract sky continuum otherwise
    mainloop: DO WHILE(j<=DimZ)

       IF(sky(j)>skysel*avg_sky) THEN !..we are entering a sky-line region

          zbin=zbin+1
          IF(zbin>1000) STOP "increase splitting array size!"

          wl(zbin)%zmin=j

          skyline_loop:DO WHILE(sky(j)>skysel*avg_sky) !..move until we exit the region or until we find a proper minimum in the unsmoothed sky
             j=j+1
             IF(sky_or(j-1)<skydeblend*avg_sky.and.sky_or(j)<skydeblend*avg_sky.and.& !..this is to avoid cutting a line at a high-flux edge
                  sky_or(MAX(1,j-1))>sky_or(j).and.sky_or(MIN(SIZE(sky_or),j+1))>sky_or(j).and.&  !..this is the minimum check
                  sky(MIN(SIZE(sky),j+1))>skysel*avg_sky) & !..this is to avoid cutting away one pixel only 
                  EXIT skyline_loop
             IF(j==DimZ) THEN
                wl(zbin)%zmax=j
                EXIT mainloop
             END IF
          END DO skyline_loop
          j=j-1

          wl(zbin)%zmax=j

       ELSE

          IF(ANY(Cube(:,:,j)/=UNDEF.and.local_sourcemask/=UNDEF).and.subtract_sky_continuum) THEN
             !..subtract the sky continuum
             SELECT CASE(TRIM(skytype))
             CASE("mean")
                this_sky_val=Mean(PACK(Cube(:,:,j),MASK=Cube(:,:,j)/=UNDEF.and.local_sourcemask==0))
             CASE("median") 
                this_sky_val=Median(PACK(Cube(:,:,j),MASK=Cube(:,:,j)/=UNDEF.and.local_sourcemask==0))
             CASE("avgsigclip") 
                CALL SigmaClip(PACK(Cube(:,:,j),MASK=Cube(:,:,j)/=UNDEF.and.local_sourcemask==0), this_sky_val, clipmedian, clipsigma)
             CASE("medsigclip")
                CALL SigmaClip(PACK(Cube(:,:,j),MASK=Cube(:,:,j)/=UNDEF.and.local_sourcemask==0), clipmean, this_sky_val, clipsigma)
             CASE("nosub")
                this_sky_val=0.
             CASE default
                STOP "selected skytype is not available!"
             END SELECT
             avg_var_over_sky=avg_var_over_sky+(clipsigma**2)/this_sky_val
             cbin_count=cbin_count+1
             WHERE(Cube(:,:,j)/=UNDEF) Cube(:,:,j)=Cube(:,:,j)-this_sky_val
          END IF


       END IF

       percent=INT(REAL(j)/DimZ*10)
       IF(percent/=old_percent) WRITE(*,'(i4,$)') percent*10 
       old_percent=percent

       j=j+1

    END DO mainloop

    IF(cbin_count/=0) THEN
       avg_var_over_sky=avg_var_over_sky/cbin_count
       print *, " "
       print *, "avg_var_over_sky=",avg_var_over_sky," calculated using n=",cbin_count, " wavelenght layers"
    ELSE
       avg_var_over_sky=1
    END IF

    print *, " "
    print *, " "
    print *, "sky lines to be processed:"
    print *, " "
    print *, "zbin  zmin     zmax    sky@zmin   sky@zmax   maxsky"
    DO j=1,zbin
       WRITE(*,'(i3,1x,i6,1x,i6,1x,3(f9.3,1x))') j,wl(j)%zmin, wl(j)%zmax, sky_or(wl(j)%zmin), sky_or(wl(j)%zmax), MAXVAL(sky_or(wl(j)%zmin:wl(j)%zmax)) 
    END DO
    print *, " "

    wl_size=zbin

    DEALLOCATE(sky,sky_or,sky_sigma)

    ALLOCATE(SkySub0(DimX,DimY,DimZ))

    IF(TRIM(hSNCube)=="??") THEN

       !..remove 0th order sky lines
       print *, "removing 0th order sky lines..."
       SkySub0=Cube
       
       DO line=1,wl_size

          !..select initial and final spectral pixel
          z1=wl(line)%zmin
          z2=wl(line)%zmax
          
          DO z=z1,z2
             SELECT CASE(TRIM(skytype))
             CASE("mean")
                this_sky_val=Mean(PACK(Cube(:,:,z),MASK=Cube(:,:,z)/=UNDEF.and.local_sourcemask==0))
             CASE("median") 
                this_sky_val=Median(PACK(Cube(:,:,z),MASK=Cube(:,:,z)/=UNDEF.and.local_sourcemask==0))
             CASE("avgsigclip") 
                CALL SigmaClip(PACK(Cube(:,:,z),MASK=Cube(:,:,z)/=UNDEF.and.local_sourcemask==0), this_sky_val, clipmedian, clipsigma,ClipVal=[-clipval,clipval])
             CASE("medsigclip")
                CALL SigmaClip(PACK(Cube(:,:,z),MASK=Cube(:,:,z)/=UNDEF.and.local_sourcemask==0), clipmean, this_sky_val, clipsigma,ClipVal=[-clipval,clipval])
             CASE("nosub")
                this_sky_val=0.
             CASE default
                STOP "selected skytype is not available!"
             END SELECT
             WHERE(SkySub0(:,:,z)/=UNDEF) SkySub0(:,:,z)=SkySub0(:,:,z)-this_sky_val
             !print *, z,this_sky_val
          END DO

       END DO
       print *, "done"
       print *, " "
       
    ELSE !..read the user provided higher SN and sky-subtracted cube for continuum source removal

       CALL ReadLocalCube(hSNCube,LocalCube=SkySub0)

    END IF

    IF(nSkipObjId>0.and.TRIM(hSNCube)/="??") THEN !..replace values of hSNCube with 0th order sky-sub from individual cube

       ALLOCATE(sky(SIZE(Cube,Dim=3)))
       sky=0.

       !..first get sky
       DO line=1,wl_size

          !..select initial and final spectral pixel
          z1=wl(line)%zmin
          z2=wl(line)%zmax
          
          DO z=z1,z2
             SELECT CASE(TRIM(skytype))
             CASE("mean")
                this_sky_val=Mean(PACK(Cube(:,:,z),MASK=Cube(:,:,z)/=UNDEF.and.local_sourcemask==0))
             CASE("median") 
                this_sky_val=Median(PACK(Cube(:,:,z),MASK=Cube(:,:,z)/=UNDEF.and.local_sourcemask==0))
             CASE("avgsigclip") 
                CALL SigmaClip(PACK(Cube(:,:,z),MASK=Cube(:,:,z)/=UNDEF.and.local_sourcemask==0), this_sky_val, clipmedian, clipsigma,ClipVal=[-clipval,clipval])
             CASE("medsigclip")
                CALL SigmaClip(PACK(Cube(:,:,z),MASK=Cube(:,:,z)/=UNDEF.and.local_sourcemask==0), clipmean, this_sky_val, clipsigma,ClipVal=[-clipval,clipval])
             CASE("nosub")
                this_sky_val=0.
             CASE default
                STOP "selected skytype is not available!"
             END SELECT
             sky(z)=this_sky_val
           END DO

       END DO

       !..replace SkySub0 where necessary
       DO y=1,SIZE(Cube,DIM=2)
          DO x=1,SIZE(Cube,DIM=1)
             IF(ANY(SkipObjId(:)==MaskCube(x,y,1))) THEN
                WHERE(Cube(x,y,:)/=UNDEF) 
                   SkySub0(x,y,:)=Cube(x,y,:)-sky(:)
                ELSEWHERE
                   SkySub0(x,y,:)=UNDEF
                END WHERE
             END IF
          END DO
       END DO
 
       DEALLOCATE(sky)
      
    END IF


   DEALLOCATE(local_sourcemask)

  END SUBROUTINE SplitCube

!------------------------------------------------------------
!
! Finds and remove continuum sources using 0th order sky subtraction
! or a higher SN cube
! 

SUBROUTINE RemoveContSources

  IMPLICIT NONE
  REAL, ALLOCATABLE :: BBimage(:,:), SkySub0_sigma(:)
  INTEGER :: n2skip, ntot, z1, z2, z1_, z2_, percent, old_percent
  REAL    :: this_cont, this_sigma, this_median, this_mean

     
  ALLOCATE(SourceContinuum(DimX,DimY,DimZ))
  SourceContinuum=0.

  !..create a global sourcemask if a Maskfile is not provided
  !..(in the latter case, full_sourcemask is created in SplitCube)
  IF(TRIM(MaskFile)=="??") THEN

     !..produce a broad band image for continuum source detection 
     ALLOCATE(BBimage(DimX,DimY))
     DO y=1,DimY
        DO x=1,DimX
           IF(ANY(SkySub0(x,y,:)/=UNDEF)) THEN
              BBimage(x,y)=Mean(PACK(SkySub0(x,y,:),MASK=SkySub0(x,y,:)/=UNDEF))
           ELSE
              BBimage(x,y)=UNDEF
           END IF
        END DO
     END DO

     CALL SigmaClip(PACK(BBimage,MASK=BBimage/=UNDEF),clipmean,clipmedian,clipsigma)

     ALLOCATE(full_sourcemask(DimX,DimY))
     full_sourcemask(:,:)=UNDEF
     WHERE(BBimage(:,:)/=UNDEF.and.BBimage(:,:)>clipmean+sclipWL*clipsigma) 
        full_sourcemask(:,:)=1.0
     ELSEWHERE(BBimage(:,:)/=UNDEF)
        full_sourcemask(:,:)=0.0
     END WHERE

     DEALLOCATE(BBimage)
     
  ENDIF
     
  IF(TRIM(sourcecube)=="??") THEN

     WRITE(*,'(a,$)') "%"

     !..loop over sky lines and estimate continuum
     n2skip=0
     ntot=0
     old_percent=0
     DO line=1,wl_size 
        
        percent=INT(REAL(line)/wl_size*10)
        IF(percent/=old_percent) WRITE(*,'(i4,$)') percent*10 
        old_percent=percent

        z1=wl(line)%zmin
        z2=wl(line)%zmax
        
        !..get an estimate of the residual sky variation from the sky subtracted cube
        !..this is needed later to check for line features
        IF(ALLOCATED(SkySub0_sigma)) DEALLOCATE(SkySub0_sigma)
        ALLOCATE(SkySub0_sigma(z1:z2))
        DO z=z1,z2
           CALL SigmaClip(PACK(SkySub0(:,:,z),MASK=SkySub0(:,:,z)/=UNDEF.and.full_sourcemask==0.),clipmean, clipmedian, SkySub0_sigma(z))
        END DO

        !..loop over spatial positions
        DO y=1,DimY
           spatial_loop: DO x=1,DimX

              ntot=ntot+1

              IF(full_sourcemask(x,y)/=1) CYCLE

              !..loop over spectral pixels around npixcont to estimate continuum
              DO z=z1,z2

                                  
                 z1_=MAX(1,z-npixcont)
                 z2_=MIN(DimZ,z+npixcont)
                 
                 CALL SigmaClip(PACK(SkySub0(x,y,z1_:z2_),MASK=SkySub0(x,y,z1_:z2_)/=UNDEF),this_mean,this_median,this_sigma,ClipVal=[-5.,5.])
              
                 IF(this_median>0) THEN
                    SourceContinuum(x,y,z)=this_median
                 ELSE
                    SourceContinuum(x,y,z)=0.
                 END IF

              END DO

              !..check for line feature presence in this spectrum comparing to the rest of the layers, if requested
              !..NB: line needs to be significant with respect to the continuum (as defined by source_line_to_cont_ratio)
              !..    and with respect to the sky line (as defined by source_line_to_sky_ratio)
              !..in this case, set SourceContinuum to UNDEF, so that spectra will be skipped

              IF(nSkipObjId>0) THEN !..check if we are skipping any object 
                 IF(ANY(SkipObjId(:)==MaskCube(x,y,1))) CYCLE spatial_loop
              END IF

              IF(linecheck) THEN
                 IF(ANY(SkySub0(x,y,z1:z2)-SourceContinuum(x,y,z1:z2)>source_line_significance*SkySub0_sigma(z1:z2).and.&
                      SkySub0(x,y,z1:z2)>source_line_to_cont_ratio*SourceContinuum(x,y,z1:z2).and.&
                      SkySub0(x,y,z1:z2)-SourceContinuum(x,y,z1:z2)>source_line_to_sky_ratio*Cube(x,y,z1:z2).and.&
                      Cube(x,y,z1:z2)>0.)) THEN
                    SourceContinuum(x,y,z1:z2)=UNDEF
                    n2skip=n2skip+1
                 ENDIF
              END IF

           END DO spatial_loop
        END DO

     END DO

     print *, " "

     IF(n2skip/=0) print *, "# of spectra that will be skipped because of a possible line feature = ", n2skip, " over ", ntot


  ELSE  !..read the SourceContinuum from the cube provided

     CALL ReadLocalCube(InpFile=sourcecube,LocalCube=SourceContinuum)

  END IF

  !CALL WriteLocalCube(SourceContinuum,"SourceCont.fits")
  !STOP


END SUBROUTINE RemoveContSources



!---------------------------------------------------------------
!
! Performs an overall shift of the flux (without propagating the variance)
!
  SUBROUTINE SpcShift(x, y, zshift)

    IMPLICIT NONE
    INTEGER(kind=4), INTENT(IN) :: x,y
    REAL(kind=4), INTENT(OUT)   :: zshift
    INTEGER(kind=4) :: i,j,k, this_bin
    REAL(kind=4) :: this_sky(z1:z2), sky_or(z1:z2), tot_var(zbins), this_shift, cfact, this_sky_s(z1:z2), this_sky_u(z1:z2)
    REAL(kind=4) :: this_var(z1:z2), var_or(z1:z2), skynorm

    sky_or=Cube(x,y,z1:z2)
    var_or=Var(x,y,z1:z2)
    this_var=1.e30
    zshift=0.
    this_sky=sky_or
    this_var=var_or
    this_sky_s=Cube_sources(x,y,z1:z2)
    IF(shiftcheck) this_sky_u=UnitCube(x,y,z1:z2)

    skynorm=1./SUM(sky_or)

    DO i=1,zbins
    
       this_shift=-maxshift+(i-1)*zstep

       IF(this_shift<0) THEN
          this_sky(z2)=sky_or(z2)-sky_or(z2)*abs(this_shift)
          DO k=z2-1,z1+1,-1
             this_sky(k)=sky_or(k)-sky_or(k)*abs(this_shift)+sky_or(k+1)*abs(this_shift)
          END DO
          this_sky(z1)=sky_or(z1)+sky_or(z1+1)*abs(this_shift)
       ELSE
          this_sky(z1)=sky_or(z1)-sky_or(z1)*this_shift
          DO k=z1+1,z2-1,1
             this_sky(k)=sky_or(k)-sky_or(k)*this_shift+sky_or(k-1)*this_shift
          END DO
          this_sky(z2)=sky_or(z2)+sky_or(z2-1)*this_shift
       END IF

       tot_var(i)=SUM(((this_sky*skynorm-skyref*refnorm)**2)*skyweight)

    END DO
    
    this_bin=MINLOC(tot_var,DIM=1)
    IF(this_bin==1.or.this_bin==zbins) THEN
       IF(.not.applyfineshift) THEN
          notconv=notconv+1
          ShiftCube(x,y,z1:z2)=UNDEF
       END IF
       RETURN
    END IF

    zshift=-maxshift+(this_bin-1)*zstep
    this_shift=zshift

    !..replace sky using shift calculated on continuum subtracted cube
    !..and apply to continuum-unsubtracted cube as well if there is a source at this location

    IF(this_shift<0) THEN
       this_sky(z2)=sky_or(z2)-sky_or(z2)*abs(this_shift)
       DO k=z2-1,z1+1,-1
          this_sky(k)=sky_or(k)-sky_or(k)*abs(this_shift)+sky_or(k+1)*abs(this_shift)
       END DO
       this_sky(z1)=sky_or(z1)+sky_or(z1+1)*abs(this_shift)
    ELSE
       this_sky(z1)=sky_or(z1)-sky_or(z1)*this_shift
       DO k=z1+1,z2-1,1
          this_sky(k)=sky_or(k)-sky_or(k)*this_shift+sky_or(k-1)*this_shift
       END DO
       this_sky(z2)=sky_or(z2)+sky_or(z2-1)*this_shift
    END IF

    IF(sourcemask(x,y)==1.) THEN
       IF(this_shift<0) THEN
          this_sky_s(z2)=Cube_sources(x,y,z2)-Cube_sources(x,y,z2)*abs(this_shift)*SkyRatio(x,y,z2)
          DO k=z2-1,z1+1,-1
             this_sky_s(k)=Cube_sources(x,y,k)-Cube_sources(x,y,k)*abs(this_shift)*SkyRatio(x,y,k)+Cube_sources(x,y,k+1)*abs(this_shift)*SkyRatio(x,y,k+1)
          END DO
          this_sky_s(z1)=Cube_sources(x,y,z1)+Cube_sources(x,y,z1+1)*abs(this_shift)*SkyRatio(x,y,z1+1)
       ELSE
          this_sky_s(z1)=Cube_sources(x,y,z1)-Cube_sources(x,y,z1)*this_shift*SkyRatio(x,y,z1)
          DO k=z1+1,z2-1,1
             this_sky_s(k)=Cube_sources(x,y,k)-Cube_sources(x,y,k)*this_shift*SkyRatio(x,y,k)+Cube_sources(x,y,k-1)*this_shift*SkyRatio(x,y,k-1)
          END DO
          this_sky_s(z2)=Cube_sources(x,y,z2)+Cube_sources(x,y,z2-1)*this_shift*SkyRatio(x,y,z2-1)
       END IF
    END IF

    IF(shiftcheck) THEN
       IF(this_shift<0) THEN
          this_sky_u(z2)=UnitCube(x,y,z2)-UnitCube(x,y,z2)*abs(this_shift)
          DO k=z2-1,z1+1,-1
             this_sky_u(k)=UnitCube(x,y,k)-UnitCube(x,y,k)*abs(this_shift)+UnitCube(x,y,k+1)*abs(this_shift)
          END DO
          this_sky_u(z1)=UnitCube(x,y,z1)+UnitCube(x,y,z1+1)*abs(this_shift)
       ELSE
          this_sky_u(z1)=UnitCube(x,y,z1)-UnitCube(x,y,z1)*this_shift
          DO k=z1+1,z2-1,1
             this_sky_u(k)=UnitCube(x,y,k)-UnitCube(x,y,k)*this_shift+UnitCube(x,y,k-1)*this_shift
          END DO
          this_sky_u(z2)=UnitCube(x,y,z2)+UnitCube(x,y,z2-1)*this_shift
       END IF
    END IF
    
    Cube(x,y,z1:z2)=this_sky(z1:z2)
    !Var(x,y,z1+INT(maxshift)+1:z2-INT(maxshift)-1)=this_var(z1+INT(maxshift)+1:z2-INT(maxshift)-1)
    IF(sourcemask(x,y)==1) Cube_sources(x,y,z1:z2)=this_sky_s(z1:z2)
    IF(shiftcheck) UnitCube(x,y,z1:z2)=this_sky_u(z1:z2)
    IF(zshift<0) THEN
       WHERE(ShiftCube(x,y,z1+1:z2)/=UNDEF) ShiftCube(x,y,z1+1:z2)=ShiftCube(x,y,z1+1:z2)+zshift
    ELSE
       WHERE(ShiftCube(x,y,z1:z2-1)/=UNDEF) ShiftCube(x,y,z1:z2-1)=ShiftCube(x,y,z1:z2-1)+zshift
    END IF

    IF(ANY(abs(ShiftCube(x,y,z1:z2))>maxtotalshift)) THEN
       notconv=notconv+1
       ShiftCube(x,y,z1:z2)=UNDEF
    ELSE
       didconv=didconv+1
    END IF


  END SUBROUTINE SpcShift

!------------------------------------------------------  
!
! Performs a pixel by pixel flux shift using groups of 3 pixel as a shift-base.
! This routine propagate also the variance. Shifting is performed starting from
! the pixel with the largest variance and continues until minimum variance is 
! found for the overall spectrum (i.e., when no more shift are requested).
! Shifting is applied to the SourceCube and UnitCube as well.
!
SUBROUTINE SpcFineShift(x, y, zmin, zmax)

    IMPLICIT NONE
    INTEGER(kind=4), INTENT(IN) :: x, y, zmin, zmax
    INTEGER(kind=4) :: i,j,k, this_bin, this_it, maxiter, step
    REAL(kind=4) :: this_sky(zmin:zmax), sky_or(zmin:zmax), tot_var, this_shift, cfact, &
         varmin, savesky(zmin:zmax), skyref_norm(zmin:zmax), skynorm, saveshift(zmin:zmax), convshift, &
         this_sky_s(zmin:zmax), sky_or_s(zmin:zmax), savesky_s(zmin:zmax), mask_var(zmin:zmax), &
         this_sky_u(zmin:zmax), sky_or_u(zmin:zmax), savesky_u(zmin:zmax), this_var(zmin:zmax), var_or(zmin:zmax), savevar(zmin:zmax)
    REAL(kind=4),PARAMETER :: UNDEF_SHIFT=999.0
    LOGICAL :: onsource


    onsource=(sourcemask(x,y)==1.0)

    maxiter=15*(zmax-zmin+1)
    convshift=2*zfinestep

    sky_or(zmin:zmax)=Cube(x,y,zmin:zmax)
    var_or(zmin:zmax)=Var(x,y,zmin:zmax)
    skynorm=1./SUM(sky_or)

    tot_var=1.e30
    varmin=1.e30
    zshift=0.
    this_sky=sky_or
    this_var=var_or

    saveshift=UNDEF_SHIFT

    IF(onsource) THEN
       sky_or_s(zmin:zmax)=Cube_sources(x,y,zmin:zmax)
       this_sky_s=sky_or_s
    END IF
    IF(shiftcheck) THEN
       sky_or_u(zmin:zmax)=UnitCube(x,y,zmin:zmax)
       this_sky_u=sky_or_u
    END IF

!---------------------
!    skynorm=1.
!    refnorm=1.
!--------------------


    iterloop:DO this_it=1,maxiter

       !..check convergence
       IF(COUNT(abs(saveshift(zmin:zmax))>convshift.and.&
            .not.skipmask(x,y,zmin:zmax))<1) THEN !..convergence reached
          WHERE(ShiftCube(x,y,zmin:zmax)==UNDEF) ShiftCube(x,y,zmin:zmax)=0.
          EXIT iterloop
       END IF


       !..find pixel with largest variance
       k=MAXLOC(((skyref(zmin:zmax)*refnorm-sky_or(zmin:zmax)*skynorm)**2)*skyweight(zmin:zmax),&
            MASK=abs(saveshift(zmin:zmax))>convshift &
            .and..not.skipmask(x,y,zmin:zmax) &
            ,DIM=1)+zmin-1

       !..apply shift

       varmin=1.e30
       this_sky=sky_or
       this_var=var_or
       IF(onsource) this_sky_s=sky_or_s
       IF(shiftcheck) this_sky_u=sky_or_u

       !..start shift iteration to find "minimum-variance-shift"
       fineloop: DO i=1,zfinebins

          this_shift=-maxfineshift+(i-1)*zfinestep

          IF(this_shift>=0) THEN
             step=1
          ELSE
             step=-1
          END IF

          !..shift sky
          IF(k/=zmin.and.k/=zmax) THEN
             this_sky(k+1)=sky_or(k+1)+step*sky_or(k-(step-1)/2)*abs(this_shift)
             this_sky(k)=sky_or(k)-sky_or(k)*abs(this_shift)+sky_or(k-step)*abs(this_shift)
             this_sky(k-1)=sky_or(k-1)-step*sky_or(k-(step+1)/2)*abs(this_shift)
          ELSEIF(k==zmin) THEN 
             this_sky(k)=sky_or(k)-step*sky_or(k-(step-1)/2)*abs(this_shift)
             this_sky(k+1)=sky_or(k+1)+step*sky_or(k-(step-1)/2)*abs(this_shift)
          ELSE
             this_sky(k)=sky_or(k)+step*sky_or(k-(step+1)/2)*abs(this_shift)
             this_sky(k-1)=sky_or(k-1)-step*sky_or(k-(step+1)/2)*abs(this_shift)         
          END IF

          !..shift variance
          IF(k/=zmin.and.k/=zmax) THEN
             this_var(k+1)=var_or(k+1)+var_or(k-(step-1)/2)*abs(this_shift)**2
             this_var(k)=var_or(k)+var_or(k)*abs(this_shift)**2+var_or(k-step)*abs(this_shift)**2
             this_var(k-1)=var_or(k-1)+var_or(k-(step+1)/2)*abs(this_shift)**2
          ELSEIF(k==zmin) THEN 
             this_var(k)=var_or(k)+var_or(k-(step-1)/2)*abs(this_shift)**2
             this_var(k+1)=var_or(k+1)+var_or(k-(step-1)/2)*abs(this_shift)**2
          ELSE
             this_var(k)=var_or(k)+var_or(k-(step+1)/2)*abs(this_shift)**2
             this_var(k-1)=var_or(k-1)+var_or(k-(step+1)/2)*abs(this_shift)**2         
          END IF
       
          !..shift sky in SourceCube
          IF(onsource) THEN
             IF(k/=zmin.and.k/=zmax) THEN
                this_sky_s(k+1)=sky_or_s(k+1)+step*sky_or_s(k-(step-1)/2)*abs(this_shift)*SkyRatio(x,y,k-(step-1)/2)
                this_sky_s(k)=sky_or_s(k)-sky_or_s(k)*abs(this_shift)*SkyRatio(x,y,k)+sky_or_s(k-step)*abs(this_shift)*SkyRatio(x,y,k-step)
                this_sky_s(k-1)=sky_or_s(k-1)-step*sky_or_s(k-(step+1)/2)*abs(this_shift)*SkyRatio(x,y,k-(step+1)/2)
             ELSEIF(k==z1) THEN 
                this_sky_s(k)=sky_or_s(k)-step*sky_or_s(k-(step-1)/2)*abs(this_shift)*SkyRatio(x,y,k-(step-1)/2)
                this_sky_s(k+1)=sky_or_s(k+1)+step*sky_or_s(k-(step-1)/2)*abs(this_shift)*SkyRatio(x,y,k-(step-1)/2)
             ELSE
                this_sky_s(k)=sky_or_s(k)+step*sky_or_s(k-(step+1)/2)*abs(this_shift)*SkyRatio(x,y,k-(step+1)/2)
                this_sky_s(k-1)=sky_or_s(k-1)-step*sky_or_s(k-(step+1)/2)*abs(this_shift)*SkyRatio(x,y,k-(step+1)/2)         
             END IF
          END IF

          !..shift sky in UnitCube
          IF(shiftcheck) THEN
            IF(k/=zmin.and.k/=zmax) THEN
                this_sky_u(k+1)=sky_or_u(k+1)+step*sky_or_u(k-(step-1)/2)*abs(this_shift)
                this_sky_u(k)=sky_or_u(k)-sky_or_u(k)*abs(this_shift)+sky_or_u(k-step)*abs(this_shift)
                this_sky_u(k-1)=sky_or_u(k-1)-step*sky_or_u(k-(step+1)/2)*abs(this_shift)
             ELSEIF(k==zmin) THEN 
                this_sky_u(k)=sky_or_u(k)-step*sky_or_u(k-(step-1)/2)*abs(this_shift)
                this_sky_u(k+1)=sky_or_u(k+1)+step*sky_or_u(k-(step-1)/2)*abs(this_shift)
             ELSE
                this_sky_u(k)=sky_or_u(k)+step*sky_or_u(k-(step+1)/2)*abs(this_shift)
                this_sky_u(k-1)=sky_or_u(k-1)-step*sky_or_u(k-(step+1)/2)*abs(this_shift)         
             END IF
          END IF

          tot_var=SUM(((this_sky*skynorm-skyref*refnorm)**2)*skyweight)

          IF(tot_var<varmin) THEN
             varmin=tot_var
             savesky=this_sky
             savevar=this_var
             IF(onsource) savesky_s=this_sky_s
             IF(shiftcheck) savesky_u=this_sky_u
             saveshift(k)=this_shift
             this_bin=i
          END IF

       END DO fineloop

       !..check if a minimum has been found
       IF(this_bin==1.or.this_bin==zfinebins) THEN
          ShiftCube(x,y,zmin:zmax)=UNDEF
          notconv=notconv+1 
          RETURN 
       END IF
       
       sky_or=savesky
       var_or=savevar
       IF(onsource) sky_or_s=savesky_s
       IF(shiftcheck) sky_or_u=savesky_u

       IF(ALL(abs(saveshift)<=convshift)) EXIT iterloop

       !..add performed shift to ShiftCube
       IF(saveshift(k)/=UNDEF_SHIFT) THEN
          IF(saveshift(k)<0) THEN
             WHERE(ShiftCube(x,y,MAX(k,zmin+1):MIN(zmax,k+1))/=UNDEF) ShiftCube(x,y,MAX(k,zmin+1):MIN(zmax,k+1))=ShiftCube(x,y,MAX(k,zmin+1):MIN(zmax,k+1))+saveshift(k)
          ELSE
             WHERE(ShiftCube(x,y,MAX(zmin,k-1):MIN(k,zmax-1))/=UNDEF) ShiftCube(x,y,MAX(zmin,k-1):MIN(k,zmax-1))=ShiftCube(x,y,MAX(zmin,k-1):MIN(k,zmax-1))+saveshift(k)
          END IF
       END IF

    END DO iterloop

    !..check if we have shifted too much
    IF(ANY(abs(ShiftCube(x,y,zmin:zmax))>maxtotalshift)) THEN
       notconv=notconv+1
       ShiftCube(x,y,zmin:zmax)=UNDEF
    END IF

    !..check if we have excedeed max number of iterations
    IF(this_it>=maxiter) THEN
       ShiftCube(x,y,zmin:zmax)=UNDEF
       notconv=notconv+1
    ELSE
       didconv=didconv+1
       !..replace sky
       Cube(x,y,zmin:zmax)=sky_or(zmin:zmax)
       Var(x,y,zmin:zmax)=var_or(zmin:zmax)
       IF(onsource) Cube_sources(x,y,zmin:zmax)=sky_or_s(zmin:zmax)
       IF(shiftcheck) UnitCube(x,y,zmin:zmax)=sky_or_u(zmin:zmax)
    END IF


  END SUBROUTINE SpcFineShift

END PROGRAM CubeSharp
  

  

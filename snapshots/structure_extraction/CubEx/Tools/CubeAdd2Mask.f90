PROGRAM CubeAdd2Mask

 USE StatLib
 USE CubeLib
 IMPLICIT NONE
 CHARACTER(len=500) :: maskname, datacube, pixtable, mapname, revert, longstring, idfilename 
 INTEGER(kind=4) :: DimX, DimY, DimZ, xbin, ybin, xbin_, ybin_, x, y, i, j, ii, jj, IDn 
 INTEGER(kind=4), PARAMETER :: os_fact=1
 INTEGER :: MUSE_ORIGIN_SHIFT_IFU=6, MUSE_ORIGIN_SHIFT_XSLICE=24
 REAL(kind=4)    :: ifustack(96)
 REAL(kind=4)    :: xfrac, yfrac, pixfrac
 INTEGER(kind=4), ALLOCATABLE :: origin(:,:), IFUmap(:,:), dq(:,:), SliceMap(:,:), slice(:,:), ifu(:,:)
 REAL(kind=4), ALLOCATABLE :: data(:,:), np(:,:), pos(:,:),IFUmap_(:,:,:), SliceMap_(:,:,:), Mask(:,:), IFUSliceMap(:,:)
 REAL(kind=8), ALLOCATABLE :: xpos(:,:), ypos(:,:), dcos_ypos(:,:), dsin_ypos(:,:), dcos_xpos(:,:)
 LOGICAL :: stackMask(24,4), saveold, interactive, ex, allslices(24,4)
 INTEGER :: unit,status, blocksize, rwstatus, naxes(2), nfound, group, ierr, stack, nmask

 CALL ReadCommandLine

 !..read previous SliceEdgeMask
 VERBOSITY=1
 print *, "reading: ", TRIM(maskname)
 CALL ReadCube(maskname)
 DimX=SIZE(Cube, DIM=1)
 DimY=SIZE(Cube, DIM=2)
 ALLOCATE(Mask(DimX,DimY))
 Mask(:,:)=Cube(:,:,1)
 DEALLOCATE(Cube)

 IF(TRIM(mapname)/="??") THEN !..use IFUSliceMap

    CALL ReadMap

 ELSE

    !..read data and WCS info from datacube
    VERBOSITY=2
    CALL ReadCube(datacube,n_ext=1)
    DimZ=SIZE(Cube, DIM=3)
    VERBOSITY=1
    print *, "MinMax Cube=", MINVAL(Cube,MASK=Cube/=UNDEF), MAXVAL(Cube,MASK=Cube/=UNDEF)

    !..read PIXTABLE 
    CALL ReadPixTable
    
    !..produce IFU and SliceMap
    CALL MapPixTable

 END IF


 IF(interactive) THEN

    CALL RunInteractiveMode

 ELSE

    IF(TRIM(longstring)=="all") THEN !..mask all IFU edges
       
       stackmask=.true.
       
    ELSEIF( TRIM(longstring) /= "??") THEN

       !..check if the entry is a file
       INQUIRE(file=TRIM(longstring),Exist=ex)
       IF(ex) THEN
          OPEN(11,file=longstring,action="read")
          READ(11,'(a)') longstring
          CLOSE(11)
       END IF

       
       !..split the longstring and assign IFU.stack masks
       stackmask=.false. !..default
       allslices=.false.
       !..finds how many entry there are in the string
       DO i=1,96  !.96 is the maximum!
          ierr=0
          READ(longstring,*,iostat=ierr) ifustack(1:i)
          IF(ierr/=0) THEN
             nmask=i-1
             EXIT
          END IF
       END DO
       !..read entries
       READ(longstring,*) ifustack(1:nmask)
       !..allocate to mask
       DO i=1,nmask
          stack=NINT(10*(abs(ifustack(i))-INT(abs(ifustack(i)))))
          IF(stack==0) THEN !..mask all stacks
             stackmask(INT(abs(ifustack(i))),:)=.true.
             IF(ifustack(i)<0.) allslices(INT(abs(ifustack(i))),:)=.true.
          ELSE !..mask individual stack
             stackmask(INT(abs(ifustack(i))),stack)=.true.
             IF(ifustack(i)<0.) allslices(INT(abs(ifustack(i))),stack)=.true.
          END IF
       END DO
       
    END IF
    
 END IF

 print *, "stackmask table per ifu and stack (muse-like orientation):"
 DO i=24,1,-1
    WRITE(*,'(i2.2,5x,4(L,3x))') i,stackmask(i,1),stackmask(i,2),stackmask(i,3),stackmask(i,4)
 END DO


 !..add selected slice values to the mask
 print *, "updating: ", TRIM(maskname)
 CALL AddToMask

  print *, "done!"

CONTAINS

!----------------------------------------------

 SUBROUTINE ReadCommandLine

   IMPLICIT NONE
   CHARACTER(len=500) :: ExeName, arg, opt, cubebase, cmd
   INTEGER :: i, nmask, ierr, narg, is
   LOGICAL :: ex

   longstring="??"
   maskname="??"
   mapname="??"
   pixtable="??"
   saveold=.true.
   revert="??"
   interactive=.false.
   idfilename="??"   
   IDn=-1

   CALL GetArg(0,ExeName)
   narg=iargc()
   IF(narg==0) THEN
       print *, " "
       WRITE(*,'(2a)')"                 CubeAdd2Mask (part of CubEx package) ",TRIM(version)
       WRITE(*,'(a)') " "
       WRITE(*,'(a)') " Scope: add to a image mask the edge slices of selected ifu (stacks) given a cube and pixtable"
       WRITE(*,'(a)') "        AND/OR add to the mask an object from a Objects_Id CubEx output"
       WRITE(*,'(a)') " "
       WRITE(*,'(a)') " by Sebastiano Cantalupo (cantalupo@phys.ethz.ch)"
       print *, " "
       WRITE(*,'(a)') "minimal usage: CubeAdd2Mask -cube <name> -add2mask <string>"
       WRITE(*,'(a)') "           or: CubeAdd2Mask -revert <maskname>"
       WRITE(*,'(a)') "  options:"
       WRITE(*,'(a)') "  -cube               <string>          : cube name (NO DEFAULT), not needed if -revert <maskname> option is used."
       WRITE(*,'(a)') "  -i                  <bol>             : if .true., runs the interactive mode using a default image name <cubename>.IM.fits (default=.false.)"
       WRITE(*,'(a)') "                                          masked ifu and stacks are saved in a file called <cubename>_maskedIFUStack.list"
       WRITE(*,'(a)') "  -add2mask           <string>          : IFU and stack list, IFU are in the range 1-24, STACKS are in the range .1-.4  Example: IFU 5, stack 2 is '5.2'"
       WRITE(*,'(a)') "                                          if stack number is not provided or .0 is used all stacks for that IFU will be removed (e.g., '5', '5.' or '5.0')"
       WRITE(*,'(a)') "                                          use 'all' to mask all ifu stack edges. List may be also provided in a file, e.g. '-add2mask ifuslice.list' with"
       WRITE(*,'(a)') "                                          the same format (all entries on a single line)."
       WRITE(*,'(a)') "                                          Use a negative number to mask EVERY slice in the IFU or stack, not just the edges, e.g. -5.2 will mask every slice"
       WRITE(*,'(a)') "                                          in the stack 2 of IFU 5. NB: this option is NOT available in interactive mode."
       WRITE(*,'(a)') "  -mapname            <string>          : associated IFUSliceMap, if a name is not provided it looks for a file called <cube_name>_IFUSliceMap.fits"
       WRITE(*,'(a)') "  -maskname           <string>          : associated mask name to modify (DEFAULT=<cube_name>_SliceEdgeMask.fits)"
       WRITE(*,'(a)') "  -pixtable           <string>          : associated pixtable name, this is NEEDED if a ifu map is not provided or not present in the current folder"
       WRITE(*,'(a)') "  -saveold            <bol>             : if .true. (default), saves the old mask as '.save_<original_mask_name>', this is needed for the '-revert' option"
       WRITE(*,'(a)') "  -revert             <name>            : if a name is provided, revert the mask to the previous state (uses the .save_<original_mask_name>) and stops."
       WRITE(*,'(a)') " "
       WRITE(*,'(a)') ' example:  CubeAdd2Mask -cube DATACUBE_FINAL_WCS_FIX_0008.fits -add2mask " 1 2 15.3 15.4 23.1 23.2 22.1 22.2 21.1 21.2 21.3" -idcube DATACUBE_FINAL_WCS_FIX_0008.Objects_Id.fits -id 8'
       STOP
    END IF
    
    
    !..read from command line
    DO i=1,narg,2
       CALL getarg(i,opt)
       CALL getarg(i+1,arg)
       SELECT CASE(TRIM(opt))
       CASE('-pixtable')     ; READ(arg,'(a)') pixtable
       CASE('-cube')         ; READ(arg,'(a)') datacube
       CASE('-maskname')     ; READ(arg,'(a)') maskname
       CASE('-add2mask')     ; READ(arg,'(a)') longstring
       CASE('-mapname')      ; READ(arg,'(a)') mapname
       CASE('-revert')       ; READ(arg,'(a)') revert
       CASE('-saveold')      ; READ(arg,*) saveold
       CASE('-i')            ; READ(arg,*) interactive
       CASE('-idcube')       ; READ(arg,'(a)') idfilename   
       CASE('-id')           ; READ(arg,*) IDn  
       CASE default
          print *, "command line argument ",TRIM(opt), " not recognized!"
          STOP
       END SELECT
    END DO

    !..revert/save previous version
    IF(TRIM(revert)/="??") THEN
       print *, "reverting: ", TRIM(revert), " to the previous version..."
       cmd="cp -f .save_"//TRIM(revert)//"  "//TRIM(revert)
       CALL SYSTEM(cmd)
       print *, "done"
       STOP
    END IF

    IF(TRIM(mapname)=="??") THEN !..check if a IFUSliceMap with a standard name is present in the current folder
       is=INDEX(TRIM(datacube),".fits")
       IF(is==0) is=LEN_TRIM(datacube)+1
       cubebase=TRIM(datacube(1:is-1))
       mapname=TRIM(cubebase)//"_IFUSliceMap.fits"  
       INQUIRE(File=mapname, Exist=ex)
       IF(.not.ex) THEN
          IF(TRIM(pixtable)=="??") THEN
             print *, "standard IFUSliceMap file: ", TRIM(mapname), " not found in current folder"
             print *, "please provide the name of the file with the -mapname option or provide a pixtable name with the -pixtable option."
             STOP
          END IF
       END IF
    END IF


    IF(TRIM(maskname)=="??") THEN !..apply a default name
       is=INDEX(TRIM(datacube),".fits")
       IF(is==0) is=LEN_TRIM(datacube)+1
       cubebase=TRIM(datacube(1:is-1))
       maskname=TRIM(cubebase)//"_SliceEdgeMask.fits"
    END IF

    IF(saveold) THEN
       cmd="cp -f "//TRIM(maskname)//" .save_"//TRIM(maskname)
       CALL SYSTEM(cmd)
    END IF


  END SUBROUTINE ReadCommandLine

!--------------------------------------------

  SUBROUTINE RunInteractiveMode

    IMPLICIT NONE
    CHARACTER(len=500) :: cmd, imagename, cubebase, ifustring, stackstring, string, ifustack_listname
    INTEGER :: is, ierr, this_IFU, this_stack, this_slice
    REAL    :: x, y
    LOGICAL :: ex

    !..creates a ds9 analysis file
    OPEN(11,file=".analysis.ds9",action="write")
    WRITE(11,*) "Print Coordinates"
    WRITE(11,*) "*.fits"
    WRITE(11,*) "bind s"
    WRITE(11,*) 'echo "$x $y" >> .coords.txt | $null'
    CLOSE(11)

    !..remove previous coords.txt file
    cmd="rm -f .coords.txt"
    CALL SYSTEM(cmd)

    !..get a image file name
    is=INDEX(TRIM(datacube),".fits")
    IF(is==0) is=LEN_TRIM(datacube)+1
    cubebase=TRIM(datacube(1:is-1))
    imagename=TRIM(cubebase)//".IM.fits"  
    ifustack_listname=TRIM(cubebase)//"_maskedIFUStacks.list"
    INQUIRE(File=imagename, Exist=ex)
    IF(.not.ex) THEN
       !..create an image with Cube2Im
       cmd="Cube2Im "//TRIM(datacube)
       CALL SYSTEM(cmd)
    END IF

    print *, "type 's' to select a IFU and stack on the ds9 window. Close the ds9 to exit"

    !..run ds9
    cmd='ds9 '//TRIM(imagename)//" -analysis load .analysis.ds9 -scale mode zscale -cmap b -zoom to fit"
    CALL SYSTEM(cmd)

    !..retrieve coordinates and convert to ifu and stacks
    stackmask=.false.
    print *, " recovered ifu and stack list: "
    OPEN(11,file=".coords.txt",action="read")
    OPEN(12,file=ifustack_listname,action="write")
    DO
       READ(11,*,iostat=ierr) x, y
       IF(ierr/=0) EXIT
       this_IFU=IFUMap(NINT(x),NINT(y))
       IF(this_IFU==0.or.this_IFU==UNDEF) STOP "one of the selected IFU/stack is out of range!"
       this_slice=SliceMap(NINT(x),NINT(y))
       IF(this_slice==0.or.this_slice==UNDEF) STOP "one of the selected IFU/stack is out of range!"
       IF(this_slice<=12) THEN
          this_stack=1
       ELSEIF(this_slice<=24) THEN
          this_stack=2
       ELSEIF(this_slice<=36) THEN
          this_stack=3
       ELSE
          this_stack=4
       END IF
       stackmask(this_IFU,this_stack)=.true.
       WRITE(ifustring,"(i2)") this_IFU
       WRITE(stackstring,"(i1)") this_stack
       string=TRIM(ifustring)//"."//TRIM(stackstring)
       WRITE(*,'(a5,1x,$)') TRIM(string) 
       WRITE(12,'(a5,1x,$)') TRIM(string) 
    END DO
    CLOSE(12)
    print *, " "
    print *, "ifu and stack list saved here: ", TRIM(ifustack_listname)
    print *, " "


  END SUBROUTINE RunInteractiveMode




!----------------------------------------------

  SUBROUTINE ReadMap

    IMPLICIT NONE
    

    CALL ReadCube(mapname)
    IF(ALLOCATED(IFUMap)) DEALLOCATE(IFUMap,SliceMap)
    ALLOCATE(IFUMap(DimX,DimY),SliceMap(DimX,DimY))
    IFUMap=NINT(Cube(:,:,1)*0.01)
    SliceMap=NINT(MOD(Cube(:,:,1),100.))

    !print *, "min max IFU:  ", MINVAL(IFUMap,MASK=IFUMap/=0), MAXVAL(IFUMap)
    !print *, "min max Slice:", MINVAL(SliceMap,MASK=SliceMap/=0), MAXVAL(SliceMap)

  END SUBROUTINE ReadMap


!-----------------------------------------------

  SUBROUTINE ReadPixTable

    IMPLICIT NONE
    REAL(kind=4) :: crpix1, crpix2, xc_cube, yc_cube, aRA, aDEC, xoff, yoff
    REAL(kind=8) :: xoff_d, yoff_d, xmin1, ymin1, xmin2, ymin2
    REAL(kind=8) ::  xmin, xmax, ymin, ymax, dp, xcen, ycen, xd, yd
    REAL(KIND=8), PARAMETER :: pi=3.14159265358979d0
    REAL(KIND=8), PARAMETER :: Rad2Deg=180.0d0/pi, Deg2Rad=pi/180.d0
    REAL(kind=8) :: rx, ry, ra_p, dec_p, phi_p, dx, dy, x, y, r_theta, phi, theta, det, CD_inv(4)
    CHARACTER(len=250) :: comment, fname
    LOGICAL :: anyf


    !..reset status
    status=0

    !..get an unused unit
    CALL ftgiou(unit,status)
    IF(status/=0) STOP "problem finding unused unit for fits open"

    print *, "reading: ", TRIM(pixtable)

    !..get declination from primary header
    fname=TRIM(pixtable)//"[0]"
    CALL ftopen(unit,fname,0,blocksize,status)
    IF(status/=0) THEN
       print *, "Unable to read file: ", TRIM(fname)
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

    ALLOCATE(origin(naxes(1),naxes(2)),&
         dq(naxes(1),naxes(2)),&
         IFU(naxes(1),naxes(2)),&
         slice(naxes(1),naxes(2)))
         !xslice(naxes(1),naxes(2)))

    !CALL ReadPixTableVar("data")
    CALL ReadPixTableVar("origin")
    CALL ReadPixTableVar("dq")
    !ALLOCATE(lambda(naxes(1),naxes(2)))
    !CALL ReadPixTableVar("lambda")
    
    CALL ftclos(unit, status)

    !..convert origin to IFU
    print *, "obtaining IFU information..."
    DO j=1,SIZE(origin)
       IFU(1,j)=IAND(ISHFT(origin(1,j),-MUSE_ORIGIN_SHIFT_IFU),X'1f')
       slice(1,j)=IAND(origin(1,j),X'3f')
       !xslice(1,j)=IAND(ISHFT(origin(1,j),-MUSE_ORIGIN_SHIFT_XSLICE),X'7f')
    END DO
    DEALLOCATE(origin)
    !min_xslice=MINVAL(xslice)
    !max_xslice=MAXVAL(xslice)
    !print *, "original min max xslice=",min_xslice,max_xslice

    !..increase image size and resolution by a factor os_fact, if requested
    DimX=DimX*os_fact
    DimY=DimY*os_fact
    xoff=xoff*os_fact
    yoff=yoff*os_fact

  END SUBROUTINE ReadPixTable

!-----------------------------------------

  SUBROUTINE ReadPixTableVar(extname)

    IMPLICIT NONE
    CHARACTER(len=*) :: extname
    INTEGER :: extver, ANY_HDU, status
    LOGICAL :: anyf

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
    !CASE("lambda")
    !   CALL ftgpve(unit,group,1,naxes(2),UNDEF,lambda,anyf,status)
    CASE("origin")
       CALL ftgpvj(unit,group,1,naxes(2),UNDEF,origin,anyf,status)
    CASE("dq")
       CALL ftgpvj(unit,group,1,naxes(2),UNDEF,dq,anyf,status)
    END SELECT
    IF(status/=0) STOP "problem reading current extension"

  END SUBROUTINE ReadPixTableVar

!---------------------------------------

  SUBROUTINE MapPixTable

    IMPLICIT NONE
    REAL(kind=8) ::  xmin, xmax, ymin, ymax, dp, xcen, ycen, xd, yd

    IF(ALLOCATED(image)) DEALLOCATE(image)
    ALLOCATE(image(DimX,DImY),IFUMap_(DimX,DimY,24),SliceMap_(DimX,DimY,48),&
         IFUMap(DimX,DimY),SliceMap(DimX,DimY))

    !..produce temporary IFU, Slice and XSlice 3D maps using ra (xpos), dec (ypos) from pixeltable
    xmin=0
    DO j=1,SIZE(xpos)
       
       !..remove bad pixels
       IF(dq(1,j)/=0) CYCLE
       
       !..select pixels in the right wavelength range
       !IF(lambda(1,j)<this_zmin.or.lambda(1,j)>this_zmax) CYCLE
       
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
          IF(ANY(Cube(ii,jj,:)/=UNDEF)) THEN
             image(i,j)=Mean(PACK(Cube(ii,jj,:),MASK=Cube(ii,jj,:)/=UNDEF))
          ELSE
             image(i,j)=UNDEF
          END IF
       END DO
    END DO
    
    print *, "minmax image=", MINVAL(image,MASK=(image/=UNDEF)), MAXVAL(image)

    !..reconstruct IFU, Slice and XSlice 2D maps
    print *, "producing IFU, Slice and XSlice maps.."
    DO j=1,DimY
       DO i=1,DimX
          IF(ANY(IFUMap_(i,j,:)/=0.)) IFUMap(i,j)=MAXLOC(IFUMap_(i,j,:),DIM=1)
          IF(ANY(SliceMap_(i,j,:)/=0.)) SliceMap(i,j)=MAXLOC(SliceMap_(i,j,:),DIM=1)
       END DO
    END DO
    WHERE(image==UNDEF) 
       IFUMap=0
       SliceMap=0
    END WHERE

  END SUBROUTINE MapPixTable

!----------------------------------------------------------

  SUBROUTINE AddToMask

    IMPLICIT NONE
    INTEGER :: top_slice(4), bottom_slice(4), i_ifu, i_stack, first_slice(4), last_slice(4)

    !..define top and bottom slice for each stack
    top_slice(:)=[10,22,34,46]
    bottom_slice(:)=[3,15,27,39]

    !..define first and last slice for each stack
    first_slice(:)=[1,13,25,37]
    last_slice(:)=[12,24,36,48]

 !   IF(TRIM(longstring) /= "??") THEN

    !..masking slices
    print *, "masking slices..."
    DO i_ifu=1,24
       DO i_stack=1,4
          IF(stackmask(i_ifu,i_stack)) THEN
             IF(allslices(i_ifu,i_stack)) THEN !..mask every slice
                WHERE(IFUMap==i_ifu.and.(SliceMap>=first_slice(i_stack).and.SliceMap<=last_slice(i_stack))) Mask=0.
             ELSE !..mask edge slices only
                WHERE(IFUMap==i_ifu.and.(SliceMap==top_slice(i_stack).or.SliceMap==bottom_slice(i_stack))) Mask=0.
             END IF
          END IF
       END DO
    END DO

 !   END IF

    !..if requested, add object to the mask
    IF(TRIM(idfilename)/= "??") THEN
       CALL ReadCube(idfilename)
       WHERE(Cube(:,:,1)==IDn) Mask=0.0
    END IF

    !..write image
    IF(.not.allocated(image)) allocate(image(Size(Mask,DIM=1),Size(Mask,DIM=2)))
    image=Mask
    CALL WriteImage(maskname)


  END SUBROUTINE AddToMask


END PROGRAM CubeAdd2Mask
 

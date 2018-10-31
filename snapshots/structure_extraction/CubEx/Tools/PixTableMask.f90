PROGRAM PixTableMask

  USE CubeLib
  IMPLICIT NONE
  CHARACTER(len=500) :: pixtable, ifu_list
  LOGICAL :: mask(24,4)

  CALL ReadParameters

  CALL UpdatePixTable

CONTAINS

  SUBROUTINE ReadParameters

    IMPLICIT NONE
    CHARACTER(len=500) :: ExeName, arg, opt, longstring
    INTEGER :: i, nmask, ierr, narg, stack
    REAL :: ifustack(96)

    CALL GetArg(0,ExeName)
    narg=iargc()
    IF(narg==0) THEN
       print *, " "
       WRITE(*,'(2a)')"       PixTableMask (part of CubEx package) ",TRIM(version)
       WRITE(*,'(a)') "  remove from pixtable the edge slices of selected ifu (stacks) "
       WRITE(*,'(a)') " "
       WRITE(*,'(a)') " by Sebastiano Cantalupo (cantalupo@phys.ethz.ch)"
       print *, " "
       WRITE(*,'(a)') "usage: PixTableMask -pixtable <name> -mask <string>"
       WRITE(*,'(a)') "  options:"
       WRITE(*,'(a)') "  -pixtable           <string>          : pixtable name to modify (NO DEFAULT)"
       WRITE(*,'(a)') "  -mask               <string>          : IFU and stack list, IFU are in the range 1-24, STACKS are in the range .1-.4. Example: IFU 5, stack 2 is '5.2'"
       WRITE(*,'(a)') "                                          if stack number is not provided or .0 is used all stacks for that IFU will be removed (e.g., '5.' or '5.0')"
       WRITE(*,'(a)') " "
       WRITE(*,'(a)') ' example: PixTableMask -old pixtable.fits -new pixtable_NEW.fits -mask "1. 2. 3.1 3.2 4. 5. 6.1 6.2 6.3 7.0 8. 14. 21.1" '
       STOP
    END IF
    
    
    !..read from command line
    DO i=1,narg,2
       CALL getarg(i,opt)
       CALL getarg(i+1,arg)
       SELECT CASE(TRIM(opt))
       CASE('-pixtable')     ; READ(arg,'(a)') pixtable
       CASE('-mask')         ; READ(arg,'(a)') longstring
       CASE default
          print *, "command line argument ",TRIM(opt), " not recognized!"
          STOP
       END SELECT
    END DO
    
    !..split the longstring and assign IFU.stack masks
    mask=.false. !..default
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
       stack=NINT(10*(ifustack(i)-INT(ifustack(i))))
       IF(stack==0) THEN !..mask all stacks
          mask(INT(ifustack(i)),:)=.true.
       ELSE !..mask individual stack
          mask(INT(ifustack(i)),stack)=.true.
       END IF
    END DO

    print *, "mask table per ifu (muse-like orientation):"
    DO i=24,1,-1
       WRITE(*,'(i2.2,5x,4(L,3x))') i,mask(i,1),mask(i,2),mask(i,3),mask(i,4)
    END DO


  END SUBROUTINE ReadParameters

!------------------------------------------------------------

  SUBROUTINE UpdatePixTable

    IMPLICIT NONE
    INTEGER :: rwstatus, unit, naxes(2), nfound, status, ANY_HDU, extver, group, i_ifu, i_stack, j
    INTEGER(kind=4), ALLOCATABLE :: dq(:,:), origin(:,:), IFU(:,:), slice(:,:)
    INTEGER :: top_slice(4), bottom_slice(4)
    INTEGER, PARAMETER :: MUSE_ORIGIN_SHIFT_IFU=6, MUSE_ORIGIN_SHIFT_XSLICE=24, flag=2**20
    LOGICAL :: anyf

    !..define top and bottom slice for each stack
    top_slice(:)=[10,22,34,46]
    bottom_slice(:)=[3,15,27,39]

    !..open file
    print *, "reading pixtable: ", TRIM(pixtable)
    rwstatus=1   ! 1=readwrite
    CALL ftdopn(unit,pixtable,rwstatus,status)
    IF(status/=0) STOP "problem reading file"

    !..get cube size
    CALL ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
    print *, "naxes=",naxes

    !..allocate variables
    ALLOCATE(origin(naxes(1),naxes(2)),&
         dq(naxes(1),naxes(2)),&
         IFU(naxes(1),naxes(2)),&
         slice(naxes(1),naxes(2)))

    !..read extensions
    extver=0
    status=0
    ANY_HDU=-1
    group=1

    CALL ftmnhd(unit, ANY_HDU, "origin", extver, status)
    CALL ftgpvj(unit,group,1,naxes(2),UNDEF,origin,anyf,status)
    IF(status/=0) STOP "problem reading extension 'origin'"

    CALL ftmnhd(unit, ANY_HDU, "dq", extver, status)
    CALL ftgpvj(unit,group,1,naxes(2),UNDEF,dq,anyf,status)
    IF(status/=0) STOP "problem reading extension 'dq'"

    !..convert origin to IFU
    print *, "obtaining IFU information..."
    DO j=1,SIZE(origin)
       IFU(1,j)=IAND(ISHFT(origin(1,j),-MUSE_ORIGIN_SHIFT_IFU),X'1f')
       slice(1,j)=IAND(origin(1,j),X'3f')
    END DO
    DEALLOCATE(origin)

    !..masking slices
    print *, "masking slices..."
    DO i_ifu=1,24
       DO i_stack=1,4
          IF(mask(i_ifu,i_stack)) THEN
             WHERE(IFU(1,:)==i_ifu.and.(slice(1,:)==top_slice(i_stack).or.slice(1,:)==bottom_slice(i_stack))) &
                  dq(1,:)=dq(1,:)+flag
          END IF
       END DO
    END DO

    !..update dq
    print *, "updating file..."
    CALL ftpprj(unit,group,1,naxes(2),dq,status)
    IF(status/=0) STOP "problem writing new dq extension"

    CALL ftclos(unit,status)
    IF(status/=0) STOP "problem closing file"
    
    print *, "done."

  END SUBROUTINE UpdatePixTable


END PROGRAM PixTableMask

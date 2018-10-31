! compile with: ifort pix2wcs.f90 -L/usr/local/scisoft/lib -lcfitsio -o pix2wcs.x
!
! converts pixel position to ra-dec given a WCS fits header 
! 
! author:SC
! based on parts of wcslib.py by Jeremy Brewer (http://wcs2kml.googlecode.com/svn-history/r25/trunk/python/wcslib.py)
!

PROGRAM main

  
  IMPLICIT NONE
  INTEGER(kind=4) :: status, unit, readwrite, ierr, i, narg
  CHARACTER(len=350) :: fname, comment, arg, WCSlabels(11)=["CRPIX1","CRPIX2","CRPIX3","CRVAL1","CRVAL2","CRVAL3","CD1_1 ","CD1_2 ","CD2_1 ","CD2_2 ","CD3_3 "]
  REAL(KIND=8) :: pix_x, pix_y, WCS(11), ra_p, dec_p, rx, ry, phi_p, dx, dy, x, y, r_theta, phi, theta, ra, dec, ra_(3), dec_(3), pix_ra, pix_dec
  REAL(KIND=8), PARAMETER :: pi=3.14159265358979d0
  REAL(KIND=8), PARAMETER :: Rad2Deg=180.0d0/pi, Deg2Rad=pi/180.0
  CHARACTER(len=300):: ra_s, dec_s, inplist, outlist
  LOGICAL :: native_spherical

  CALL ReadArgs

! -- get WCS info from file 

  status=0

  !..get an unused unit
  CALL ftgiou(unit,status)

  IF(status/=0) STOP "problem with ftgiou"

  readwrite=0
  CALL ftdopn(unit,fname,readwrite,status)

  !..get WCS values
  DO i=1,SIZE(WCS)
     CALL ftgkyd(unit,TRIM(WCSlabels(i)),WCS(i),comment,status)
  END DO

  call ftclos(unit,status)

!--- transform coordinates

  IF(pix_x==-1.d0.and.pix_ra==-1.d0) THEN !..reads/write value from/to a list

     OPEN(1,file=inplist,action="read") 
     OPEN(11,file=outlist,action="write")
     ierr=0
     DO 
        READ(1,*,iostat=ierr) pix_x,pix_y
        IF(ierr/=0) EXIT
        CALL transform
        WRITE(11,*) TRIM(ra_s), " ", TRIM(dec_s)
     END DO
     CLOSE(1)
     CLOSE(11)
     
  ELSE !..do a single calculation and produce terminal output

     IF(pix_x/=-1.d0) THEN

        CALL transform
  
        !..terminal output:

        print *, " "
        print *, " ra & dec (deg):      ", ra, dec
        print *, " "
        print *, " ra & dec (sexagesimal): ", TRIM(ra_s), "           ", TRIM(dec_s)
        print *, " "

     ELSE

        CALL RaDec_to_pixel(pix_ra, pix_dec, pix_x, pix_y, native_spherical=.false.)

        print *, " "
        WRITE(*,'(a,f6.2,2x,f6.2)') "pix_x & pix_y = ", pix_x, pix_y

     END IF


  END IF

CONTAINS 

!------------------------------------

  SUBROUTINE ReadArgs

    IMPLICIT NONE

!..get command line arguments
    narg=iargc()
    IF(narg<3) THEN
       print *, "converts pixel position to ra-dec given a WCS fits header and VICEVERSA"
       print *, "usage: pix2wcs.x <fitsfile> pixel_x pixel_y"
       print *, "                   OR                      "
       print *, "       pix2wcs.x <fitsfile> <pixel_input_list> <ra_dec_output_list>"
       print *, " "
       print *, "if pixel_x and pixel_y are in sexagesimal, transforms ra-dec to pixels"
       STOP
    END IF

    pix_x=-1.d0
    pix_y=-1.d0
    pix_ra=-1.d0
    pix_dec=-1.d0

    CALL getarg(1,fname)
    CALL getarg(2,arg)

    IF(INDEX(arg,":")==0) THEN !..read pixels or inplist/outlist

       READ(arg,*,iostat=ierr) pix_x
       IF(ierr/=0) READ(arg,*) inplist
       CALL getarg(3,arg)
       READ(arg,*,iostat=ierr) pix_y
       IF(ierr/=0) READ(arg,*) outlist

    ELSE

       !print *, "read ra & dec"

       !..read ra dec and convert to degrees
       !..splitting strings into pieces

       READ(arg(1:INDEX(arg,":")-1),*) ra_(1)
       READ(arg(INDEX(arg,":")+1:INDEX(arg,":",BACK=.true.)-1),*) ra_(2)
       READ(arg(INDEX(arg,":",BACK=.true.)+1:),*) ra_(3)

       !print *, ra_(:)

       !..convert to degrees
       pix_ra=(ra_(1)+ra_(2)/60.+ra_(3)/3600.)/24.*360.0

       !print *, pix_ra

       CALL getarg(3,arg)

       READ(arg(1:INDEX(arg,":")-1),*) dec_(1)
       READ(arg(INDEX(arg,":")+1:INDEX(arg,":",BACK=.true.)-1),*) dec_(2)
       READ(arg(INDEX(arg,":",BACK=.true.)+1:),*) dec_(3)

       !print *, dec_(:)

       !..convert to degrees
       pix_dec=SIGN((abs(dec_(1))+dec_(2)/60.+dec_(3)/3600.),dec_(1))

       !print *, pix_dec

    END IF


  END SUBROUTINE ReadArgs

!-----------------------------------

  SUBROUTINE transform

    IMPLICIT NONE
 

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

    !..get ra & dec:
    ra=ra_p + datan2(-dcos(theta)*dsin(phi-phi_p), &
         dsin(theta)*dcos(dec_p)-dcos(theta)*dsin(dec_p)*dcos(phi-phi_p))
    dec = dasin(dsin(theta)*dsin(dec_p)+dcos(theta)*dcos(dec_p)*dcos(phi - phi_p))

    ! convert ra, dec to degrees
    ra=ra*Rad2Deg
    IF(ra<0.d0) ra=ra+360.
    dec=dec*Rad2Deg

    !..convert ra dec to sexagesimal
    ra=ra*24.d0/360.d0
    ra_(1)=ra ; ra_(2)=(ra_(1)-INT(ra_(1)))*60. ; ra_(3)=(ra_(2)-INT(ra_(2)))*60.
    dec_(1)=dec ; dec_(2)=abs(dec_(1)-INT(dec_(1)))*60. ; dec_(3)=(dec_(2)-INT(dec_(2)))*60.
    !..readjust critical rounding up cases:
    IF(ra_(3)-INT(ra_(3))>=0.9995) ra_(3)=INT(ra_(3))+1.d0
    IF(dec_(3)-INT(dec_(3))>=0.995) dec_(3)=INT(dec_(3))+1.d0
    !..produce 'well-formatted' strings
    WRITE(ra_s,'(i2.2,a,i2.2,a,i2.2,f0.3)') INT(ra_(1)),":",INT(ra_(2)),":",INT(ra_(3)),ra_(3)-INT(ra_(3))
    WRITE(dec_s,'(i3.2,a,i2.2,a,i2.2,f0.2)') INT(dec_(1)),":",INT(dec_(2)),":",INT(dec_(3)),dec_(3)-INT(dec_(3))

  END SUBROUTINE transform


!---------------------------------------------------------------

  SUBROUTINE RaDec_to_pixel(pix_ra, pix_dec, pix_x, pix_y, native_spherical)

    IMPLICIT NONE
    REAL(kind=8), INTENT(OUT)  :: pix_x, pix_y
    REAL(kind=8), INTENT(IN)   :: pix_ra, pix_dec 
    LOGICAL, INTENT(IN)        :: native_spherical
    REAL(KIND=8), PARAMETER :: pi=3.14159265358979d0
    REAL(KIND=8), PARAMETER :: Rad2Deg=180.0d0/pi, Deg2Rad=pi/180.d0
    REAL(kind=8) :: rx, ry, ra_p, dec_p, phi_p, dx, dy, x, y, r_theta, phi, theta, det, CD_inv(4), ra, dec


!..compute the inverse of the CD matrix
    det=WCS(7)*WCS(10)-WCS(8)*WCS(9)
    IF(det==0) STOP "CD matrix has 0 determinant!"
    CD_inv(1)=WCS(10)/det
    CD_inv(2)=WCS(8)/det
    CD_inv(3)=WCS(9)/det
    CD_inv(4)=WCS(7)/det

    IF(.not.native_spherical) THEN !..convert to native spherical coordinates
       
       !..convert angles to radians
       ra=pix_ra*Deg2Rad
       dec=pix_dec*Deg2Rad

       !..convert WCS to radians
       ra_p=WCS(4)*Deg2Rad
       dec_p=WCS(5)*Deg2Rad

       ! native longitude of the celestial pole in radians
       phi_p = pi
       IF(dec_p >= pi/2.d0) phi_p = 0.0

       !..convert to native spherical
       phi = phi_p + datan2(-dcos(dec) * dsin(ra - ra_p), &
            dsin(dec) * dcos(dec_p) - &
            dcos(dec) * dsin(dec_p) * &
            dcos(ra - ra_p))
       theta = dasin(dsin(dec) * dsin(dec_p) + &
                          dcos(dec) * dcos(dec_p) * &
                          dcos(ra - ra_p))

    ELSE

       phi=pix_ra
       theta=pix_dec

    END IF

!..compute the intermediate world coordinates
    rx=(180.d0/pi)*dsin(phi)/dtan(theta)
    ry=-(180.d0/pi)*dcos(phi)/dtan(theta)

!..compute pixel position
    pix_x=CD_inv(1)*rx-CD_inv(2)*ry+WCS(1)
    pix_y=-CD_inv(3)*rx+CD_inv(4)*ry+WCS(2)

  END SUBROUTINE RaDec_to_pixel


END PROGRAM main

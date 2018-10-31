SUBROUTINE WriteCatalogue(ThisCatalogue)

  USE Globalmodule
  IMPLICIT NONE
  CHARACTER(len=250), INTENT(IN) :: ThisCatalogue
  INTEGER(kind=4) :: i, ii, col_n(23), col_n_WCS(29)
  CHARACTER(len=300) :: xcen_ra, xcen_dec, lcen_ra, lcen_dec, xcen_lambda, lcen_lambda

  IF(COUNT(Obj(:)%Id>0)==0) RETURN !..no objects for the catalogue

  IF(ALL(WCS==0.d0)) THEN !..produce standard output without ra,dec and lambda

     DO i=1,23; col_n(i)=i; END DO
    
        OPEN(1,file=ThisCatalogue,action="write")

        IF(printheader) THEN
           WRITE(1,'(a)') "###################################################"
           WRITE(1,'(a)') "# col 1:     Id "
           WRITE(1,'(a)') "# col 2:     number of connected spaxels above detection threshold "
           WRITE(1,'(a)') "# col 3-5:   geometrical centroid position (units=spaxels) "
           WRITE(1,'(a)') "# col 6-8:   flux-weighted-centroid position (units=spaxels) "
           WRITE(1,'(a)') "# col 9-11:  minimum coordinates of the box containing the object (units=spaxels) "
           WRITE(1,'(a)') "# col 12-14: maximum coordinates of the box containing the object (units=spaxels) "
           WRITE(1,'(a)') "# col 15:    projected iso-area (pixels)"
           WRITE(1,'(a)') "# col 16:    projected size in wavelenght dimension (pixels)"
           WRITE(1,'(a)') "# col 17-18: IsoFlux and IsoErr (units of the original datacube)"
           WRITE(1,'(a)') "# col 19-20: Flux and Err within the box defined by the coordinates in col 9-14"
           WRITE(1,'(a)') "# col 21-22: Flux and Err in the Cylindrical apertures defined in the parameter file"
           WRITE(1,'(a)') "# col 23:    Associated object Id (=0 if no association found/requested)"
           WRITE(1,'(a)') "# "
           WRITE(1,'(a,i7,1x,i8,1x,3(i7,1x),1x,3(i7,1x),6(i7,x),1x,2(i7,1x),6(i12,1x),i7)') "#",(col_n(i),i=1,22)
           WRITE(1,'(a)') "# "
        END IF
        ii=0
        DO i=1,SIZE(Obj)
           IF(Obj(i)%Id>0) THEN
              ii=ii+1
              WRITE(1,'(2(i8,x),1x,3(f7.2,1x),1x,3(f7.2,1x),6(i7,x),1x,2(i7,1x),6(es12.4,1x),i7)') Obj(i)%Id, Obj(i)%NSpax, Obj(i)%xcen(1:3), Obj(i)%lcen(1:3), Obj(i)%boxmin, Obj(i)%boxmax, &
                   Obj(i)%area, Obj(i)%dz, Obj(i)%IsoFlux, Obj(i)%IsoErr, Obj(i)%BoxFlux, Obj(i)%BoxErr, Obj(i)%AperFlux, Obj(i)%AperErr, Obj(i)%Assoc
           END IF
        END DO
        CLOSE(1)

     ELSE  !..compute and add ra, dec, and lambda in the catalogue

      DO i=1,29; col_n_WCS(i)=i; END DO
    
        OPEN(1,file=ThisCatalogue,action="write")

        IF(printheader) THEN
           WRITE(1,'(a)') "###################################################"
           WRITE(1,'(a)') "# col 1:     Id "
           WRITE(1,'(a)') "# col 2:     number of connected spaxels above detection threshold "
           WRITE(1,'(a)') "# col 3-5:   geometrical centroid position (units=spaxels) "
           WRITE(1,'(a)') "# col 6-8:   flux-weighted-centroid position (units=spaxels) "
           WRITE(1,'(a)') "# col 9-11:  minimum coordinates of the box containing the object (units=spaxels) "
           WRITE(1,'(a)') "# col 12-14: maximum coordinates of the box containing the object (units=spaxels) "
           WRITE(1,'(a)') "# col 15:    projected iso-area (pixels)"
           WRITE(1,'(a)') "# col 16:    projected size in wavelenght dimension (pixels)"
           WRITE(1,'(a)') "# col 17-18: IsoFlux and IsoErr (units of the original datacube)"
           WRITE(1,'(a)') "# col 19-20: Flux and Err within the box defined by the coordinates in col 9-14"
           WRITE(1,'(a)') "# col 21-22: Flux and Err in the Cylindrical apertures defined in the parameter file"
           WRITE(1,'(a)') "# col 23-25: geometrical centroid in units of ra, dec and lambda"
           WRITE(1,'(a)') "# col 26-28: flux-weighted-centroid in units of ra, dec and lambda"
           WRITE(1,'(a)') "# col 29:    Associated object Id (=0 if no association found/requested)"
           WRITE(1,'(a)') "# "
           WRITE(1,'(a,i7,1x,i8,1x,3(i7,1x),1x,3(i7,1x),6(i7,x),1x,2(i7,1x),6(i12,1x),2(2(i13,1x),i8,1x),i7)') "#",(col_n_WCS(i),i=1,29)
           WRITE(1,'(a)') "# "
        END IF
        ii=0
        DO i=1,SIZE(Obj)

           IF(Obj(i)%Id>0) THEN

              ii=ii+1
              
              ! --- get ra & dec
              CALL pix2wcs(REAL(Obj(i)%xcen(1),KIND=8), REAL(Obj(i)%xcen(2),KIND=8), xcen_ra, xcen_dec)
              CALL pix2wcs(REAL(Obj(i)%lcen(1),KIND=8), REAL(Obj(i)%lcen(2),KIND=8), lcen_ra, lcen_dec)
              !..get lambda
              WRITE(xcen_lambda,'(f8.2)') (Obj(i)%xcen(3)-WCS(3))*WCS(11)+WCS(6)
              WRITE(lcen_lambda,'(f8.2)') (Obj(i)%lcen(3)-WCS(3))*WCS(11)+WCS(6)


              WRITE(1,'(2(i8,x),1x,3(f7.2,1x),1x,3(f7.2,1x),6(i7,x),1x,2(i7,1x),6(es12.4,1x),2(2(a13,a1),a8,a1),i7)') Obj(i)%Id, Obj(i)%NSpax, Obj(i)%xcen(1:3), Obj(i)%lcen(1:3), Obj(i)%boxmin, Obj(i)%boxmax, &
                   Obj(i)%area, Obj(i)%dz, Obj(i)%IsoFlux, Obj(i)%IsoErr, Obj(i)%BoxFlux, Obj(i)%BoxErr, Obj(i)%AperFlux, Obj(i)%AperErr, &
                   TRIM(xcen_ra), " ", TRIM(xcen_dec), " ", TRIM(xcen_lambda), " ",TRIM(lcen_ra)," ", TRIM(lcen_dec), " ",TRIM(lcen_lambda), &
                   " ", INT(Obj(i)%Assoc)


           END IF
        END DO
        CLOSE(1)       

     ENDIF


CONTAINS

!----------------------------------------------------------
! converts pixel position to ra-dec given a WCS fits header 
! based on parts of wcslib.py by Jeremy Brewer (http://wcs2kml.googlecode.com/svn-history/r25/trunk/python/wcslib.py)
!
!
  SUBROUTINE pix2wcs(pix_x, pix_y, ra_s, dec_s)

    IMPLICIT NONE
    REAL(kind=8), INTENT(IN) :: pix_x, pix_y
    CHARACTER(len=300), INTENT(OUT) :: ra_s, dec_s
    REAL(KIND=8) :: ra_p, dec_p, rx, ry, phi_p, dx, dy, x, y, r_theta, phi, theta, ra, dec, ra_(3), dec_(3)
    REAL(KIND=8), PARAMETER :: pi=3.14159265358979d0
    REAL(KIND=8), PARAMETER :: Rad2Deg=180.0d0/pi, Deg2Rad=pi/180.0


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

  
  END SUBROUTINE pix2wcs

END SUBROUTINE WriteCatalogue



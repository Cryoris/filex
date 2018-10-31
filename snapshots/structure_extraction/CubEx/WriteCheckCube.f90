SUBROUTINE WriteCheckCube

  USE Globalmodule
  IMPLICIT NONE
  CHARACTER(len=250) :: bovfile, cmd
  CHARACTER(len=8) :: varname
  REAL, ALLOCATABLE :: CheckCubeData(:,:,:)
  INTEGER :: i,j,k, status, naxis, naxes(3), group, blocksize, bitpix, unit, nelements, fpixel, ii, thisId
  LOGICAL :: ex, extend, simple

  ALLOCATE(CheckCubeData(1:DimX,1:DimY,1:DimZ))

  DO ii=1,NCheckCubes

     print *, "Writing ChekCube file:", TRIM(CheckCube(ii))
     print *, "format=", TRIM(CheckCubeFMT)

!..produce checkimage from cube excluding ghost zones
     SELECT CASE(TRIM(CheckCubeType(ii)))
     CASE("Objects")    
        WHERE(Mask(1:DimX,1:DimY,1:DimZ)/=0) 
           CheckCubeData=Cube(1:DimX,1:DimY,1:DimZ)
        ELSEWHERE
           CheckCubeData=UNDEF
        END WHERE
        varname="Flux"
     CASE("Residuals")    
        WHERE(Mask(1:DimX,1:DimY,1:DimZ)==0) 
           CheckCubeData=Cube(1:DimX,1:DimY,1:DimZ)
        ELSEWHERE
           CheckCubeData=UNDEF
        END WHERE
        varname="Flux"
     CASE("Objects_F")    
        WHERE(Mask(1:DimX,1:DimY,1:DimZ)/=0) 
           CheckCubeData=Cube(1:DimX,1:DimY,1:DimZ)
        ELSEWHERE
           CheckCubeData=UNDEF
        END WHERE
        varname="Flux"
     CASE("SNR")
        WHERE(Var(1:DimX,1:DimY,1:DimZ)/=UNDEF.and.Var(1:DimX,1:DimY,1:DimZ)/=0.)
           CheckCubeData=Cube(1:DimX,1:DimY,1:DimZ)/sqrt(Var(1:DimX,1:DimY,1:DimZ))
        ELSEWHERE
           CheckCubeData=0.
        END WHERE
        varname="SNR"
     CASE("SNR_F")
        WHERE(VarF(1:DimX,1:DimY,1:DimZ)/=UNDEF.and.VarF(1:DimX,1:DimY,1:DimZ)/=0.)
           CheckCubeData=CubeF(1:DimX,1:DimY,1:DimZ)/sqrt(VarF(1:DimX,1:DimY,1:DimZ))
        ELSEWHERE
           CheckCubeData=0.
        END WHERE
        varname="SNR"
     CASE("Objects_SNR")
        WHERE(Mask(1:DimX,1:DimY,1:DimZ)/=0)
           CheckCubeData=Cube(1:DimX,1:DimY,1:DimZ)/sqrt(Var(1:DimX,1:DimY,1:DimZ))
        ELSEWHERE
           CheckCubeData=Objects_SNR_Floor
        END WHERE
        varname="SNR"
     CASE("Objects_SNR_F")
        WHERE(Mask(1:DimX,1:DimY,1:DimZ)/=0)
           CheckCubeData=CubeF(1:DimX,1:DimY,1:DimZ)/sqrt(VarF(1:DimX,1:DimY,1:DimZ))
        ELSEWHERE
           CheckCubeData=Objects_SNR_Floor
        END WHERE
        varname="SNR"
     CASE("Objects_Id")
        print *, "HaloLabel", HaloLabel
        DO k=1,DimZ
           DO j=1,DimY
              DO i=1,DimX
                 IF(mask(i,j,k)>0) THEN
                    CheckCubeData(i,j,k)=LabelToId(mask(i,j,k))
                 ELSE
                    CheckCubeData(i,j,k)=0.
                 END IF
              END DO
           END DO
        END DO
        varname="Id"
     CASE("Objects_Id_Assoc")
        DO k=1,DimZ
           DO j=1,DimY
              DO i=1,DimX
                 IF(mask(i,j,k)>0) THEN
                    thisId=LabelToId(mask(i,j,k))
                    IF(Obj(thisId)%Assoc/=0) THEN
                       CheckCubeData(i,j,k)=Obj(thisId)%Assoc
                    ELSE
                       CheckCubeData(i,j,k)=thisId
                    END IF
                 ELSE
                    CheckCubeData(i,j,k)=0.
                 END IF
              END DO
           END DO
        END DO
        varname="Id"
     CASE("CubeF")
        CheckCubeData=CubeF(1:DimX,1:DimY,1:DimZ)
     CASE("Cube_F")
        CheckCubeData=CubeF(1:DimX,1:DimY,1:DimZ)
     CASE("VarF")
        CheckCubeData=VarF(1:DimX,1:DimY,1:DimZ)
     CASE("TEST")   !..reserved for testing purposes
        CheckCubeData=VarF(1:DimX,1:DimY,1:DimZ)/sqrt(Var(1:DimX,1:DimY,1:DimZ))
     CASE default
        print *, "CheckCubeType not recognized! CheckCube will not be produced..."
        RETURN
     END SELECT

     IF(TRIM(CheckCubeFMT)=="fits") THEN

        status=0

        !..get an unused unit
        CALL ftgiou(unit,status)
     
        IF(status/=0) STOP "problem with ftgiou"
     
        INQUIRE(file=CheckCube(ii),EXIST=ex)
        IF(ex) THEN !..open, delete previous data and re-initialize the file 
           CALL ftdopn(unit,CheckCube(ii),1,status)  
           IF(status/=0) STOP "problem with ftdopn"
           CALL ftdelt(unit, status)
           IF(status/=0) STOP "problem with ftdelt"
        ENDIF
        blocksize=1
        CALL ftinit(unit,CheckCube(ii),blocksize,status)
        IF(status/=0) STOP "problem with ftinit"

        !  Initialize parameters about the FITS image.
        simple=.true.
        extend=.true.
        bitpix=-32
        naxis=3
        naxes=[DimX,DimY,DimZ] ! excluding ghost zones

        !  Write the required header keywords to the file
        call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
        IF(status/=0) STOP "problem with ftphpr"
        
        !..Write the name of the variable
        CALL ftpkys(unit,"OBJECT",TRIM(varname),' ',status)

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
        group=1
        fpixel=1
        nelements=PRODUCT(naxes)
        call ftppre(unit,group,fpixel,nelements,CheckCubeData,status)
        IF(status/=0) STOP "problem with ftppre"

        call ftclos(unit, status)
        IF(status/=0) STOP "problem with ftclos!"

     ELSEIF(TRIM(CheckCubeFMT)=="bov") THEN
           
        !..write output
        OPEN(1,file=CheckCube(ii),action='write',form='unformatted')
        WRITE(1) CheckCubeData(1:DimX,1:DimY,1:DimZ)
        CLOSE(1)

        !..write bov header
        bovfile=TRIM(CheckCube(ii))//'.bov'
        OPEN(1,file=bovfile,action='write',form='formatted')
        WRITE(1,*) 'TIME: 0'
        WRITE(1,*) 'DATA_FILE: ',TRIM(CheckCube(ii))
        WRITE(1,*) 'DATA_SIZE: ',DimX, DimY, DimZ
        WRITE(1,*) 'DATA_FORMAT: FLOAT'
        WRITE(1,*) 'VARIABLE: CheckCube'
        WRITE(1,*) 'DATA_ENDIAN: LITTLE'
        WRITE(1,*) 'CENTERING: zonal'
        WRITE(1,*) 'BRICK_ORIGIN: 0 0 0'
        WRITE(1,*) 'BRICK_SIZE: ',DimX, DimY, DimZ
        WRITE(1,*) 'BYTE_OFFSET: 4'
        CLOSE(1)

     END IF

  END DO

  DEALLOCATE(CheckCubeData)

END SUBROUTINE WriteCheckCube

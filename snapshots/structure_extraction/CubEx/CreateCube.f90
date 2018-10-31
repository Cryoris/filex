SUBROUTINE CreateCube

  USE Globalmodule
  IMPLICIT NONE
  INTEGER(kind=4) :: nx=300, nobj=5, i, pos(3), ii
  REAL(kind=4)    :: R(3), size_(3)
  CHARACTER(len=300) :: outfile, bovfile
  
!..dimension without ghost zones
  DimX=nx
  DimY=nx
  DimZ=nx


  ALLOCATE(Cube(0:nx+1,0:nx+1,0:nx+1),Var(0:nx+1,0:nx+1,0:nx+1),Mask(0:nx+1,0:nx+1,0:nx+1))


  !..create a random variance cube with a median of 1
  CALL RANDOM_NUMBER(Var)
  !Var=Var+0.5

  Var=1

  !..create a random background cube with a median of 1
  CALL RANDOM_NUMBER(Cube)
  !Cube=Cube+0.5 
  
  Cube=0

  !..create a set of objects at random locations with random sizes
  DO i=1,nobj

     CALL RANDOM_NUMBER(R)
     pos=INT(R*nx)

     WHERE(pos<10) pos=10
     WHERE(pos>nx-10) pos=nx-10

     CALL RANDOM_NUMBER(size_)
     DO ii=1,3
        size_(ii)=MAX(3.,size_(ii)*10.)
     END DO

     print *, "pos=", pos
     print *, "boxmin=", INT(pos(:)-size_(:)/2)
     print *, "boxmax=", INT(pos(:)+size_(:)/2)

     Cube(pos(1)-INT(size_(1))/2:pos(1)+INT(size_(1))/2,pos(2)-INT(size_(2))/2:pos(2)+INT(size_(2))/2,pos(3)-INT(size_(3))/2:pos(3)+INT(size_(3))/2)=&
          5
         !Cube(pos(1)-size_(1)/2:pos(1)+size_(1)/2,pos(2)-size_(2)/2:pos(2)+size_(2)/2,pos(3)-size_(3)/2:pos(3)+size_(3)/2)*3.

  END DO


!... write cube

  outfile="CubeTest"
           
  !..write output
  OPEN(1,file=outfile,action='write',form='unformatted')
  WRITE(1) Cube(1:nx,1:nx,1:nx)/Var(1:nx,1:nx,1:nx)
  CLOSE(1)

  !..write bov header
  bovfile=TRIM(outfile)//'.bov'
  OPEN(1,file=bovfile,action='write',form='formatted')
  WRITE(1,*) 'TIME: 0'
  WRITE(1,*) 'DATA_FILE: ',TRIM(outfile)
  WRITE(1,*) 'DATA_SIZE: ',nx,nx,nx
  WRITE(1,*) 'DATA_FORMAT: FLOAT'
  WRITE(1,*) 'VARIABLE: FLUX'
  WRITE(1,*) 'DATA_ENDIAN: LITTLE'
  WRITE(1,*) 'CENTERING: zonal'
  WRITE(1,*) 'BRICK_ORIGIN: 0 0 0'
  WRITE(1,*) 'BRICK_SIZE: ',nx,nx,nx
  WRITE(1,*) 'BYTE_OFFSET: 4'
  CLOSE(1)

 


END SUBROUTINE CreateCube

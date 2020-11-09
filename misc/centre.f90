PROGRAM centre
  IMPLICIT NONE
  INTEGER :: i, j, natom, nstart, io, ncon, k
  CHARACTER :: pdbfile*120, lign80*80, dummy*120
  REAL(8) :: xm,ym,zm,xmax,ymax,zmax,xmin,ymin,zmin
  INTEGER, ALLOCATABLE,DIMENSION(:) :: atnum,resnum,resvals
  REAL(8), ALLOCATABLE,DIMENSION(:) ::x,y,z,occ,bfac,vals,mass
  CHARACTER, ALLOCATABLE, DIMENSION(:) :: atom*6,name*4,res*3,chain*1,elem*2,chag*2,stlin*80,con*80
  LOGICAL :: getoption,alig


  IF(.not.getoption('-pdb',.true.,pdbfile)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF

  IF(getoption('-fit',.true.,dummy)) THEN
     alig=.true.
  ELSE
     alig=.false.
  END IF


  OPEN(file=pdbfile,form='FORMATTED',status='OLD',unit=2345)
  
  i=0
  j=0
  k=0
  DO
     READ(2345,'(A)',IOSTAT=io) lign80
     IF (io > 0) THEN
        WRITE(*,*) 'Check input.  Something was wrong'
        STOP
        ! End of file
     ELSE IF (io < 0) THEN
        EXIT
        ! Count number of atoms
     ELSE
        IF ((lign80(1:4).eq.'ATOM').or. &
             (lign80(1:6).eq.'HETATM')) THEN
           i=i+1
        ELSE IF ((lign80(1:3).ne.'END').and.(lign80(1:3).ne.'TER').and. &
             (lign80(1:6).ne.'CONECT').and.(lign80(1:6).ne.'MASTER')) THEN
           j=j+1
        ELSE IF (lign80(1:6).eq.'CONECT') THEN
           k=k+1
        END IF
     END IF
  END DO

  REWIND(2345) 

  natom=i
  nstart=j
  ncon=k

  ALLOCATE(x(natom))
  ALLOCATE(y(natom))
  ALLOCATE(z(natom))
  ALLOCATE(atom(natom))
  ALLOCATE(atnum(natom))
  ALLOCATE(name(natom))
  ALLOCATE(res(natom))
  ALLOCATE(chain(natom))
  ALLOCATE(resnum(natom))
  ALLOCATE(occ(natom))
  ALLOCATE(bfac(natom))
  ALLOCATE(elem(natom))
  ALLOCATE(chag(natom))
  ALLOCATE(mass(natom))
  ALLOCATE(stlin(nstart))
  ALLOCATE(con(ncon))

  i=0
  j=0
  k=0
  DO
     READ(2345,'(A)',IOSTAT=io) lign80
     IF (io > 0) THEN
        WRITE(*,*) 'Check input.  Something was wrong'
        STOP
        ! End of file
     ELSE IF (io < 0) THEN
        EXIT
        ! Count number of atoms
     ELSE
        IF ((lign80(1:4).eq.'ATOM').or. &
             (lign80(1:6).eq.'HETATM')) THEN
           i=i+1
           READ(lign80,'(A6,I5,1X,A4,1X,A3,1X,A1,I4,1X,3X,3F8.3,2F6.2,10X,2A2)') &
                atom(i),atnum(i),name(i),res(i),chain(i),resnum(i),x(i),y(i),z(i),&
                occ(i),bfac(i),elem(i),chag(i)
        ELSE IF ((lign80(1:3).ne.'END').and.(lign80(1:3).ne.'TER').and. &
             (lign80(1:6).ne.'CONECT').and.(lign80(1:6).ne.'MASTER')) THEN
           j=j+1
           READ(lign80,'(A)') stlin(j)
        ELSE IF (lign80(1:6).eq.'CONECT') THEN
           k=k+1
           READ(lign80,'(A)') con(k)
        END IF
     END IF
  END DO

  CLOSE(2345)

  DO i=1,natom
     IF (elem(i).eq." H") THEN
        mass(i)=1.008
     ELSE IF (elem(i).eq." C") THEN
        mass(i)=12.011
     ELSE IF (elem(i).eq." N") THEN
        mass(i)=14.007
     ELSE IF (elem(i).eq." O") THEN
        mass(i)=15.999
     ELSE IF (elem(i).eq." P") THEN
        mass(i)=30.974
     ELSE IF (elem(i).eq." S") THEN
        mass(i)=32.066
     ELSE
        WRITE(6,*) 'Unknown atom type please add to list.  ', elem(i)
        mass(i)=1
     END IF
  END DO

  xm=0
  ym=0
  zm=0
  DO i=1,natom
     xm=xm+x(i)
     ym=ym+y(i)
     zm=zm+z(i)
  END DO
  xm=xm/natom
  ym=ym/natom
  zm=zm/natom

  DO i=1,natom
     x(i)=x(i)-xm
     y(i)=y(i)-ym
     z(i)=z(i)-zm
  END DO  

  xmax=0
  ymax=0
  zmax=0
  xmin=0
  ymin=0
  zmin=0

  DO i=1,natom
     xmax=max(x(i),xmax)
     ymax=max(y(i),ymax)
     zmax=max(z(i),zmax)
     xmin=min(x(i),xmin)
     ymin=min(y(i),ymin)
     zmin=min(z(i),zmin)
  END DO  

  IF (alig) THEN
     xm=(xmin+xmax)/2
     ym=(ymin+ymax)/2
     zm=(zmin+zmax)/2
     
     DO i=1,natom
        x(i)=x(i)-xm
        y(i)=y(i)-ym
        z(i)=z(i)-zm
     END DO

  END IF
  xmax=0
  ymax=0
  zmax=0
  xmin=0
  ymin=0
  zmin=0
  DO i=1,natom
     xmax=max(x(i),xmax)
     ymax=max(y(i),ymax)
     zmax=max(z(i),zmax)
     xmin=min(x(i),xmin)
     ymin=min(y(i),ymin)
     zmin=min(z(i),zmin)
  END DO

  WRITE(*,'(A,3(1X,F7.3))') 'Max coordinates', xmax, ymax, zmax
  WRITE(*,'(A,3(1X,F7.3))') 'Min coordinates', xmin, ymin, zmin















  OPEN(file='centred.pdb',form='FORMATTED',unit=2346)

  DO i=1,nstart
     WRITE(2346,'(A)') stlin(i)
  END DO
  
  DO i=1,natom
     WRITE(2346,'(A6,I5,1X,A4,1X,A3,1X,A1,I4,1X,3X,3F8.3,2F6.2,10X,2A2)') &
                atom(i),atnum(i),name(i),res(i),chain(i),resnum(i),x(i),y(i),z(i),&
                occ(i),bfac(i),elem(i),chag(i)
  END DO

  DO i=1,ncon
     WRITE(2346,'(A)') con(i)
  END DO

  CLOSE(2346)

END PROGRAM centre

SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
 
  WRITE(iunit,'(A)')"                                                                         "

END SUBROUTINE helptext


FUNCTION getoption(flag,getval,cvalue)
  IMPLICIT NONE
  CHARACTER(*) :: flag,cvalue
  CHARACTER(160) :: arg
  LOGICAL :: getoption,getval
  INTEGER :: l,i,j
  
  getoption=.false.
  i=0
  DO j=1,iargc()
     CALL getarg(j,arg)
     IF (arg.eq.flag) i=j
  END DO
  IF (i.gt.0) THEN
     getoption=.true.
  END IF
  IF (getval) call getarg(i+1,cvalue)
END FUNCTION getoption

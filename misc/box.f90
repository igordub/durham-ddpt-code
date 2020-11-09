PROGRAM box
  IMPLICIT NONE
  INTEGER :: i, j, npro, nwat, io,k, numcon
  CHARACTER :: pdbfile*120, lign80*80, dummy*120,st1*30,st2*30,st3*30,en1*15,en2*15,en3*15,group*7
  REAL(8) :: xm,ym,zm,xmax,ymax,zmax,xmin,ymin,zmin,x1,y1,z1,x2,y2,z2,x3,y3,z3,dim
  INTEGER, ALLOCATABLE,DIMENSION(:) :: atnum,resnum,resvals
  REAL(8), ALLOCATABLE,DIMENSION(:) ::x,y,z,occ,bfac,vals,mass
  CHARACTER, ALLOCATABLE, DIMENSION(:) :: atom*6,name*4,res*3,chain*1,elem*2,chag*2,stlin*80,con*80
  LOGICAL :: getoption,alig

  IF(.not.getoption('-pdb',.true.,pdbfile)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF

  IF(getoption('-dim',.true.,dummy)) THEN
     READ(dummy,*)dim
  ELSE
     dim=100
     WRITE(*,*) "dimesion set to 100"
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
           IF (lign80(18:20).ne."WAT") THEN
              i=i+1
           ELSE
              j=j+1
           END IF
        ELSE IF (lign80(1:6).eq.'CONECT') THEN
           k=k+1
        END IF
     END IF
  END DO

  REWIND(2345) 

  npro=i
  nwat=j
  numcon=k

  ALLOCATE(x(npro))
  ALLOCATE(y(npro))
  ALLOCATE(z(npro))
  ALLOCATE(atom(npro))
  ALLOCATE(atnum(npro))
  ALLOCATE(name(npro))
  ALLOCATE(res(npro))
  ALLOCATE(chain(npro))
  ALLOCATE(resnum(npro))
  ALLOCATE(occ(npro))
  ALLOCATE(bfac(npro))
  ALLOCATE(elem(npro))
  ALLOCATE(chag(npro))
  ALLOCATE(mass(npro))
  ALLOCATE(con(numcon))

  i=0
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
           IF (lign80(18:20).ne."WAT") THEN
              i=i+1
              READ(lign80,'(A6,I5,1X,A4,1X,A3,1X,A1,I5,3X,3F8.3,2F6.2,10X,2A2)') &
                   atom(i),atnum(i),name(i),res(i),chain(i),resnum(i),x(i),y(i),z(i),&
                   occ(i),bfac(i),elem(i),chag(i)
           END IF
        ELSE IF (lign80(1:6).eq.'CONECT') THEN
           k=k+1
           READ(lign80,'(A)') con(k)
        END IF
     END IF
  END DO

  DO i=1,npro
     IF (index(name(i),'C').gt.0) THEN
        elem(i)=" C"
     ELSE IF ((index(name(i),'NA').gt.0).or.(index(name(i),'NA').gt.0)) THEN
        elem(i)="Na"   
     ELSE IF (index(name(i),'N').gt.0) THEN
        elem(i)=" N"
     ELSE IF (index(name(i),'P').gt.0) THEN
        elem(i)=" P"
     ELSE IF (index(name(i),'S').gt.0) THEN
        elem(i)=" S"
     ELSE IF (index(name(i),'O').gt.0) THEN
        elem(i)=" O"   
     ELSE IF (index(name(i),'H').gt.0) THEN
        elem(i)=" H"
     ELSE 
        WRITE(*,*) 'name', name(i), 'not found. Please add.'
     END IF
  END DO

  REWIND(2345)

  xmax=0
  ymax=0
  zmax=0
  xmin=0
  ymin=0
  zmin=0
  DO i=1,npro
     xmax=max(x(i),xmax)
     ymax=max(y(i),ymax)
     zmax=max(z(i),zmax)
     xmin=min(x(i),xmin)
     ymin=min(y(i),ymin)
     zmin=min(z(i),zmin)
  END DO


  OPEN(file='out.pdb',form='FORMATTED',unit=2346)  
  DO i=1,npro
     IF (chain(i).eq."A") THEN
        group="chainA "
     ELSE IF (chain(i).eq."B") THEN
        group="chainB "
     ELSE
        group="chainC "
     END IF
     
     WRITE(2346,'(A6,I5,1X,A4,1X,A3,1X,A1,I4,1X,3X,3F8.3,2F6.2,10X,2A2)') &
                atom(i),atnum(i),name(i),res(i),chain(i),resnum(i),x(i),y(i),z(i),&
                occ(i),bfac(i),elem(i),chag(i)
  END DO

  DEALLOCATE(x)
  DEALLOCATE(y)
  DEALLOCATE(z)
  DEALLOCATE(atom)
  DEALLOCATE(atnum)
  DEALLOCATE(name)
  DEALLOCATE(res)
  DEALLOCATE(chain)
  DEALLOCATE(resnum)
  DEALLOCATE(occ)
  DEALLOCATE(bfac)
  DEALLOCATE(elem)
  DEALLOCATE(chag)
  DEALLOCATE(mass)

  i=0
  dim=dim/2
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
           IF (lign80(18:20).eq."WAT") THEN
              READ(lign80,'(A30,3F8.3,A15)') st1,x1,y1,z1,en1
              READ(2345,'(A30,3F8.3,A15)') st2,x2,y2,z2,en2
              READ(2345,'(A30,3F8.3,A15)') st3,x3,y3,z3,en3

              IF ((dabs(x1).gt.dim).or.(dabs(y1).gt.dim).or.(dabs(z1).gt.dim).or. &
                   (dabs(x2).gt.dim).or.(dabs(y2).gt.dim).or.(dabs(z2).gt.dim).or. &
                   (dabs(x3).gt.dim).or.(dabs(y3).gt.dim).or.(dabs(z3).gt.dim)) THEN
                 WRITE(*,*) "Removed ", st1
              ELSE
                 i=i+1
 !                WRITE(group,'(I7)') i
                 group='sol    '
                 WRITE(2346,'(A30,3F8.3,A15,7X,2A2)') st1,x1,y1,z1,en1," O","  " 
                 WRITE(2346,'(A30,3F8.3,A15,7X,2A2)') st2,x2,y2,z2,en2," H","  "
                 WRITE(2346,'(A30,3F8.3,A15,7X,2A2)') st3,x3,y3,z3,en3," H","  "
              END IF
           END IF
        END IF
     END IF
  END DO

  DO k=1,numcon
      WRITE(2346,'(A)') con(k)
  END DO

  WRITE(2346,'(A)') 'END'


  WRITE(*,'(A,3(1X,F7.3))') 'Max coordinates', xmax, ymax, zmax
  WRITE(*,'(A,3(1X,F7.3))') 'Min coordinates', xmin, ymin, zmin

  CLOSE(2345)
  CLOSE(2346)

END PROGRAM box

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

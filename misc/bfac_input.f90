PROGRAM bfac_input
  IMPLICIT NONE
  LOGICAL :: getoption, hetatm, caonly, resmass, atmass, cutvect,cutassign,forceres,ukn
  CHARACTER :: dummy*120, pdbfile*120, lign80*80, nomatm*120, nomres*120,pdbfile2*120,fill(100)*4,fill2(100)*4
  INTEGER :: natom, io, i, j, natomold, nident, ncusres,ii,ijjjj,iseed,jat,jj,ll,nntr,nnzero, pos(100),pos2(100), &
       length1,length2,natom2, k
  REAL(8) :: cutoff, rave, rdev, rmin, rmax,anmp,ddf,rkh,dist,dist2,dmax,dmin,distave,drms,kij, &
       shift,rx,ry,rz,trace,random, sum
  LOGICAL,ALLOCATABLE,DIMENSION(:) :: otherres
  INTEGER,ALLOCATABLE,DIMENSION(:) :: atnum,resnum,resnumold, resone, restwo
  REAL(8),ALLOCATABLE,DIMENSION(:) :: x,y,z,occ,bfac,vals,mass,massold,cutvalue,acutoff,kijcust
  CHARACTER,ALLOCATABLE,DIMENSION(:) :: atom*6,name*4,res*3,chain*1,elem*2,chag*2,chainold*1,ident*4
  INTEGER,ALLOCATABLE,DIMENSION(:) :: atnum2,resnum2
  REAL(8),ALLOCATABLE,DIMENSION(:) :: x2,y2,z2,occ2,bfac2,vals2,mass2
  CHARACTER,ALLOCATABLE,DIMENSION(:) :: atom2*6,name2*4,res2*3,chain2*1,elem2*2,chag2*2

  IF (getoption('-help',.false.,dummy)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF

  iseed=27041961
  shift=1e-8
!----------------------------------------------------------------------------------------
! Read in options

  IF (.not.getoption('-pdb',.true.,pdbfile)) THEN
     WRITE(0,'(A)') "Need to input a pdb file with -pdb"
     CALL helptext(0)
     CALL exit(0)
  END IF

  IF (.not.getoption('-pdb2',.true.,pdbfile2)) THEN
     WRITE(0,'(A)') "Need to input a pdb file with -pdb"
     CALL helptext(0)
     CALL exit(0)
  END IF

!------------------------------------------------------------------------------------
! Read in pdb file
  OPEN(file=pdbfile,form="FORMATTED",status="OLD",unit=2356)

  i=0
  DO
     READ(2356,'(A)',IOSTAT=io) lign80
     IF (io > 0) THEN
        WRITE(*,*) 'Check input.  Something was wrong'
        STOP
        ! End of file
     ELSE IF (io < 0) THEN
        EXIT
        ! Count number of atoms
     ELSE
        IF ((lign80(1:4).eq.'ATOM').or. &
             ((lign80(1:6).eq.'HETATM').and.hetatm)) THEN
           i=i+1
        END IF
     END IF
  END DO
  
  natom=i

  REWIND(2356)

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

  REWIND(2356)
 
  i=0
  DO
     READ(2356,'(A)',IOSTAT=io) lign80
     IF (io > 0) THEN
        WRITE(*,*) 'Check input.  Something was wrong'
        STOP
        ! End of file
     ELSE IF (io < 0) THEN
        EXIT
        ! Count number of atoms
     ELSE
        IF ((lign80(1:4).eq.'ATOM').or. &
             ((lign80(1:6).eq.'HETATM').and.hetatm)) THEN
           i=i+1
           READ(lign80,'(A6,I5,1X,A4,1X,A3,1X,A1,I5,3X,3F8.3,2F6.2,10X,2A2)') &
                atom(i),atnum(i),name(i),res(i),chain(i),resnum(i),x(i),y(i),z(i),&
                occ(i),bfac(i),elem(i),chag(i)
        END IF
     END IF
  END DO

  CLOSE(2356)

  WRITE(6,'(A,I5,A)') "Read in pdbfile, found", natom, " atoms"

!------------------------------------------------------------------------------------------------------------------

! Read in pdb file
  OPEN(file=pdbfile2,form="FORMATTED",status="OLD",unit=2356)

  i=0
  DO
     READ(2356,'(A)',IOSTAT=io) lign80
     IF (io > 0) THEN
        WRITE(*,*) 'Check input.  Something was wrong'
        STOP
        ! End of file
     ELSE IF (io < 0) THEN
        EXIT
        ! Count number of atoms
     ELSE
        IF ((lign80(1:4).eq.'ATOM').or. &
             ((lign80(1:6).eq.'HETATM').and.hetatm)) THEN
           i=i+1
        END IF
     END IF
  END DO
  
  natom2=i

  REWIND(2356)

  ALLOCATE(x2(natom2))
  ALLOCATE(y2(natom2))
  ALLOCATE(z2(natom2))
  ALLOCATE(atom2(natom2))
  ALLOCATE(atnum2(natom2))
  ALLOCATE(name2(natom2))
  ALLOCATE(res2(natom2))
  ALLOCATE(chain2(natom2))
  ALLOCATE(resnum2(natom2))
  ALLOCATE(occ2(natom2))
  ALLOCATE(bfac2(natom2))
  ALLOCATE(elem2(natom2))
  ALLOCATE(chag2(natom2))
  ALLOCATE(mass2(natom2))

  REWIND(2356)
 
  i=0
  DO
     READ(2356,'(A)',IOSTAT=io) lign80
     IF (io > 0) THEN
        WRITE(*,*) 'Check input.  Something was wrong'
        STOP
        ! End of file
     ELSE IF (io < 0) THEN
        EXIT
        ! Count number of atoms
     ELSE
        IF ((lign80(1:4).eq.'ATOM').or. &
             ((lign80(1:6).eq.'HETATM').and.hetatm)) THEN
           i=i+1
           READ(lign80,'(A6,I5,1X,A4,1X,A3,1X,A1,I5,3X,3F8.3,2F6.2,10X,2A2)') &
                atom2(i),atnum2(i),name2(i),res2(i),chain2(i),resnum2(i),x2(i),y2(i),z2(i),&
                occ2(i),bfac2(i),elem2(i),chag2(i)
        END IF
     END IF
  END DO

  CLOSE(2356)

  WRITE(6,'(A,I5,A)') "Read in pdbfile2, found", natom2, " atoms"

!-----------------------------------------------------------------------

  DO j=1,1000
     
     k=0
     DO i=1,natom
        IF (resnum(i).eq.j) THEN
           k=k+1
           fill(k)=name(i)
           pos(k)=i
        END IF
     END DO
     length1=k

     k=0
     DO i=1,natom2
        IF (resnum2(i).eq.j) THEN
           k=k+1
           fill2(k)=name2(i)
           pos2(k)=i
        END IF
     END DO
     length2=k

     DO i=1,length1
        DO k=1,length2
           IF (fill(i).eq.fill2(k)) THEN
              bfac(pos(i))=bfac2(pos2(k))
           END IF
        END DO
     END DO

     sum=0
     DO i=1,length2
        sum=sum+bfac2(pos2(i))
     END DO
     sum=sum/real(length2)

     DO i=1,length1
        IF (bfac(pos(i)).eq.0) THEN
           bfac(pos(i))=sum
        END IF
     END DO


  END DO

  
  OPEN(file='out.pdb',form='FORMATTED',unit=8765)

  DO i=1,natom
     WRITE(8765,'(A6,I5,1X,A4,1X,A3,1X,A1,I4,1X,3X,3F8.3,2F6.2,10X,2A2)') &
                atom(i),atnum(i),name(i),res(i),chain(i),resnum(i),x(i),y(i),z(i),&
                occ(i),bfac(i),elem(i),chag(i)
  END DO

  CLOSE(8765)

! End of reading in pdb file
!----------------------------------------------------------------------------------------- 



END PROGRAM bfac_input


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

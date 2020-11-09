PROGRAM allo_chain
  IMPLICIT NONE
  LOGICAL :: getoption, hetatm, caonly, resmass, atmass, cutvect,cutassign,forceres,ukn
  CHARACTER :: dummy*120, pdbfile*120, lign80*80, nomatm*120, res2*4, nomres*120
  INTEGER :: natom, io, i, j, natomold, nident, ncusres,ii,ijjjj,iseed,jat,jj,ll,nntr,nnzero
  REAL(8) :: cutoff, rave, rdev, rmin, rmax,anmp,ddf,rkh,dist,dist2,dmax,dmin,distave,drms,kij, &
       shift,rx,ry,rz,trace,random
  LOGICAL,ALLOCATABLE,DIMENSION(:) :: otherres
  INTEGER,ALLOCATABLE,DIMENSION(:) :: atnum,resnum,resnumold, resone, restwo
  REAL(8),ALLOCATABLE,DIMENSION(:) :: x,y,z,occ,bfac,vals,mass,massold,cutvalue,acutoff,kijcust
  CHARACTER,ALLOCATABLE,DIMENSION(:) :: atom*6,name*4,res*3,chain*1,elem*2,chag*2,chainold*1,ident*4
  REAL(8),ALLOCATABLE,DIMENSION(:,:) :: a

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



  DO i=1,natom
     IF (resnum(i).le.200) THEN
        chain(i)='A'
!     ELSE IF (resnum(i).le.507) THEN
!        chain(i)='B'
!     ELSE IF (resnum(i).le.508) THEN
!        chain(i)='A'
     ELSE
        chain(i)='B'
     END IF
  END DO

  
  OPEN(file='chain.pdb',form='FORMATTED',unit=8765)

  DO i=1,natom
     IF (resnum(i).lt.10000) THEN
        WRITE(8765,'(A6,I5,1X,A4,1X,A3,1X,A1,I4,1X,3X,3F8.3,2F6.2,10X,2A2)') &
                atom(i),atnum(i),name(i),res(i),chain(i),resnum(i),x(i),y(i),z(i),&
                occ(i),bfac(i),elem(i),chag(i)
     ELSE
        WRITE(8765,'(A6,I5,1X,A4,1X,A3,1X,A1,I5,3X,3F8.3,2F6.2,10X,2A2)') &
                atom(i),atnum(i),name(i),res(i),chain(i),resnum(i),x(i),y(i),z(i),&
                occ(i),bfac(i),elem(i),chag(i)
     END IF
  END DO

  CLOSE(8765)

! End of reading in pdb file
!----------------------------------------------------------------------------------------- 



END PROGRAM allo_chain


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

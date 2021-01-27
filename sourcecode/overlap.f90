! calculates the overlap of modes
! Tom Rodgers 13 04 2011
PROGRAM overlap
  USE read_files
  USE nrutil
  IMPLICIT NONE
  CHARACTER :: lign80*80, pdbfile*120, hetchoice*20, file1*120, file2*120, pdbfile2*40, dummy*120
  INTEGER :: endvec, startvec, ndim, vector, i, j, k, io, natom, m, n
  INTEGER :: nat, resid, residold, st, en, startvec2, endvec2, natom2, num1, num2
  REAL(DP), ALLOCATABLE,DIMENSION(:) :: crosscor, eigenval1, eigenval2
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: eigen,  eigen2, over, cumover
  REAL(DP) :: sum, dot, lengthi, lengthj
  LOGICAL :: hetatm,getoption

  IF (getoption('-help',.false.,dummy)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF

  IF (getoption('-d1',.false.,dummy)) THEN
     ndim=1
  ELSE
     ndim=3
  END IF

  ! First eigenvector file
  IF(.not.getoption('-i',.true.,file1)) THEN
     WRITE(0,*) 'Specify input file with -i'
     CALL helptext(0)
     CALL exit(1)
  END IF

  ! Start vector
  IF(getoption('-s1',.true.,dummy)) THEN
     READ(dummy,*) startvec
  ELSE
     startvec=1
     WRITE(0,*) 'No start vector was specified, calculating on vector 1'
  END IF
 
  ! End vector
  IF(getoption('-e1',.true.,dummy)) THEN
     READ(dummy,*) endvec
  ELSE
     endvec=startvec
  END IF

! Second eigenvector file
  IF(.not.getoption('-i2',.true.,file2)) THEN
     file2=file1
  END IF

  ! Second start vector
  IF(getoption('-s2',.true.,dummy)) THEN
     READ(dummy,*) startvec2
  ELSE
     startvec2=1
     WRITE(0,*) 'No start vector was specified, calculating on vector 1'
  END IF

  ! Second end vector
  IF(getoption('-e2',.true.,dummy)) THEN
     READ(dummy,*) endvec2
  ELSE
     endvec2=startvec2
  END IF

  IF (getoption('-het',.false.,dummy)) THEN
     hetatm=.true.
  ELSE
     hetatm=.false.
  END IF

  CALL read_eigenfacs(file1,natom,startvec,endvec,ndim,eigenval1,eigen,num1)

  CALL read_eigenfacs(file2,natom2,startvec2,endvec2,ndim,eigenval2,eigen2,num2)

! Make sure the number of atoms are the same
  IF (natom2.lt.natom) THEN
     WRITE(*,*) "Second vector has less atoms than the first,"
     WRITE(*,*) "only reading the first ", natom2, " of the first"
     WRITE(*,*) "molecule"
     natom=natom2
!     pdbfile=pdbfile2
  ELSE IF (natom2.gt.natom) THEN
     WRITE(*,*) "Second vector has more atoms than the first,"
     WRITE(*,*) "only reading the first ", natom, " of the second"
     WRITE(*,*) "molecule"
  ELSE
  END IF

  ALLOCATE(over(num1,num2))
  ALLOCATE(cumover(num1,num2))

  DO m=1,num1
     DO n=1,num2    
        sum=0
        DO i=1,natom
           dot=0
           lengthi=0
           lengthj=0
           DO k=1,ndim
              dot=dot+eigen((i-1)*ndim+k,m)*eigen2((i-1)*ndim+k,n)
              lengthi=lengthi+eigen((i-1)*ndim+k,m)**2
              lengthj=lengthj+eigen2((i-1)*ndim+k,n)**2
           END DO
           sum=sum+dot/(lengthi**0.5*lengthj**0.5)
        END DO 
        over(m,n)=dabs(sum)/REAL(natom)
     END DO
  END DO


  cumover=0
  DO i=1,num1
     DO j=1,num2
        DO m=1,i
           DO n=1,j
           cumover(i,j)=cumover(i,j)+over(m,n)**2
           END DO
        END DO
     END DO
  END DO
  
  cumover=cumover**0.5

  OPEN(file="summary.dat",form="FORMATTED",unit=5897) 
  OPEN(file="summary_cum.dat",form="FORMATTED",unit=5898)
  WRITE(5897,'(A1)') '#'
  WRITE(5898,'(A1)') '#'
  DO i=1,num1
     DO j=1,num2
        WRITE(5897,'(1X,I5,1X,I5,1X,G11.4)') i+startvec-1, j+startvec2-1, over(i,j)
        WRITE(5898,'(1X,I5,1X,I5,1X,G11.4)') i+startvec-1, j+startvec2-1, cumover(i,j)
     END DO
     WRITE(5897,'(1X)')
     WRITE(5898,'(1X)')
  END DO
  CLOSE(5897)
  CLOSE(5898)

  DEALLOCATE(eigen)
  DEALLOCATE(eigen2)
  DEALLOCATE(over)
  DEALLOCATE(cumover)

END PROGRAM overlap

SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
  WRITE(iunit,'(A)')"                           O  V  E  R  L  A  P                           "
  WRITE(iunit,'(A)')"                               VERSION 1.0                               "
  WRITE(iunit,'(A)')"                                                                         "  
  WRITE(iunit,'(A)')"                               Written by:                               "
  WRITE(iunit,'(A)')"                      Tom Rodgers and David Burnell                      "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program calculates the overlap between modes from either one set of "
  WRITE(iunit,'(A)')"normal modes or the overlap between modes of two different sets. It      "
  WRITE(iunit,'(A)')"requires eigenfacs formated input files which can be produced from the   "
  WRITE(iunit,'(A)')"included ENM code, GNM code, GroNM, GroED, or AmberED packages.          "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"overlap.dat contains the overlap of each atom and sum.dat contains the   "
  WRITE(iunit,'(A)')"average values for the whole system.                                     "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"     Usage:                                                              "
  WRITE(iunit,'(A)')"             overlap -i infile -s startvec [-e endvec]                   "
  WRITE(iunit,'(A)')"                     [-i2 infile2] [-s2 startvec2] [-d1] [-het]          "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"Option    Type      Value     Description                                "
  WRITE(iunit,'(A)')"------------------------------------------------------------             "
  WRITE(iunit,'(A)')"  -i      Input               First eigenvector file                     "
  WRITE(iunit,'(A)')"  -s1     Input               First included eigenvector number          "
  WRITE(iunit,'(A)')"  -e2     Input,Opt   -s      Last included eigenvector number           "
  WRITE(iunit,'(A)')"  -i2     Input,Opt   -i      Second eigenvector file                    "
  WRITE(iunit,'(A)')"  -s2     Input,Opt   1       First included eigenvector number from     "
  WRITE(iunit,'(A)')"                              second eigenvector file                    "
  WRITE(iunit,'(A)')"  -e2     Input,Opt   -s2     Last included eigenvector number from     "
  WRITE(iunit,'(A)')"                              second eigenvector file                    "
  WRITE(iunit,'(A)')"  -d1      Opt                Allows 1d vectors with flag                "
  WRITE(iunit,'(A)')"  -het     Opt                Incudes reading HETATM records with flag   "
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

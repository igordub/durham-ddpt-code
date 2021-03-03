! Program to calculate the eigenvalues and vectors
! Tom Rodgers 31 03 2011
PROGRAM diagstd
  USE nrtype
  USE utils
  USE write_files
  IMPLICIT NONE
  CHARACTER :: filename*64, lign80*80, slin*30, elin*26, dummy*120,outfile*120
  LOGICAL :: hetatm, getoption, qexist, covar
  REAL(DP),ALLOCATABLE,DIMENSION(:) :: d,e
  INTEGER,ALLOCATABLE,DIMENSION(:) :: g
  REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: a,b
  REAL(DP) :: f
  INTEGER :: i, j, io, natom, k

  IF (getoption('-help',.false.,dummy)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF
  
  IF(.not.getoption('-i',.true.,filename)) THEN
     filename="matrix.sdijf"
  END IF

  IF (getoption('-covar',.false.,dummy)) THEN
     covar=.true.
  ELSE
     covar=.false.
  END IF

  INQUIRE(file=filename,exist=qexist)
  IF (.not.qexist) THEN
     WRITE(6,'(A)') filename, "does not exist"
     STOP "Specify with -i flag if not matrix.sdijk"
  END IF

  outfile='matrix.eigenfacs'

!----------------------------------------------------------
  
  OPEN(file=filename,form='FORMATTED',status='OLD',unit=6897)
  
  natom=0

  DO
     READ(6897,'(1X,I9,1X,I9,1X,G19.11)',IOSTAT=io) i,j,f
     IF (io > 0) THEN
        WRITE(*,*) 'Check input.  Something was wrong'
        STOP
        ! End of file
     ELSE IF (io < 0) THEN
        EXIT
     ELSE
        natom=max(natom,i,j)
     ENDIF
  END DO
  
  natom=natom/3

  ALLOCATE(a(3*natom,3*natom))
  
  REWIND(6897)

  a=0

  DO
     READ(6897,'(A)',IOSTAT=io) lign80
     IF (io > 0) THEN
        WRITE(*,*) 'Check input.  Something was wrong'
        STOP
        ! End of file
     ELSE IF (io < 0) THEN
        EXIT
     ELSE
        READ(lign80,'(1X,I9,1X,I9,1X,G19.11)') i,j,a(i,j)
        a(j,i)=a(i,j)
     ENDIF
  END DO
  
  CLOSE(6897)

  ALLOCATE(d(3*natom))
  ALLOCATE(e(3*natom)) 
  
! Break matrix into tri-diagonal matrix
  CALL tred2(a,d,e) 
! a is replaced by orthoganal matrix Q
! d is the diagonal elements of the tri-diagonal matrix
! e is the off diagonal matrix

! Calculates the eigenvalues and eigenvectors
  CALL tqli(d,e,a)
! d is replaced as the eigenvalues
! a becomes the eigenvectors a(:,k) is for d(k)
! e is destroyed

! sort d and a in order of eigenvalue

  DEALLOCATE(e)


  IF (covar) THEN
     CALL sort_max(a,d,3*natom)
     DO i=1,3*natom
        IF (d(i).lt.1E-08) d(i)=1E-8
     END DO
     d=0.5961752/d
  ELSE
     CALL sort_min(a,d,3*natom)
  END IF
  
  DO i=1,3*natom
     IF (d(i).lt.1E-012) d(i)=1E-12
  END DO
  
  CALL write_eigenfacs(outfile,a,d,3*natom,3*natom)

DEALLOCATE(a)
DEALLOCATE(d)

END PROGRAM diagstd

SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
  WRITE(iunit,'(A)')"                           D  I  A  G  S  T  D                           "
  WRITE(iunit,'(A)')"                               VERSION 1.0                               "
  WRITE(iunit,'(A)')"                                                                         "  
  WRITE(iunit,'(A)')"                               Written by:                               "
  WRITE(iunit,'(A)')"                      Tom Rodgers and David Burnell                      "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program calculates the normal modes for a ENM based on the reduced  "
  WRITE(iunit,'(A)')"Hermitian matrix file, matrix.sdijf, produced from the ENM program.      "
  WRITE(iunit,'(A)')"A matrix.eigenfacs file is produced which can be further analysed using  "
  WRITE(iunit,'(A)')"the programs in this toolbox.                                            "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"     Usage:                                                              "
  WRITE(iunit,'(A)')"             diagstd [-i matrix]                                         "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"Option    Type       Value           Description                         "
  WRITE(iunit,'(A)')"------------------------------------------------------------             "
  WRITE(iunit,'(A)')"  -i      Input,Opt  matrix.sdijf    reduced Hermitian matrix file       "
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

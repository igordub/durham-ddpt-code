! calculates the cross correlation for each residue over the selected number
! of modes
! Tom Rodgers 08 04 2011
PROGRAM cross_cor
  USE read_files
  USE nrutil
  IMPLICIT NONE
  CHARACTER :: lign80*80, pdbfile*120, hetchoice*20, filename*120, dummy*120
  INTEGER :: skipvec, numvec, endvec, startvec, ndim, vector, i, j, k, io, natom, m
  INTEGER :: nat, resid, residold, st, en
  INTEGER, ALLOCATABLE,DIMENSION(:) :: natvec
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: eigenval
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: eigen, crosscor
  REAL(DP) :: sum, dot, lengthi, lengthj
  LOGICAL :: hetatm, getoption, qexist,nopdb
  INTEGER, ALLOCATABLE,DIMENSION(:) :: atnum,resnum
  REAL(DP), ALLOCATABLE,DIMENSION(:) ::x,y,z,occ,bfac
  CHARACTER, ALLOCATABLE, DIMENSION(:) :: atom*6,name*4,res*3,chain*1,elem*2,chag*2

  IF (getoption('-help',.false.,dummy)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF
  
  IF(.not.getoption('-i',.true.,filename)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF

  INQUIRE(file=filename,exist=qexist)
  IF (.not.qexist) THEN
     STOP "eigenfacs file specified does not exist"
  END IF

  nopdb=.false.
  IF(.not.getoption('-pdb',.true.,pdbfile)) THEN
     nopdb=.true.
  END IF

  INQUIRE(file=pdbfile,exist=qexist)
  IF (.not.qexist) THEN
     nopdb=.true.
  END IF
  
  IF (getoption('-d1',.false.,dummy)) THEN
     ndim=1
  ELSE
     ndim=3
  END IF
  
  ! Start vector
  IF(getoption('-s',.true.,dummy)) THEN
     READ(dummy,*) startvec
  ELSE
     startvec=7
  END IF
  skipvec=startvec-1

  ! End vector
  IF(getoption('-e',.true.,dummy)) THEN
     READ(dummy,*) endvec
  ELSE
     endvec=31
  END IF  
  numvec=endvec-skipvec

  IF (getoption('-het',.false.,dummy)) THEN
     hetatm=.true.
  ELSE
     hetatm=.false.
  END IF

  CALL read_eigenfacs(filename,natom,startvec,endvec,ndim,eigenval,eigen,i)

  print*, natom

  ALLOCATE(crosscor(natom,natom))

  DO i=1,natom
     DO j=1,natom
        sum=0
        DO m=1,numvec
           dot=0
           lengthi=0
           lengthj=0
           DO k=1,ndim
              dot=dot+eigen((i-1)*ndim+k,m)*eigen((j-1)*ndim+k,m)
              lengthi=lengthi+eigen((i-1)*ndim+k,m)**2
              lengthj=lengthj+eigen((j-1)*ndim+k,m)**2
           END DO
           sum=sum+dot/(lengthi**0.5*lengthj**0.5)
        END DO
        crosscor(i,j)=sum/numvec
     END DO
  END DO


  IF (nopdb) THEN
     ALLOCATE(resnum(natom))
     DO i=1,natom
        resnum(i)=i
     END DO
  ELSE
     CALL read_pdb(pdbfile,hetatm,natom,atom,atnum,name,res,chain,resnum,x,y,z, &
          occ,bfac,elem,chag) 
  END IF
  
!!$  OPEN(file=pdbfile,form="FORMATTED",status="OLD",unit=2356)
!!$
!!$  DO
!!$     READ(2356,'(A)',IOSTAT=io) lign80
!!$     IF (io > 0) THEN
!!$        WRITE(*,*) 'Check input.  Something was wrong'
!!$        STOP
!!$        ! End of file
!!$     ELSE IF (io < 0) THEN
!!$        EXIT
!!$        ! Count up number of atoms
!!$     ELSE
!!$        IF (lign80(1:4).eq.'ATOM') THEN
!!$           READ(lign80,'(22X,I4)') residold
!!$           EXIT
!!$        END IF
!!$     END IF
!!$  END DO
!!$
!!$  REWIND(2356)
!!$
!!$  nat=0
!!$  resid=0
!!$  i=0
!!$  j=1
!!$  ALLOCATE(natvec(natom))
!!$
!!$  DO
!!$     READ(2356,'(A)',IOSTAT=io) lign80
!!$     IF (io > 0) THEN
!!$        WRITE(*,*) 'Check input.  Something was wrong'
!!$        STOP
!!$        ! End of file
!!$     ELSE IF (io < 0) THEN
!!$        IF (lign80(1:4).eq.'ATOM') THEN
!!$           nat=nat+1
!!$           natvec(nat)=i
!!$        END IF
!!$        EXIT
!!$        ! Count up number of atoms
!!$     ELSE
!!$        IF (lign80(1:4).eq.'ATOM') THEN
!!$           READ(lign80,'(22X,I4)') resid
!!$           IF (resid.eq.residold) THEN
!!$              i=i+1
!!$           ELSE
!!$              nat=nat+1
!!$              natvec(nat)=i
!!$              i=1
!!$           END IF
!!$           residold=resid
!!$           j=0
!!$        ELSE IF ((lign80(1:6).eq.'HETATM').and.hetatm) THEN
!!$           j=j+1
!!$           IF (j.eq.1) THEN
!!$              nat=nat+1
!!$              natvec(nat)=i
!!$              i=1
!!$              residold=resid
!!$           END IF
!!$           nat=nat+1
!!$           natvec(nat)=1
!!$        ELSE
!!$        END IF
!!$     END IF
!!$  END DO
!!$  
!!$  CLOSE(2356)
!!$
!!$  IF (nat.lt.natom) THEN
!!$     WRITE(*,*) "Adjusting cross-correlation in terms of residues"
!!$
!!$     DO j=1,natom
!!$        sum=0
!!$        DO k=1,natvec(1)
!!$           sum=sum+crosscor(k,j)
!!$        END DO
!!$        crosscor(1,j)=sum/natvec(1)
!!$     END DO
!!$
!!$     DO i=2,nat
!!$        st=0
!!$        DO j=1,i-1
!!$           st=st+natvec(j)
!!$        END DO
!!$        en=st+natvec(i)
!!$        st=st+1
!!$        
!!$        DO j=1,natom
!!$           sum=0
!!$           DO k=st,en
!!$              sum=sum+crosscor(k,j)
!!$           END DO
!!$           crosscor(i,j)=sum/natvec(i)
!!$        END DO
!!$
!!$     END DO
!!$
!!$     DO j=1,nat
!!$        sum=0
!!$        DO k=1,natvec(1)
!!$           sum=sum+crosscor(j,k)
!!$        END DO
!!$        crosscor(j,1)=sum/natvec(1)
!!$     END DO
!!$
!!$     DO i=2,nat
!!$        st=0
!!$        DO j=1,i-1
!!$           st=st+natvec(j)
!!$        END DO
!!$        en=st+natvec(i)
!!$        st=st+1
!!$        
!!$        DO j=1,nat
!!$           sum=0
!!$           DO k=st,en
!!$              sum=sum+crosscor(j,k)
!!$           END DO
!!$           crosscor(j,i)=sum/natvec(i)
!!$        END DO
!!$
!!$     END DO
!!$  
!!$  END IF
  

!  OPEN(file="atom",form="FORMATTED",unit=5898)
!  WRITE(5898,'(I4)') nat
!  CLOSE(5898)


  OPEN(file="crosscor.dat",form="FORMATTED",unit=5897)

  DO i=1,natom
     DO j=1,natom
        WRITE(5897,'(1X,I5,1X,I5,1X,G11.4)') resnum(i), resnum(j), crosscor(i,j)
     END DO
     WRITE(5897,'(1X)')
  END DO

  CLOSE(5897)
  
!  DEALLOCATE(eigen)
  DEALLOCATE(crosscor)
!  DEALLOCATE(natvec)

END PROGRAM cross_cor

SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
  WRITE(iunit,'(A)')"                           C  R  O  S  C  O  R                           "
  WRITE(iunit,'(A)')"                               VERSION 1.0                               "
  WRITE(iunit,'(A)')"                                                                         "  
  WRITE(iunit,'(A)')"                               Written by:                               "
  WRITE(iunit,'(A)')"                      Tom Rodgers and David Burnell                      "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program calculates the average residue cross correlation from an    "
  WRITE(iunit,'(A)')"eigenfacs formated file which can be produced from the included ENM      "
  WRITE(iunit,'(A)')"code, GNM code, or the GroAMED packages.                                 "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"The output file, crosscor.dat can be plotted in gnuplot with the command:"
  WRITE(iunit,'(A)')'splot "crosscor.dat" using 1:2:3 notitle                                 '  
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"     Usage:                                                              "
  WRITE(iunit,'(A)')"             cross_cor -i infile [-s startvec]              "
  WRITE(iunit,'(A)')"                       [-e endvec] [-d1] [-het]                          "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"Option    Type       Value             Description                       "
  WRITE(iunit,'(A)')"------------------------------------------------------------             "
  WRITE(iunit,'(A)')"  -i       Input                       First eigenvector file            "
!  WRITE(iunit,'(A)')" -pdb      Input                       First eigenvector file            "
  WRITE(iunit,'(A)')"  -s      Input,Opt   7                First included eigenvector number "
  WRITE(iunit,'(A)')"  -e      Input,Opt   31               Last included eigenvector number  "
  WRITE(iunit,'(A)')"  -d1       Opt                        Allows 1d vectors with flag       "
  WRITE(iunit,'(A)')"  -het      Opt                        Incudes reading HETATM records with flag"
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

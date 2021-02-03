! Program to calculate the GNM
! Tom Rodgers 31 03 2011
PROGRAM GNM
  USE nrtype
  USE utils
  USE read_files
  IMPLICIT NONE
  CHARACTER :: filename*120, linerd*80, slin*30, elin*26, dummy*120, fname*120
  LOGICAL :: hetatm, getoption
  INTEGER, ALLOCATABLE,DIMENSION(:) :: atnum,resnum,g
  REAL(DP), ALLOCATABLE,DIMENSION(:) ::x,y,z,occ,bfac,d,e
  CHARACTER, ALLOCATABLE, DIMENSION(:) :: atom*6,name*4,res*3,chain*1,elem*2,chag*2
  REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: a,b
  REAL(DP) :: dist, cutoff, f, xot, yot, zot, force
  INTEGER :: i, j, io, natom, k

  IF (getoption('-help',.false.,dummy)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF
  
  IF(.not.getoption('-pdb',.true.,filename)) THEN
     WRITE(0,'(A)') "Need to input a pdb file with -pdb"
     CALL helptext(0)
     CALL exit(0)
  END IF

  IF(getoption('-c',.true.,dummy)) THEN
     READ(dummy,*) cutoff
  ELSE
     cutoff=12
  END IF

  IF(getoption('-f',.true.,dummy)) THEN
     READ(dummy,*) force
  ELSE
     force=1
  END IF

  IF (getoption('-het',.false.,dummy)) THEN
     hetatm=.true.
  ELSE
     hetatm=.false.
  END IF
  


  CALL read_pdb(filename,hetatm,natom,atom,atnum,name,res,chain,resnum,x,y,z, &
                  occ,bfac,elem,chag)

  j=0
  DO i=1,natom
     IF ((name(i).eq." CA ").or.(atom(i).eq."HETATM".and.hetatm)) THEN
        j=j+1
        atom(j)=atom(i)
        atnum(j)=atnum(i)
        name(j)=name(i)
        res(j)=res(i)
        chain(j)=chain(i)
        resnum(j)=resnum(i)
        x(j)=x(i)
        y(j)=y(i)
        z(j)=z(i)
        occ(j)=occ(i)
        bfac(j)=bfac(i)
        elem(j)=elem(i)
        chag(j)=chag(i)
     END IF
  END DO
  natom=j


  OPEN(file="CAonly.pdb",form="FORMATTED",unit=9873)
  OPEN(file="flat.pdb",form="FORMATTED",unit=9874)
  
  DO i=1,natom
     WRITE(9873,'(A6,I5,1X,A4,1X,A3,1X,A1,I4,1X,3X,3F8.3,2F6.2,10X,2A2)') &
          atom(i),atnum(i),name(i),res(i),chain(i),resnum(i),x(i),0.d0,z(i),&
          occ(i),bfac(i),elem(i),chag(i)
     WRITE(9874,'(A6,I5,1X,A4,1X,A3,1X,A1,I4,1X,3X,3F8.3,2F6.2,10X,2A2)') &
          atom(i),atnum(i),name(i),res(i),chain(i),resnum(i),x(i),0.d0,z(i),&
          occ(i),bfac(i),elem(i),chag(i)
  END DO

  CLOSE(9873)
  CLOSE(9873)

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

  ALLOCATE(a(natom,natom))

  OPEN(file="GNM.vmd",form="FORMATTED",unit=5978)

  WRITE(5978,'(A)') '#!/usr/local/bin/vmd'
  WRITE(5978,'(A)') '# script for VMD (Visual Molecular Dynamics)'
  WRITE(5978,'(A)') '# Goal: visualizing the elastic network'
  WRITE(5978,'(A)') '# Type: vmd -e this-file'
  WRITE(5978,'(A)') 'color Display {Background} white'
  WRITE(5978,'(A)') 'mol new'
  WRITE(5978,'(A)') 'draw color black'


  a=0
  DO i=1,natom
     DO j=1,natom
        IF (i.ne.j) THEN
           dist=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
           
           IF (dist.le.cutoff) THEN
              a(i,j)=-1.000
              WRITE(5978,'(A,3F12.4,A,3F12.4,A)') 'draw line {',x(i),0.0,z(i),'} {',x(j),0.0,z(j),'}'
           END IF
        END IF
     END DO
  END DO

  WRITE(5978,'(A)') 'mol load pdb flat.pdb'
  CLOSE(5978)


  ALLOCATE(b(natom,natom))

  b=a

  DO i=1,natom
     a(i,i)=-sum(b(i,:))
  END DO

  DEALLOCATE(b)
  
  a=force*a

  OPEN(file="pdbmat.sdijf",unit=6897)
  
  DO i=1,natom
     DO j=i,natom
        IF (a(i,j).ne.0) THEN
           WRITE(6897,'(1X,I9,1X,I9,1X,G19.11)') i,j,a(i,j)
        END IF
     END DO
END DO
  
CLOSE(6897)

  ALLOCATE(d(natom))
  ALLOCATE(e(natom)) 
  

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

  CALL sort_min(a,d,natom)

  fname='matrix.eigenfacs'

  OPEN(file=fname,form="FORMATTED",unit=4679)
  
  DO k=1,natom
     
     WRITE(4679,'(1X,A6,I5,7X,A5,1PG12.4)') "VECTOR", k, "VALUE", d(k) 
     WRITE(4679,'(A36)') " -----------------------------------"
     
     DO i=1,natom
        WRITE(4679,'(1PG12.4)') a(i,k)
     END DO
     
  END DO

  CLOSE(4679)

DEALLOCATE(x)
DEALLOCATE(y)
DEALLOCATE(z)
DEALLOCATE(a)
DEALLOCATE(d)

END PROGRAM GNM

SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
  WRITE(iunit,'(A)')"                           G  N  M  P  R  O  D                           "
  WRITE(iunit,'(A)')"                               VERSION 1.0                               "
  WRITE(iunit,'(A)')"                                                                         "  
  WRITE(iunit,'(A)')"                               Written by:                               "
  WRITE(iunit,'(A)')"                      Tom Rodgers and David Burnell                      "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program calculates the normal modes for a GNM based on the pdb      "
  WRITE(iunit,'(A)')"inputed with -pdb. A matrix.eigenfacs file is produced which can be      "
  WRITE(iunit,'(A)')"further analysed using the programs in this toolbox.                     "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program also outputs a map of the spring connections, GNM.vmd which "
  WRITE(iunit,'(A)')"can be viewed using VMD, a representation of the pdb in 2d, flat.pdb, and"
  WRITE(iunit,'(A)')"a pdb of the Ca atoms only, CAonly.pdb.                                  "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"     Usage:                                                              "
  WRITE(iunit,'(A)')"             GNM -pdb pdbfile [-c cut-off] [-f force] [-het]             "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"Option    Type       Value           Description                         "
  WRITE(iunit,'(A)')"------------------------------------------------------------             "
  WRITE(iunit,'(A)')" -pdb     Input                      pdb file name                       "
  WRITE(iunit,'(A)')"  -c      Input,Opt  12              Cut-off for spring connectivity     "
  WRITE(iunit,'(A)')"  -f      Input,Opt  1               Spring constant                     "
  WRITE(iunit,'(A)')" -het      Opt                       Include HETATM references with flag "
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

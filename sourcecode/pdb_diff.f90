PROGRAM pdb_diff
  USE nrtype
  USE read_files
  USE write_files
  IMPLICIT NONE
  INTEGER, ALLOCATABLE,DIMENSION(:) :: atnum,resnum
  REAL(DP), ALLOCATABLE,DIMENSION(:) ::x1,y1,z1,x2,y2,z2,occ,bfac,eigenval,mass
  REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: eigenvec
  REAL(DP) :: tot
  CHARACTER, ALLOCATABLE, DIMENSION(:) :: atom*6,name*4,res*3,chain*1,elem*2,chag*2
  INTEGER :: i, natom
  CHARACTER :: pdbfile1*120, dummy*120, pdbfile2*120, filename*120
  LOGICAL :: hetatm, getoption

  IF (getoption('-help',.false.,dummy)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF

  IF (.not.getoption('-pdb1',.true.,pdbfile1)) THEN
     WRITE(0,'(A)') "Need to input a pdb file with -pdb1"
     CALL helptext(0)
     CALL exit(0)
  END IF

  IF (.not.getoption('-pdb2',.true.,pdbfile2)) THEN
     WRITE(0,'(A)') "Need to input a pdb file with -pdb2"
     CALL helptext(0)
     CALL exit(0)
  END IF

  IF(getoption('-het',.false.,dummy)) THEN
     hetatm=.true.
  ELSE
     hetatm=.false.
  END IF 

  CALL read_pdb(pdbfile1,hetatm,natom,atom,atnum,name,res,chain,resnum,x1,y1,z1, &
       occ,bfac,elem,chag)

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
  
  CALL read_pdb(pdbfile2,hetatm,natom,atom,atnum,name,res,chain,resnum,x2,y2,z2, &
       occ,bfac,elem,chag)

  ALLOCATE(eigenvec(3*natom,1))
  ALLOCATE(eigenval(1))

  DO i=1,natom
     eigenvec(3*(i-1)+1,1)=x1(i)-x2(i)
     eigenvec(3*(i-1)+2,1)=y1(i)-y2(i)
     eigenvec(3*(i-1)+3,1)=z1(i)-z2(i)
  END DO

  tot=0
  DO i=1,3*natom
     tot=tot+eigenvec(i,1)**2.d0
  END DO
  tot=tot**0.5d0

  eigenvec=eigenvec/tot
  eigenval(1)=1/tot**2.d0

  filename='pdb.eigenfacs'
  CALL write_eigenfacs(filename,eigenvec,eigenval,1,3*natom)

END PROGRAM pdb_diff

SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
  WRITE(iunit,'(A)')"                           P  D  B  D  I  F  F                           "
  WRITE(iunit,'(A)')"                               VERSION 1.0                               "
  WRITE(iunit,'(A)')"                                                                         "  
  WRITE(iunit,'(A)')"                               Written by:                               "
  WRITE(iunit,'(A)')"                      Tom Rodgers and David Burnell                      "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program projects the eigenvectors onto the pdb file                 "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"     Usage:                                                              "
  WRITE(iunit,'(A)')"             genENM -pdb pdbfile [-i eigenfacsfile] [-s start]           "
  WRITE(iunit,'(A)')"                    [-e end] [-scale n] [-het]                           " 
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"Option    Type       Value                Description                    "
  WRITE(iunit,'(A)')"------------------------------------------------------------             "
  WRITE(iunit,'(A)')" -pdb     Input                           pdb file name                  "
  WRITE(iunit,'(A)')"  -i    Input,Opt    matrix.eigenfacs     eigenvectors file              "
  WRITE(iunit,'(A)')"  -s    Input,Opt    7                    first eigenvector included     "
  WRITE(iunit,'(A)')"  -e    Input,Opt    31                   last eigenvector included      "
  WRITE(iunit,'(A)')"-scale  Input,Opt    1                    scale factor                   "
  WRITE(iunit,'(A)')" -het      Opt                            include HETATMs                "
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

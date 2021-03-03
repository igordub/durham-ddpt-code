PROGRAM domains
  USE nrutil
  USE read_files
  USE write_files
  USE utils
  IMPLICIT NONE
  INTEGER, ALLOCATABLE,DIMENSION(:) :: atnum, resnum, sortarray, domb, atnumdom
  REAL(DP), ALLOCATABLE,DIMENSION(:) ::x, y, z, occ, bfac, mass
  CHARACTER, ALLOCATABLE, DIMENSION(:) :: atom*6, name*4, res*3, chain*1, &
       elem*2, chag*2, cdom*1, chaindom*1
  INTEGER :: io, i, natom, ndom, j, num, natom1
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: dom
  CHARACTER :: pdbfile*120, lign80*80, dummy*120, fname*120, pdb2*120, &
       feigname*120, filename*120
  LOGICAL :: hetatm, getoption, rev, qexist
  REAL(8),ALLOCATABLE,DIMENSION(:,:) :: eigenvec, eigenvechold
  REAL(8),ALLOCATABLE,DIMENSION(:) :: eigenval
    
  IF (getoption('-help',.false.,dummy)) THEN
     CALL helptext(0)
     CALL EXIT(0)
  END IF
  
  IF (.not.getoption('-pdb',.true.,pdbfile)) THEN
     WRITE(6,'(A)') "Need to assign pdb file with -pdb flag"
     CALL helptext(0)
     CALL EXIT(0) 
  END IF
  
  INQUIRE(file=pdbfile,exist=qexist)
  IF (.not.qexist) THEN
     STOP "pdb file specified does not exist"
  END IF

  IF (.not.getoption('-pdb2',.true.,pdb2)) THEN
     pdb2='domains.pdb'
  END IF

  IF (.not.getoption('-i',.true.,filename)) THEN
     filename='matrix.eigenfacs'
  END IF

  IF (getoption('-het',.false.,dummy)) THEN
     hetatm=.true.
  ELSE
     hetatm=.false.
  END IF

  IF (getoption('-rev',.false.,dummy)) THEN
     rev=.true.
  ELSE
     rev=.false.
  END IF

  CALL read_pdb(pdbfile,hetatm,natom,atom,atnum,name,res,chain,resnum,x,y,z, &
       occ,bfac,elem,chag)
  
  ALLOCATE(mass(natom))
  mass=0

  IF (rev) THEN
     DEALLOCATE(atom)
     DEALLOCATE(name)
     DEALLOCATE(res)
     DEALLOCATE(resnum)
     DEALLOCATE(x)
     DEALLOCATE(y)
     DEALLOCATE(z)
     DEALLOCATE(occ)
     DEALLOCATE(bfac)
     DEALLOCATE(elem)
     DEALLOCATE(chag)
     
     INQUIRE(file=pdb2,exist=qexist)
     IF (.not.qexist) THEN
        WRITE(6,'(A)') pdb2, "does not exist"
        STOP "specify with -pdb2 flag"
     END IF

     CALL read_pdb(pdb2,hetatm,natom,atom,atnumdom,name,res,chaindom,resnum,x,y,z, &
          occ,bfac,elem,chag)
     natom1=natom
     INQUIRE(file=filename,exist=qexist)
     IF (.not.qexist) THEN
        WRITE(6,'(A)') filename, "does not exist"
        STOP "specify with -i flag"
     END IF

     CALL read_eigenfacs(filename,natom,1,10000,3,eigenval,eigenvec,num)
     
     IF (natom1.lt.natom) THEN
          WRITE(6,'(A)') "eigenfac file contains more atoms than the pdb file"
          WRITE(6,'(A)') "Maybe -het was used or the vectors are not Ca and pdb is"
          STOP
       ELSE IF (natom1.gt.natom) THEN
          WRITE(6,'(A)') "eigenfac file contains less atoms than the pdb file"
          WRITE(6,'(A)') "Maybe the vectors are Ca only and pdb isn't"
          STOP
       END IF

     ALLOCATE(eigenvechold(natom*3,num))
     eigenvechold=eigenvec
     
     DO i=1,natom
        DO j=1,natom
           IF (atnum(i).eq.atnumdom(j).and.chain(i).eq.chaindom(j)) THEN
              eigenvec((i-1)*3+1,:)=eigenvechold((j-1)*3+1,:)
              eigenvec((i-1)*3+2,:)=eigenvechold((j-1)*3+2,:)
              eigenvec((i-1)*3+3,:)=eigenvechold((j-1)*3+3,:)
           END IF
        END DO
     END DO

     feigname='matrix.eigenfacssort'

     CALL write_eigenfacs(feigname,eigenvec,eigenval,num,natom*3)

  ELSE
     INQUIRE(file='input.domain',exist=qexist)
     IF (.not.qexist) THEN
        STOP "Need an input.domain file to allocate domains"
     END IF
     CALL get_domain(ndom,dom,cdom,domb)

     ALLOCATE(sortarray(natom))
     CALL sort_array(sortarray,natom,ndom,chain,cdom,resnum,dom)
     
     CALL run_sort_c(atom,natom,sortarray)
     CALL run_sort_i(atnum,natom,sortarray)
     CALL run_sort_c(name,natom,sortarray)
     CALL run_sort_c(res,natom,sortarray)
     CALL run_sort_c(chain,natom,sortarray)
     CALL run_sort_i(resnum,natom,sortarray)
     CALL run_sort_r(x,natom,sortarray)
     CALL run_sort_r(y,natom,sortarray)
     CALL run_sort_r(z,natom,sortarray)
     CALL run_sort_r(occ,natom,sortarray)
     CALL run_sort_r(bfac,natom,sortarray)
     CALL run_sort_c(elem,natom,sortarray)
     CALL run_sort_c(chag,natom,sortarray)
     
     fname='domains.pdb'
     CALL write_pdb(fname,atom,atnum,name,res,chain,resnum,x,y,z,&
          occ,bfac,mass,elem,chag,natom,.false.,hetatm)

  END IF

END PROGRAM domains

SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
  WRITE(iunit,'(A)')"                           D  O  M  A  I  N  S                           "
  WRITE(iunit,'(A)')"                               VERSION 1.0                               "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"                               Written by:                               "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program rewrites a pdbfile in order of the domains listed in        "
  WRITE(iunit,'(A)')"input.domain for use with rtb domains.                                   "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program also rewrites the domain eigenfacs file back into the       "
  WRITE(iunit,'(A)')"pdb order.                                                               " 
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"     Usage:                                                              "
  WRITE(iunit,'(A)')"             domains -pdb pdbfile [-i matrix]                            "
  WRITE(iunit,'(A)') "                                                                        "
  WRITE(iunit,'(A)')"Option    Type       Value             Description                       "
  WRITE(iunit,'(A)')"------------------------------------------------------------             "
  WRITE(iunit,'(A)')" -pdb      Input                        pdb file                         "
  WRITE(iunit,'(A)')" -pdb2   Input,Opt  domains.pdb       Domain ordered pdb file            "
  WRITE(iunit,'(A)')"  -i     Input,Opt  matrix.eigenfacs  Domain ordered eigenfacs file      "
  WRITE(iunit,'(A)')" -het       Opt                       Include reading HETATM records     "
  WRITE(iunit,'(A)')" -rev       Opt                       Sort eigenfacs file instead of pdb "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')" input.domain format:                                                    "
  WRITE(iunit,'(A)')"   WRITE('input.domain','(I4,1X,A1,1X,I4,1X,I4)') domb, cdom, doms, dome "
  WRITE(iunit,'(A)')"   atomname is in pdb format                                             "  
  WRITE(iunit,'(A)')"   e.g.                                                                  "  
  WRITE(iunit,'(A)')"   |   1 A   10   15                                                     "
  WRITE(iunit,'(A)')"   |   2 B   20   30                                                     "
  WRITE(iunit,'(A)')"                                                                         "

END SUBROUTINE helptext

FUNCTION getoption(flag,getval,cvalue)
  IMPLICIT NONE
  CHARACTER(*):: flag,cvalue
  CHARACTER(160):: arg
  LOGICAL ::getoption,getval
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
  IF (getval) CALL getarg(i+1,cvalue)
END FUNCTION getoption


PROGRAM full_to_ca
  USE nrutil
  USE read_files
  USE write_files
  IMPLICIT NONE
 
    INTEGER, ALLOCATABLE,DIMENSION(:) :: atnum,resnum, inclu
    REAL(DP), ALLOCATABLE,DIMENSION(:) ::x,y,z,occ,bfac,mass, eigenval
    REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: eigenvec
    CHARACTER, ALLOCATABLE, DIMENSION(:) :: atom*6,name*4,res*3,chain*1,elem*2,chag*2
    INTEGER :: io, i, natom, j, ndim, vector, k, startvec, endvec,num, natom1
    CHARACTER :: pdbfile*120,filename*120,outfile*120, dummy*120, filebfac*120, lign80*80
    LOGICAL :: hetatm, eigenvecs, getoption, bfacs, qexist
    REAL(DP) :: ei

    IF (.not.getoption('-pdb',.true.,pdbfile)) THEN
       WRITE(0,'(A)') "Need to input a pdb file with -pdb"
       CALL helptext(0)
       CALL exit(0)
    END IF
    
    INQUIRE(file=pdbfile,exist=qexist)
    IF (.not.qexist) THEN
       STOP "pdb file specified does not exist"
    END IF
    
    IF (getoption('-het',.false.,dummy)) THEN
       hetatm=.true.
    ELSE
       hetatm=.false.
    END IF

    IF (getoption('-i',.true.,filename)) THEN
       eigenvecs=.true.
       INQUIRE(file=filename,exist=qexist)
       IF (.not.qexist) THEN
          STOP "eigenfacs file specified does not exist"
       END IF
    ELSE
       eigenvecs=.false.
    END IF

    IF (getoption('-bfac',.true.,filebfac)) THEN
       bfacs=.true.
       INQUIRE(file=filebfac,exist=qexist)
       IF (.not.qexist) THEN
          STOP "bfactor file specified does not exist"
       END IF
    ELSE
       bfacs=.false.
    END IF

    CALL read_pdb(pdbfile,hetatm,natom,atom,atnum,name,res,chain,resnum,x,y,z, &
         occ,bfac,elem,chag)
    
    outfile='CAonly.pdb'

    ALLOCATE(mass(natom))
    mass=0
    
    CALL write_pdb(outfile,atom,atnum,name,res,chain,resnum,x,y,z,&
         occ,bfac,mass,elem,chag,natom,.true.,hetatm)
    natom1=natom

    ALLOCATE(inclu(natom))

    DO i=1,natom
       IF ((index(name(i),'CA').gt.0) .or. (atom(i).eq.'HETATM'.and.hetatm))THEN
          inclu(i)=1
       ELSE
          inclu(i)=0
       END IF
    END DO

    IF (eigenvecs) THEN
        
       CALL read_eigenfacs(filename,natom,1,3*natom,3,eigenval,eigenvec,num)
 
       IF (natom1.lt.natom) THEN
          WRITE(6,'(A)') "eigenfac file contains more atoms than the pdb file"
          WRITE(6,'(A)') "Maybe -het was used or the vectors are not Ca and pdb is"
          STOP
       ELSE IF (natom1.gt.natom) THEN
          WRITE(6,'(A)') "eigenfac file contains less atoms than the pdb file"
          WRITE(6,'(A)') "Maybe the vectors are Ca only and pdb isn't"
          STOP
       END IF

       j=0
       DO i=1,natom
          IF (inclu(i).eq.1) THEN
             j=j+1
             eigenvec((j-1)*3+1,:)=eigenvec((i-1)*3+1,:)
             eigenvec((j-1)*3+2,:)=eigenvec((i-1)*3+2,:)
             eigenvec((j-1)*3+3,:)=eigenvec((i-1)*3+3,:)
          END IF
       END DO
       
       outfile='matrix.eigenfacsca'
       CALL write_eigenfacs(outfile,eigenvec,eigenval,num,3*j)
       
    END IF

    IF (bfacs) THEN
       
       OPEN(file=filebfac,form='FORMATTED',status='OLD',unit=8134)
       OPEN(file='mode.bfacca',form='FORMATTED',unit=8135)

       READ(8134,'(A)') lign80
       WRITE(8135,'(A)') lign80
       READ(8134,'(A)') lign80
       WRITE(8135,'(A)') lign80

       DO i=1,natom
          READ(8134,'(A)') lign80
          IF (inclu(i).eq.1) THEN
             WRITE(8135,'(A)') lign80
          END IF
       END DO
       
       CLOSE(8134)

    END IF

END PROGRAM full_to_ca

SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
  WRITE(iunit,'(A)')"                           F  U  L  L  2  C  A                           "
  WRITE(iunit,'(A)')"                               VERSION 1.0                               "
  WRITE(iunit,'(A)')"                                                                         "  
  WRITE(iunit,'(A)')"                               Written by:                               "
  WRITE(iunit,'(A)')"                      Tom Rodgers and David Burnell                      "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program converts any imput files to the versions with only the      "
  WRITE(iunit,'(A)')"alpha Carbon atoms included.                                             "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"     Usage:                                                              "
  WRITE(iunit,'(A)')"             genENM -pdb pdbfile [-c cut-off] [-f force] [-ccust cfile]  "
  WRITE(iunit,'(A)')"                    [-fcust ffile] [-mass] [-ca] [-res] [-het] [-hin h]  "
  WRITE(iunit,'(A)')"                    [-an p]                                              " 
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"Option    Type       Value       Description                             "
  WRITE(iunit,'(A)')"------------------------------------------------------------             "
  WRITE(iunit,'(A)')" -pdb     Input                  pdb file name                           "
  WRITE(iunit,'(A)')"  -i    Input,Opt                Convert the eigenfacs style file        "
  WRITE(iunit,'(A)')" -bfac  Input,Opt                Convert the mode.bfactors file          "
  WRITE(iunit,'(A)')" -het      Opt                   Includes reading HETATM records         "
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

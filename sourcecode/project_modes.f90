PROGRAM project_modes
  USE nrtype
  USE read_files
  USE write_files
  IMPLICIT NONE
  INTEGER, ALLOCATABLE,DIMENSION(:) :: atnum,resnum
  REAL(DP), ALLOCATABLE,DIMENSION(:) ::x,y,z,occ,bfac,eigenval,xm,ym,zm, mass
  REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: eigenvec
  CHARACTER, ALLOCATABLE, DIMENSION(:) :: atom*6,name*4,res*3,chain*1,elem*2,chag*2
  INTEGER :: i, natom, startvec, endvec, num, j
  CHARACTER :: pdbfile*120, dummy*120, filename*120, fname*120
  LOGICAL :: hetatm, getoption
  REAL(DP) :: scale
  
  IF (getoption('-help',.false.,dummy)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF

  IF (.not.getoption('-pdb',.true.,pdbfile)) THEN
     WRITE(0,'(A)') "Need to input a pdb file with -pdb"
     CALL helptext(0)
     CALL exit(0)
  END IF

  IF (.not.getoption('-i',.true.,filename)) THEN
     filename='matrix.eigenfacs'
  END IF

  ! Start vector
  IF(getoption('-s',.true.,dummy)) THEN
     READ(dummy,*) startvec
  ELSE
     startvec=7
  END IF

  ! End vector
  IF(getoption('-e',.true.,dummy)) THEN
     READ(dummy,*) endvec
  ELSE
     endvec=31
  END IF  

  IF(getoption('-scale',.true.,dummy)) THEN
     READ(dummy,*) scale
  ELSE
     scale=1
  END IF  

  IF(getoption('-het',.false.,dummy)) THEN
     hetatm=.true.
  ELSE
     hetatm=.false.
  END IF 

  CALL read_pdb(pdbfile,hetatm,natom,atom,atnum,name,res,chain,resnum,x,y,z, &
       occ,bfac,elem,chag) 

  CALL read_eigenfacs(filename,natom,startvec,endvec,3,eigenval,eigenvec,num)

  ALLOCATE(xm(natom))
  ALLOCATE(ym(natom))
  ALLOCATE(zm(natom))
  ALLOCATE(mass(natom))
  mass=0.0d0

  DO i=1,num
     DO j=1,natom
        xm(j)=x(j)+(1/eigenval(i)**0.5d0)*scale*eigenvec((j-1)*3+1,i)
        ym(j)=y(j)+(1/eigenval(i)**0.5d0)*scale*eigenvec((j-1)*3+2,i)
        zm(j)=z(j)+(1/eigenval(i)**0.5d0)*scale*eigenvec((j-1)*3+3,i)
     END DO
     WRITE(fname,'(A5,I3.3,A4)') "Mode_", i-1+startvec,".pdb"
     CALL write_pdb(fname,atom,atnum,name,res,chain,resnum,xm,ym,zm,&
           occ,bfac,mass,elem,chag,natom,.false.,hetatm)
  END DO

END PROGRAM project_modes

SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
  WRITE(iunit,'(A)')"                           P  R  O  J  E  C  T                           "
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

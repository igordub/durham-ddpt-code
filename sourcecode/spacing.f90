PROGRAM spacing
  USE nrtype
  USE read_files
  IMPLICIT NONE
  INTEGER, ALLOCATABLE,DIMENSION(:) :: atnum,resnum,resnumold
  REAL(DP), ALLOCATABLE,DIMENSION(:) ::x,y,z,occ,bfac
  REAL(DP) :: masstol, xav, yav, zav, bfacav
  CHARACTER, ALLOCATABLE, DIMENSION(:) :: atom*6,name*4,res*3,chain*1,chainold*1,elem*2,chag*2
  INTEGER :: io, i, natom, natomold, j, fat
  CHARACTER :: pdbfile*120, lign80*80, dummy*120
  LOGICAL :: hetatm,getoption,caonly,lig1,qexist,firstat
  
  IF (getoption('-help',.false.,dummy)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF
  
  IF (.not.getoption('-pdb',.true.,pdbfile)) THEN
     WRITE(0,'(A)') "Need to input a pdb file with -pdb"
     CALL helptext(0)
     CALL exit(0)
  END IF
  
  INQUIRE(file=pdbfile,exist=qexist)
  IF (.not.qexist) THEN
     STOP "pdb file does not exist"
  END IF

  IF (getoption('-het',.false.,dummy)) THEN
     hetatm=.true.
  ELSE
     hetatm=.false.
  END IF
  
  IF (getoption('-ca',.false.,dummy)) THEN
     caonly=.true.
  ELSE
     caonly=.false.
  END IF
  
  IF (getoption('-lig1',.false.,dummy)) THEN
     lig1=.true.
  ELSE
     lig1=.false.
  END IF 

  CALL read_pdb(pdbfile,hetatm,natom,atom,atnum,name,res,chain,resnum,x,y,z, &
       occ,bfac,elem,chag)

  IF (caonly) THEN
     firstat=.true.
     ALLOCATE(resnumold(natom))
     ALLOCATE(chainold(natom))
     j=0
     DO i=1,natom
        resnumold(i)=resnum(i)
        chainold(i)=chain(i)
        IF (name(i).eq." CA ") THEN
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
        ELSE IF (atom(i).eq."HETATM".and.hetatm) THEN
           ! Make each HETATM residue 1 average bead
           IF (lig1) THEN
              IF (atom(i-1).ne.'HETATM') firstat=.true.
              IF (firstat) THEN
                 j=j+1
                 xav=x(i)
                 yav=y(i)
                 zav=z(i)
                 bfacav=bfac(i)
                 masstol=1
                 fat=i
                 firstat=.false.
              ELSE
                 IF (resnum(i).ne.resnum(i-1)) THEN
                    atom(j)=atom(i-1)
                    atnum(j)=atnum(i-1)
                    name(j)=' CA '
                    res(j)=res(i-1)
                    chain(j)=chain(i-1)
                    resnum(j)=resnum(i-1)
                    x(j)=xav/masstol
                    y(j)=yav/masstol
                    z(j)=zav/masstol
                    occ(j)=1.d0
                    bfac(j)=bfacav/masstol
                    elem(j)=' C'
                    chag(j)='  '
                                        
                    j=j+1
                    xav=x(i)
                    yav=y(i)
                    zav=z(i)
                    bfacav=bfac(i)
                    masstol=1
                    fat=i
                 ELSE
                    xav=xav+x(i)
                    yav=yav+y(i)
                    zav=zav+z(i)
                    bfacav=bfacav+bfac(i)
                    masstol=masstol+1
                 END IF
              END IF
              IF (i.eq.natom) THEN
                 atom(j)=atom(i-1)
                 atnum(j)=atnum(i-1)
                 name(j)=' CA '
                 res(j)=res(i-1)
                 chain(j)=chain(i-1)
                 resnum(j)=resnum(i-1)
                 x(j)=xav/masstol
                 y(j)=yav/masstol
                 z(j)=zav/masstol
                 occ(j)=1.d0
                 bfac(j)=bfacav/masstol
                 elem(j)=' C'
                 chag(j)='  '
              END IF
              ! print all HETATMs
           ELSE
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
           ! catch the print of the average HETATM if HETATMs are not all at end
        ELSE
           IF (atom(i-1).eq."HETATM".and.hetatm.and.lig1) THEN
              atom(j)=atom(i-1)
              atnum(j)=atnum(i-1)
              name(j)=' CA '
              res(j)=res(i-1)
              chain(j)=chain(i-1)
              resnum(j)=resnum(i-1)
              x(j)=xav/masstol
              y(j)=yav/masstol
              z(j)=zav/masstol
              occ(j)=1.d0
              bfac(j)=bfacav/masstol
              elem(j)=' C'
              chag(j)='  '
           END IF
        END IF
     END DO
     
     natomold=natom
     natom=j
     
     WRITE(6,'(A,I5,A)') "Adjusted to Ca atoms only, found", natom, " atoms"

     
     DEALLOCATE(resnumold)
     DEALLOCATE(chainold)
     
  END IF
  
  OPEN(file='dist.dat',form='FORMATTED',unit=6434)
  DO i=1,natom
     DO j=1,natom
        WRITE(6434,'(1X,I5,1X,I5,1X,G11.4)') resnum(i), resnum(j), ((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)**0.5
     END DO
     WRITE(6434,'(1X)')
  END DO
  CLOSE(6434)



END PROGRAM spacing


SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
  WRITE(iunit,'(A)')"                           S  P  A  C  I  N  G                           "
  WRITE(iunit,'(A)')"                               VERSION 1.0                               "
  WRITE(iunit,'(A)')"                                                                         "  
  WRITE(iunit,'(A)')"                               Written by:                               "
  WRITE(iunit,'(A)')"                      Tom Rodgers and David Burnell                      "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program produces a file, dist.dat for plotting the spacing between  "
  WRITE(iunit,'(A)')"atoms.                                                                   "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"     Usage:                                                              "
  WRITE(iunit,'(A)')"             spacing -pdb pdbfile [-ca] [-het] [-lig1]                   "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"Option    Type       Value       Description                             "
  WRITE(iunit,'(A)')"------------------------------------------------------------             "
  WRITE(iunit,'(A)')" -pdb     Input                  pdb file name                           "
  WRITE(iunit,'(A)')" -ca       Opt                   Calculates Ca model only                "
  WRITE(iunit,'(A)')" -lig1     Opt                   Assigns each residue in the HETATMs one "
  WRITE(iunit,'(A)')"                                 average point, only works with -het and "
  WRITE(iunit,'(A)')"                                 -ca flags as well                       "
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

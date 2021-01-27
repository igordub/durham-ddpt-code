PROGRAM map
  USE nrutil
  USE read_files
  USE utils
  IMPLICIT NONE
  INTEGER, ALLOCATABLE,DIMENSION(:) :: atnum,resnum
  INTEGER, ALLOCATABLE,DIMENSION(:,:) :: g
  REAL(DP), ALLOCATABLE,DIMENSION(:) :: x,y,z,occ,bfac,eigenval,coor,mass
  REAL(DP), ALLOCATABLE,DIMENSION(:,:) :: eigenvec,veclow,b
  REAL(DP) :: dist,total,totold,cutoff,k,totmass,diff,scale,f,c
  CHARACTER, ALLOCATABLE, DIMENSION(:) :: atom*6,name*4,res*3,chain*1,elem*2,chag*2
  INTEGER :: natom,startvec,endvec,num,i,j,ii,jump
  CHARACTER :: pdbfile*120,eigenfile*120,dummy*120
  LOGICAL :: hetatm,qexist,getoption

  IF (getoption('-help',.false.,dummy)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF

  IF(.not.getoption('-pdb',.true.,pdbfile)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF

  INQUIRE(file=pdbfile,exist=qexist)
  IF (.not.qexist) THEN
     WRITE(6,'(A)') pdbfile, "does not exist"
     STOP "Specify with -pdb flag"
  END IF
  
  IF(.not.getoption('-i',.true.,eigenfile)) THEN
     eigenfile="matrix.eigenfacs"
  END IF

  INQUIRE(file=eigenfile,exist=qexist)
  IF (.not.qexist) THEN
     WRITE(6,'(A)') eigenfile, "does not exist"
     STOP "Specify with -i flag if not matrix.eigenfacs"
  END IF

  IF(getoption('-s',.true.,dummy)) THEN
     READ(dummy,*) startvec
  ELSE
     startvec=7
  END IF 

  IF(getoption('-e',.true.,dummy)) THEN
     READ(dummy,*) endvec
  ELSE
     endvec=31
  END IF 

  IF(getoption('-f',.true.,dummy)) THEN
     READ(dummy,*) k
  ELSE
     k=1.d0
  END IF 

  IF(getoption('-c',.true.,dummy)) THEN
     READ(dummy,*) cutoff
  ELSE
     cutoff=8.d0
  END IF 

  CALL read_pdb(pdbfile,hetatm,natom,atom,atnum,name,res,chain,resnum,x,y,z, &
       occ,bfac,elem,chag)

  ALLOCATE(coor(3*natom))
  DO i=1,natom
     coor(3*(i-1)+1)=x(i)
     coor(3*(i-1)+2)=y(i)
     coor(3*(i-1)+3)=z(i)
  END DO
  DEALLOCATE(x)
  DEALLOCATE(y)
  DEALLOCATE(z)

  CALL read_eigenfacs(eigenfile,natom,startvec,endvec,3,eigenval,eigenvec,num)

  ALLOCATE(mass(natom))
  CALL res_mass(natom,res,mass)

  totmass=0.d0
  DO i=1,natom
     totmass=totmass+mass(i)
  END DO

  total=0.d0
  DO i=1,natom
     DO j=i,natom
        dist=0.d0
        dist=dist+(coor(3*(i-1)+1)-coor(3*(j-1)+1))**2
        dist=dist+(coor(3*(i-1)+2)-coor(3*(j-1)+2))**2
        dist=dist+(coor(3*(i-1)+3)-coor(3*(j-1)+3))**2

        IF (dist.le.cutoff**2) THEN
           total=total+k*dist
        END IF
     END DO
  END DO
  
  print*, total
  
  ALLOCATE(veclow(natom,natom))
  veclow=0.d0
  DO j=1,natom
     DO ii=1,natom
        DO i=1,num
           veclow(j,ii)=veclow(j,ii)+((((coor(3*(j-1)+1)+eigenvec(3*(j-1)+1,i)-(coor(3*(ii-1)+1)+eigenvec(3*(ii-1)+1,i)))**2 + &
                (coor(3*(j-1)+2)+eigenvec(3*(j-1)+2,i)-(coor(3*(ii-1)+2)+eigenvec(3*(ii-1)+2,i)))**2 + &
                (coor(3*(j-1)+3)+eigenvec(3*(j-1)+3,i)-(coor(3*(ii-1)+3)+eigenvec(3*(ii-1)+3,i)))**2)**0.5 - &
                ((coor(3*(j-1)+1)-coor(3*(ii-1)+1))**2 + &
                (coor(3*(j-1)+2)-coor(3*(ii-1)+2))**2 + &
                (coor(3*(j-1)+3)-coor(3*(ii-1)+3))**2)**0.5)**2)/eigenval(i)
        END DO
     END DO
  END DO
  veclow=veclow**0.5

  jump=resnum(1)-1

  OPEN(file="movement.dat",form="FORMATTED",unit=5897) 
!  WRITE(5897,'(A1)') '#'
  DO i=1,natom
     DO j=1,natom
        WRITE(5897,'(1X,I5,1X,I5,1X,G11.4)') i+jump, j+jump, veclow(i,j)
     END DO
     WRITE(5897,'(1X)')
  END DO
  CLOSE(5897)




END PROGRAM map

SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
  WRITE(iunit,'(A)')"                           M  O  V  E  I  N  G                           "
  WRITE(iunit,'(A)')"                               VERSION 1.0                               "
  WRITE(iunit,'(A)')"                                                                         "  
  WRITE(iunit,'(A)')"                               Written by:                               "
  WRITE(iunit,'(A)')"                      Tom Rodgers and David Burnell                      "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program calculates the change in the distance between atoms between "
  WRITE(iunit,'(A)')"the equilibrium value and the value after applying the eigenvectors.     "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"The output file, movement.dat, can be plotted in gnuplot with the        "
  WRITE(iunit,'(A)')"command:                                                                 "
  WRITE(iunit,'(A)')'splot "moving.dat" using 1:2:3 notitle                                   '
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"     Usage:                                                              "
  WRITE(iunit,'(A)')"             movement -pdb pdbfile -i eigenfacs file [-s start]          "
  WRITE(iunit,'(A)')"                      [-e end]                                           "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"Option    Type       Value             Description                       "
  WRITE(iunit,'(A)')"------------------------------------------------------------             "
  WRITE(iunit,'(A)')" -pdb     Input                        pdb file name                     "
  WRITE(iunit,'(A)')"  -i      Input                        Eigenvector file name             "
  WRITE(iunit,'(A)')"  -s      Input,Opt  7                 First eigenvalue included         "
  WRITE(iunit,'(A)')"  -e      Input,Opt  31                Last eigenvalue included          "
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

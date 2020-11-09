PROGRAM centre
  USE nrtype
  USE utils
  USE read_files
  USE write_files
  USE trajutils
  IMPLICIT NONE
  INTEGER :: i, j, natom, nstart, io, ncon, k
  CHARACTER :: pdbfile*120, dummy*120, pdbfile2*120
  REAL(DP) :: xm,ym,zm,xmax,ymax,zmax,xmin,ymin,zmin
  INTEGER, ALLOCATABLE,DIMENSION(:) :: atnum,resnum,resvals
  REAL(DP), ALLOCATABLE,DIMENSION(:) ::x,y,z,occ,bfac,vals,mass,full
  CHARACTER, ALLOCATABLE, DIMENSION(:) :: atom*6,name*4,res*3,chain*1,elem*2,chag*2
  INTEGER, ALLOCATABLE,DIMENSION(:) :: atnum2,resnum2,resvals2
  REAL(DP), ALLOCATABLE,DIMENSION(:) ::x2,y2,z2,occ2,bfac2,vals2,full2
  CHARACTER, ALLOCATABLE, DIMENSION(:) :: atom2*6,name2*4,res2*3,chain2*1,elem2*2,chag2*2
  LOGICAL :: getoption,alig, hetatm


  IF(.not.getoption('-pdb',.true.,pdbfile)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF

  IF(getoption('-fit',.false.,dummy)) THEN
     alig=.true.
     IF(.not.getoption('-pdb2',.true.,pdbfile2)) THEN
        CALL helptext(0)
        CALL exit(0)
     END IF
  ELSE
     alig=.false.
  END IF

  IF(getoption('-het',.false.,dummy)) THEN
     hetatm=.true.
  ELSE
     hetatm=.false.
  END IF


  CALL read_pdb(pdbfile,hetatm,natom,atom,atnum,name,res,chain,resnum,x,y,z, &
       occ,bfac,elem,chag)
 
  ALLOCATE(mass(natom))
  
  CALL atom_mass(natom,elem,mass)
  
  xm=0
  ym=0
  zm=0
  DO i=1,natom
     xm=xm+x(i)
     ym=ym+y(i)
     zm=zm+z(i)
  END DO
  xm=xm/natom
  ym=ym/natom
  zm=zm/natom

  DO i=1,natom
     x(i)=x(i)-xm
     y(i)=y(i)-ym
     z(i)=z(i)-zm
  END DO  

  xmax=0
  ymax=0
  zmax=0
  xmin=0
  ymin=0
  zmin=0

  DO i=1,natom
     xmax=max(x(i),xmax)
     ymax=max(y(i),ymax)
     zmax=max(z(i),zmax)
     xmin=min(x(i),xmin)
     ymin=min(y(i),ymin)
     zmin=min(z(i),zmin)
  END DO  

  IF (alig) THEN
     xm=(xmin+xmax)/2
     ym=(ymin+ymax)/2
     zm=(zmin+zmax)/2
     
     DO i=1,natom
        x(i)=x(i)-xm
        y(i)=y(i)-ym
        z(i)=z(i)-zm
     END DO

  END IF
  xmax=0
  ymax=0
  zmax=0
  xmin=0
  ymin=0
  zmin=0
  DO i=1,natom
     xmax=max(x(i),xmax)
     ymax=max(y(i),ymax)
     zmax=max(z(i),zmax)
     xmin=min(x(i),xmin)
     ymin=min(y(i),ymin)
     zmin=min(z(i),zmin)
  END DO

  WRITE(*,'(A,3(1X,F7.3))') 'Max coordinates', xmax, ymax, zmax
  WRITE(*,'(A,3(1X,F7.3))') 'Min coordinates', xmin, ymin, zmin


  IF (alig) THEN

     CALL read_pdb(pdbfile2,hetatm,natom2,atom2,atnum2,name2,res2,chain2,resnum2,x2,y2, &
          z2,occ2,bfac2,elem2,chag2)

     ALLOCATE(full(natom*3))
     ALLOCATE(full2(natom*3))

     DO i=1,natom
        full((i-1)*3+1)=x(i)
        full((i-1)*3+2)=y(i)
        full((i-1)*3+3)=z(i)
        full2((i-1)*3+1)=x2(i)
        full2((i-1)*3+2)=y2(i)
        full2((i-1)*3+3)=z2(i)
     END DO

     CALL trajfit(full,full2,natom,1,0)

     DO i=1,natom
        x(i)=full((i-1)*3+1)
        y(i)=full((i-1)*3+2)
        z(i)=full((i-1)*3+3)
     END DO

  END IF

  fname='centred.pdb'
  CALL write_pdb(fname,atom,atnum,name,res,chain,resnum,x,y,z,&
       occ,bfac,mass,elem,chag,natom,.false.,.true.)

END PROGRAM centre

SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
 
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

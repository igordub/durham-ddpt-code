PROGRAM rms
  USE read_files
  USE nrutil
  USE utils
  IMPLICIT NONE
  CHARACTER :: lign80*80, dummy*120, filename*120, pdbfile*120
  LOGICAL :: getoption, hetatm, atmass, resmass
  INTEGER :: ndim, startvec, endvec, natom, j, i, io, numvec, vector, k, ncalpha, nexp, num
  REAL(DP) :: length, sum, coll, alpha, corr, kb, T, moyp, moyx, myscale, rmsp, rmsx, &
       rsmall
  REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: eigen, rmsval
  REAL(DP),ALLOCATABLE,DIMENSION(:) :: rmsmean, cbfac, eigenval, mass
  INTEGER, ALLOCATABLE,DIMENSION(:) :: atnum,resnum
  REAL(DP), ALLOCATABLE,DIMENSION(:) ::x,y,z,occ,bfac
  CHARACTER, ALLOCATABLE, DIMENSION(:) :: atom*6,name*4,res*3,chain*1,elem*2,chag*2

  IF (getoption('-help',.false.,dummy)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF

  IF(.not.getoption('-i',.true.,filename)) THEN
     filename='matrix.eigenfacs'
  END IF

  IF(.not.getoption('-pdb',.true.,pdbfile)) THEN
     pdbfile='CAonly.pdb'
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

  ! End vector
  IF(getoption('-e',.true.,dummy)) THEN
     READ(dummy,*) endvec
  ELSE
     endvec=31
  END IF  
  numvec=endvec-startvec+1

  IF (getoption('-het',.false.,dummy)) THEN
     hetatm=.true.
  ELSE
     hetatm=.false.
  END IF

  IF (getoption('-res',.false.,dummy)) THEN
     resmass=.true.
  ELSE
     resmass=.false.
  END IF
  
  IF (getoption('-mass',.false.,dummy)) THEN
     atmass=.true.
  ELSE
     atmass=.false.
  END IF
  
  CALL read_eigenfacs(filename,natom,startvec,endvec,ndim,eigenval,eigen,numvec)

  ALLOCATE(rmsval(natom,numvec))
  ALLOCATE(rmsmean(natom))
  ALLOCATE(cbfac(natom))

  DO i=1,natom
     DO j=1,numvec
        length=0
        DO k=1,ndim
           length=length+eigen((i-1)*ndim+k,j)**2
        END DO
        rmsval(i,j)=length**0.5
     END DO
  END DO

  CALL read_pdb(pdbfile,hetatm,natom,atom,atnum,name,res,chain,resnum,x,y,z, &
       occ,bfac,elem,chag)

  ALLOCATE(mass(natom))
! Assigning atom masses
  IF (atmass) THEN
     IF (resmass) THEN
        CALL res_mass(natom,res,mass)
        WRITE(6,'(A)') "Assigned atoms their residue mass"
     ELSE
        CALL atom_mass(natom,elem,mass)
        WRITE(6,'(A)') "Assigned atoms their true mass"
     END IF
  ELSE
     mass=1
     WRITE(6,'(A)') "Atom masses all set to 1 amu"
  END IF

  OPEN(file='rms.val',form='FORMATTED',unit=2357)
  
  WRITE(2357,'(A)',advance='NO') '#   Num At Res Resid Chain'
  DO i=1,numvec
     WRITE(2357,'(1X,A4,I3)',advance='NO') 'Mode', i+startvec-1
  END DO
  WRITE(2357,'(1X)')
  
  DO i=1,natom
     WRITE(2357,'(A6,I5,1X,A4,1X,A3,1X,A1,I4)',advance='NO') atom(i),atnum(i),name(i),res(i),chain(i),resnum(i)
     DO j=1,numvec
        WRITE(2357,'(1X,F7.4)',advance='NO') rmsval(i,j)
     END DO
     WRITE(2357,'(1X)')
  END DO

  CLOSE(2357)

  OPEN(file="collectivity.val",form="FORMATTED",unit=3858)
  WRITE(3858,'(A)') "#mode collectivity"

  DO i=1,numvec
     sum=0
     DO j=1,natom
        sum=sum+rmsval(j,i)**2
     END DO

     alpha=1/sum
     sum=0
     DO j=1,natom
        sum=sum+alpha*rmsval(j,i)**2*LOG(alpha*rmsval(j,i)**2)
     END DO
     
     coll=EXP(-sum)/natom
     WRITE(3858,'(I3,1X,F7.4)') i+startvec-1,coll
  END DO

  CLOSE(3858)
  
  kb=925.10966 ! A^2 u cm^-2 K^-1
  T=300 ! K
  rsmall=1e-04

  num=1

DO k=1,14

  nexp=0
  ncalpha=0
  DO i=num,num+371
     IF (name(i).eq.' CA ') THEN
        ncalpha=ncalpha+1
        IF (bfac(i).lt.21.d0) THEN
           ncalpha=ncalpha-1
        END IF 
        IF (bfac(i).lt.0.d0) THEN
           nexp=nexp+1
        END IF
     END IF
  END DO

  corr=0.d0
  moyp=0.d0
  moyx=0.d0
  myscale=0.d0
  rmsp=0.d0
  rmsx=0.d0
  DO i=num,num+371
     cbfac(i)=0.d0
     DO j=1,numvec
        cbfac(i)=cbfac(i)+rmsval(i,j)**2/(108.591**2*eigenval(j))  ! e^2/w^2
     END DO
     cbfac(i)=2*TWOPI**2*kb*T*cbfac(i)/(3*mass(i))
     
     IF (nexp.eq.0.and.name(i).eq.' CA '.and.bfac(i).ge.21.d0 ) THEN
        corr=corr+cbfac(i)*bfac(i)
        moyp=moyp+cbfac(i)
        moyx=moyx+bfac(i)
        myscale=myscale+bfac(i)/cbfac(i)
        rmsp=rmsp+cbfac(i)*cbfac(i)
        rmsx=rmsx+bfac(i)*bfac(i) 
     END IF
  END DO

  OPEN(file='mode.bfactors',form='FORMATTED',unit=6745)

 
  IF (ncalpha.gt.0) THEN
     moyp=moyp/dfloat(ncalpha)
     moyx=moyx/dfloat(ncalpha)
     myscale=myscale/dfloat(ncalpha)
     rmsp=rmsp/dfloat(ncalpha)-moyp*moyp
     rmsx=rmsx/dfloat(ncalpha)-moyx*moyx
     rmsp=sqrt(rmsp) 
     rmsx=sqrt(rmsx) 
 
     IF (nexp.eq.0) THEN
        IF (rmsx.gt.rsmall) THEN
           corr=corr/dfloat(ncalpha)-moyp*moyx
           IF (rmsp.gt.rsmall) corr=corr/(rmsp*rmsx)
           WRITE(6,'(/A,F10.3,A,I6,A)') ' Correlation= ',corr,' for ',ncalpha,' C-alpha atoms.'
           WRITE(6745,'(A,F10.3,A,I6,A)')'#Correlation= ',corr,' for ',ncalpha,' C-alpha atoms.'
        ELSE
           WRITE(6,'(/A/)')' Experimental B-factors are nearly constant !'
           WRITE(6745,'(A)')'#Experimental B-factors are nearly constant !'
        END IF
     END IF
 
     WRITE(6,'(A,F10.3,A,F7.2)')' <B-predict>= ',moyp,' +/- ',rmsp
     WRITE(6,'(A,F10.3,A,F7.2)')' <B-crystal>= ',moyx,' +/- ',rmsx
     WRITE(6,'(A,F10.3)')       ' Shiftng-fct= ',moyx-moyp
     IF (rmsp.gt.rsmall) THEN
        WRITE(6,'(A,F10.3)')    ' Scaling-fct= ',rmsx/rmsp
        WRITE(6,'(A,F10.3)')    ' Just Scale = ',myscale 
        WRITE(6,'(/A)')         ' Predicted, Scaled and Experimental B-factors are saved.'
     END IF
 
     WRITE(6745,'(A)')'#     Number Name Res   Num   Predicted   FullScaled    Scaled    Experimental'

     DO i=num,num+371
        WRITE(6745,'(A6,I5,1X,A4,1X,A3,1X,A1,I4,1X,4F12.4)') atom(i),atnum(i),name(i),res(i),chain(i),resnum(i), &
             cbfac(i),(cbfac(i)-moyp)*rmsx/rmsp+moyx,cbfac(i)*myscale,bfac(i)
     END DO

  END IF

  num=num+372

END DO


  CLOSE(6745)

  DEALLOCATE(eigen)
  DEALLOCATE(rmsval)
  DEALLOCATE(rmsmean)

END PROGRAM rms

SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
  WRITE(iunit,'(A)')"                           R  M  S  /  C  O  L                           "
  WRITE(iunit,'(A)')"                               VERSION 1.0                               "
  WRITE(iunit,'(A)')"                                                                         "  
  WRITE(iunit,'(A)')"                               Written by:                               "
  WRITE(iunit,'(A)')"                      Tom Rodgers and David Burnell                      "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program calculates the RMSD of each of the atoms for each mode,     "
  WRITE(iunit,'(A)')"output in rms.val, from an eigenfacs formated file which can be produced "
  WRITE(iunit,'(A)')"from the included ENM code, GNM code, GroNM, GroED, or AmberED packages. "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program also outputs the collectivity of each mode, output in       "
  WRITE(iunit,'(A)')"collectivity.val.                                                        "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program also calculates the predicted bfactors and compares them    "
  WRITE(iunit,'(A)')"with the experimental values given in the original pdbfile, output in    "
  WRITE(iunit,'(A)')"mode.bfactors.                                                           "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"     Usage:                                                              "
  WRITE(iunit,'(A)')"             rms [-i infile] [-pdb pdbfile] [-s startvec]                "
  WRITE(iunit,'(A)')"                 [-e endvec] [-d1] [-het]                                "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"Option    Type       Value             Description                       "
  WRITE(iunit,'(A)')"------------------------------------------------------------             "
  WRITE(iunit,'(A)')"  -i      Input,Opt  matrix.eigenfacs  First eigenvector file            "
  WRITE(iunit,'(A)')" -pdb     Input,Opt  CAonly.pdb        pdb file                          "
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

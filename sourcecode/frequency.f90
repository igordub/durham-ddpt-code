!  program to calculate the frequency of the normal modes in cm-1
!  Tom Rodgers 31 04 2011
PROGRAM frequency
  USE read_files
  USE nrutil
  IMPLICIT NONE
  REAL(DP):: Z,en,S,Gsch,Ssch,T
  REAL(DP),ALLOCATABLE,DIMENSION(:) :: eigen, freq
  REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: eigenvec
  INTEGER:: num, i,io, startvec, endvec
  CHARACTER:: eigenfile*120,lign80*80,dummy*120
  LOGICAL:: getoption, qexist

  IF (getoption('-help',.false.,dummy)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF
  
  IF(.not.getoption('-i',.true.,eigenfile)) THEN
     eigenfile="matrix.eigenfacs"
  END IF

  INQUIRE(file=eigenfile,exist=qexist)
  IF (.not.qexist) THEN
     WRITE(6,'(A)') eigenfile, "does not exist"
     STOP "Specify with -i flag if not matrix.eigenfacs"
  END IF
  
  IF(getoption('-t',.true.,dummy)) THEN
     READ(dummy,*) T
  ELSE
     T=298.0d0
  END IF

  IF(getoption('-s',.true.,dummy)) THEN
     READ(dummy,*) startvec
  ELSE
     startvec=1
  END IF 

  IF(getoption('-e',.true.,dummy)) THEN
     READ(dummy,*) endvec
  ELSE
     endvec=1000
  END IF 

  CALL read_eigenfacs(eigenfile,i,startvec,endvec,1,eigen,eigenvec,num)

  ALLOCATE(freq(num))
  
  DO i=1,num
     freq(i) = 108.591*dsqrt(eigen(i))
  END DO

  OPEN(file="mode.frequencies",unit=9874)
  WRITE(9874,'(A)')"# Normal mode frequencies - can be read into grace and plotted as a histogram"
  WRITE(9874,'(A)')"# in cm-1 assuming that Hessian was in CHARMM units"

  OPEN(file="mode.energy",unit=9875)
  WRITE(9875,'(A)')"#  Mode      Z         G/kT      Gs/kT      S/k      Ss/k   "

  DO i=1,num
     WRITE(9874,'(F12.6)') freq(i)
     Z=1/(1-dexp(-1.438775057*freq(i)/T))
     en=-dlog(Z)
     S=(1.438775057*freq(i)/T)/(dexp(1.438775057*freq(i)/T)-1)-dlog(1-dexp(-1.438775057*freq(i)/T))
     Gsch=-0.5d0*(dlog(1+(dexp(1.0d0))**2/(1.438775057*freq(i)/T)**2)+2*(1.438775057*freq(i)/T)/ &
          dexp(1.0d0)*datan(dexp(1.0d0)/(1.438775057*freq(i)/T))-2)
     Ssch=0.5d0*dlog((1+(dexp(1.0d0))**2/(1.438775057*freq(i)/T)**2))
     WRITE(9875,'(1X,I5,1X,F18.16,4(F18.16))') i, Z, en, Gsch, S, Ssch
  END DO
 
  CLOSE(9874)
  CLOSE(9875)

END PROGRAM frequency


SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
  WRITE(iunit,'(A)')"                           F  R  E  Q  /  E  N                           "
  WRITE(iunit,'(A)')"                               VERSION 1.0                               "
  WRITE(iunit,'(A)')"                                                                         "  
  WRITE(iunit,'(A)')"                               Written by:                               "
  WRITE(iunit,'(A)')"                      Tom Rodgers and David Burnell                      "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program calculates the normal mode frequencies, output in           "
  WRITE(iunit,'(A)')"mode.frequencies, from an eigenfacs formated file which can be produced  "
  WRITE(iunit,'(A)')"from the included ENM code, GNM code, and the GroAMED packages.          "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program also outputs the energy parameters for the normal modes,    "
  WRITE(iunit,'(A)')"output in mode.energy. The value of the partition function, Z; the       "
  WRITE(iunit,'(A)')"full harmonic free energy, G; the Schlitter free energy, Gs; the full    "
  WRITE(iunit,'(A)')"harmanic entropy, S; and the Schlitter entropy, Ss, are given for each   "
  WRITE(iunit,'(A)')"mode.                                                                    "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"     Usage:                                                              "
  WRITE(iunit,'(A)')"             frequency [-i infile] [-t temperature]                      "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"Option    Type       Value             Description                       "
  WRITE(iunit,'(A)')"------------------------------------------------------------             "
  WRITE(iunit,'(A)')"  -i      Input,Opt  matrix.eigenfacs  Eigenvector file name             "
  WRITE(iunit,'(A)')"  -t      Input,Opt  298               Temperature of entropy calculation"
  WRITE(iunit,'(A)')"  -s      Input,Opt  1                 First eigenvalue included         "
  WRITE(iunit,'(A)')"  -e      Input,Opt  1000              Last eigenvalue included          "
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


PROGRAM plot_on_pdb
  USE read_files
  USE nrutil
  USE write_files
  IMPLICIT NONE
  CHARACTER :: lign80*80, pdbfile*120,listfile*120,dummy*120,fname*120
  INTEGER :: skipvec, numvec, endvec, startvec, ndim, vector, i, j, k, io, natom, m
  INTEGER :: nat, resid, residold, st, en, numvals, adjust,offset
  INTEGER, ALLOCATABLE,DIMENSION(:) :: atnum,resnum,resvals
  REAL(DP), ALLOCATABLE,DIMENSION(:) ::x,y,z,occ,bfac,vals,mass
  CHARACTER, ALLOCATABLE, DIMENSION(:) :: atom*6,name*4,res*3,chain*1,elem*2,chag*2
  REAL(DP) :: sum, dot, lengthi, lengthj,dx,dy,dz,cutoff,scale
  LOGICAL :: hetatm,getoption


  IF (getoption('-help',.false.,dummy)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF

  IF (.not.getoption('-i',.true.,listfile)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF

  IF (.not.getoption('-pdb',.true.,pdbfile)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF

  IF (getoption('-het',.false.,dummy)) THEN
     hetatm=.true.
  ELSE
     hetatm=.false.
  END IF

  IF (getoption('-o',.true.,dummy)) THEN
     READ(dummy,*) offset
  ELSE
     offset=0
  END IF

    IF (getoption('-scale',.true.,dummy)) THEN
     READ(dummy,*) scale
  ELSE
     scale=1.d0
  END IF


  CALL read_pdb(pdbfile,hetatm,natom,atom,atnum,name,res,chain,resnum,x,y,z, &
                  occ,bfac,elem,chag)


  OPEN(file=listfile,form="FORMATTED",status="OLD",unit=2357)

  i=0
  DO
     READ(2357,'(A)',IOSTAT=io) lign80
     IF (io > 0) THEN
        WRITE(*,*) 'Check input.  Something was wrong'
        STOP
        ! End of file
     ELSE IF (io < 0) THEN
        EXIT
        ! Count number of atoms
     ELSE
       i=i+1
     END IF
  END DO

  numvals=i

  ALLOCATE(resvals(numvals))
  ALLOCATE(vals(numvals))
  
  REWIND(2357)

  DO i=1,numvals
     READ(2357,*) resvals(i), vals(i)
  END DO

  CLOSE(2357)

!  DO i=1,natom
!     IF (chain(i).ne.chain(1)) THEN
!        EXIT
!     END IF
!  END DO

!  adjust=resnum(i)-resnum(1)+offset

!  DO j=i,natom
!     IF (chain(j).ne.chain(i)) THEN
!        WRITE(*,*) "More than 2 different chains"
!        EXIT
!     END IF
!     resnum(j)=resnum(j)-adjust
!  END DO

  DO i=1,natom
     m=0
     DO j=1,numvals
        IF (resnum(i).eq.resvals(j)) THEN
           bfac(i)=vals(j)
           m=1
        END IF
     END DO
     IF (m.eq.0) THEN
        bfac(i)=0
        WRITE(*,*) "Res number ",resnum(i)," is not listed."
     END IF
  END DO

  ALLOCATE(mass(natom))
  mass=0
  
  fname='plot.pdb'

  CALL write_pdb(fname,atom,atnum,name,res,chain,resnum,x,y,z,&
       occ,scale*bfac,mass,elem,chag,natom,.false.,hetatm)

  DEALLOCATE(x)
  DEALLOCATE(y)
  DEALLOCATE(z)
  DEALLOCATE(resvals)
  DEALLOCATE(vals)
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


END PROGRAM plot_on_pdb


SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
  WRITE(iunit,'(A)')"                           P  L  O  T  P  D  B                           "
  WRITE(iunit,'(A)')"                               VERSION 1.0                               "
  WRITE(iunit,'(A)')"                                                                         "  
  WRITE(iunit,'(A)')"                               Written by:                               "
  WRITE(iunit,'(A)')"                      Tom Rodgers and David Burnell                      "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program copies values for each residue into the bfactor position    "
  WRITE(iunit,'(A)')"in a pdb, if Ca only pdb then the values will be assigned to each Ca     "
  WRITE(iunit,'(A)')"atom, if a full pdf then the values will be assigned to each residue. It "
  WRITE(iunit,'(A)')"adjusts the second chain assuming that the same residue starts, the      "
  WRITE(iunit,'(A)')"offset can be adjusted with -o num.                                      "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"overlap.dat contains the overlap of each atom and sum.dat contains the   "
  WRITE(iunit,'(A)')"average values for the whole system.                                     "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"     Usage:                                                              "
  WRITE(iunit,'(A)')"             plotpdb -i infile -pdb pdbfile [-het] [-o offset]           "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"Option    Type      Value     Description                                "
  WRITE(iunit,'(A)')"------------------------------------------------------------             "
  WRITE(iunit,'(A)')"  -i      Input               Res data file  |resnumber value| no header "
  WRITE(iunit,'(A)')"  -pdb    Input               pdbfile                                    "
  WRITE(iunit,'(A)')"  -het    Opt                 Incudes reading HETATM records with flag   "
  WRITE(iunit,'(A)')"  -o      Input,Opt 0         Value of the offset if needed              "
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

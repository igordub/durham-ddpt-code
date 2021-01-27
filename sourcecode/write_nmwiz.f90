! Program to produce NMWiz files to use VMD plugin
! Tom Rodgers 31 04 2011
PROGRAM write_nmwiz
  USE read_files
  USE nrutil
  IMPLICIT NONE
  INTEGER :: io, natom, i, serial(5), char(2), mode, num, modes, resid
  INTEGER, ALLOCATABLE, DIMENSION(:) :: resSeq, atnum
  REAL(DP) :: temp, xm, ym, zm, eigen
  REAL(DP),ALLOCATABLE, DIMENSION(:) :: x, y, z, bfactor,occ
  CHARACTER(6),ALLOCATABLE, DIMENSION(:) :: atomname
  CHARACTER(4),ALLOCATABLE, DIMENSION(:) :: name
  CHARACTER(3),ALLOCATABLE, DIMENSION(:) :: resName
  CHARACTER(1),ALLOCATABLE, DIMENSION(:) :: chainID
  CHARACTER(2),ALLOCATABLE, DIMENSION(:) :: elem,chag
  CHARACTER :: altLoc, filename*120, linerd*80, iCode, vect*7, val*5, heta*5, &
               dummy*120, eigenfile*120, bfacfile*120, lign80*80
  LOGICAL :: hetatm, getoption
  
  IF (getoption('-help',.false.,dummy)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF
  
  IF(.not.getoption('-i',.true.,eigenfile)) THEN
     eigenfile='matrix.eigenfacs'
  END IF
  
  IF(.not.getoption('-pdb',.true.,filename)) THEN
     filename='CAonly.pdb'
  END IF

  IF (getoption('-het',.false.,dummy)) THEN
     hetatm=.true.
  ELSE
     hetatm=.false.
  END IF

  IF(getoption('-n',.true.,dummy)) THEN
     READ(dummy,*) modes
  ELSE
     modes=25
  END IF

  IF(.not.getoption('-b',.true.,bfacfile)) THEN
     bfacfile='mode.bfactors'
  END IF


  ! Open the eigenvector file
  OPEN(file=eigenfile,form="FORMATTED",status="OLD",unit=4539)
  
  ! Open the nmwiz file for writting
  OPEN(file='nmwiz.nmd',form="FORMATTED",unit=4540)
  
  ! Open the bfactors file
  OPEN(file=bfacfile,form="FORMATTED",status="OLD",unit=4541)
  
  CALL read_pdb(filename,hetatm,natom,atomname,atnum,name,resName,chainID,resSeq,x,y,z, &
       occ,bfactor,elem,chag)
   
  ! Write out autoload header
  WRITE(4540, '(A20)') 'nmwiz_load nmwiz.nmd'
  
  ! Write out the atomnames
  WRITE(4540,'(A9,1X)',advance="NO") 'atomnames'
  DO i=1,natom
     WRITE(4540,'(A4,1X)',advance="NO") name(i)
  END DO
  WRITE(4540,'(1X)')
  
  ! Write out the resnames
  WRITE(4540,'(A8,1X)',advance="NO") "resnames"
  DO i=1,natom
     WRITE(4540,'(A3,1X)',advance="NO") resName(i)
  END DO
  WRITE(4540,'(1X)')
  
  ! Write out the chainids
  WRITE(4540,'(A8,1X)',advance="NO") "chainids"
  DO i=1,natom
     WRITE(4540,'(A1,1X)',advance="NO") chainID(i)
  END DO
  WRITE(4540,'(1X)')
  
  ! Write out resids
  WRITE(4540,'(A6,1X)',advance="NO") "resids"
  DO i=1,natom
     WRITE(4540,'(I5,1X)',advance="NO") resSeq(i)
  END DO
  WRITE(4540,'(1X)')

  i=0
  DO
     READ(4541,'(A)',IOSTAT=io) lign80
     IF (io > 0) THEN
        WRITE(*,*) 'Check input.  Something was wrong'
        STOP
        ! End of file
     ELSE IF (io < 0) THEN
        EXIT
        ! Count up number of atoms
     ELSE
        IF ((index(lign80,'#').gt.10).or.(index(lign80,'#').eq.0)) THEN
           i=i+1
           READ(lign80,'(43X,F12.4)') bfactor(i)
        END IF
     END IF
  END DO
  
  ! Write out the bfactors
  WRITE(4540,'(A8,1X)',advance="NO") "bfactors"
  DO i=1,natom
     WRITE(4540,'(F9.4,1X)',advance="NO") bfactor(i)
  END DO
  WRITE(4540,'(1X)')
  
  ! Write out the coordinates
  WRITE(4540,'(A11,1X)',advance="NO") "coordinates"
  DO i=1,natom
     WRITE(4540,'(F8.3,1X)',advance="NO") x(i)
     WRITE(4540,'(F8.3,1X)',advance="NO") y(i)
     WRITE(4540,'(F8.3,1X)',advance="NO") z(i)
  END DO
  WRITE(4540,'(1X)')
  
  ! Write out the modes
  DO mode=1,modes
     
     ! Read eigenvalue
     READ(4539,'(24X,G12.4)') eigen
     ! Skip over ------- line
     READ(4539,*)
     ! Scale factor is:
     eigen=dsqrt(1/eigen)
     
     ! Write the mode heading
     WRITE(4540,'(A4,1X,I6,1X,F12.3,1X)',advance="NO") "mode", mode, eigen
     
     DO i=1,natom
        ! Read in 
        READ(4539,*) xm, ym, zm
        ! Print the modes to the output
        WRITE(4540,'(F9.6,1X)',advance="NO") xm
        WRITE(4540,'(F9.6,1X)',advance="NO") ym
        WRITE(4540,'(F9.6,1X)',advance="NO") zm
     END DO
     WRITE(4540,'(1X)')

  END DO
  
  CLOSE(4539)
  CLOSE(4540)
  CLOSE(4541)
  
  DEALLOCATE(x)
  DEALLOCATE(y)
  DEALLOCATE(z)
  DEALLOCATE(atomname)
  DEALLOCATE(name)
  DEALLOCATE(resName)
  DEALLOCATE(chainID)
  DEALLOCATE(resSeq)
  DEALLOCATE(bfactor)
  
  WRITE(*,*) 'nmwiz> writing nmwiz file'
  WRITE(*,*) 'nmwiz> run with vmd -e nmwiz.nmd'
  
END PROGRAM write_nmwiz

SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
  WRITE(iunit,'(A)')"                           N  M  W  I  Z  W  T                           "
  WRITE(iunit,'(A)')"                               VERSION 1.0                               "
  WRITE(iunit,'(A)')"                                                                         "  
  WRITE(iunit,'(A)')"                               Written by:                               "
  WRITE(iunit,'(A)')"                      Tom Rodgers and David Burnell                      "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program creates a file suitable for use with the VMD plugin NMWiz,  "
  WRITE(iunit,'(A)')"output as nmwix.nmd, from an eigenfacs formated file which can be        "
  WRITE(iunit,'(A)')"produced from the included ENM code, GNM code, GroNM, GroED, or AmberED  "
  WRITE(iunit,'(A)')"packages. It also requires the matching pdb and the mode.bfactors file   "
  WRITE(iunit,'(A)')"produced with the rms program.                                           "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"     Usage:                                                              "
  WRITE(iunit,'(A)')"             w_nmwiz [-i infile] [-pdb pdbfile] [-n numberofmodes]       "
  WRITE(iunit,'(A)')"                     [-b bfactorfile] [-het]                             "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"Option    Type       Value             Description                       "
  WRITE(iunit,'(A)')"------------------------------------------------------------             "
  WRITE(iunit,'(A)')"  -i      Input,Opt  matrix.eigenfacs  Eigenvector file                  "
  WRITE(iunit,'(A)')" -pdb     Input,Opt  CAonly.pdb        PDB filename                      "
  WRITE(iunit,'(A)')"  -n      Input,Opt   25               Number of modes written           "
  WRITE(iunit,'(A)')"  -b      Input,Opt  mode.bfactors     Calculated bfactor file           "
  WRITE(iunit,'(A)')"  -het     Opt                Incudes reading HETATM records with flag   "
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

PROGRAM diagrtb
  USE nrutil
  USE read_files
  USE utils
  USE rtb_util
  IMPLICIT NONE
  INTEGER :: nbrt    
  CHARACTER :: cformat*12, cstatus*12, feig*120, &
       field*128, fmtx*128, fpdb*128, &
       namfich*128, motinp*128, &
       string*128, sstr*4, dummy*120
  LOGICAL ::   alive, qexist, qinter, qstop, hetatm, getoption, atmass, eigen, mat, resmass, covar, fidet, is_numeric
  INTEGER ::   i,  k, klist, ksep, &
       lmot, lnompdb, natblocs, natom, nb, nddres, nddres2, ndim, &
       nmots, nrbl, nunit, nvec, undat, uninp, ii, unrtb
  REAL(DP) :: rave, rdev, rmax, rmin
  INTEGER,ALLOCATABLE,DIMENSION(:) :: blength
  CHARACTER,ALLOCATABLE,DIMENSION(:) :: ssu*1

  INTEGER, ALLOCATABLE,DIMENSION(:) :: atnum,resnum
  REAL(DP), ALLOCATABLE,DIMENSION(:) ::x,y,z,occ,bfac,corlin,mass,masshold
  CHARACTER, ALLOCATABLE, DIMENSION(:) :: atom*6,name*4,res*3,chain*1,elem*2,chag*2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nbrt=6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  WRITE(6,'(/A)') 'Starting D  I  A  G  R  T  B'

! Read in options from tags
  
  IF (getoption('-help',.false.,dummy)) THEN
     CALL helptext(0)
     CALL EXIT(0)
  END IF
  
! Number of residues per block
  IF (getoption('-r',.true.,dummy)) THEN
     IF (is_numeric(dummy(1:1))) THEN
        READ(dummy,*) nrbl
     ELSE
        STOP "The value of r following -r must be a positive real number"
     END IF
  ELSE
     nrbl=1
  END IF

! Substructuring options
! NONE: Only the number of residues per block is taken into account
!   (Chain breaks interrupt blocks in RESIdue case)
! SECO: Whole secondary structures are put in blocks  
!   Structures are defined in the coordinate file, by the
!   HELIX and  SHEET keyworks.
!   Non-specified portions of the sequence are split
!   into blocks, according to the NRBL value.
! DOMN: Whole specified domains are put into blocks
  
  IF (getoption('-str',.true.,dummy)) THEN
     READ(dummy,'(A4)') sstr
  ELSE
     sstr = 'NONE'
  END IF

! Atomic masses:
  IF (getoption('-mass',.false.,dummy)) THEN
     atmass = .true.
  ELSE
     atmass = .false.
  END IF

! Residue masses:
  IF (getoption('-res',.false.,dummy)) THEN
     resmass = .true.
  ELSE
     resmass = .false.
  END IF

! Number of eigenvectors to compute:
  IF (getoption('-e',.true.,dummy)) THEN
     IF (is_numeric(dummy(1:1))) THEN
        READ(dummy,*) nvec
     ELSE
        STOP "The number of vectors following -e must be a positive real number"
     END IF 
  ELSE
     nvec=56
  END IF

  IF (getoption('-covar',.false.,dummy)) THEN
     covar=.true.
  ELSE
     covar=.false.
  END IF
  
! PDB default filename with the coordinates of the system.
  IF (.not.getoption('-pdb',.true.,fpdb)) THEN
     CALL helptext(0)
     CALL EXIT(0)         
  END IF

  INQUIRE(file=fpdb,exist=qexist)  
  IF (.not.qexist) THEN
     STOP 'pdb file specified does not exist'
  END IF
  
! Matrix default filename with the non-zero elements.
  IF (.not.getoption('-i',.true.,fmtx)) THEN
     fmtx = 'matrix.sdijf'     
  END IF

  INQUIRE(file=fmtx,exist=qexist)  
  IF (.not.qexist) THEN
     WRITE(6,'(A)') fmtx, "does not exist"
     STOP 'specify with -i if not matrix.sdijf'
  END IF

! Eigenvector output filename.
  feig = 'matrix.eigenfacs'
  
! READ HETATMs
  IF (getoption('-het',.false.,dummy)) THEN
     hetatm=.true.
  ELSE
     hetatm=.false.
  END IF

! Keep eigenvalues for RTBs
  IF (getoption('-eig',.false.,dummy)) THEN
     eigen=.true.
  ELSE
     eigen=.false.
  END IF

! Print RTB matrix
  IF (getoption('-mat',.false.,dummy)) THEN
     mat=.true.
  ELSE
     mat=.false.
  END IF

  IF (getoption('-det',.false.,dummy)) THEN
     fidet=.true.
  ELSE
     fidet=.false.
  END IF

! Read in values from pdb file
  CALL read_pdb(fpdb,hetatm,natom,atom,atnum,name,res,chain,resnum,x,y,z, &
                  occ,bfac,elem,chag)
  DEALLOCATE(atom)
  DEALLOCATE(atnum)
  DEALLOCATE(name)
  DEALLOCATE(occ)
  DEALLOCATE(bfac)
  DEALLOCATE(chag) 

! Find blocks from pdb structure
  CALL BLOCPDB(fpdb, sstr, nrbl, resnum, blength, chain, natblocs, natom, nb, hetatm)
  DEALLOCATE(resnum)
  DEALLOCATE(chain)

  ! Allocate new holders for the position and the mass
  ALLOCATE(corlin(3*natom))
  ALLOCATE(masshold(natom))
  ALLOCATE(mass(3*natom))
  ! Assign mass values
  IF (atmass) THEN
     IF (resmass) THEN
        CALL res_mass(natom,res,masshold)
        WRITE(6,'(A)') "Adjusted Ca mass to residue mass"   
     ELSE
        CALL atom_mass(natom,elem,masshold)
        WRITE(6,'(/A)')'Masses were assigned based on element types.'
     END IF
     WRITE(6,'(/A)')'Mass statistics: '
     CALL vecstat(masshold,natom,rmin,rmax,rave,rdev)
     WRITE(6,'(4(A,F12.6))')' <m>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax
  ELSE
     masshold=1.0d0
     WRITE(6,'(A)')'All masses set to unity.'
  END IF

! Assign position and mass holders
  DO i=1,natom
     corlin(3*(i-1)+1)=x(i)
     corlin(3*(i-1)+2)=y(i)
     corlin(3*(i-1)+3)=z(i)
     mass(3*(i-1)+1)=masshold(i)
     mass(3*(i-1)+2)=masshold(i)
     mass(3*(i-1)+3)=masshold(i)
  END DO

!  DEALLOCATE(masshold)
!  DEALLOCATE(x)
!  DEALLOCATE(y)
!  DEALLOCATE(z)
!  DEALLOCATE(elem)
!  DEALLOCATE(res)

  ndim = 6*nb
!   
 IF (ndim.lt.nvec) THEN
     nvec = ndim
     WRITE(6,'(/I10,A)') nvec,' eigenvectors, only, can be determined.'
  END IF

! Read in the matrix  
  CALL PREPMAT(fmtx, natblocs, natom, nb, blength)

  DEALLOCATE(blength)

! Translate into rotational translational blocks  
  CALL RTB(natom, nb, nbrt, 3*natblocs, (3*natblocs)**2,corlin,mass)

  DEALLOCATE(corlin)
  DEALLOCATE(mass)

! Find eigenvalues and vectors  
  CALL DIAGSTD (ndim, nvec, covar, mat, fidet)

! Translate back to atomistic  
  CALL RTBTOMODES (natom, nb, nbrt, 3*natblocs, nvec, feig)

! Remove temporary files
  OPEN(10,file='diagrtb_work.blocs')
  CLOSE(10,status='delete')
  OPEN(10,file='diagrtb_work.matblocs')
  CLOSE(10,status='delete')
  OPEN(10,file='diagrtb_work.sdijb')
  CLOSE(10,status='delete')
  IF (.not.eigen) THEN
     OPEN(10,file='diagrtb_work.eigenfacs')
     CLOSE(10,status='delete')
  END IF
  
  WRITE(6,'(/A)') 'End of D  I  A  G  R  T  B'
  
END PROGRAM diagrtb

!----------------------------------------------------------------------

SUBROUTINE DIAGSTD(ndim,nvecout,covar,mat,fidet)
  USE nrutil
  USE utils
  USE write_files
  IMPLICIT NONE

  LOGICAL :: mat, covar, fidet
  INTEGER :: i, ii, j, jj, k, natom, nbig, nbelow, ndim,  &
       nord, nredond, ntrace, nvec, nvecout, nunit, nzero,  &
       rdunit, unmess, unmodes, io
  REAL(DP) :: matrd, trace, det, FindDet
  REAL(DP),ALLOCATABLE,DIMENSION(:) :: work, ev
  REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: amat, b
  CHARACTER :: cformat*20, cstatus*20, nomfich*40, fname*120
     
  unmess=6
  nvec=ndim
  
  nunit=10
  rdunit=nunit
  nunit=nunit+1
  
  OPEN(file='diagrtb_work.sdijb',form='UNFORMATTED',status='OLD',unit=rdunit)
  
  k=0
  nord=0
  DO
     READ(rdunit,IOSTAT=io) i,j
     IF (io < 0) THEN
        EXIT
     END IF
     
     k=k+1
     IF (i.le.0.or.j.le.0) THEN
        WRITE(unmess,'(/A,I9,2(A,I6))') 'Line: ',k,' I= ',i,' J= ',j
        STOP '*Wrong matrix*'
     END IF
     IF (i.gt.nord) nord=i
     IF (j.gt.nord) nord=j
  END DO
  
  WRITE(unmess,'(A,I9)')'Projected matrix order =',nord
  WRITE(unmess,'(A,I9)')'Nb of non-zero elements:',k
  
  IF (nord.gt.ndim) THEN
     WRITE(unmess,'(/A)') 'Matrix can not be read.'
     IF (nord.gt.ndim) WRITE(unmess,'(2(A,I9))') &
          ' Nord=  ',nord,' > Ndim=  ',ndim
     STOP '*Wrong matrix*'
  END IF
  
  rewind(rdunit)
  
  nredond=0
  ntrace=0
  trace=0.d0
  nbig=0

  ALLOCATE(amat(nord,nord))
  amat=0.d0
  
  DO jj=1,k
     
     READ(rdunit,IOSTAT=io) i,j,matrd
     IF (io > 0) THEN
        WRITE(unmess,'(/A,I6)') 'Reading ligne ',k
        WRITE(unmess,'(2I6,F16.8)') ' i, j, matrd= ',i,j,matrd
        STOP '*Wrong matrix*'
     END IF
     
     IF (dabs(matrd).gt.0.d0) THEN
        amat(i,j)=matrd
        amat(j,i)=matrd
        IF (i.eq.j) THEN 
           trace=trace+matrd
           ntrace=ntrace+1
        END IF
        IF (matrd.gt.1E+10) THEN
           nbig=nbig+1
           IF (nbig.lt.10) THEN
              WRITE(unmess,'(A,2I12,A,G12.3)')'Element: ',i,j,' = ',matrd
           ELSE 
              IF (nbig.eq.10) WRITE(unmess,*) '...'
           END IF
        END IF
     ELSE
        nredond=nredond+1
     END IF
  END DO
  
  IF (nredond.gt.0) &
       WRITE(unmess,'(A,I9)') 'Nb of matrix elements found twice: ',nredond
  IF (nbig.gt.0) &
       WRITE(unmess,'(A,I9)') 'Nb of elements    > 1E+10 :',nbig
  
  WRITE(6,'(A,F12.4)') 'Projected matrix trace = ',trace
  
  IF (nord-ntrace.ne.0) &
       WRITE(unmess,'(A,I11)') 'Nb on zero elements on the trace:',nord-ntrace
  
!     Diagonalisation:
  
  WRITE(unmess,'(A)')'Diagonalization.'
  
  IF (nvec.gt.nord) nvec=nord
  WRITE(unmess,'(/I7,A)')nvec,' eigenvectors are computed.'
  IF (nvecout.gt.nvec) nvecout=nvec
  
  WRITE(unmess,'(I7,A)') nvecout,' of them to be saved.'
  
  IF (mat) THEN
     OPEN(file='matrix.rtb',form='FORMATTED',unit=7354)
     DO i=1,nord
        DO j=i,nord
           IF (amat(i,j).ne.0) THEN
              WRITE(7354,'(2I10,1PG20.12)') i, j, amat(i,j)
           END IF
        END DO
     END DO 
  END IF


  IF (fidet) THEN
     WRITE(6,'(/A)')'Calculating Covariance matrix determinate'
     ALLOCATE(b(nord,nord))
     b=amat
     det=FindDet(b,nord)
     DEALLOCATE(b)
     WRITE(6,*) "Covarience matrix determinate = ", det
  END IF


  ALLOCATE(work(nord))
  ALLOCATE(ev(nord))

  CALL TRED2(amat,ev,work)
  CALL TQLI(ev,work,amat)

  DEALLOCATE(work)

  IF (covar) THEN
     CALL sort_max(amat,ev,nord)
     ev=0.5961752/ev
  ELSE
     CALL sort_min(amat,ev,nord)
  END IF


! Write eigenfacs file for RTBs
  fname='diagrtb_work.eigenfacs'
  CALL write_eigenfacs(fname,amat(:,1:nvecout),ev(1:nvecout),nvecout,nord)
  
  DEALLOCATE(ev)
  DEALLOCATE(amat)

END SUBROUTINE DIAGSTD

SUBROUTINE PREPMAT (fmtx, natblocs, natom, nb, tabb)
  USE nrutil
  IMPLICIT NONE
  !     Purpose:
  !     =======
  
  !     --------------------------------------------
  !     Matrix preparation: it is split into blocks.
  !     --------------------------------------------
  !
  !     Arguments
  !     =========
  !
  !     natblocs  : maximum number of atoms found in a block.
  !     natom     : number of atoms in the system.
  !     nb        : number of blocks the system is split into.
  !     qtrace    : flag to remember that a given diagonal element is known.
  !     indexi    : i-index of non-zero elements in block hw.
  !     indexj    : j-index of non-zero elements in block hw.
  !     tabb      : block lengths.
  !     a         : band-matrix corresponding to a given block.
  !     hw        : block being filled.
  ! 
  !-----------------------------------------------------------------------  
  CHARACTER :: fmtx*120
  INTEGER :: natblocs, natom, nb 
  LOGICAL,ALLOCATABLE,DIMENSION(:) :: qtrace
  INTEGER :: tabb(nb)
  REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: a
  REAL(DP),ALLOCATABLE,DIMENSION(:) :: hw
  LOGICAL :: qinit
  INTEGER :: bloc, i, ii, iii, imax, j, jbloc, k, kk, &
       nblocs, nbnnul, nempty, ni1, ni2, nj1, nj2, & 
       nlign, nn, nread, ntotn1, tbloc, tt, io
  REAL(DP) :: toto, trace
  INTEGER,ALLOCATABLE,DIMENSION(:) :: indexi, indexj

  ALLOCATE(indexi(9*natblocs*natblocs))
  ALLOCATE(indexj(9*natblocs*natblocs))
  ALLOCATE(a(3*natblocs,3*natom))
  ALLOCATE(hw(9*natblocs*natblocs))

  nblocs=nb
  
  open(unit=50,file=fmtx,status='old',form='formatted')
    
  open(unit=51,file='diagrtb_work.matblocs',status='unknown', &
       form='unformatted')
  
! Creation of diagrtb_work.matblocs
  
  WRITE(6,'(A)')'Rewriting the matrix'
  
  WRITE(51) natom,NB,(TABB(k),k=1,NB)
  nlign=1
  nempty=0
  
  NI1=1
  NTOTN1=0
  nn=1
  NI2=TABB(nn)
  qinit=.true.
  TBLOC=0
  
  imax=-1
  trace=0.d0
  ALLOCATE(qtrace(3*natom))
  qtrace=.false.

  nread=1
  DO
     
     IF (qinit) THEN
        DO i=1,TABB(nn)
           kk=NTOTN1+1
           DO j=kk,3*natom
              a(i,j)=0.d0
           END DO
        END DO
     END IF
     
     READ(50,*,IOSTAT=io) i,j,toto
     IF (io > 0) THEN
        WRITE(6,'(/A,I9)') 'Reading matrix line ',nread
        STOP '*Wrong matrix*'
     ELSE IF (io < 0) THEN
        EXIT
     END IF
     
     nread=nread+1
     
     IF (i.gt.3*natom.or.j.gt.3*natom) THEN
        WRITE(6,'(/A,I6,A,I6,A/A,I6,A)') 'Matrix element ',i,' j= ',j, &
             ' found.',' More than ',natom,' atoms in matrix file !'
        STOP '*Wrong matrix or wrong coordinate file*'
     END IF
     
     IF (i.eq.j) THEN
        IF (.not.qtrace(i)) THEN
           trace=trace+toto
           qtrace(i)=.true. 
        END IF
     END IF
     
     IF (i.gt.imax) imax=i
     IF (j.gt.imax) imax=j
     
     IF (i.le.NI2) THEN
        iii=i-NTOTN1
        a(iii,j)=toto
        qinit=.false.
     END IF
     
     IF (i.gt.NI2) THEN
        
        backspace 50
        NJ1=NI1
        NJ2=NI2
          
        tt=0
        TBLOC=TBLOC+1
        BLOC=TBLOC
        
        DO JBLOC=1,NB
           nbnnul=0
           ii=0
           tt=tt+1
           
           IF (tt.eq.1) THEN
              DO i=1,TABB(nn)
                 DO j=kk,NI2   
                    a(j-NTOTN1,i+NTOTN1)=a(i,j)
                 END DO
              END DO
              DO i=1,TABB(nn)
                 DO j=NJ1,NJ2
                    IF (a(i,j).ne.0.d0) THEN
                       ii=ii+1
                       HW(ii)=0.d0
                       nbnnul=ii
                       HW(ii)=a(i,j)
                       indexi(ii)=i+NTOTN1
                       indexj(ii)=j
                    END IF
                 END DO
              END DO
              IF (nbnnul.gt.9*natblocs*natblocs) THEN
                 WRITE(6,'(/A/2(A,I12))') 'Too many matrix elements in a single block: ', &
                      ' Nbnnul= ',nbnnul,' Max= ',9*natblocs*natblocs
                 STOP '*Too large block*'
              END IF
              IF (nbnnul.le.0) nempty=nempty+1
              nlign=nlign+1
              WRITE(51) nbnnul,(indexi(ii),indexj(ii),HW(ii),ii=1,nbnnul) 
           END IF

           IF (tt.gt.1) THEN
              BLOC=BLOC+1
              NJ1=NJ2+1
              NJ2=NJ2+TABB(BLOC)
              DO i=1,TABB(nn)
                 DO j=NJ1,NJ2
                    IF (a(i,j).ne.0.d0) THEN
                       ii=ii+1
                       HW(ii)=0.d0
                       nbnnul=ii
                       HW(ii)=a(i,j)
                       indexi(ii)=i+NTOTN1
                       indexj(ii)=j
                    END IF
                 END DO
              END DO
              IF (nbnnul.gt.9*natblocs*natblocs) THEN
                 WRITE(6,'(/A/2(A,I12))') 'Too many matrix elements in a single block: ', &
                      ' Nbnnul= ',nbnnul,' Max= ',9*natblocs*natblocs
                 STOP '*Too large block*'
              END IF
              IF (nbnnul.le.0) nempty=nempty+1
              nlign=nlign+1
              WRITE(51) nbnnul,(indexi(ii),indexj(ii),HW(ii),ii=1,nbnnul)
           END IF
        END DO
        
        qinit=.true.
        NTOTN1=NTOTN1+TABB(nn)
        NB=NB-1
        nn=nn+1
        NI1=NI2+1
        NI2=NI2+TABB(nn)
     END IF
     
  END DO
  
  WRITE(6,'(/I13,A)') nread,' matrix lines read.'
  
  WRITE(6,'(/A,I10)') 'Matrix order    = ',imax
  WRITE(6,'(/A,F15.4)') 'Matrix trace    = ',trace
  
  IF (imax/3.ne.natom) THEN
     WRITE(6,'(/I6,A/A,I6,A)') imax, &
          ' coordinates found in this matrix, ', &
          ' instead of ',natom,'*3, as expected. Not the right one ?'
     STOP '*Wrong matrix or wrong coordinate file*'
  END IF
  
  WRITE(6,'(A,2I7,F15.4)') 'Last element read: ',i,j,toto
  
  iii=i-NTOTN1
  a(iii,j)=toto
  NJ1=NI1
  NJ2=NI2
  ii=0
  DO i=1,TABB(nn)
     DO j=kk,NJ2
        a(j-NTOTN1,i+NTOTN1)=a(i,j)
     END DO
  END DO
  DO i=1,TABB(nn)
     DO j=NJ1,NJ2
        IF (a(i,j).ne.0.d0) THEN
           ii=ii+1
           HW(ii)=0.d0
           nbnnul=ii
           HW(ii)=a(i,j)
           indexi(ii)=i+NTOTN1
           indexj(ii)=j
        END IF
     END DO
  END DO
  IF (nbnnul.gt.9*natblocs*natblocs) THEN
     WRITE(6,'(/A/2(A,I12))') 'Too many matrix elements in a single block: ', &
          ' Nbnnul= ',nbnnul,' Max= ',9*natblocs*natblocs
     STOP '*Too large block*'
  END IF
  IF (nbnnul.le.0) nempty=nempty+1
  
  nlign=nlign+1
  WRITE(51) nbnnul,(indexi(ii),indexj(ii),HW(ii),ii=1,nbnnul)
  
  nb=nblocs 
  
  WRITE(6,'(I13,A)') nlign,' lines saved.'
  WRITE(6,'(I13,A)') nempty,' empty lines.'
  
  IF (nlign.ne.nb*(nb+1)/2+1) THEN
     WRITE(6,'(/I6,A)') nb*(nb+1)/2+1,' lines are expected on output !'
     STOP '*Unexpected big problem*'
  END IF
  
  WRITE(6,'(A)') 'Number of lines on output is as expected.'
  
  close(50)
  close(51)
  
  DEALLOCATE(qtrace)
  DEALLOCATE(indexi)
  DEALLOCATE(indexj)
  DEALLOCATE(a)
  DEALLOCATE(hw)

END SUBROUTINE PREPMAT
!----------------------------------------------------------------------
SUBROUTINE RTB ( natom, nb, nbrt, nddres, nddres2, corlin, amass)
  USE nrutil
  IMPLICIT NONE
  !     Purpose:
  !     =======
  !
  !     ------------------
  !     Matrix projection.
  !     ------------------
  !
  !     Arguments
  !     =========
  !
  !     natom   : number of atoms in the system.
  !     nb      : number of blocks the system is split into.
  !     nbrt    : number of degrees of freedom of each block.
  !     nddres  : size of largest block.
  !     nddres2 : nddres**2
  !     indexi  : i-index of non-zero elements in block hw.
  !     indexj  : j-index of non-zero elements in block hw.
  !     n       : block lengths.
  !     amass   : masses.
  !     corlin  : coordinates.
  !     d       : tranlation and rotation vectors of each block.
  !     h       : block being considered.
  !     hd      : working matrix for projection of a given block.
  !     hh      : working matrix for projection of a given block.
  !     mat     : projected matrix.
  !     rt      : tranlation and rotation vectors of all blocks.
  !     s       : working matrix for projection of a given block.
  !
  !-----------------------------------------------------------------------
  INTEGER :: natom, nb, nbrt, nddres, nddres2
  INTEGER :: indexi(nddres2), indexj(nddres2), n(nb)
  REAL(DP) :: amass(3*natom), corlin(3*natom)
  REAL(DP),ALLOCATABLE,DIMENSION(:) :: h
  REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: d, hd, hh, mat, s 
  REAL(DP),ALLOCATABLE,DIMENSION(:,:,:) :: rt
  LOGICAL :: qok
  INTEGER :: i, ibloc, ii, indxi, indxj, j, jbloc, k, kk, kr, &
       mmm, natcur, nbcur, nbnnul, nempty, ni, nj, &
       nlign, nnskip, nnuls, nskip, ntot, ntoti, ntotj, io
  REAL(DP) :: amasstot, trace, xg, yg, zg
  
  ALLOCATE(h(nddres2))
  ALLOCATE(d(nbrt,nddres))
  ALLOCATE(hd(nddres,nbrt))
  ALLOCATE(hh(nddres,nddres))
  ALLOCATE(mat(6*nb,6*nb))
  ALLOCATE(s(nbrt,nbrt))

  OPEN(unit=51,file='diagrtb_work.matblocs',status='old', &
       form='unformatted')
  
  OPEN(unit=52,file='diagrtb_work.sdijb',status='unknown', &
       form='unformatted')
  OPEN(unit=35,file='diagrtb_work.blocs',status='unknown', &
       form='unformatted')
  
  qok=.true.
  AMASSTOT = 0.D0 
  DO I=1,3*natom
     AMASSTOT = AMASSTOT + AMASS(I)
     IF (amass(i).lt.0) qok=.false.
  END DO
  AMASSTOT = AMASSTOT / 3.D0 
  
  IF (.not.qok) THEN
     WRITE(6,'(A)') 'Masses with negative values read in temporary file.'
     STOP '*Unexpected big problem*'
  END IF
  
  WRITE(6,'(A,F16.4)') 'Total mass = ',AMASSTOT
  !
  natcur=natom
  nbcur=nb
  READ(51) natom,NB,(n(k),k=1,NB)
  
  WRITE(6,'(A,I6)') 'Number of atoms found in matrix:',natom
  
  IF (natom.ne.natcur) THEN
     WRITE(6,'(A)') 'Pdb file and matrix are not consistent any more !'
     STOP '*Unexpected big problem*'
  END IF
  
  WRITE(6,'(A,I6)') 'Number of blocks =',nbcur

  IF (nb.ne.nbcur) THEN
     WRITE(6,'(/I6,A)') nbcur,' blocks... up to now !'
     STOP '*Unexpected big problem*'
  END IF
  WRITE(35) natom,NB,(n(k),k=1,NB)
  
  natcur=0
  DO i=1,nb
     natcur=natcur+n(i)
  END DO
  IF (natcur.ne.3*natom) THEN
     WRITE(6,'(A,I6,A,I6)') 'Sum of block lengths: ',natcur, &
          ' instead of ',3*natom
     STOP '*Unexpected big problem*'
  END IF
  
  nlign=1
  nempty=0
  nskip=0
  
  WRITE(6,'(A)')'Projection started'

  kk=0
  NTOT=0
  ALLOCATE(rt(nbrt,nddres,nb))

  DO IBLOC=1,NB  

     nbnnul=0
     
     READ(51,IOSTAT=io)  &
          nbnnul,(indexi(ii),indexj(ii),H(ii),ii=1,nbnnul)
     IF (io > 0) THEN
        WRITE(6,'(A,I3,A,I3)') 'Error while reading I-bloc: ',ibloc, &
             ' nbnnul= ',nbnnul
        WRITE(6,'(I6,A)') nlign,' lines read.'
        WRITE(6,'(I6,A)') nempty,' empty lines found.'
        WRITE(6,'(I6,A)') nskip,' lines skipped.'
        STOP '*Matrix i/o error*'
     ELSE IF (io < 0) THEN
        WRITE(6,'(A,I3,A,I3)') 'End-of-file while reading I-bloc: ',ibloc, &
             ' nbnnul= ',nbnnul
        WRITE(6,'(I6,A)') nlign,' lines read.'
        WRITE(6,'(I6,A)') nempty,' empty lines found.'
        WRITE(6,'(I6,A)') nskip,' lines skipped.'
        STOP '*Matrix i/o error*'
     ELSE
        nlign=nlign+1
     END IF
     
     IF (nbnnul.gt.nddres2) THEN
        WRITE(6,'(I6,A,I6,A,I6)') nbnnul,' elements in bloc ',ibloc, &
             ' Maximum allowed: ',nddres2
        STOP '*Too large block*'
     END IF
     IF (nbnnul.le.0) THEN
        nempty=nempty+1
     ELSE
        DO i=1,nbnnul
           IF (indexi(i).le.NTOT.or.indexj(i).le.NTOT) THEN
              WRITE(6,'(5(A,I6))') 'Ntot=',ntot,' but i= ',indexi(i),' j= ', &
                   indexj(i),' for element ',i,' of bloc ',ibloc 
              STOP '*Unexpected big problem*'
           END IF
        END DO
     END IF
     
     kk=kk+1
     DO kr=1,NB-kk
        READ(51,IOSTAT=io) nnskip
        IF (io > 0) THEN
           WRITE(6,'(3(A,I3))') 'Error after reading I-bloc: ',ibloc,  &
                ' NB= ',NB,' kk= ',kk 
           WRITE(6,'(I6,A)') nlign,' lines read.'
           WRITE(6,'(I6,A)') nempty,' empty lines found.'
           WRITE(6,'(I6,A)') nskip,' lines skipped.'
           STOP '*Matrix i/o error*'
        ELSE IF (io < 0) THEN
           WRITE(6,'(A,I3,A,I3)') 'Error while reading I-bloc: ',ibloc, &
                ' nbnnul= ',nbnnul
           WRITE(6,'(I6,A)') nlign,' lines read.'
           WRITE(6,'(I6,A)') nempty,' empty lines found.'
           WRITE(6,'(I6,A)') nskip,' lines skipped.'
           STOP '*Matrix i/o error*'
        ELSE
           IF (nnskip.le.0) nempty=nempty+1
           nskip=nskip+1
        END IF
     END DO
     
     IF (n(ibloc).gt.nddres.or.n(ibloc).le.0) THEN
        WRITE(6,'(A,I6,A,I6,A,I6)') 'N(ibloc)= ',n(ibloc), &
             ' for bloc ',ibloc, &
             ' Maximum allowed: ',nddres
        STOP '*Too large block*'
     END IF
     
     DO i=1,N(IBLOC)
        DO j=1,N(IBLOC)
           HH(i,j)=0.d0
        END DO
     END DO
     
     DO ii=1,nbnnul
        i=indexi(ii)-NTOT
        j=indexj(ii)-NTOT
        HH(i,j)=H(ii)
     END DO
     
     AMASSTOT = 0.D0
     DO I=1,N(IBLOC)
        AMASSTOT = AMASSTOT + AMASS(NTOT+I)
     END DO
     AMASSTOT=AMASSTOT/3.d0
     
     !
     XG = 0.D0
     YG = 0.D0  
     ZG = 0.D0
     
     DO I=0,( (N(IBLOC)/3)-1 )
        XG = XG + AMASS(NTOT+3*I+1)*corlin(NTOT+3*I+1)
        YG = YG + AMASS(NTOT+3*I+2)*corlin(NTOT+3*I+2)
        ZG = ZG + AMASS(NTOT+3*I+3)*corlin(NTOT+3*I+3)
     END DO
     
     XG = XG / AMASSTOT
     YG = YG / AMASSTOT
     ZG = ZG / AMASSTOT
     !
     DO I=0,( (N(IBLOC)/3)-1)      
        corlin(NTOT+3*I+1)= corlin(NTOT+3*I+1) - XG
        corlin(NTOT+3*I+2)= corlin(NTOT+3*I+2) - YG
        corlin(NTOT+3*I+3)= corlin(NTOT+3*I+3) - ZG
     END DO
     
     DO I=1,6
        DO J=1,N(IBLOC)
           D(I,J) = 0.D0 
        END DO
     END DO
     !
     DO I=1, N(IBLOC)/3
        II=NTOT 
        D(1,1+3*(I-1))= DSQRT(AMASS(II+1+3*(I-1)))
        D(2,2+3*(I-1))= DSQRT(AMASS(II+2+3*(I-1)))
        D(3,3+3*(I-1))= DSQRT(AMASS(II+3+3*(I-1)))
        D(4,2+3*(I-1))=-DSQRT(AMASS(II+2+3*(I-1)))*corlin(II+3+3*(I-1))
        D(4,3+3*(I-1))= DSQRT(AMASS(II+3+3*(I-1)))*corlin(II+2+3*(I-1))
        D(5,1+3*(I-1))= DSQRT(AMASS(II+1+3*(I-1)))*corlin(II+3+3*(I-1))
        D(5,3+3*(I-1))=-DSQRT(AMASS(II+3+3*(I-1)))*corlin(II+1+3*(I-1))
        D(6,1+3*(I-1))=-DSQRT(AMASS(II+1+3*(I-1)))*corlin(II+2+3*(I-1))
        D(6,2+3*(I-1))= DSQRT(AMASS(II+2+3*(I-1)))*corlin(II+1+3*(I-1))
     END DO
     !
     
     MMM=6
     
     CALL SCHMIDT(MMM,N(IBLOC),D ) 
     !
     DO i=1,6
        DO j=1,N(IBLOC)
           rt(i,j,ibloc)=D(i,j)
           WRITE(35) i,j,ibloc,rt(i,j,ibloc)
        END DO
     END DO
     
     NTOT=NTOT+N(IBLOC)
  END DO
  
  CLOSE(51)
  
  OPEN(unit=51,file='diagrtb_work.matblocs',status='old', &
       form='unformatted')
  
  READ(51) natom,NB,(N(k),k=1,NB)
  
  ni=0
  nj=0
  indxi=0
  indxj=0
  NTOTI=0
  NTOTJ=0
  
  DO IBLOC=1,NB
     DO JBLOC=IBLOC,NB
        
        READ(51,IOSTAT=io)  &
             nbnnul,(indexi(ii),indexj(ii),H(ii),ii=1,nbnnul)
        IF (io > 0) THEN
           WRITE(6,'(A,I3,A,I3)') 'Error while reading J-bloc: ',jbloc, &
                ' for band: ',ibloc 
           WRITE(6,'(I6,A)') nlign,' lines read.'
           WRITE(6,'(I6,A)') nempty,' empty lines found.'
           WRITE(6,'(I6,A)') nskip,' lines skipped.'
           STOP '*Matrix i/o error*'
        ELSE IF (io < 0) THEN
           WRITE(6,'(A,I3,A,I3)') 'Error while reading J-bloc: ',jbloc, &
                ' for band: ',ibloc
           WRITE(6,'(I6,A)') nlign,' lines read.'
           WRITE(6,'(I6,A)') nempty,' empty lines found.'
           WRITE(6,'(I6,A)') nskip,' lines skipped.'
           STOP '*Matrix i/o error*'
        END IF
        
        
        DO i=1,N(IBLOC)
           DO j=1,N(JBLOC)
              HH(i,j)=0.d0
           END DO
        END DO
        
        DO ii=1,nbnnul
           i=indexi(ii)-NTOTI
           j=indexj(ii)-NTOTJ
           HH(i,j)=H(ii)
        END DO
        
        
        DO j=1,N(IBLOC)
           DO i=1,6
              HD(j,i)=0.d0
              DO k=1,N(JBLOC)
                 HD(j,i)=HD(j,i)+HH(j,k)*rt(i,k,jbloc)   
              END DO
           END DO
        END DO
        
        DO i=1,6
           DO j=1,6
              s(i,j)=0.d0
              DO k=1,N(IBLOC)
                 s(i,j)=s(i,j)+rt(i,k,ibloc)*HD(k,j)
              END DO
              ni=i+indxi
              nj=j+indxj
              mat(ni,nj)=0.d0
              mat(ni,nj)=s(i,j)
           END DO
        END DO
        
        
        NTOTJ=NTOTJ+N(JBLOC)
        
        indxj=indxj+6
        
     END DO
     
     indxi=indxi+6
     indxj=indxi 
     NTOTI=NTOTI+N(IBLOC)
     NTOTJ=NTOTI
     
  END DO
  
  DEALLOCATE(rt)

  WRITE(6,'(A)') 'Projected matrix is being saved'
  
  trace=0.d0
  nnuls=0
  DO i=1,6*nb
     DO j=i,6*nb
        IF (mat(i,j).ne.0.d0) THEN
           nnuls=nnuls+1
           WRITE(52) i,j,mat(i,j)
           IF (i.eq.j) trace=trace+mat(i,j)
        END IF
     END DO
  END DO
  
  CLOSE(35)
  CLOSE(51)
  CLOSE(52)

  DEALLOCATE(h)
  DEALLOCATE(d)
  DEALLOCATE(hd)
  DEALLOCATE(hh)
  DEALLOCATE(mat)
  DEALLOCATE(s)
  
  WRITE(6,'(A,F12.4)') 'Projected matrix trace = ',trace
  
  WRITE(6,'(I10,A)') nnuls,' non-zero elements.'
    
END SUBROUTINE RTB

! -----------------------------------------------
SUBROUTINE SCHMIDT(M,N,C)
  ! ----------------------------------------------- 
  USE nrutil
  IMPLICIT NONE
  !
  ! GRAM-SCHMIDT
  !
  INTEGER ::   I,J,K,M,N
  REAL(DP) :: aaa, anorm, C(6,*), REC(6,6)
  !
  ANORM = 0.D0
  DO I=1,N
     ANORM = ANORM + C(1,I)*C(1,I)
  END DO
  ANORM = 1.D0 / ( DSQRT(ANORM) )
  !
  DO I=1,N
     C(1,I) = ANORM * C(1,I)
  END DO
  ! 
  DO  I=2,M      
     DO J=1,I-1
        REC(J,I) = 0.D0
        DO K=1,N
           REC(J,I) = REC(J,I) + C(J,K)*C(I,K)
        END DO
     END DO
     DO K=1,N
        AAA = 0.D0
        DO J=1,I-1
           AAA = AAA + C(J,K)*REC(J,I)
        END DO
        C(I,K) = C(I,K) - AAA
     END DO
     !
     ANORM = 0.D0
     DO K=1,N
        ANORM = ANORM + C(I,K)*C(I,K)
     END DO
     ANORM = 1.D0 / ( DSQRT(ANORM) )
     !
     DO K=1,N
        C(I,K) = ANORM*C(I,K)
     END DO
     !
  END DO

END SUBROUTINE SCHMIDT
!----------------------------------------------------------------------
SUBROUTINE RTBTOMODES (natom, nb, nbrt, nddres, nvec, feig)
  USE nrutil
  USE write_files
  USE read_files
  IMPLICIT NONE
  !     Purpose:
  !     =======
  !
  !     -----------------------------------------------------------
  !     From eigenvectors in block-rotation-translation coordinates
  !     to eigenvectors in cartesian coordinates.
  !     -----------------------------------------------------------
  !
  !     Input : matrix file 'diagrtb_work.eigenfacs'
  !     Output: matrix file feig
  !
  !     Arguments
  !     =========
  !
  !     natom     : maximum number of atoms in the system.
  !     nresmax   : maximum number of residues and of blocks.
  !     ntab      : working array.
  !     resi      : residue number each atom belongs to.
  !     tabb      : block lengths.
  !     amass     : coordinate masses.
  !     corlin    : coordinates.
  !     massat    : atomic masses.
  !     xat       : atomic x-coordinate.
  !     yat       : atomic y-coordinate.
  !     zat       : atomic z-coordinate.
  !     natblocs  : maximum number of atoms found in blocks.
  !     natom     : number of atoms in the system.
  !     nb        : number of blocks the system is split into.
  !
  !-----------------------------------------------------------------------
  INTEGER :: natom, nb, nbrt, nddres, nvec
  REAL(DP),ALLOCATABLE,DIMENSION(:,:,:) :: rt
  REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: matvec1, covec
  REAL(DP),ALLOCATABLE,DIMENSION(:) :: freq1, norm
  CHARACTER :: feig*120, cformat*40, cstatus*40, nomeig1*40, nomeig2*40 ,filename*120
  LOGICAL :: qinter
  INTEGER :: i, ibloc, ii, ivec, j, k, natcur, nbcur, nbread, nddl1, io, &
       ntot, nunit, uneig1, uneig2, natomhold
  REAL(DP) :: dvec
  INTEGER,ALLOCATABLE,DIMENSION(:) :: n

  nunit=10
  OPEN(unit=35,file='diagrtb_work.blocs',status='old', &
       form='unformatted')
  
  natcur=natom
  nbcur=nb

  READ(35) natom,nb
  ALLOCATE(n(nb))
  REWIND(35)
  READ(35) natom, nb, (n(i),i=1,nb)
 
  WRITE(6,'(/A,I6)') 'Number of atoms in temporary block-file = ',natom
  
  IF (natom.ne.natcur) THEN
     WRITE(6,'(/I6,A)') natcur,' atoms... up to now.'
     STOP '*Unexpected big problem*'
  END IF
  
  WRITE(6,'(A,I6)') 'Number of blocs = ',nb
  
  IF (nb.ne.nbcur) THEN
     WRITE(6,'(/I6,A)') nbcur,' blocks... up to now.'
     STOP '*Unexpected big problem*'
  END IF
  
  filename='diagrtb_work.eigenfacs'
  CALL read_eigenfacs(filename,natomhold,1,6*natom,3,freq1,matvec1,nvec)

  nddl1=natomhold*3

  WRITE(6,'(/I5,A,I6,A)') nvec,' vectors, with ',nddl1,' coordinates in vector file.'
  
  IF (nddl1.ne.6*nb) THEN
     WRITE(6,'(/A,I6,A)') 'Vectors of length: ',6*nb, &
          ' were expected. Block and vector files are not consistent.'
     STOP '*Unexpected big problem*'
  END IF
  !
  ALLOCATE(norm(nvec))

!  WRITE(6,'(/A)') ' Norm of eigenvectors in projected coordinates (one expected):'  
!  DO ivec=1,nvec
!     norm(ivec)=0.d0
!     DO i=1,nddl1
!        dvec=matvec1(i,ivec)**2
!        norm(ivec)=norm(ivec)+dvec
!     END DO
!     norm(ivec)=dsqrt(norm(ivec))
!  END DO
!  WRITE(6,'(5F8.5)') (norm(ivec),ivec=1,nvec)
  
  !
  WRITE(6,'(/A)') 'RTB block-file is being read.'
  !
  ALLOCATE(rt(nbrt,nddres,nb))
  rt=0.d0
  nbread=0
  !
  DO 
     READ(35,IOSTAT=io) i,j,ibloc,rt(i,j,ibloc)
     IF (io < 0) THEN
        EXIT
     END IF
     nbread=nbread+1
     IF (i.gt.nbrt.or.i.le.0) THEN
        WRITE(6,'(I6,A,I1,A,I4)') i, &
             ' rigid body ddl instead of ',nbrt,' for block ',ibloc
        STOP '*Unexpected big problem*'
     END IF
     IF (j.gt.nddres.or.j.le.0) THEN
        WRITE(6,'(I6,A,I6,A,I6)') j, &
             ' ddl for block ',ibloc,'. Maximum allowed is: ',nddres
        STOP '*Unexpected big problem*'
     END IF
     IF (ibloc.gt.nb.or.ibloc.le.0) THEN
        WRITE(6,'(I4,A,I4)') ibloc, &
             ' blocks at least. Maximum allowed is: ',nb
        STOP '*Unexpected big problem*'
     END IF
  END DO
  
  WRITE(6,'(I9,A)')nbread,' lines found in RTB file.'
  
  ALLOCATE(covec(3*natom,nvec))
  covec=0.d0
  
  DO ivec=1,nvec
     ntot=0
     ii=0
     DO ibloc=1,nb
        DO j=1,n(ibloc)
           covec(j+ntot,ivec)=0.d0
           DO i=1,6
              covec(j+ntot,ivec)=covec(j+ntot,ivec)+ &
                   rt(i,j,ibloc)*matvec1(i+ii,ivec)
           END DO
        END DO
        ii=ii+6
        ntot=ntot+n(ibloc)
     END DO
  END DO
  
  DEALLOCATE(n)
  DEALLOCATE(rt)
  DEALLOCATE(matvec1)

!  WRITE(6,'(/A)')  &
!       ' Norm of eigenvectors in cartesian coordinates (one expected):'
!  DO ivec=1,nvec
!     norm(ivec)=0.d0
!     DO i=1,natom
!        ii=3*i-2
!        dvec=covec(ii,ivec)**2+covec(ii+1,ivec)**2 &
!             +covec(ii+2,ivec)**2
!        norm(ivec)=norm(ivec)+dvec
!     END DO
!     norm(ivec)=dsqrt(norm(ivec))
!  END DO
!  WRITE(6,'(5F8.5)') (norm(ivec),ivec=1,nvec)
!  
!  WRITE(6,'(/A)')  &
!       ' Orthogonality of first eigenvectors (zero expected):'
!  DO i=2,min(10,nvec)
!     DO j=1,i-1
!        norm(j)=0.d0
!        DO k=1,natom
!           ii=3*k-2
!           norm(j)=norm(j)+ &
!                covec(ii,i)*covec(ii,j)+covec(ii+1,i)*covec(ii+1,j)+ &
!                covec(ii+2,i)*covec(ii+2,j)
!        END DO
!     END DO
!     WRITE(6,'(A,I3,A,10F6.3)')  &
!          ' Vector ',i,':',(norm(j),j=1,i-1)
!  END DO
  DEALLOCATE(norm)

  CALL write_eigenfacs(feig,covec,freq1,nvec,3*natom)
  
  DEALLOCATE(covec)
  DEALLOCATE(freq1)

  CLOSE(35)
  
  WRITE(6,'(/I6,A)') nvec, ' eigenvectors saved.'

END SUBROUTINE RTBTOMODES
!
!-----------------------------------------------------------------------------

SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
  WRITE(iunit,'(A)')"                           D  I  A  G  R  T  B                           "
  WRITE(iunit,'(A)')"                               VERSION 1.0                               "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"                               Written by:                               "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program calculates the normal modes for a ENM using the rotational  " 
  WRITE(iunit,'(A)')"translational block method based on the reduced Hermitian matrix file,   " 
  WRITE(iunit,'(A)')"matrix.sdijf, produced from the ENM program.                             "
  WRITE(iunit,'(A)')"The ENM can be split by residues, by the secondary structure (SECO) using"
  WRITE(iunit,'(A)')"the HELIX and SHEET comands in the pdb, or by custom domains (DOMN) read "
  WRITE(iunit,'(A)')"from input.domain file.                                                  "
  WRITE(iunit,'(A)')"A matrix.eigenfacs file is produced which can be further analysed using  "
  WRITE(iunit,'(A)')"the programs in this toolbox.                                            "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"     Usage:                                                              "
  WRITE(iunit,'(A)')"             diagrtb -pdb pdbfile [-i matrix] [-r num] [-str option]     "
  WRITE(iunit,'(A)')"                   [-mass] [-het] [-e num] [-high] [-eig] [-mat]         "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"Option    Type       Value           Description                         "
  WRITE(iunit,'(A)')"------------------------------------------------------------             "
  WRITE(iunit,'(A)')" -pdb     Input                      pdb file                            "
  WRITE(iunit,'(A)')"  -i      Input,Opt  matrix.sdijf    Reduced Hermitian matrix file       "
  WRITE(iunit,'(A)')"  -r      Input,Opt  1               Number of residues per block        "
  WRITE(iunit,'(A)')" -str     Input,Opt  NONE            Substructuring options, NONE, SECO, "
  WRITE(iunit,'(A)')"                                     or DOMN                             "
  WRITE(iunit,'(A)')" -mass      Opt                      Use true atom masses                "
  WRITE(iunit,'(A)')" -res       Opt                      Use atom residue mass               " 
  WRITE(iunit,'(A)')" -het       Opt                      Also read HETATM records            "
  WRITE(iunit,'(A)')"  -e      Input,Opt  56              Number of eigenvectors output       "
  WRITE(iunit,'(A)')" -eig       Opt                      Print eigenvalues for RTBs          "
  WRITE(iunit,'(A)')" -mat       Opt                      Print RTB matrix                    "
  WRITE(iunit,'(A)')"-covar      Opt                      Run with covar matrix input         "
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

FUNCTION FindDet(matrix, n)
  USE nrtype
  IMPLICIT NONE
  REAL(DP), DIMENSION(n,n) :: matrix
  INTEGER, INTENT(IN) :: n
  REAL(DP) :: m, temp, FindDet
  INTEGER :: i, j, k, l
  LOGICAL :: DetExists = .TRUE.
  l = 1
  !Convert to upper triangular form
  DO k = 1, n-1
     IF (matrix(k,k) == 0) THEN
        DetExists = .FALSE.
        DO i = k+1, n
           IF (matrix(i,k) /= 0) THEN
              DO j = 1, n
                 temp = matrix(i,j)
                 matrix(i,j)= matrix(k,j)
                 matrix(k,j) = temp
              END DO
              DetExists = .TRUE.
              l=-l
              EXIT
           ENDIF
        END DO
        IF (DetExists .EQV. .FALSE.) THEN
           FindDet = 0
           return
        END IF
     ENDIF
     DO j = k+1, n
        m = matrix(j,k)/matrix(k,k)
        DO i = k+1, n
           matrix(j,i) = matrix(j,i) - m*matrix(k,i)
        END DO
     END DO
  END DO
  
  !Calculate determinant by finding product of diagonal elements
  FindDet = l
  DO i = 1, n
     FindDet = FindDet * matrix(i,i)
  END DO
  
END FUNCTION FindDet

FUNCTION is_numeric(string)
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) :: string
LOGICAL :: is_numeric
REAL :: x
INTEGER :: e,n
CHARACTER(len=12) :: fmt

n=LEN_TRIM(string)
WRITE(fmt,'("(F",I0,".0)")') n
READ(string,fmt,IOSTAT=e) x
is_numeric = e == 0
END FUNCTION is_numeric

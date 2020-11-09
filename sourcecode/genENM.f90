PROGRAM genENM
  USE nrutil
  USE utils
  USE read_files
  USE write_files
  USE foul
  IMPLICIT NONE
  LOGICAL :: getoption, hetatm, caonly, resmass, atmass, cutvect,cutassign, &
       forceres,ukn,usesec,qexist,lig1,firstat, custb, dna, hin, is_numeric 
  CHARACTER :: dummy*120, pdbfile*120, lign80*80, nomatm*120, res2*4, nomres*120, &
       foutname*120,custbfile*120
  INTEGER :: natom, io, i, j, natomold, nident, ncusres,ii,ijjjj,iseed,jat,jj, &
       ll,nntr,nnzero,nhel,nshe,k,l,fat,cbs,cbnum,nseconds
  REAL(DP) :: cutoff, cutoffdef, rave, rdev, rmin, rmax,anmp,ddf,rkh,dist,dist2,dmax,dmin, &
       distave,drms,kij,kset,shift,rx,ry,rz,trace,random,masstol,xav,yav,zav,bfacav,entot
  LOGICAL,ALLOCATABLE,DIMENSION(:) :: otherres
  INTEGER,ALLOCATABLE,DIMENSION(:) :: atnum,resnum,resnumold, resone, restwo
  REAL(DP),ALLOCATABLE,DIMENSION(:) :: x,y,z,occ,bfac,vals,mass,massold,cutvalue, &
       acutoff,kijcust
  CHARACTER,ALLOCATABLE,DIMENSION(:) :: atom*6,name*4,res*3,chain*1,elem*2,chag*2, &
       chainold*1,ident*4,cshe*1,chel*1,chainone*1,chaintwo*1,seconds*80
  CHARACTER,ALLOCATABLE,DIMENSION(:,:) :: custchain
  REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: a
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: she, hel, custres

  IF (getoption('-help',.false.,dummy)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF

  iseed=27041961
  shift=1e-8
!-----------------------------------------------------------------------------------
! Read in options

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

  IF (getoption('-dna',.false.,dummy)) THEN
     dna=.true.
  ELSE
     dna=.false.
  END IF

  IF (getoption('-ca',.false.,dummy)) THEN
     caonly=.true.
  ELSE
     caonly=.false.
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

  IF (getoption('-c',.true.,dummy)) THEN
     READ(dummy,*) cutoffdef
  ELSE
     cutoffdef=12
     WRITE(6,'(A/A)') "Cut off set to 12 A","Assign with -c if a different value is needed"
  END IF 

  IF (getoption('-ccust',.true.,nomatm)) THEN
     cutvect=.true.
     INQUIRE(file=nomatm,exist=qexist)
     IF (.not.qexist) THEN
        STOP "ccust file specified does not exist"
     END IF 
  ELSE
     cutvect=.false.
  END IF

  IF (getoption('-f',.true.,dummy)) THEN
     READ(dummy,*) kset
  ELSE
     kset=1
     WRITE(6,'(A/A)') "Spring constant set to 1 kcal/mol/A^2","Assign with -f if a different value is needed"
  END IF 

  IF (getoption('-fcust',.true.,nomres)) THEN
     forceres=.true.
     INQUIRE(file=nomres,exist=qexist)
     IF (.not.qexist) THEN
        STOP "fcust file specified does not exist"
     END IF
  ELSE
     forceres=.false.
  END IF

  IF (getoption('-spcust',.true.,custbfile)) THEN
     custb=.true.
     INQUIRE(file=custbfile,exist=qexist)
     IF (.not.qexist) THEN
        STOP "spcust file specified does not exist"
     END IF
  ELSE
     custb=.false.
  END IF

  IF (getoption('-hine',.true.,dummy)) THEN
     IF (is_numeric(dummy(1:1))) THEN
        READ(dummy,*) rkh
     ELSE
        STOP "The value of h following -hine must be a positive real number"
     END IF
  ELSE
     rkh=-1
  END IF 
  
  IF (getoption('-hin',.false.,dummy)) THEN
     hin=.true.
  ELSE
     hin=.false.
  END IF 

  IF (getoption('-an',.true.,dummy)) THEN
     IF (is_numeric(dummy(1:1))) THEN
        READ(dummy,*) anmp
     ELSE
        STOP "The value of p following -an must be a positive real number"
     END IF
  ELSE
     anmp=-1
  END IF 

  IF (getoption('-sec',.false.,dummy)) THEN
     usesec=.true.
  ELSE
     usesec=.false.
  END IF 

  IF (getoption('-lig1',.false.,dummy)) THEN
     lig1=.true.
  ELSE
     lig1=.false.
  END IF 

!------------------------------------------------------------------------------------
! Read in pdb file
  CALL read_pdb(pdbfile,hetatm,natom,atom,atnum,name,res,chain,resnum,x,y,z, &
       occ,bfac,elem,chag)

  WRITE(6,'(A,I5,A)') "Read in pdbfile, found", natom, " atoms"
! End of reading in pdb file
!----------------------------------------------------------------------------------------- 

  ALLOCATE(mass(natom))
! Assigning atom masses
  IF (atmass) THEN
     CALL atom_mass(natom,elem,mass)
     WRITE(6,'(A)') "Assigned atoms their true mass"
  ELSE
     mass=1
     WRITE(6,'(A)') "Atom masses all set to 1 amu"
  END IF
!-----------------------------------------------------------------------------------------
! Adjust to ca atoms only if needed

  IF (caonly) THEN
     firstat=.true.
     ALLOCATE(resnumold(natom))
     ALLOCATE(massold(natom))
     ALLOCATE(chainold(natom))
     j=0
     DO i=1,natom
        resnumold(i)=resnum(i)
        massold(i)=mass(i)
        chainold(i)=chain(i)
        IF ((index(name(i),'CA').gt.0) .or.(name(i).eq." P  ".and.index(res(i),'D').gt.0.and.dna) &
             .or.(name(i).eq." C4'".and.index(res(i),'D').gt.0.and.dna) &
             .or.(name(i).eq." C2 ".and.index(res(i),'D').gt.0.and.dna))THEN
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
           mass(j)=mass(i)
        ELSE IF (atom(i).eq."HETATM".and.hetatm) THEN
! Make each HETATM residue 1 average bead
           IF (lig1) THEN
              IF (atom(i-1).ne.'HETATM') firstat=.true.
              IF (firstat) THEN
                 j=j+1
                 xav=x(i)*mass(i)
                 yav=y(i)*mass(i)
                 zav=z(i)*mass(i)
                 bfacav=bfac(i)*mass(i)
                 masstol=mass(i)
                 fat=i
                 firstat=.false.
              ELSE
                 IF (resnum(i).ne.resnum(i-1).or.chain(i).ne.chain(i-1)) THEN
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
                    mass(j)=masstol/REAL(i-fat)
                    
                    j=j+1
                    xav=x(i)*mass(i)
                    yav=y(i)*mass(i)
                    zav=z(i)*mass(i)
                    bfacav=bfac(i)*mass(i)
                    masstol=mass(i)
                    fat=i
                 ELSE
                    xav=xav+x(i)*mass(i)
                    yav=yav+y(i)*mass(i)
                    zav=zav+z(i)*mass(i)
                    bfacav=bfacav+bfac(i)*mass(i)
                    masstol=masstol+mass(i)
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
                 mass(j)=masstol/REAL(i+1-fat)
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
              mass(j)=mass(i)
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
              mass(j)=masstol/REAL(i-fat)
           END IF
        END IF
     END DO
     
     natomold=natom
     natom=j

     WRITE(6,'(A,I5,A)') "Adjusted to Ca atoms only, found", natom, " atoms"

     IF (resmass) THEN
        IF (atmass) THEN
           CALL res_mass(natom,res,mass)
        ELSE
           DO i=1,natom
              mass(i)=0
              DO j=1,natomold
                 IF ((resnum(i).eq.resnumold(j)).and.(chain(i).eq.chainold(j))) THEN
                    mass(i)=mass(i)+mass(j)
                 END IF
              END DO
           END DO
        END IF
        
        WRITE(6,'(A)') "Adjusted Ca mass to residue mass"
        
     END IF

     DEALLOCATE(resnumold)
     DEALLOCATE(massold)
     DEALLOCATE(chainold)

  END IF

! End of adjusting to ca atoms only
!-----------------------------------------------------------------------------------------

 IF (usesec) THEN
    CALL get_second(pdbfile,nhel,hel,chel,nshe,she,cshe)
    j=0
    k=1
    l=1
     DO i=1,natom
        ! Check if the residue is part of a HELIX
        IF ((resnum(i).ge.hel(k,1).and.resnum(i).le.hel(k,2)) &
             .or.(atom(i).eq."HETATM".and.hetatm)) THEN
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
           mass(j)=mass(i)
        ELSEIF (resnum(i).gt.hel(k,2)) THEN
           IF (k.lt.nhel) THEN
              k=k+1
              IF (resnum(i).ge.hel(k,1).and.resnum(i).le.hel(k,2)) THEN
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
                 mass(j)=mass(i)
              END IF
           END IF
        END IF
        ! Or part of a SHEET
        IF(resnum(i).ge.she(l,1).and.resnum(i).le.she(l,2)) THEN
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
           mass(j)=mass(i)
        ELSEIF (resnum(i).gt.she(l,2)) THEN
           IF (l.lt.nshe) THEN
              l=l+1
              IF (resnum(i).ge.she(l,1).and.resnum(i).le.she(l,2)) THEN
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
                 mass(j)=mass(i)
              END IF
           END IF
        END IF

     END DO

     natom=j

     WRITE(6,'(A,I5,A)') "Adjusted to Secondary structure atoms only, using", natom, " atoms"

! Write secondary structure pdb
     foutname='sec_strut.pdb'
     CALL write_pdb(foutname,atom,atnum,name,res,chain,resnum,x,y,z,&
           occ,bfac,mass,elem,chag,natom,.false.,hetatm)

  END IF

  CALL read_second(pdbfile,nseconds,seconds)


! Write Caonly.pdb
  foutname='CAonly.pdb'
  CALL write_pdb_sec(foutname,atom,atnum,name,res,chain,resnum,x,y,z,&
           occ,bfac,mass,elem,chag,natom,.true.,hetatm,nseconds,seconds)

!---------------------------------------------------------------------------------------
! Find custom atom cutoffs if required
  IF (cutvect) THEN
     OPEN(file=nomatm,form='FORMATTED',status='OLD',unit=5698)

     nident=0
     DO
        READ(5698,'(A)',IOSTAT=io)  lign80
        IF (io > 0) THEN
           WRITE(*,*) 'Check input.  Something was wrong'
           STOP
        ELSEIF (io < 0) THEN
           EXIT
        ELSE
           nident = nident + 1
        ENDIF
     ENDDO
     
     REWIND(5698)
     
     ALLOCATE(ident(nident))
     ALLOCATE(cutvalue(nident))
     ALLOCATE(acutoff(natom))

     nident=0
     WRITE(6,'(A)') "Custom cutoff parameters are:"
     DO
        READ(5698,'(A)',IOSTAT=io)  lign80
        IF (io > 0) THEN
           WRITE(*,*) 'Check input.  Something was wrong'
           STOP
        ELSEIF (io < 0) THEN
           EXIT
        ELSE
           nident = nident + 1
           READ(lign80,'(A4,1X,F7.3)') ident(nident), cutvalue(nident)
           WRITE(6,'(A4,1X,F7.3)') ident(nident), cutvalue(nident)
        ENDIF
     ENDDO

     DO i=1,natom
        cutassign=.false.
        DO j=1,nident
           IF (name(i).eq.ident(j)) THEN
              acutoff(i)=cutvalue(j)
              cutassign=.true.
              EXIT
           END IF
        END DO
        IF (.not.cutassign) THEN
           acutoff(i)=cutoff
        END IF
     END DO


     CLOSE(5698)
     
     WRITE(6,'(A)') "Read in custom cutoff parameters"

     DEALLOCATE(ident)
     DEALLOCATE(cutvalue)

  END IF

! End of custom atom cutoffs
!-----------------------------------------------------------------------------------------  

! Find custom interaction pairs if required
  IF (forceres) THEN
     OPEN(file=nomres,form='FORMATTED',status='OLD',unit=5699)

     ncusres=0

     DO
        READ(5699,'(A)',IOSTAT=io) lign80
        IF (io > 0) THEN
            WRITE(*,*) 'resforce > check input'
            STOP
         ELSE IF (io < 0) THEN
            EXIT
         ELSE
            ncusres=ncusres+1
         END IF
      END DO
      
      ALLOCATE(resone(ncusres))
      ALLOCATE(restwo(ncusres))
      ALLOCATE(chainone(ncusres))
      ALLOCATE(chaintwo(ncusres))
      ALLOCATE(kijcust(ncusres))
      ALLOCATE(otherres(ncusres))
      
      REWIND(5699)

      ncusres=0

      WRITE(6,'(A)') "Custom residue interactions are:"
      DO
        READ(5699,'(A)',IOSTAT=io) lign80
        IF (io > 0) THEN
            WRITE(*,*) 'resforce > check input'
            STOP
         ELSE IF (io < 0) THEN
            EXIT
         ELSE
            ncusres=ncusres+1
            READ(lign80,'(1X,I4,1X,A1,1X,A4,1X,A1,1X,F8.3)') resone(ncusres), chainone(ncusres), res2, &
                 chaintwo(ncusres), kijcust(ncusres)
            WRITE(6,'(1X,I4,1X,A1,1X,A4,1X,A1,1X,F8.3)') resone(ncusres), chainone(ncusres), res2, &
                 chaintwo(ncusres), kijcust(ncusres)

            IF (index(res2,"*").gt.0) THEN
               otherres(ncusres)=.false.
               restwo(ncusres)=0
            ELSE
               otherres(ncusres)=.true.
               READ(res2,'(I4)') restwo(ncusres)
            END IF
            
         END IF

      END DO

     CLOSE(5699)
     
     WRITE(6,'(A)') "Read in custom interaction parameters"
     
  END IF

! End of custom interaction pairs
!-----------------------------------------------------------------------------------------
! Read in custom bonds
  IF (custb) THEN
     OPEN(file=custbfile,form='FORMATTED',status='OLD',unit=5799)

     cbs=0

     DO
        READ(5799,'(A)',IOSTAT=io) lign80
        IF (io > 0) THEN
            WRITE(*,*) 'custbonds > check input'
            STOP
         ELSE IF (io < 0) THEN
            EXIT
         ELSE
            cbs=cbs+1
         END IF
      END DO
      
      ALLOCATE(custres(2*cbs,2))
      ALLOCATE(custchain(2*cbs,2))
      
      REWIND(5799)

      cbs=0

      WRITE(6,'(A)') "Custom springs are:"
      DO
        READ(5799,'(A)',IOSTAT=io) lign80
        IF (io > 0) THEN
            WRITE(*,*) 'resforce > check input'
            STOP
         ELSE IF (io < 0) THEN
            EXIT
         ELSE
            cbs=cbs+1
            READ(lign80,'(2(1X,I4,1X,A1))') custres(cbs,1), custchain(cbs,1),custres(cbs,2), custchain(cbs,2)
            WRITE(6,'(2(1X,I4,1X,A1))') custres(cbs,1), custchain(cbs,1),custres(cbs,2), custchain(cbs,2)
            cbs=cbs+1
            custres(cbs,1)=custres(cbs-1,2)
            custchain(cbs,1)=custchain(cbs-1,2)
            custres(cbs,2)=custres(cbs-1,1)
            custchain(cbs,2)=custchain(cbs-1,1)
            WRITE(6,'(2(1X,I4,1X,A1))') custres(cbs,1), custchain(cbs,1),custres(cbs,2), custchain(cbs,2)
         END IF

      END DO

     CLOSE(5799)
     
     WRITE(6,'(A)') "Read in custom springs"
  END IF

! End of reading custom bonds
!-----------------------------------------------------------------------------------------

! Some stats
  WRITE(6,'(/A)')' Coordinate statistics: '

  CALL vecstat(x,natom,rmin,rmax,rave,rdev)
  WRITE(6,'(4(A,F12.6))')' <x>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax

  CALL vecstat(y,natom,rmin,rmax,rave,rdev)
  WRITE(6,'(4(A,F12.6))')' <y>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax

  CALL vecstat(z,natom,rmin,rmax,rave,rdev)
  WRITE(6,'(4(A,F12.6))')' <z>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax

  WRITE(6,'(/A)')' Mass statistics: '
  CALL vecstat(mass,natom,rmin,rmax,rave,rdev)
  WRITE(6,'(4(A,F12.6))')' <m>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax

!-----------------------------------------------------------------------------------------

  ALLOCATE(a(3,3*natom))

  OPEN(file='ENM.vmd',form='FORMATTED',unit=7432)

  WRITE(7432,'(A)') '#!/usr/local/bin/vmd'
  WRITE(7432,'(A)') '# script for VMD (Visual Molecular Dynamics)'
  WRITE(7432,'(A)') '# Goal: visualizing the elastic network'
  WRITE(7432,'(A)') '# Type: vmd -e this-file'
  WRITE(7432,'(A)') 'color Display {Background} white'
  WRITE(7432,'(A)') 'mol new'
  WRITE(7432,'(A)') 'draw color black'

  OPEN(file='matrix.sdijf',form='FORMATTED',unit=9432)

  trace=0.d0
  dmin=0.d0
  dmax=0.d0
  distave=0.d0
  drms=0.d0
  nnzero=0
  nntr=0
  ll=0
  entot=0.d0

  DO i=1,natom
     ii=3*i-2

     DO j=1,3*natom
        a(1,j)=0.d0
        a(2,j)=0.d0
        a(3,j)=0.d0
     END DO
     
     DO j=1,natom
        IF (i.ne.j) THEN
           jj=3*j-2
           
           !------------------------------------------------------------------------
           ! Set spring constant to custom between residues if required
           kij=kset
           IF (forceres) THEN
              DO ijjjj=1,ncusres
                 IF (((resone(ijjjj).eq.resnum(i)).and.(chainone(ijjjj).eq.chain(i))) &
                    .or.((resone(ijjjj).eq.resnum(j)).and.(chainone(ijjjj).eq.chain(j)))) THEN
                    IF (otherres(ijjjj)) THEN
                       IF (((restwo(ijjjj).eq.resnum(i)).and.(chaintwo(ijjjj).eq.chain(i))) &
                            .or.((restwo(ijjjj).eq.resnum(j)).and.(chaintwo(ijjjj).eq.chain(j)))) THEN
                          kij=kijcust(ijjjj)
                       END IF
                    ELSE
                       kij=kijcust(ijjjj)
                    END IF
                 END IF
              END DO
           END IF
           !-------------------------------------------------------------------

           rx=x(i)-x(j)
           ry=y(i)-y(j)
           rz=z(i)-z(j)
           dist2=rx*rx + ry*ry + rz*rz
           dist=sqrt(dist2)
           
           ! If a cutoff vector has been defined then the atom-atom cutoff 
           ! is taken here as the larger of the 2 values
           cutoff=cutoffdef
           IF (cutvect) THEN
              cutoff=max(acutoff(i),acutoff(j))
           END IF
           
           IF (hin) THEN
              IF (dist.lt.4.0) THEN
                 kij=kij*(205.54*dist-571.21)
              ELSE 
                 kij=kij*(305920/dist**6)
              END IF
           ! Adjust spring constant for Hinsen exp model
           ELSEIF (rkh.gt.0.d0) THEN
              kij=kij*exp(-(dist/rkh)**2.d0)   
           ! Adjust spring constant for anisotropic network model
           ELSEIF (anmp.gt.0.d0) THEN
              kij=kij/(dist**anmp)
           END IF

           ! Set custom bonds
           IF (custb) THEN
              DO cbnum=1,cbs
                 IF (resnum(i).eq.custres(cbnum,1).and.chain(i).eq.custchain(cbnum,1).and. &
                      resnum(j).eq.custres(cbnum,2).and.chain(j).eq.custchain(cbnum,2)) THEN
                    cutoff=9999.d9
                 END IF
              END DO
           END IF

           !--------------------------------------------------------------------------
           ! Calculation of the element harmonic potential
           IF (dist.le.cutoff) THEN 
              
              entot=entot+kij*dist2

              ll=ll+1
              IF (j.gt.i) THEN
                 WRITE(7432,'(A,3F12.4,A,3F12.4,A)') 'draw line {',x(i),y(i),z(i),'} {',x(j),y(j),z(j),'}'
              END IF
              
              IF (ll.eq.1.or.dist.lt.dmin) dmin=dist
              IF (ll.eq.1.or.dist.gt.dmax) dmax=dist
              
              distave=distave+dist
              drms=drms+dist2
              
              ! Diagonal elements of blocks i and j:
              !-----------------------------------
              ddf=kij/dist2
              a(1,ii)=a(1,ii)+rx*rx*ddf
              a(1,jj)=a(1,jj)-rx*rx*ddf
              a(2,ii+1)=a(2,ii+1)+ry*ry*ddf
              a(2,jj+1)=a(2,jj+1)-ry*ry*ddf
              a(3,ii+2)=a(3,ii+2)+rz*rz*ddf
              a(3,jj+2)=a(3,jj+2)-rz*rz*ddf
              
              ! Extra-diagonal elements of the two blocks:
              !---------------------------------------
              a(1,ii+1)=a(1,ii+1)+rx*ry*ddf
              a(2,ii)=a(2,ii)+rx*ry*ddf
              a(1,jj+1)=a(1,jj+1)-rx*ry*ddf
              a(2,jj)=a(2,jj)-rx*ry*ddf
              a(1,ii+2)=a(1,ii+2)+rx*rz*ddf
              a(3,ii)=a(3,ii)+rx*rz*ddf
              a(1,jj+2)=a(1,jj+2)-rx*rz*ddf
              a(3,jj)=a(3,jj)-rx*rz*ddf
              a(2,ii+2)=a(2,ii+2)+ry*rz*ddf
              a(3,ii+1)=a(3,ii+1)+ry*rz*ddf
              a(2,jj+2)=a(2,jj+2)-ry*rz*ddf
              a(3,jj+1)=a(3,jj+1)-ry*rz*ddf
           endif
        endif
     enddo
     
     ! Only the upper half of the matrix calculated
     
     ! Level-shift, to avoid the zero numerical values at 
     ! the time of the diagonalisation (minimization is perfect, 
     ! by definition). The chance is to raise the degeneration 
     ! of the six null eigenvalues, and to differentiate rotations 
     ! and translations.
     
     a(1,ii)  =a(1,ii)   + shift*random(iseed)
     a(2,ii+1)=a(2,ii+1) + shift*random(iseed)
     a(3,ii+2)=a(3,ii+2) + shift*random(iseed)
     
     DO j=ii,3*natom
        jat=(j-1)/3+1
        IF (a(1,j).ne.0.d0) THEN
           nnzero=nnzero+1
           WRITE(9432,'(2I10,1PG20.12)') ii,j,a(1,j)/sqrt(mass(i)*mass(jat))
        END IF
     END DO
     
     DO j=ii+1,3*natom
        jat=(j-1)/3+1
        IF (a(2,j).ne.0.d0) THEN
           nnzero=nnzero+1
           WRITE(9432,'(2I10,1PG20.12)')ii+1,j,a(2,j)/sqrt(mass(i)*mass(jat))
        END IF
     END DO
     
     DO j=ii+2,3*natom
        jat=(j-1)/3+1
        IF (a(3,j).ne.0.d0) THEN
           nnzero=nnzero+1
           WRITE(9432,'(2I10,1PG20.12)') ii+2,j,a(3,j)/sqrt(mass(i)*mass(jat))
        END IF
     END DO
     
     nntr=nntr+1
     
     trace=trace+(a(1,ii)+a(2,ii+1)+a(3,ii+2))/mass(i)
  END DO
  
  CLOSE(9432)
  
  WRITE(7432,'(2A)') 'mol load pdb ',pdbfile
  CLOSE(7432)
  
  WRITE(6,'(/A,F8.4,A)')' The matrix is ', 100.d0*dfloat(nnzero)/dfloat(3*natom*(3*natom+1)/2),' % Filled.'
  WRITE(6,'(I12,A)') nnzero,'  non-zero elements.'
  distave=distave/float(ll)
  drms=drms/float(ll)-distave**2.d0
  IF (drms.gt.0.d0) drms=sqrt(drms)
  WRITE(6,'(/A,F9.2,A,F9.2/(A,F9.2))') ' Average dist.  = ',distave,' +/- ',drms, &
       ' Maximum dist.  = ',dmax,' Minimum dist.  = ',dmin
  
  WRITE(6,'(/A,1PG13.6)') ' Matrix trace   = ',trace

  DEALLOCATE(a)

  WRITE(6,'(/A)')' Hessian matrix ready.'
  
  WRITE(6,'(A,1PG13.6,A)')' Potential energy:', entot/2, 'kcal mol-1'

  DO i=1,natom
     IF (occ(i).lt.1.d0) THEN
        CALL write_formatted('One or more residues have an occupancy of less than 1,', 'bright red')
        CALL write_formatted('this often means that the residues have alternative positions,', 'bright red')
        CALL write_formatted('please check your pdb and re-run if needed.', 'bright red')
        EXIT
     END IF
  END DO


END PROGRAM genENM

FUNCTION RANDOM(ISEED)
!C-----------------------------------------------------------------------
!C     RANDOM NUMBER GENERATOR: UNIFORM DISTRIBUTION (0,1)
!C     ISEED: SEED FOR GENERATOR. ON THE FIRST CALL THIS HAS TO
!C     HAVE A VALUE IN THE EXCLUSIVE RANGE (1, 2147483647)
!C     AND WILL BE REPLACED BY A NEW VALUE TO BE USED IN
!C     FOLLOWING CALL.
!C
!C     REF: Lewis, P.A.W., Goodman, A.S. & Miller, J.M. (1969)
!C     "Pseudo-random number generator for the System/360", IBM
!C     Systems Journal 8, 136.
!C
!C     This is a "high-quality" machine independent generator.
!C     INTEGERS are supposed to be 32 bits or more.
!C     The same algorithm is used as the basic IMSL generator.
!C
!C     Author: Lennart Nilsson
!C
  USE nrutil
  IMPLICIT NONE
  INTEGER :: ISEED
  REAL(DP) :: DSEED,DIVIS,DENOM,MULTIP
  REAL(DP) :: random
  DATA DIVIS/2147483647.D0/
  DATA DENOM /2147483711.D0/
  DATA MULTIP/16807.D0/
!C
  IF(ISEED.LE.1) ISEED=314159
  DSEED=MULTIP*ISEED
  DSEED=MOD(DSEED,DIVIS)
  RANDOM=DSEED/DENOM
  ISEED=DSEED
  !C
END FUNCTION RANDOM

SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
  WRITE(iunit,'(A)')"                           G  e  n  E  N  M  M                           "
  WRITE(iunit,'(A)')"                               VERSION 1.0                               "
  WRITE(iunit,'(A)')"                                                                         "  
  WRITE(iunit,'(A)')"                               Written by:                               "
  WRITE(iunit,'(A)')"                      Tom Rodgers and David Burnell                      "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program produces the Hermitian matrix for a ENM based on the pdb    "
  WRITE(iunit,'(A)')"inputed with -pdb. A matrix.sdijf file is produced which can be          "
  WRITE(iunit,'(A)')"diagonalised using DIAGSTD to calculate the eigenvalues and vectors.     "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program also outputs a map of the spring connections, ENM.vmd which "
  WRITE(iunit,'(A)')"can be viewed using VMD, and a pdb of the Ca atoms only, CAonly.pdb      "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"     Usage:                                                              "
  WRITE(iunit,'(A)')"             genENM -pdb pdbfile [-c cut-off] [-f force] [-ccust cfile]  "
  WRITE(iunit,'(A)')"                    [-fcust ffile] [-spcust spfile] [-mass] [-ca]        "
  WRITE(iunit,'(A)')"                    [-lig1] [-res] [-het] [-hin] [-hine h] [-an p]       " 
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"Option    Type       Value       Description                             "
  WRITE(iunit,'(A)')"------------------------------------------------------------             "
  WRITE(iunit,'(A)')" -pdb     Input                  pdb file name                           "
  WRITE(iunit,'(A)')"  -c      Input,Opt  12          Cut-off for spring connectivity         "
  WRITE(iunit,'(A)')"  -f      Input,Opt  1           Spring constant                         "
  WRITE(iunit,'(A)')"-ccust    Input,Opt              Allows custom cutoffs read from cfile   "
  WRITE(iunit,'(A)')"-fcust    Input,Opt              Allows custom spring constants read     "
  WRITE(iunit,'(A)')"                                 from ffile                              "
  WRITE(iunit,'(A)')"-spcust   Input,Opt              Allows custom springs to be included    "
  WRITE(iunit,'(A)')"-mass      Opt                   Uses actual atom masses                 "
  WRITE(iunit,'(A)')" -ca       Opt                   Calculates Ca model only                "
  WRITE(iunit,'(A)')" -lig1     Opt                   Assigns each residue in the HETATMs one "
  WRITE(iunit,'(A)')"                                 average point, only works with -het and "
  WRITE(iunit,'(A)')"                                 -ca flags as well                       "
  WRITE(iunit,'(A)')" -res      Opt                   Assigns ca atom the whole residue mass, "
  WRITE(iunit,'(A)')"                                 only works with -ca flag as well        "
  WRITE(iunit,'(A)')" -het      Opt                   Includes reading HETATM records         "
  WRITE(iunit,'(A)')" -hin      Opt                   Uses Hinsen fitted interactions         "
  WRITE(iunit,'(A)')" -hine     Input,Opt             Uses Hinsen exponential interactions    "
  WRITE(iunit,'(A)')"                                 instead, h is the decay factor          "
  WRITE(iunit,'(A)')"  -an      Input,Opt             Uses anisotropic interactions instead,  "
  WRITE(iunit,'(A)')"                                 p is the order of the power decay       "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')" cfile format:                                                           "
  WRITE(iunit,'(A)')"       WRITE(cfile,'(A4,1X,F7.3)') atomname, cutoff                      "
  WRITE(iunit,'(A)')"   atomname is in pdb format                                             "  
  WRITE(iunit,'(A)')"   e.g.                                                                  "  
  WRITE(iunit,'(A)')"   | C     8.000                                                         "
  WRITE(iunit,'(A)')"   | N     8.000                                                         "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')" ffile format:                                                           "
  WRITE(iunit,'(A)')"       WRITE(ffile,'(1X,I4,1X,A1,1X,I4,1X,A1,1X,F8.3)')                  "
  WRITE(iunit,'(A)')"                                         res1, chain1, res2, chain2, kij "
  WRITE(iunit,'(A)')"   res is the residue numbers, chain is the chain code.                  "
  WRITE(iunit,'(A)')"         res2 can be * which means all connections with res1             "
  WRITE(iunit,'(A)')"   e.g.                                                                  "  
  WRITE(iunit,'(A)')"   |   34 A   24 A      3.0                                              "
  WRITE(iunit,'(A)')"   |  329 A    * B      3.0                                              "
  WRITE(iunit,'(A)')"   N.B. this will only set connections if within cutoff                  "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')" spfile format:                                                          "
  WRITE(iunit,'(A)')"       WRITE(spfile,'(2(1X,I4,1X,A1))') res1, chain1, res2, chain2       "
  WRITE(iunit,'(A)')"   res is the residue numbers. chain is the chain associated with res    "
  WRITE(iunit,'(A)')"   e.g.                                                                  "  
  WRITE(iunit,'(A)')"   |   34 A   24 B                                                       "
  WRITE(iunit,'(A)')"   |  329 C  124 D                                                       "
  WRITE(iunit,'(A)')"   N.B. this will only set connections regardless of cutoff              "
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


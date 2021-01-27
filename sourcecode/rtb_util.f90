MODULE rtb_util

CONTAINS
  
  SUBROUTINE BLOCPDB ( fpdb, sstr, nbb, resi, &
       tabb, chain, natbmax, natom, nb, hetatom)
    USE nrutil
    USE read_files
    USE utils
    IMPLICIT NONE
    
    CHARACTER :: fpdb*(*), sstr*(*)
    INTEGER :: natbmax, natom, nb, nbb
    INTEGER :: numbloc(natom), resi(natom)
    INTEGER,ALLOCATABLE,DIMENSION(:) :: tabb
    REAL(DP),PARAMETER :: unknown=9999.999999 
    CHARACTER :: chain(natom)*1
    LOGICAL :: hetatom
    
    CHARACTER :: aa1*4, aa2*4, atname*4, cformat*12, cstatus*12,  &
         lign80*80, namfil*20, ssu1*1, ssu2*1, ssuch*2, li*6, elem(natom)*2
    LOGICAL ::   qerror, qexist, qinterr, qok
    INTEGER ::   at, i, ii, j, k, lnompdb, natbloc, natbmin, ndat,  &
         noff, nres, nresb, nssu, num1, num2, nunit,  &
         unrtb, io
    REAL(DP) :: tot
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: she, hel, dom
    INTEGER,ALLOCATABLE,DIMENSION(:) :: domb
    INTEGER :: nhel,nshe,ndom
    CHARACTER,ALLOCATABLE,DIMENSION(:) :: cshe*1, chel*1, cdom*1
    
    nb=0
    DO i=1,natom
       numbloc(i)=0
    END DO
    
    ! Assign based on the secondary structure in pdbfile  
    IF (sstr.eq.'SECO') THEN
       WRITE(6,'(/A)')'Substructuring:'
       WRITE(6,'(A)')'HELIX or SHEET keyword expected in coordinate file.' 
       CALL get_second(fpdb,nhel,hel,chel,nshe,she,cshe)
       
       ! loop over the HELIXs
       DO k=1,nhel
          nb=nb+1
          natbloc=0
          DO i=1,natom
             IF ((chain(i).eq.'*'.or.chel(k).eq.chain(i)).and. &
                  resi(i).ge.hel(k,1).and.resi(i).le.hel(k,2)) THEN
                IF (numbloc(i).gt.0) THEN
                   WRITE(6,'(A,I6,A)') 'Atom ',i, &
                        ' is already assigned block ',numbloc(i)
                END IF
                natbloc=natbloc+1
                numbloc(i)=nb
             END IF
          END DO
          IF (natbloc.eq.0) THEN
             WRITE(6,'(A)') 'Substructuring information ignored.'
             nb=nb-1
          END IF
       END DO
       ! loop over the SHEETs
       DO k=1,nshe
          nb=nb+1
          natbloc=0
          DO i=1,natom
             IF ((chain(i).eq.'*'.or.cshe(k).eq.chain(i)).and. &
                  resi(i).ge.she(k,1).and.resi(i).le.she(k,2)) THEN
                IF (numbloc(i).gt.0) THEN
                   WRITE(6,'(A,I6,A)') 'Atom ',i, &
                        ' is already assigned block ',numbloc(i)
                END IF
                natbloc=natbloc+1
                numbloc(i)=nb
             END IF
          END DO
          IF (natbloc.eq.0) THEN
             WRITE(6,'(A)') 'Substructuring information ignored.'
             nb=nb-1
          END IF
       END DO
       ! Split the loops into blocks:      
       IF (nb.le.0) THEN
          WRITE(6,'(/A)') 'Substructuring information not found. NONE assumed.'
       ELSE
          WRITE(6,'(I6,A)') nb, ' secondary structure elements, each in its own block.'   
          WRITE(6,'(I6,A)') nbb, ' residue(s) per block, otherwise.'
       END IF
       ! Check for any groups not included     
       nresb=1
       DO i=1,natom
          IF (numbloc(i).eq.0) THEN
             IF (i.gt.1.and.resi(i).ne.resi(i-1)) THEN
                nres=nres+1
                nresb=nresb+1
             END IF
             
             qok=.true.
             IF (i.gt.1)  &
                  qok=(resi(i).eq.resi(i-1).and.chain(i).eq.chain(i-1)).or. &
                  (resi(i)-resi(i-1).eq.1.and.chain(i).eq.chain(i-1))      
             
             IF (nresb.gt.nbb.or..not.qok) THEN
                IF (natbloc.eq.0) THEN
                   WRITE(6,'(A)') 'Substructuring information ignored.'
                   nb=nb-1
                END IF
                natbloc=1
                nresb=1
                nb=nb+1
             ELSE
                natbloc=natbloc+1
             END IF
             numbloc(i)=nb
          END IF
       END DO
       
       ! Assign based on custom domains
    ELSE IF (sstr.eq.'DOMN') THEN
       
       WRITE(6,'(/A)')'Substructuring:'
       WRITE(6,'(A)')'DOMAINs expected in input.domain file.' 
       ! Find number of DOMAINs
       INQUIRE(file='input.domain',exist=qexist)  
       IF (.not.qexist) THEN
          STOP 'An input.domain file is needed to divide into domains'
       END IF
       CALL get_domain(ndom,dom,cdom,domb)
       
       ! loop over the DOMAINs
       DO k=1,ndom
          IF (k.eq.1) THEN
             nb=nb+1
             natbloc=0
          END IF
          
          IF ((k.gt.1).and.(domb(k).ne.domb(k-1))) THEN
             nb=nb+1
             natbloc=0
          END IF
          
          DO i=1,natom
             IF ((chain(i).eq.'*'.or.cdom(k).eq.chain(i)).and. &
                  resi(i).ge.dom(k,1).and.resi(i).le.dom(k,2)) THEN
                IF (numbloc(i).gt.0) THEN
                   WRITE(6,'(A,I6,A)') 'Atom ',i, &
                        ' is already assigned block ',numbloc(i)
                END IF
                natbloc=natbloc+1
                numbloc(i)=nb
             END IF
          END DO
          IF (natbloc.eq.0) THEN
             WRITE(6,'(A)') 'Substructuring information ignored.'
             nb=nb-1
          END IF
       END DO
       
       ! Split the loops into blocks:      
       IF (nb.le.0) THEN
          WRITE(6,'(/A)') 'Substructuring information not found. NONE assumed.'
       ELSE
          WRITE(6,'(I6,A)') nb, ' secondary structure elements, each in its own block.'   
          WRITE(6,'(I6,A)') nbb, ' residue(s) per block, otherwise.'
       END IF
       ! Check for any groups not included     
       nresb=1
       DO i=1,natom
          IF (numbloc(i).eq.0) THEN
             IF (i.gt.1.and.resi(i).ne.resi(i-1)) THEN
                nres=nres+1
                nresb=nresb+1
             END IF
             
             qok=.true.
             IF (i.gt.1)  &
                  qok=(resi(i).eq.resi(i-1).and.chain(i).eq.chain(i-1)).or. &
                  (resi(i)-resi(i-1).eq.1.and.chain(i).eq.chain(i-1))      
             
             IF (nresb.gt.nbb.or..not.qok) THEN
                IF (natbloc.eq.0) THEN
                   WRITE(6,'(A)') 'Substructuring information ignored.'
                   nb=nb-1
                END IF
                natbloc=1
                nresb=1
                nb=nb+1
             ELSE
                natbloc=natbloc+1
             END IF
             numbloc(i)=nb
          END IF
       END DO
       
       ! Or assign based on the number of residues
       ! (Residues should be written in order one after the other)
    ELSE
       
       WRITE(6,'(I6,A)') nbb,' residue(s) per block.'
       nb=1
       nres=1
       nresb=1
       natbloc=0
       
       DO i=1,natom
          IF (i.gt.1.and.resi(i).ne.resi(i-1)) THEN
             nres=nres+1
             nresb=nresb+1
          END IF
          
          qok=.true.
          IF (i.gt.1)  &
               qok=(resi(i).eq.resi(i-1).and.chain(i).eq.chain(i-1)).or. &
               (resi(i)-resi(i-1).eq.1.and.chain(i).eq.chain(i-1))
          
          IF (nresb.gt.nbb.or..not.qok) THEN
             IF (natbloc.eq.0) THEN
                WRITE(6,'(A)') 'Substructuring information ignored.'
                nb=nb-1
             END IF
             natbloc=1
             nresb=1
             nb=nb+1
          ELSE
             natbloc=natbloc+1
          END IF
          numbloc(i)=nb
       END DO
       IF (natbloc.eq.0) THEN
          WRITE(6,'(A)') 'Substructuring information ignored.'
          nb=nb-1
       END IF
       
       WRITE(6,'(I6,A)') nres,' residues.'
    END IF
    
    ! Length of each block:    
    
    ALLOCATE(tabb(nb))
    k=1
    natbloc=1

    DO i=2,natom
       IF (i.eq.natom.or.numbloc(i).ne.numbloc(i-1)) THEN
          IF (i.eq.natom) natbloc=natbloc+1
          
          IF (natbloc.lt.3) THEN
             WRITE(6,'(/I2,A,I4,A)') natbloc,' atoms in block ',k, &
                  ' i.e., less than what is required for a rigid body.'
             
             IF (i.eq.natom) THEN
                WRITE(6,'(A,I6,2A,I6)') 'Last atom in this block is the ',i, &
                     'th, in residue ',chain(i),resi(i)
                WRITE(6,'(A)')'It is merged with the previous one.'
                
                k=k-1
                TABB(k)=TABB(k)+3*natbloc
                
                IF (k.eq.1.or.TABB(k).lt.natbmin) natbmin=TABB(k)
                IF (k.eq.1.or.TABB(k).gt.natbmax) natbmax=TABB(k) 
                
                WRITE(6,'(I6,A,I6)') TABB(k)/3,' atoms in block ',k 
                
                k=k+1
                natbloc=1
             ELSE
                WRITE(6,'(A,I6,2A,I6)') 'Last atom in this block is the ',i-1, &
                     'th, in residue ',chain(i-1),resi(i-1)
                WRITE(6,'(A)') 'It will be merged with next block.'
                natbloc=natbloc+1
             END IF
          ELSE
             TABB(k)=3*natbloc
             
             IF (k.eq.1.or.TABB(k).lt.natbmin) natbmin=TABB(k)
             IF (k.eq.1.or.TABB(k).gt.natbmax) natbmax=TABB(k) 
             
             IF (i.ne.natom) THEN
                k=k+1
                natbloc=1
             END IF
          END IF
       ELSE
          natbloc=natbloc+1
       END IF
    END DO
    nb=k
    natbmax=natbmax/3
    natbmin=natbmin/3
    
    WRITE(6,'(/I6,A)') nb,' blocks.'
    
    WRITE(6,'(/A,I6,A)')'At most, ',natbmax,' atoms in each of them.'
    WRITE(6,'(A,I6,A)')'At least,',natbmin,' atoms in each of them.'
    
  END SUBROUTINE BLOCPDB
  
END MODULE rtb_util

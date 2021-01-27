MODULE read_files
! read_pdb - reads pdb files
! read_second - reads the secondary structures from pdb for copy
! get_second - reads the secondary structures from a pdbfile
! get_domain - reads an input.domain file
! read_amb_prm - reads amber parameter files
! read_amb_trj - reads amber trajectory files
! read_eigenfacs - reads eigenfacs format files
! read_gro - reads gromacs .gro file
! read_trr_out - read gromacs gmxdump version of .trr file
! read_CONFIG - read dl_poly CONFIG file
! read_HISTORY - read dl_poly HISTORY file
!
CONTAINS
! Read pdb files
  SUBROUTINE read_pdb(pdbfile,hetatm,natom,atom,atnum,name,res,chain,resnum,x,y,z, &
                  occ,bfac,elem,chag)
    USE nrtype
    IMPLICIT NONE
    INTEGER, ALLOCATABLE,DIMENSION(:) :: atnum,resnum
    REAL(DP), ALLOCATABLE,DIMENSION(:) ::x,y,z,occ,bfac
    CHARACTER, ALLOCATABLE, DIMENSION(:) :: atom*6,name*4,res*3,chain*1,elem*2,chag*2
    INTEGER :: io, i, natom
    CHARACTER :: pdbfile*120, lign80*80
    LOGICAL :: hetatm
    
    ! Find number of atoms
    OPEN(file=pdbfile,form="FORMATTED",status="OLD",unit=2356)
    i=0
    DO
       READ(2356,'(A)',IOSTAT=io) lign80
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
          ! Count number of atoms
       ELSE
          IF ((lign80(1:4).eq.'ATOM').or. &
               ((lign80(1:6).eq.'HETATM').and.hetatm)) THEN
             i=i+1
          END IF
       END IF
    END DO
    
    natom=i
    
    REWIND(2356)
    
    ALLOCATE(x(natom))
    ALLOCATE(y(natom))
    ALLOCATE(z(natom))
    ALLOCATE(atom(natom))
    ALLOCATE(atnum(natom))
    ALLOCATE(name(natom))
    ALLOCATE(res(natom))
    ALLOCATE(chain(natom))
    ALLOCATE(resnum(natom))
    ALLOCATE(occ(natom))
    ALLOCATE(bfac(natom))
    ALLOCATE(elem(natom))
    ALLOCATE(chag(natom))  
    
    ! Read in atom information
    i=0
    DO
       READ(2356,'(A)',IOSTAT=io) lign80
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
          ! Count number of atoms
       ELSE
          IF ((lign80(1:4).eq.'ATOM').or. &
               ((lign80(1:6).eq.'HETATM').and.hetatm)) THEN
             i=i+1
             READ(lign80,'(A6,I5,1X,A4,1X,A3,1X,A1,I4,1X,3X,3F8.3,2F6.2,10X,2A2)') &
                  atom(i),atnum(i),name(i),res(i),chain(i),resnum(i),x(i),y(i),z(i),&
                  occ(i),bfac(i),elem(i),chag(i)
          END IF
       END IF
    END DO
    
    CLOSE(2356)
    
  END SUBROUTINE read_pdb

! Gets the secondard structure from pdbfile
  SUBROUTINE read_second(pdbfile,nshe,she)
    IMPLICIT NONE
    INTEGER :: io, j, nshe
    CHARACTER :: pdbfile*120, lign80*80
    CHARACTER,ALLOCATABLE,DIMENSION(:) :: she*80

    ! Find number of atoms
    OPEN(file=pdbfile,form="FORMATTED",status="OLD",unit=2356)
    nshe=0
    DO
       READ(2356,'(A)',IOSTAT=io) lign80
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
          ! Count number of atoms
       ELSE
          IF ((lign80(1:5).eq.'HELIX').or.(lign80(1:5).eq.'SHEET')) THEN
             nshe=nshe+1
          END IF
       END IF
    END DO
    
    REWIND(2356)

    ALLOCATE(she(nshe))
    
    ! Read in atom information
    j=0
    DO
       READ(2356,'(A)',IOSTAT=io) lign80
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
          ! Count number of atoms
       ELSE
          IF ((lign80(1:5).eq.'HELIX').or.(lign80(1:5).eq.'SHEET')) THEN
             j=j+1
             READ(lign80,'(A80)') she(j)
          END IF
       END IF
    END DO

    CLOSE(2356)

  END SUBROUTINE read_second

! Gets the secondard structure from pdbfile
  SUBROUTINE get_second(pdbfile,nhel,hel,chel,nshe,she,cshe)
    IMPLICIT NONE
    INTEGER :: io, i, j, nhel,nshe,serNum,class,len,strand,f
    CHARACTER :: pdbfile*120, lign80*80, name*6, helixID*3, ResName*3, chain*1, init*1,comm*29,sheetID*3
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: she, hel, b, g
    CHARACTER,ALLOCATABLE,DIMENSION(:) :: cshe*1, chel*1, v*1

    ! Find number of atoms
    OPEN(file=pdbfile,form="FORMATTED",status="OLD",unit=2356)
    nhel=0
    nshe=0
    DO
       READ(2356,'(A)',IOSTAT=io) lign80
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
          ! Count number of atoms
       ELSE
          IF (lign80(1:5).eq.'HELIX') THEN
             nhel=nhel+1
          ELSEIF (lign80(1:5).eq.'SHEET') THEN
             nshe=nshe+1
          END IF
       END IF
    END DO
    
    REWIND(2356)

    ALLOCATE(hel(nhel,2))
    ALLOCATE(she(nshe,2))
    ALLOCATE(cshe(nshe))
    ALLOCATE(chel(nhel))
    

    ! Read in atom information
    i=0
    j=0
    DO
       READ(2356,'(A)',IOSTAT=io) lign80
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
          ! Count number of atoms
       ELSE
          IF (lign80(1:5).eq.'HELIX') THEN
             i=i+1
             READ(lign80,'(A6,1X,I3,1X,A3,1X,A3,1X,A1,1X,I4,A1,1X,A3,1X,A1,1X,I4,A1,I2,A29,I4)') &
                  name,serNum,helixID,ResName,chel(i),hel(i,1),init,ResName,chain,hel(i,2),init,class,comm,len
          ELSEIF (lign80(1:5).eq.'SHEET') THEN
             j=j+1
             READ(lign80,'(A6,1X,I3,1X,A3,I2,1X,A3,1X,A1,I4,A1,1X,A3,1X,A1,I4)') &
                  name,strand,sheetID,len,ResName,cshe(j),she(j,1),init,ResName,chain,she(j,2)
          END IF
       END IF
    END DO
    
    CLOSE(2356)
    
    ! Sort the structures as SHEET groups are not in order
    
    ALLOCATE(b(nhel,2))
    ALLOCATE(g(nhel,1))
    ALLOCATE(v(nhel))
    b=hel
    v=chel
    DO j=1,nhel
       f=MINVAL(hel(:,1))
       
       DO i=1,nhel
          IF (f.eq.hel(i,1)) THEN
             g(j,1)=i
             hel(i,1)=9000
             EXIT
          END IF
       END DO
    END DO

    DO i=1,nhel
       hel(i,:)=b(g(i,1),:)
       chel(i)=v(g(i,1))
    END DO

    DEALLOCATE(b)
    DEALLOCATE(g)
    DEALLOCATE(v)
    ALLOCATE(b(nshe,2))
    ALLOCATE(g(nshe,1))
    ALLOCATE(v(nshe))
    b=she
    v=cshe
    DO j=1,nshe
       f=MINVAL(she(:,1))
       
       DO i=1,nshe
          IF (f.eq.she(i,1)) THEN
             g(j,1)=i
             she(i,1)=9000
             EXIT
          END IF
       END DO
    END DO

    DO i=1,nshe
       she(i,:)=b(g(i,1),:)
       cshe(i)=v(g(i,1))
    END DO

    DEALLOCATE(b)
    DEALLOCATE(g)
    DEALLOCATE(v)


  END SUBROUTINE get_second

  SUBROUTINE get_domain(ndom,dom,cdom,domb)
    IMPLICIT NONE
    INTEGER :: io, i, j, ndom, f
    CHARACTER :: lign80*80
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: dom, b
    INTEGER,ALLOCATABLE,DIMENSION(:) :: g, domb, w
    CHARACTER,ALLOCATABLE,DIMENSION(:) :: cdom*1, v*1

OPEN(file='input.domain',form="FORMATTED",status="OLD",unit=2356)
    ndom=0
    DO
       READ(2356,'(A)',IOSTAT=io) lign80
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
          ! Count number of atoms
       ELSE
          ndom=ndom+1
       END IF
    END DO
    
    REWIND(2356)

    ALLOCATE(dom(ndom,2))
    ALLOCATE(cdom(ndom))
    ALLOCATE(domb(ndom))

    ! Read in DOMAIN information
    i=0
    DO
       READ(2356,'(A)',IOSTAT=io) lign80
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
          ! Count number of atoms
       ELSE
             i=i+1
             READ(lign80,'(I4,1X,A1,1X,I4,1X,I4)')domb(i),cdom(i),dom(i,1),dom(i,2)
       END IF
    END DO
    
    CLOSE(2356)

    ALLOCATE(b(ndom,2))
    ALLOCATE(g(ndom))
    ALLOCATE(v(ndom))
    ALLOCATE(w(ndom))
    b=dom
    v=cdom
    w=domb
    
    DO j=1,ndom
       f=MINVAL(domb)
       
       DO i=1,ndom
          IF (f.eq.domb(i)) THEN
             g(j)=i
             domb(i)=9000
             EXIT
          END IF
       END DO
    END DO

    DO i=1,ndom
       dom(i,:)=b(g(i),:)
       cdom(i)=v(g(i))
       domb(i)=w(g(i))
    END DO

    DEALLOCATE(b)
    DEALLOCATE(g)
    DEALLOCATE(v)
    DEALLOCATE(w)


  END SUBROUTINE get_domain

! Read the amber paramaters file
  SUBROUTINE read_amb_par(prmfile,natomq,natom,atom,atnum,nameprm,nameres, &
       chain,resnumamb,occ,bfac,mass,elem,chag,natomsol,sol)
    USE nrtype
    IMPLICIT NONE
    CHARACTER :: prmfile*120, lign80*80
    INTEGER :: natomq, natom,nres,i,io, j,len,natomsol
    CHARACTER,ALLOCATABLE,DIMENSION(:) :: nameprm*4, nameprmtemp*4,nameres*3, namerestemp*4,atom*6, &
         chain*1,elem*2,chag*2
    CHARACTER :: name*4
    REAL(DP),ALLOCATABLE,DIMENSION(:) :: mass,occ,bfac
    INTEGER,ALLOCATABLE,DIMENSION(:) :: respoint,resnumamb,atnum
    LOGICAL :: sol

! prmfile  - The amber parameter filename
! natomq   - Full number of atoms
! natom - Number of atoms that are not SOL (i.e. the ones we want)
! nres     - Number of residues
! nameprm  - Atom names
! nameres  - Residue names
! mass     - Atom masses

    OPEN(file=prmfile,form='FORMATTED',status='OLD',unit=5349)
    
    ! Find the number of atoms in the trajectory
    DO
       READ(5349,'(A)',IOSTAT=io) lign80
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
          ! Count number of atoms
       ELSE
          IF (lign80(2:14).eq."FLAG POINTERS") THEN
             READ(5349,*)
             READ(5349,'(I8)') natomq
             READ(5349,'(2I8)')i,nres
             EXIT
          END IF
       END IF
    END DO
    
    ! Read the atom names
    ALLOCATE(nameprmtemp(natomq))
    DO
       READ(5349,'(A)',IOSTAT=io) lign80
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
          ! Count number of atoms
       ELSE
          IF (lign80(2:15).eq."FLAG ATOM_NAME") THEN
             READ(5349,*)
             DO
                READ(5349,'(20A4)',IOSTAT=io) (nameprmtemp(i), i=1,natomq)
                IF (io > 0) THEN
                   WRITE(*,*) 'Check input.  Something was wrong'
                   STOP
                   ! End of file
                ELSE IF (io < 0) THEN
                   EXIT
                ELSE
                   EXIT
                END IF
             END DO
             EXIT
          END IF
       END IF
    END DO
    
    ! Find the number of atoms before the SOL as we are not interested in this
    
    DO i=1,natomq
       IF (i.le.natomq-2) THEN
          IF (nameprmtemp(i).eq."O   ".and.nameprmtemp(i+1).eq."H1  ".and.nameprmtemp(i+2).eq."H2  ") THEN
             natom=i-1
             EXIT
          END IF
       ELSE
          natom=natomq
       END IF
    END DO
    natomsol=natom
    IF (sol) THEN
       natom=natomq
    END IF
    
    ALLOCATE(nameprm(natom))
    nameprm=nameprmtemp(1:natom)
    DEALLOCATE(nameprmtemp)

    ! fix atom names to the right aligned
    
    DO i=1,natom
       len=index(nameprm(i),' ')
       len=len-1
       name=nameprm(i)(1:len)
       IF (len.eq.1) THEN
          WRITE(nameprm(i),'(A,A,2A)')' ',name(1:len),'  '
       ELSE IF (len.eq.2) THEN
          WRITE(nameprm(i),'(A,2A,A)')' ',name(1:len),' '
       ELSE IF (len.eq.3) THEN
          WRITE(nameprm(i),'(A,3A)')' ',name(1:len)
       END IF
    END DO
   
    WRITE(6,'(A,I5,A)') "Trajectory contains ",natomq, " atoms,"
    WRITE(6,'(A,I5,A)') "of which ",natomsol, " atoms are not SOL."
    
    ! Read in the atom masses
    ALLOCATE(mass(natom))
    DO
       READ(5349,'(A)',IOSTAT=io) lign80
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
          ! Count number of atoms
       ELSE
          IF (lign80(2:10).eq."FLAG MASS") THEN
             READ(5349,*)
             DO
                READ(5349,'(5E16.8)',IOSTAT=io) (mass(i), i=1,natom)
                IF (io > 0) THEN
                   WRITE(*,*) 'Check input.  Something was wrong'
                   STOP
                   ! End of file
                ELSE IF (io < 0) THEN
                   EXIT
                   ! Count number of atoms
                ELSE
                   EXIT
                END IF
             END DO
             EXIT
          END IF
       END IF
    END DO
    
    WRITE(6,'(A)') "Atom masses read from parameter file."

    ! Read in residue names 
    ALLOCATE(namerestemp(nres))
    DO
       READ(5349,'(A)',IOSTAT=io) lign80
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
          ! Count number of atoms
       ELSE
          IF (lign80(2:19).eq."FLAG RESIDUE_LABEL") THEN
             READ(5349,*)
             DO
                READ(5349,'(20A4)',IOSTAT=io) (namerestemp(i), i=1,nres)
                IF (io > 0) THEN
                   WRITE(*,*) 'Check input.  Something was wrong'
                   STOP
                   ! End of file
                ELSE IF (io < 0) THEN
                   EXIT
                   ! Count number of atoms
                ELSE
                   EXIT
                END IF
             END DO
             EXIT
          END IF
       END IF
    END DO
    !Strip of SOL as we are not interested in these
    DO i=1,nres
       IF (namerestemp(i).eq."WAT ") THEN
          EXIT
       END IF
    END DO
    nres=i-1
  
    WRITE(6,'(A)') "Atom residues read from parameter file."

    ! Find atoms in each residue
    ALLOCATE(respoint(nres))
    DO
       READ(5349,'(A)',IOSTAT=io) lign80
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
          ! Count number of atoms
       ELSE
          IF (lign80(2:21).eq."FLAG RESIDUE_POINTER") THEN
             READ(5349,*)
             DO
                READ(5349,'(10I8)',IOSTAT=io) (respoint(i), i=1,nres)
                IF (io > 0) THEN
                   WRITE(*,*) 'Check input.  Something was wrong'
                   STOP
                   ! End of file
                ELSE IF (io < 0) THEN
                   EXIT
                   ! Count number of atoms
                ELSE
                   EXIT
                END IF
             END DO
             EXIT
          END IF
       END IF
    END DO

    CLOSE(5349)

   WRITE(6,'(A)') "Atom residues numbers from parameter file." 

    ALLOCATE(resnumamb(natom))
    ALLOCATE(nameres(natom))
    j=1
    DO i=1,natom
       resnumamb(i)=j
       READ(namerestemp(j),'(A3)') nameres(i)
       IF (j.lt.nres) THEN
          IF (i.lt.respoint(j+1)) THEN
          ELSE
             j=j+1
             resnumamb(i)=j
             READ(namerestemp(j),'(A3)') nameres(i)
          END IF
       END IF
    END DO
    
    DEALLOCATE(namerestemp)
    DEALLOCATE(respoint)

    ALLOCATE(atom(natom))
    ALLOCATE(atnum(natom))
    ALLOCATE(chain(natom))
    ALLOCATE(elem(natom))
    ALLOCATE(chag(natom))    
    ALLOCATE(occ(natom))
    ALLOCATE(bfac(natom))

    DO i=1,natom
       atom(i)="ATOM  "
       atnum(i)=i
       chain(i)="A"
       chag(i)="  "
       occ(i)=1.d0
       bfac(i)=0.d0
       IF ((index(nameprm(i),'Cl').gt.0).or.(index(nameprm(i),'CL').gt.0)) THEN
          elem(i)="Cl"   
       ELSE IF (index(nameprm(i),'C').gt.0) THEN
          elem(i)=" C"
       ELSE IF ((index(nameprm(i),'NA').gt.0).or.(index(nameprm(i),'Na').gt.0)) THEN
          elem(i)="Na"   
       ELSE IF (index(nameprm(i),'N').gt.0) THEN
          elem(i)=" N"
       ELSE IF (index(nameprm(i),'P').gt.0) THEN
          elem(i)=" P"
       ELSE IF (index(nameprm(i),'S').gt.0) THEN
          elem(i)=" S"
       ELSE IF (index(nameprm(i),'O').gt.0) THEN
          elem(i)=" O"   
       ELSE IF (index(nameprm(i),'H').gt.0) THEN
          elem(i)=" H"
       ELSE IF ((index(nameprm(i),'Fe').gt.0).or.(index(nameprm(i),'FE').gt.0)) THEN
          elem(i)="Fe"   
       ELSE 
          WRITE(*,*) 'name', nameprm(i), 'not found. Please add.'
       END IF
    END DO

  END SUBROUTINE read_amb_par
  
!-----------------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------------  
  SUBROUTINE read_amb_trj(trajfile,natomq,natom,x,xav,start,end,nframes,box)
    USE nrutil
    IMPLICIT NONE
    CHARACTER :: trajfile*120
    INTEGER :: natomq, natom, i, j, k, start, end, io, nframes
    REAL(DP) :: box1,box2,box3
    REAL(DP),ALLOCATABLE,DIMENSION(:) :: xav,xi,yi,zi
    REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: x,box

    ALLOCATE(xi(natomq))
    ALLOCATE(yi(natomq))
    ALLOCATE(zi(natomq))
  
    ALLOCATE(xav(natom*3))
    j=0
    k=0
    xav=0
    
    OPEN(file=trajfile,form="FORMATTED",status='OLD',unit=6999)
    READ(6999,*)
    DO
       READ(6999,'(10F8.3)',IOSTAT=io) (xi(i), yi(i), zi(i), i=1,natomq)
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
          ! Count number of atoms
       ELSE
          READ(6999,*)
          j=j+1
          IF ((j.ge.start).and.(j.le.end)) THEN 
             k=k+1
             DO i=1,natom
                xav(3*(i-1)+1)=xav(3*(i-1)+1)+xi(i)
                xav(3*(i-1)+2)=xav(3*(i-1)+2)+yi(i)
                xav(3*(i-1)+3)=xav(3*(i-1)+3)+zi(i)
             END DO
          END IF
       END IF
    END DO
    nframes=k
    
    xav=xav/nframes
    
    WRITE(6,'(A,I5,A)') "Found ", j, " frames in the trajectory,"
    WRITE(6,'(A,I5,A)') "using ", nframes, " frames."
    WRITE(6,'(A)') "Calculated average coordinates"  
    
    WRITE(6,'(A)') "Reading trajectory frames"
    
    ALLOCATE(x(3*natom,nframes))
    ALLOCATE(box(9,nframes))
    box=0

    j=0
    k=0
    WRITE(6,'(A11)',ADVANCE='NO') "Read frame "
    
    REWIND(6999)
    READ(6999,*)
    DO
       READ(6999,'(10F8.3)',IOSTAT=io) (xi(i), yi(i), zi(i), i=1,natomq)
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
          ! Count number of atoms
       ELSE
          READ(6999,*) box1, box2, box3
          j=j+1
          IF ((j.ge.start).and.(j.le.end)) THEN 
             k=k+1
             DO i=1,natom
                x(3*(i-1)+1,k)=xi(i)
                x(3*(i-1)+2,k)=yi(i)
                x(3*(i-1)+3,k)=zi(i)
             END DO
             box(1,k)=box1
             box(5,k)=box2
             box(9,k)=box3

             IF (nframes.lt.10) THEN
                WRITE(6,'(1X,I5)',ADVANCE='NO') k
             ELSE
                IF (int(k/(nframes/10))*(nframes/10).eq.k) THEN
                   WRITE(6,'(1X,I5)',ADVANCE='NO') k
                END IF
             END IF
          END IF
       END IF
    END DO

    WRITE(6,'(A)') ""
   
    DEALLOCATE(xi)
    DEALLOCATE(yi)
    DEALLOCATE(zi)
    
  END SUBROUTINE read_amb_trj
  
!-----------------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------------  
  SUBROUTINE read_amb_inpcrd(trajfile,natomq,natom,x,xav,start,end,nframes)
    USE nrutil
    IMPLICIT NONE
    CHARACTER :: trajfile*120
    INTEGER :: natomq, natom, i, j, k, start, end, io, nframes
    REAL(DP),ALLOCATABLE,DIMENSION(:) :: xav,xi,yi,zi
    REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: x

    ALLOCATE(xi(natomq))
    ALLOCATE(yi(natomq))
    ALLOCATE(zi(natomq))
  
    ALLOCATE(xav(natom*3))
    j=0
    k=0
    xav=0
    
    OPEN(file=trajfile,form="FORMATTED",status='OLD',unit=6999)
    READ(6999,*)
    READ(6999,*) natomq

    DO
       READ(6999,*,IOSTAT=io) (xi(i), yi(i), zi(i), i=1,natomq) !'(10F8.3)'
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
          ! Count number of atoms
       ELSE
          READ(6999,*)
          j=j+1
          IF ((j.ge.start).and.(j.le.end)) THEN 
             k=k+1
             DO i=1,natom
                xav(3*(i-1)+1)=xav(3*(i-1)+1)+xi(i)
                xav(3*(i-1)+2)=xav(3*(i-1)+2)+yi(i)
                xav(3*(i-1)+3)=xav(3*(i-1)+3)+zi(i)
             END DO
          END IF
       END IF
    END DO
    nframes=k
    
    xav=xav/nframes
    
    WRITE(6,'(A,I5,A)') "Found ", j, " frames in the trajectory,"
    WRITE(6,'(A,I5,A)') "using ", nframes, " frames."
    WRITE(6,'(A)') "Calculated average coordinates"  
    
    WRITE(6,'(A)') "Reading trajectory frames"
    
    ALLOCATE(x(3*natom,nframes))

    j=0
    k=0
    WRITE(6,'(A11)',ADVANCE='NO') "Read frame "
    
    REWIND(6999)
    READ(6999,*)
    DO
       READ(6999,*,IOSTAT=io) (xi(i), yi(i), zi(i), i=1,natomq) !'(10F8.3)'
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
          ! Count number of atoms
       ELSE
          READ(6999,*)
          j=j+1
          IF ((j.ge.start).and.(j.le.end)) THEN 
             k=k+1
             DO i=1,natom
                x(3*(i-1)+1,k)=xi(i)
                x(3*(i-1)+2,k)=yi(i)
                x(3*(i-1)+3,k)=zi(i)
             END DO
             
             IF (nframes.lt.10) THEN
                WRITE(6,'(1X,I5)',ADVANCE='NO') k
             ELSE
                IF (int(k/(nframes/10))*(nframes/10).eq.k) THEN
                   WRITE(6,'(1X,I5)',ADVANCE='NO') k
                END IF
             END IF
          END IF
       END IF
    END DO

    WRITE(6,'(A)') ""
   
    DEALLOCATE(xi)
    DEALLOCATE(yi)
    DEALLOCATE(zi)
    
  END SUBROUTINE read_amb_inpcrd
  
!-----------------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------------
  SUBROUTINE read_eigenfacs(filename,natom,startvec,endvec,ndim,eigenval,eigenvec,num)
    USE nrutil
    IMPLICIT NONE
    CHARACTER :: filename*120, lign80*80
    INTEGER :: natom, j, io, ndim, vector, i, k, startvec, endvec,num
    REAL(DP) :: ei
    REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: eigenvec
    REAL(DP),ALLOCATABLE,DIMENSION(:) :: eigenval

  OPEN(file=filename,form="FORMATTED",status="OLD",unit=4565)

  natom=0
  j=0
  ! Count number of atoms
  DO
     READ(4565,'(A)',IOSTAT=io) lign80
     IF (io > 0) THEN
        WRITE(*,*) 'Check input.  Something was wrong'
        STOP
        ! End of file
     ELSE IF (io < 0) THEN
        EXIT
        ! Count up number of atoms
     ELSE
        IF (lign80(1:7).eq.' VECTOR') THEN
           READ(lign80,'(7X,I5)') vector
           IF ((vector.ge.startvec).and.(vector.le.endvec)) THEN
              j=j+1   
           END IF
           natom=0
        ELSE IF (lign80(1:7).eq.' ------') THEN
        ELSE
           natom = natom + 1
        END IF
     END IF
  END DO
  
  REWIND(4565)

  ALLOCATE(eigenvec(natom*ndim,j))
  ALLOCATE(eigenval(j))

  j=0
  k=0
  DO
     
     READ(4565,'(A)',IOSTAT=io) lign80
     IF (io > 0) THEN
        WRITE(*,*) 'Check input.  Something was wrong'
        STOP
        ! End of file
     ELSE IF (io < 0) THEN
        EXIT
     ELSE
        IF (lign80(1:7).eq.' VECTOR') THEN
           READ(lign80,'(7X,I5,12X,G12.4)') vector,ei
           i=0
           IF ((vector.ge.startvec).and.(vector.le.endvec)) THEN
              j=vector-startvec+1
              eigenval(j)=ei
           END IF
        ELSE IF (lign80(1:7).eq.' ------') THEN
        ELSE 
           IF ((vector.ge.startvec).and.(vector.le.endvec)) THEN
              j=vector-startvec+1
              READ(lign80,'(3(1X,G11.4))') (eigenvec(i+k,j),k=1,ndim)
              i=i+ndim
           END IF
        END IF
     END IF
     
  END DO

  num=j

  CLOSE(4565)

  END SUBROUTINE read_eigenfacs
  

  SUBROUTINE read_gro(grofile,natomq,natom,atom,atnum,name,res, &
       chain,resnum,occ,bfac,mass,elem,chag)
    USE nrtype
    IMPLICIT NONE
    CHARACTER :: grofile*120, lign80*80
    INTEGER :: natomq, natom,nres,i,io, j
    CHARACTER,ALLOCATABLE,DIMENSION(:) :: name*4,res*3,atom*6, &
         chain*1,elem*2,chag*2
    REAL(DP),ALLOCATABLE,DIMENSION(:) :: mass,occ,bfac
    INTEGER,ALLOCATABLE,DIMENSION(:) :: respoint,resnum,atnum

    OPEN(file=grofile,form='FORMATTED',status='OLD',unit=9324)

    READ(9324,*) !skip header line

    READ(9324,'(I6)') natomq

    i=0
    DO j=1,natomq
       READ(9324,'(A)',IOSTAT=io) lign80
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
          ! Count number of atoms
       ELSE
          IF (lign80(6:8).ne.'SOL'.and.lign80(6:7).ne.'NA'.and.lign80(6:7).ne.'CL') THEN
             i=i+1
          END IF
       END IF
    END DO

    natom=i

    ALLOCATE(atom(natom))
    ALLOCATE(atnum(natom))
    ALLOCATE(name(natom))
    ALLOCATE(res(natom))
    ALLOCATE(chain(natom))
    ALLOCATE(resnum(natom))
    ALLOCATE(occ(natom))
    ALLOCATE(bfac(natom))
    ALLOCATE(mass(natom))
    ALLOCATE(elem(natom))
    ALLOCATE(chag(natom))


    REWIND(9324)
    READ(9324,*) !skip header line
    READ(9324,*) !skip number of atoms line
    
    i=0
    DO j=1,natomq
       READ(9324,'(A)',IOSTAT=io) lign80
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
          ! Count number of atoms
       ELSE
          IF (lign80(6:8).ne.'SOL'.and.lign80(6:7).ne.'NA'.and.lign80(6:7).ne.'CL') THEN
             i=i+1
             READ(lign80,'(I5,A3,3X,A4,I5)') resnum(i),res(i),name(i),atnum(i)
          ELSE
             EXIT
          END IF
       END IF
    END DO
    
    CLOSE(9324)
    
    WRITE(6,'(A,I5,A)') "Trajectory contains ",natomq, " atoms,"
    WRITE(6,'(A,I5,A)') "of which ",natom, " atoms are not SOL."
    
    DO i=1,natom
       atom(i)="ATOM  "
       chain(i)="A"
       elem(i)="  "
       chag(i)="  "
       occ(i)=1.d0
       bfac(i)=0.d0
       mass(i)=1.d0
    END DO

  END SUBROUTINE read_gro

  SUBROUTINE read_trr_out(trajfile,natomq,natom,x,xav,start,end,skip,nframes,box)
    USE nrutil
    IMPLICIT NONE
    CHARACTER :: trajfile*120, lign80*80
    INTEGER :: natomq, natom, i, j, k, start, end, io, nframes,sval,eval,skip,l,nframorg, zzz
    REAL(DP),ALLOCATABLE,DIMENSION(:) :: xav,xi,yi,zi
    REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: x,box
    
    OPEN(file=trajfile,form="FORMATTED",status="OLD",unit=4568)
    
    ! Get the number of atoms and eigenvectors
    
    i=0
    j=0
    DO
       READ(4568,'(A)',IOSTAT=io) lign80
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
       ELSE     
          IF (INDEX(lign80,'frame').gt.0) THEN
             i=i+1
             IF (i.ge.start.and.i.le.end) THEN
                j=j+1
             END IF
          ELSE IF (INDEX(lign80,'natoms').gt.0) THEN
             READ(lign80,'(11X,I9)') natomq
          END IF
       END IF
    END DO
    
    nframorg=j
    nframes=int(j/skip)

    WRITE(6,'(A,I5,A)') "Found ", i, " frames in the trajectory,"
    WRITE(6,'(A,I5,A,I5,A)') "Start frame ", start, ", End frame", min(i,end),","
    WRITE(6,'(A,I5,A)') "Number of frames in range is ", nframorg, " frames,"
    WRITE(6,'(A,I5,A)') "Taking every ", skip, " frames,"
    WRITE(6,'(A,I5,A)') "using ", nframes, " frames."
    
    WRITE(6,'(A)') "Reading trajectory frames"

    ALLOCATE(x(3*natom,nframes))
    ALLOCATE(xav(3*natom))
    ALLOCATE(box(9,nframes))
    
    REWIND(4568)
    
    ! Read in the trajectory
    j=0
    i=0
    k=0
    l=0
    WRITE(6,'(A11)',ADVANCE='NO') "Read frame "
    DO
       READ(4568,'(A)',IOSTAT=io) lign80
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
       ELSE
          IF (INDEX(lign80,'frame').gt.0) THEN
             j=j+1
             i=0
             zzz=0
             IF (j.ge.start.and.j.le.end) THEN
                k=k+1
                IF (nframes.gt.10) THEN
                   IF (int(k/(nframorg/10))*(nframorg/10).eq.k) THEN
                      WRITE(6,'(1X,I5)',ADVANCE='NO') k
                   END IF
                ELSE
                   WRITE(6,'(1X,I5)',ADVANCE='NO') k
                END IF
                IF (int(k/skip)*skip.eq.k) THEN
                   l=l+1
                END IF
             END IF
          ELSE IF (INDEX(lign80,'natoms').gt.0) THEN
          ELSE IF (INDEX(lign80,'box[').gt.0) THEN
             IF (j.ge.start.and.j.le.end) THEN
                IF (int(k/skip)*skip.eq.k) THEN
                   zzz=zzz+1
                   sval=INDEX(lign80,'{')
                   eval=INDEX(lign80,'}')
                   READ(lign80(sval+1:eval-1),*) box(3*(zzz-1)+1,l), &
                           box(3*(zzz-1)+2,l), box(3*(zzz-1)+3,l) 
                END IF
             END IF
          ELSE IF (INDEX(lign80,'x3').gt.0) THEN
          ELSE
             IF (j.ge.start.and.j.le.end) THEN
                IF (int(k/skip)*skip.eq.k) THEN
                   i=i+1
                   IF (i.le.natom) THEN
                      sval=INDEX(lign80,'{')
                      eval=INDEX(lign80,'}')
                      READ(lign80(sval+1:eval-1),*) x(3*(i-1)+1,l), &
                           x(3*(i-1)+2,l), x(3*(i-1)+3,l)                  
                   END IF
                END IF
             END IF
          END IF
       END IF
    END DO
    WRITE(6,'(A)') ""
    CLOSE(4568)
    
    ! Convert from nm to A so it matches standard format
    x=x*10.d0
    box=box*10.d0

    xav=0
    DO i=1,nframes
       xav=xav+x(:,i)
    END DO
    
    xav=xav/real(nframes)

    WRITE(6,'(A)') "Calculated average coordinates"
    
    
  END SUBROUTINE read_trr_out

  SUBROUTINE read_CONFIG(dlfile,natomq,natom,atom,atnum,name,res, &
       chain,resnum,occ,bfac,mass,elem,chag)
    USE nrtype
    IMPLICIT NONE
    CHARACTER :: dlfile*120, lign80*80,temp*9
    CHARACTER*1 :: alp(26)
    INTEGER :: natomq, natom,nres,i,io, j, levcfg, imcon,sval,val
    CHARACTER,ALLOCATABLE,DIMENSION(:) :: name*4,res*3,atom*6, &
         chain*1,elem*2,chag*2
    REAL(DP),ALLOCATABLE,DIMENSION(:) :: mass,occ,bfac
    INTEGER,ALLOCATABLE,DIMENSION(:) :: respoint,resnum,atnum
    
    alp(1)='A';alp(2)='B';alp(3)='C';alp(4)='D';alp(5)='E';alp(6)='F'
    alp(7)='G';alp(8)='H';alp(9)='I';alp(10)='J';alp(11)='K';alp(12)='L'
    alp(13)='M';alp(14)='N';alp(15)='O';alp(16)='P';alp(17)='Q';alp(18)='R'
    alp(19)='S';alp(20)='T';alp(21)='U';alp(22)='V';alp(23)='W';alp(24)='X'
    alp(25)='Y';alp(26)='Z'

    OPEN(file=dlfile,form='FORMATTED',status='OLD',unit=6524)
    
    READ(6524,'(A)') lign80 ! Header line
    READ(6524,*) levcfg, imcon ! format info

    IF (imcon.gt.0) THEN ! if periodic bondaries skip over these lines
       READ(6524,'(A)') lign80
       READ(6524,'(A)') lign80
       READ(6524,'(A)') lign80
    END IF
    
    i=0
    j=0
    
    DO
       READ(6524,'(A)',IOSTAT=io) lign80
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
       ELSE  
          IF (lign80(1:2).ne."OW".and.lign80(1:2).ne."HW".and. &
              lign80(1:2).ne."Na") THEN
             i=i+1
          END IF
          j=j+1
          READ(6524,'(A)') lign80 ! skip coordinates
          IF (levcfg.gt.0) THEN
             READ(6524,'(A)') lign80 ! skip velocities
          END IF
          IF (levcfg.gt.1) THEN
             READ(6524,'(A)') lign80 ! skip forces
          END IF
       END IF
    END DO

    natom=i
    natomq=j

    REWIND(6524)

    ALLOCATE(atom(natom))
    ALLOCATE(atnum(natom))
    ALLOCATE(name(natom))
    ALLOCATE(res(natom))
    ALLOCATE(chain(natom))
    ALLOCATE(resnum(natom))
    ALLOCATE(occ(natom))
    ALLOCATE(bfac(natom))
    ALLOCATE(mass(natom))
    ALLOCATE(elem(natom))
    ALLOCATE(chag(natom))

    READ(6524,'(A)') lign80 ! Header line
    READ(6524,'(A)') lign80 ! format info

    IF (imcon.gt.0) THEN ! if periodic bondaries skip over these lines
       READ(6524,'(A)') lign80
       READ(6524,'(A)') lign80
       READ(6524,'(A)') lign80
    END IF

    i=0
    
    DO
       READ(6524,'(A)',IOSTAT=io) lign80
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
       ELSE     
          IF (lign80(1:2).ne."OW".and.lign80(1:2).ne."HW".and. &
              lign80(1:2).ne."Na") THEN
             i=i+1
             READ(lign80,'(A4,10X,I7,5X,A9)') name(i),atnum(i),temp

             sval=10
             DO j=1,26
                val=INDEX(temp,alp(j))
                IF (val.gt.0) THEN
                   sval=MIN(sval,val)
                END IF
             END DO
             res(i)=temp(sval:sval+3)

             IF (sval.lt.10) THEN
                READ(temp(1:sval-1),*) resnum(i)
             END IF

          END IF
          READ(6524,'(A)') lign80
          IF (levcfg.gt.0) THEN
             READ(6524,'(A)') lign80 ! skip velocities
          END IF
          IF (levcfg.gt.1) THEN
             READ(6524,'(A)') lign80 ! skip forces
          END IF
       END IF
    END DO

    CLOSE(6524)

    WRITE(6,'(A,I5,A)') "Trajectory contains ",natomq, " atoms,"
    WRITE(6,'(A,I5,A)') "of which ",natom, " atoms are not SOL."
    
    DO i=1,natom
       atom(i)="ATOM  "
       chain(i)="A"
       elem(i)="  "
       chag(i)="  "
       occ(i)=1.d0
       bfac(i)=0.d0
       mass(i)=1.d0
    END DO

  END SUBROUTINE read_CONFIG
  
  SUBROUTINE read_HISTORY(trajfile,natomq,natom,x,xav,start,end,nframes)
    USE nrtype
    IMPLICIT NONE
    CHARACTER :: trajfile*120, lign80*80, name*4
    INTEGER :: levcfg, imcon, megatm, frame, records, i, j, natom, natomq, nframes, start, end, atnum, io, k
    REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: x
    REAL(DP),ALLOCATABLE,DIMENSION(:) :: xav
    REAL(DP) :: mass,charge,val
    
    OPEN(file=trajfile,form='FORMATTED',status='OLD',unit=6524)
    
    READ(6524,'(A)') lign80 ! Header line
    READ(6524,*) levcfg, imcon, megatm, frame, records ! format info

    READ(6524,'(A)') lign80
    IF (imcon.gt.0) THEN ! if periodic bondaries skip over these lines
       READ(6524,'(A)') lign80
       READ(6524,'(A)') lign80
       READ(6524,'(A)') lign80
    END IF
    
    i=0
    j=0
    
    DO
       READ(6524,'(A)',IOSTAT=io) lign80
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
       ELSE
          IF (lign80(1:8).eq."timestep") THEN 
             EXIT
          ELSE IF (lign80(1:2).ne."OW".and.lign80(1:2).ne."HW".and. &
               lign80(1:2).ne."Na") THEN
             i=i+1
          END IF
          j=j+1
          READ(6524,'(A)') lign80 ! skip coordinates
          IF (levcfg.gt.0) THEN
             READ(6524,'(A)') lign80 ! skip velocities
          END IF
          IF (levcfg.gt.1) THEN
             READ(6524,'(A)') lign80 ! skip forces
          END IF
       END IF
    END DO

    natom=i
    natomq=j
    nframes=frame

    REWIND(6524)

    ALLOCATE(x(natom*3,nframes))
    ALLOCATE(xav(natom*3))
    xav=0.d0

    READ(6524,'(A)') lign80 ! Header line
    READ(6524,'(A)') lign80 ! format info

    i=0
    j=0
    
    DO
       READ(6524,'(A)',IOSTAT=io) lign80
       IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          STOP
          ! End of file
       ELSE IF (io < 0) THEN
          EXIT
       ELSE
          IF (lign80(1:8).eq."timestep") THEN 
             j=j+1
             i=0
             IF (imcon.gt.0) THEN ! if periodic bondaries skip over these lines
                READ(6524,'(A)') lign80
                READ(6524,'(A)') lign80
                READ(6524,'(A)') lign80
             END IF
          ELSE IF (lign80(1:2).ne."OW".and.lign80(1:2).ne."HW".and. &
               lign80(1:2).ne."Na") THEN
             i=i+1
             READ(lign80,'(A4,10X,I7,F9.3,3X,F9.3,2X,F10.4)') name,atnum,mass,charge,val
             READ(6524,*) x(3*(i-1)+1,j),x(3*(i-1)+2,j),x(3*(i-1)+3,j) 
             IF (levcfg.gt.0) THEN
                READ(6524,'(A)') lign80 ! skip velocities
             END IF
             IF (levcfg.gt.1) THEN
                READ(6524,'(A)') lign80 ! skip forces
             END IF
             xav(3*(i-1)+1)=xav(3*(i-1)+1)+x(3*(i-1)+1,j)
             xav(3*(i-1)+2)=xav(3*(i-1)+2)+x(3*(i-1)+2,j)
             xav(3*(i-1)+3)=xav(3*(i-1)+3)+x(3*(i-1)+3,j)
          ELSE
             READ(6524,'(A)') lign80
             IF (levcfg.gt.0) THEN
                READ(6524,'(A)') lign80 ! skip velocities
             END IF
             IF (levcfg.gt.1) THEN
                READ(6524,'(A)') lign80 ! skip forces
             END IF
          END IF
       END IF
    END DO

    CLOSE(6524)

    xav=xav/nframes

    WRITE(6,'(A,I5,A)') "Trajectory contains ",natomq, " atoms,"
    WRITE(6,'(A,I5,A)') "of which ",natom, " atoms are not SOL."
    WRITE(6,'(A,I5,A)') "Trajectory contains ",nframes, " frames."


  END SUBROUTINE read_HISTORY
END MODULE read_files

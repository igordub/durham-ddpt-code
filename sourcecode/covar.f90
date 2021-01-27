PROGRAM covar
  USE nrtype
  USE read_files
  USE read_pdbs
  USE utils
  USE pbc_tools
  USE trajutils
  USE write_files
  IMPLICIT NONE
  CHARACTER :: prmfile*120, lign80*80, trajfile*120,dummy*120,let*1,form*120,filenamehold*120,outfile*120, &
       reffile*120,chainfile*120
  INTEGER :: natomq, natom,nres,i,io, j,start,start1,end1,end,skip,nframes,nframeshold, &
       k,btype,len,ncorners,ierr,len2,nfiles,nt,len3,ii,refnatom,endvec,startvec,ncha, num
  CHARACTER,ALLOCATABLE,DIMENSION(:) :: name*4,res*3,atom*6,chain*1,elem*2,chag*2,refname*4,refres*3,refatom*6, &
       refchain*1,refelem*2,refchag*2,chalab*1
  REAL(DP),ALLOCATABLE,DIMENSION(:) :: mass,occ,bfac,xav,volume,rootmass,d,e,x,y,z,refmass,refocc,refbfac, &
       coorref,refx,refy,refz,ave
  REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: coor,coorhold,coor_old,box,boxhold,box_old,a,b,fitval
  REAL(DP),ALLOCATABLE,DIMENSION(:,:,:) :: corners
  REAL(DP) :: det,x1(3),y1(3),z1(3),val,FindDet,dot
  INTEGER,ALLOCATABLE,DIMENSION(:) :: resnum,atnum,refresnum,refatnum
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: chares
  LOGICAL :: getoption,file_exists,caonly,fidet,mat,resmass,hetatm,pbc,ref,qexist,setchain

  IF (getoption('-help',.false.,dummy)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF

  IF(.not.getoption('-p',.true.,prmfile)) THEN
     CALL exit(0)
  END IF

  INQUIRE(file=prmfile,exist=qexist)
  IF (.not.qexist) THEN
     CALL helptext(0)
     STOP "paramter file does not exist"
  END IF

  IF(.not.getoption('-i',.true.,filenamehold)) THEN
     CALL exit(0)
  END IF

  IF(getoption('-s',.true.,dummy)) THEN
     READ(dummy,*) start
  ELSE
     start=1
  END IF
  
  IF(getoption('-e',.true.,dummy)) THEN
     READ(dummy,*) end
  ELSE
     end=100000
  END IF
  
  IF(getoption('-box',.true.,dummy)) THEN
     READ(dummy,*) btype
  ELSE
     btype=3
  END IF  

  IF(.not.getoption('-form',.true.,form)) THEN
     form="amber"
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

  IF (getoption('-het',.false.,dummy)) THEN
     hetatm=.true.
  ELSE
     hetatm=.false.
  END IF

  IF (getoption('-det',.false.,dummy)) THEN
     fidet=.true.
  ELSE
     fidet=.false.
  END IF

  IF (getoption('-mat',.false.,dummy)) THEN
     mat=.true.
  ELSE
     mat=.false.
  END IF 

  IF (getoption('-skip',.true.,dummy)) THEN
     READ(dummy,*) skip
  ELSE
     skip=1
  END IF
  
  IF (getoption('-pbc',.false.,dummy)) THEN
     pbc=.true.
  ELSE
     pbc=.false.
  END IF

  IF(getoption('-chain',.true.,chainfile)) THEN
     setchain=.true.
  ELSE
     setchain=.false.
  END IF

  IF (getoption('-ref',.true.,reffile)) THEN
     ref=.true.
     INQUIRE(file=reffile,exist=qexist)
     IF (.not.qexist) THEN
        WRITE(*,'(A)') "Reference file does not exist"
     END IF
  ELSE
     ref=.false.
  END IF

  IF (getoption('-sv',.true.,dummy)) THEN
     READ(dummy,*) startvec
  ELSE
     startvec=1
  END IF  

  IF (getoption('-ev',.true.,dummy)) THEN
     READ(dummy,*) endvec
  ELSE
     endvec=3
  END IF


!--------------------------------------------------------------------------------------

  IF (form.eq."amber") THEN ! amber input
     CALL read_amb_par(prmfile,natomq,natom,atom,atnum,name,res, &
          chain,resnum,occ,bfac,mass,elem,chag,i,.false.)
    
  ELSEIF (form.eq."gromacs") THEN ! gromacs input
     CALL read_gro(prmfile,natomq,natom,atom,atnum,name,res, &
          chain,resnum,occ,bfac,mass,elem,chag)

  ELSEIF (form.eq."dlpoly") THEN ! dl_poly input
     
     
     
  ELSEIF (form.eq."pdb") THEN ! pdb input
     
     
     
  ELSE
     len=INDEX(form,' ')
     WRITE(6,'(3A)') "a -form parameter of ", form(1:len-1), " is not known"
     CALL exit(0)
  END IF


  IF (setchain) THEN
     INQUIRE(file=chainfile,exist=qexist)
     IF (.not.qexist) THEN
        CALL helptext(0)
        STOP "chain file does not exist"
     ELSE
               
        OPEN(file=chainfile,form='FORMATTED',status='OLD',unit=7219)
        ncha=0
        DO
           READ(7219,'(A)',IOSTAT=io) lign80
           IF (io > 0) THEN
              WRITE(*,*) 'Check input.  Something was wrong'
              STOP
              ! End of file
           ELSE IF (io < 0) THEN
              EXIT
              ! Count number of atoms
           ELSE
              ncha=ncha+1
           END IF
        END DO
        
        REWIND(7219)
        
        ALLOCATE(chares(ncha,2))
        ALLOCATE(chalab(ncha))
        
        i=0
        DO
           READ(7219,'(A)',IOSTAT=io) lign80
           IF (io > 0) THEN
              WRITE(*,*) 'Check input.  Something was wrong'
              STOP
              ! End of file
           ELSE IF (io < 0) THEN
              EXIT
           ELSE
              i=i+1
              READ(lign80,'(A1,1X,I4,1X,I4)')chalab(i),chares(i,1),chares(i,2)
           END IF
        END DO
        
        CLOSE(7291)
        
        DO i=1,natom
           DO j=1,ncha
              IF (atnum(i).ge.chares(j,1).and.atnum(i).le.chares(j,2)) THEN
                 chain(i)=chalab(j)
              END IF
           END DO
        END DO
        
     END IF
     
  END IF


  ! Read trajectories
  len=index(filenamehold,'*')
  len2=index(filenamehold,' ')
  IF (len.gt.0) THEN
     nfiles=999
     nt=3
  ELSE
     nfiles=0
     nt=1
  END IF

  nframes=0
  nframeshold=0
  start1=start
  end1=end
  ALLOCATE(coor(3*natom,nframes))
  ALLOCATE(box(3*natom,nframes))
  loop_files:  DO i=0,nfiles
     DO j=1,nt
        IF (len.gt.0) THEN
           IF (j.eq.1) THEN
              IF (i.lt.10) THEN
                 WRITE(trajfile,'(A,I1,A)')filenamehold(1:len-1),i,filenamehold(len+1:len2-1)
              ELSE IF (i.lt.100) THEN
                 WRITE(trajfile,'(A,I2,A)')filenamehold(1:len-1),i,filenamehold(len+1:len2-1)
              ELSE
                 WRITE(trajfile,'(A,I3,A)')filenamehold(1:len-1),i,filenamehold(len+1:len2-1)
              END IF
           ELSE IF (j.eq.2.and.i.lt.10) THEN
              WRITE(trajfile,'(A,I2.2,A)')filenamehold(1:len-1),i,filenamehold(len+1:len2-1)
           ELSE IF (j.eq.3.and.i.lt.100) THEN
              WRITE(trajfile,'(A,I3.3,A)')filenamehold(1:len-1),i,filenamehold(len+1:len2-1)
           ELSE
              CYCLE
           END IF
           
           file_exists=.false.
           INQUIRE(FILE=trajfile,EXIST=file_exists)
           
           IF (file_exists) THEN
              len3=index(trajfile,' ')
              WRITE(*,'(A,A)') 'Reading trajectory file: ',trajfile(1:len3-1)
              
              
              start1=start1-nframes
              IF (start1.lt.0) start1=1
              end1=end1-nframeshold
              IF (end1.le.0) EXIT loop_files
              
              IF (form.eq."amber") THEN ! amber input
                 CALL read_amb_trj(trajfile,natomq,natom,coorhold,xav,start1,end1,nframeshold,boxhold)
              ELSEIF (form.eq."gromacs") THEN ! gromacs input
                 CALL read_trr_out(trajfile,natomq,natom,coorhold,xav,start1,end1,skip,nframeshold,boxhold)
              ELSEIF (form.eq."dlpoly") THEN ! dl_poly input
                      
              ELSEIF (form.eq."pdb") THEN ! pdb input   
                 
              END IF
              
              ALLOCATE(coor_old(3*natom,nframes))
              coor_old=coor
              DEALLOCATE(coor)
              
              ALLOCATE(box_old(9,nframes))
              box_old=box
              DEALLOCATE(box)
              
              nframes=nframes+nframeshold
              
              ALLOCATE(coor(3*natom,nframes))
              coor(:,1:nframes-nframeshold)=coor_old
              coor(:,nframes-nframeshold+1:nframes)=coorhold
              DEALLOCATE(coorhold)
              DEALLOCATE(coor_old)
              
              ALLOCATE(box(9,nframes))
              box(:,1:nframes-nframeshold)=box_old
              box(:,nframes-nframeshold+1:nframes)=boxhold
              DEALLOCATE(boxhold)
              DEALLOCATE(box_old)
              
              DEALLOCATE(xav)
           END IF
           
        ELSE
           trajfile=filenamehold
           
           DEALLOCATE(coor)
           DEALLOCATE(box)
           file_exists=.false.
           INQUIRE(FILE=trajfile,EXIST=file_exists)
           
           IF (file_exists) THEN
              len3=index(trajfile,' ')
              WRITE(*,'(A,A)') 'Reading trajectory file: ',trajfile(1:len3-1)
              
              IF (form.eq."amber") THEN ! amber input
                 CALL read_amb_trj(trajfile,natomq,natom,coor,xav,start,end,nframes,box)
              ELSEIF (form.eq."gromacs") THEN ! gromacs input
                 CALL read_trr_out(trajfile,natomq,natom,coor,xav,start,end,skip,nframes,box)   
              ELSEIF (form.eq."dlpoly") THEN ! dl_poly input
                      
              ELSEIF (form.eq."pdb") THEN ! pdb input  
                 CALL read_pdb2(trajfile,hetatm,natom,atom,atnum,name,res,chain,resnum,coor, &
                      occ,bfac,elem,chag,start,end,skip,nframes)
                 ALLOCATE(mass(natom))
                 ALLOCATE(xav(natom*3)) 
              END IF
              
              DEALLOCATE(xav)
           END IF
           
        END IF
        
        
     END DO
  END DO loop_files
  
  IF (pbc) THEN
     
     ! RE-CENTRE BOX at (0,0,0)
     DO j=1,nframes
        DO i=1,natom
           coor(3*(i-1)+1,j)=coor(3*(i-1)+1,j)-0.5d0*box(1,j)-0.5d0*box(4,j)-0.5d0*box(7,j)
           coor(3*(i-1)+2,j)=coor(3*(i-1)+2,j)-0.5d0*box(2,j)-0.5d0*box(5,j)-0.5d0*box(8,j)
           coor(3*(i-1)+3,j)=coor(3*(i-1)+3,j)-0.5d0*box(3,j)-0.5d0*box(6,j)-0.5d0*box(9,j)
        END DO
        CALL pbcwrap(box(:,j),coor(:,j),natom,btype)
     END DO
     
     DO j=1,nframes
        CALL pbcrepare(box(:,j),coor(:,j),natom,btype,chain)
     END DO
     
     ! Make protein whole
     DO k=1,nframes
        CALL pro_fix(box(:,k),coor(:,k),natom,chain,res,btype)
     END DO
     
     ! Center box on protein
     DO k=1,nframes
        CALL pro_cent(box(:,k),coor(:,k),natom,natom,chain,res,btype)
     END DO
     
  END IF

  !----------------------------------------------------------------------------------------------------
  ! Adjusting if CA only
  IF (caonly) THEN
     j=0
     DO i=1,natom
        IF (index(name(i),'CA').gt.0) THEN
           j=j+1
           coor(3*(j-1)+1,:)=coor(3*(i-1)+1,:)
           coor(3*(j-1)+2,:)=coor(3*(i-1)+2,:)
           coor(3*(j-1)+3,:)=coor(3*(i-1)+3,:)
           atom(j)=atom(i)
           atnum(j)=atnum(i)
           name(j)=name(i)
           res(j)=res(i)
           chain(j)=chain(i)
           resnum(j)=resnum(i)
           occ(j)=occ(i)
           bfac(j)=bfac(i)
           elem(j)=elem(i)
           chag(j)=chag(i)
           mass(j)=mass(i)
        END IF
     END DO
     natom=j
     
     WRITE(6,'(A,I5,A)') "Taken only Ca atoms: ", natom," atoms." 
     
     IF (resmass) THEN
        CALL res_mass(natom,res,mass)
        WRITE(6,'(A)') "Adjusted Ca mass to residue mass"
     END IF
  END IF
  
  ALLOCATE(xav(3*natom))
  xav=0.d0
  DO k=1,nframes
     xav=xav+coor(1:3*natom,k)
  END DO
  xav=xav/real(nframes)

  IF (ref) THEN
     CALL read_pdb(reffile,hetatm,refnatom,refatom,refatnum,refname,refres,refchain,refresnum,refx,refy,refz, &
          refocc,refbfac,refelem,refchag)
     ALLOCATE(coorref(3*refnatom))
     DO i=1,refnatom
        coorref(3*(i-1)+1)=refx(i)
        coorref(3*(i-1)+2)=refy(i)
        coorref(3*(i-1)+3)=refz(i)
     END DO

     DEALLOCATE(refatom)
     DEALLOCATE(refatnum)
     DEALLOCATE(refres)
     DEALLOCATE(refchain)
     DEALLOCATE(refresnum)
     DEALLOCATE(refocc)
     DEALLOCATE(refbfac)
     DEALLOCATE(refelem)
     DEALLOCATE(refchag)
     DEALLOCATE(refx)
     DEALLOCATE(refy)
     DEALLOCATE(refz)

     IF (caonly) THEN
        j=0
        DO i=1,refnatom
           IF (index(refname(i),'CA').gt.0) THEN
              j=j+1
              coorref(3*(j-1)+1)=coorref(3*(i-1)+1)
              coorref(3*(j-1)+2)=coorref(3*(i-1)+2)
              coorref(3*(j-1)+3)=coorref(3*(i-1)+3)
           END IF
        END DO
        refnatom=j
     END IF

     CALL trajfit(xav(1:3*natom),coorref(1:3*natom),natom,1,ierr)

     DEALLOCATE(coorref)
     DEALLOCATE(refname)

  END IF
  
  CALL trajfit(coor(1:3*natom,:),xav(1:3*natom),natom,nframes,ierr)  

!-------------------------------------------------------------------------------------------
! Write average coordinate pdb file
  outfile='average.pdb'
  ALLOCATE(x(natom))
  ALLOCATE(y(natom))
  ALLOCATE(z(natom))
  DO i=1,natom
     x(i)=xav(3*(i-1)+1)
     y(i)=xav(3*(i-1)+2)
     z(i)=xav(3*(i-1)+3)
  END DO

  CALL write_pdb(outfile,atom,atnum,name,res,chain,resnum,x,y,z,&
       occ,bfac,mass,elem,chag,natom,caonly,hetatm)

  DEALLOCATE(x)
  DEALLOCATE(y)
  DEALLOCATE(z)

  WRITE(6,'(A)') "Wrote average structure pdb"
  
! Finished writting pdb file
!-------------------------------------------------------------------------------------------
  ALLOCATE(rootmass(natom))
  DO i=1,natom
     rootmass(i)=mass(i)**0.5
  END DO
  
  DO i=1,nframes
     coor(1:3*natom,i)=coor(1:3*natom,i)-xav(1:3*natom)
     DO j=1,natom
        coor(3*(j-1)+1,i)=coor(3*(j-1)+1,i)*rootmass(j)
        coor(3*(j-1)+2,i)=coor(3*(j-1)+2,i)*rootmass(j)
        coor(3*(j-1)+3,i)=coor(3*(j-1)+3,i)*rootmass(j)
     END DO
  END DO

  DEALLOCATE(xav)
  DEALLOCATE(rootmass)
!----------------------------------------------------------------------------------------
! Write the covarience matrix

  WRITE(6,'(A)') "Writting covarience matrix"

  WRITE(6,'(A)',ADVANCE='NO') "Wrote atom "

  ALLOCATE(a(3*natom,3*natom))

  DO i=1,natom
     ii=3*i-2 
     DO j=ii,3*natom
        val=0
       DO k=1,nframes
           val=val+coor(ii,k)*coor(j,k)
        END DO
        val=val/real(nframes)
        a(ii,j)=val
        a(j,ii)=val
     END DO
     
     DO j=ii+1,3*natom
        val=0
        DO k=1,nframes
           val=val+coor(ii+1,k)*coor(j,k)
        END DO
        val=val/real(nframes)
        a(ii+1,j)=val
        a(j,ii+1)=val
     END DO
     
     DO j=ii+2,3*natom
        val=0
        DO k=1,nframes
           val=val+coor(ii+2,k)*coor(j,k)
        END DO
        val=val/real(nframes)
        a(ii+1,j)=val
        a(j,ii+1)=val
     END DO
     
     IF (int(i/(int(natom/10)))*(int(natom/10)).eq.i) THEN
        WRITE(6,'(1X,I5)',ADVANCE='NO') i
     END IF
  END DO

  WRITE(6,'(/A)')'Covariance matrix ready.'

  WRITE(6,'(/A)')'Writing atom rmsd data - rmsd.dat'

  OPEN(file='rmsd.dat',form='FORMATTED',unit=9439)
  WRITE(9439,'(A)') "#   atom         rmsd    "
  DO i=1,natom
     WRITE(9439,'(I10,1PG20.12)') i, &
          ((a(3*(i-1)+1,3*(i-1)+1)+a(3*(i-1)+2,3*(i-1)+2)+a(3*(i-1)+3,3*(i-1)+3))/mass(i))**0.5
  END DO
  CLOSE(9439)
     
! Write out matrix if needed
  IF (mat) THEN
     WRITE(6,'(/A)')'Writting Covariance matrix in matrix.sdijf'
     OPEN(file='matrix.sdijf',form='FORMATTED',unit=9432)
     DO i=1,3*natom
        DO j=i,3*natom  
           IF (a(i,j).ne.0) THEN
              WRITE(9432,'(2I10,1PG20.12)') i,j,a(i,j)
           END IF
        END DO
     END DO
     CLOSE(9432)
  END IF

  ALLOCATE(d(3*natom))
  ALLOCATE(e(3*natom)) 

  IF (fidet) THEN
     WRITE(6,'(/A)')'Calculating Covariance matrix determinate'
     ALLOCATE(b(3*natom,3*natom))
     b=0.152324397*298*a+1.d0
     det=FindDet(b,3*natom)
     DEALLOCATE(b)
     WRITE(6,'(A,G8.4)') "Covarience matrix determinate = ", det
     WRITE(6,'(A,G8.4)') "Schlitter entropy / k = ", 0.5d0*log(det)
  END IF
  
! Break matrix into tri-diagonal matrix
  CALL tred2(a,d,e) 
! a is replaced by orthoganal matrix Q
! d is the diagonal elements of the tri-diagonal matrix
! e is the off diagonal matrix

! Calculates the eigenvalues and eigenvectors
  CALL tqli(d,e,a)
! d is replaced as the eigenvalues
! a becomes the eigenvectors a(:,k) is for d(k)
! e is destroyed

! sort d and a in order of eigenvalue

  DEALLOCATE(e)

  WRITE(6,'(/A)')'Writting matrix.eigenfacs file.'

  CALL sort_max(a,d,3*natom)
  
  DO i=1,3*natom
     IF (d(i).lt.1E-08) d(i)=1E-8
  END DO

  outfile='matrix.eigenfacs'
  CALL write_eigenfacs(outfile,a,0.5961752/d,3*natom,3*natom)

  DEALLOCATE(d)
  
  DEALLOCATE(mass)

  WRITE(6,'(/A)')'Writting trajectory eigenfit file.'

  num=endvec-startvec+1

  ALLOCATE(fitval(nframes,num))

  DO i=1,nframes
     DO j=startvec,endvec
        dot=0
        DO k=1,3*natom
           dot=dot+coor(k,i)*a(k,j)
        END DO
        fitval(i,j-startvec+1)=dot
     END DO
  END DO

  ALLOCATE(ave(num))

  DO i=1,num
     ave(i)=0
     DO j=start,end
        ave(i)=ave(i)+fitval(j,i)
     END DO
     ave(i)=ave(i)/real(end-start+1)
  END DO

  DO j=1,nframes
     fitval(j,:)=fitval(j,:)-ave
  END DO

  OPEN(file='fiteigen.dat',form='FORMATTED',unit=2354)

  WRITE(2354,'(A)',ADVANCE='NO') "#"
  DO i=1,num  
    WRITE(2354,'(1X,F11.5)',ADVANCE='NO') ave(i)
  END DO
  WRITE(2354,'(A)') ""

  DO i=1,nframes
    WRITE(2354,'(I6)',ADVANCE='NO') i
      DO j=1,num
        WRITE(2354,'(1X,F11.6)',ADVANCE='NO') fitval(i,j)
      END DO
    WRITE(2354,'(A)') ""
  END DO
  CLOSE(2354)


!  CALL draw_box(corners,box,nframes,volume,btype,ncorners)
!
!  ! WRITE ATOMS OUT
!  OPEN(file='out.pdb',form='FORMATTED',unit=6589)
!  WRITE(6589,'(A)') "HEADER    gro with corners of the periodic box"
!  
!  DO j=1,nframes
!     DO i=1,natom
!        WRITE(6589,'(A6,I5,1X,A4,1X,A3,1X,A1,I4,1X,3X,3F8.3,2F6.2,1X,F8.3,1X,2A2)') &
!             atom(i),atnum(i),name(i),res(i),chain(i),resnum(i),coor(3*(i-1)+1,j),coor(3*(i-1)+2,j),coor(3*(i-1)+3,j),&
!             occ(i),bfac(i),mass(i),elem(i),chag(i)
!     END DO
!     ! WRITE PBC  
!     DO i=1,ncorners
!        WRITE(6589,'(A6,I5,1X,A4,1X,A3,1X,A1,I4,1X,3X,3F8.3,2F6.2,1X,F8.3,1X,2A2)') &
!             'ATOM  ',natom+i,'  C ','COR','A',1,corners(i,1,j),corners(i,2,j),corners(i,3,j),&
!             0.0,0.0,0.0,'  ','  '
!     END DO
!     
!     WRITE(6589,'(A6,2I5)') "CONECT",natom+1,natom+2
!     WRITE(6589,'(A6,2I5)') "CONECT",natom+2,natom+3
!     WRITE(6589,'(A6,2I5)') "CONECT",natom+3,natom+4
!     WRITE(6589,'(A6,2I5)') "CONECT",natom+1,natom+4
!     
!     WRITE(6589,'(A6,2I5)') "CONECT",natom+5,natom+6
!     WRITE(6589,'(A6,2I5)') "CONECT",natom+6,natom+7
!     WRITE(6589,'(A6,2I5)') "CONECT",natom+7,natom+8
!     WRITE(6589,'(A6,2I5)') "CONECT",natom+5,natom+8
!
!     IF (btype.le.3) THEN
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+1,natom+5
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+2,natom+6
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+3,natom+7
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+4,natom+8
!        
!     ELSE IF (btype.eq.4) THEN
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+10,natom+11
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+12,natom+13
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+14,natom+15
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+9,natom+16
!        
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+9,natom+17
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+10,natom+17
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+11,natom+18
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+12,natom+18
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+13,natom+19
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+14,natom+19
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+15,natom+20
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+16,natom+20
!        
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+1,natom+17
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+2,natom+18
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+3,natom+19
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+4,natom+20
!
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+9,natom+21
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+10,natom+21
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+11,natom+22
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+12,natom+22
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+13,natom+23
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+14,natom+23
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+15,natom+24
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+16,natom+24
!
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+5,natom+21
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+6,natom+22
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+7,natom+23
!        WRITE(6589,'(A6,2I5)') "CONECT",natom+8,natom+24
!        
!     END IF
!
!     WRITE(6589,'(A)') "END MDL"
!  END DO
!  
!  CLOSE(6589)
  
END PROGRAM covar

SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
  WRITE(iunit,'(A)')"                              C  O  V  A  R                              "
  WRITE(iunit,'(A)')"                               VERSION 1.0                               "
  WRITE(iunit,'(A)')"                                                                         "  
  WRITE(iunit,'(A)')"                               Written by:                               "
  WRITE(iunit,'(A)')"                      Tom Rodgers and David Burnell                      "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program calculates the covarience matrix from a trajectory and      "
  WRITE(iunit,'(A)')"outputs an eigenfacs type file of the eigenvalues and eigenvectors which "
  WRITE(iunit,'(A)')"can be further analysed using the rest of this toolbox. The input        "
  WRITE(iunit,'(A)')"trajectories that can be used are:                                       "
  WRITE(iunit,'(A)')"    .prmtop parameter file and .mdcrd trajectory files from amber        "
  WRITE(iunit,'(A)')"    .gro parameter file and gmx_dump written trajectory file from gromacs"
  WRITE(iunit,'(A)')"    .CONFIG parameter file and .HISTORY trajectory file from dl_poly     "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"Output files:                                                            "
  WRITE(iunit,'(A)')"     average.pdb - average coordinates of the trajectory                 "
  WRITE(iunit,'(A)')"     rmsd.dat - the rmsd displacement of the atoms from the average      "
  WRITE(iunit,'(A)')"     matrix.sdijf - the covarience matrix, only with -mat option         "
  WRITE(iunit,'(A)')"     fiteigen.dat - the trajectory coordinates dot product with the      "
  WRITE(iunit,'(A)')"                    eigenvectors                                         "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"     Usage:                                                              "
  WRITE(iunit,'(A)')"             covar -form filetype -p parameterfile -i trajectoryfile     "
  WRITE(iunit,'(A)')"                   [-ca] [-res] [-s startframe] [-e endframe]            "
  WRITE(iunit,'(A)')"                   [-skip n] [-mat] [-chain chain.dat]                   "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"Option    Type      Value     Description                                "
  WRITE(iunit,'(A)')"------------------------------------------------------------             "
  WRITE(iunit,'(A)')" -form    Input               Input file types: amber, gromacs, dl_poly  "
  WRITE(iunit,'(A)')"  -p      Input               Parameter file: .prmtop, .gro, CONFIG      "
  WRITE(iunit,'(A)')"  -i      Input               Trajectory file: .mdcrd, gmx_dump, HISTORY "
  WRITE(iunit,'(A)')"  -ca     Opt                 Only uses Ca atoms                         "
  WRITE(iunit,'(A)')" -res     Opt                 Takes the mass of the residue for the Ca   "
  WRITE(iunit,'(A)')"                              mass in the analysis, only works with -ca  "
  WRITE(iunit,'(A)')"  -s      Input,Opt 1         First frame included from the trajectory   "
  WRITE(iunit,'(A)')"  -e      Input,Opt 100000    Last frame included from the trajectory    "
  WRITE(iunit,'(A)')" -skip    Input,Opt 1         Take every n frames                        "
  WRITE(iunit,'(A)')" -mat     Opt                 Output the covar matrix, matrix.sdijf      "
  WRITE(iunit,'(A)')" -chain   Input,Opt           File to list the atoms in each custom user "
  WRITE(iunit,'(A)')"                              chain                                      "
  WRITE(iunit,'(A)')"  -sv     Input,Opt 1         First eigenvector included in plot         "
  WRITE(iunit,'(A)')"  -ev     Input,Opt 3         Last eigenvector included in plot          "
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

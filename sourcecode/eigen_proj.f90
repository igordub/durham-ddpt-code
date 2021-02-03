PROGRAM eigen_proj
  USE nrtype
  USE read_files
  USE utils
  USE pbc_tools
  USE trajutils
  USE write_files
  IMPLICIT NONE
  CHARACTER :: prmfile*120, lign80*80, trajfile*120,dummy*120,let*1,form*120,filenamehold*120, &
       outfile*120,pdbfile*120,eigenfile*120
  INTEGER :: natomq, natom,nres,i,io, j,start,start1,end1,end,skip,nframes,nframeshold, &
  k,btype,len,ncorners,ierr,len2,nfiles,nt,len3,ii,startvec,endvec,num
  CHARACTER,ALLOCATABLE,DIMENSION(:) :: name*4,res*3,atom*6,chain*1,elem*2,chag*2
  REAL(DP),ALLOCATABLE,DIMENSION(:) :: mass,occ,bfac,xav,ave,x,y,z,eigenval
  REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: coor,coorhold,coor_old,box,boxhold,box_old,eigen, fitval
  REAL(DP) :: dot,x1(3),y1(3),z1(3),val
  INTEGER,ALLOCATABLE,DIMENSION(:) :: resnum,atnum
  LOGICAL :: getoption,file_exists,caonly,fidet,mat,resmass,hetatm,pbc

  IF (getoption('-help',.false.,dummy)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF

  IF(.not.getoption('-p',.true.,prmfile)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF

  IF(.not.getoption('-i',.true.,filenamehold)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF

  IF (.not.getoption('-pdb',.true.,pdbfile)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF

  IF (.not.getoption('-eigen',.true.,eigenfile)) THEN
     CALL helptext(0)
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

  IF (getoption('-pbc',.false.,dummy)) THEN
     pbc=.true.
  ELSE
     pbc=.false.
  END IF

!--------------------------------------------------------------------------------------

  IF (form.eq."amber") THEN ! amber input
     CALL read_amb_par(prmfile,natomq,natom,atom,atnum,name,res, &
          chain,resnum,occ,bfac,mass,elem,chag,i,.false.)
     DO i=1,natom
        IF (resnum(i).le.200) THEN
           chain(i)='A'
        ELSE
           chain(i)='B'
        END IF
     END DO
  ELSEIF (form.eq."gromacs") THEN ! gromacs input
     CALL read_gro(prmfile,natomq,natom,atom,atnum,name,res, &
          chain,resnum,occ,bfac,mass,elem,chag)
     ! ALLOCATE CHAIN
     let='A'
     j=0
!     DO i=1,natom
!        IF ((index(res(i),'SOL').eq.0).and.(index(res(i),'NA').eq.0).and.(index(res(i),'CL').eq.0)) THEN
!           j=j+1
!           IF (resnum(i).lt.resnum(i-1)) THEN
!              let='B'
!           END IF
!           chain(i)=let
!        END IF
!     END DO
     DO i=1,natom
        IF (resnum(i).le.200) THEN
           chain(i)='A'
        ELSE
           chain(i)='B'
        END IF
     END DO 
  ELSE
     len=INDEX(form,' ')
     WRITE(6,'(3A)') "a -form parameter of ", form(1:len-1), " is not known"
     CALL exit(0)
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

  ALLOCATE(xav(3*natom))
  xav=0.d0
  DO k=1,nframes
     xav=xav+coor(:,k)
  END DO
  xav=xav/real(nframes)

  !----------------------------------------------------------------------------------------------------
! Adjusting if CA only
  IF (caonly) THEN
     j=0
     DO i=1,natom
        IF (index(name(i),'CA').gt.0) THEN
           j=j+1
           xav(3*(j-1)+1)=xav(3*(i-1)+1)
           xav(3*(j-1)+2)=xav(3*(i-1)+2)
           xav(3*(j-1)+3)=xav(3*(i-1)+3)
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
     
  END IF

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

  CALL read_pdb(pdbfile,hetatm,natom,atom,atnum,name,res,chain,resnum,x,y,z, &
       occ,bfac,elem,chag)
  
  IF (caonly) THEN
     j=0
     DO i=1,natom
        IF (index(name(i),'CA').gt.0) THEN
           j=j+1
           xav(3*(j-1)+1)=x(i)
           xav(3*(j-1)+2)=y(i)
           xav(3*(j-1)+3)=z(i)
        END IF
     END DO
     natom=j
  END IF

  CALL trajfit(coor(1:3*natom,:),xav(1:3*natom),natom,nframes,ierr)

  CALL read_eigenfacs(eigenfile,natom,startvec,endvec,3,eigenval,eigen,num)
  
  IF (caonly) THEN
     j=0
     DO i=1,natom
        IF (index(name(i),'CA').gt.0) THEN
           j=j+1
           eigen(3*(j-1)+1,:)=eigen(3*(i-1)+1,:)
           eigen(3*(j-1)+2,:)=eigen(3*(i-1)+2,:)
           eigen(3*(j-1)+3,:)=eigen(3*(i-1)+3,:)
           atom(j)=atom(i)
           atnum(j)=atnum(i)
           name(j)=name(i)
           res(j)=res(i)
           chain(j)=chain(i)
           resnum(j)=resnum(i)
           occ(j)=occ(i)
           bfac(j)=bfac(i)
           elem(j)=elem(i)
           chag(j)=elem(i)
        END IF
     END DO
     natom=j
  END IF

  ALLOCATE(fitval(nframes,num))

  DO i=1,nframes
     DO j=1,num
        dot=0
        DO k=1,3*natom
           dot=dot+coor(k,i)*eigen(k,j)
        END DO
        fitval(i,j)=dot
     END DO
  END DO

  ALLOCATE(ave(num))

  DO i=1,num
     ave(i)=0
     DO j=1,nframes
        ave(i)=ave(i)+fitval(j,i)
     END DO
     ave(i)=ave(i)/real(nframes)
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

  DEALLOCATE(fitval)
  DEALLOCATE(coor)
  DEALLOCATE(eigen)

END PROGRAM eigen_proj

SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
  WRITE(iunit,'(A)')"                           E  G  N  P  R  O  J                           "
  WRITE(iunit,'(A)')"                               VERSION 1.0                               "
  WRITE(iunit,'(A)')"                                                                         "  
  WRITE(iunit,'(A)')"                               Written by:                               "
  WRITE(iunit,'(A)')"                      Tom Rodgers and David Burnell                      "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program performs a projection on the eigenspace defined by          "
  WRITE(iunit,'(A)')"precalculated eigenvectors and outputs a plotable file over the number of"
  WRITE(iunit,'(A)')"specified trajectories.                                                  "
  WRITE(iunit,'(A)')"The input trajectories that can be used are:                             "
  WRITE(iunit,'(A)')"    .prmtop parameter file and .mdcrd trajectory files from amber        "
  WRITE(iunit,'(A)')"    .gro parameter file and gmx_dump written trajectory file from gromacs"
  WRITE(iunit,'(A)')"    .CONFIG parameter file and .HISTORY trajectory file from dl_poly     "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"     Usage:                                                              "
  WRITE(iunit,'(A)')"           eigen_proj -form filetype -p parameterfile -i trajectoryfile  "
  WRITE(iunit,'(A)')"                   -pdb pdbfile -eigen eigenfacs [-ca] [-s startframe]   "
  WRITE(iunit,'(A)')"                   [-e endframe] [-sv startvector] [-ev endvector]       "
  WRITE(iunit,'(A)')"                   [-skip n]                                             "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"Option    Type      Value     Description                                "
  WRITE(iunit,'(A)')"------------------------------------------------------------             "
  WRITE(iunit,'(A)')" -form    Input               Input file types: amber, gromacs, dl_poly  "
  WRITE(iunit,'(A)')" -eigen   Input               Eigenvectors eigenfacs file name           "
  WRITE(iunit,'(A)')" -pdb     Input               Average pdb structure matching eigenfacs   "
  WRITE(iunit,'(A)')"  -p      Input               Parameter file: .prmtop, .gro, CONFIG      "
  WRITE(iunit,'(A)')"  -i      Input               Trajectory file: .mdcrd, gmx_dump, HISTORY "
  WRITE(iunit,'(A)')"  -ca     Opt                 Only uses Ca atoms                         "
  WRITE(iunit,'(A)')"  -s      Input,Opt 1         First frame included from the trajectory   "
  WRITE(iunit,'(A)')"  -e      Input,Opt 100000    Last frame included from the trajectory    "
  WRITE(iunit,'(A)')" -skip    Input,Opt 1         Take every n frames                        " 
  WRITE(iunit,'(A)')"  -sv     Input,Opt 1         First eigenvector included for analysis    "
  WRITE(iunit,'(A)')"  -ev     Input,Opt 3         Last eigenvector included for analysis     "
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

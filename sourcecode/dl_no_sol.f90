PROGRAM dl_no_sol
  USE nrutil
  USE read_files
  IMPLICIT NONE
  CHARACTER :: dlfile*120,trajfile*120 
  INTEGER :: natomq, natom,i,j,start,end,nframes,fram
  CHARACTER,ALLOCATABLE,DIMENSION(:) :: name*4,res*3,atom*6, &
       chain*1,elem*2,chag*2
  REAL(DP),ALLOCATABLE,DIMENSION(:) :: mass,occ,bfac,xav,charge,dis
  REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: x
  INTEGER,ALLOCATABLE,DIMENSION(:) :: resnum,atnum
  REAL(DP) :: val,step

  val=100.d0
  step=0.002
  fram=200
  dlfile='CONFIG'
  trajfile='HISTORY'

  CALL read_CONFIG(dlfile,natomq,natom,atom,atnum,name,res, &
       chain,resnum,occ,bfac,mass,elem,chag)

  CALL read_HISTORY(trajfile,natomq,natom,x,xav,start,end,nframes)

  ALLOCATE(charge(natom))
  charge=0.d0
  ALLOCATE(dis(natom))
  dis=0.d0


  OPEN(file='CONFIG_no',form='FORMATTED',unit=9874)
  OPEN(file='HISTORY_no',form='FORMATTED',unit=9875)
  
  WRITE(9874,'(A)') 'No solvent CONFIG file'
  WRITE(9874,'(1X,I9,1X,I9,1X,I9)') 0, 1, natom
  WRITE(9874,'(3F20.10)') val, 0.d0, 0.d0
  WRITE(9874,'(3F20.10)') 0.d0, val, 0.d0
  WRITE(9874,'(3F20.10)') 0.d0, 0.d0, val

  DO i=1,natom
     WRITE(9874,'(A4,10X,I7,5X,I6,A3)') name(i), atnum(i),resnum(i),res(i)
     WRITE(9874,'(F16.8,4X,F16.8,4X,F16.8)') xav(3*(i-1)+1),xav(3*(i-1)+2),xav(3*(i-1)+3)
  END DO
  
  CLOSE(9874)

  WRITE(9875,'(A)') 'No solvent HISTORY file'
  WRITE(9875,'(1X,I9,1X,I9,1X,I9,1X,I20,1X,I20)') 0, 1, natom,nframes, (natom*2+4)*nframes+2

  DO j=1,nframes
     
     WRITE(9875,'(A8,1X,I9,1X,I9,1X,I9,1X,I9,1X,F11.6,1X,F11.6)') 'timestep', (j-1)*fram+1, natom, 0, 1, &
          step, ((real(j)-1.d0)*real(fram)+1.d0)*step
     WRITE(9875,'(3F20.10)') val, 0.d0, 0.d0
     WRITE(9875,'(3F20.10)') 0.d0, val, 0.d0
     WRITE(9875,'(3F20.10)') 0.d0, 0.d0, val
     
     DO i=1,natom
        WRITE(9875,'(A4,10X,I7,F9.3,3X,F9.3,2X,F10.4)') name(i), atnum(i), mass(i), charge(i), dis(i)
        WRITE(9875,'(F16.8,4X,F16.8,4X,F16.8)') x(3*(i-1)+1,j),x(3*(i-1)+2,j),x(3*(i-1)+3,j)
     END DO
  END DO
  
  CLOSE(9875)

END PROGRAM dl_no_sol

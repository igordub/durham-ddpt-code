PROGRAM bfac_diff
  IMPLICIT NONE
  CHARACTER :: lign80*80
  REAL(8),ALLOCATABLE,DIMENSION(:) :: bfactororig, bfactor
  INTEGER,ALLOCATABLE,DIMENSION(:) :: resnum
  INTEGER :: i, num, io

  OPEN(file='mode.bfactorsorg',status='OLD',form='FORMATTED',unit=4541)
  
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
        END IF
     END IF
  END DO
  
  num=i
  ALLOCATE(bfactororig(num))
  ALLOCATE(resnum(num))
  ALLOCATE(bfactor(num))
  
  REWIND(4541)
  
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
           READ(lign80,'(22X,I4,1X,F12.4)') resnum(i), bfactororig(i)
        END IF
     END IF
  END DO
 
  CLOSE(4541)
  OPEN(file='mode.bfactorsorg',status='OLD',form='FORMATTED',unit=4542)

  i=0
  DO
     READ(4542,'(A)',IOSTAT=io) lign80
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
           READ(lign80,'(22X,I4,1X,F12.4)') resnum(i), bfactor(i)
        END IF
     END IF
  END DO

  CLOSE(4542)
  OPEN(file='mode.bfactordiff',form='FORMATTED',unit=4543)

  DO i=1,num
     WRITE(4543,'(1X,I4,1X,F12.4)') resnum(i), bfactor(i)-bfactororig(i)
  END DO

  CLOSE(4543)
  DEALLOCATE(bfactororig)
  DEALLOCATE(resnum)
  DEALLOCATE(bfactor)

END PROGRAM bfac_diff

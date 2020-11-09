PROGRAM total_en
  IMPLICIT NONE
  REAL(8) :: Z, G1, G2, S1, S2, Gsum, Ssum
  INTEGER :: n, io

  OPEN(file='mode.energy',status='OLD',form='FORMATTED',unit=5698)

  READ(5698,'(A)')
  Gsum=0
  Ssum=0

  DO
     READ(5698,*,IOSTAT=io) n, Z, G1, G2, S1, S2
     IF (io > 0) THEN
        WRITE(*,*) 'Check input.  Something was wrong'
        STOP
        ! End of file
     ELSE IF (io < 0) THEN
        EXIT
        ! Count number of atoms
     ELSE
        Gsum=Gsum+G1
        Ssum=Ssum+S1
     END IF
  END DO
  
  CLOSE(5698)
  
  OPEN(file='total.energy',form='FORMATTED',unit=5699)

  WRITE(5699,*) Gsum, Ssum

  CLOSE(5699)


END PROGRAM total_en

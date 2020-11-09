PROGRAM matrix
  IMPLICIT NONE
  INTEGER :: i, num, num1, j, io
  CHARACTER :: filename*120, lign80*80
  REAL(8),ALLOCATABLE,DIMENSION(:,:) :: a
  LOGICAL :: getoption

  IF(.not.getoption('-i',.true.,filename)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF

  num=0

  OPEN(file=filename,form="FORMATTED",status="OLD",unit=2356)
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
             READ(lign80,'(I10)') num1
             num=max(num,num1)
       END IF
    END DO
  
    ALLOCATE(a(num,num))
    a=0.d0
    
    REWIND(2356)

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
          READ(lign80,'(2I10,1PG20.12)') i, j, a(i,j)
          a(j,i)=a(i,j)
       END IF
    END DO  

    CLOSE(2356)
    OPEN(file='mat.mat',form='FORMATTED',unit=2357)

    DO i=1,num
       WRITE(2357,*) a(:,i)
    END DO

    CLOSE(2357)

END PROGRAM matrix

SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
 
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

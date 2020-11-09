PROGRAM mult
  IMPLICIT NONE
  REAL(8) :: num, tot
  INTEGER :: i,count
  LOGICAL :: getoption
  CHARACTER :: dummy*120

  IF(getoption('-n',.true.,dummy)) THEN
     READ(dummy,*) count
  ELSE
     count=25
  END IF 

  OPEN(file='eigenvalues',form='FORMATTED',status='OLD',unit=2659)

  DO i=1,6

     READ(2659,*)

  END DO

  tot=1

  DO i=1,count

     READ(2659,*) num
     tot=tot*num
     
  END DO

  CLOSE(2659)

  OPEN(file='total',form='FORMATTED',unit=8975)
  WRITE(8975,*) tot
  CLOSE(8975)


END PROGRAM mult
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

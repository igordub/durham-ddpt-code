PROGRAM mult2
  IMPLICIT NONE
  REAL(8) :: num(3), tot
  INTEGER :: i


  OPEN(file='totals',form='FORMATTED',status='OLD',unit=2659)


  tot=1

  DO i=1,3

     READ(2659,*) num(i)
     
  END DO

  tot=num(3)*num(1)/num(2)/num(2)

  CLOSE(2659)

  OPEN(file='C',form='FORMATTED',unit=8975)
  WRITE(8975,*) tot
  CLOSE(8975)


END PROGRAM mult2

! Trajutils: general utilities that operate on MD trajectory streams
! Copyright (C) Charlie Laughton 2006
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License along
!   with this program; if not, write to the Free Software Foundation, Inc.,
!   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! Charlie Laughton
! School of Pharmacy
! University of Nottingham
! NG7 2RD
! UK
! charles.laughton@nottingham.ac.uk
module trajutils
!
!  a number of basic trajectory file manipulation utilities
!
!  subroutine trajfit(x,ref,natoms,nframes,ierr)
!    fits the frames in x to coordinates in ref
!    input: x(3,natoms,nframes),ref(3,natoms)
!    output: fitted coordinates in x
!
!  subroutine trajavg(x,xavg,natoms,nframes,ierr)
!    calculates average coordinates
!  
  contains
  
  subroutine trajfit(x,ref,natoms,nframes,ierr)
    USE nrutil
  implicit none
!
!  least-squares fits the frames in x to the structure ref
!
  real(DP)    ::  x(3,natoms,nframes),ref(3,natoms)
  real(DP)    ::  r(3,3),v(3)
  logical ::  findmove
  integer ::  natoms,nframes,ierr,i,j
  real(DP)    ::  xt,yt,zt,rmsd

  findmove=.true.
  do i=1,nframes
    call matfit(natoms,ref,x(1,1,i),r,v,rmsd,findmove)
    do j=1,natoms
      xt=r(1,1)*x(1,j,i)+r(1,2)*x(2,j,i)+r(1,3)*x(3,j,i)+v(1)
      yt=r(2,1)*x(1,j,i)+r(2,2)*x(2,j,i)+r(2,3)*x(3,j,i)+v(2)
      zt=r(3,1)*x(1,j,i)+r(3,2)*x(2,j,i)+r(3,3)*x(3,j,i)+v(3)
      x(1,j,i)=xt
      x(2,j,i)=yt
      x(3,j,i)=zt
    end do
  end do

  end subroutine trajfit

  subroutine trajavg(x,xavg,natoms,nframes,ierr)
    USE nrutil
  implicit none
!
!  average coordinates in the trajectory
!
  real(DP)     ::  x(3,natoms,nframes), xavg(3,natoms)
  integer  ::  natoms,nframes,ierr,i,j
  real(DP)     ::  rframes

  rframes=1.0/nframes
  xavg=0.0
  do i=1,nframes
    do j=1,natoms
      xavg(1,j)=xavg(1,j)+x(1,j,i)
      xavg(2,j)=xavg(2,j)+x(2,j,i)
      xavg(3,j)=xavg(3,j)+x(3,j,i)
    end do
  end do
  xavg=xavg*rframes
  end subroutine trajavg
end module trajutils


SUBROUTINE MATFIT (N, XA, XB, R, V, RMSE, ENTRY )
  USE nrutil
  implicit real*8 (a-h,o-z)
  INTEGER :: N
  REAL(DP) :: XA(3,N),XB(3,N),R(3,3),V(3),RMSE
  LOGICAL :: ENTRY
  !
  !     SUBROUTINE TO FIT THE COORD SET XA(3,N) TO THE SET XB(3,N)
  !     IN THE SENSE OF:
  !            XA= R*XB +V
  !     R IS A UNITARY 3.3 RIGHT HANDED ROTATION MATRIX
  !     AND V IS THE OFFSET VECTOR. THIS IS AN EXACT SOLUTION
  !
  !     IF ENTRY IS LOGICALLY FALSE ONLY THE RMS COORDINATE ERROR
  !     WILL BE RETURNED (BUT QUICKLY)
  !
  !    THIS SUBROUTINE IS A COMBINATION OF MCLACHLAN'S AND KABSCH'S
  !    TECHNIQUES. SEE
  !     KABSCH, W. ACTA CRYST A34, 827,1978
  !     MCLACHAN, A.D., J. MOL. BIO. NNN, NNNN 1978
  !     WRITTEN BY S.J. REMINGTON 11/78.
  !
  !     THIS SUBROUTINE USES THE IBM SSP EIGENVALUE ROUTINE 'EIGEN'
  !
  DIMENSION CMA(3),CMB(3),UMAT(3,3)
  XN=N
  XASQ=0.0
  XBSQ=0.0
  XNI=1.0/XN
  !
  !     ACCUMULATE UNCORRECTED (FOR C.M.) SUMS AND SQUARES
  !
  DO I=1,3
     CMA(I) = 0.0
     CMB(I) = 0.0
     DO J=1,3
        UMAT(I,J) = 0.0
     END DO
     !
     DO J=1,N
        DO K=1,3
           UMAT(I,K) = UMAT(I,K) + XA(I,J)*XB(K,J)
        END DO
        T = XA(I,J)
        CMA(I) = CMA(I) + T
        XASQ = XASQ + T*T
        T = XB(I,J)
        CMB(I) = CMB(I) + T
        XBSQ = XBSQ + T*T
     END DO
  END DO
  !
  !     SUBTRACT CM OFFSETS
  !
  DO I=1,3
     XASQ = XASQ - CMA(I)*CMA(I)*XNI
     XBSQ = XBSQ - CMB(I)*CMB(I)*XNI
     DO J=1,3
        UMAT(I,J) = (UMAT(I,J) - CMA(I)*CMB(J)*XNI)*XNI
     END DO
  END DO
  !
  !     FIT IT
  !
  CALL QKFIT( UMAT, RTSUM, R, ENTRY )
  RMSE =(XASQ + XBSQ)*XNI - 2.0*RTSUM
  IF( RMSE .LT. 0.0 ) RMSE = 0.0
  RMSE = SQRT (RMSE)
  !
  !      CALCULATE OFFSET IF ENTRY=.TRUE.
  !
  IF(.NOT. ENTRY) RETURN
  !
  DO I=1,3
     T = 0.0
     DO J=1,3
        T = T + R(I,J)*CMB(J)
     END DO
     V(I) = ( CMA(I) - T )*XNI
  END DO
  
END SUBROUTINE matfit

SUBROUTINE QKFIT (UMAT, RTSUM, R, ENTRY )
  USE nrutil
  IMPLICIT REAL(DP) (A-H,O-Z)
  DIMENSION UMAT(3,3), R(3,3)
  REAL(DP) :: r
  LOGICAL :: ENTRY
  !
  !     THE 'EIGENVALUE ONLY' ENTRY WAS
  !     ADAPTED FROM PROGRAM BY A.D. MCLACHAN 7/78
  !
  DIMENSION USQMAT(3,3),ROOT(3),A(3,3),B(3,3),UTR(6)
  !
  EQUIVALENCE (AAM,USQMAT(1,1)),(BAM,USQMAT(2,2)),(CAM,USQMAT(3,3))
  EQUIVALENCE (FAM,USQMAT(2,3)),(GAM,USQMAT(1,3)),(HAM,USQMAT(1,2))
  EQUIVALENCE( A(1,1), USQMAT(1,1) ), (UTR(1), B(1,1))
  !
  DATA EPS/1.0D-12/
!  DATA PI/3.14159265358979/
  DATA ONE/1.0/,ZERO/0.0/,TWO/2.0/,THREE/3.0/
  DATA HALF/0.5/,THIRD/0.333333333/,FORTHR/1.333333333/
  DATA USQMAT(2,1),USQMAT(3,1),USQMAT(3,2)/3*0.0/
  ISIG=1
  !
  !      IF ENTRY IS .TRUE. GET OUT THE ROTATION MATRIX
  !
  IF (ENTRY) THEN
     
     !     FORM USQ = (UT).U    (IN UPPER TRIANGULAR SYMMETRIC STORAGE MODE)
     !
     DO I=1,3
        DO J=I,3
           T = 0.0
           DO K=1,3
              T = T + UMAT(K,I)*UMAT(K,J)
           END DO
           IA = I + (J*J-J)/2
           UTR(IA) = T
        END DO
     END DO
     !%    WRITE(6,999) UTR
     !
     !     CALCULATE EIGENVALUES AND VECTORS
     !
     CALL EIGEN (UTR, A, 3, 0)
     CALL ESORT(UTR,A,3,0)
     !%    WRITE(6,999) UTR
     !
     ROOT(1) = UTR(1)
     ROOT(2) = UTR(3)
     ROOT(3) = UTR(6)
     !%    WRITE(6,999) ROOT
     !%    WRITE(6,999) A
     !
     !     SET A3 = A1 CROSS A2
     !     ROOTS ARE IN ORDER R(1) >= R(2) >= R(3) >= 0
     !
     A(1,3) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
     A(2,3) = A(3,1)*A(1,2) - A(1,1)*A(3,2)
     A(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)
     !%    WRITE(6,999) A
     !
     !     VECTOR SET B=U.A
     !
     DO I=1,3
        DO J=1,3
           T = 0.0
           DO K=1,3
              T = T + UMAT(J,K)*A(K,I)
           END DO
           B(J,I) = T
        END DO
     END DO
     !
     !      NORMALIZE B1 AND B2 AND CALCULATE B3 = B1 CROSS B2
     !
     B1 = SQRT( B(1,1)*B(1,1) + B(2,1)*B(2,1) + B(3,1)*B(3,1) )
     B2 = SQRT( B(1,2)*B(1,2) + B(2,2)*B(2,2) + B(3,2)*B(3,2) )
     DO I=1,3
        B(I,1) = B(I,1)/B1
        B(I,2) = B(I,2)/B2
     END DO
     !
     !      CHECK FOR LEFT HANDED ROTATION
     !
     B13 = B(2,1)*B(3,2) - B(3,1)*B(2,2)
     B23 = B(3,1)*B(1,2) - B(1,1)*B(3,2)
     B33 = B(1,1)*B(2,2) - B(2,1)*B(1,2)
     !
     S = B13*B(1,3) + B23*B(2,3) + B33*B(3,3)
     IF (S .LT. 0.0) ISIG = -1
     B(1,3) = B13
     B(2,3) = B23
     B(3,3) = B33
     !%    WRITE(6,999) B
     !
     !     CALCULATE ROTATION MATRIX R
     !
     DO I=1,3
        DO J=1,3
           T = 0.0
           DO K=1,3
              T = T + B(I,K)*A(J,K)
           END DO
           R(I,J) = T
        END DO
     END DO
     !
     !     RMS ERROR
     !
     DO I=1,3
        IF (ROOT(I) .LT. 0.0) ROOT(I) = 0.0
        ROOT(I) = SQRT( ROOT(I) )
     END DO
     !
     !     CHANGE SIGN OF EVAL #3 IF LEFT HANDED
     !
     IF (ISIG .LT. 0) ROOT(3)=-ROOT(3)
     RTSUM = ROOT(3) + ROOT(2) + ROOT(1)
     
  ELSE
     !
     !     CALC DET OF UMAT
     !
     DU11=(UMAT(2,2)*UMAT(3,3))-(UMAT(2,3)*UMAT(3,2))
     DU21=(UMAT(2,3)*UMAT(3,1))-(UMAT(2,1)*UMAT(3,3))
     DU31=(UMAT(2,1)*UMAT(3,2))-(UMAT(2,2)*UMAT(3,1))
     DETU=(UMAT(1,1)*DU11)+(UMAT(1,2)*DU21)+(UMAT(1,3)*DU31)
     !%    WRITE(6,999) DETU
     !%999 FORMAT(/(3F12.5))
     IF(DETU.LT. 0.0) ISIG=-1
     !
     !     FORM USQMAT AS POSITIVE SEMI DEFINITE MATRIX
     !
     DO J=1,3
        DO I=1,J
           USQMAT(I,J)=(UMAT(1,I)*UMAT(1,J))+(UMAT(2,I)*UMAT(2,J))+ &
                (UMAT(3,I)*UMAT(3,J))
        END DO
     END DO
     !%    WRITE(6,999) USQMAT
     !
     !     REDUCE AVG OF DIAGONAL TERMS TO ZERO
     !
     DIGAV=(AAM+BAM+CAM)*THIRD
     !%    WRITE(6,999) DIGAV
     AAM=AAM-DIGAV
     BAM=BAM-DIGAV
     CAM=CAM-DIGAV
     !
     !     SETUP COEFFS OF SECULAR EQUATION OF MATRIX WITH TRACE ZERO
     !
     CC=(FAM*FAM)+(GAM*GAM)+(HAM*HAM)-(AAM*BAM)-(BAM*CAM)-(CAM*AAM)
     DD=(AAM*BAM*CAM)+TWO*(FAM*GAM*HAM)-(AAM*FAM*FAM) &
          -(BAM*GAM*GAM)-(CAM*HAM*HAM)
     !
     !     THE SECULAR EQN IS Y**3-CC*Y-DD=0  AND DD IS DET(USQMAT)
     !     REDUCE THIS TO THE FORM COS**3-(3/4)COS-
     !     (1/4)COS3THETA = 0
     !     WITH SOLUTIONS COSTHETA.  SO Y=QQ*COSTHETA
     !
     IF(CC.LE.EPS) THEN!
        !     SPECIAL FOR TRIPLY DEGENERATE
        !
        ROOT(1)=0.0
        ROOT(2)=0.0
        ROOT(3)=0.0
     ELSE
        
        QQ=SQRT(FORTHR*CC)
        COS3TH=(THREE*DD)/(CC*QQ)
        IF(ABS(COS3TH).GT.ONE) COS3TH=SIGN(ONE,COS3TH)
        !
        !     FUNCTION ARCOS
        !
        IF(COS3TH.NE.0.0) THEN
           ARGSQ=COS3TH*COS3TH
           THETA=ATAN(SQRT(1.0-ARGSQ)/COS3TH)
           IF(COS3TH.LT.0.0) THETA=PI-THETA
        ELSE
           THETA= 1.570796327
        END IF
        !
        !     ROOTS IN ORDER OF SIZE GO 1,2,3 1 LARGEST
        !
        THETA=THETA*THIRD
        ROOT(1)=QQ*COS(THETA)
        DIFF=HALF*SQRT(THREE*(QQ*QQ-ROOT(1)*ROOT(1)))
        ROOT(2)=-ROOT(1)*HALF+DIFF
        ROOT(3)=-ROOT(1)*HALF-DIFF
     END IF
     
     !     ADD ON DIGAV AND TAKE SQRT
     DO J=1,3
        RT=ROOT(J)+DIGAV
        IF(RT.LT.EPS) RT=0.0
        ROOT(J)=SQRT(RT)
     END DO
     !%    WRITE(6,999) ROOT
     !     IF DETU IS <0 CHANGE SIGN OF ROOT(3)
     IF(ISIG.EQ.-1) ROOT(3)=-ROOT(3)
     RTSUM=ROOT(1)+ROOT(2)+ROOT(3)
     !%    WRITE(6,999) RTSUM
  END IF
  !
END SUBROUTINE QKFIT
!---- SUBROUTINE TO COMPUTE EIGENVALUES & EIGENVECTORS OF A REAL SYMMETRIC
!---- MATRIX, STOLEN FROM IBM SSP MANUAL (SEE P165)
!---- DESCRIPTION OF PARAMETERS -
!---- A - ORIGINAL MATRIX STORED COLUMNWISE AS UPPER TRIANGLE ONLY,
!---- I.E. "STORAGE MODE" = 1.  EIGENVALUES ARE WRITTEN INTO DIAGONAL
!---- ELEMENTS OF A  I.E.  A(1)  A(3)  A(6)  FOR A 3*3 MATRIX.
!---- R - RESULTANT MATRIX OF EIGENVECTORS STORED COLUMNWISE IN SAME
!---- ORDER AS EIGENVALUES.
!---- N - ORDER OF MATRICES A & R.
!---- MV = 0 TO COMPUTE EIGENVALUES & EIGENVECTORS.
SUBROUTINE EIGEN(A,R,N,MV)
  USE nrutil
  REAL(DP) :: A(1),R(1),ANORM,ANRMX,THR,X,Y,SINX,SINX2,COSX,COSX2,SINCS,RANGE
  !---- FOR DOUBLE PRECISION, SQRT IN STATEMENTS 40,68,75&78 MUST BE DSQRT,
  !---- ABS IN 62 MUST BE DABS AND 1.E-6 IN 5 MUST BE 1.D-12 .
  !5	RANGE=1.E-6
  RANGE=1.D-12
  IF((MV-1).ne.0) THEN
     IQ=-N
     DO J=1,N
        IQ=IQ+N
        DO I=1,N
           IJ=IQ+I
           R(IJ)=0.
           IF((I-J).eq.0) THEN
              R(IJ)=1.
           END IF
        END DO
     END DO
     !---- INITIAL AND FINAL NORMS (ANORM & ANRMX)
  END IF

  ANORM=0.
  DO I=1,N
     DO J=I,N
        IF((I-J).ne.0) THEN
           IA=I+(J*J-J)/2
           ANORM=ANORM+A(IA)**2
        END IF
     END DO
  END DO

  IF ((ANORM).gt.0) THEN
     !40	ANORM=SQRT(2.*ANORM)
     ANORM=DSQRT(2.*ANORM)
     ANRMX=ANORM*RANGE/N
     !---- INITIALIZE INDICATORS AND COMPUTE THRESHOLD
     IND=0
     THR=ANORM
     DO
        THR=THR/N
        DO
           L=1
           DO
              M=L+1
              !---- COMPUTE SIN & COS
              DO
                 MQ=(M*M-M)/2
                 LQ=(L*L-L)/2
                 LM=L+MQ
                 IF((ABS(A(LM))-THR).ge.0) THEN 
                    IND=1
                    LL=L+LQ
                    MM=M+MQ
                    X=.5*(A(LL)-A(MM))
                    Y=-A(LM)/DSQRT(A(LM)**2+X*X)
                    IF((X).lt.0) THEN
                       Y=-Y
                    END IF
                    SINX=Y/DSQRT(2.*(1.+(DSQRT(1.-Y*Y))))
                    SINX2=SINX**2
                    COSX=DSQRT(1.-SINX2)
                    COSX2=COSX**2
                    SINCS=SINX*COSX
                    !---- ROTATE L & M COLUMNS
                    ILQ=N*(L-1)
                    IMQ=N*(M-1)
                    DO I=1,N
                       IQ=(I*I-I)/2
                       IF((I-L).ne.0) THEN
                          IF((I-M).gt.0) THEN
                             IM=M+IQ
                          ELSE IF ((I-M).lt.0) THEN
                             IM=I+MQ
                          END IF
                          IF ((I-M).ne.0) THEN
                             IF((I-L).lt.0) THEN
                                IL=I+LQ
                             ELSE
                                IL=L+IQ
                             END IF
                             X=A(IL)*COSX-A(IM)*SINX
                             A(IM)=A(IL)*SINX+A(IM)*COSX
                             A(IL)=X
                          END IF
                       END IF
                       IF((MV-1).ne.0) THEN
                          ILR=ILQ+I
                          IMR=IMQ+I
                          X=R(ILR)*COSX-R(IMR)*SINX
                          R(IMR)=R(ILR)*SINX+R(IMR)*COSX
                          R(ILR)=X
                       END IF
                    END DO
                    X=2.*A(LM)*SINCS
                    Y=A(LL)*COSX2+A(MM)*SINX2-X
                    X=A(LL)*SINX2+A(MM)*COSX2+X
                    A(LM)=(A(LL)-A(MM))*SINCS+A(LM)*(COSX2-SINX2)
                    A(LL)=Y
                    A(MM)=X
                    !---- TESTS FOR COMPLETION
                    !---- TEST FOR M = LAST COLUMN   
                 END IF
                 IF((M-N).eq.0) THEN 
                    EXIT
                 ELSE
                    M=M+1
                 END IF
              END DO
              !---- TEST FOR L =PENULTIMATE COLUMN
              IF((L-(N-1)).eq.0) THEN
                 EXIT
              ELSE
                 L=L+1
              END IF
           END DO
           IF((IND-1).eq.0) THEN
              IND=0
           ELSE
              EXIT
           END IF
        END DO
        !---- COMPARE THRESHOLD WITH FINAL NORM
        IF((THR-ANRMX).le.0) THEN
           EXIT
        END IF
     END DO
     
  END IF
  
END SUBROUTINE EIGEN

SUBROUTINE ESORT(A,R,N,MV)
  IMPLICIT DOUBLE PRECISION(A-H,O-Z)
  DIMENSION A(1),R(1)
  IQ=-N
  DO I=1,N
     IQ=IQ+N
     LL=I+(I*I-I)/2
     JQ=N*(I-2)
     DO J=I,N
        JQ=JQ+N
        MM=J+(J*J-J)/2
        IF((A(LL)-A(MM)).lt.0) THEN
           X=A(LL)
           A(LL)=A(MM)
           A(MM)=X
           IF((MV-1).ne.0) THEN
              DO K=1,N
                 ILR=IQ+K
                 IMR=JQ+K
                 X=R(ILR)
                 R(ILR)=R(IMR)
                 R(IMR)=X
              END DO
           END IF
        END IF
     END DO
  END DO
  
END SUBROUTINE ESORT

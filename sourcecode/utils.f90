!CONTAINS:
! utils module
! pythag function
!
MODULE utils
! atom_mass - assigns the atom mass based on elem
! res_mass - assigns the residue mass based on res
! tred2 - turns symetric matrix into a tri-diagonal matrix
! tqli - finds the eigenvalues and vectors of a tri-diagonal matrix
! sort_max - sorts vectors by associted values, HIGH to LOW
! sort_min - sorts vectors by associted values, LOW to HIGH
! sort_array - generates a custom sorting array
! run_sort - sorts based on a sorting array
! sort_matfile - sorts an input matrix file based on a sorting array 
! vecstat - gives min, max, and deviation of values  
!
CONTAINS
  SUBROUTINE atom_mass(natom,elem,mass)
    USE nrtype
    IMPLICIT NONE
    CHARACTER,ALLOCATABLE,DIMENSION(:) :: elem*2
    INTEGER :: natom,i
    REAL(DP),ALLOCATABLE,DIMENSION(:) :: mass
    
    DO i=1,natom
       IF(elem(i).eq." H") THEN
          mass(i)=1.008
       ELSE IF(elem(i).eq." C") THEN
          mass(i)=12.011
       ELSE IF(elem(i).eq." N") THEN
          mass(i)=14.007
       ELSE IF(elem(i).eq." O") THEN
          mass(i)=15.999
       ELSE IF(elem(i).eq." P") THEN
          mass(i)=30.974
       ELSE IF(elem(i).eq." S") THEN
          mass(i)=32.066
       ELSE IF(elem(i).eq."ZN") THEN
          mass(i)=65.39
       ELSE
          WRITE(6,'(2A)')"Unknown atom type please add to list: ", elem(i)
          mass(i)=1
       END IF
    END DO
    
  END SUBROUTINE atom_mass

  SUBROUTINE res_mass(natom,res,mass)
    USE nrtype
    IMPLICIT NONE
    CHARACTER,ALLOCATABLE,DIMENSION(:) :: res*3, resin*3
    CHARACTER :: lign80*80
    INTEGER :: natom,i,io,j,num
    REAL(DP),ALLOCATABLE,DIMENSION(:) :: mass, massin
    LOGICAL :: file_exists

    DO i=1,natom
       IF (res(i).eq.'ARG') THEN
          mass(i)=174.2017-18.015
       ELSE IF ((res(i).eq.'HIS').or.(res(i).eq.'HIE')) THEN
          mass(i)=155.1552-18.015
       ELSE IF (res(i).eq.'LYS') THEN
          mass(i)=146.1882-18.015
       ELSE IF (res(i).eq.'ASP') THEN
          mass(i)=133.1032-18.015
       ELSE IF (res(i).eq.'GLU') THEN
          mass(i)=147.1299-18.015
       ELSE IF (res(i).eq.'SER') THEN
          mass(i)=105.0930-18.015
       ELSE IF (res(i).eq.'THR') THEN
          mass(i)=119.1197-18.015
       ELSE IF (res(i).eq.'ASN') THEN
          mass(i)=132.1184-18.015
       ELSE IF (res(i).eq.'GLN') THEN
          mass(i)=146.1451-18.015
       ELSE IF (res(i).eq.'CYS') THEN
          mass(i)=121.1590-18.015
       ELSE IF (res(i).eq.'GLY') THEN
          mass(i)=75.0669-18.015
       ELSE IF (res(i).eq.'PRO') THEN
          mass(i)=115.1310-18.015
       ELSE IF (res(i).eq.'ALA') THEN
          mass(i)=89.0935-18.015
       ELSE IF (res(i).eq.'VAL') THEN
          mass(i)=117.1469-18.015
       ELSE IF (res(i).eq.'ILE') THEN
          mass(i)=131.1736-18.015
       ELSE IF (res(i).eq.'LEU') THEN
          mass(i)=131.1736-18.015
       ELSE IF (res(i).eq.'MET') THEN
          mass(i)=149.2124-18.015
       ELSE IF (res(i).eq.'PHE') THEN
          mass(i)=165.1900-18.015
       ELSE IF (res(i).eq.'TYR') THEN
          mass(i)=181.1894-18.015
       ELSE IF (res(i).eq.'TRP') THEN
          mass(i)=204.2262-18.015
       ELSE IF (res(i).eq.' DT') THEN
          mass(i)=42.037
       ELSE IF (res(i).eq.' DC') THEN
          mass(i)=37.033
       ELSE IF (res(i).eq.' DG') THEN
          mass(i)=50.377
       ELSE IF (res(i).eq.' DA') THEN
          mass(i)=45.043
       ELSE
          file_exists=.false.
          INQUIRE(FILE='resmass.dat',EXIST=file_exists)
          
          IF (file_exists) THEN
             OPEN(file='resmass.dat',form='FORMATTED',status='OLD',unit=2975)
             
             j=0
             DO
                READ(2975,'(A)',IOSTAT=io) lign80
                IF (io > 0) THEN
                   WRITE(*,*) 'Check input.  Something was wrong'
                   STOP
                   ! End of file
                ELSE IF (io < 0) THEN
                   EXIT
                   ! Count number of atoms
                ELSE
                   IF ((index(lign80,'#').gt.10).or.(index(lign80,'#').eq.0)) THEN
                      j=j+1
                   END IF
                END IF
             END DO
             
             num=j
             
             REWIND(2975)
             ALLOCATE(resin(num))
             ALLOCATE(massin(num))
             
             j=0
             DO
                READ(2975,'(A)',IOSTAT=io) lign80
                IF (io > 0) THEN
                   WRITE(*,*) 'Check input.  Something was wrong'
                   STOP
                   ! End of file
                ELSE IF (io < 0) THEN
                   EXIT
                   ! Count number of atoms
                ELSE
                   IF ((index(lign80,'#').gt.10).or.(index(lign80,'#').eq.0)) THEN
                      j=j+1
                      READ(lign80,*) resin(j),massin(j)
                   END IF
                END IF
             END DO

             CLOSE(2975)

             DO j=1,num
                IF (res(i).eq.resin(j)) THEN
                   mass(i)=massin(j)
                END IF
             END DO

             DEALLOCATE(resin)
             DEALLOCATE(massin)
             
          ELSE
             WRITE(6,'(3A)')"Unknown residue type please add to resmass.dat: ", res(i), " mass value"
          END IF
       END IF
    END DO
    
  END SUBROUTINE res_mass


  SUBROUTINE tred2(a,d,e)
    ! From Numerical Recipies for Fortran 90
    USE nrtype
    USE nrutil, ONLY : assert_eq,outerprod
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
    REAL(DP), DIMENSION(:), INTENT(OUT) :: d,e
    INTEGER(I4B) :: i,j,l,n,si
    REAL(DP) :: f,g,h,hh,scale
    REAL(DP), DIMENSION(size(a,1)) :: gg
    LOGICAL(LGT), SAVE :: yesvec=.true.
    
    n=assert_eq(size(a,1),size(a,2),size(d),size(e),'tred2')
    
    do i=n,2,-1
       l=i-1
       h=0.0
       if (l > 1) then
          scale=sum(abs(a(i,1:l)))
          if (scale == 0.0) then
             e(i)=a(i,l)
          else
             a(i,1:l)=a(i,1:l)/scale
             h=sum(a(i,1:l)**2)
             f=a(i,l)
             g=-sign(sqrt(h),f)
             e(i)=scale*g
             h=h-f*g
             a(i,l)=f-g
             if (yesvec) a(1:l,i)=a(i,1:l)/h
             do j=1,l
                e(j)=(dot_product(a(j,1:j),a(i,1:j)) &
                     +dot_product(a(j+1:l,j),a(i,j+1:l)))/h
             end do
             f=dot_product(e(1:l),a(i,1:l))
             hh=f/(h+h)
             e(1:l)=e(1:l)-hh*a(i,1:l)
             do j=1,l
                a(j,1:j)=a(j,1:j)-a(i,j)*e(1:j)-e(j)*a(i,1:j)
             end do
          end if
       else
          e(i)=a(i,l)
       end if
       d(i)=h
    end do
    if (yesvec) d(1)=0.0
    e(1)=0.0
    do i=1,n
       if (yesvec) then
          l=i-1
          if (d(i) /= 0.0) then
             gg(1:l)=matmul(a(i,1:l),a(1:l,1:l))
             a(1:l,1:l)=a(1:l,1:l)-outerprod(a(1:l,i),gg(1:l))
          end if
          d(i)=a(i,i)
          a(i,i)=1.0
          a(i,1:l)=0.0
          a(1:l,i)=0.0
       else
          d(i)=a(i,i)
       end if
    end do
  END SUBROUTINE tred2
  
  
  SUBROUTINE tqli(d,e,z)
    ! From Numerical Recipies for Fortran 90
    USE nrtype
    USE nrutil, ONLY : assert_eq,nrerror
    USE nr, ONLY : pythag
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: d,e
    REAL(DP), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: z
    INTEGER(I4B) :: i,iter,l,m,n,ndum
    REAL(DP) :: b,c,dd,f,g,p,r,s
    REAL(DP), DIMENSION(size(e)) :: ff
    n=assert_eq(size(d),size(e),'tqli: n')
    if (present(z)) ndum=assert_eq(n,size(z,1),size(z,2),'tqli: ndum')
    e(:)=eoshift(e(:),1)
    do l=1,n
       iter=0
       iterate: do
          do m=l,n-1
             dd=abs(d(m))+abs(d(m+1))
             if (abs(e(m))+dd == dd) exit
          end do
          if (m == l) exit iterate
          if (iter == 30) call nrerror('too many iterations in tqli')
          iter=iter+1
          g=(d(l+1)-d(l))/(2.0_sp*e(l))
          r=pythag(g,1.0_dp)
          g=d(m)-d(l)+e(l)/(g+sign(r,g))
          s=1.0
          c=1.0
          p=0.0
          do i=m-1,l,-1
             f=s*e(i)
             b=c*e(i)
             r=pythag(f,g)
             e(i+1)=r
             if (r == 0.0) then
                d(i+1)=d(i+1)-p
                e(m)=0.0
                cycle iterate
             end if
             s=f/r
             c=g/r
             g=d(i+1)-p
             r=(d(i)-g)*s+2.0_sp*c*b
             p=s*r
             d(i+1)=g+p
             g=c*r-b
             if (present(z)) then
                ff(1:n)=z(1:n,i+1)
                z(1:n,i+1)=s*z(1:n,i)+c*ff(1:n)
                z(1:n,i)=c*z(1:n,i)-s*ff(1:n)
             end if
          end do
          d(l)=d(l)-p
          e(l)=g
          e(m)=0.0
       end do iterate
    end do
  END SUBROUTINE tqli
 

  SUBROUTINE sort_max(a,d,n)
    USE nrtype
    IMPLICIT NONE
    INTEGER :: i, j, n
    REAL(DP) :: d(n)
    REAL(DP) :: a(n,n)
    REAL(DP),ALLOCATABLE,DIMENSION(:) :: e
    REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: b
    INTEGER,ALLOCATABLE,DIMENSION(:) :: g
    REAL(DP) :: f
    
    
    ALLOCATE(b(n,n))
    ALLOCATE(e(n))
    ALLOCATE(g(n))
    
    e=d
    b=a
    DO j=1,n
       f=MAXVAL(d)
       
       DO i=1,n
          IF (f==d(i)) THEN
             g(j)=i
             d(i)=-9d9
             EXIT
          END IF
       END DO
    END DO
    
    DO i=1,n
       d(i)=e(g(i))
       a(1:n,i)=b(1:n,g(i))
    END DO
    
    DEALLOCATE(b)
    DEALLOCATE(e)
    DEALLOCATE(g)


  END SUBROUTINE sort_max

  SUBROUTINE sort_min(a,d,n)
    USE nrtype
    IMPLICIT NONE
    INTEGER :: i, j, n
    REAL(DP) :: d(n)
    REAL(DP) :: a(n,n)
    REAL(DP),ALLOCATABLE,DIMENSION(:) :: e
    REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: b
    INTEGER,ALLOCATABLE,DIMENSION(:) :: g
    REAL(DP) :: f
    
    
    ALLOCATE(b(n,n))
    ALLOCATE(e(n))
    ALLOCATE(g(n))
    
    e=d
    b=a
    DO j=1,n
       f=MINVAL(d)
       
       DO i=1,n
          IF (f==d(i)) THEN
             g(j)=i
             d(i)=9d9
             EXIT
          END IF
       END DO
    END DO
    
    DO i=1,n
       d(i)=e(g(i))
       a(1:n,i)=b(1:n,g(i))
    END DO
    
    DEALLOCATE(b)
    DEALLOCATE(e)
    DEALLOCATE(g)

  END SUBROUTINE sort_min
!
!
  SUBROUTINE sort_array(sortarray,natom,n,chain,cdom,resnum,dom)
    USE nrutil
    IMPLICIT NONE
    INTEGER :: natom, n, i, j, k
    INTEGER :: sortarray(natom)
    LOGICAL,ALLOCATABLE,DIMENSION(:) :: added
    INTEGER :: resnum(natom), dom(n,2)
    CHARACTER :: chain(natom)*1, cdom(n)*1

    sortarray=0
    ALLOCATE(added(natom))
    added=.false.

    j=0
    DO k=1,n
       DO i=1,natom
          IF (chain(i).eq.cdom(k).and.resnum(i).ge.dom(k,1).and.resnum(i).le.dom(k,2)) THEN
             j=j+1
             sortarray(j)=i
             added(i)=.true.
          END IF
       END DO
    END DO

    IF (j.lt.natom) THEN
       WRITE(6,'(A)') 'Not all atoms were given domains, these will be treated sperately.'
       
       DO i=1,natom
          IF (.not.added(i)) THEN
             j=j+1
             sortarray(j)=i
          END IF
       END DO
    END IF


  END SUBROUTINE sort_array
!
!
  SUBROUTINE run_sort_i(a,n,sortarray)
    IMPLICIT NONE
    INTEGER :: n, i
    INTEGER :: a(n), sortarray(n)
    INTEGER,ALLOCATABLE,DIMENSION(:) :: asort

    ALLOCATE(asort(n))
    
    asort=a

    DO i=1,n
       a(i)=asort(sortarray(i))
    END DO
    
    DEALLOCATE(asort)
    
  END SUBROUTINE run_sort_i
  SUBROUTINE run_sort_r(a,n,sortarray)
    USE nrutil
    IMPLICIT NONE
    INTEGER :: n, i
    INTEGER :: sortarray(n)
    REAL(DP) :: a(n)
    REAL(DP),ALLOCATABLE,DIMENSION(:) :: asort

    ALLOCATE(asort(n))
    
    asort=a

    DO i=1,n
       a(i)=asort(sortarray(i))
    END DO
    
    DEALLOCATE(asort)
    
  END SUBROUTINE run_sort_r
   SUBROUTINE run_sort_c(a,n,sortarray)
    IMPLICIT NONE
    INTEGER :: n, i
    INTEGER :: sortarray(n)
    CHARACTER*(*) :: a(n)
    CHARACTER,ALLOCATABLE,DIMENSION(:) :: asort*500

    ALLOCATE(asort(n))
    
    DO i=1,n
       READ(a(i),'(A)') asort(i)
    END DO

    DO i=1,n
       READ(asort(sortarray(i)),'(A)') a(i)
    END DO
    
    DEALLOCATE(asort)
    
  END SUBROUTINE run_sort_c 
!
!
  SUBROUTINE vecstat(vect,natom,rmin,rmax,rave,rdev)
    USE nrutil
    IMPLICIT NONE
    INTEGER :: natom
    REAL(DP) :: rave, rdev, rmax, rmin, vect(natom)
    
    INTEGER :: i
    
    rave=0.d0
    rdev=0.d0
    rmin=-9999.d0
    rmax=9999.d0
    
    DO i=1,natom
       IF (vect(i).gt.rmax.or.i.eq.1) rmax=vect(i)
       IF (vect(i).lt.rmin.or.i.eq.1) rmin=vect(i)
       rave=rave+vect(i)
       rdev=rdev+vect(i)**2.0
    END DO
    
    rave=rave/dfloat(natom) 
    rdev=rdev/dfloat(natom)-rave*rave
    IF (rdev.gt.0.d0) rdev=sqrt(rdev)
    
  END SUBROUTINE vecstat
 
END MODULE utils

FUNCTION pythag_dp(a,b)
  USE nrtype
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: a,b
  REAL(DP) :: pythag_dp
  REAL(DP) :: absa,absb
  absa=abs(a)
  absb=abs(b)
  if (absa > absb) then
     pythag_dp=absa*sqrt(1.0_dp+(absb/absa)**2)
  else
     if (absb == 0.0) then
        pythag_dp=0.0
     else
        pythag_dp=absb*sqrt(1.0_dp+(absa/absb)**2)
     end if
  end if
END FUNCTION pythag_dp

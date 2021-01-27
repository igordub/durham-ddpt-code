MODULE pbc_tools
! pbcwrap - wrap around periodic box
! pbcrepare - repare protein across periodic boundaries
! pro_cent - center protein
! pro_fix - fix protein multimers
! cross - calculates crossproduct of a 3x3 matrix
! invert - inverts a 3x3 matrix
! draw_box - draws the box
CONTAINS
  SUBROUTINE pbcwrap(bound,coor,natom,type)
    USE nrutil
    IMPLICIT NONE
    INTEGER :: natom, i, type
    REAL(DP) :: coor(3*natom),bound(9),rbound(9),det,a,b,c,rt2
    
    IF (type.eq.1) THEN ! Standard cubic boundary conditions
       
       a=1.0d0/bound(1)
       
       DO i=1,natom
          coor(3*(i-1)+1)=coor(3*(i-1)+1)-bound(1)*Anint(a*coor(3*(i-1)+1))
          coor(3*(i-1)+2)=coor(3*(i-1)+2)-bound(1)*Anint(a*coor(3*(i-1)+2))
          coor(3*(i-1)+3)=coor(3*(i-1)+3)-bound(1)*Anint(a*coor(3*(i-1)+3))
       END DO
       
    ELSEIF (type.eq.2) THEN ! Cuboid boundary conditions
       
       a=1.0d0/bound(1)
       b=1.0d0/bound(5)
       c=1.0d0/bound(9)
       
       DO i=1,natom
          coor(3*(i-1)+1)=coor(3*(i-1)+1)-bound(1)*Anint(a*coor(3*(i-1)+1))
          coor(3*(i-1)+2)=coor(3*(i-1)+2)-bound(5)*Anint(b*coor(3*(i-1)+2))
          coor(3*(i-1)+3)=coor(3*(i-1)+3)-bound(9)*Anint(c*coor(3*(i-1)+3))
       END DO
       
    ELSEIF (type.eq.3) THEN ! Parrallelpiped boundary conditions
       
       CALL invert(bound,rbound,det)
       
       DO i=1,natom
          a=rbound(1)*coor(3*(i-1)+1)+rbound(4)*coor(3*(i-1)+2)+rbound(7)*coor(3*(i-1)+3)
          b=rbound(2)*coor(3*(i-1)+1)+rbound(5)*coor(3*(i-1)+2)+rbound(8)*coor(3*(i-1)+3)
          c=rbound(3)*coor(3*(i-1)+1)+rbound(6)*coor(3*(i-1)+2)+rbound(9)*coor(3*(i-1)+3)
          
          a=a-Anint(a)
          b=b-Anint(b)
          c=c-Anint(c)
          
          coor(3*(i-1)+1)=bound(1)*a+bound(4)*b+bound(7)*c
          coor(3*(i-1)+2)=bound(2)*a+bound(5)*b+bound(8)*c
          coor(3*(i-1)+3)=bound(3)*a+bound(6)*b+bound(9)*c
       END DO
       
    ELSEIF (type.eq.4) THEN ! Truncated octohedron boundary conditions
       
       IF (.not.(Abs(bound(1)-bound(5)).lt.1.0E-6.and.Abs(bound(5)-bound(9)).lt.1.0E-6)) STOP "bound problems"
       
       a=1.0d0/bound(1)
       
       DO i=1,natom
          coor(3*(i-1)+1)=coor(3*(i-1)+1)-bound(1)*Anint(a*coor(3*(i-1)+1))
          coor(3*(i-1)+2)=coor(3*(i-1)+2)-bound(1)*Anint(a*coor(3*(i-1)+2))
          coor(3*(i-1)+3)=coor(3*(i-1)+3)-bound(1)*Anint(a*coor(3*(i-1)+3))
          
          IF ((Abs(coor(3*(i-1)+1))+Abs(coor(3*(i-1)+2))+Abs(coor(3*(i-1)+3))) .ge. 0.75d0*bound(1)) THEN
             coor(3*(i-1)+1)=coor(3*(i-1)+1)-0.5d0*Sign(bound(1),coor(3*(i-1)+1))
             coor(3*(i-1)+2)=coor(3*(i-1)+2)-0.5d0*Sign(bound(1),coor(3*(i-1)+2))
             coor(3*(i-1)+3)=coor(3*(i-1)+3)-0.5d0*Sign(bound(1),coor(3*(i-1)+3))
          END IF
       END DO
       
    ELSEIF (type.eq.5) THEN ! Rhombic dodecahedral boundary conditions
       rt2=1.4142135662373095d0
       
       IF (.not.(Abs(bound(1)-bound(5)).lt.1.0E-6.and.Abs(bound(9)-bound(1)*rt2).lt.1.0E-6)) STOP "bound problems"
       
       a=1.0d0/bound(1)
       b=1.0d0/bound(9)
       
       DO i=1,natom
          coor(3*(i-1)+1)=coor(3*(i-1)+1)-bound(1)*Anint(a*coor(3*(i-1)+1))
          coor(3*(i-1)+2)=coor(3*(i-1)+2)-bound(1)*Anint(a*coor(3*(i-1)+2))
          coor(3*(i-1)+3)=coor(3*(i-1)+3)-bound(9)*Anint(b*coor(3*(i-1)+3))
          
          IF ((Abs(coor(3*(i-1)+1))+Abs(coor(3*(i-1)+2))+Abs(rt2*coor(3*(i-1)+3))).ge.bound(1)) THEN
             coor(3*(i-1)+1)=coor(3*(i-1)+1)-0.5d0*Sign(bound(1),coor(3*(i-1)+1))
             coor(3*(i-1)+2)=coor(3*(i-1)+2)-0.5d0*Sign(bound(1),coor(3*(i-1)+2))
             coor(3*(i-1)+3)=coor(3*(i-1)+3)-0.5d0*Sign(bound(9),coor(3*(i-1)+3))
          END IF
       END DO
       
    END IF
    
  END SUBROUTINE pbcwrap
  
  SUBROUTINE pbcrepare(bound,coor,natom,type,chain)
    USE nrutil
    IMPLICIT NONE
    INTEGER :: natom, i, type
    CHARACTER :: chain(natom)*1
    REAL(DP) :: coor(3*natom),bound(9),rbound(9),det,a,b,c,rt2,dx,dy,dz
    
    IF (type.eq.1) THEN ! Standard cubic boundary conditions
       
       a=1.0d0/bound(1)
       
       DO i=2,natom
          IF (chain(i).eq.chain(i-1)) THEN
             dx=coor(3*(i-1)+1)-coor(3*(i-2)+1)
             dy=coor(3*(i-1)+2)-coor(3*(i-2)+2)
             dz=coor(3*(i-1)+3)-coor(3*(i-2)+3)
             
             dx=dx-bound(1)*Anint(a*dx)
             dy=dy-bound(1)*Anint(a*dy)
             dz=dz-bound(1)*Anint(a*dz)
             
             coor(3*(i-1)+1)=coor(3*(i-2)+1)+dx
             coor(3*(i-1)+2)=coor(3*(i-2)+2)+dy
             coor(3*(i-1)+3)=coor(3*(i-2)+3)+dz 
          END IF
       END DO
       
    ELSEIF (type.eq.2) THEN ! Cuboid boundary conditions
       
       a=1.0d0/bound(1)
       b=1.0d0/bound(5)
       c=1.0d0/bound(9)
       
       DO i=2,natom
          IF (chain(i).eq.chain(i-1)) THEN
             dx=coor(3*(i-1)+1)-coor(3*(i-2)+1)
             dy=coor(3*(i-1)+2)-coor(3*(i-2)+2)
             dz=coor(3*(i-1)+3)-coor(3*(i-2)+3)
             
             dx=dx-bound(1)*Anint(a*dx)
             dy=dy-bound(5)*Anint(b*dy)
             dz=dz-bound(9)*Anint(c*dz)

             coor(3*(i-1)+1)=coor(3*(i-2)+1)+dx
             coor(3*(i-1)+2)=coor(3*(i-2)+2)+dy
             coor(3*(i-1)+3)=coor(3*(i-2)+3)+dz
          END IF
       END DO
       
    ELSEIF (type.eq.3) THEN ! Parrallelpiped boundary conditions
       
       CALL invert(bound,rbound,det)
       
       DO i=2,natom
          IF (chain(i).eq.chain(i-1)) THEN
             dx=coor(3*(i-1)+1)-coor(3*(i-2)+1)
             dy=coor(3*(i-1)+2)-coor(3*(i-2)+2)
             dz=coor(3*(i-1)+3)-coor(3*(i-2)+3)
             
             a=rbound(1)*dx+rbound(4)*dy+rbound(7)*dz
             b=rbound(2)*dx+rbound(5)*dy+rbound(8)*dz
             c=rbound(3)*dx+rbound(6)*dy+rbound(9)*dz
             
             a=a-Anint(a)
             b=b-Anint(b)
             c=c-Anint(c)
             
             dx=bound(1)*a+bound(4)*b+bound(7)*c
             dy=bound(2)*a+bound(5)*b+bound(8)*c
             dz=bound(3)*a+bound(6)*b+bound(9)*c

             coor(3*(i-1)+1)=coor(3*(i-2)+1)+dx
             coor(3*(i-1)+2)=coor(3*(i-2)+2)+dy
             coor(3*(i-1)+3)=coor(3*(i-2)+3)+dz
          END IF
       END DO
       
    ELSEIF (type.eq.4) THEN ! Truncated octohedron boundary conditions
       
       IF (.not.(Abs(bound(1)-bound(5)).lt.1.0E-6.and.Abs(bound(5)-bound(9)).lt.1.0E-6)) STOP "bound problems"
       
       a=1.0d0/bound(1)
       
       DO i=2,natom
          IF (chain(i).eq.chain(i-1)) THEN
             dx=coor(3*(i-1)+1)-coor(3*(i-2)+1)
             dy=coor(3*(i-1)+2)-coor(3*(i-2)+2)
             dz=coor(3*(i-1)+3)-coor(3*(i-2)+3)
             
             dx=dx-bound(1)*Anint(a*dx)
             dy=dy-bound(1)*Anint(a*dy)
             dz=dz-bound(1)*Anint(a*dz)
             
             IF ((Abs(dx)+Abs(dy)+Abs(dz)).ge.0.75d0*bound(1)) THEN
                dx=dx-0.5d0*Sign(bound(1),dx)
                dy=dy-0.5d0*Sign(bound(1),dy)
                dz=dz-0.5d0*Sign(bound(1),dz)
             END IF

             coor(3*(i-1)+1)=coor(3*(i-2)+1)+dx 
             coor(3*(i-1)+2)=coor(3*(i-2)+2)+dy
             coor(3*(i-1)+3)=coor(3*(i-2)+3)+dz
          END IF
       END DO
       
    ELSEIF (type.eq.5) THEN ! Rhombic dodecahedral boundary conditions
       rt2=1.4142135662373095d0
       
       IF (.not.(Abs(bound(1)-bound(5)).lt.1.0E-6.and.Abs(bound(9)-bound(1)*rt2).lt.1.0E-6)) STOP "bound problems"
       
       a=1.0d0/bound(1)
       b=1.0d0/bound(9)
       
       DO i=2,natom
          IF (chain(i).eq.chain(i-1)) THEN
             dx=coor(3*(i-1)+1)-coor(3*(i-2)+1)
             dy=coor(3*(i-1)+2)-coor(3*(i-2)+2)
             dz=coor(3*(i-1)+3)-coor(3*(i-2)+3)

             dx=dx-bound(1)*Anint(a*dx)
             dy=dy-bound(1)*Anint(a*dy)
             dz=dz-bound(9)*Anint(b*dz)
             
             IF ((Abs(dx)+Abs(dy)+Abs(dz)).ge.bound(1)) THEN
                dx=dx-0.5d0*Sign(bound(1),dx)
                dy=dy-0.5d0*Sign(bound(1),dy)
                dz=dz-0.5d0*Sign(bound(9),dz)
             END IF

             coor(3*(i-1)+1)=coor(3*(i-2)+1)+dx
             coor(3*(i-1)+2)=coor(3*(i-2)+2)+dy
             coor(3*(i-1)+3)=coor(3*(i-2)+3)+dz
          END IF
       END DO
          
    END IF

  END SUBROUTINE pbcrepare
  
  SUBROUTINE pro_cent(box,coor,natom,nprot,chain,res,btype)
    USE nrutil
    INTEGER :: natom, j, i,btype, nprot
    CHARACTER :: chain(natom)*1,res(natom)*3
    REAL(DP) :: coor(3*natom),box(9),a,b,c
    
    j=0
    a=0
    b=0
    c=0
    
    DO i=1,nprot
       a=a+coor(3*(i-1)+1)
       b=b+coor(3*(i-1)+2)
       c=c+coor(3*(i-1)+3)
    END DO
    a=a/real(nprot)
    b=b/real(nprot)
    c=c/real(nprot)
    
    ! Center the box on the protein center
    DO i=1,natom
       coor(3*(i-1)+1)=coor(3*(i-1)+1)-a
       coor(3*(i-1)+2)=coor(3*(i-1)+2)-b
       coor(3*(i-1)+3)=coor(3*(i-1)+3)-c
    END DO
    
  END SUBROUTINE pro_cent
  
  SUBROUTINE pro_fix(box,coor,natom,chain,res,btype)
    USE nrutil
    IMPLICIT NONE
    INTEGER :: natom, j, i, nchain, k,btype
    CHARACTER :: chain(natom)*1,res(natom)*3
    REAL(DP) :: coor(3*natom),box(9),rbox(9),a,b,c,dx(3),dx1(3),det
    REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: mass,mov
    
    j=1
    DO i=2,natom
!       IF ((index(res(i),'SOL').eq.0).and.(index(res(i),'NA').eq.0).and.(index(res(i),'CL').eq.0)) THEN
          IF (chain(i).ne.chain(i-1)) THEN
             j=j+1
          END IF
!       END IF
    END DO
    
    nchain=j
    
    ALLOCATE(mass(3,nchain))
    ALLOCATE(mov(3,nchain))
    mov=0.d0
    
    k=1
    j=1
    a=coor(1)
    b=coor(2)
    c=coor(3)
    
    DO i=2,natom
!       IF ((index(res(i),'SOL').eq.0).and.(index(res(i),'NA').eq.0).and.(index(res(i),'CL').eq.0)) THEN
          IF (chain(i).ne.chain(i-1)) THEN
             mass(1,k)=a/real(j)
             mass(2,k)=b/real(j)
             mass(3,k)=c/real(j)
             k=k+1
             a=coor(3*(i-1)+1)
             b=coor(3*(i-1)+2)
             c=coor(3*(i-1)+3)
             j=1
          ELSE
             j=j+1
             a=a+coor(3*(i-1)+1)
             b=b+coor(3*(i-1)+2)
             c=c+coor(3*(i-1)+3)
          END IF
!       END IF
    END DO
    mass(1,k)=a/real(j)
    mass(2,k)=b/real(j)
    mass(3,k)=c/real(j)

    DO i=2,nchain
       dx(1)=mass(1,1)-mass(1,i)
       dx(2)=mass(2,1)-mass(2,i)
       dx(3)=mass(3,1)-mass(3,i)
       
       dx1=dx
       
       CALL pbcwrap(box,dx1,1,btype)
       
       mov(:,i)=dx-dx1

    END DO

    j=1
    DO i=1,natom
!       IF ((index(res(i),'SOL').eq.0).and.(index(res(i),'NA').eq.0).and.(index(res(i),'CL').eq.0)) THEN
          IF (chain(i).ne.chain(i-1)) THEN
             j=j+1
          END IF
          coor(3*(i-1)+1)=coor(3*(i-1)+1)-mov(1,j)
          coor(3*(i-1)+2)=coor(3*(i-1)+2)-mov(2,j)
          coor(3*(i-1)+3)=coor(3*(i-1)+3)-mov(3,j)
!       END IF
    END DO
    
  END SUBROUTINE pro_fix
  
  FUNCTION cross(a, b)
    USE nrutil
    IMPLICIT NONE
    REAL(8), DIMENSION(3) :: cross
    REAL(8), DIMENSION(3), INTENT(IN) :: a, b
    
    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  END FUNCTION cross
  
  SUBROUTINE invert(a,b,d)
    USE nrutil
    IMPLICIT NONE
    
    REAL(8) :: a(9)
    REAL(8) :: b(9)
    REAL(8) :: d,r
    
    ! calculate adjoint matrix
    
    b(1)=a(5)*a(9)-a(6)*a(8)
    b(2)=a(3)*a(8)-a(2)*a(9)
    b(3)=a(2)*a(6)-a(3)*a(5)
    b(4)=a(6)*a(7)-a(4)*a(9)
    b(5)=a(1)*a(9)-a(3)*a(7)
    b(6)=a(3)*a(4)-a(1)*a(6)
    b(7)=a(4)*a(8)-a(5)*a(7)
    b(8)=a(2)*a(7)-a(1)*a(8)
    b(9)=a(1)*a(5)-a(2)*a(4)
    
    ! calculate determinant
    
    d=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
    r=0.0d0
    IF (Abs(d).gt.0.0d0) r=1.0d0/d
    
    ! complete inverse matrix
    
    b(1)=r*b(1)
    b(2)=r*b(2)
    b(3)=r*b(3)
    b(4)=r*b(4)
    b(5)=r*b(5)
    b(6)=r*b(6)
    b(7)=r*b(7)
    b(8)=r*b(8)
    b(9)=r*b(9)
    
  END SUBROUTINE invert
  
  SUBROUTINE draw_box(corners,box,nframes,volume,type,ncorners)
    USE nrutil
    IMPLICIT NONE
    INTEGER :: nframes, i, type, ncorners
    REAL(DP),ALLOCATABLE,DIMENSION(:,:,:) :: corners,corners1
    REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: x,y,z
    REAL(DP),ALLOCATABLE,DIMENSION(:) :: volume, a
    REAL(DP) :: box(9,nframes),x1(3),y1(3),z1(3),rot(3,3)


    IF (type.le.3) THEN
       
       ncorners=8
       
       ALLOCATE(corners(ncorners,3,nframes))
       ALLOCATE(x(3,nframes))
       ALLOCATE(y(3,nframes))
       ALLOCATE(z(3,nframes))
       x(1,:)=box(1,:)
       x(2,:)=box(2,:)
       x(3,:)=box(3,:)
       y(1,:)=box(4,:)
       y(2,:)=box(5,:)
       y(3,:)=box(6,:)
       z(1,:)=box(7,:)
       z(2,:)=box(8,:)
       z(3,:)=box(9,:)
       
       ALLOCATE(volume(nframes))
       DO i=1,nframes
          x1=x(:,i)
          y1=y(:,i)
          z1=z(:,i)
          volume(i)=abs(dot_product(x1,cross(y1,z1)))
       END DO
       
       corners(1,:,:)=-0.5d0*x-0.5d0*y-0.5d0*z
       corners(2,:,:)=0.5d0*x-0.5d0*y-0.5d0*z
       corners(3,:,:)=0.5d0*x+0.5d0*y-0.5d0*z
       corners(4,:,:)=-0.5d0*x+0.5d0*y-0.5d0*z
       corners(5,:,:)=-0.5d0*x-0.5d0*y+0.5d0*z
       corners(6,:,:)=0.5d0*x-0.5d0*y+0.5d0*z
       corners(7,:,:)=0.5d0*x+0.5d0*y+0.5d0*z
       corners(8,:,:)=-0.5d0*x+0.5d0*y+0.5d0*z
       
       DEALLOCATE(x)
       DEALLOCATE(y)
       DEALLOCATE(z)
       
       ELSEIF (type.eq.4) THEN

          ncorners=24

          ALLOCATE(corners(ncorners,3,nframes))    
          ALLOCATE(corners1(ncorners,3,nframes)) 
          ALLOCATE(a(nframes))
          a=box(1,:)/(2*2**0.5)

          ALLOCATE(volume(nframes))
          volume=8.d0*2.d0**0.5*a**3

          corners=0
          ! top square
          corners(1,1,:)=-0.5d0*a
          corners(1,2,:)=-0.5d0*a
          corners(1,3,:)=0.5d0*box(1,:)
          corners(2,1,:)=-0.5d0*a
          corners(2,2,:)=0.5d0*a
          corners(2,3,:)=0.5d0*box(1,:)
          corners(3,1,:)=0.5d0*a
          corners(3,2,:)=0.5d0*a
          corners(3,3,:)=0.5d0*box(1,:)
          corners(4,1,:)=0.5d0*a
          corners(4,2,:)=-0.5d0*a
          corners(4,3,:)=0.5d0*box(1,:)
          ! bottom square
          corners(5,1,:)=-0.5d0*a
          corners(5,2,:)=-0.5d0*a
          corners(5,3,:)=-0.5d0*box(1,:)
          corners(6,1,:)=-0.5d0*a
          corners(6,2,:)=0.5d0*a
          corners(6,3,:)=-0.5d0*box(1,:)
          corners(7,1,:)=0.5d0*a
          corners(7,2,:)=0.5d0*a
          corners(7,3,:)=-0.5d0*box(1,:)
          corners(8,1,:)=0.5d0*a
          corners(8,2,:)=-0.5d0*a
          corners(8,3,:)=-0.5d0*box(1,:)
          ! middle octagon
          corners(9,1,:)=-0.5d0*a
          corners(9,2,:)=-1.5d0*a
          corners(9,3,:)=0.d0
          corners(10,1,:)=-1.5d0*a
          corners(10,2,:)=-0.5d0*a
          corners(10,3,:)=0.d0
          corners(11,1,:)=-1.5d0*a
          corners(11,2,:)=0.5d0*a
          corners(11,3,:)=0.d0
          corners(12,1,:)=-0.5d0*a
          corners(12,2,:)=1.5d0*a
          corners(12,3,:)=0.d0
          corners(13,1,:)=0.5d0*a
          corners(13,2,:)=1.5d0*a
          corners(13,3,:)=0.d0
          corners(14,1,:)=1.5d0*a
          corners(14,2,:)=0.5d0*a
          corners(14,3,:)=0.d0
          corners(15,1,:)=1.5d0*a
          corners(15,2,:)=-0.5d0*a
          corners(15,3,:)=0.d0
          corners(16,1,:)=0.5d0*a
          corners(16,2,:)=-1.5d0*a
          corners(16,3,:)=0.d0
          ! middle top square
          corners(17,1,:)=-a
          corners(17,2,:)=-a
          corners(17,3,:)=0.707106781*a
          corners(18,1,:)=-a
          corners(18,2,:)=a
          corners(18,3,:)=0.707106781*a
          corners(19,1,:)=a
          corners(19,2,:)=a
          corners(19,3,:)=0.707106781*a
          corners(20,1,:)=a
          corners(20,2,:)=-a
          corners(20,3,:)=0.707106781*a
          ! middle bottom square
          corners(21,1,:)=-a
          corners(21,2,:)=-a
          corners(21,3,:)=-0.707106781*a
          corners(22,1,:)=-a
          corners(22,2,:)=a
          corners(22,3,:)=-0.707106781*a
          corners(23,1,:)=a
          corners(23,2,:)=a
          corners(23,3,:)=-0.707106781*a
          corners(24,1,:)=a
          corners(24,2,:)=-a
          corners(24,3,:)=-0.707106781*a

          corners1(:,3,:)=corners(:,2,:)
          corners1(:,2,:)=corners(:,3,:)
          corners1(:,1,:)=corners(:,1,:)
          corners=corners1

          DEALLOCATE(corners1)

          rot=0.d0
          rot(1,1)=0.707106781
          rot(1,3)=0.707106781
          rot(2,2)=1.d0
          rot(3,1)=-0.707106781
          rot(3,3)=0.707106781

          DO i=1,nframes
             corners(:,:,i)=matmul(corners(:,:,i),rot)
          END DO
          
       ELSEIF (type.eq.5) THEN



       END IF

  END SUBROUTINE draw_box

END MODULE pbc_tools

MODULE write_files
! write_pdb - write a pdbfile
! write_pdb_sec - writes a pdbfile with secondary structure labels
! write_eigenfacs - write the eigenfacs type file
! write_gro - write a .gro file
!
CONTAINS

  SUBROUTINE write_pdb(fname,atom,atnum,name,res,chain,resnum,x,y,z,&
           occ,bfac,mass,elem,chag,natom,caonly,hetatm)
    USE nrutil
    IMPLICIT NONE
    INTEGER :: natom,i
    CHARACTER :: fname*120
    LOGICAL :: caonly,hetatm
    INTEGER :: atnum(natom),resnum(natom)
    REAL(DP) :: x(natom),y(natom),z(natom),occ(natom),bfac(natom),mass(natom)
    CHARACTER :: atom(natom)*6,name(natom)*4,res(natom)*3,chain(natom)*1, &
         elem(natom)*2,chag(natom)*2
    
    OPEN(file=fname,form='FORMATTED',unit=6589)
    WRITE(6589,'(A)') "HEADER    GenENMM auto generated pdb file"
    
    IF (caonly) THEN
       DO i=1,natom
          IF (index(name(i),"CA").gt.0.or.(atom(i).eq."HETATM".and.hetatm).or.(name(i).eq." P  ".and.index(res(i),'D').gt.0) &
               .or.(name(i).eq." C4'".and.index(res(i),'D').gt.0) &
               .or.(name(i).eq." C2 ".and.index(res(i),'D').gt.0)) THEN
             WRITE(6589,'(A6,I5,1X,A4,1X,A3,1X,A1,I4,1X,3X,3F8.3,2F6.2,1X,F8.3,1X,2A2)') &
                  atom(i),atnum(i),name(i),res(i),chain(i),resnum(i),x(i),y(i),z(i),&
                  occ(i),bfac(i),mass(i),elem(i),chag(i)
          END IF
          IF (i.lt.natom) THEN
             IF (chain(i).ne.chain(i+1)) THEN
                WRITE(6589,'(A6,I5,6X,A3,1X,A1,I4)') "TER   ",atnum(i)+1,res(i),chain(i),resnum(i)
             END IF
          ELSE
             WRITE(6589,'(A6,I5,6X,A3,1X,A1,I4)') "TER   ",atnum(i)+1,res(i),chain(i),resnum(i)
          END IF
       END DO

    ELSE
       DO i=1,natom
          IF ((atom(i).eq."ATOM  ").or.(atom(i).eq."HETATM".and.hetatm)) THEN
             WRITE(6589,'(A6,I5,1X,A4,1X,A3,1X,A1,I4,1X,3X,3F8.3,2F6.2,1X,F8.3,1X,2A2)') &
                  atom(i),atnum(i),name(i),res(i),chain(i),resnum(i),x(i),y(i),z(i),&
                  occ(i),bfac(i),mass(i),elem(i),chag(i)
          END IF
          IF (i.lt.natom) THEN
             IF (chain(i).ne.chain(i+1)) THEN
                WRITE(6589,'(A6,I5,6X,A3,1X,A1,I4)') "TER   ",atnum(i)+1,res(i),chain(i),resnum(i)
             END IF
          ELSE
             WRITE(6589,'(A6,I5,6X,A3,1X,A1,I4)') "TER   ",atnum(i)+1,res(i),chain(i),resnum(i)
          END IF
       END DO

    END IF
       
    WRITE(6589,'(A)') "END"
    
    CLOSE(6589)
    
  END SUBROUTINE write_pdb
!
!
  SUBROUTINE write_pdb_sec(fname,atom,atnum,name,res,chain,resnum,x,y,z,&
           occ,bfac,mass,elem,chag,natom,caonly,hetatm,nseconds,seconds)
    USE nrutil
    IMPLICIT NONE
    INTEGER :: natom,i,nseconds
    CHARACTER :: fname*120
    LOGICAL :: caonly,hetatm
    INTEGER :: atnum(natom),resnum(natom)
    REAL(DP) :: x(natom),y(natom),z(natom),occ(natom),bfac(natom),mass(natom)
    CHARACTER :: atom(natom)*6,name(natom)*4,res(natom)*3,chain(natom)*1, &
         elem(natom)*2,chag(natom)*2,seconds(nseconds)*80
    
    OPEN(file=fname,form='FORMATTED',unit=6589)
    WRITE(6589,'(A)') "HEADER    GenENMM auto generated pdb file"
    
    DO i=1,nseconds
       WRITE(6589,'(A)') seconds(i)
    END DO

    IF (caonly) THEN
       DO i=1,natom
          IF (index(name(i),"CA").gt.0.or.(atom(i).eq."HETATM".and.hetatm).or.(name(i).eq." P  ".and.index(res(i),'D').gt.0) &
               .or.(name(i).eq." C4'".and.index(res(i),'D').gt.0) &
               .or.(name(i).eq." C2 ".and.index(res(i),'D').gt.0)) THEN
             WRITE(6589,'(A6,I5,1X,A4,1X,A3,1X,A1,I4,1X,3X,3F8.3,2F6.2,1X,F8.3,1X,2A2)') &
                  atom(i),atnum(i),name(i),res(i),chain(i),resnum(i),x(i),y(i),z(i),&
                  occ(i),bfac(i),mass(i),elem(i),chag(i)
          END IF
          IF (i.lt.natom) THEN
             IF (chain(i).ne.chain(i+1)) THEN
                WRITE(6589,'(A6,I5,6X,A3,1X,A1,I4)') "TER   ",atnum(i)+1,res(i),chain(i),resnum(i)
             END IF
          ELSE
             WRITE(6589,'(A6,I5,6X,A3,1X,A1,I4)') "TER   ",atnum(i)+1,res(i),chain(i),resnum(i)
          END IF
       END DO

    ELSE
       DO i=1,natom
          IF ((atom(i).eq."ATOM  ").or.(atom(i).eq."HETATM".and.hetatm)) THEN
             WRITE(6589,'(A6,I5,1X,A4,1X,A3,1X,A1,I4,1X,3X,3F8.3,2F6.2,1X,F8.3,1X,2A2)') &
                  atom(i),atnum(i),name(i),res(i),chain(i),resnum(i),x(i),y(i),z(i),&
                  occ(i),bfac(i),mass(i),elem(i),chag(i)
          END IF
          IF (i.lt.natom) THEN
             IF (chain(i).ne.chain(i+1)) THEN
                WRITE(6589,'(A6,I5,6X,A3,1X,A1,I4)') "TER   ",atnum(i)+1,res(i),chain(i),resnum(i)
             END IF
          ELSE
             WRITE(6589,'(A6,I5,6X,A3,1X,A1,I4)') "TER   ",atnum(i)+1,res(i),chain(i),resnum(i)
          END IF
       END DO

    END IF
       
    WRITE(6589,'(A)') "END"
    
    CLOSE(6589)
    
  END SUBROUTINE write_pdb_sec
!
!
  SUBROUTINE write_eigenfacs(fname,eigenvecs,eigenvals,len,natom3)
    USE nrutil
    IMPLICIT NONE
    CHARACTER :: fname*120
    INTEGER :: len, natom3, k, i
    REAL(DP) :: eigenvals(len),eigenvecs(natom3,len)
    
    OPEN(file=fname,form="FORMATTED",unit=4679)
    
    DO k=1,len
       
       WRITE(4679,'(1X,A6,I5,7X,A5,1PG12.4)') "VECTOR", k, "VALUE", eigenvals(k)
       WRITE(4679,'(A36)') " -----------------------------------"
       
       WRITE(4679,'(3(1PG12.4))') (eigenvecs(i,k),i=1,natom3)
       
    END DO
    
    CLOSE(4679)
    
  END SUBROUTINE write_eigenfacs

   SUBROUTINE write_gro(grofile,natom,num,atnum,name,res,resnum,traj)
    USE nrtype
    IMPLICIT NONE
    CHARACTER :: grofile*120
    INTEGER :: natom, num, i, j
    CHARACTER :: name(natom)*4,res(natom)*3
    REAL(DP) :: traj(3*natom,num)
    INTEGER :: resnum(natom),atnum(natom)

    OPEN(file=grofile,form='FORMATTED',unit=9324)

    DO j=1,num
       
       WRITE(9324,'(A,I5)') "Self Generated gro file. Step = ", num

       WRITE(9324,'(I6)') natom
       
       DO i=1,natom
          WRITE(9324,'(I5,A3,3X,A4,I5,3F8.3)') resnum(i),res(i),name(i),atnum(i), &
               traj(3*(i-1)+1,j),traj(3*(i-1)+2,j),traj(3*(i-1)+3,j)
       END DO

    END DO

    CLOSE(9324)

  END SUBROUTINE write_gro

END MODULE write_files

MODULE read_pdbs
! read_pdb - reads pdb file
!
CONTAINS
! Read pdb files
  SUBROUTINE read_pdb2(pdbfile,hetatm,natom,atom,atnum,name,res,chain,resnum,x, &
                  occ,bfac,elem,chag,start,end,skip,nframes)
    USE nrtype
    IMPLICIT NONE
    INTEGER, ALLOCATABLE,DIMENSION(:) :: atnum,resnum
    REAL(DP), ALLOCATABLE,DIMENSION(:) :: occ,bfac
    REAL(DP), ALLOCATABLE,DIMENSION(:,:) :: x
    CHARACTER, ALLOCATABLE, DIMENSION(:) :: atom*6,name*4,res*3,chain*1,elem*2,chag*2
    INTEGER :: io, i, natom, start, end, skip, j, k, nframes, l, nframorg
    CHARACTER :: pdbfile*120, lign80*80
    LOGICAL :: hetatm
    
    ! Find number of atoms
    OPEN(file=pdbfile,form="FORMATTED",status="OLD",unit=2356)
    i=0
    j=0
    k=0
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
          IF (lign80(1:3).eq.'END') THEN
             j=j+1
             IF (j.ge.start.and.j.le.end) THEN
                k=k+1
             END IF
          ELSE IF((lign80(1:4).eq.'ATOM').or. &
               ((lign80(1:6).eq.'HETATM').and.hetatm)) THEN
             i=i+1
          END IF
       END IF
    END DO
    
    natom=i/j
    nframorg=k
    nframes=int(k/skip)

    WRITE(6,'(A,I5,A)') "Found ", natom, " atoms in the trajectory,"
    WRITE(6,'(A,I5,A)') "Found ", j, " frames in the trajectory,"
    WRITE(6,'(A,I5,A,I5,A)') "Start frame ", start, ", End frame", min(j,end),","
    WRITE(6,'(A,I5,A)') "Number of frames in range is ", nframorg, " frames,"
    WRITE(6,'(A,I5,A)') "Taking every ", skip, " frames,"
    WRITE(6,'(A,I5,A)') "using ", nframes, " frames."
    
    REWIND(2356)
    
    ALLOCATE(x(3*natom,nframes))
    ALLOCATE(atom(natom))
    ALLOCATE(atnum(natom))
    ALLOCATE(name(natom))
    ALLOCATE(res(natom))
    ALLOCATE(chain(natom))
    ALLOCATE(resnum(natom))
    ALLOCATE(occ(natom))
    ALLOCATE(bfac(natom))
    ALLOCATE(elem(natom))
    ALLOCATE(chag(natom))  
 
  
    ! Read in atom information
    i=0

    IF (start.eq.1) THEN
       j=1
       k=1
       l=1
    ELSE
       j=0
       k=0
       l=0
    END IF
    
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
          IF (lign80(1:3).eq.'END') THEN
             j=j+1
             i=0
             IF (j.ge.start.and.j.le.end) THEN
                k=k+1
                IF (nframes.gt.10) THEN
                   IF (int(k/(nframes/10))*(nframes/10).eq.k) THEN
                      WRITE(6,'(1X,I5)',ADVANCE='NO') k
                   END IF
                ELSE
                   WRITE(6,'(1X,I5)',ADVANCE='NO') k
                END IF
                IF (int(k/skip)*skip.eq.k) THEN
 	           l=l+1
                END IF
             END IF
          ELSE IF ((lign80(1:4).eq.'ATOM').or. &
               ((lign80(1:6).eq.'HETATM').and.hetatm)) THEN
             IF (j.ge.start.and.j.le.end) THEN
                IF (int(k/skip)*skip.eq.k) THEN
                   i=i+1
                   READ(lign80,'(A6,I5,1X,A4,1X,A3,1X,A1,I4,1X,3X,3F8.3,2F6.2,10X,2A2)') &
                        atom(i),atnum(i),name(i),res(i),chain(i),resnum(i),x(3*(i-1)+1,l),x(3*(i-1)+2,l),x(3*(i-1)+3,l),&
                        occ(i),bfac(i),elem(i),chag(i)
                END IF
             END IF

          END IF
       END IF
    END DO
    
    CLOSE(2356)
    WRITE(6,'(1X)')    
  END SUBROUTINE read_pdb2
END MODULE read_pdbs

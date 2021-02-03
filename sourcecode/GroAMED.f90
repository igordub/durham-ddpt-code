PROGRAM GroAMED
  IMPLICIT NONE
  CHARACTER :: lign80*80, pdbfilename*64, atname*4, attype*3, lig*3, pdborig*64,fitfind*5,cafind*5, &
               dummy*120, filename*120, freqfile*120
  INTEGER :: io, fnum, sval, eval, natom, i, val, k, j, res, eigennum
  REAL, ALLOCATABLE, DIMENSION(:,:) :: eigenvec
  REAL, ALLOCATABLE, DIMENSION(:) :: eigenval, inc, x, y, z
  REAL :: eigenv, mass, vals(256)
  LOGICAL :: fit, ca, getoption, pca, gro, qexist

  IF (getoption('-help',.false.,dummy)) THEN
     CALL helptext(0)
     CALL exit(0)
  END IF
  
  IF(.not.getoption('-i',.true.,filename)) THEN
     WRITE(0,*) 'Specify input file with -i'
     CALL helptext(0)
     CALL exit(1)
  END IF

  INQUIRE(file=filename,exist=qexist)
  IF (.not.qexist) THEN
     STOP "eigenvec file specified doesn't exist"
  END IF
  
  IF(.not.getoption('-gropdb',.true.,pdbfilename)) THEN
     WRITE(0,*) 'Specify gromacs pdb file with -gropdb'
     CALL helptext(0)
     CALL exit(1)
  END IF

  INQUIRE(file=pdbfilename,exist=qexist)
  IF (.not.qexist) THEN
     STOP "The MD pdb file specified doesn't exist"
  END IF

  IF(.not.getoption('-pdb',.true.,pdborig)) THEN
     WRITE(0,*) 'Specify original pdb file with -pdb'
     CALL helptext(0)
     CALL exit(1)
  END IF

  INQUIRE(file=pdborig,exist=qexist)
  IF (.not.qexist) THEN
     STOP "The original pdb file specified doesn't exist"
  END IF

  IF(.not.getoption('-eigen',.true.,freqfile)) THEN
   freqfile='eigenval.xvg'  
  END IF
  
  IF(.not.getoption('-lig',.true.,lig)) THEN
     lig='NONE'
  END IF

  IF(getoption('-nofit',.false.,dummy)) THEN
     fit=.FALSE.
  ELSE
     fit=.TRUE.
  END IF

  IF(getoption('-ca',.false.,dummy)) THEN
     ca=.TRUE.
  ELSE
     ca=.FALSE.
  END IF

  IF(getoption('-mass',.true.,dummy)) THEN
     READ(dummy,*) mass
  ELSE
     mass=99.466
  END IF

  IF(getoption('-nm',.false.,dummy)) THEN
     pca=.false.
  ELSE
     pca=.true.
  END IF

  IF(getoption('-amb',.false.,dummy)) THEN
     gro=.false.
  ELSE
     gro=.true.
  END IF




  IF (gro) THEN
     
     
     IF (pca) THEN
        
        
        OPEN(file=filename,form="FORMATTED",status="OLD",unit=4568)
        
        ! Get the number of atoms and eigenvectors
        
        IF (fit) THEN
           fnum=-2
        ELSE
           fnum=-1
        END IF
        
        DO
           READ(4568,'(A)',IOSTAT=io) lign80
           IF (io > 0) THEN
              WRITE(*,*) 'Check input.  Something was wrong'
              STOP
              ! End of file
           ELSE IF (io < 0) THEN
              EXIT
           ELSE     
              IF (INDEX(lign80,'frame').gt.0) THEN
                 fnum=fnum+1
              ELSE IF (INDEX(lign80,'natoms').gt.0) THEN
                 READ(lign80,'(11X,I9)') natom
              END IF
           END IF
        END DO
        
        ALLOCATE(eigenvec(fnum*natom,3))
        ALLOCATE(eigenval(fnum))
        ALLOCATE(inc(natom))
        ALLOCATE(x(natom))
        ALLOCATE(y(natom))
        ALLOCATE(z(natom))
        
        REWIND(4568)
        
        ! Read in the eigenvectors
        IF (fit) THEN
           fnum=-2
        ELSE
           fnum=-1
        END IF
        i=0
        j=0
        DO
           READ(4568,'(A)',IOSTAT=io) lign80
           IF (io > 0) THEN
              WRITE(*,*) 'Check input.  Something was wrong'
              STOP
              ! End of file
           ELSE IF (io < 0) THEN
              EXIT
           ELSE
              IF (INDEX(lign80,'frame').gt.0) THEN
                 fnum=fnum+1
              ELSE IF (INDEX(lign80,'natoms').gt.0) THEN
              ELSE IF (INDEX(lign80,'box').gt.0) THEN
              ELSE IF (INDEX(lign80,'x3').gt.0) THEN
              ELSE
                 IF (fnum.gt.0) THEN
                    i=i+1
                    sval=INDEX(lign80,'{')
                    eval=INDEX(lign80,'}')
                    READ(lign80(sval+1:eval-1),*) eigenvec(i,1), eigenvec(i,2), eigenvec(i,3)
                 ELSE IF (fnum.eq.0) THEN
                    j=j+1
                    sval=INDEX(lign80,'{')
                    eval=INDEX(lign80,'}')
                    READ(lign80(sval+1:eval-1),*) x(j), y(j), z(j)
                 END IF
              END IF
           END IF
        END DO
        
        CLOSE(4568)
        
        OPEN(file=freqfile,form="FORMATTED",status="OLD",unit=4567)
        
        DO
           READ(4567,'(A)',IOSTAT=io) lign80
           IF (io > 0) THEN
              WRITE(*,*) 'Check input.  Something was wrong'
              STOP
              ! End of file
           ELSE IF (io < 0) THEN
              EXIT
           ELSE
              IF (lign80(1:1).eq.'#') THEN
              ELSE IF (lign80(1:1).eq.'@') THEN
              ELSE IF (lign80(1:1).eq.'&') THEN
              ELSE
                 READ(lign80,*) val, eigenv
                 IF (val.le.fnum) eigenval(val)=eigenv
              END IF
           END IF
        END DO
        
        CLOSE(4567)
        
        IF (ca) THEN
           
           OPEN(file="matrix.eigenfacs",form="FORMATTED",unit=4568)
           
           DO k=1,fnum
              
              WRITE(4568,'(1X,A6,I5,7X,A5,G12.4)') "VECTOR", k, "VALUE", 0.005961752/(eigenval(k)*mass/12.011) 
              WRITE(4568,'(A36)') " -----------------------------------"
              
              DO i=1,natom
                 WRITE(4568,'(3(1X,G11.4))') 10*eigenvec((k-1)*natom+i,1),&
                      10*eigenvec((k-1)*natom+i,2),10*eigenvec((k-1)*natom+i,3)
              END DO
              
           END DO
           
           CLOSE(4568)
           
           i=0
           OPEN(file=pdborig,form="FORMATTED",status="OLD",unit=4599)
           OPEN(file="CAonly.pdb",form="FORMATTED",unit=4600)
           DO
              READ(4599,'(A)',IOSTAT=io) lign80
              IF (io > 0) THEN
                 WRITE(*,*) 'Check input.  Something was wrong'
                 STOP
                 ! End of file
              ELSE IF (io < 0) THEN
                 EXIT
              ELSE
                 IF ((lign80(1:4).eq.'ATOM').or.(lign80(1:6).eq.'HETATM')) THEN
                    READ(lign80,'(12X,A4,1X,A3)') atname,attype
                    IF ((atname.eq.' CA ').or.(attype.eq.lig)) THEN
                       i=i+1
                       WRITE(4600,'(A30,3F8.3,A26)') lign80(1:30),10*x(i),10*y(i),10*z(i),lign80(55:80)
                    ELSE
                    END IF
                 END IF
              END IF
           END DO
           
           
           CLOSE(4599)
           CLOSE(4600)
           
        ELSE
           
           i=0
           
           OPEN(file=pdbfilename,form="FORMATTED",unit=4569)
           DO
              READ(4569,'(A)',IOSTAT=io) lign80
              IF (io > 0) THEN
                 WRITE(*,*) 'Check input.  Something was wrong'
                 STOP
                 ! End of file
              ELSE IF (io < 0) THEN
                 EXIT
              ELSE
                 IF (lign80(1:4).eq.'ATOM') THEN
                    READ(lign80,'(12X,A4,1X,A3,3X,I3)') atname,attype,res
                    i=i+1
                    IF ((atname.eq.' CA ').or.((attype.eq.lig).and..not.(atname(2:2).eq."H"))) THEN
                       inc(i)=1
                       !                 IF (res.eq.2) inc(i)=0
                    ELSE
                       inc(i)=0
                    END IF
                 END IF
              END IF
           END DO
           
           CLOSE(4569)
           
           OPEN(file="matrix.eigenfacs",form="FORMATTED",unit=4568)
           
           DO k=1,fnum
              
              WRITE(4568,'(1X,A6,I5,7X,A5,G12.4)') "VECTOR", k, "VALUE", 0.005961752/eigenval(k) 
              WRITE(4568,'(A36)') " -----------------------------------"
              
              DO i=1,natom
                 IF (inc(i).eq.1) THEN
                    WRITE(4568,'(3(1X,G11.4))') 10*eigenvec((k-1)*natom+i,1),&
                         10*eigenvec((k-1)*natom+i,2),10*eigenvec((k-1)*natom+i,3)
                 END IF
              END DO
              
           END DO
           
           
           i=0
           DO j=1,natom
              IF (inc(j).eq.1) THEN
                 i=i+1
                 x(i)=x(j)
                 y(i)=y(j)
                 z(i)=z(j)
              END IF
           END DO
           i=0
           OPEN(file=pdborig,form="FORMATTED",status="OLD",unit=4599)
           OPEN(file="CAonly.pdb",form="FORMATTED",unit=4600)
           DO
              READ(4599,'(A)',IOSTAT=io) lign80
              IF (io > 0) THEN
                 WRITE(*,*) 'Check input.  Something was wrong'
                 STOP
                 ! End of file
              ELSE IF (io < 0) THEN
                 EXIT
              ELSE
                 IF ((lign80(1:4).eq.'ATOM').or.(lign80(1:6).eq.'HETATM')) THEN
                    READ(lign80,'(12X,A4,1X,A3)') atname,attype
                    IF ((atname.eq.' CA ').or.(attype.eq.lig)) THEN
                       i=i+1
                       WRITE(4600,'(A30,3F8.3,A26)') lign80(1:30),10*x(i),10*y(i),10*z(i),lign80(55:80)
                    ELSE
                    END IF
                 END IF
              END IF
           END DO
           
           
           CLOSE(4599)
           CLOSE(4600)
           
        END IF
        
        DEALLOCATE(eigenvec)
        DEALLOCATE(eigenval)
        DEALLOCATE(inc)
        DEALLOCATE(x)
        DEALLOCATE(y)
        DEALLOCATE(z)
        
     ELSE
        
        OPEN(file=filename,form="FORMATTED",status="OLD",unit=4568)
        
        ! Get the number of atoms and eigenvectors
        fnum=-1
        DO
           READ(4568,'(A)',IOSTAT=io) lign80
           IF (io > 0) THEN
              WRITE(*,*) 'Check input.  Something was wrong'
              STOP
              ! End of file
           ELSE IF (io < 0) THEN
              EXIT
           ELSE
              IF (INDEX(lign80,'frame').gt.0) THEN
                 fnum=fnum+1
              ELSE IF (INDEX(lign80,'natoms').gt.0) THEN
                 READ(lign80,'(11X,I9)') natom
              END IF
           END IF
        END DO
        
        ALLOCATE(eigenvec(fnum*natom,3))
        ALLOCATE(eigenval(fnum))
        ALLOCATE(inc(natom))
        ALLOCATE(x(natom))
        ALLOCATE(y(natom))
        ALLOCATE(z(natom))
        
        REWIND(4568)
        
        ! Read in the eigenvectors
        fnum=-1
        i=0
        j=0
        DO
           READ(4568,'(A)',IOSTAT=io) lign80
           IF (io > 0) THEN
              WRITE(*,*) 'Check input.  Something was wrong'
              STOP
              ! End of file
           ELSE IF (io < 0) THEN
              EXIT
           ELSE
              IF (INDEX(lign80,'frame').gt.0) THEN
                 fnum=fnum+1
              ELSE IF (INDEX(lign80,'natoms').gt.0) THEN
              ELSE IF (INDEX(lign80,'box').gt.0) THEN
              ELSE IF (INDEX(lign80,'x3').gt.0) THEN
              ELSE
                 IF (fnum.gt.0) THEN
                    i=i+1
                    sval=INDEX(lign80,'{')
                    eval=INDEX(lign80,'}')
                    READ(lign80(sval+1:eval-1),*) eigenvec(i,1), eigenvec(i,2), eigenvec(i,3)
                 ELSE IF (fnum.eq.0) THEN
                    j=j+1
                    sval=INDEX(lign80,'{')
                    eval=INDEX(lign80,'}')
                    READ(lign80(sval+1:eval-1),*) x(j), y(j), z(j)
                 END IF
              END IF
           END IF
        END DO
        
        CLOSE(4568)
        
        OPEN(file=freqfile,form="FORMATTED",status="OLD",unit=4567)
        
        DO
           READ(4567,'(A)',IOSTAT=io) lign80
           IF (io > 0) THEN
              WRITE(*,*) 'Check input.  Something was wrong'
              STOP
              ! End of file
           ELSE IF (io < 0) THEN
              EXIT
           ELSE
              IF (lign80(1:1).eq.'#') THEN
              ELSE IF (lign80(1:1).eq.'@') THEN
              ELSE IF (lign80(1:1).eq.'&') THEN
              ELSE
                 READ(lign80,*) val, eigenval(val)
                 IF (eigenval(val).lt.1E-6) eigenval(val)=1E-6
                 eigenval(val)=(eigenval(val)/108.591)**2
              END IF
           END IF
        END DO
        
        CLOSE(4567)
        
        
        i=0
        
        OPEN(file=pdbfilename,form="FORMATTED",unit=4569)
        DO
           READ(4569,'(A)',IOSTAT=io) lign80
           IF (io > 0) THEN
              WRITE(*,*) 'Check input.  Something was wrong'
              STOP
              ! End of file
           ELSE IF (io < 0) THEN
              EXIT
           ELSE
              IF (lign80(1:4).eq.'ATOM') THEN
                 READ(lign80,'(12X,A4,1X,A3,3X,I3)') atname,attype,res
                 i=i+1
                 IF ((atname.eq.' CA ').or.((attype.eq.lig).and..not.(atname(2:2).eq."H"))) THEN
                    inc(i)=1
                    !              IF (res.eq.2) inc(i)=0
                 ELSE
                    inc(i)=0
                 END IF
              END IF
           END IF
        END DO
        
        CLOSE(4569)
        
        OPEN(file="matrix.eigenfacs",form="FORMATTED",unit=4568)
        
        DO k=1,fnum
           
           WRITE(4568,'(1X,A6,I5,7X,A5,G12.4)') "VECTOR", k, "VALUE", eigenval(k) 
           WRITE(4568,'(A36)') " -----------------------------------"
           
           DO i=1,natom
              IF (inc(i).eq.1) THEN
                 WRITE(4568,'(3(1X,G11.4))') 10*eigenvec((k-1)*natom+i,1),&
                      10*eigenvec((k-1)*natom+i,2),10*eigenvec((k-1)*natom+i,3)
              END IF
           END DO
           
        END DO
        
        
        i=0
        DO j=1,natom
           IF (inc(j).eq.1) THEN
              i=i+1
              x(i)=x(j)
              y(i)=y(j)
              z(i)=z(j)
           END IF
        END DO
        i=0
        OPEN(file=pdborig,form="FORMATTED",status="OLD",unit=4599)
        OPEN(file="CAonly.pdb",form="FORMATTED",unit=4600)
        DO
           READ(4599,'(A)',IOSTAT=io) lign80
           IF (io > 0) THEN
              WRITE(*,*) 'Check input.  Something was wrong'
              STOP
              ! End of file
           ELSE IF (io < 0) THEN
              EXIT
           ELSE
              IF ((lign80(1:4).eq.'ATOM').or.(lign80(1:6).eq.'HETATM')) THEN
                 READ(lign80,'(12X,A4,1X,A3)') atname,attype
                 IF ((atname.eq.' CA ').or.(attype.eq.lig)) THEN
                    i=i+1
                    WRITE(4600,'(A30,3F8.3,A26)') lign80(1:30),10*x(i),10*y(i),10*z(i),lign80(55:80)
                 ELSE
                 END IF
              END IF
           END IF
        END DO
        
        
        CLOSE(4599)
        CLOSE(4600)
        
        DEALLOCATE(eigenvec)
        DEALLOCATE(eigenval)
        DEALLOCATE(inc)
        DEALLOCATE(x)
        DEALLOCATE(y)
        DEALLOCATE(z)
        
     END IF
     
  ELSE
     
     OPEN(file=filename,form="FORMATTED",status="OLD",unit=4568)
     
     ! Get the number of atoms and eigenvectors
     fnum=-1
     
     DO
        READ(4568,'(A)',IOSTAT=io) lign80
        IF (io > 0) THEN
           WRITE(*,*) 'Check input.  Something was wrong'
           STOP
           ! End of file
        ELSE IF (io < 0) THEN
           EXIT
        ELSE
           IF (INDEX(lign80,'Eigenvector').gt.0) THEN
              fnum=fnum+1
              READ(4568,*) natom
           ELSE IF (INDEX(lign80,'****').gt.0) THEN
              fnum=fnum+1
           END IF
        END IF
     END DO
     
     natom=natom/3
     
     ALLOCATE(eigenvec(natom*3,fnum))
     ALLOCATE(eigenval(fnum))
     ALLOCATE(inc(natom))
     ALLOCATE(x(natom))
     ALLOCATE(y(natom))
     ALLOCATE(z(natom))
     
     REWIND(4568)
     
     ! Read in the eigenvectors
     fnum=-1
     i=0
     DO
        READ(4568,'(A)',IOSTAT=io) lign80
        IF (io > 0) THEN
           WRITE(*,*) 'Check input.  Something was wrong'
           STOP
           ! End of file
        ELSE IF (io < 0) THEN
           EXIT
        ELSE
           IF (INDEX(lign80,'Eigenvector').gt.0) THEN
              fnum=fnum+1
              i=0
              READ(4568,*) 
           ELSE IF (INDEX(lign80,'****').gt.0) THEN
              fnum=fnum+1
              READ(4568,*) eigennum, eigenval(eigennum)
              i=0
           ELSE
              IF (fnum.gt.0) THEN
                 vals=1000
                 READ(lign80,*,IOSTAT=io) (vals(k),k=1,256)
                 IF (io > 0) THEN
                    WRITE(*,*) 'Check input.  Something was wrong'
                    STOP
                    ! End of file
                 ELSE IF (io < 0) THEN
                 ELSE
                 END IF
                 
                 DO k=1,256
                    IF (vals(k).lt.500) THEN
                       i=i+1
                       eigenvec(i,fnum)=vals(k)
                    END IF
                 END DO
              END IF
           END IF
        END IF
     END DO
     
     CLOSE(4568)
     
     eigenval=(eigenval/108.591)**2
     
     i=0
     
     OPEN(file=pdbfilename,form="FORMATTED",unit=4569)
     DO
        READ(4569,'(A)',IOSTAT=io) lign80
        IF (io > 0) THEN
           WRITE(*,*) 'Check input.  Something was wrong'
           STOP
           ! End of file
        ELSE IF (io < 0) THEN
           EXIT
        ELSE
           IF (lign80(1:4).eq.'ATOM') THEN
              READ(lign80,'(12X,A4,1X,A3,10X,3F8.3)') atname,attype,x(i+1),y(i+1),z(i+1)
              i=i+1
              IF ((atname.eq.' CA ').or.((attype.eq.lig).and..not.(atname(2:2).eq."H"))) THEN
                 inc(i)=1
              ELSE
                 inc(i)=0
              END IF
           END IF
        END IF
     END DO
     
     CLOSE(4569)
     
     OPEN(file="matrix.eigenfacs",form="FORMATTED",unit=4510)
     
     DO k=1,fnum
        
        WRITE(4510,'(1X,A6,I5,7X,A5,G12.4)') "VECTOR", k, "VALUE", eigenval(k) 
        WRITE(4510,'(A36)') " -----------------------------------"
        
        DO i=1,natom
           IF (inc(i).eq.1) THEN
              WRITE(4510,'(3(1X,G11.4))') 10*eigenvec(3*(i-1)+1,k),&
                   10*eigenvec(3*(i-1)+2,k),10*eigenvec(3*(i-1)+3,k)
           END IF
        END DO
        
     END DO
     
     CLOSE(4510)
     
     
     i=0
     DO j=1,natom
        IF (inc(j).eq.1) THEN
           i=i+1
           x(i)=x(j)
           y(i)=y(j)
           z(i)=z(j)
        END IF
     END DO
     i=0
     OPEN(file=pdborig,form="FORMATTED",status="OLD",unit=4599)
     OPEN(file="CAonly.pdb",form="FORMATTED",unit=4600)
     DO
        READ(4599,'(A)',IOSTAT=io) lign80
        IF (io > 0) THEN
           WRITE(*,*) 'Check input.  Something was wrong'
           STOP
           ! End of file
        ELSE IF (io < 0) THEN
           EXIT
        ELSE
           IF ((lign80(1:4).eq.'ATOM').or.(lign80(1:6).eq.'HETATM')) THEN
              READ(lign80,'(12X,A4,1X,A3)') atname,attype
              IF ((atname.eq.' CA ').or.(attype.eq.lig)) THEN
                 i=i+1
                 WRITE(4600,'(A30,3F8.3,A26)') lign80(1:30),x(i),y(i),z(i),lign80(55:80)
              ELSE
              END IF
           END IF
        END IF
     END DO
     
     
     CLOSE(4599)
     CLOSE(4600)
     
     DEALLOCATE(eigenvec)
     DEALLOCATE(eigenval)
     DEALLOCATE(inc)
     DEALLOCATE(x)
     DEALLOCATE(y)
     DEALLOCATE(z)
     
     
  END IF
  
END PROGRAM GroAMED


SUBROUTINE helptext(iunit)
  IMPLICIT NONE
  INTEGER :: iunit
  WRITE(iunit,'(A)')"                           G  R  O  A  M  E  D                           "
  WRITE(iunit,'(A)')"                               VERSION 1.0                               "
  WRITE(iunit,'(A)')"                                                                         "  
  WRITE(iunit,'(A)')"                               Written by:                               "
  WRITE(iunit,'(A)')"                      Tom Rodgers and David Burnell                      "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"This program converts the output from gromacs g_covar, gromacs g_nmeig,  "
  WRITE(iunit,'(A)')"and amber ptraj commands into files that are compatable with this        "
  WRITE(iunit,'(A)')"toolbox. For gromacs, the eigenvec.trr file needs to be converted into a "
  WRITE(iunit,'(A)')"readable format using gromacs:                                           "
  WRITE(iunit,'(A)')"       gmxdump -f eigenvec.trr > eigenvec.out                            "
  WRITE(iunit,'(A)')"This file can then be input using the -i flag.                           "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"The average pdb file produced by gromacs is also needed with the original"
  WRITE(iunit,'(A)')"pdb file that was used. The eigenval.xvg produced by gromacs is also     "
  WRITE(iunit,'(A)')"used.                                                                    "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"A matrix.eigenfacs file is outputted with a pdb file of the Ca atoms,    "
  WRITE(iunit,'(A)')"CAonly.pdb                                                               "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"     Usage:                                                              "
  WRITE(iunit,'(A)')"           GroAMED -i infile -mdpdb mdpdb -pdb originalpdb               "
  WRITE(iunit,'(A)')"                   [-eigen eigenfile] [-lig ligand] [-nofit] [-ca]       "
  WRITE(iunit,'(A)')"                   [-mass resmass]                                       "
  WRITE(iunit,'(A)')"                                                                         "
  WRITE(iunit,'(A)')"Option    Type      Value        Description                             "
  WRITE(iunit,'(A)')"------------------------------------------------------------             "
  WRITE(iunit,'(A)')"  -i      Input                  CA ouputfile                            "
  WRITE(iunit,'(A)')" -mdpdb   Input                  MD output pdbfile                       "
  WRITE(iunit,'(A)')" -pdb     Input                  Original pdb file                       "
  WRITE(iunit,'(A)')" -eigen   Input,Opt eigenval.xvg Generated eigenvalue file               "
  WRITE(iunit,'(A)')" -lig     Input,Opt  NONE        Name of any ligand included             "
  WRITE(iunit,'(A)')" -nofit    Opt                   Use if no fitting was used in gromacs   "
  WRITE(iunit,'(A)')"  -ca      Opt                   Use if only Ca atoms were fit in gromacs"
  WRITE(iunit,'(A)')" -mass    Input,Opt  99.466      Average residue mass for rescaling if   "
  WRITE(iunit,'(A)')"                                 -ca is used                             "
  WRITE(iunit,'(A)')"  -nm      Opt                   Use normal mode data instead of PCA data"
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

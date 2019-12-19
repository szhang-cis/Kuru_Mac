 SUBROUTINE pcte05 (flag,iload)
 !
 !  Write mesh information for TECPLOT for SOL3D elements
 !
 USE data_db
 IMPLICIT NONE
 LOGICAL, INTENT(IN) :: flag  !.FALSE. mesh only, .TRUE. results
 INTEGER(kind=4), INTENT(IN) :: iload

 TYPE( sol3d ), POINTER :: e
 INTEGER :: ngv,i,n,iset,l,iel,nnode,nlgv,ln(8)
 CHARACTER (len=15) :: etype
 CHARACTER (len=32) :: sname
 LOGICAL :: one                       !flag
 REAL(kind=4), POINTER :: gv(:,:)     !global variables
 REAL(kind=4), ALLOCATABLE :: lgv(:)  !local gauss variables

 INTERFACE
   INCLUDE 'tecgetgv.h'
 END INTERFACE

 ! tasks to do
 ! 1-PRINT mesh (flag=.FALSE.)
 !   a-print nodal coordinates (coord or coors if stage displements) of Truss sets only
 !   b-print connectivities of each truss set  using VARSHARELIST=([1,2,3]) (ndime)
 ! 2-PRINT results (flag=.TRUE.)
 !   variables are of three 'Groups'
 !         I-Global nodal variables: displacements, velocities, accelerations, pressure and temperatures (1:ngvar)
 !        II-Nodal Smoothed variables (truss_nvarn) for this type of elements (n1:n1-1+truss_nvarn), n1 must be determined
 !       III-Gaussian variables (truss_nvarg) for this type of elements (n2:n2-1+truss_nvarg), n2 must be determined
 !   a-print nodal displacements, velocities, etc. and smoothed variables at the first zone
 !   b-print Gaussian variables at each zone using VARLOCATION VARLOCATION=([?,?]=CELLCENTERED)
 !   c-it must be determined what to include  VARSHARELIST    (type I and II)
 !   d-it must be determined what to include  PASSIVEVARLIST  (thouse no included in type I, II and III)

 CALL tecgetgv(flag,gv,sol3d_nodes,sol3d_nps,ngv)  !get global variables
 IF( flag) CALL tecgenvarshls(4)   !Generate Variables Share List and Passive Variables list
 nlgv = t_ev(2,4)                  !number of Local Gauss Variables

 e => sol3d_head                 !point to first element set
 one = .TRUE.                    !true for first set of sol3d elemenents
 DO iset=1,sol3d_sets            !for each element set
   nnode = MIN(e%nnode,4)        !lower order nodes
   IF( e%nnode == 4)THEN
     etype = 'FETETRAHEDROM  '
   ELSE                          !Hexahedra or Prism
     etype = 'FEBRICK        '
   END IF

   l = LEN_TRIM(e%sname)
   sname(1:l+2) = '"'//e%sname(1:l)//'"'     !set name with added apostrophes
   WRITE(11,"('ZONE T=',a,' Nodes=',i6,' Elements=',i6,' ZONETYPE=',a15,    &
              ' DATAPACKING=BLOCK  SOLUTIONTIME=',e12.4 )") &   !,' STRANDID=',i4
               sname(1:l+2),sol3d_nps, e%nelem, etype, ttime    !, iload

   IF( flag )THEN !SOLUTION
     IF( t_npv > 0)WRITE(11,"(a)")TRIM(ps_var)  !print list of passive variables
     IF( nlgv > 0) WRITE(11,"(a)")TRIM(lc_var)  !print variables location
     IF( one )THEN                  !only print for first set
       DO i=1,ngv    !print global variables
         WRITE(11,"(10e15.6)")(gv(i,n),n=1,sol3d_nps)
         WRITE(11,"()")            !leave a line between data blocks
       END DO
       DO i=1,sol3d_nvarn !print smoothed variables
         WRITE(11,"(10e15.6)")(sol3d_vargs(i,n),n=1,sol3d_nps)
         WRITE(11,"()")            !leave a line between data blocks
       END DO

     ELSE
       WRITE(11,"(a)")TRIM(sh_var)  !print lis of common variables
     END IF
     ! print Gaussian Variables
     IF( nlgv > 0 )THEN
       ALLOCATE( lgv(e%nelem) )   !get memory for average element values
       DO i=1,nlgv                ! for each variabls
         DO iel=1,e%nelem         ! compute cell-center values
           lgv(iel) = SUM(e%elvar(i,:,iel))/e%ngaus
         END DO
         WRITE(11,"(10e15.6)") (lgv(iel),iel=1,e%nelem)  !print gauss values
         WRITE(11,"()")            !leave a line between data blocks
       END DO
       DEALLOCATE( lgv )
     END IF

   ELSE           !GRID
     IF( one )THEN               !only print for first set
       DO i=1,ngv
         WRITE(11,"(10e15.6)")(gv(i,n),n=1,sol3d_nps)
         WRITE(11,"()")            !leave a line between data blocks
       END DO
     ELSE
       WRITE(11,*)TRIM(sh_var)
     END IF
     DO iel = 1,e%nelem            !for each element in the set
       IF( e%nnode == 6 .OR. e%nnode == 15 )THEN   !6/15-node prism
         ln = (/ e%lnods(1:3,iel),e%lnods(3:6,iel),e%lnods(6,iel) /)
         WRITE(11,"(10i8)") sol3d_nodes(ln,1)  !print local numeration
       ELSE ! hexahedra
         WRITE(11,"(10i8)") sol3d_nodes(e%lnods(1:8,iel),1)  !print local numeration
       END IF
     END DO
   END IF
   one = .FALSE.
   e => e%next
 END DO

 RETURN

 END SUBROUTINE pcte05

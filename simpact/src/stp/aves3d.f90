SUBROUTINE aves3d( )
  !
  ! interpolate smoothed values for quadratic elements
  ! 3-D solid elements
  !
  USE data_db
  IMPLICIT NONE
  
  ! local variables
  TYPE (sol3d), POINTER :: eset
  INTEGER :: ielem,i,n,np,nn,nelem,nnode,nvert,lnods(20),ii
  INTEGER, POINTER :: nodes(:,:)
  REAL (kind=8), POINTER :: vargs(:,:),accpn(:)
  INTEGER (kind=4) :: id(2,12,3) = reshape((/                    &
       1,2, 2,3, 3,1, 1,4, 2,4, 3,4, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, & !10-nodes Tetrahedra
       1,2, 2,3, 3,1, 1,4, 2,5, 3,6, 4,5, 5,6, 6,4, 0,0, 0,0, 0,0, & !15-nodes Prism
       1,2, 2,3, 3,4, 4,1, 1,5, 2,6, 3,7, 4,8, 5,6, 6,7, 7,8, 8,5 /), & !20-nodes Hexahedra
       (/2,12,3/))
  
  
  vargs => sol3d_vargs         !pointers to use short names
  accpn => sol3d_accpn
  nodes => sol3d_nodes
  
  eset => sol3d_head           !point to first set
  
  DO
     IF( .NOT.ASSOCIATED(eset) )EXIT
     nelem = eset%nelem      !number of elements in the set
     nnode = eset%nnode      !number of nodes per element
     IF( nnode > 8 )THEN     !only for quadratic elements
        nvert = NINT(nnode/2.5) !
        SELECT CASE (nnode)
        CASE(10)    ! nvert = 4
           ii = 1
        CASE(15)    ! nvert = 6
           ii = 2
        CASE(20)    ! nvert = 8
           ii = 3
        END SELECT
        
        DO ielem=1,nelem                             !for each element
           lnods(1:nnode) = eset%lnods(:,ielem)       !element connectivities
           lnods(1:nnode) = nodes(lnods(1:nnode),1)   !local nodes
           DO i=nvert+1,nnode                     !for each mid side node
              n = lnods(i)                         !associated node
              IF( accpn(n) /= 0) CYCLE
              np = lnods(id(1,i-nvert,ii))                    !previous node
              nn = lnods(id(2,i-nvert,ii))                    !next  node
              vargs(:,n) = ( vargs(:,np) + vargs(:,nn) )/2d0  !average value
           END DO
        END DO
     END IF
     eset => eset%next
  END DO
  
  RETURN
END SUBROUTINE aves3d

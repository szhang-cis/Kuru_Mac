 SUBROUTINE pcgi05 ( )
 !
 !  Write mesh information for GiD for 3-D-Solids
 !
 USE data_db
 IMPLICIT NONE
 TYPE( sol3d ), POINTER :: e  !element set
 REAL(kind=8) :: a,b,c
 INTEGER :: iset,iel,n,k,g
 LOGICAL :: gauss
 CHARACTER (len=13) :: etype
 CHARACTER (len=32) :: sname
 CHARACTER (len=34) :: gpname
 REAL(kind=8), ALLOCATABLE :: p(:),w(:)

 a = 1d0/SQRT(3D0)      !TTT gauss point position
 b = 1d0/3d0            !center of a triangle
 c = 1d0/6d0            !

 gauss =  sol3d_nvarg > 0

 e => sol3d_head                 !point to first element set

 DO iset=1,sol3d_sets            !for each element set
   k = LEN_TRIM(e%sname)
   sname = '"'//e%sname(1:k)//'"'
   gpname = '"GP'//sname(2:32)
   k = k + 4
   CALL prtcoo(e%nnode,3,sname,sol3d_nodes(1,2),iset)        !write coordinates and headers

   DO iel = 1,e%nelem            !for each element in the set
     WRITE(11,"(22i8)") iel,(label(e%lnods(n,iel)),n=1,e%nnode),e%matno(iel) !e%set
   END DO
   IF( gauss )THEN
!-----------------------------------------------------------------------------
     SELECT CASE (e%nnode)
     CASE ( 4 )
       etype = 'Tetrahedra   '
     CASE ( 6,15 )
       etype = 'Prism        '
     CASE (8, 20)
       etype = 'Hexahedra    '
     END SELECT

     WRITE(13,"('GaussPoints  ',a,' Elemtype  ',a13,3x,a)") &
                gpname(1:k),etype,sname

     SELECT CASE (e%nnode)
     CASE ( 4 )   !Tetrahedra
       WRITE(13,"('  Number Of Gauss Points:',i3,/,        &
               &  '  Natural Coordinates: internal')")e%ngaus

     CASE ( 6,15 )  !linear or Quadratic Prism
       SELECT CASE (e%ngaus)
       CASE ( 1, 6 )
         WRITE(13,"('  Number Of Gauss Points:',i3,/,        &
                 &  '  Natural Coordinates: internal')") 1 !e%ngaus
       CASE (2)
         IF( given ) THEN
           WRITE(13,"('  Number Of Gauss Points:',i3,/,        &
                   &  '  Natural Coordinates: Given')")e%ngaus
           IF( e%etype == 16 )THEN          !Prism
             WRITE(13,"(3e20.12)") b, b,-a
             WRITE(13,"(3e20.12)") b, b, a
           ELSE                             !Sprism (extrapolated to faces)
             WRITE(13,"(3e20.12)") b, b,-1d0
             WRITE(13,"(3e20.12)") b, b, 1d0
           END IF
         ELSE
           ! extend to six
           WRITE(13,"('  Number Of Gauss Points:',i3,/,        &
                   &  '  Natural Coordinates: internal')")e%ngaus*3
         END IF
       END SELECT

     CASE (8, 20)

       SELECT CASE (e%ngaus)
       CASE ( 1, 8 )
         WRITE(13,"('  Number Of Gauss Points:',i3,/,        &
                 &  '  Natural Coordinates: internal')")e%ngaus

       CASE (2)
         IF( given ) THEN
           WRITE(13,"('  Number Of Gauss Points:',i3,/,        &
                   &  '  Natural Coordinates: Given')")e%ngaus
           WRITE(13,"(3e20.12)") 0.0, 0.0, -1.0
           WRITE(13,"(3e20.12)") 0.0, 0.0,  1.0
         ELSE
           ! extend to 8
           WRITE(13,"('  Number Of Gauss Points:',i3,/,        &
                   &  '  Natural Coordinates: internal')") 8
         END IF

       CASE ( 4 )
         IF( given ) THEN
           WRITE(13,"('  Number Of Gauss Points:',i3,/,        &
                   &  '  Natural Coordinates: Given')")e%ngaus
           WRITE(13,"(3e20.12)") a,-a,-a
           WRITE(13,"(3e20.12)")-a, a,-a
           WRITE(13,"(3e20.12)")-a,-a, a
           WRITE(13,"(3e20.12)") a, a, a
         ELSE
           ! extend to 8
           WRITE(13,"('  Number Of Gauss Points:',i3,/,        &
                   &  '  Natural Coordinates: internal')")e%ngaus*2
         END IF
       END SELECT
     END SELECT
     WRITE(13,"('End GaussPoints')")

   END IF
   e => e%next
 END DO

 RETURN

 END SUBROUTINE pcgi05

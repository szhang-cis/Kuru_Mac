 SUBROUTINE readcs ( )
 !
 !   read surface data to compute surface contact press
 !
 USE cont_db
 USE data_db, ONLY : ndime
 IMPLICIT NONE

 CHARACTER (len=30) :: sname
 LOGICAL :: spress,swrink,cpress
 INTEGER (kind=4) :: ncnod,nsegm
 TYPE (surf_db), POINTER :: surf

 INTEGER (kind=4) :: i,j

 nsurf = 0
 CALL ini_srf(shead,stail)

 DO
   READ(44) sname,spress,swrink,cpress
   IF( .NOT.spress .AND. .NOT.swrink ) EXIT
   READ(44) ncnod,nsegm
   CALL new_surf(surf)
   ALLOCATE( surf%lcnod(ncnod), surf%lcseg(ndime,nsegm) )
   READ(44) (surf%lcnod(i),i=1,ncnod)
   READ(44)((surf%lcseg(j,i),j=1,ndime),i=1,nsegm)
   surf%sname = sname
   surf%press = spress
   surf%cpress= cpress
   surf%wrink = swrink
   surf%ncnod = ncnod
   surf%nsegm = nsegm
   nsurf = nsurf+1
   CALL add_srf(surf,shead,stail)
   IF( spress )press = .TRUE.
   IF( swrink )wrink = .TRUE.
 END DO

 RETURN
 END SUBROUTINE readcs

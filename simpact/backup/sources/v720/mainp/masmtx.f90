SUBROUTINE masmtx ()

 ! *** calculates consistent (non-diagonal) mass matrix

  USE cms_db
  USE ctrl_db, ONLY: ndime, ndofn, npoin, neq, numct, mscal, ascal, lumped
  USE outp_db, ONLY: sumat,iwrit
  USE npo_db, ONLY :  emass,ymass,label,ifact,ifpre
  USE kinc_db, ONLY : ndumn,npsdf,nesdf,ftsdf,maxav,maxa
  USE ndp_db, ONLY : nndp
  USE sms_db, ONLY : selective_mass_scaling, fibre_mass

  IMPLICIT NONE

  INTEGER (kind=4) :: i,j,k,ipoin,chnode,dpos(ndofn)
  REAL(kind=8), ALLOCATABLE :: xcmas(:) !ndofn*(ndofn+1)/2)

  INTERFACE
    INCLUDE 'elemnt.h'
    INCLUDE 'contac.h'
    INCLUDE 'colsol.h'
    INCLUDE 'lump_mass.h'
    INCLUDE 'reduce_mass.h'
  END INTERFACE

  ! compute dependencies of rigid bodies for a correct assemble for non-diagonal matrices
  IF(.NOT.lumped .AND. ndumn > 0) CALL rigbdc (nndp)   !nodal dependencies of slave DOFs
  emass = 0d0 !initializes mass matrix
  ymass = 0d0

  IF (nconm > 0) THEN  !add concentrated masses
    IF(iwrit == 1) WRITE(lures,"(5x,'Concentrated Masses')",ERR=9999)
    ALLOCATE ( xcmas(ndofn*(ndofn+1)/2) )
    xcmas = 0d0
    k = 1
    DO j=1,ndofn
      dpos(j) = k
      k = k + ndofn - j + 1
    END DO

    DO i = 1,nconm
      ipoin = nodcms(i)
      ipoin = chnode(ipoin)
      xcmas(dpos) = cmass(1:ndofn,i)
      IF(iwrit == 1)WRITE(lures,"(i8,6e12.4)",ERR=9999) label(ipoin),cmass(1:ndofn,i)
      emass(1:ndofn,ipoin) = emass(1:ndofn,ipoin) + cmass(1:ndofn,i)
      IF(.NOT.lumped)CALL ensmat(ndofn,ifpre(1,ipoin),xcmas(1),ymass(1))
    END DO
  END IF

  sumat = 0d0

  CALL elemnt ('LUMASS')  !compute element masses
  IF( numct > 0) CALL contac ('LUMASS',0)  !contact surface mass
  IF( lumped ) THEN
    CALL elemnt ('RIGBDY')  !center or mass of rigid bodies
    IF( ndumn > 0) CALL rigbdc (nndp)   !nodal dependencies of slave DOFs
  END IF

  WRITE(lures, "(//'Total Mass of the System (Without Point Masses):', &
               & e15.7)",ERR=9999) sumat

  IF(iwrit == 1) THEN
    IF(ndime == 2)WRITE(lures,"(//,'  Lumped Mass Matrix' /,          &
                              &    ' Node',6X,'X',11X,'Y',11X,'a')",ERR=9999)
    IF(ndime == 3)WRITE(lures,"(//,'  Lumped Mass Matrix' /,' Node  ',  &
                & 6X,'X',11X,'Y',11X,'z',10x,'RX',10X,'RY',10X,'RZ')",ERR=9999)
    DO ipoin=1,npoin
      WRITE(lures,"(i7,8e12.4)",ERR=9999) label(ipoin),emass(1:ndofn,ipoin)
    END DO
  END IF
  IF (ascal) CALL scale_mass( mscal ) ! auto-scale mass matrix factor evaluation

  IF( ABS(mscal-1d0) > 1e-8)THEN
    DO ipoin=1,npoin
      emass(:,ipoin) = emass(:,ipoin)*mscal !scale mass
    END DO
  END IF

  IF( lumped )THEN
    CALL ensmal(ndofn*npoin,neq,ifpre(1,1),emass(1,1),ymass(1),npsdf(1),nesdf(1),ftsdf(1))
    IF( selective_mass_scaling )CALL fibre_mass( emass )  ! where ???
    DO i=1,neq
      ymass(i) = 1d0/ymass(i)
    END DO
  ELSE
    IF( ABS(mscal-1d0) > 1e-8)THEN
      DO i=1,neq
        ymass(i) = ymass(i)*mscal !scale mass
      END DO
    END IF
    print *,'reducing mass matrix size'
    CALL reduce_mass(ymass(:),maxav(:),neq,maxa)
    !DO i=1,neq
    !  j = maxav(i)
    !  k = maxav(i+1)-1
    !  WRITE(58,"(i7,(7e15.4))")i,ymass(j:k)
    !END DO
    print *,'Factoring mass matrix'
    CALL colsol(ymass(:),maxav(:),neq,1,58,0,0)
    !DO i=1,neq
    !  j = maxav(i)
    !  WRITE(58,"(i7,(7e15.4))")i,ymass(j)
    !END DO
    print *,'mass matrix factorized'
  END IF


RETURN
9999 CALL runen2('')
END SUBROUTINE masmtx

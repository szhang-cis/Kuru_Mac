SUBROUTINE lumass ()

  ! *** calculates lumped mass for each element & assembles global mass matrix

  USE cms_db
  USE ctrl_db, ONLY: ndime, ndofn, npoin, neq, numct, mscal, ascal
  USE outp_db, ONLY: sumat,iwrit
  USE npo_db, ONLY :  emass,ymass,label,ifact,ifpre
  USE kinc_db, ONLY : ndumn,npsdf,nesdf,ftsdf
  USE ndp_db, ONLY : nndp
  USE sms_db, ONLY : selective_mass_scaling, fibre_mass

  IMPLICIT NONE

  INTEGER (kind=4) :: i,ipoin,chnode
  REAL    (kind=8) :: xcmas(ndofn)

  INTERFACE
    INCLUDE 'elemnt.h'
    INCLUDE 'contac.h'
  END INTERFACE

  emass = 0d0 !initializes mass matrix

  IF (nconm > 0) THEN  !add concentrated masses
    IF(iwrit == 1) WRITE(lures,"(5x,'Concentrated Masses')",ERR=9999)

    DO i = 1,nconm
      ipoin = nodcms(i)
      ipoin = chnode(ipoin)
      xcmas(1:ndofn) = cmass(1:ndofn,i)
      IF(iwrit == 1)WRITE(lures,"(i8,6f10.3)",ERR=9999) label(ipoin),xcmas
      emass(1:ndofn,ipoin) = emass(1:ndofn,ipoin) + xcmas
    END DO
  END IF

  sumat = 0d0

  CALL elemnt ('LUMASS')  !compute element masses
  IF( numct > 0) CALL contac ('LUMASS',0)  !contact surface mass
  CALL elemnt ('RIGBDY')  !center or mass of rigid bodies
  IF( ndumn > 0) CALL rigbdc (nndp)

  WRITE(lures, "(//'Total Mass of the System (Without Point Masses):', &
               & e15.7)",ERR=9999) sumat

  IF(iwrit == 1) THEN
    IF(ndime == 2)WRITE(lures,"(//,'  Lumped Mass Matrix' /,          &
                              &    ' Node',6X,'X',11X,'Y',11X,'a')",ERR=9999)
    IF(ndime == 3)WRITE(lures,"(//,'  Lumped Mass Matrix' /,' Node',  &
                & 6X,'X',11X,'Y',11X,'z',10x,'RX',10X,'RY',10X,'RZ')",ERR=9999)
    DO ipoin=1,npoin
      WRITE(lures,"(i5,6e12.4)",ERR=9999) label(ipoin),emass(1:ndofn,ipoin)
    END DO
  END IF
  IF (ascal) CALL scale_mass( mscal ) ! auto-scale mass matrix factor evaluation
  emass = emass*mscal !scale mass

!           Assembles mass for active DOFs
  ymass = 0d0
  CALL ensmal(ndofn*npoin,neq,ifpre(1,1),emass(1,1),ymass(1),npsdf(1),nesdf(1),ftsdf(1))
  IF( selective_mass_scaling )CALL fibre_mass( emass )
  DO i=1,neq
    ymass(i) = 1d0/ymass(i)
  END DO
RETURN
9999 CALL runen2('')
END SUBROUTINE lumass

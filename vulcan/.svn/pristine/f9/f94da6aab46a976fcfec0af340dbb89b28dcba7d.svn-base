      SUBROUTINE CEDM20(DESIG,DMATX,DVOLU,EHIST,KGAUS,NDIME,
     .                  NNODE,NSTR1,PRESG,PROPS,SGTOT,SIGMA,STRAN,
     .                  ELCOD,RADUS,THICK)
C***********************************************************************
C
C****THIS ROUTINE EVALUATES TOTAL STRESSES FOR ELASTO/DAMAGED MATERIALS
C
C        NCRIT = 21        Isotropic damage
C        NCRIT = 22        Tensile   damage ( stress based )
C        NCRIT = 23        Model III ( AUSTRIAN PAPER )
C
C    Input variables:
C      
C      DESIG(NSTR1)        Working space
C      DMATX(NSTRS,NSTRS)  Elastic constitutive matrix
C      DVOLU(NGAUS)        Volume of the integration points
C      ELCOD(NDIME,NNODE)  Nodal coordinates
C      PRESG(NSTR1)        Previous total stresses
C      PROPS(NPROP)        Material properties
C      STRAN(NSTR1)        Current total strains
C      SIGMA(NSTR1)        Elastic stress predictor
C
C    Output variables:
C
C      SGTOT(NSTR1)        Current total stresses
C
C    Input/Output variables:
C
C      EHIST(1:6)          Internal variables for visco-elasticity 
C      EHIST(10)           Current "equivalent energy"        
C      EHIST(11)           Current damage internal variable
C      EHIST(12)           Current damage threshold
C      EHIST(13:15)        Characteristic lengths
C      EHIST(17)           Softening parameter
C      EHIST(19)           Time of onset of damage
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'auxl_om.f'
      DIMENSION DESIG(*), DMATX(NSTRS,*), EHIST(*), ELCOD(NDIME,*), 
     .          PRESG(*), PROPS(*),       SIGMA(*), SGTOT(*),
     .          STRAN(*), DVOLU(*)
C
C***COMPUTE TOTAL ELASTIC STRESS
C
      DO 10 ISTRE=1,NSTR1
   10 SIGMA(ISTRE)=0.0D00
      DO 20 ISTRE=1,NSTRS
      DO 20 JSTRE=1,NSTRS
   20 SIGMA(ISTRE)=SIGMA(ISTRE)+DMATX(ISTRE,JSTRE)*STRAN(JSTRE)
C
C   DEAL WITH DAMAGE
C
      IF(NCRIT.LT.24) THEN
        CALL SETDMG(DESIG,DMATX,DVOLU,ELCOD,NDIME,NNODE,
     .              NSTR1,PRESG,PROPS,SGTOT,SIGMA,STRAN,RADUS,
     .              EHIST(1),EHIST(10),EHIST(11),EHIST(12),
     .              EHIST(13),EHIST(17),EHIST(19),THICK)
      ENDIF
C
      RETURN
      END

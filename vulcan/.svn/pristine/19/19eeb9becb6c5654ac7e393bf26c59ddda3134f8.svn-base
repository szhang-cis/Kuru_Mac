      SUBROUTINE FIXCRK(BETAM,STNCU,STSVL,ALENG,NCRAK,NDIME,NSTR1,
     .                  PROPS,SGTOT,SIGMA,STRAN)
C***********************************************************************
C
C**** THIS ROUTINE DEALS WITH THE EXISTING ACTIVE CRACKS AND EVALUATES 
C     THE CURRENT STRESSES
C
C    Input variables:
C      
C      PROPS(NPROP)       Material properties
C      SIGMA(NSTR1)       'Trial' total stresses
C      STRAN(NSTR1)       Current total strains
C
C    Output variables:
C
C      SGTOT(NSTR1)       Current total stresses
C
C    Input/Output variables:
C
C      BETAM(NDIME,NDIME) Transformation matrix
C      STNCU(3)           Current strains in cracks
C      STSVL(3)           Secant stiffness moduli for cracks
C      ALENG(3)           Characteristic lengths for cracks        
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'auxl_om.f'
C
      DIMENSION ALENG(*), BETAM(NDIME,*), 
     .          PROPS(*), SGTOT(*),       SIGMA(*), 
     .          STNCU(*), STRAN(*),       STSVL(*)
      DIMENSION INDEX(3), TRANS(6,6)
      DATA MCRAK/3/
      DATA INDEX/1,2,4/
C
C***MAKE SGTOT=SIGMA
C
      DO 10 ISTRE=1,NSTR1
   10 SGTOT(ISTRE)=SIGMA(ISTRE)
C
C***IF NO CRACKS ARE PRESENT RETURN
C
      IF(NCRAK.EQ.0) RETURN
C
C***DEAL WITH EXISTING CRACKS
C
      YOUNG=PROPS( 2)
      FTULT=PROPS(32)
      GFVAL=PROPS(33)
      CONS1=       GFVAL/FTULT
      CONS2=0.5D00*FTULT/YOUNG
C
C   (A) LOOP ON NO. OF ACTIVE CRACKS
C
      DO 20 ICRAK=1,NCRAK 
C
C   (B) COMPUTE SOFTENING CONSTANT FOR THIS CRACK
C
       ALONG=ALENG(ICRAK)
       CONSS=CONS1/ALONG-CONS2
       IF(CONSS.LT.0.0) CONSS=0.1E-08
C
C   (C) ADJUST STRESS NORMAL TO THE CRACK
C
       ISTRE=INDEX(ICRAK)
C
       CALL STSEFF(CONSS,FTULT, STNCU(ICRAK),STSVL(ICRAK),
     .             STRAN(ISTRE),SGTOT(ISTRE),YOUNG)
C
   20 CONTINUE
C
C   (D) TRANSFORM STRESSES INTO GLOBAL SYSTEM
C
      CALL TRASIG(SGTOT,BETAM,NSTRE,NDIME,'L_TO_G')
C
      RETURN
      END

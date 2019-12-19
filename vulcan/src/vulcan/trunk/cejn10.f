      SUBROUTINE CEJN10(DMATX,ELCOD,EFFST,EPSTN,KGAUS,NDIME,NNODE,
     .                  NSTR1,PRESG,PROPS,RMAT1,SGTOT,SIGMA,STRAN,   
     .                  XLENG,conju,teno3,teno4,dstra)
C***********************************************************************
C
C****THIS ROUTINE EVALUATES TOTAL STRESSES FOR "SMEARED" JOINT MATERIALS
C                     NCRIT= 10
C
C    Input variables:
C      
C      DMATX(NSTRS,*)      Material matrix
C      ELCOD(NDIME,NNODE)  Coordinates of the element nodes
C      PRESG(NSTR1)        Previous total stresses
C      PROPS(NPROP)        Material properties
C      RMAT1(NDIME,NDIME)  Director cosines of material system
C      SIGMA(NSTR1)        'Trial' total stresses
C      STRAN(NSTR1)        Current total strains
C
C    Output variables:
C
C      SGTOT(NSTR1)        Current total stresses
C
C    Input/Output variables:
C
C      EFFST               Recorded plastic effective strain
C      EPSTN               Recorded plastic effective stress
C      XLENG               Recorded characteristic length
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'auxl_om.f'
C
      DIMENSION DMATX(NSTRS,*), ELCOD(NDIME,*), PRESG(*), PROPS(*), 
     .          RMAT1(NDIME,*), SGTOT(*),       SIGMA(*), STRAN(*),
     .          DSTRA(*)
C
      YOUNG=PROPS(2)
      OPENG=PROPS(32)
C
C   (A) TRANSFORM STRAINS INTO LOCAL SYSTEM 
C       index  = 2 deformacion  ;  ifl = 1  glob >> loca
C
      call tenvec(ndime,nstrs,STRAN,rmat1, 2 ,  1 )
      call tenvec(ndime,nstrs,DSTRA,rmat1, 2 ,  1 )
c      CALL TRASTR(STRAN,RMAT1,NSTR1,NDIME,'G_TO_L')
C
C   (B) TRANSFORM STRESSES INTO MATERIAL SYSTEM
C       index  = 1 tension  ;  ifl = 1  glob >> loca
C
      call tenvec(ndime,nstrs,SIGMA,rmat1, 1 ,  1 )
      call tenvec(ndime,nstrs,PRESG,rmat1, 1 ,  1 )
C      CALL TRASIG(SIGMA,RMAT1,NSTR1,NDIME,'G_TO_L')
C      CALL TRASIG(PRESG,RMAT1,NSTR1,NDIME,'G_TO_L')
C
C   (C) FIND CHARACTERISTIC LENGTH
C 
      IF(XLENG.EQ.0.0D00)
     .  CALL LENJNT(XLENG,ELCOD,RMAT1,NDIME,NNODE)
C
C   (D) MAKE SGTOT=SIGMA
C
      DO 10 ISTRE=1,NSTR1
   10 SGTOT(ISTRE)=SIGMA(ISTRE)
C
C   (E) ADJUST STRESS NORMAL TO THE JOINT
C
       CALL STSNML(XLENG,OPENG,STRAN(1),SGTOT(1),YOUNG,conju,sgtot(2),
     .             props(33),teno3,teno4,props(34),sgtot(4),STRAN(4),
     .             props(35),dstra(1))
C
C   (E) ADJUST SHEAR STRESSES
C
      CALL PLASJT(NSTR1,PROPS,SGTOT)
C
C   (F) TRANSFORM STRESSES INTO GLOBAL SYSTEM
C       index  = 1 tension  ;  ifl = 1  glob >> loca
C
      call tenvec(ndime,nstrs,SGTOT,rmat1, 1 ,  2 )
C      CALL TRASIG(SGTOT,RMAT1,NSTR1,NDIME,'L_TO_G')
C
      RETURN
      END

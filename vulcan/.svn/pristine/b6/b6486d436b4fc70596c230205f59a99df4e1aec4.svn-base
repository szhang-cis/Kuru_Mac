      SUBROUTINE STIF31(PROPS,ESTIF,TENOD,ELDIS,EHIST)
C*****************************************************************
C
C**** THIS ROUTINE COMPUTES THE STIFFNESS MATRIX FOR:
C                ( ELEMENT TYPE NO. 31 )
C
C*****************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION PROPS(NPROP),       ESTIF(NKOVA)
      DIMENSION TENOD(NNODL)
      DIMENSION ELDIS(NDOFN,NNODL), EHIST(NHIST,*)
      DIMENSION VECTN(3,3),         COEFN(3)
C
C**** ASSIGN PROPERTIES
C
      IF(NDOFN.EQ.2) THEN
C
C     VECTN(1,1)= NX    VECTN(1,2)= NY
C     VECTN(2,1)= TX    VECTN(2,2)= TY
C
C     COEFN(1)= CN (EHIST(2,1))     COEFN(2)= CT (EHIST(2,2))
C
       VECTN(1,1)=PROPS(1)
       VECTN(1,2)=PROPS(2)
       VECTN(2,1)=PROPS(3)
       VECTN(2,2)=PROPS(4)
       COEFN(1)=EHIST(2,1)
       COEFN(2)=EHIST(2,2)
C
       AAA=COEFN(1)*VECTN(1,1)*VECTN(1,1)+COEFN(2)*VECTN(2,1)*VECTN(2,1)
       BBB=COEFN(1)*VECTN(1,1)*VECTN(1,2)+COEFN(2)*VECTN(2,1)*VECTN(2,2)
       CCC=COEFN(1)*VECTN(1,2)*VECTN(1,2)+COEFN(2)*VECTN(2,2)*VECTN(2,2)
C
       ESTIF(1) = AAA
       ESTIF(2) = BBB
       ESTIF(3) =-AAA
       ESTIF(4) =-BBB
       ESTIF(9) = CCC
       ESTIF(10)=-BBB
       ESTIF(11)=-CCC
       ESTIF(16)= AAA
       ESTIF(17)= BBB
       ESTIF(22)= CCC
C
      ELSE
       CALL RUNEND('STIF31: 3D CASE NOT IMPLEMENTED YET')
      ENDIF
C
      RETURN
      END

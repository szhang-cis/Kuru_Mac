C
C------------------------------------------------------------------------
C
      SUBROUTINE CONTAC(A,IGRAD,ICODE,NDOFN,NDOFC,NPOIN,NTOTV,IFFIX,
     .                    FORZL,RKMAL,DISIL,DISPL,DISTL,DISIT,DISTO,
     .                    DISPR,IPOIN,REFRL,DGAPL,IITER,KCODE,
     .                                                           IND)
C***********************************************************************
C
C****ESTA RUTINA TRATA EL PROBLEMA DE CONTACTO Y DESPEGUE
C
C                IND =1    :NUEVO DESPEGUE
C                IND =2    :SIGE EL DESPEGUE
C                IND =3    :NUEVO CONTACTO
C                IND =4    :SIGUE EL CONTACTO
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IFFIX(*),FORZL(3),RKMAL(3,3),DISIL(3),DISPL(3),DISTL(3),
     .          DISIT(*),DISTO(*),DISPR(*),A(3,3),REFRL(3,3),DGAPL(3)
C
       DO ILOCA=1,NDOFN
        DGAPL(ILOCA)  =0.0
       ENDDO
C
       IF(KCODE.LT.0)RETURN
C
c       IF(IITER.GT.1)RETURN
C
       IF(IND.EQ.1)THEN
         ICODE=ICODE+1
         IFFIX(IGRAD)=0
         REFRL(1,3)=0.0
         FORZL(1)  =0.0
       ENDIF
C
       IF(IND.EQ.2)THEN
         IFFIX(IGRAD)=0
         REFRL(1,3)=0.0
         FORZL(1)  =0.0
       ENDIF
C
       IF(IND.EQ.3)THEN
        ICODE=ICODE+1
        IFFIX(IGRAD)=1
        DGAPL(1)  =-DISTL(1)
        REFRL(1,3)=0.0
c        FORZL(1)  =0.0
       ENDIF
C
       IF(IND.EQ.4)THEN
        IFFIX(IGRAD)=1
       ENDIF
C
       RETURN
       END

      SUBROUTINE TRSKEW(NPOIN,NTOTV,NDOFN,DISIT,INFRI,COFRI,INDDD,
     .                  NSKEW)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS A VECTOR TRANSFORMATION (TRANSFORM DISIT)
C
C     INDDD=1   GLOBAL > LOCAL
C     INDDD=2   LOCAL > GLOBAL
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION DISIT(NTOTV),
     .          INFRI(NPOIN), COFRI(NSKEW,NDOFN,*)         ! ndime=ndofn
C
      DIMENSION DISIL(3),     A(3,3)
C
C**** LOOP OVER NPOIN
C
      DO IPOIN=1,NPOIN
       IOEF3=INFRI(IPOIN)
       IF(IOEF3.NE.0) THEN
C
C**** ASSIGN COFRI TO MATRIX A
C
        DO IDOFN=1,NDOFN
         DO JDOFN=1,NDOFN
          A(IDOFN,JDOFN)=COFRI(IOEF3,IDOFN,JDOFN)
         ENDDO
        ENDDO
C
        IF(INDDD.EQ.1)THEN                          ! global ----> local
         DO ILOCA=1,NDOFN
          ITOTV=(IPOIN-1)*NDOFN+ILOCA
          DISIL(ILOCA)=0.0
C
          DO JLOCA=1,NDOFN
           JTOTV=(IPOIN-1)*NDOFN+JLOCA
           DISIL(ILOCA)=DISIL(ILOCA)+A(ILOCA,JLOCA)*DISIT(JTOTV)
          ENDDO
         ENDDO
C
C**** LOAD DISIL TO DISIT
C
         DO ILOCA=1,NDOFN
          ITOTV=(IPOIN-1)*NDOFN+ILOCA
          DISIT(ITOTV)=DISIL(ILOCA)
         ENDDO
C
        ELSE                                        ! local ----> global
C
         DO ILOCA=1,NDOFN
          ITOTV=(IPOIN-1)*NDOFN+ILOCA
          DISIL(ILOCA)=0.0
C
          DO JLOCA=1,NDOFN
           JTOTV=(IPOIN-1)*NDOFN+JLOCA
           DISIL(ILOCA)=DISIL(ILOCA)+A(JLOCA,ILOCA)*DISIT(JTOTV)
          ENDDO
         ENDDO
C
C**** LOAD DISIL TO DISIT
C
         DO ILOCA=1,NDOFN
          ITOTV=(IPOIN-1)*NDOFN+ILOCA
          DISIT(ITOTV)=DISIL(ILOCA)
         ENDDO
C
        ENDIF               ! INDDD.EQ.1
       ENDIF                ! IOEF3.NE.0
      ENDDO                 ! IPOIN=1,NPOIN
C
      RETURN
      END

      SUBROUTINE GAUCEK(DVOLU,ITYPE,ITASK,NGAUL,NNODL,ICEKE)
C***********************************************************************
C
C**** THIS ROUTINE CHECKS GAUSSIAN VOLUME FOR AXISYMMETRIC PROBLEMS FOR 
C     OUTPUT OPERATIONS.
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION DVOLU(*)
C
      ICEKE=0
      TOLVO=1.0D-08
C
C**** SOLID ELEMENTS
C
      IF(ITYPE.EQ.5.OR.ITYPE.EQ.30) THEN
       IF((ITASK.GE.12.AND.ITASK.LE.15).OR.(ITASK.EQ.17).OR.
     .                                     (ITASK.EQ.19)) THEN
        DVOLL=0.0D0                                       ! total volume
        DO IGAUL=1,NGAUL
         DVOLL=DVOLL+DVOLU(IGAUL)
        ENDDO
C
        IAVER=0
        DO IGAUL=1,NGAUL
         IF(DVOLU(IGAUL).LT.(TOLVO*DVOLL/NGAUL)) IAVER=1
        ENDDO
C
        IF(IAVER.EQ.1) THEN
         IF(NGAUL.EQ.3) THEN              ! nnodl=3,6,7 - simplification
          DO IGAUL=1,NGAUL
           IF(DVOLU(IGAUL).LT.(TOLVO*DVOLL/NGAUL)) THEN
            DVOLU(IGAUL)=DVOLL*1.0D0/3.0D0
           ELSE
            DVOLU(IGAUL)=DVOLL*1.0D0/3.0D0
           ENDIF
          ENDDO
         ENDIF
         IF(NGAUL.EQ.4) THEN
          IF(NNODL.EQ.4.OR.NNODL.EQ.8.OR.NNODL.EQ.9) THEN
           DO IGAUL=1,NGAUL
            IF(DVOLU(IGAUL).LT.(TOLVO*DVOLL/NGAUL)) THEN
             DVOLU(IGAUL)=DVOLL/8.0D0
            ELSE
             DVOLU(IGAUL)=DVOLL*3.0D0/8.0D0
            ENDIF
           ENDDO
          ENDIF
          IF(NNODL.EQ.6.OR.NNODL.EQ.7) THEN
           ICEKE=1
          ENDIF
         ENDIF
         IF(NGAUL.EQ.6) THEN                                 ! nnodl=6,7
          ICEKE=2
         ENDIF
         IF(NGAUL.EQ.9) THEN                                 ! nnodl=8,9
          DO IGAUL=1,NGAUL
           IF(DVOLU(IGAUL).LT.(TOLVO*DVOLL/NGAUL)) THEN
            DVOLU(IGAUL)=DVOLL/27.0D0
           ENDIF
          ENDDO
         ENDIF
        ENDIF
       ENDIF
      ENDIF
C
C**** BOUNDARY ELEMENTS
C
      IF(ITYPE.EQ.32.OR.ITYPE.EQ.33.OR.
     .   ITYPE.EQ.101.OR.ITYPE.EQ.104) THEN
       DVOLL=0.0D0                                        ! total volume
       DO IGAUL=1,NGAUL
        DVOLL=DVOLL+DVOLU(IGAUL)
       ENDDO
C
       IAVER=0
       DO IGAUL=1,NGAUL
        IF(DVOLU(IGAUL).LT.(TOLVO*DVOLL/NGAUL)) IAVER=1
       ENDDO
C
       IF(IAVER.EQ.1) THEN
        IF(NGAUL.EQ.2) THEN
         DO IGAUL=1,NGAUL
          IF(DVOLU(IGAUL).LT.(TOLVO*DVOLL/NGAUL)) THEN
           DVOLU(IGAUL)=DVOLL/4.0D0
          ELSE
           DVOLU(IGAUL)=DVOLL*3.0D0/4.0D0
          ENDIF
         ENDDO
        ENDIF
        IF(NGAUL.EQ.3) THEN
         ICEKE=3
        ENDIF
       ENDIF
      ENDIF
C
      RETURN
      END

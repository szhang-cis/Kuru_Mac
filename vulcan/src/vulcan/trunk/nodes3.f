      SUBROUTINE NODES3(COORD,LNODS,NDIME,IELEM,NNODE)
C***********************************************************************
C
C**** THIS ROUTINE INTERPOLATES THE MID-SIDE NODES OF STRAIGHT
C     SIDES OF ELEMENTS FOR 3D
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION COORD(NDIME,*), LNODS(*)
      DIMENSION DERIV(3,20),  ELCOD(3,20), SHAPE(20)
C
      IF(NNODE.LT.10) RETURN
C
      IF(NNODE.EQ.10)                NNOD1=5
      IF(NNODE.EQ.20.OR.NNODE.EQ.27) NNOD1=9
C
C**** COMPUTE THE NUMBER OF THE MIDSIDE NODE
C
      DO NODEB=NNOD1,NNODE
        NODMD=LNODS(NODEB)
C
C**** CHECK IF THE COORDINATES OF THE MIDSIDE NODE HAVE BEEN INPUT
C
        TOTAL=DABS(COORD(1,NODMD))+DABS(COORD(2,NODMD))
     .       +DABS(COORD(3,NODMD))
        IF(TOTAL.LE.0.0) THEN
C
C**** COMPUTE THE NUMBERS OF THE CORNER NODES OF THE CURRENT EDGE
C
          IF(NODEB.LE.20) THEN
            IF(NNOD1.EQ.5) THEN
              IF(NODEB.LE.6) THEN
                NODEA=NODEB-NNOD1+1
                NODEC=NODEA+1
                IF(NODEC.GT.3) NODEC=1
              ELSE
                NODEA=NODEB-NNOD1-2
                NODEC=4
              ENDIF
            ENDIF
            IF(NNOD1.EQ.9) THEN
              IF(NODEB.LE.12) THEN
                NODEA=NODEB-NNOD1+1
                NODEC=NODEA+1
                IF(NODEC.GT.4) NODEC=1
              ELSE IF(NODEB.LE.16) THEN
                NODEA=NODEB-NNOD1-3
                NODEC=NODEA+4
              ELSE IF(NODEB.LE.20) THEN
                NODEA=NODEB-NNOD1-3
                NODEC=NODEA+1
                IF(NODEC.GT.8) NODEC=5
              ENDIF
            ENDIF
            NODST=LNODS(NODEA)
            NODFN=LNODS(NODEC)
C
C**** INTERPOLATE BY A STRAIGHT LINE
C  
            COORD(1,NODMD)=(COORD(1,NODST)+COORD(1,NODFN))*0.5
            COORD(2,NODMD)=(COORD(2,NODST)+COORD(2,NODFN))*0.5
            COORD(3,NODMD)=(COORD(3,NODST)+COORD(3,NODFN))*0.5
C
          ELSE
            IF(NODEB.EQ.21) THEN
              EXISP= 0.
              ETASP= 0.
              EZETA=-1.
            ENDIF
            IF(NODEB.EQ.22) THEN
              EXISP= 0.
              ETASP=-1.
              EZETA= 0.
            ENDIF
            IF(NODEB.EQ.23) THEN
              EXISP= 1.
              ETASP= 0.
              EZETA= 0.
            ENDIF
            IF(NODEB.EQ.24) THEN
              EXISP= 0.
              ETASP= 1.
              EZETA= 0.
            ENDIF
            IF(NODEB.EQ.25) THEN
              EXISP=-1.
              ETASP= 0.
              EZETA= 0.
            ENDIF
            IF(NODEB.EQ.26) THEN
              EXISP= 0.
              ETASP= 0.
              EZETA= 1.
            ENDIF
            IF(NODEB.EQ.27) THEN
              EXISP= 0.
              ETASP= 0.
              EZETA= 0.
            ENDIF
C
C**** EVALUATE THE COORDINATES OF THE ELEMENT NODALS POINT
C
            DO INODE=1,20
              KLOCA=LNODS(INODE)
              DO IDIME=1,NDIME
                ELCOD(IDIME,INODE)=COORD(IDIME,KLOCA)
              ENDDO
            ENDDO
C
C**** COMPUTE THE SHAPE FUNCTION
C
            CALL SHAFUN(DERIV,EXISP,ETASP,EZETA,NDIME,   20,    0,    0,
     .                  SHAPE)
C
C**** COMPUTE THE GLOBAL COORDINATES
C
            KLOCA=LNODS(NODEB)
            DO IDIME=1,NDIME
              COORD(IDIME,KLOCA)=0.0
              DO INODE=1,NNODE
                COORD(IDIME,KLOCA)=COORD(IDIME,KLOCA)+SHAPE(INODE)*
     .                             ELCOD(IDIME,INODE)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDDO
C
      RETURN
      END

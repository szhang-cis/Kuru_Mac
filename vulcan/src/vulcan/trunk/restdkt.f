      SUBROUTINE RESTDKT(NUNIT,IFLAG,MATRI,LENGH,LENRC,NRCUR)
C***********************************************************************
C
C**** THIS ROUTINE WRITES (IFLAG=1) OR READS (IFLAG=2)
C     ONE MATRIX ('MATRI' WHICH LENGTH IS 'LENGH') TO OR FROM
C     'NRCUR'   CURRENT RECORD POINTER OF 'NUNIT' FILE
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
      DIMENSION MATRI(LENGH)
C
C**** LENGTH OF RECORD IN WORDS
C
      GO TO (101,102,102,101,101,101,101,101), NMACHI
C
  101 CONTINUE
      LENWD=LENRC/4
      GO TO 1000
C
  102 CONTINUE
      LENWD=LENRC
      GO TO 1000
C
 1000 CONTINUE
C
      NREST=LENGH-(LENGH/LENWD)*LENWD       ! # of resting words
      NFULL=LENGH-NREST                     ! # of full records
C
      IF(IFLAG.EQ.2) THEN     ! READ FROM DISK
       DO ILENG=1,NFULL,LENWD
        CALL READDK(NUNIT,MATRI(ILENG),LENWD,NRCUR)    ! full records
        NRCUR=NRCUR+1
       ENDDO
       IF(NREST.GT.0) THEN
        CALL READDK(NUNIT,MATRI(NFULL+1),NREST,NRCUR)  ! resting words
        NRCUR=NRCUR+1
       ENDIF
      ELSE                    ! WRITE TO DISK
       DO ILENG=1,NFULL,LENWD
        CALL WRITDK(NUNIT,MATRI(ILENG),LENWD,NRCUR)    ! full records
        NRCUR=NRCUR+1
       ENDDO
       IF(NREST.GT.0) THEN
        CALL WRITDK(NUNIT,MATRI(NFULL+1),NREST,NRCUR)  ! resting words
        NRCUR=NRCUR+1
       ENDIF
      ENDIF
C
      RETURN
      END

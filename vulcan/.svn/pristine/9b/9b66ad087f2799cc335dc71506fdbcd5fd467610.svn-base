      SUBROUTINE RSSETP
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES  THE RECORD ADDRES FOR DIFFERENT ARRAYS
C     STORED IN THE PROCESSING AREA
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
C
C**** SET IDATP POINTERS
C
      DO INDEX=1,12
       IDATP(INDEX,1)=0
       IDATP(INDEX,2)=0
      ENDDO
C
C**** EVALUATE LENGTH OF EACH ARRAY IN DOUBLE PRECISION WORDS
C
      NMEMA6M=1-NMEMO6M
      NMEMA7M=1-NMEMO7M
C
      IDATP( 1,1)=NDATA                                   ! ELDAT
      IDATP( 2,1)=NPREV                                   ! ELPRE
      IDATP( 3,1)=NSTAT                                   ! ELVAR
      IDATP( 4,1)=NPREV                                   ! ELPRE
      IDATP( 5,1)=NSTAT                                   ! ELVAR
      IDATP( 6,1)=NEVAC*NEVAC*NMEMA7M                     ! CSTIF
      IDATP( 7,1)=NKOVA*NMEMA6M                           ! ESTIF
      IF(KDYNA.NE.0) IDATP( 8,1)=NKOVA*NMEMA6M            ! WSTIF
      IF(KPORE.EQ.2) THEN
       IDATP( 9,1)=NKOND                                  ! PSTIF
       IDATP(10,1)=NKOND                                  ! QSTIF
       IDATP(11,1)=NEVAB*NNODE                            ! HSTIF
      ENDIF
      IDATP(12,1)=NTOTV                                   ! DISTO
C
      GO TO (101,102,102,101,101,101,101,101), NMACHIM
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
C**** EVALUATE LENGTH OF EACH ARRAY IN RECORDS
C
      DO INDEX=1,12
       IF(IDATP(INDEX,1).NE.0)
     .  IDATP(INDEX,2)=(2*IDATP(INDEX,1)-1)/LENWD+1
      ENDDO
C
C**** INITIALIZE POINTERS FOR POSITION IN DOUBLE PRECISION WORDS
C
      IDATP( 1,3)=0
      DO INDEX=2,7
       IDATP(INDEX,3)=IDATP(INDEX-1,3)+IDATP(INDEX-1,1)
      ENDDO
C
      IF((KPROB.EQ.1).AND.(KPORE.NE.2.AND.KDYNA.NE.1)) THEN
       IDATP(8,3)=IDATP(7,3)
       IDATP(7,3)=IDATP(6,3)
      ELSE
       IDATP(8,3)=IDATP(7,3)+IDATP(7,1)
      ENDIF
C
      IDATP(9,3)=IDATP(8,3)+IDATP(8,1)
C
      IDATP(10,3)=IDATP( 9,3)+IDATP( 9,1)
      IDATP(11,3)=IDATP(10,3)+IDATP(10,1)
      IDATP(12,3)=IDATP(11,3)+IDATP(11,1)
C
C**** INITIALIZE POINTERS FOR POSITION IN RECORDS
C
      IDATP( 1,4)=0
      DO INDEX=2,7
       IDATP(INDEX,4)=IDATP(INDEX-1,4)+IDATP(INDEX-1,2)
      ENDDO
C
      IF((KPROB.EQ.1).AND.(KPORE.NE.2.AND.KDYNA.NE.1)) THEN
       IDATP(8,4)=IDATP(7,4)
       IDATP(7,4)=IDATP(6,4)
      ELSE
       IDATP(8,4)=IDATP(7,4)+IDATP(7,2)
      ENDIF
C
      IDATP(9,4)=IDATP(8,4)+IDATP(8,2)
C
      IDATP(10,4)=IDATP( 9,4)+IDATP( 9,2)
      IDATP(11,4)=IDATP(10,4)+IDATP(10,2)
      IDATP(12,4)=IDATP(11,4)+IDATP(11,2)
C
C**** LAST RECORD OF THE BLOCK
C
      NLENP=IDATP(12,4)
C
C**** LAST RECORD OF THE FILE
C
      NRECP=NLENP*NELEM
C
C**** INITIALIZE FLAGS FOR ARRAYS ALLOCATION
C
      DO INDEX=1,11
       IDATP(INDEX,5)=-1
      ENDDO
      IDATP(12,5)=0
C
      RETURN
      END

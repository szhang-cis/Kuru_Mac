      SUBROUTINE OUTGAP(MATNO,PROEL,TGAPS,PREAS,LNODS,ACCPN)
C***********************************************************************
C
C**** THIS ROUTINE WRITES REACTIONS TO OUTPUT FILE
C
C     KPRI11=0 DO NOT WRITE NORMAL GAP
C            1 WRITE NORMAL GAP
C     KFEMV=0  DO NOT WRITE NORMAL GAP TO POSTPROCESSOR FILE ! NOT YET
C           1  WRITE NORMAL GAP TO POSTPROCESSOR FILE        ! NOT YET
C
C     Note: PREAS & TGAPS are printed even for uncoupled analysis
C           (see addpri.f, para_om.f, forcin.f, frin04*.f & frin32.f)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inpo_om.f'
C
      DIMENSION MATNO(NELEM), PROEL(NPREL,NGRUP), TGAPS(NPOIN),
     .          PREAS(NPOIN), LNODS(NNODE,NELEM)
      DIMENSION ACCPN(NPOIN)
C
C**** CONTROLS THE EXISTENCE OF GAP ELEMENTS
C
      IF(KGAPC.EQ.0) RETURN
C
      IF(KPRI11.EQ.0) RETURN              ! OUTPUT NORMAL GAP & PRESSURE
C
C**** COMPUTES NUMBER OF POINTS WITHOUT CONTACT (AL METHOD)
C
      NPOI1=NPOIN-NPOIC
C
C**** INITIALIZATION
C
      FACTA=1.D+00
      DO IPOIN=1,NPOIN
       ACCPN(IPOIN)=0.0D+00
      ENDDO
C
C**** LOOP OVER THE ELEMENTS (ONLY CONSIDERS CONTACT ELEMENTS)
C
      DO 100 IELEM=1,NELEM
       LGRUP=MATNO(IELEM)
       LMATS=INT(PROEL(1,LGRUP))
       LTYPE=INT(PROEL(5,LGRUP))
       NNODL=INT(PROEL(2,LGRUP))
C
       IF(LTYPE.EQ.2) GO TO 100       ! SKIP REINFORCING ELEMENTS
       IF(LTYPE.EQ.3) GO TO 100       ! SKIP LINK ELEMENTS
       IF(LTYPE.EQ.4) THEN
        DO INODE=1,NNODL
         LPOIN=LNODS(INODE,IELEM)
         IF(LPOIN.NE.0) ACCPN(LPOIN)=1.0D0
        ENDDO
       ENDIF
       IF(LTYPE.EQ.30) GO TO 100      ! SKIP SOLID ELEMENTS
       IF(LTYPE.EQ.32) THEN
        DO INODE=1,NNODL
         LPOIN=LNODS(INODE,IELEM)
         IF(LPOIN.NE.0) ACCPN(LPOIN)=ACCPN(LPOIN)+1.0D0
        ENDDO
       ENDIF
       IF(LTYPE.EQ.33) GO TO 100      ! SKIP LINK ELEMENTS
C
  100 CONTINUE
C
      DO IPOIN=1,NPOIN
       IF(ACCPN(IPOIN).LT.1.0D0) ACCPN(IPOIN)=1.0D+00
      ENDDO
C
C**** COMPUTES AVERAGE NORMAL GAP & PRESSURE
C
      DO IPOIN=1,NPOI1
       TGAPS(IPOIN)=TGAPS(IPOIN)/ACCPN(IPOIN)
       PREAS(IPOIN)=PREAS(IPOIN)/ACCPN(IPOIN)
      ENDDO
C
C**** WRITES NODAL NORMAL GAP & PRESSURE
C
      IF(IPRCO.GT.0) THEN
       WRITE(LURES,900)
       WRITE(LURES,960)
       DO IPOIN=1,NPOIN
        DO IELEM=1,NELEM
         LGRUP=MATNO(IELEM)
         LTYPE=INT(PROEL(5,LGRUP))
         IF(LTYPE.EQ.4.OR.LTYPE.EQ.32) THEN
          NNODL=INT(PROEL(2,LGRUP))
          DO INODL=1,NNODL
           IF(LNODS(INODL,IELEM).EQ.IPOIN) THEN
            WRITE(LURES,800) IPOIN,TGAPS(IPOIN),PREAS(IPOIN)
            GO TO 200
           ENDIF
          ENDDO                                  ! inodl=1,nnodl
         ENDIF                                   ! ltype.eq.4 or 32
        ENDDO                                    ! ielem=1,nelem
  200   CONTINUE
       ENDDO                                     ! ipoin=1,npoin
      ENDIF                                      ! iprco.gt.0
C
C**** WRITES TO POSTPROCESSOR FILE
C
      isie=1                           ! it should be implemented
      if(isie.eq.0) then
      IF(KFEMV.EQ.1)
     . WRITE(LUPOS) (SNGL(TGAPS(IPOIN)),IPOIN=1,NPOIN)
      endif
      RETURN
C
  800 FORMAT(5X,I5,E15.6,2X,E15.6)
  900 FORMAT(/,6X,'NORMAL GAP & PRESSURE')
  960 FORMAT(6X,'NODE',3X,'NORMAL GAP       NORMAL PRESSURE')
C
      END

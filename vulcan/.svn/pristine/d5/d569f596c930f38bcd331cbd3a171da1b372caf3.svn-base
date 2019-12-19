      SUBROUTINE DAINVA(NCRIT,PROPS,NTYPE,STRES,ANGFI,CAPAP,
     .                  PMEAN,VARJ2,VARJ3,DEVIA,STEFF,SMAXS,THETA,
     .                  DFTEQ,DFAFI,YIELD,NSTR1,ALFAT,BETAT,
     .                  GAMAT,RETEN,A1COE,A2COE,A3COE,BACKS,CCERO,
     .                  POROS,IDAMG,DAMAG,TAUAM,TAUAP,
     .                  CCEROM,CCEROP,GAMMAM,GAMMAP,LARGE)
C***********************************************************************
C
C**** THIS SUBROUTINE EVALUATES THE STRESS INVARIANTS AND THE CURRENT
C     VALUE OF THE YIELD FUNCTION
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION STRES(*), PROPS(*), DEVIA(6), BACKS(*), DBACK(6)
      DIMENSION STRESM(6), STRESP(6), STRSP(3), CGAMAM(6,6), CGAMAP(6,6)
C
      IF(LARGE.NE.0) RETURN      ! not implemented yet for large strains
C
      TWOPI=6.283185307179586
      ZM=1.0E-15
C
      IF(IDAMG.EQ.0) RETURN
C
      IF(IDAMG.LE.20) THEN   ! damage models governed by pl. or viscopl.
       RETURN
      ELSE                   ! damage models governed by a damage crit.
C
       GO TO (131,132) (IDAMG-20)
C
  131  CONTINUE              ! concrete damage model
       DO ISTR1=1,NSTR1
        STRESM(ISTR1)=0.0D0
        STRESP(ISTR1)=0.0D0
        DO JSTR1=1,NSTR1
         CGAMAM(ISTR1,JSTR1)=0.0D0
        ENDDO
       ENDDO
       CGAMAM(1,1)=1.0D0
       CGAMAP(1,1)=1.0D0
       IF(NTYPE.NE.5) THEN
        CGAMAM(2,2)=1.0D0
        CGAMAM(1,2)=-GAMMAM
        CGAMAM(2,1)=-GAMMAM
        CGAMAP(2,2)=1.0D0
        CGAMAP(1,2)=-GAMMAP
        CGAMAP(2,1)=-GAMMAP
        IF(NTYPE.NE.1) THEN
         CGAMAM(4,4)=1.0D0
         CGAMAM(1,4)=-GAMMAM
         CGAMAM(2,4)=-GAMMAM
         CGAMAM(4,1)=-GAMMAM
         CGAMAM(4,2)=-GAMMAM
         CGAMAP(4,4)=1.0D0
         CGAMAP(1,4)=-GAMMAP
         CGAMAP(2,4)=-GAMMAP
         CGAMAP(4,1)=-GAMMAP
         CGAMAP(4,2)=-GAMMAP
        ENDIF
       ENDIF
C
       IF(NTYPE.LE.3) THEN
        CALL PRIDIRN(NTYPE,STRES,STRSP)        ! eigenvalues not ordered
        DO I=1,2
        ANGLE=STRSP(3)
        IF(I.EQ.2) ANGLE=STRSP(3)+TWOPI/4.0D00
         CA=DCOS(ANGLE)
         SA=DSIN(ANGLE)
         IF(STRSP(I).LT.0.0) THEN
          STRESM(1)=STRESM(1)+STRSP(I)*CA*CA
          STRESM(2)=STRESM(2)+STRSP(I)*SA*SA
          STRESM(3)=STRESM(3)+STRSP(I)*SA*CA
         ELSE
          STRESP(1)=STRESP(1)+STRSP(I)*CA*CA
          STRESP(2)=STRESP(2)+STRSP(I)*SA*SA
          STRESP(3)=STRESP(3)+STRSP(I)*SA*CA
         ENDIF
        ENDDO
        IF(NTYPE.NE.1) THEN
         IF(STRES(4).LT.0.0) THEN
          STRESM(4)=STRES(4)
         ELSE
          STRESP(4)=STRES(4)
         ENDIF
        ENDIF
C
C       DO ISTR1=1,NSTR1                                  ! verification
C        IF(DABS(STRESM(ISTR1)+STRESP(ISTR1)-STRES(ISTR1)).GT.ZM)
C    .    CALL RUNEND('DAINVA: ERROR IN THE COMPUTATION OF sigma-+')
C       ENDDO
C
        TAUAM=0.0
        TAUAP=0.0
        DO ISTR1=1,NSTR1
         DO JSTR1=1,NSTR1
          TAUAM=TAUAM+STRESM(ISTR1)*CGAMAM(ISTR1,JSTR1)*STRES(JSTR1)
          TAUAP=TAUAP+STRESP(ISTR1)*CGAMAP(ISTR1,JSTR1)*STRES(JSTR1)
         ENDDO
        ENDDO
C       IF((TAUAM+TAUAP).LT.0.0)                          ! verification
C    .   CALL RUNEND('WARNING: TAUAM+TAUAP LT 0.0')
        IF(TAUAM.LT.0.0) TAUAM=-TAUAM
        IF(TAUAP.LT.0.0) TAUAP=-TAUAP
        TAUAM=DSQRT(TAUAM)/CCEROM
        TAUAP=DSQRT(TAUAP)/CCEROP
       ENDIF
       IF(NTYPE.EQ.4) THEN
        call runend('DAINVA: NTYPE=4 not implemented')
        CALL PRINSI(STRES,STRSP)
       ENDIF
       IF(NTYPE.EQ.5) THEN
        call runend('DAINVA: NTYPE=5 not implemented')
       ENDIF
       RETURN
C
  132  call runend('idamg=22 not implemented')
       RETURN
C
      ENDIF                ! idamg.le.20
C
      RETURN
      END

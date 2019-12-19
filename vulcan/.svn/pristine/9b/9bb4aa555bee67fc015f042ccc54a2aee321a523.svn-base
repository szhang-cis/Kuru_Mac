      SUBROUTINE LDGR05S(DVOLUS,PROPSS,SHAPES,WORK1S,ELDISS,VELCMS,
     .                   EHISTS,ILDGRS,WHAPES,TEINIS,FPCHLS,DVOLIS,
     .                   ADVEMS,CARTDS,XJACMS)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE EQUIVALENT "NODAL FORCES" DUE TO 
C     THE INTERNAL HEAT ( ELEMENT NO. 5 )
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_oms.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'   ! thermal-mechanical
      INCLUDE 'nued_om.f'   ! thermal-microstructural
      INCLUDE 'nuef_om.f'   ! thermal-flow
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      COMMON/SUVOHEAS/SVHEAS
C
      DIMENSION DVOLUS(*),        PROPSS(*),
     .          SHAPES(NNODLS,*), WORK1S(*),
     .          ELDISS(NDOFCS,*), VELCMS(*),
     .          EHISTS(NHISTS,*)
      DIMENSION WHAPES(NNODLS,*), TEINIS(NDOFCS,*),
     .          FPCHLS(NFPCH,*),  DVOLIS(*)
      DIMENSION ADVEMS(NDIMES,*), CARTDS(NDIMES,NNODLS,*),
     .          XJACMS(NDIMES,*)
      DIMENSION TSTRAV(6),        TSTREV(6)
C
C**** ENTER LOOPS FOR AREA NUMERICAL INTEGRATION
C
      DO 100 IGAUSS=1,NGAULS
C
      DTEMPS=0.0
      TGAUSS=0.0
      TGAUIS=0.0
      PSEUDOS=0.0
      DO INODLS=1,NNODLS
       IF(KDYNAS.EQ.1)
     . DTEMPS=DTEMPS+SHAPES(INODLS,IGAUSS)*VELCMS(INODLS)
       TGAUSS=TGAUSS+SHAPES(INODLS,IGAUSS)*ELDISS(1,INODLS)
       TGAUIS=TGAUIS+SHAPES(INODLS,IGAUSS)*TEINIS(1,INODLS)
       IF(IFILLS.EQ.1) THEN
        IF(IMICR.EQ.0) THEN
         IPSEUS=2*NNUPTS+1
        ELSE
         IPSEUS=2*NNUPTS+NNUPO+1
        ENDIF
        PSEUDOS=PSEUDOS+SHAPES(INODLS,IGAUSS)*FPCHLS(IPSEUS,INODLS)
       ENDIF
      END DO
C
C**** UPDATES THE DENSITY
C
      ILAHES=-1
      ISINRS=2
C
      IF(NMEMO3S.EQ.0) THEN
       CALL CAPCOFS( BASMM, BASCC,PROPSS,TGAUSS,
     .              DTEMPS,SOUR1S,SOUR2S,ILAHES,ISINRS,
     .              DSOURS,COUTDS, BASMI,TGAUIS,PSEUDOS)
      ELSE
       IX=NDIMETOS-1
       IF(IMICR.EQ.1) THEN
        BASMM=EHISTS(1,IGAUSS)      ! density computed in micr05t
        BASMI=EHISTS(6+IX,IGAUSS)   ! density computed in micr05t
       ELSE
        IF(ICONVS.EQ.0) THEN
         IDECIS=0
         IF(NMEMO10S.EQ.1) IDECIS=1              ! density changes
         IF(ILDGRS.EQ.0) THEN
          IF(KDYNAS.EQ.0) THEN
           IF(IDECIS.EQ.0) THEN
            CALL CAPCOFS(EHISTS(1,IGAUSS), BASCC,PROPSS,TGAUSS,
     .                   DTEMPS,SOUR1S,SOUR2S,ILAHES,ISINRS,
     .                 DSOURS,COUTDS,EHISTS(6+IX,IGAUSS),TGAUIS,PSEUDOS)
           ENDIF
          ENDIF
         ELSE
          IF(KDYNAS.EQ.0) THEN
           IF(IDECIS.EQ.0) THEN
            CALL CAPCOFS(EHISTS(1,IGAUSS), BASCC,PROPSS,TGAUSS,
     .                   DTEMPS,SOUR1S,SOUR2S,ILAHES,ISINRS,
     .                 DSOURS,COUTDS,EHISTS(6+IX,IGAUSS),TGAUIS,PSEUDOS)
           ENDIF
          ELSE
           IF(IDECIS.EQ.0) THEN ! warning: EHISTT(1) & EHISTT(6+IX) will
C                               ! recomputed in frdy05t.f
            CALL CAPCOFS(EHISTS(1,IGAUSS), BASCC,PROPSS,TGAUSS,
     .                   DTEMPS,SOUR1S,SOUR2S,ILAHES,ISINRS,
     .                 DSOURS,COUTDS,EHISTS(6+IX,IGAUSS),TGAUIS,PSEUDOS)
           ENDIF
          ENDIF
         ENDIF                       ! ildgrs.eq.0
        ENDIF                        ! iconvs.eq.0
        BASMM=EHISTS(1,IGAUSS)       ! computed in frin05s.f if IDECIS=1
C   ! or in fric05s.f if ICONVS=1 or in frdy05s.f if KDYNAS=1 & ILDGRS=0
        BASMI=EHISTS(6+IX,IGAUST)    ! computed in frin05t.f if IDECIS=1
C   ! or in fric05s.f if ICONVS=1 or in frdy05s.f if KDYNAS=1 & ILDGRS=0
       ENDIF
      ENDIF
C
C**** COMPUTES DENSITY
C
      BASMM=BASMI
C
      IF(ITERME.GT.0) THEN                      ! bidirectional coupling

       call runends('iterme > 0 not implemented in ldgr05s')

c      IF(ITERMD.GT.0) THEN                     ! deformed shape
c       IF(LARGET.EQ.1) THEN                    ! TLF
c        BASMM=BASMI
c       ENDIF
c       IF(LARGET.EQ.2) THEN                    ! ULF
c        detjj=1.0      ! determinant of F (from mechanical computation)
c        BASMM=BASMI/DETJJ
c       ENDIF
c       IF(LARGET.EQ.3) THEN                    ! Eulerian
c        BASMM=BASMI*DVOLIT(IGAUST)/DVOLUT(IGAUST)
c       ENDIF
c      ENDIF
      ENDIF
C
      DENSET=BASMM
C
C**** REFERENCE SOURCE
C
      IF(ILDGRT.EQ.0) THEN
       SOURCT=SVHEAT
      ELSE
       SOURCT=0.0                  ! this should be improved !
C
       IF(ITERMEF.GT.0) THEN

        call runends('itermef > 0 not implemented in ldgr05s')

c       ITERMPF=0                  ! better as input
c       IF(ITERMPF.GT.0) THEN
C
c        VISCOT=VISC1F
c        IF(IFILL.EQ.1) VISCOT=PSEUDO*VISC1F+(1.0-PSEUDO)*VISC2F
C
C**** CALCULATE THE RATE-OF-DEFORMATION TENSOR
C
c        CALL PROMA2(XJACMT,ADVEMT,CARTDT(1,1,IGAUST),NDIMET,NNODLT)
C
c        GO TO (1,2,3), NDIMET     ! NTYPET should be used
C
c   1    NSTR1V=1
c        TSTRAV(1)=XJACMT(1,1)
c        TSTREV(1)=2.0*VISCOT*TSTRAV(1)
c        GOTO 6
c   2    NSTR1V=3
c        TSTRAV(1)=XJACMT(1,1)
c        TSTRAV(2)=XJACMT(2,2)
c        TSTRAV(3)=XJACMT(1,2)+XJACMT(2,1)
c        TSTREV(1)=2.0*VISCOT*TSTRAV(1)
c        TSTREV(2)=2.0*VISCOT*TSTRAV(2)
c        TSTREV(3)=    VISCOT*TSTRAV(3)
c        GOTO 6
c   3    NSTR1V=6
c        TSTRAV(1)=XJACMT(1,1)
c        TSTRAV(2)=XJACMT(2,2)
c        TSTRAV(3)=XJACMT(1,2)+XJACMT(2,1)
c        TSTRAV(4)=XJACMT(3,3)
c        TSTRAV(5)=XJACMT(1,3)+XJACMT(3,1)
c        TSTRAV(6)=XJACMT(2,3)+XJACMT(3,2)
c        TSTREV(1)=2.0*VISCOT*TSTRAV(1)
c        TSTREV(2)=2.0*VISCOT*TSTRAV(2)
c        TSTREV(3)=    VISCOT*TSTRAV(3)
c        TSTREV(4)=2.0*VISCOT*TSTRAV(4)
c        TSTREV(5)=    VISCOT*TSTRAV(5)
c        TSTREV(6)=    VISCOT*TSTRAV(6)
c        GOTO 6
c   6    CONTINUE
C
C**** COMPUTES THE VISCOUS ENERGY TERM
C
c        VISCET=0.0
c        DO ISTR1V=1,NSTR1V
c         VISCET=VISCET+TSTREV(ISTR1V)*TSTRAV(ISTR1V)
c        ENDDO
c        SOURCT=VISCET
c       ENDIF                      ! itermpf.gt.0
       ENDIF                       ! itermef.gt.0
c      SOURCT=-SOURCT              ! treated as a residual heat
      ENDIF                        ! ildgrt.eq.1
C
C**** EVALUATES THE SOURCE TERM
C
      CALL SOURCTS(DENSES,SOURCS,QINTES)
C
C**** CALCULATE LOADS AND ASSOCIATE WITH ELEMENT NODAL POINTS
C
      IEVABS=0
      DO 130 INODES=1,NNODLS
      SHAVOS=WHAPES(INODES,IGAUSS)*DVOLUS(IGAUSS)
C
      IEVABS=IEVABS+1
      WORK1S(IEVABS)=WORK1S(IEVABS)+QINTES*SHAVOS
  130 CONTINUE
C
  100 CONTINUE
C
      RETURN
      END

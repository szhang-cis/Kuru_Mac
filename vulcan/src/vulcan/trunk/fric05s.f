      SUBROUTINE FRIC05S(CARTDS,DVOLUS,EHISTS,ELCODS,ELDISS,EPMTXS,
     .                   GPCODS,LNODSS,PROPSS,RMAT1S,SHAPES,STRANS,
     .                   STRSGS,STRA0S,STRS0S,
     .                   BMSIGS,BMATXS,DESIGS,DMATXS,DSTRAS,PRESGS,
     .                   SGTOTS,SIGMAS,TSTRAS,XJACMS,ELELTS,VELCMS,
     .                   WHAPES,ADVEMS,TEINIS,FPCHLS)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE INTERNAL RESISTING HEATS DUE TO
C     CONVECTION EFFECTS
C     ( FOR ELEMENT NO. 5 )
C
C***********************************************************************
C
C     Index of variables:
C
C     EHISTT(   1) = Density
C     EHISTT(   2) = Specific Heat coefficient
C     EHISTT(   3) = Isotropic Conductivity or Conduct. x (orthot. mat.)
C     EHISTT(   4) = Conductivity y (orthotropic material)
C     EHISTT(   5) = Conductivity z (orthotropic material)
C     EHISTT(3:11) = Conductivity for the fully anisotropc mat. (3D)
C     EHISTT(4+IX) = L*Phase-change function
C     EHISTT(5+IX) = L*Phase-change function rate
C     EHISTT(6+IX) = Initial density
C     EHISTT(7+IX) = Coupling coefficient
C     EHISTT(8+IX) = Temperature derivative of phase-change function
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
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DIMENSION CARTDS(NDIMES,NNODLS,*), DVOLUS(*),
     .          EHISTS(NHISTS,*),        ELCODS(NDIMES,*),
     .          ELDISS(NDOFCS,*),
     .          EPMTXS(*),               GPCODS(NDIMES,*), 
     .          LNODSS(*),               PROPSS(*),
     .          RMAT1S(NDIMES,*),        SHAPES(NNODLS,*),
     .          STRANS(NSTR1S,*),        STRSGS(NSTR1S,*),
     .          STRA0S(NSTR1S,*),        STRS0S(NSTR1S,*)
      DIMENSION BMATXS(NSTR1S,*),        BMSIGS(*),
     .          DESIGS(*),               DMATXS(NSTR1S,*),
     .          DSTRAS(*),               PRESGS(*),
     .          SGTOTS(*),               SIGMAS(*),
     .          TSTRAS(*),               XJACMS(NDIMES,*),
     .          ELELTS(*),               VELCMS(*)
      DIMENSION WHAPES(NNODLS,*),        ADVEMS(NDIMES,*),
     .          TEINIS(NDOFCS,*),        FPCHLS(NFPCH,*)
      DIMENSION VELOCS(3)
C
C**** LOOP OVER INTEGRATION POINTS
C
      DO IGAUSS=1,NGAULS
C
C**** COMPUTE TEMPERATURE AT GAUSSIAN POINT
C
       DTEMPS=0.0
       TGAUSS=0.0
       TGAUIS=0.0
       PSEUDOS=0.0
       DO INODLS=1,NNODLS
        IF(KDYNAS.EQ.1)
     .  DTEMPS=DTEMPS+SHAPES(INODLS,IGAUSS)*VELCMS(INODLS)
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
       DO IDIMES=1,NDIMES
        VELOCS(IDIMES)=0.0
        DO INODLS=1,NNODLS
         VELOCS(IDIMES)=VELOCS(IDIMES)+SHAPES(INODLS,IGAUSS)*
     .                  ADVEMS(IDIMES,INODLS)
        ENDDO
       END DO
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + alpha dt)
C
       IF(KDYNAS.EQ.1) THEN
        IF(KINTES.EQ.1)       ! Euler method
     .   TGAUSS=TALFAS*TGAUSS+(1.0-TALFAS)*(TGAUSS-DTEMPS*DTIMES)
       ENDIF
C
C**** UPDATE THE DENSITY & CAPACITY COEFFICIENT
C
       ILAHES=0
       ISINRS=2
C
       IF(NMEMO3S.EQ.0) THEN
        CALL CAPCOFS(COEFMS,COEFCS,PROPSS,TGAUSS,
     .               DTEMPS,SOUR1S,SOUR2S,ILAHES,ISINRS,
     .               DSOURS,COUTDS, BASMI,TGAUIS,PSEUDOS)
       ELSE
        IX=NDIMETOS-1
        IF(IMICR.EQ.1) THEN
         COEFMS=EHISTS(1,IGAUSS)  ! COEFMS=Density computed in micr05s.f
         BASMI= EHISTS(6+IX,IGAUSS) ! Init. density comp. in micr05s.f
         COEFCS=EHISTS(2,IGAUSS)  ! COEFCS=Capacity comp. in micr05s.f
        ELSE
         IDECIS=0
         IF(NMEMO10S.EQ.1) IDECIS=1              ! density changes
         IF(IDECIS.EQ.0) THEN
          CALL CAPCOFS(EHISTS(1,IGAUSS),EHISTS(2,IGAUSS),PROPSS,TGAUSS,
     .                 DTEMPS,SOUR1S,SOUR2S,ILAHES,ISINRS,
     .                 DSOURS,EHISTS(7+IX,IGAUSS),EHISTS(6+IX,IGAUSS),
     .                                                   TGAUIS,PSEUDOS)
         ELSE
          CALL CAPCOFS(COEFMS,EHISTS(2,IGAUSS),PROPSS,TGAUSS,
     .                 DTEMPS,SOUR1S,SOUR2S,ILAHES,ISINRS,
     .                 DSOURS,EHISTS(7+IX,IGAUSS),BASMI,TGAUIS,PSEUDOS)
         ENDIF
         COEFMS=EHISTS(1,IGAUSS)     ! computed in frin05s.f if IDECIS=1
         BASMI= EHISTS(6+IX,IGAUSS)  ! computed in frin05s.f if IDECIS=1
         COEFCS=EHISTS(2,IGAUSS)     ! COEFCT=Capacity
        ENDIF
       ENDIF
C
C**** COMPUTES DENSITY
C
       COEFMT=BASMI
C
C**** COMPUTE DENSITY*CAPACITY*VELOCITY
C
       DO IDIMES=1,NDIMES               ! nstr1s=ndimes
        DMATXS(IDIMES,1)=COEFMS*COEFCS*VELOCS(IDIMES)
       ENDDO
C
C**** CALCULATE THE CONVECTIVE TERM
C
       CALL LIHEATCS(DSTRAS,TSTRAS,DMATXS,CONVES,NDIMES,NSTR1S,
     .               NNODLS,NDOFCS,CARTDS(1,1,IGAUSS),
     .               VELCMS,ELDISS,KDYNAS)
C
C**** CALCULATES THE "CONVECTIVE FORCES"
C
       CALL EQHEATC(BMSIGS,WHAPES(1,IGAUSS),DVOLUS(IGAUSS),CONVES,
     .              NNODLS,NEVABS)
C
      END DO ! IGAUSS=1,NGAULS
C
      RETURN
      END

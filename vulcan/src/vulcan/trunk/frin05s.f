      SUBROUTINE FRIN05S(CARTDS,DVOLUS,EHISTS,ELCODS,ELDISS,EPMTXS,
     .                   GPCODS,LNODSS,PROPSS,RMAT1S,SHAPES,STRANS,
     .                   STRSGS,STRA0S,STRS0S,
     .                   BMSIGS,BMATXS,DESIGS,DMATXS,DSTRAS,PRESGS,
     .                   SGTOTS,SIGMAS,TSTRAS,XJACMS,ELELTS,VELCMS,
     .                   WARTDS,TEINIS,FPCHLS,DVOLIS)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE RESISTING HEATS 
C               ( FOR ELEMENT NO. 5 )
C
C***********************************************************************
C
C     Index of variables:
C
C     IX=NDIMETO-1
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
      DIMENSION WARTDS(NDIMES,NNODLS,*), TEINIS(NDOFCS,*),
     .          FPCHLS(NFPCH,*),         DVOLIS(*)
      DIMENSION COEFKS(9)
C
C**** INITIALIZE ELEMENT CONTRIBUTION OF ELELTT (K*T)
C
      DO IEVABS=1,NNODLS
       ELELTS(IEVABS)=0.0D00
      END DO
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
         PSEUDOS=PSEUDOS+SHAPEs(INODLS,IGAUSS)*FPCHLS(IPSEUS,INODLS)
        ENDIF
       END DO
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + alpha dt)
C
       IF(KDYNAS.EQ.1) THEN
        IF(KINTES.EQ.1)              ! Euler method
     .   TGAUSS=TALFAS*TGAUSS+(1.0-TALFAS)*(TGAUSS-DTEMPS*DTIMES)
       ENDIF
C
C**** UPDATES CONDUCTIVITY COEFFICIENT (COEFKT=Conductivitity)
C
       ILAHES=-1
       ISINRS=2
C
       IF(NMEMO3S.EQ.0) THEN
        IDECIS=0
        IF(NMEMO10S.EQ.1) IDECIS=1                ! density changes
        IF(IDECIS.EQ.1)
     .   CALL CAPCOFS(BASMM,BASCC,PROPSS,TGAUSS,
     .                DTEMPS,SOUR1S,SOUR2S,ILAHES,ISINRS,
     .                DSOURS,COUTDS, BASMI,TGAUIS,PSEUDOS)
        CALL CONCOFS(COEFKS,PROPSS,TGAUSS,PSEUDOS)
       ELSE
        IX=NDIMETOS-1
        IF(IMICR.EQ.1) THEN

         call runends('imicr=1 in frin05s.f')

         BASMM= EHISTS(1,IGAUSS)                 ! computed in micr05t.f
         BASMI= EHISTS(6+IX,IGAUSS)              ! computed in micr05t.f
         DO I=1,NDIMETOS
          COEFKS(I)=EHISTS(3+I-1,IGAUSS)         ! computed in micr05t.f
         ENDDO
        ELSE
         IDECIS=0
         IF(NMEMO10S.EQ.1) IDECIS=1              ! density changes
         IF(IDECIS.EQ.1)
     .    CALL CAPCOFS(EHISTS(1,IGAUSS),BASCC,PROPSS,TGAUSS,
     .                 DTEMPS,SOUR1S,SOUR2S,ILAHES,ISINRS,
     .                 DSOURS,COUTDS,EHISTS(6+IX,IGAUSS),TGAUIS,PSEUDOS)
         CALL CONCOFS(EHISTS(3,IGAUSS),PROPSS,TGAUSS,PSEUDOS)
         IF(IDECIS.EQ.1) THEN
          BASMM= EHISTS(1,IGAUSS)
          BASMI= EHISTS(6+IX,IGAUSS)
         ENDIF
         DO I=1,NDIMETOS
          COEFKS(I)=EHISTS(3+I-1,IGAUSS)
         ENDDO
        ENDIF
       ENDIF
C
C**** COMPUTE CONDUCTIVITY TENSOR     (only isotropic case available !!)
C
       IF(NMEMO10S.EQ.1) THEN
        IF(BASMM.GT.0.0) THEN
         DO I=1,NDIMETOS
          COEFKS(I)=COEFKS(I)*(BASMI/BASMM)
         ENDDO
        ENDIF
       ENDIF
C
       IF(ITERME.GT.0) THEN                     ! bidirectional coupling

        call runends('iterme > 0 in frin05s.f')

        IF(ITERMD.GT.0) THEN                    ! deformed shape
         IF(LARGES.EQ.1) THEN                   ! TLF
          DO I=1,NDIMETOS
           detjj=1.0    ! determinant of F (from mechanical computation)
           BASMM=BASMI/DETJJ
           COEFKS(I)=COEFKS(I)*(BASMI/BASMM)
          ENDDO
         ENDIF
         IF(LARGES.EQ.2) THEN                   ! ULF
          DO I=1,NDIMETOS
           COEFKS(I)=COEFKS(I)
          ENDDO
         ENDIF
         IF(LARGES.EQ.3) THEN                   ! Eulerian
          DO I=1,NDIMETOS
           COEFKS(I)=COEFKS(I)
           IF(IFREKS.EQ.2)
     .      COEFKS(I)=COEFKS(I)*DVOLIS(IGAUSS)/DVOLUS(IGAUSS)
          ENDDO
         ENDIF
        ENDIF
       ENDIF
C
       I=1
       DO IDIMES=1,NDIMES                ! nstr1t=ndimet
        DO JDIMES=1,NDIMES
         DMATXS(IDIMES,JDIMES)=0.0
         IF(ISOTRS.EQ.0.OR.ISOTRS.EQ.1) THEN
          IF(IDIMES.EQ.JDIMES) THEN
           IF(ISOTRS.EQ.1) I=IDIMES
           DMATXS(IDIMES,JDIMES)=COEFKS(I)
          ENDIF
          IF(ISOTRS.EQ.2) THEN
           CALL RUNENDS('ERROR: FULLY ANIST. MAT. NOT IMPLEMENTED')
          ENDIF
         ENDIF
        ENDDO
       ENDDO
C
C**** CALCULATE THE TOTAL HEAT FLUX AT GAUSSIAN POINTS
C
       CALL LIHEATS(CARTDS(1,1,IGAUSS),DEPSVS,SGTOTS,
     .              DSTRAS,
     .              ELDISS,GPCODS(1,IGAUSS),LARGES,NDIMES,NDOFNS,NNODLS,
     .              NSTR1S,
     .              PROPSS,SHAPES(1,IGAUSS),STRANS(1,IGAUSS),
     .              STRA0S(1,IGAUSS),
     .              TSTRAS,XJACMS,NDOFCS,SIGMAS,VELCMS,DMATXS,KDYNAS)
C
C**** STORE TOTAL HEAT FLUX IN THE PROPER ARRAY
C
       IF(NMEMO4S.EQ.1) THEN
        DO ISTR1S=1,NSTR1S
         STRSGS(ISTR1S,IGAUSS)=-SGTOTS(ISTR1S)
         STRANS(ISTR1S,IGAUSS)= TSTRAS(ISTR1S)
        ENDDO
       ENDIF
C
C**** INTEGRATE THE HEATS INTO THE INTERNAL "HEAT FORCES"
C
       CALL EQHEATT(BMATXS,BMSIGS,CARTDS(1,1,IGAUSS),DVOLUS(IGAUSS),
     .              ELDISS,
     .              GPCODS(1,IGAUSS),LARGES,NDIMES,NDOFCS,NDOFNS,
     .              NEVABS,NNODLS,
     .              NSTR1S,NTYPES,SHAPES(1,IGAUSS),SGTOTS,XJACMS,
     .              SIGMAS,ELELTS,WARTDS(1,1,IGAUSS))
C 
      END DO      ! igauss=1,ngauls
C
      RETURN
      END

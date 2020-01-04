      SUBROUTINE SEI101T(DVOLUT,ELCODT,GPCODT,LNODST,PROPST,
     .                   SHAPET,THICKT,ELDIST,
     .                   DVOL1T,GPCO1T,SHAP1T,
     .                   DERIVT,POSGPT,WEIGPT,XJACMT,
     .                   ELCO1T,ELDI1T,DVOLIT,ELCOIT,DISPTT,ITASKT)
C***********************************************************************
C
C**** THIS ROUTINE SETS UP SOME NEEDED CONSTANT MATRICES FOR
C     FUTURE USE ( ELEMENT NO. 5 )
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      COMMON/JACOBSTA/IERORT,KERORT
C
      DIMENSION DVOLUT(*),
     .          ELCODT(NDIMET,*),        GPCODT(NDIMET,*),
     .          LNODST(*),               PROPST(*),
     .          SHAPET(NNODLT,*)
      DIMENSION DERIVT(NDIMLT,NNODLT,*), POSGPT(NDIMLT,*),
     .          WEIGPT(*),               XJACMT(NDIMET,*),
     .          ELCO1T(NDIMET,*)
      DIMENSION DVOL1T(*),
     .          GPCO1T(NDIMET,*),        SHAP1T(NNODLT,*)
      DIMENSION ELDIST(NDOFCT,*),        ELDI1T(NDOFCT,*)
      DIMENSION DVOLIT(*),               ELCOIT(NDIMET,*),
c    .          DISPTT(NDOFCM,*),        DISRTT(NDOFCM,*)
     .          DISPTT(NDOFCM,*)
C
C*** CHECK THE RIGHT USE OF THIS ROUTINE
C
      IF(NDIMET.NE.1) THEN
       IF(NDIMET.LE.NDIMLT)
     .  CALL RUNENDT('ERROR: NDIMET LE NDIMLT IN SEI101T')
      ENDIF
C
C**** OBTAINS MATERIAL PROPERTIES
C
      CALL IDENPHT(PROPST,    2)
C
      IF(NMEMO1.EQ.0) THEN
       DO IDIMET=1,NDIMET
        DO INODLT=1,NNODLT
         ELCO1T(IDIMET,INODLT)=ELCODT(IDIMET,INODLT)
        ENDDO
       ENDDO
      ENDIF
C
      IF(NMEMO5.EQ.0) THEN
       DO IDOFCT=1,NDOFCT
        DO INODLT=1,NNODLT
         ELDI1T(IDOFCT,INODLT)=ELDIST(IDOFCT,INODLT)
        ENDDO
       ENDDO
      ENDIF
C
C**** DEALS WITH DEFORMED SHAPE
C
C     Notes:
C
C     LARGET=1 => ELCO1T: material coord. & ELCOIT=spatial coord.
C     LARGET=2 => ELCO1T: updated coord. & ELCOIT=spatial coord.
C     LARGET=3 => ELCO1T: spatial coord. & ELCOIT=material coord.
C
C     DISPTT: total displacements
C     DISRTT: incremental displacements
C
      IF(LARGET.NE.0) THEN      ! large strains-displacements

       IF(LARGET.EQ.1) THEN
        DO IDIMET=1,NDIMET
         DO INODLT=1,NNODLT
          ELCOIT(IDIMET,INODLT)=ELCO1T(IDIMET,INODLT)+
     .                          DISPTT(IDIMET,INODLT)
         ENDDO
        ENDDO
       ENDIF                    ! larget.eq.1
       IF(LARGET.EQ.2) THEN
        DO IDIMET=1,NDIMET
         DO INODLT=1,NNODLT
          ELCOIT(IDIMET,INODLT)=ELCO1T(IDIMET,INODLT)+
     .                          DISPTT(IDIMET,INODLT)
         ENDDO
        ENDDO
        DO IDIMET=1,NDIMET
         DO INODLT=1,NNODLT
          CALL RUNENDT('ERROR: LARGE=2 NOT IMPLEMENTED - SETI05T')
c         ELCO1T(IDIMET,INODLT)=ELCOIT(IDIMET,INODLT)-
c    .                          DISRTT(IDIMET,INODLT)
         ENDDO
        ENDDO
       ENDIF                    ! larget.eq.2
       IF(LARGET.EQ.3) THEN
        DO IDIMET=1,NDIMET
         DO INODLT=1,NNODLT
          ELCOIT(IDIMET,INODLT)=ELCO1T(IDIMET,INODLT)
         ENDDO
        ENDDO
        DO IDIMET=1,NDIMET
         DO INODLT=1,NNODLT
          ELCO1T(IDIMET,INODLT)=ELCOIT(IDIMET,INODLT)+
     .                          DISPTT(IDIMET,INODLT)
         ENDDO
        ENDDO
       ENDIF                    ! larget.eq.3
      ELSE                      ! small strains-displacements
       DO IDIMET=1,NDIMET
        DO INODLT=1,NNODLT
         ELCOIT(IDIMET,INODLT)=ELCO1T(IDIMET,INODLT)
        ENDDO
       ENDDO
      ENDIF                     ! larget.ne.0
C
      IF(NMEMO2.EQ.0) THEN
       IF(NRULET.GT.5)
     .  CALL RUNENDT('ERROR: NRULET GT 5 WITH NMEMO2=0')
       IF(ITASKT.EQ.1) THEN
        CALL SET101T(DVOLUT,ELCO1T,GPCODT,LNODST,PROPST,
     .               SHAPET,THICKT,
     .               DERIVT,POSGPT,WEIGPT,XJACMT)
       ELSE
        IF(LARGET.NE.0)                         ! just to compute DVOLIT
     .   CALL SET101T(DVOLIT,ELCOIT,GPCO1T,LNODST,PROPST,
     .                SHAP1T,THICKT,
     .                DERIVT,POSGPT,WEIGPT,XJACMT)
        DO IGAUST=1,NGAULT
         DVOL1T(IGAUST)=DVOLUT(IGAUST)
         DO INODLT=1,NNODLT
          SHAP1T(INODLT,IGAUST)=SHAPET(INODLT,IGAUST)
          DO IDIMET=1,NDIMET
	   GPCO1T(IDIMET,IGAUST)=GPCODT(IDIMET,IGAUST)
	  ENDDO
	 ENDDO
        ENDDO
       ENDIF                ! itaskt.eq.1
      ELSE
C
C**** NMEMO2=1 >> THERMAL PROBLEM INDEPENDENT OF ITASK
C
C     ITASK=1: checks negative volumes
C     ITASK>1: compute all the necessary arrays
C
       IF(LARGET.NE.0)                          ! just to compute DVOLIT
     . CALL SET101T(DVOLIT,ELCOIT,GPCO1T,LNODST,PROPST,   
     .              SHAP1T,THICKT,   
     .              DERIVT,POSGPT,WEIGPT,XJACMT)   
       CALL SET101T(DVOL1T,ELCO1T,GPCO1T,LNODST,PROPST,
     .              SHAP1T,THICKT,
     .              DERIVT,POSGPT,WEIGPT,XJACMT)
      ENDIF                 ! nmemo2.eq.1
C
      RETURN
      END
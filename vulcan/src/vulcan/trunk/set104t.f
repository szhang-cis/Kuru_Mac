      SUBROUTINE SET104T(DVOLUT,ELCODT,GPCODT,LNODST,PROPST,SHAPET,
     .                   THICKT,DERIVT,POSGPT,WEIGPT,XJACMT,VNORLT,
     .                   STRSGT,ITASKT)
C***********************************************************************
C
C**** THIS ROUTINE SETS UP SOME NEEDED CONSTANT MATRICES FOR
C     FUTURE USE ( ELEMENT NO. 104 )
C
C     Note: SHAPET & DERIVT must be dimensioned with
C           NNODNT due to non-coincident contact mesh (NOCOLT=1)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'     ! thermal-mechanical
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      COMMON/JACOBSTA/IERORT,KERORT
C
      DIMENSION DERIVT(NDIMLT,NNODNT,*), DVOLUT(*),
     .          ELCODT(NDIMET,NNODLT),
     .          GPCODT(NDIMET,*),        LNODST(*), PROPST(*),
     .          SHAPET(NNODNT,*)
      DIMENSION POSGPT(NDIMLT,*),        WEIGPT(*),
     .          XJACMT(NDIMET,NDIMLT)
      DIMENSION VNORLT(NDIMET,*)
      DIMENSION STRSGT(NSTR1T,*)        ! useful for non-coincident mesh
C
      DIMENSION EXIDI(3), ETADI(3), PVECT(3,9)
C
      TWOPIT=6.283185307179586D0
C
C**** INITIALIZATION
C
      NNOBOT=NNODLT/2
      IF(NOCOLT.EQ.1) NNOBOT=NNODLT                ! non-coincident mesh
C
C**** EVALUATE QUANTITIES ASSOCIATED WITH THE GAUSSIAN POINTS
C
      if(ndimet.eq.1) then
       posgpt(1,1)=1.0D0
       weigpt(1)=1.0D0
       shapet(1,1)=1.0D0
       derivt(1,1,1)=1.0D0
       dvolut(1)=thickt
       return
      endif
C
C**** IDENTIFY SAMPLING POINTS AND WEIGHTS
C
      CALL RULEPW(NDIMLT,NGAULT,NRULET,POSGPT,WEIGPT)
C
C**** LOOP ON INTEGRATION POINTS
C
      DO IGAUST=1,NGAULT
C
C**** COMPUTE SHAPE FUNCTIONS AND DERIVATIVES
C
       EXISPT=POSGPT(1,IGAUST)
       ETASPT=0.0D0
       IF(NDIMLT.EQ.2) ETASPT=POSGPT(2,IGAUST)
C
       CALL SHAFUN(DERIVT(1,1,IGAUST),EXISPT,ETASPT,EZETAT,NDIMLT,
     .             NNOBOT,     1,     0,SHAPET(1,IGAUST))
C
       IF(NOCOLT.EQ.1) THEN                        ! non-coincident mesh
        IAUXZT=1
        EXISPT=STRSGT(IAUXZT,IGAUST)           ! unnecessary for itask=1
        ETASPT=0.0D0
        IF(NDIMLT.EQ.2) ETASPT=STRSGT(IAUXZT+1,IGAUST)
C
        CALL SHAFUN(DERIVT(1,1,IGAUST),EXISPT,ETASPT,EZETAT,NDIMLT,
     .              NNOBOT,     1,NNOBOT,SHAPET(1,IGAUST))
       ENDIF
C
C**** COMPUTE JACOBIAN ( NOTE CARTESIAN DERIVATIVES ARE NOT COMPUTED )
C
       CALL JACBOUT(DERIVT(1,1,IGAUST),DETJMT,ELCODT,GPCODT(1,IGAUST),
     .              IELEMT,NDIMET,NDIMLT,NNOBOT,SHAPET(1,IGAUST),
     .              XJACMT,LUREST,LUPRIT)
       IF(IERORT.NE.0) RETURN
C
C**** INTEGRATION WEIGHTS
C
                         DVOLUT(IGAUST)=WEIGPT(IGAUST)*DETJMT
       IF(NTYPET.EQ.3)   DVOLUT(IGAUST)=DVOLUT(IGAUST)*TWOPIT*
     .                                  GPCODT(1,IGAUST)
       IF(NTYPET.EQ.5)   DVOLUT(IGAUST)=DVOLUT(IGAUST)*THICKT ! 1-D
C
      END DO
C
C**** CHECKS GAUSSIAN VOLUME FOR AXISYMMETRIC PROBLEMS FOR OUTPUT OPER.
C
      IF(NTYPET.EQ.3) THEN
       CALL GAUCEK(DVOLUT,   104,     1,NGAULT,NNODLT,ICEKET)
       IF(ICEKET.NE.0) CALL RUNENDT('ERROR:WRONG NGAUS IN gaucek.f-104')
      ENDIF
C
C**** COMPUTES OUTWARD NORMAL
C
      INDEXX=0
      IF(ITERME.GT.0) THEN                       ! bidirectional coupled
       IF(ITERMG.GT.0) THEN                      ! gap dependency
        IF(NMEMO11.EQ.1) THEN
         IF(INDEXX.EQ.0) INDEXX=1
        ENDIF
       ENDIF
      ENDIF
      IF(NOCOLT.EQ.1) THEN                       ! non-coincident mesh
       IF(LARGET.EQ.0) THEN
        IF(ITASKT.EQ.1) THEN
         IF(INDEXX.EQ.0) INDEXX=1
        ENDIF
       ELSE
        IF(INDEXX.EQ.0) INDEXX=1
       ENDIF
      ENDIF
C
      IF(INDEXX.EQ.1) THEN
C
C**** LOOP OVER INTEGRATION POINTS
C
       DO 100 IGAUST=1,NGAULT
C
C**** COMPUTES OUTWARD UNIT NORMAL TO BODY 1 (should be INOTE ?)
C
C     Note: the outward unit normal vector of a body (1 in this case) is
C           obtained if the correct (i.e., anti-clockwise) node
C           connectivity is defined for the contact element looking at
C           it from the body (this is valid for 2D & 3D problems).
C
       IF(NDIMET.EQ.1) THEN                    ! 1D CASE
        EXIDI(1)=EXIDI(1)+ELCODT(IDIMET,1)*DERIVT(1,1,IGAUST)
        EXINM=EXIDI(1)*EXIDI(1)
        PNORM=EXINM
        PVECT(1,IGAUST)=-EXIDI(1)                     ! normal to body 2
        PVECT(1,IGAUST)= EXIDI(1)                     ! normal to body 1
C
       ELSE IF(NDIMET.EQ.2) THEN               ! 2D CASE
        EXINM=0.0D0
        DO IDIMET=1,NDIMET
         EXIDI(IDIMET)=0.0D0
         DO INODET=1,NNOBOT
          EXIDI(IDIMET)=EXIDI(IDIMET)+
     .                  ELCODT(IDIMET,INODET)*DERIVT(1,INODET,IGAUST)
         ENDDO
         EXINM=EXINM+EXIDI(IDIMET)*EXIDI(IDIMET)
        ENDDO
        PNORM=EXINM
        PVECT(1,IGAUST)= EXIDI(2)                     ! normal to body 2
        PVECT(2,IGAUST)=-EXIDI(1)
        PVECT(1,IGAUST)=-EXIDI(2)                     ! normal to body 1
        PVECT(2,IGAUST)= EXIDI(1)
C
       ELSE                                    ! 3D CASE
        DO IDIMET=1,NDIMET          ! compute two vectors on the surface
         EXIDI(IDIMET)=0.0D0
         ETADI(IDIMET)=0.0D0
         DO INODET=1,NNOBOT
          EXIDI(IDIMET)=EXIDI(IDIMET)+DERIVT(1,INODET,IGAUST)*
     .                  ELCODT(IDIMET,INODET)
          ETADI(IDIMET)=ETADI(IDIMET)+DERIVT(2,INODET,IGAUST)*
     .                  ELCODT(IDIMET,INODET)
         ENDDO
        ENDDO
        PNORM=0.0D0
        DO IDIMA=1,NDIMET           ! compute normal as cross product
         IDIMB=IDIMA+1-IDIMA/3*3
         IDIMC=IDIMB+1-IDIMB/3*3
         PVECT(IDIMA,IGAUST)=EXIDI(IDIMB)*ETADI(IDIMC)-
     .                      ETADI(IDIMB)*EXIDI(IDIMC) ! normal to body 2
         PVECT(IDIMA,IGAUST)=-PVECT(IDIMA,IGAUST)     ! normal to body 1
         PNORM=PNORM+PVECT(IDIMA,IGAUST)*PVECT(IDIMA,IGAUST)
        ENDDO
       ENDIF
C
       PNORM=DSQRT(PNORM)           ! compute normalized normal
       DO IDIMET=1,NDIMET
        PVECT(IDIMET,IGAUST)=PVECT(IDIMET,IGAUST)/PNORM
       ENDDO
C
  100  CONTINUE
C
C**** COMPUTES VNORL
C
       DO IDIMET=1,NDIMET
        DO IGAUST=1,NGAULT
         VNORLT(IDIMET,IGAUST)=PVECT(IDIMET,IGAUST)
        ENDDO
       ENDDO
      ENDIF                         ! indexx.eq.1
C
      RETURN
      END

      SUBROUTINE FRIN32(PROPS,LNODS,ELDIS,BMSIG,DVOLU,SHAPE,EHIST,TENOD,
     .                  ELCOD,PWOEL,PREAL,TGAPL,VNORL,VTANL,DISIL,DISPL)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE CONTACT FORCES FOR ELEMENT 32
C
C     Notes:
C     For contact and linking problems, NDIME=NDOFN
C
C     SHAPE, ELDIS & TENOD must be dimensioned with NNODN due to
C     non-coincident contact mesh (NOCOL=1)
C
C     For frictional/debounding problems:
C     EHIST(1:NDIME-1,IGAUS)=frictional tangential force
C     EHIST(NDIME:NDIME+NDIME*NDIME-1,IGAUS)=frictional const. tensor
C
C     For mode I debounding model:
C     EHIST(NDIME+NDIME*NDIME,IGAUS)=debounding parameter
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
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      COMMON/TUNING4/RITEN,RITEF,TRUPL,TRUPM,TOLGA,TOLGAM
C
      DIMENSION PROPS(*),       ELDIS(NDOFC,*),
     .          BMSIG(*),       DVOLU(*),
     .          SHAPE(NNODN,*), EHIST(NHIST,*)
      DIMENSION TENOD(*),       ELCOD(NDIME,*),
     .          PWOEL(*),       PREAL(*),
     .          TGAPL(*),       LNODS(*)
      DIMENSION VNORL(NDIME,*), DGAPS(3),
     .          VTANL(NDIME,NDIME-1,*),
     .          DISIL(NDOFC,*), DISPL(NDOFC,*)
      DIMENSION DGAPT(2),       PREST(2),
     .          RIGI1(3,3),     RIGI2(2,3), RIGI3(2,2)
      REAL(KIND=8), ALLOCATABLE, SAVE :: HISTO(:,:,:),HISTOTEMP(:,:,:)
      LOGICAL                  , SAVE :: INIT=.FALSE.

      
      IF(IELEM.EQ.329) WRITE(LURES,*) "INCIALES,histo,histotemp:"
     .                ,HISTO(329,1,1),HISTOTEMP(329,1,1)
      IF(.NOT. INIT) THEN
         ALLOCATE(HISTO(NELEM,NGAUL,NDIME))
         ALLOCATE(HISTOTEMP(NELEM,NGAUL,NDIME-1))
         DO IDIME=1,NDIME-1
           HISTO(:,:,IDIME)=0
           HISTOTEMP(:,:,IDIME)=0
         END DO
         INIT=.TRUE.
      ENDIF
      
      IF (IITER.EQ.-00) THEN
        DO IDIME=1,NDIME-1
           HISTO(IELEM,:,IDIME)=HISTOTEMP(IELEM,:,IDIME)
        END DO
      END IF
      
      DO IDIME=1,NDIME-1
        HISTOTEMP(IELEM,:,IDIME)=0
      END DO
C
C**** DETERMINES BETWEEN CONTACT OR LINKING
C
      ICOLI=INT(PROPS(1))
      IF(ICOLI.EQ.1) THEN
C
       NNOBO=NNODL/2
       IF(NOCOL.EQ.1) NNOBO=NNODL                  ! non-coincident mesh
       NEVBO=NNOBO*NDOFN
C
C**** PROPERTIES FOR THE CONTACT PROBLEM
C
C     PROPS(2)=RIGIN:NORNAL STIFFNESS
C     PROPS(5)=TEMPL:LIQUIDUS TEMPERATURE
C     PROPS(6)=INOTE:DENOTES WHICH FACE BELONGS TO PIECE
C
       RIGIN=PROPS(2)
       TEMPL=PROPS(5)
       INOTE=INT(PROPS(6))
       ICOMO=INT(PROPS(7))
       IITEN=INT(RITEN)     ! 0: tun. not considered; 1: tun. considered
C
       IF(ICONC.EQ.0) THEN
        IF(IITEN.EQ.0) THEN
         INP14=INT(PROPS(14))
         IF(INP14.EQ.1) THEN
          TOLGAM=PROPS(16)
         ENDIF
        ENDIF
       ENDIF
C
       IFRIC=INT(PROPS(51))
       IF(IFRIC.EQ.1) THEN
        IFRIM=INT(PROPS(52))
        IF(IFRIM.EQ.1) THEN
         IFRCO=INT(PROPS(53))
         IF(IFRCO.EQ.1) THEN
          RIGIT=PROPS(54)
          FRIMU=PROPS(55)
         ENDIF
        ENDIF
        IF(IFRIM.EQ.2) THEN
         RIGIT=PROPS(54)
        ENDIF
       ENDIF
C
C**** INITIALISES MECHANICAL COUPLING TERM (only for friction problems)
C
       IF(ITERME.GT.0) THEN                     ! bidirectional coupling
        DO INODL=1,NNODL
         PWOEL(INODL)=0.0D0
        END DO
        IGAFR=0
        IF(IFRIC.EQ.1) THEN
         IF(IFRIM.EQ.1) THEN
          IF(ITERMF.GT.0) THEN                  ! frictional heating
           IGAFR=INT(PROPS(101))
           AGAFP=0.0D0
           IF(IGAFR.EQ.1) AGAFP=PROPS(102)      ! partition coefficient
          ENDIF
         ENDIF
        ENDIF
       ENDIF
C
C**** LOOP OVER INTEGRATION POINTS
C
       DO 100 IGAUS=1,NGAUL
C
       IF(NOCOL.EQ.1) THEN                              ! checks contact
        TGAUS=0.0D0                   ! initialization for nodal gn & pn
        DGAPN=0.0D0
        PRESN=0.0D0
        TGAPL(IGAUS)=0.0D0
        PREAL(IGAUS)=0.0D0
        IAUXY=(IGAUS-1)*NNODL
        IF(LNODS(NNOBO).EQ.LNODS(NNODL+NNOBO+IAUXY)) GO TO 100
       ENDIF
C
C**** COMPUTES NORMAL GAP
C
C     DGAPN  > 0: CONTACT PROBLEM
C     DGAPN <= 0: NO CONTACT PROBLEM
C
       DGAPN=0.0D0
       DO IDOFN=1,NDOFN
        DO INODL=1,NNOBO
         ELDI1=ELDIS(IDOFN,INODL)
         IF(NOCOL.EQ.0) THEN     
          ELDI2=ELDIS(IDOFN,INODL+NNOBO)
          DGAPN=DGAPN+SHAPE(INODL,IGAUS)*(ELDI1-ELDI2)*
     .                                                VNORL(IDOFN,IGAUS)
         ELSE
          ELDI2=ELDIS(IDOFN,INODL+NNOBO+IAUXY)
          DGAPN=DGAPN+(SHAPE(INODL,IGAUS)*ELDI1-
     .                 SHAPE(INODL+NNOBO,IGAUS)*ELDI2)*
     .                                                VNORL(IDOFN,IGAUS)
         ENDIF
        END DO
       END DO
C
C**** COMPUTES TANGENTIAL GAP INCREMENT
C
       IF(IFRIC.EQ.1) THEN
        DO IDIME=1,NDIME-1
         DGAPT(IDIME)=0.0D0
         DO IDOFN=1,NDOFN
          DO INODL=1,NNOBO
           ELDI1=DISPL(IDOFN,INODL)+DISIL(IDOFN,INODL)
           IF(NOCOL.EQ.0) THEN
            ELDI2=DISPL(IDOFN,INODL+NNOBO)+DISIL(IDOFN,INODL+NNOBO)
            DGAPT(IDIME)=DGAPT(IDIME)+SHAPE(INODL,IGAUS)*(ELDI1-ELDI2)*
     .                                          VTANL(IDOFN,IDIME,IGAUS)
           ELSE
            ELDI2=DISPL(IDOFN,INODL+NNOBO+IAUXY)+
     .            DISIL(IDOFN,INODL+NNOBO+IAUXY)
            DGAPT(IDIME)=DGAPT(IDIME)+(SHAPE(INODL,IGAUS)*ELDI1-
     .                                 SHAPE(INODL+NNOBO,IGAUS)*ELDI2)*
     .                                          VTANL(IDOFN,IDIME,IGAUS)
           ENDIF
          END DO
         END DO
        END DO
       ENDIF
C
C**** COMPUTES TEMPERATURE
C
       TGAUS=0.0D0
       IF(ITERME.GE.0) THEN              ! coupled problems
        IF(INOTE.EQ.1) THEN
         DO INODL=1,NNOBO
          TGAUS=TGAUS+SHAPE(INODL,IGAUS)*TENOD(INODL)
         END DO
        ELSE
         DO INODL=1,NNOBO
          IF(NOCOL.EQ.0) THEN
           TGAUS=TGAUS+SHAPE(INODL,IGAUS)*TENOD(INODL+NNOBO)
          ELSE
           CALL RUNEND('ERROR: INOTE=2 NOT IMPLEMENTED YET (frin32.f)')
          ENDIF
         END DO
        ENDIF
       ENDIF
C
C**** COMPUTES NORMAL & TANGENTIAL PRESSURE IN THE LOCAL SYSTEM AND
C     COMPUTES FRICTIONAL CONSTITUTIVE TENSOR (see mass32.f)
C
       IF(DGAPN.GT.TOLGAM) DGAPN=TOLGAM
       PRESN=RIGIN*DGAPN                 ! ICOMO=0            ! line
       IF(ICOMO.EQ.1) THEN
        PRESN=RIGIN*DGAPN*DGAPN                               ! parabola
        RIGIN=2.0D0*RIGIN*DABS(DGAPN)
       ENDIF
       IF(ICOMO.EQ.2) THEN
        DGABA=PROPS(8)
        PRESN=RIGIN*(DGAPN-DGABA/2.0D0)                       ! line
        IF(DGAPN.GT.0.0D0.AND.DGAPN.LT.DGABA) THEN
         PRESN=RIGIN/(2.0D0*DGABA)*DGAPN*DGAPN                ! parabola
         RIGIN=RIGIN/DGABA*DGAPN
        ENDIF
       ENDIF
       IF(ICOMO.EQ.3) THEN
        RIGI0=PROPS(9)
        PRESN=RIGIN*DGAPN*DGAPN+RIGI0*DGAPN                   ! parabola
        RIGIN=2.0D0*RIGIN*DABS(DGAPN)+RIGI0
       ENDIF
       IF(ICOMO.EQ.4) THEN
        DGABA=PROPS(8)
        RIGI0=PROPS(9)
        PRESN=RIGIN*(DGAPN-DGABA/2.0D0)+RIGI0/2.0D0*DGABA     ! line
        IF(DGAPN.GT.0.0D0.AND.DGAPN.LT.DGABA) THEN
         PRESN=(RIGIN-RIGI0)/(2.0D0*DGABA)*DGAPN*DGAPN+       ! parabola
     .         RIGI0*DGAPN
         RIGIN=(RIGIN-RIGI0)/DGABA*DGAPN+RIGI0
        ENDIF
       ENDIF
       IF(ICOMO.EQ.5) THEN                     ! mode I debounding model
        IF(DGAPN.GT.0.0D0) THEN                ! compression
         RIGIN=PROPS(2)                        ! K_n_compression
         PRESN=RIGIN*DGAPN
        ELSE                                   ! tension
         RIGIN=PROPS(17)*PROPS(17)/(2.0D0*PROPS(9)*PROPS(8))  ! K_n
         GAPCN=2.0D0*PROPS(8)/PROPS(17)                       ! u^c_n
         GAPNB=PROPS(9)*GAPCN                                 ! u_n-bar
         DAMAN=EHIST(NDIME+NDIME*NDIME,IGAUS)                 ! d_n
         DACRT=PROPS(18)                       ! critical debounding
         IF(DAMAN.LT.DACRT) THEN
          IF(DABS(DGAPN).GT.GAPNB) THEN        ! debounding
           IF(DABS(DGAPN).LE.GAPCN) THEN       ! beginning of debounding
            DAMAX=(DABS(DGAPN)-GAPNB)/DABS(DGAPN)*GAPCN/(GAPCN-GAPNB)
            IF(DAMAX.GT.DAMAN) DAMAN=DAMAX
           ELSE                                ! end of debounding
            DAMAN=DACRT
           ENDIF
          ENDIF
         ENDIF
         PRESN=(1.0D0-DAMAN)*RIGIN*DGAPN
         RIGIN=(1.0D0-DAMAN)*RIGIN
         IF(IITER.GT.0) EHIST(NDIME+NDIME*NDIME,IGAUS)=DAMAN
        ENDIF
       ENDIF
C
       IF(IFRIC.EQ.1) THEN
        DO IDIME=1,NDIME-1
         PREST(IDIME)=EHIST(IDIME,IGAUS)
         IF(DGAPN.GT.0.0D0) THEN
          PREST(IDIME)=PREST(IDIME)+RIGIT
     .                      *(DGAPT(IDIME)-HISTO(IELEM,IGAUS,IDIME))
         ELSE
          PREST(IDIME)=0.0D0     ! tangential pressure is not temp.-dep.
         ENDIF
        END DO
        IF(IELEM.EQ.329)
     .   WRITE(LURES,*)  "INITIAL ITER, PREST(1),IGAUS",PREST(1),IGAUS
C
        PRETN=0.0D0
        DGATN=0.0D0
        DO IDIME=1,NDIME-1
         PRETN=PRETN+PREST(IDIME)*PREST(IDIME)
         DGATN=DGATN+DGAPT(IDIME)*DGAPT(IDIME)
        ENDDO
        IF(PRETN.GT.0.0D0) PRETN=DSQRT(PRETN)
        IF(DGATN.GT.0.0D0) DGATN=DSQRT(DGATN)
C
        DO IDIME=1,NDIME-1
         DO JDIME=1,NDIME-1
          RIGI3(IDIME,JDIME)=0.0D0
          IF(IDIME.EQ.JDIME) RIGI3(IDIME,JDIME)=1.0D0
         ENDDO
        ENDDO
C
        DO IDIME=1,NDIME         ! frictional (stick) const. tensor
         DO JDIME=1,NDIME
          RIGI1(IDIME,JDIME)=0.0D0
          IF(DGAPN.GE.0.0D0.AND.PRETN.GE.0.0D0) THEN  ! GE instead of GT
           DO KDIME=1,NDIME-1
            DO LDIME=1,NDIME-1
             RIGI1(IDIME,JDIME)=RIGI1(IDIME,JDIME)+
     .                                   RIGIT*VTANL(IDIME,KDIME,IGAUS)*
     .                      RIGI3(KDIME,LDIME)*VTANL(JDIME,KDIME,IGAUS)
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ENDDO
C
        IF(ITIME*ISTEP.EQ.1) THEN    ! to avoid initial ill-conditioning
         IF(IITER.EQ.0) THEN
          KDIME=NDIME-1
          DO IDIME=1,NDIME
           DO JDIME=1,NDIME
            KDIME=KDIME+1
            EHIST(KDIME,IGAUS)=RIGI1(IDIME,JDIME)
           ENDDO
          ENDDO
         ENDIF
        ENDIF
       ENDIF
C
C**** TEMPERATURE-DEPENDENT NORMAL PRESSURE
C
       IF(TGAUS.LT.TEMPL) THEN
        IF(DGAPN.LE.0.0D0) THEN
         PRESN=0.0D0
        ENDIF
       ENDIF
C
C**** CORRECTS & STORES TANGENTIAL PRESSURE & FRICTIONAL CONSTITUTIVE
C     TENSOR AS A HISTORY VARIABLES
C
       IF(IITER.GT.0) THEN
        IF(IFRIC.EQ.1) THEN
         IF(IFRIM.EQ.1) THEN
          AFRIC=PRETN-FRIMU*PRESN
          IF(AFRIC.GT.0.0D0.AND.PRESN.GT.0.0D0.AND.PRETN.GT.0.0D0) THEN
           IF(IELEM.EQ.329)
     .      WRITE(LURES,*) "SLIP,DGATN,PRETN,FRIMU*PRESN"
     .                       ,DGATN,PRETN,FRIMU*PRESN
           DO IDIME=1,NDIME-1
            HISTOTEMP(IELEM,IGAUS,IDIME)=DGAPT(IDIME)
            PREST(IDIME)=FRIMU*PRESN*PREST(IDIME)/PRETN
c           PREST(IDIME)=FRIMU*PRESN*DGAPT(IDIME)/DGATN  ! other possib.
           END DO
C
           PRETN=1.0D0           ! redefines norm of tangential pressure
           IF(FRIMU.GT.0.0D0) PRETN=FRIMU*PRESN
C
           DO IDIME=1,NDIME-1    ! frictional (slip) const. tensor
            DO JDIME=1,NDIME
             RIGI2(IDIME,JDIME)=PREST(IDIME)/PRETN*VNORL(JDIME,IGAUS)
c            RIGI2(IDIME,JDIME)=DGAPT(IDIME)/DGATN*VNORL(JDIME,IGAUS)
            END DO
           END DO
           DO IDIME=1,NDIME-1
            DO JDIME=1,NDIME-1
             RIGI3(IDIME,JDIME)=PREST(IDIME)*PREST(JDIME)/(PRETN*PRETN)
            ENDDO
           ENDDO
           DO IDIME=1,NDIME
            DO JDIME=1,NDIME
             DO KDIME=1,NDIME-1
              DO LDIME=1,NDIME-1
               RIGI1(IDIME,JDIME)=RIGI1(IDIME,JDIME)-
     .                                   RIGIT*VTANL(IDIME,KDIME,IGAUS)*
     .                      RIGI3(KDIME,LDIME)*VTANL(JDIME,KDIME,IGAUS)
              ENDDO
             ENDDO
             DO KDIME=1,NDIME-1
              RIGI1(IDIME,JDIME)=RIGI1(IDIME,JDIME)+
     .           FRIMU*RIGIN*VTANL(IDIME,KDIME,IGAUS)*RIGI2(KDIME,JDIME)
             ENDDO
            ENDDO
           ENDDO
          ENDIF
C
          DO IDIME=1,NDIME-1
           EHIST(IDIME,IGAUS)=PREST(IDIME)
          END DO
          KDIME=NDIME-1
          DO IDIME=1,NDIME
           DO JDIME=1,NDIME
            KDIME=KDIME+1
            EHIST(KDIME,IGAUS)=RIGI1(IDIME,JDIME)
           ENDDO
          ENDDO
         ENDIF           ! ifrim.eq.1
        ENDIF            ! ifric.eq.1
        IF(IELEM.EQ.329) THEN
         WRITE(LURES,*) "PREST(1),IGAUS, FINAL ITER", PREST(1),IGAUS
        END IF
       ENDIF             ! iiter.gt.0
C
C**** EVALUATES ELEMENT CONTRIBUT. OF CONTACT FORCE IN THE GLOBAL SYSTEM
C     (AS AN EXTERNAL FORCE >> -, AS A RESIDUAL FORCE >> +)
C
       DO INODE=1,NNOBO
        SAREA=SHAPE(INODE,IGAUS)*DVOLU(IGAUS)
        IEVAB=(INODE-1)*NDIME
C$DIR SCALAR
        DO IDIME=1,NDIME
         PVECT=PRESN*VNORL(IDIME,IGAUS)
         IF(IFRIC.EQ.1) THEN
          DO JDIME=1,NDIME-1
           PVECT=PVECT+PREST(JDIME)*VTANL(IDIME,JDIME,IGAUS)
          END DO
         ENDIF
         IEVAB=IEVAB+1
         BMSIG(IEVAB)=BMSIG(IEVAB)+PVECT*SAREA
        ENDDO
       ENDDO
C
C**** COMPUTES FRICTIONAL HEATING
C
       IF(IFRIC.EQ.1) THEN
        IF(IFRIM.EQ.1) THEN
         IF(ITERME.GT.0) THEN   ! bidirectional coupling
          IF(ITERMF.GT.0) THEN  ! frictional heating
           IF(IGAFR.EQ.1) THEN
            COUTF=0.0D0
            DO IDIME=1,NDIME-1
             COUTF=COUTF+PREST(IDIME)*DGAPT(IDIME)/DTIME*AGAFP
            END DO
            IF(COUTF.LT.0.0D0)
     .       CALL RUNEND('ERROR: FRICTIONAL DISSIPATION < 0')
            DO INODE=1,NNOBO
             SAREA=SHAPE(INODE,IGAUS)*DVOLU(IGAUS)
             PWOEL(INODE)=PWOEL(INODE)+COUTF*SAREA/COUFAC
            END DO
           ENDIF
          ENDIF
         ENDIF
        ENDIF
       ENDIF
C
C**** COMPUTES NODAL NORMAL GAP & PRESSURE TO CONVECTION-RADIATION
C     COEFFICIENT OR JUST TO PRINT (see outgap.f) ONLY 
C
C     Assumption: GAUSS POINTS COINCIDE WITH NODES
C
C     Note: THE COMPUTATION OF MORE CORRECT NODAL NORMAL GAP & PRESSURE
C           FOR NON-COINCIDENT MESHES (NOCOL=1) CAN BE CARRIED OUT
C           THROUGH A SMOOTHING OPERATION (via, for example, smoi32.f;
C           with DGAPN & PRESN being internal variables; this is too
C           complicate to implement!)
C
       IF(NOCOL.EQ.0) THEN                      ! coincident mesh
        INODL=IGAUS
        TGAPL(INODL)=0.0D0                      ! contact
        TGAPL(INODL+NNOBO)=0.0D0
        IF(DGAPN.LT.0.0D0) THEN                 ! normal gap
         TGAPL(INODL)=-DGAPN
         TGAPL(INODL+NNOBO)=-DGAPN
        ENDIF
        PREAL(INODL)=PRESN
        PREAL(INODL+NNOBO)=PRESN
       ELSE                                     ! non-coincident mesh
        TGAPL(IGAUS)=0.0D0                      ! contact
        IF(DGAPN.LT.0.0D0) TGAPL(IGAUS)=-DGAPN  ! normal gap
        PREAL(IGAUS)=PRESN
       ENDIF                         ! nocoi.eq.0
C
  100  CONTINUE
       IF (IELEM.EQ.329) THEN
         WRITE(LURES,*) "FINAL,HISTO,HISTOTEMP",HISTO(IELEM,1,1)
     .                                     ,HISTOTEMP(IELEM,1,1)
         WRITE(LURES,*) " "
       END IF
C
       IF(NOCOL.EQ.0) THEN
        DO IEVAB=1,NEVBO
         BMSIG(IEVAB+NEVBO)=-BMSIG(IEVAB)
        END DO
       ENDIF                         ! nocoi.eq.0
C
      ELSE              ! icoli=2 (linking)
C
       NNOBO=NNODL/2
       NEVBO=NNOBO*NDOFN
C
C**** PROPS(1+IDOFN)=PENALTIES PARAMETERS IN X-Y-Z DIRECTIONS
C     PROPS(5)=LIQUIDUS TEMPERATURE
C
       TEMPL=PROPS(5)
       INOTE=INT(PROPS(6))
C
C**** INITIALISE MECHANICAL COUPLING TERM
C
       IF(ITERME.GT.0) THEN   ! bidirectional coupling
        DO INODL=1,NNODL
         PWOEL(INODL)=0.0D0
        END DO
       ENDIF
C
C**** LOOP OVER INTEGRATION POINTS
C
       DO 101 IGAUS=1,NGAUL
C
C**** COMPUTES TEMPERATURE
C
       TGAUS=0.0D0
       IF(ITERME.GE.0) THEN              ! coupled problems
        IF(INOTE.EQ.1) THEN
         DO INODL=1,NNOBO
          TGAUS=TGAUS+SHAPE(INODL,IGAUS)*TENOD(INODL)
         END DO
        ELSE
         DO INODL=1,NNOBO
          TGAUS=TGAUS+SHAPE(INODL,IGAUS)*TENOD(INODL+NNOBO)
         END DO
        ENDIF
       ENDIF
C
C**** COMPUTES GAPS 
C
       DO IDOFN=1,NDOFN
        DGAPS(IDOFN)=0.0D+00
       ENDDO
C
       IF(TGAUS.GT.TEMPL) THEN
        DO IDOFN=1,NDOFN
         DO INODL=1,NNOBO
          ELDI1=ELDIS(IDOFN,INODL)
          ELDI2=ELDIS(IDOFN,INODL+NNOBO)
          DGAPS(IDOFN)=DGAPS(IDOFN)+SHAPE(INODL,IGAUS)*(ELDI1-ELDI2)
         ENDDO
        ENDDO
       ENDIF            ! tgaus.gt.templ
C
C**** EVALUATE ELEMENT CONTRIBUTION OF CONTACT FORCE
C     (AS A EXTERNAL FORCE >> -, AS A RESIDUAL FORCE >> +)
C
       DO INODE=1,NNOBO
        SAREA=SHAPE(INODE,IGAUS)*DVOLU(IGAUS)
        IEVAB=(INODE-1)*NDIME
C$DIR SCALAR
        DO IDIME=1,NDIME
         PVECT=PROPS(1+IDIME)*DGAPS(IDIME)
         IEVAB=IEVAB+1
         BMSIG(IEVAB)=BMSIG(IEVAB)+PVECT*SAREA
        ENDDO
       ENDDO
C
  101  CONTINUE
C
       DO IEVAB=1,NEVBO
        BMSIG(IEVAB+NEVBO)=-BMSIG(IEVAB)
       END DO
C
      ENDIF                  ! icoli.eq.1
C
      RETURN
      END

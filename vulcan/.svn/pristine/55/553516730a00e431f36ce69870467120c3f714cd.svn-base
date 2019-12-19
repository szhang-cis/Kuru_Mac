      SUBROUTINE FRIN32B(PROPS,ELDIS,BMSIG,DVOLU,SHAPE,EHIST,TENOD,
     .                   ELCOD,PWOEL,PREAL,TGAPL,VNORL)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE CONTACT FORCES FOR ELEMENT 32
C
C     TOTAL AUGMENTED LAGRANGIAN CONTACT METHOD
C
C     Note: for contact and linking problems, NDIME=NDOFN
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
      DIMENSION PROPS(*),       ELDIS(NDOFC,*),
     .          BMSIG(*),       DVOLU(*),
     .          SHAPE(NNODL,*), EHIST(NHIST,*)
      DIMENSION TENOD(*),       ELCOD(NDIME,*),
     .          PWOEL(*),       PREAL(*),
     .          TGAPL(*)
      DIMENSION VNORL(NDIME,*), DGAPS(3)
C
      NNOBO=NNODL/3
      NEVBO=NNOBO*NDOFN
C
C**** PROPERTIES FOR THE CONTACT PROBLEM
C
C     PROPS(2)=RIGIN            : NORNAL STIFFNESS
C     PROPS(3)=IMODEL           : CONTACT MODEL
C     PROPS(5)=TEMPL            : LIQUIDUS TEMPERATURE
C     PROPS(6)=INOTE            : DENOTES WHICH FACE BELONGS TO PIECE
C
      RIGIN=PROPS(2)
      IMODEL=INT(PROPS(3))  ! 0: symmetric matrix; 1: unsymmetric matrix
      TEMPL=PROPS(5)
      INOTE=INT(PROPS(6))
      ICOMO=INT(PROPS(7))
C
C**** INITIALISE MECHANICAL COUPLING TERM (only for friction problems)
C
      IF(ITERME.GT.0) THEN   ! bidirectional coupling
       DO INODL=1,NNODL
        PWOEL(INODL)=0.0
       END DO
      ENDIF
C
C**** LOOP OVER INTEGRATION POINTS
C
      DO 100 IGAUS=1,NGAUL
C
C**** ENSURES A UNIT OUTWARD NORMAL (only necessary when NMEMO2=0)
C
      VNORC=0.0
      DO IDIME=1,NDIME
       VNORA=VNORL(IDIME,IGAUS)
       VNORC=VNORC+VNORA*VNORA
      ENDDO
      VNORC=DSQRT(VNORC)
      DO IDIME=1,NDIME
       VNORA=VNORL(IDIME,IGAUS)
       VNORL(IDIME,IGAUS)=VNORA/VNORC
      ENDDO
C
C**** COMPUTES NORMAL GAP
C
C     DGAPN  > 0: CONTACT PROBLEM
C     DGAPN <= 0: NO CONTACT PROBLEM
C
      DGAPN=0.0
      DO IDOFN=1,NDOFN
       DO INODL=1,NNOBO
        ELDI1=ELDIS(IDOFN,INODL)
        ELDI2=ELDIS(IDOFN,INODL+NNOBO)
        DGAPN=DGAPN+SHAPE(INODL,IGAUS)*(ELDI1-ELDI2)*VNORL(IDOFN,IGAUS)
       END DO
      END DO
C
C**** COMPUTES TEMPERATURE
C
      TGAUS=0.0
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
C**** COMPUTES NORMAL PRESSURE IN THE GLOBAL SYSTEM
C
      PRESN=0.0
      DO INODL=1,NNOBO
       ELDIP=ELDIS(1,INODL+2*NNOBO)
       PRESN=PRESN+SHAPE(INODL,IGAUS)*ELDIP
      END DO
C
      IF(ICOMO.EQ.1) THEN
       DGABA=PROPS(8)
       IF(DGAPN.GT.0.0.AND.DGAPN.LT.DGABA)
     .  PRESN=RIGIN/(2.0*DGABA)*DGAPN*DGAPN
      ENDIF
C
C**** TEMPERATURE-DEPENDENT NORMAL PRESSURE (APPR. KUHN-TUCKER COND.)
C
      IF(TGAUS.LT.TEMPL) THEN
       IF(DGAPN.LE.0.0) THEN
        DGAPN=0.0
        PRESN=0.0
       ENDIF
       IF(PRESN.GE.0.0) THEN
        DGAPN=0.0
        PRESN=0.0
       ENDIF
      ENDIF
C
C**** EVALUATE ELEMENT CONTRIBUTION OF CONTACT FORCE
C     (AS A EXTERNAL FORCE >> -, AS A RESIDUAL FORCE >> +)
C
      DO INODE=1,NNOBO
       SAREA=SHAPE(INODE,IGAUS)*DVOLU(IGAUS)
       IEVAB=(INODE-1)*NDIME
       JEVAB=2*NNOBO*NDIME+IEVAB
C$DIR SCALAR
       DO IDIME=1,NDIME
        PVECT=-PRESN*VNORL(IDIME,IGAUS)
C
        QVECT=0.0
        IF(IDIME.EQ.1) THEN               ! normal component
         IF(IAUGM.EQ.0) THEN
          if(imodel.eq.0) then
           QVECT=-dgapn-presn/rigin
          else
           QVECT= dgapn+presn/rigin
          endif
         ELSE
          if(imodel.eq.0) then
           QVECT=-dgapn
          else
           QVECT= dgapn
          endif
         ENDIF
        ENDIF
C
        IEVAB=IEVAB+1
        JEVAB=JEVAB+1
        BMSIG(IEVAB)=BMSIG(IEVAB)+PVECT*SAREA
        BMSIG(JEVAB)=BMSIG(JEVAB)+QVECT*SAREA
       ENDDO
      ENDDO
C
  100 CONTINUE
C
      DO IEVAB=1,NEVBO
       BMSIG(IEVAB+NEVBO)=-BMSIG(IEVAB)
      END DO
C
C**** COMPUTES NODAL NORMAL GAP & PRESSURE TO CONVECTION-RADIATION
C     COEFFICIENT (ASSUMPTION: VNORL IS THE SAME FOR EVERY GAUSS POINT)
C
      IF(ITERME.GT.0) THEN       ! bidirectional coupling
       DO INODL=1,NNOBO
        DGAPN=0.0
        DO IDOFN=1,NDOFN
         ELDI1=ELDIS(IDOFN,INODL)
         ELDI2=ELDIS(IDOFN,INODL+NNOBO)
         ELDIP=ELDIS(IDOFN,INODL+2*NNOBO)
         DGAPN=DGAPN+(ELDI1-ELDI2)*VNORL(IDOFN,1)
        END DO
C
        PRESN=0.0
        ELDIP=ELDIS(1,INODL+2*NNOBO)
        PRESN=PRESN+SHAPE(INODL,IGAUS)*ELDIP
C
        IF(DGAPN.GE.0.0) THEN      ! contact
         DGAPN=0.0
         PRESN=-PRESN              ! absolute value
        ELSE                       ! a normal gap is produced
         PRESN=0.0
         DGAPN=-DGAPN              ! absolute value
        ENDIF
C
        TGAPL(INODL)=DGAPN
        PREAL(INODL)=PRESN
       END DO
      ENDIF
C
      RETURN
      END

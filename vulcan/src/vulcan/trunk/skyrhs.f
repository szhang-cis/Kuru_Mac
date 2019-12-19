      SUBROUTINE SKYRHS(CSTIF,DISIT,REFOR,LNODS,NDOFN,NELEM,NEVAB,
     .                  NNODE,NPOIN,DISIM,FOREL,IFLAG,
     .                  NPREL,NGRUP,NPROP,NMATS,
     .                  NDATA,NPREV,NSTAT,NMATX,NTOTV,NDIME,
     .                  NDISR,NDISO,
     .                  MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .                  ELVAR,ELMAT,WORK1,TEMPN,DISPR,DTEMP,
     .                  VNORM,DISTO,COORD,INFRI,COFRI,LACTI)
C***********************************************************************
C
C**** THIS ROUTINE MULTIPLIES K * d TO COMPUTE THE RHS
C          IFLAG=0  INITIALIZE & RECOMPUTE RHS (FOR REACTIONS)
C          IFLAG=1  MODIFY RHS     (FOR IMPOSED DISPLACEMENTS)
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
      INCLUDE 'auxl_om.f'
C
      DIMENSION CSTIF(NEVAB,NEVAB), DISIM(*),       DISIT(*),
     .          REFOR(*),           FOREL(*),       LNODS(NNODE,*)
C
      DIMENSION MATNO(NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),       ELPRE(NPREV),
     .          ELVAR(NSTAT),       ELMAT(NMATX),
     .          WORK1(*)
      DIMENSION TEMPN(NPOIN,2),     DISPR(NTOTV,NDISR),
     .          DTEMP(NPOIN),       VNORM(NTOTV)
      DIMENSION DISTO(NTOTV,NDISO), COORD(NDIME,NPOIN)
      DIMENSION INFRI(*),           COFRI(*),
     .          LACTI(NELEM)
C
C**** IF NECESSARY INITIALIZE RHS
C
      IF(IFLAG.EQ.0) THEN
       DO ITOTV=1,NPOIN*NDOFN
        REFOR(ITOTV)=0.0
       ENDDO
      ENDIF
C
      DO 10 IELEM=1,NELEM
C
      NNODL=0
      DO INODE=1,NNODE
       IPOIN=LNODS(INODE,IELEM)
       IF(IPOIN.NE.0) NNODL=NNODL+1
      ENDDO
C
C**** GATHER DISIT > DISIM
C
      CALL GATHER(DISIT,NDOFN,NPOIN,DISIM,NDOFN,NNODL,LNODS(1,IELEM))
C
      JMPOS=0
      IF(IFLAG.EQ.1) THEN
       DO IEVAB=1,NEVAB
        DISIM(IEVAB)=-DISIM(IEVAB)
        IF(DISIM(IEVAB).NE.0.0) JMPOS=1
       ENDDO
      ENDIF
C
C**** IFLAG=0 COMPUTING REACTIONS -- JMPOS=1 PRESCRIBED DISPLACEMENTS
C
      IF(IFLAG.EQ.0.OR.JMPOS.EQ.1) THEN
C
C**** READ CSTIF FROM DATA BASE
C 
       IF(NMEMO7M.EQ.0) THEN
        CALL DATBAS(CSTIF,    6,    2)
       ELSE
        CALL STIFMS(ELDAT,ELPRE,ELVAR,ELMAT,LNODS,MATNO,PROEL,
     .              PROPS,WORK1,TEMPN(1,1),DISPR(1,1),DTEMP,VNORM,COORD,
     .              DISTO(1,1),INFRI,COFRI,LACTI,CSTIF)
       ENDIF
C
C**** MULTIPLY K*d
C
       DO IEVAB=1,NEVAB
        FOREL(IEVAB)=0.0
        DO JEVAB=1,NEVAB
         FOREL(IEVAB)=FOREL(IEVAB)+CSTIF(IEVAB,JEVAB)*DISIM(JEVAB)
        ENDDO
       ENDDO
C
C**** SCATTER FOREL > REFOR
C
       CALL SCATER(FOREL,NDOFN,NNODL,REFOR,NDOFN,NPOIN,LNODS(1,IELEM))
C
      ENDIF
C
   10 CONTINUE
C
      RETURN
      END

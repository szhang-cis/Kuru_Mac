      SUBROUTINE INCREM(ELDAT,ELPRE,ELVAR,ELMAT,HEADS,HTLOD,
     .                  IFFIX,PRESC,LNODS,MATNO,PROEL,PROPS,
     .                  RLOAD,RLOAH,PRESF,PREHF,FICTO,TFICT,TLOAD)
C***********************************************************************
C
C**** THIS ROUTINE INCREMENTS THE APPLIED LOADING
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION MATNO(NELEM),         LNODS(NNODE,NELEM),
     .          PROEL(NPREL,NGRUP),   PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),         ELPRE(NPREV),
     .          ELVAR(NSTAT),         ELMAT(NMATX)
      DIMENSION HEADS(NPOIN,*), 
     .          HTLOD(NHLOD,NSUBF,*), IFFIX(*),       
     .          RLOAD(NTOTV),         TLOAD(*)
      DIMENSION PRESC(NTOTV,2),       RLOAH(NTOTV,NFUNC),
     .          FICTO(NFUNC),         TFICT(NFUNC)
      DIMENSION PRESF(NNODE,NDIME,NELEM),
     .          PREHF(NNODE,NDIME,NELEM,NFUNC)
C
C**** DETERMINE THE LOAD/DISPLACEMENT INCREMENT FACTOR
C
      CALL INCFAC(HTLOD,FICTO,TFICT)
C
C**** INCREMENT THE APPLIED LOAD
C
      CALL INCFOR(HEADS,IFFIX,RLOAD,RLOAH,PRESF,PREHF,FICTO,TLOAD)
C
C**** INCREMENT THE NON-TENSIONAL STRAINS
C
      CALL INCSTN(ELDAT,ELPRE,ELVAR,ELMAT,HTLOD,LNODS,MATNO,
     .            PROEL,PROPS)
C
      RETURN
      END

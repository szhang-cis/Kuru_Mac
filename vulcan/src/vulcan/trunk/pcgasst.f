      SUBROUTINE PCGASST(CSTIF,GSTDI,LNODS,
     .                   NDOFN,NELEM,NEVAB,NNODE,NPOIN,KWIDT,
     .                   NPREL,NGRUP,NPROP,NMATS,
     .                 NDATA,NPREV,NSTAT,NMATX,NTOTV,NDIME,NTOTVM,NFPCH,
     .                   MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .                   ELVAR,ELMAT,WORK1,DISTO,COORD,DISIT,
     .                   ADVEL,TEMPI,PREAS,TGAPS,DISPL,FPCHA,LACTI)
C***********************************************************************
C
C**** THIS ROUTINE ASSEMBLES PART OF THE ELEMENTAL STIFFNESS MATRIX 
C     INTO THE GLOBAL PRECONDITIONING MATRIX
C
C.... INPUT PARAMETERS
C
C     CSTIF(NEVAB,NEVAB)      - SQUARE ARRAY CONTAINING ELEMENT MATRIX
C     KWIDT                   - FLAG FOR ASSEMBLY,
C                                   = 0   ONLY DIAGONAL
C                                   < 0   BLOCK DIAGONAL
C
C.... OUPUT PARAMETERS
C
C     GSTDI(*)                - PRECONDITIONING MATRIX
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
      INCLUDE 'auxl_omt.f'
C
      DIMENSION CSTIF(NEVAB,NEVAB), GSTDI(*), LNODS(NNODE,*)
C
      DIMENSION MATNO(NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),       ELPRE(NPREV),
     .          ELVAR(NSTAT),       ELMAT(NMATX)
      DIMENSION DISTO(NTOTV,3),     COORD(NDIME,NPOIN),
     .          DISIT(NTOTV)
      DIMENSION ADVEL(NTOTV*NDIME), TEMPI(NPOIN,2)
      DIMENSION PREAS(NPOIN),       TGAPS(NPOIN)
      DIMENSION DISPL(NTOTVM),      FPCHA(NFPCH,NPOIN),
     .          LACTI(NELEM)
C
C**** INITIALISE GLOBAL PRECONDITIONING MATRIX  
C
      NSIZE=NDOFN
      IF(KWIDT.LT.0) NSIZE=NDOFN*NDOFN
      NNDEX=NPOIN*NSIZE
      DO 10 INDEX=1,NNDEX
   10 GSTDI(INDEX)=0.0
C
C**** ASSEMBLE GLOBAL PRECONDITIONING MATRIX 
C
      DO 1000 IELEMT=1,NELEM
C
      NNODL=0
      DO INODE=1,NNODE
       IPOIN=LNODS(INODE,IELEMT)
       IF(IPOIN.NE.0) NNODL=NNODL+1
      ENDDO
C
C**** READ CSTIF FROM DATA BASE
C 
      IF(NMEMO7.EQ.0) THEN
       CALL DATBAST(CSTIF,    6,    2)
      ELSE
       CALL STIFMST(ELDAT,ELPRE,ELVAR,ELMAT,LNODS,MATNO,PROEL,
     .              PROPS,WORK1,DISTO(1,2),DISIT,COORD,
     .              ADVEL,TEMPI(1,1),PREAS,TGAPS,DISTO(1,1),
     .              TEMPI(1,2),DISPL,FPCHA,LACTI,CSTIF)
      ENDIF
C
      IF(KWIDT.EQ.0) THEN               ! Diagonal Preconditioning
       DO INODE=1,NNODL
        IPOIN=LNODS(INODE,IELEMT)
        IEVAB=(INODE-1)*NDOFN
        IPOSN=(IPOIN-1)*NSIZE
C$DIR SCALAR
        DO IDOFN=1,NDOFN
         IEVAB=IEVAB+1
         IPOSN=IPOSN+1 
         GSTDI(IPOSN)=GSTDI(IPOSN)+CSTIF(IEVAB,IEVAB)
        ENDDO
       ENDDO
      ELSE                              ! Block-Diagonal Preconditioning
       DO INODE=1,NNODL
        IPOIN=LNODS(INODE,IELEMT)
        IEVA0=(INODE-1)*NDOFN
        IPOSN=(IPOIN-1)*NSIZE
C$DIR SCALAR
        DO JEVAB=IEVA0+1,IEVA0+NDOFN
C$DIR SCALAR
         DO IEVAB=IEVA0+1,IEVA0+NDOFN
          IPOSN=IPOSN+1 
          GSTDI(IPOSN)=GSTDI(IPOSN)+CSTIF(IEVAB,JEVAB)
         ENDDO
        ENDDO
       ENDDO   
      ENDIF
C
 1000 CONTINUE
C
      RETURN
      END
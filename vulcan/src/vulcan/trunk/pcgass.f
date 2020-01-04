      SUBROUTINE PCGASS(CSTIF,GSTDI,LNODS,
     .                  NDOFN,NELEM,NEVAB,NNODE,NPOIN,KWIDT,
     .                  NPREL,NGRUP,NPROP,NMATS,
     .                  NDATA,NPREV,NSTAT,NMATX,NTOTV,NDIME,
     .                  NDISR,NDISO,
     .                  MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .                  ELVAR,ELMAT,WORK1,TEMPN,DISPR,DTEMP,
     .                  VNORM,DISTO,COORD,INFRI,COFRI,LACTI)
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
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
C
      DIMENSION CSTIF(NEVAB,NEVAB), GSTDI(*), LNODS(NNODE,*)
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
      DO 1000 IELEM=1,NELEM
C
      NNODL=0
      DO INODE=1,NNODE
       IPOIN=LNODS(INODE,IELEM)
       IF(IPOIN.NE.0) NNODL=NNODL+1
      ENDDO
C
C**** READ CSTIF FROM DATA BASE
C 
      IF(NMEMO7M.EQ.0) THEN
       CALL DATBAS(CSTIF,    6,    2)
      ELSE
       CALL STIFMS(ELDAT,ELPRE,ELVAR,ELMAT,LNODS,MATNO,PROEL,
     .             PROPS,WORK1,TEMPN(1,1),DISPR(1,1),DTEMP,VNORM,COORD,
     .             DISTO(1,1),INFRI,COFRI,LACTI,CSTIF)
      ENDIF
C
      IF(KWIDT.EQ.0) THEN               ! Diagonal Preconditioning
       DO INODE=1,NNODL
        IPOIN=LNODS(INODE,IELEM)
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
        IPOIN=LNODS(INODE,IELEM)
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
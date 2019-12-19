      SUBROUTINE LINCOM(BSBAR,DEPSV,DESIG,DMATX,DSTRA,ELDIS,GPCOD,
     .                  LARGE,NDIME,NDOFN,NNODE,NSTR1,PROPS,SHAPE,
     .                  STRAN,STRA0,TSTRA,XJACM,SGTOT)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE INCOMPRESSIBILITY CONDITION
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
      DIMENSION BSBAR(NSTR1,*), DMATX(NSTRS,*), DESIG(*), DSTRA(*),
     .          ELDIS(NDIME,*), PROPS(*),       GPCOD(*), SHAPE(*),
     .          STRAN(*),       STRA0(*),       TSTRA(*), XJACM(NDIME,*)
      DIMENSION SGTOT(*)
C
C**** CALCULATES THE INCOMPRESSIBILITY CONDITION
C
      DO ISTR1=1,NSTR1
       TSTRA(ISTR1)=0.0D0
       DO INODE=1,NNODE
        DO IDIME=1,NDIME
         IEVAB=(INODE-1)*NDIME+IDIME
         TSTRA(ISTR1)=TSTRA(ISTR1)+BSBAR(ISTR1,IEVAB)*
     .                ELDIS(IDIME,INODE)
        ENDDO
       ENDDO
      ENDDO
C
      DO ISTR1=1,NSTR1
       SGTOT(ISTR1)=0.0
      ENDDO
      DO ISTRS=1,NSTRS
       DO JSTRS=1,NSTRS
        SGTOT(ISTRS)=SGTOT(ISTRS)+DMATX(ISTRS,JSTRS)*
     .               TSTRA(JSTRS)
       ENDDO
      ENDDO
C
      RETURN
      END

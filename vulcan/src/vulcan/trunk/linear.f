      SUBROUTINE LINEAR(CARTD,DEPSV,DESIG,DMATX,DSTRA,ELDIS,GPCOD,
     .                  LARGE,NDIME,NDOFN,NNODE,NSTR1,PROPS,SHAPE,
     .                  STRAN,STRA0,TSTRA,
     .                  XJACM,XJACI,XJA3M,XJA3I,DETJM,DETJB,
     .                  ELCOD,
     .                  ITASL)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES STRAINS AND INCREMENTAL ELASTIC STRESSES
C
C     Notes:
C
C     THE JACOBIAN MATRIX IS:
C
C     THE DISPLACEMENTS GRADIENT FOR LARGE=0
C     (ELDIS=DISPLACEMENTS)
C
C
C     THE DEFORMATION GRADIENT FOR LARGE > 0
C
C     FOR LARGE=1 (XJACM = MATERIAL DISPLACEMENTS GRADIENT + IDENTITY)
C
C     FOR LARGE=2 (XJACM = UPDATED INCREMENTAL DISPLACEMENTS GRADIENT +
C                          IDENTITY)
C
C     IN BOTH CASES:
C     ELDIS = SPATIAL COORDINATES = MATERIAL COORDINATES + DISPLACEMENTS
C
C     LARGE flag independent computation of TSTRA will be attempted for
C     large strains
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
      DIMENSION CARTD(NDIME,*), DMATX(NSTRS,*), DESIG(*), DSTRA(*),
     .          ELDIS(NDIME,*), PROPS(*),       GPCOD(*), SHAPE(*),
     .          STRAN(*),       STRA0(*),       TSTRA(*),
     .          XJACM(NDIME,*), XJACI(NDIME,*)
      DIMENSION ELCOD(NDIME,*)
C
C**** CALCULATE THE JACOBIAN MATRIX
C
      CALL PROMA2(XJACM,ELDIS,CARTD,NDIME,NNODE)
C
C**** CALCULATE TOTAL STRAINS
C
      IF(LARGE.EQ.0) THEN
       GO TO (1,2,3,4,5), NTYPE
C
    1  TSTRA(1)=XJACM(1,1)
       TSTRA(2)=XJACM(2,2)
       TSTRA(3)=XJACM(1,2)+XJACM(2,1)
       TSTRA(4)=0.0             ! this component is computed in planes.f
       GOTO 6
    2  TSTRA(1)=XJACM(1,1)
       TSTRA(2)=XJACM(2,2)
       TSTRA(3)=XJACM(1,2)+XJACM(2,1)
       TSTRA(4)=0.0
       GOTO 6
    3  TSTRA(1)=XJACM(1,1)
       TSTRA(2)=XJACM(2,2)
       TSTRA(3)=XJACM(1,2)+XJACM(2,1)
       TSTRA(4)=0.0
       DO INODE=1,NNODE
        TSTRA(4)=TSTRA(4)+ELDIS(1,INODE)*SHAPE(INODE)/GPCOD(1)
       ENDDO
       GOTO 6
    4  TSTRA(1)=XJACM(1,1)
       TSTRA(2)=XJACM(2,2)
       TSTRA(3)=XJACM(1,2)+XJACM(2,1)
       TSTRA(4)=XJACM(3,3)
       TSTRA(5)=XJACM(1,3)+XJACM(3,1)
       TSTRA(6)=XJACM(2,3)+XJACM(3,2)
       GOTO 6
    5  TSTRA(1)=XJACM(1,1)
    6  CONTINUE
C
      ELSE                  ! large strains
C
       IF(LARGE.EQ.1) THEN                         ! large strains (TLF)
C
C**** CALCULATE THE DETERMINANT AND INVERS OF JACOBIAN MATRIX
C
        CALL INVMTX(XJACM,XJACI,DETJM,NDIME)
C
        XJA3M=1.0
        XJA3I=1.0
        IF(NTYPE.EQ.3) THEN
         GPCOI=0.0
         GPCOA=0.0D0                 ! average radial spatial coordinate
         DO INODE=1,NNODE
          GPCOI=GPCOI+ELDIS(1,INODE)*SHAPE(INODE)
          GPCOA=GPCOA+ELDIS(1,INODE)
         ENDDO
         GPCOA=GPCOA/NNODE
         IF(GPCOD(1).GT.(1.0D-10*GPCOA).AND.
     .         GPCOI.GT.(1.0D-10*GPCOA)) THEN         ! see GPCOA in *.f
          DETJM=DETJM*GPCOI/GPCOD(1)
          XJA3M=GPCOI/GPCOD(1)
          XJA3I=1.0/XJA3M
         ENDIF
        ENDIF
        DETJB=DETJM
C
        IF(ITASL.EQ.1) RETURN
C
        GO TO (11,12,13,14,15), NTYPE
C
   11   TSTRA(1)=0.5*(XJACM(1,1)*XJACM(1,1)+XJACM(2,1)*XJACM(2,1)-1.0)
        TSTRA(2)=0.5*(XJACM(1,2)*XJACM(1,2)+XJACM(2,2)*XJACM(2,2)-1.0)
        TSTRA(3)=     XJACM(1,1)*XJACM(1,2)+XJACM(2,1)*XJACM(2,2)
        TSTRA(4)=0.0            ! this component is computed in planes.f
        GO TO 7
   12   TSTRA(1)=0.5*(XJACM(1,1)*XJACM(1,1)+XJACM(2,1)*XJACM(2,1)-1.0)
        TSTRA(2)=0.5*(XJACM(1,2)*XJACM(1,2)+XJACM(2,2)*XJACM(2,2)-1.0)
        TSTRA(3)=     XJACM(1,1)*XJACM(1,2)+XJACM(2,1)*XJACM(2,2)
        TSTRA(4)=0.0
        GO TO 7
   13   TSTRA(1)=0.5*(XJACM(1,1)*XJACM(1,1)+XJACM(2,1)*XJACM(2,1)-1.0)
        TSTRA(2)=0.5*(XJACM(1,2)*XJACM(1,2)+XJACM(2,2)*XJACM(2,2)-1.0)
        TSTRA(3)=     XJACM(1,1)*XJACM(1,2)+XJACM(2,1)*XJACM(2,2)
        TSTRA(4)=0.0
        DO INODE=1,NNODE
         ELDCO=ELDIS(1,INODE)-ELCOD(1,INODE)       ! total displacements
         TSTRA(4)=TSTRA(4)+ELDCO*SHAPE(INODE)/GPCOD(1)
        ENDDO
        TSTRA(4)=TSTRA(4)+0.5*TSTRA(4)*TSTRA(4)
        GO TO 7
   14   TSTRA(1)=0.5*(XJACM(1,1)*XJACM(1,1)+XJACM(2,1)*XJACM(2,1)
     .               +XJACM(3,1)*XJACM(3,1)-1.0)
        TSTRA(2)=0.5*(XJACM(1,2)*XJACM(1,2)+XJACM(2,2)*XJACM(2,2)
     .               +XJACM(3,2)*XJACM(3,2)-1.0)
        TSTRA(3)=     XJACM(1,1)*XJACM(1,2)+XJACM(2,1)*XJACM(2,2)
     .               +XJACM(3,1)*XJACM(3,2)
        TSTRA(4)=0.5*(XJACM(1,3)*XJACM(1,3)+XJACM(2,3)*XJACM(2,3)
     .               +XJACM(3,3)*XJACM(3,3)-1.0)
        TSTRA(5)=     XJACM(1,1)*XJACM(1,3)+XJACM(2,1)*XJACM(2,3)
     .               +XJACM(3,1)*XJACM(3,3)
        TSTRA(6)=     XJACM(1,2)*XJACM(1,3)+XJACM(2,2)*XJACM(2,3)
     .               +XJACM(3,2)*XJACM(3,3)
        GO TO 7
   15   TSTRA(1)=0.5*(XJACM(1,1)*XJACM(1,1)-1.0)
    7   CONTINUE
       ENDIF                ! large.eq.1
C
       IF(LARGE.EQ.2) THEN                         ! large strains (ULF)
c       GO TO (21,22,23,24,25), NTYPE
        CALL RUNEND('LINEAR: LARGE=2 NOT IMPLEMENTED')
       ENDIF                ! large.eq.2
C
      ENDIF                 ! large.eq.0
C
C**** AND THE CORRESPONDING INCREMENTAL STRAINS, NEW TOTAL STRAINS
C
      IF(NMEMOM.EQ.0) THEN
       DO ISTRE=1,NSTR1
        DSTRA(ISTRE)=TSTRA(ISTRE)-STRAN(ISTRE)    ! INCR. STRAINS
        STRAN(ISTRE)=TSTRA(ISTRE)                 ! UPDATE TOTAL STRAINS
       ENDDO
      ELSE
       DO ISTRE=1,NSTR1
        DSTRA(ISTRE)=TSTRA(ISTRE)-STRAN(ISTRE)
     .                           -STRA0(ISTRE)    ! INCR. STRAINS
        STRAN(ISTRE)=TSTRA(ISTRE)-STRA0(ISTRE)    ! UPDATE TOTAL STRAINS
        TSTRA(ISTRE)=STRAN(ISTRE)               
       ENDDO
      ENDIF
C
      DEPSV=DSTRA(1)+DSTRA(2)+DSTRA(4)
C
C**** AND THE INCREMENTAL ELASTIC STRESSES 
C
      DO ISTRE=1,NSTR1
       DESIG(ISTRE)=0.0D00
      ENDDO
      IF(NMEMOM.EQ.1) THEN
       DO ISTRE=1,NSTRS
        DO JSTRE=1,NSTRS
         DESIG(ISTRE)=DESIG(ISTRE)+DMATX(ISTRE,JSTRE)*DSTRA(JSTRE)
        ENDDO
       ENDDO
      ENDIF
C
      RETURN
      END

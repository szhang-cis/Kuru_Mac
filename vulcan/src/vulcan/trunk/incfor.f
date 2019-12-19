      SUBROUTINE INCFOR(HEADS,IFFIX,RLOAD,RLOAH,PRESF,PREHF,FICTO,TLOAD)
C***********************************************************************
C
C**** THIS ROUTINE INCREMENTS THE APPLIED FORCE
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
C
      DIMENSION HEADS(NPOIN,*),     IFFIX(*),     RLOAD(NTOTV),
     .          RLOAH(NTOTV,NFUNC), FICTO(NFUNC), TLOAD(*)
      DIMENSION PRESF(NNODE,NDIME,NELEM),
     .          PREHF(NNODE,NDIME,NELEM,NFUNC)
C
      DO IDOFN=1,NDOFN
C$DIR NO_RECURRENCE
       DO IPOIN=1,NPOIN
        ITOTV=(IPOIN-1)*NDOFC+IDOFN
        IF(IFFIX(ITOTV).EQ.0) THEN    ! load in reaction is not included
         DO IFUNC=1,NFUNC
          RLOA1=RLOAH(ITOTV,IFUNC)*FICTO(IFUNC)
          TLOAD(ITOTV)=TLOAD(ITOTV)+RLOA1
         ENDDO
        ENDIF
       ENDDO
      ENDDO
C
      IF(NLDSF.EQ.1) THEN             ! deformation-dependent load
       DO IELEM=1,NELEM
        DO INODE=1,NNODE
         DO IDIME=1,NDIME
          DO IFUNC=1,NFUNC
           PRESF(INODE,IDIME,IELEM)=PRESF(INODE,IDIME,IELEM)+
     .                       PREHF(INODE,IDIME,IELEM,IFUNC)*FICTO(IFUNC)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDIF
C
      IF(KPORE.EQ.2)THEN
       DO IPOIN=1,NPOIN
        HEADS(IPOIN,3)=HEADS(IPOIN,4)
        HEADS(IPOIN,4)=0.0
        TLOAD(IPOIN*NDOFC)=0.0
       ENDDO
      ENDIF
C
      RETURN
      END

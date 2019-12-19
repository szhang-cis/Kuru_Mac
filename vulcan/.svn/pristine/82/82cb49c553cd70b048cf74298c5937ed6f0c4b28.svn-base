      SUBROUTINE RESLOD(DISIT,HEADS,IFFIX,REFOR,TLOAD)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE RESIDUAL FORCES, REACTIONS & GCURN
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION DISIT(*), HEADS(NPOIN,*), IFFIX(*), REFOR(*),
     .          TLOAD(NTOTV,*)
C
      GCURN=0.
      DO IDOFN=1,NDOFN
C$DIR NO_RECURRENCE
       DO IPOIN=1,NPOIN
        ITOTV=(IPOIN-1)*NDOFC+IDOFN
C
        IF(KINTE.EQ.0.OR.KINTE.EQ.2)
     .   REFOR(ITOTV)=TLOAD(ITOTV,1)-REFOR(ITOTV)
C
        IF(KINTE.EQ.1) THEN                  ! Wilson
         CONST=3.0D0/TBETA                   ! theta
         REFOR(ITOTV)=CONST*TLOAD(ITOTV,1)+(1.0D0-CONST)*TLOAD(ITOTV,2)-
     .                                                      REFOR(ITOTV)
        ENDIF
C
        IF(KINTE.EQ.3)                       ! Central difference
     .   REFOR(ITOTV)=TLOAD(ITOTV,2)-REFOR(ITOTV)
C
        IF(KINTE.EQ.5) THEN                  ! HHT
         CNEW3=TGAMA-1.0D0                   ! gama (alfa-HHT)
         REFOR(ITOTV)=(1.0D0+CNEW3)*TLOAD(ITOTV,1)-CNEW3*TLOAD(ITOTV,2)-
     .                                                      REFOR(ITOTV)
        ENDIF
C
        IF(KINTE.EQ.6) THEN                  ! Collocation
         CNEW3=TDELT                         ! theta
         REFOR(ITOTV)=CNEW3*TLOAD(ITOTV,1)+(1.0D0-CNEW3)*TLOAD(ITOTV,2)-
     .                                                      REFOR(ITOTV)
        ENDIF
C
        GCURN=GCURN+REFOR(ITOTV)*DISIT(ITOTV)
       ENDDO
      ENDDO
C
      IF(KPORE.EQ.2) THEN
       CONST=WLUMP*DTIME
       IF(KINTE.EQ.2) CONST=CONST/TBETA
       DO IPOIN=1,NPOIN
        TSMUS=HEADS(IPOIN,4)
        IF(KINTE.EQ.2) TSMUS=TSMUS+(1.-TBETA)*HEADS(IPOIN,3)
        ITOTV=IPOIN*NDOFC
        IF(IFFIX(ITOTV).NE.0) TLOAD(ITOTV,1)=REFOR(ITOTV)-TSMUS
        REFOR(ITOTV)=CONST*(TLOAD(ITOTV,1)+TSMUS-REFOR(ITOTV))
       ENDDO
      ENDIF
C
      RETURN
      END

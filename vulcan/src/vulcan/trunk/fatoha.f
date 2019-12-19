      SUBROUTINE FATOHA(PREYS,YIELD,CCERO,
     .                  FREDU,ANUMF,STENI,STEN1,STEN2,STMIN,STMAX,REVEF)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE REDUCTION FUNCTION FOR FATIGUE PROBLEMS
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
      INCLUDE 'auxl_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'prob_om.f'
C
      COMMON/FATIPROE/PFATI,PCYCL,PREVF,REVCO
C
      KFATI=INT(PFATI)
      IF(KFATI.EQ.0) RETURN
C
      IF(IFATI.EQ.0) RETURN
C
C**** COMPUTE FATIGUE VARIABLES
C
      NCYCL=INT(PCYCL)
      IF(NCYCL.EQ.0) THEN     ! compute reversibility factor R (IREVF=1)
       call runend('min & max (R) not implemented yet')
c                             ! assigment of steni,sten1,sten2,stmin,
c                             !              stmax,revef
c      REVEF=STMIN/STMAX
      ELSE                    ! compute f_red
C
       IREVF=INT(PREVF)
       IF(IREVF.EQ.0) THEN    ! constant R
        REVER=REVCO
       ELSE                   ! variable R
c       REVER=REVEF
       ENDIF
C
       IF(IFATM.EQ.1) THEN
        IF(REVER.LT.VFATI(6)) REVER=VFATI(6)
        IF(REVER.GT.VFATI(7)) REVER=VFATI(7)
        SLIMI=CCERO*(VFATI(2)+VFATI(3)*REVER)
        IF(YIELD.LT.SLIMI) THEN                     ! yield < sigma_lim
         SINFI=PREYS-YIELD*(PREYS-SLIMI)/SLIMI
         IF(PCYCL.LT.VFATI(1)) THEN                 ! n_cycles < limit n
          ALFAFAT=VFATI(4)*(VFATI(5)**(2.0D0*REVER))
          BETAFAT=VFATI(8)*ALFAFAT
          AUXFA1=(DLOG10(VFATI(1)))**BETAFAT
          BBFAT=-DLOG(SINFI/PREYS)/AUXFA1
          AUXFA2=(DLOG10(PCYCL))**BETAFAT
          FREDA=DEXP(-BBFAT*AUXFA2)
         ELSE                                       ! n_cycles > limit n
          FREDA=SINFI/PREYS
         ENDIF
        ELSE                                        ! yield > sigma_lim
         ALFAFAT=VFATI(4)*(VFATI(5)**(2.0D0*REVER))
         BETAFAT=VFATI(8)*ALFAFAT
         AUXFA0=(DLOG(YIELD/PREYS)/DLOG(SILIM))**(1.0D0/ALFAFAT)
         ANUMX=10.0D0**(DLOG10(VFATI(1))*AUXFA0)
         IF(ANUMF.EQ.0.0) ANUMF=ANUMX
         IF(PCYCL.LT.ANUMX) ANUMF=ANUMX
         AUXFA1=(DLOG10(ANUMF))**BETAFAT
         BBFAT=-DLOG(YIELD/PREYS)/AUXFA1
         AUXFA2=(DLOG10(PCYCL))**BETAFAT
         FREDA=DEXP(-BBFAT*AUXFA2)
        ENDIF
c       IF(FREDA.LT.FREDU) FREDU=FREDA
        fredu=freda   ! to be removed !!!!
        PREYS=PREYS*FREDU
       ENDIF                          ! ifatm.eq.1
C
       IF(IFATM.EQ.2) THEN
        call runend('ERROR: model 2 not implemented yet')
       ENDIF                          ! ifatm.eq.2
C
      ENDIF                           ! ncycl.eq.0
C

      write(lures,*) 'ccero,slimi,yield,preys',ccero,slimi,yield,preys

C
      RETURN
      END

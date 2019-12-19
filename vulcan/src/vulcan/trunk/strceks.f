      SUBROUTINE STRCEKS
C***********************************************************************
C
C**** THIS ROUTINE CHECKS THE CONTROLLING PARAMETERS FOR THE
C     CURRENT TIME INTERVAL
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
C
      DIMENSION NERORS(17)
C
      DO 10 IERORS=1,17
   10 NERORS(IERORS)=0
C
C**** CHECK ON PARAMETERS FOR STATIC-UNCOUPLED PROBLEM
C
      IF(KDYNAS.EQ.0.AND.KPORES.NE.2) THEN
       IF(KINTES.NE.0)                                    NERORS(2)=1
       IF(KOPTIS.NE.0.AND.LAUTOS.NE.0)                    NERORS(3)=1
       IF(KARCLS.EQ.4.AND.NCDISS.EQ.0)                    NERORS(9)=1
       IF(KSMUSS.GT.0)                                    NERORS(10)=1
       IF(KARCLS.NE.0.AND.LAUTOS.EQ.0)                    NERORS(12)=1
      ENDIF
C
C**** CHECK ON PARAMETERS FOR STATIC-COUPLED PROBLEM
C
      IF(KDYNAS.EQ.0.AND.KPORES.EQ.2) THEN
       IF(KINTES.EQ.0)                                    NERORS(4)=1
       IF(LAUTOS.NE.0)                                    NERORS(5)=1
       IF(KSMUSS.GT.0.AND.LINESS.NE.0)                    NERORS(6)=1 !?
       IF(LINESS.NE.0) WRITE(LUREST,800)
       IF(KARCLS.NE.0)                                    NERORS(7)=1
       IF(LACCES.NE.0)                                    NERORS(8)=1
      ENDIF
C
C**** CHECK ON PARAMETERS FOR DYNAMIC-UNCOUPLED PROBLEM
C
      IF(KDYNAS.EQ.1.AND.KPORES.NE.2) THEN
       IF(KINTES.EQ.0)                                    NERORS(4)=1
       IF(LAUTOS.NE.0)                                    NERORS(5)=1
       IF(KARCLS.NE.0)                                    NERORS(7)=1
       IF(LACCES.NE.0.AND.KPROBS.NE.4)                    NERORS(8)=1
       IF(KSMUSS.GT.0)                                    NERORS(10)=1
       IF(KINTES.NE.2.AND.KPROBS.NE.4)                    NERORS(11)=1
       IF(KOPTIS.NE.0)                                    NERORS(13)=1
      ENDIF
C
C**** CHECK ON PARAMETERS FOR DYNAMIC-COUPLED PROBLEM
C
      IF(KDYNAS.EQ.1.AND.KPORES.EQ.2) THEN
       IF(LAUTOS.NE.0)                                    NERORS(5)=1
       IF(KSMUSS.GT.0.AND.LINESS.NE.0)                    NERORS(6)=1
       IF(LINESS.NE.0) WRITE(LURESS,800)
       IF(KARCLS.GT.0)                                    NERORS(7)=1
       IF(LACCES.NE.0)                                    NERORS(8)=1
       IF(KINTES.NE.2)                                    NERORS(11)=1
       IF(KOPTIS.NE.0)                                    NERORS(13)=1
      ENDIF
C
C**** CHECK ON PARAMETERS FOR THE "EXPLICIT_SOLUTION" OPTION
C
      IF(KSOLVS.EQ.4) THEN
       IF(KDYNAS.EQ.0)                                    NERORS(14)=1
       IF(KDYNAS.EQ.1) THEN
        IF(KINTES.EQ.1) THEN
         IF(TALFAS.NE.0.0)                                NERORS(15)=1
        ELSE
                                                          NERORS(16)=1
        ENDIF
       ENDIF
       IF(ICONVS.EQ.1) THEN
        IF(IPERTS.EQ.3)                                   NERORS(17)=1
       ENDIF
      ENDIF
C
C**** CHECK IF ERRORS HAD OCCURED
C
      KOUNTS=0
      DO 20 IERORS=1,17
      IF(NERORS(IERORS).NE.1)GOTO 20
      KOUNTS=KOUNTS+1
   20 CONTINUE
C
      IF(KOUNTS.EQ.0) RETURN
C
C**** ERRORS DIAGNOSED BY ROUTINE STRCEK
C
      WRITE(LURESS,600)KOUNTS
  600 FORMAT(//,20X,I5,' INPUT ERRORS DETECTED BY STRATE.'/,
     .          25X,   ' DIAGNOSTIC FOLLOWS :',//)
C
      IF(NERORS(2).EQ.1)  WRITE(LURESS,705) 
      IF(NERORS(3).EQ.1)  WRITE(LURESS,710) 
      IF(NERORS(4).EQ.1)  WRITE(LURESS,715) 
      IF(NERORS(5).EQ.1)  WRITE(LURESS,720)
      IF(NERORS(6).EQ.1)  WRITE(LURESS,725) 
      IF(NERORS(7).EQ.1)  WRITE(LURESS,730) 
      IF(NERORS(8).EQ.1)  WRITE(LURESS,735)
      IF(NERORS(9).EQ.1)  WRITE(LURESS,740)
      IF(NERORS(10).EQ.1) WRITE(LURESS,745)
      IF(NERORS(11).EQ.1) WRITE(LURESS,750)
      IF(NERORS(12).EQ.1) WRITE(LURESS,755)
      IF(NERORS(13).EQ.1) WRITE(LURESS,760)
      IF(NERORS(14).EQ.1) WRITE(LURESS,765)
      IF(NERORS(15).EQ.1) WRITE(LURESS,770)
      IF(NERORS(16).EQ.1) WRITE(LURESS,775)
      IF(NERORS(17).EQ.1) WRITE(LURESS,776)
C
      CALL RUNENDS('    ERROR DETECTED IN STRATET      ')
C
      RETURN
C
  705 FORMAT(
     .' THE PROBLEM IS STATIC-UNCOUPLED AND YOU ARE ASKING FOR'
     .' A TIME STEPPING ALGORITHM',/)
C
  710 FORMAT(
     .' YOU ARE ASKING SIMULTANEOUSLY FOR AUTOMATIC LOAD AND TIME'
     .' INCREMENTATION',/,' >>>> INCOMPATIBLE',/)
C
  715 FORMAT(
     .' THE PROBLEM ADVANCES IN TIME BUT YOU ARE NOT ASKING FOR',/,
     .' ANY TIME STEPPING ALGORITHM OR YOU ARE ASKING FOR A ',/,
     .' NON-AVAILABLE ONE',/)
C
  720 FORMAT(
     .' AUTOMATIC LOAD INCREMENTATION IS INCOMPATIBLE WITH TIME',/,
     .' STEPPING ALGORITHMS',/)
C
  725 FORMAT(
     .' SMOOTHING AND LINE SEARCH ARE NOT AVAILABLE SIMULTANEOUSLY'/)
C
  730 FORMAT(
     .' ARC-LENGTH IS NOT COMPATIBLE WITH TIME STEPPING PROBLEMS',/)
C
  735 FORMAT(
     .' CONVERGENCE ACCELERATOR IS NOT COMPATIBLE WITH TIME STEPPING',
     .' PROBLEMS',/)
C
  740 FORMAT(
     .' YOU HAVE SELECTED DISPLACEMENT CONTROL BUT HAVE NOT SPECIFIED',
     .' A CONTROLLING D.O.F.',/)
C
  745 FORMAT(
     .' THIS IS NOT A COUPLED PROBLEM AND YOU ARE ACTIVATING THE',
     .' SMOOTHING ALGORITHM',/)
C
  750 FORMAT(
     .' THE PROBLEM IS SECOND ORDER IN TIME AND YOU ARE NOT SELECTING'/,
     .' AN APPROPRIATE TIME STEPPING ALGORITHM',/)
C
  755 FORMAT(
     .' THE ARC-LENGTH ALGORITHM REQUIRES AUTOMATIC LOAD '/,
     .' PLEASE, SPECIFY ONE OF THE AVAILABLE OPTIONS',/)
C
  760 FORMAT(
     .' FOR THE DYNAMIC PROBLEM IT IS SUPPOSED TO HAVE CONSTANT TIME'/,
     .' WITHIN EACH INTERVAL. THUS IT SHOULD SET THE TIME PARAMETERS'/,
     .' IN CARD "TIME_INCREMENTION,<PARA1>,<PARA2>" EQUAL TO ZERO'/)
C
  765 FORMAT(
     .'THE "EXPLICIT_SOLUTION" OPTION CAN ONLY BE USED FOR TRANSIENT',/,
     .' PROBLEMS',/)
C
  770 FORMAT(
     .' THE "EXPLICIT_SOLUTION" OPTION MUST BE USED WITH A CERO ALFA',/,
     .' PARAMETER FOR THE EULER METHOD',/)
C
  775 FORMAT(
     .' THE "EXPLICIT_SOLUTION" OPTION IS NOT IMPLEMENTED YET FOR',/,
     .' HUGHES TIME MARCHING SCHEMES',/)
C
  776 FORMAT(
     .' THE "EXPLICIT_SOLUTION" OPTION IS NOT RECOMMENDABLE FOR  ',/,
     .' ADVECTION PROBLEMS WITH IPERT=3',/)
C
  800 FORMAT(
     .' WARNING: THE OPTION OF LINE-SEARCH FOR COUPLED PROBLEM HAS '/,
     .' NOT BEEN YET TESTED PROPERLY. PAY ATTENTIONS ON THE RESULTS.'/)
C
      END

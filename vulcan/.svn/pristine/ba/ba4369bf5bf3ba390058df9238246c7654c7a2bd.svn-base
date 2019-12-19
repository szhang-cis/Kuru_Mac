      SUBROUTINE IDEPROS9(PROPST,IPLAT,IA1)
C***********************************************************************
C     
C**** THIS ROUTINE ORDERS THE MICROSTRUCTURAL PROPERTIES OF MODEL
C     NUMBER 9 (IPCMO=9) OF RATE PHASE-CHANGE FORMULATIONS
C     
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C     
C**** THERMAL VARIABLES
C     
      INCLUDE 'auxl_omt.f'
      INCLUDE 'prob_omt.f'
C     
      DIMENSION PROPST(*)
C     
C***  chemical composition
C     
      IA2= IA1+2                ! 2=ipcfo,ipcmo

C     carbon
      CA     = PROPST(IA2+ 1)
      AKCAL  = PROPST(IA2+ 2)   ! carbon partition coeffic.
      AKCAS  = PROPST(IA2+ 3)   ! carbon partition coeffic.
      ICALSMX= INT(PROPST(IA2+ 4)) ! fla, fl & not consider Si sgregation
      ICASSMX= INT(PROPST(IA2+ 5)) ! fla, fl & not consider Si sgregation
C     silicon
      SI     = PROPST(IA2+ 6)
      AKSIL  = PROPST(IA2+ 7)   ! silicon partition coeffic.
      AKSIS  = PROPST(IA2+ 8)   ! silicon partition coeffic.
      ISILSMX= INT(PROPST(IA2+ 9)) ! fla, fl & not consider Si sgregation
      ISISSMX= INT(PROPST(IA2+10)) ! fla, fl & not consider Si sgregation
C     phosporus
      PH     = PROPST(IA2+11)
      AKPHL  = PROPST(IA2+12)   ! phosporus partition coeffic.
      AKPHS  = PROPST(IA2+13)   ! phosporus  partition coeffic.
      IPHLSMX= INT(PROPST(IA2+14)) ! fla, fl & not consider P sgregation
      IPHSSMX= INT(PROPST(IA2+15)) ! fla, fl & not consider P sgregation
C     copper 
      CU     = PROPST(IA2+16)
      AKCUL  = PROPST(IA2+17)   ! copper partition coeffic.
      AKCUS  = PROPST(IA2+18)   ! copper  partition coeffic.
      ICULSMX= INT(PROPST(IA2+19)) ! fla, fl & not consider Cu sgregation
      ICUSSMX= INT(PROPST(IA2+20)) ! fla, fl & not consider Cu sgregation
C     manganeso (manganese) MN (la Q es por la definicion)
      QN     = PROPST(IA2+21)
      AKQNL  = PROPST(IA2+22)   ! manganese partition coeffic.
      AKQNS  = PROPST(IA2+23)   ! manganese  partition coeffic.
      IQNLSMX= INT(PROPST(IA2+24)) ! fla, fl & not consider Mn sgregation
      IQNSSMX= INT(PROPST(IA2+25)) ! fla, fl & not consider Mn sgregation
C     magnesio  (magnesium) MG (la Q es por la definicion)
      QG     = PROPST(IA2+26)
      AKQGL  = PROPST(IA2+27)   ! magnesium partition coeffic.
      AKQGS  = PROPST(IA2+28)   ! magnesium  partition coeffic.
      IQGLSMX= INT(PROPST(IA2+29)) ! fla, fl & not consider Mg sgregation
      IQGSSMX= INT(PROPST(IA2+30)) ! fla, fl & not consider Mg sgregation
C     niobio (niobium) NB (la Q es por la definicion)
      QB     = PROPST(IA2+31)
      AKQBL  = PROPST(IA2+32)   ! niobium partition coeffic.
      AKQBS  = PROPST(IA2+33)   ! niobium  partition coeffic.
      IQBLSMX= INT(PROPST(IA2+34)) ! fla, fl & not consider Nb sgregation
      IQBSSMX= INT(PROPST(IA2+35)) ! fla, fl & not consider Nb sgregation
C     estano (tin) TIN
      SN     = PROPST(IA2+36)
      AKSNL  = PROPST(IA2+37)   ! tin partition coeffic.
      AKSNS  = PROPST(IA2+38)   ! tin  partition coeffic.
      ISNLSMX= INT(PROPST(IA2+39)) ! fla, fl & not consider Tin sgregation
      ISNSSMX= INT(PROPST(IA2+40)) ! fla, fl & not consider Tin sgregation
C     cromo (chromium) Cr
      CR    = PROPST(IA2+41)
      AKCRL  = PROPST(IA2+42)   ! chromium partition coeffic.
      AKCRS  = PROPST(IA2+43)   ! chromium  partition coeffic.
      ICRLSMX= INT(PROPST(IA2+44)) ! fla, fl & not consider Cr sgregation
      ICRSSMX= INT(PROPST(IA2+45)) ! fla, fl & not consider Cr sgregation
C     molibdeno (molybdenum) Mo
      QO     = PROPST(IA2+46)
      AKQOL  = PROPST(IA2+47)   ! molybdenum partition coeffic.
      AKQOS  = PROPST(IA2+48)   ! molybdenum  partition coeffic.
      IQOLSMX= INT(PROPST(IA2+49)) ! fla, fl & not consider Mo sgregation
      IQOSSMX= INT(PROPST(IA2+50)) ! fla, fl & not consider Mo sgregation
C     molibdeno (nickel) Ni
      QI     = PROPST(IA2+51)
      AKQIL  = PROPST(IA2+52)   ! nickel partition coeffic.
      AKQIS  = PROPST(IA2+53)   ! nickel  partition coeffic.
      IQILSMX= INT(PROPST(IA2+54)) ! fla, fl & not consider Ni sgregation
      IQISSMX= INT(PROPST(IA2+55)) ! fla, fl & not consider Ni sgregation
C     
C***  graphite nucleation
C     
      INUCMX= INT(PROPST(IA2+56))
      IF(INUCMX.EQ.1.OR.INUCMX.EQ.2.OR.
     .     INUCMX.EQ.3.OR.INUCMX.EQ.4) THEN ! Su, Boeri, Rappaz & Stef.
         ANUCA= PROPST(IA2+57)
         ANUCB= PROPST(IA2+58)
         ANUCC= PROPST(IA2+59)         
      ENDIF
C     
      INUCAX= INT(PROPST(IA2+60))
C     
      IGRGMX= INT(PROPST(IA2+61))
      IF(IGRGMX.EQ.1.OR.IGRGMX.EQ.2) THEN ! Dardati & Zener
         DIFCL= PROPST(IA2+62)
         RNODA= PROPST(IA2+63)
         RNODO= PROPST(IA2+64)
         AUSGR= PROPST(IA2+65)
         DENSA= PROPST(IA2+66)
         DENSG= PROPST(IA2+67)
      ENDIF
C     
      IGRAMX= INT(PROPST(IA2+68))
      IF(IGRAMX.EQ.1.OR.IGRAMX.EQ.2) THEN ! Rappaz & Stefanescu
         DIFCA= PROPST(IA2+69)
         AMLQ = PROPST(IA2+70)
         AKCA = PROPST(IA2+71)
         COGTH= PROPST(IA2+72)
      ENDIF
C     
      IDELRN = INT(PROPST(IA2+73)) ! model to growth zone's 2 radius
      ISDASGR= INT(PROPST(IA2+74)) ! model to growth zone's 2 radius
      IMICOUP= INT(PROPST(IA2+75)) ! models to couple micro and macro
C     
      IKMICX= INT(PROPST(IA2+76))
      IKAUX = 0
      IF(IKMICX.EQ.1) THEN
         IF(IMICOUP.EQ.1.OR.IMICOUP.EQ.2.OR.IMICOUP.EQ.3) THEN
            IKAUX= 3
            BASKS= PROPST(IA2+77)
            BASKM= PROPST(IA2+78)
            BASKL= PROPST(IA2+79)
         ENDIF
      ENDIF
C     
      IFPCDT= INT(PROPST(IA2+77+IKAUX))
      IAFLOJ= INT(PROPST(IA2+78+IKAUX))
C     
C**** transfer to VPLAT array -- explanation
C     
      VPLAT(IPLAT, 6)= CA
      VPLAT(IPLAT, 7)= AKCAL
      VPLAT(IPLAT, 8)= AKCAS
      VPLAT(IPLAT, 9)= FLOAT(ICALSMX)
      VPLAT(IPLAT,10)= FLOAT(ICASSMX)
      VPLAT(IPLAT,11)= SI
      VPLAT(IPLAT,12)= AKSIL
      VPLAT(IPLAT,13)= AKSIS
      VPLAT(IPLAT,14)= FLOAT(ISILSMX)
      VPLAT(IPLAT,15)= FLOAT(ISISSMX)
      VPLAT(IPLAT,16)= PH
      VPLAT(IPLAT,17)= AKPHL 
      VPLAT(IPLAT,18)= AKPHS
      VPLAT(IPLAT,19)= FLOAT(IPHLSMX)
      VPLAT(IPLAT,20)= FLOAT(IPHSSMX)
      VPLAT(IPLAT,21)= CU
      VPLAT(IPLAT,22)= AKCUL
      VPLAT(IPLAT,23)= AKCUS
      VPLAT(IPLAT,24)= FLOAT(ICULSMX)
      VPLAT(IPLAT,25)= FLOAT(ICUSSMX)
      VPLAT(IPLAT,26)= QN 
      VPLAT(IPLAT,27)= AKQNL
      VPLAT(IPLAT,28)= AKQNS
      VPLAT(IPLAT,29)= FLOAT(IQNLSMX)
      VPLAT(IPLAT,30)= FLOAT(IQNSSMX)
      VPLAT(IPLAT,31)= QG
      VPLAT(IPLAT,32)= AKQGL
      VPLAT(IPLAT,33)= AKQGS
      VPLAT(IPLAT,34)= FLOAT(IQGLSMX)
      VPLAT(IPLAT,35)= FLOAT(IQGSSMX)
      VPLAT(IPLAT,36)= QB 
      VPLAT(IPLAT,37)= AKQBL
      VPLAT(IPLAT,38)= AKQBS
      VPLAT(IPLAT,39)= FLOAT(IQBLSMX)
      VPLAT(IPLAT,40)= FLOAT(IQBSSMX)
      VPLAT(IPLAT,41)= SN
      VPLAT(IPLAT,42)= AKSNL
      VPLAT(IPLAT,43)= AKSNS
      VPLAT(IPLAT,44)= FLOAT(ISNLSMX)
      VPLAT(IPLAT,45)= FLOAT(ISNSSMX)
      VPLAT(IPLAT,46)= CR
      VPLAT(IPLAT,47)= AKCRL
      VPLAT(IPLAT,48)= AKCRS
      VPLAT(IPLAT,49)= FLOAT(ICRLSMX)
      VPLAT(IPLAT,50)= FLOAT(ICRSSMX)
      VPLAT(IPLAT,51)= QO
      VPLAT(IPLAT,52)= AKQOL
      VPLAT(IPLAT,53)= AKQOS
      VPLAT(IPLAT,54)= FLOAT(IQOLSMX)
      VPLAT(IPLAT,55)= FLOAT(IQOSSMX)
      VPLAT(IPLAT,56)= QI
      VPLAT(IPLAT,57)= AKQIL 
      VPLAT(IPLAT,58)= AKQIS
      VPLAT(IPLAT,59)= FLOAT(IQILSMX)
      VPLAT(IPLAT,60)= FLOAT(IQISSMX)
C     
      VPLAT(IPLAT,61)= FLOAT(INUCMX)
      VPLAT(IPLAT,62)= ANUCA
      VPLAT(IPLAT,63)= ANUCB
      VPLAT(IPLAT,64)= ANUCC
C     
      VPLAT(IPLAT,65)= FLOAT(INUCAX)
C     
      VPLAT(IPLAT,66)= FLOAT(IGRGMX)
      VPLAT(IPLAT,67)= DIFCL
      VPLAT(IPLAT,68)= RNODA
      VPLAT(IPLAT,69)= RNODO
      VPLAT(IPLAT,70)= AUSGR
      VPLAT(IPLAT,71)= DENSA
      VPLAT(IPLAT,72)= DENSG
C     
      VPLAT(IPLAT,73)= FLOAT(IGRAMX)
      VPLAT(IPLAT,74)= DIFCA
      VPLAT(IPLAT,75)= AMLQ  
      VPLAT(IPLAT,76)= AKCA  
      VPLAT(IPLAT,77)= COGTH 
C     
      VPLAT(IPLAT,78)= FLOAT(IDELRN)
      VPLAT(IPLAT,79)= FLOAT(ISDASGR)
      VPLAT(IPLAT,80)= FLOAT(IMICOUP)
C     
      VPLAT(IPLAT,81)= FLOAT(IKMICX)
      IF(IKMICX.EQ.1) THEN
         IF(IMICOUP.EQ.1.OR.IMICOUP.EQ.2.OR.IMICOUP.EQ.3) THEN
            VPLAT(IPLAT,82)= BASKS
            VPLAT(IPLAT,83)= BASKM
            VPLAT(IPLAT,84)= BASKL
         ENDIF
      ENDIF
C     
      VPLAT(IPLAT,82+IKAUX)= FLOAT(IFPCDT)
      VPLAT(IPLAT,83+IKAUX)= FLOAT(IAFLOJ)
C     
      IMODE= 78+IKAUX           ! imode=total number of prop. of model 9
C     
      IA1  = IA2+IMODE
C     
      RETURN
      END

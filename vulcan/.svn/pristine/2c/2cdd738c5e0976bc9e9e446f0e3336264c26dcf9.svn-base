      SUBROUTINE PHCHFVT(REFORT,RATIS,IFLAGX)
C***********************************************************************
C
C**** THIS ROUTINE APPLIES A PHASE-CHANGE FRONT VELOCITY CONTROL
C     ALGORITHM
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'                  ! variable strategy (Oct/97)
C
      DIMENSION REFORT(NTOTVT,2)
      DIMENSION RATIS(4)
C
      save mitertx,tolertx,itrucx,truchx
C
C**** PHASE-CHANGE FRONT VELOCITY CONTROL ALGORITHM (see convert.f)
C
      IF(IFLAGX.EQ.1) THEN
       if(iitert.eq.1) then
        call solcont(ratis(1))
c
c**** read itruc & truch for every step from a variable strategy file
c
        if(itrucv.eq.1) then
         itapet=lustrt                   ! always an external file
C
         IERR1T=1
         IF(ITAPET.NE.LUDATT.AND.ITAPET.NE.LUSTRT) GOTO 1000
         IF(ITAPET.EQ.LUSTRT) THEN       ! external .str file
          IF(IOSTRT.EQ.0)
     .                 OPEN(UNIT=LUSTRT,FILE=CI1T,STATUS='OLD',ERR=1000)
          IOSTRT=1
          IERR1T=0
C
 1000     IF(IERR1T.NE.0) THEN
           IF(IERR1T.EQ.1) WRITE(LUREST,901)
           CALL RUNENDT('ERROR IN OPENING FILE ')
          ENDIF
         ENDIF      ! itapet.eq.lustrt
C
         nprint=1
         call listent('convert',nprint,itapet)
C
         isteptx=int(paramt(1))
         mitertx=int(paramt(2))
         tolertx=paramt(3)
         itrucx=int(paramt(4))
         truchx=paramt(5)
C
c        extupx=paramt(6)      ! to be checked for advective problems !!
c        if(extupx.ne.0.0d0) extup=extupx
c        isupwx=int(paramt(7))
c        if(isupwx.ne.0) isupw=isupwx
C
         if(istept.ne.isteptx)
     .    call runendt('istept ne isetptx in strategy file')
        endif
       endif
c
c**** assign values for variable strategy
c
       if(itrucv.eq.1) then
        if(mitertx.gt.0.and.tolertx.gt.0.0d0) then
         mitert=mitertx
         tolert=tolertx
         itruc=itrucx
         if(truchx.gt.0.0d0) truch=truchx
        endif
       endif
c
c**** itruc < 0 ==> only tsoli varies (not nsol_i)
c
       if(iitert.eq.1) then
        if(itruc.gt.0) then
         isoli=nsol1-itruc
         nsol1=nsol1-isoli
         nsol2=nsol2-isoli
         nsol3=nsol3-isoli
        endif
       endif
c
c**** options (better as input): a) tsoli=1.0  b) tsoli=truch
c
       tsoli=1.0d0
       if(iitert.ge.nsol1) then
        if(iitert.lt.nsol2) tsoli=tsol1+(tsol2-tsol1)/
     .                      (nsol2-nsol1)*(iitert-nsol1)
        if(iitert.ge.nsol2.and.iitert.lt.nsol3) tsoli=tsol2+
     .                      (tsol3-tsol2)/(nsol3-nsol2)*(iitert-nsol2)
        if(iitert.ge.nsol3) tsoli=tsol3
        tsoli=tsoli*truch
       endif
C
       IF(IITERT.GE.NSOL1) THEN       ! tsoli < 1.0
        DO IDOFNT=1,NDOFNT
         DO IPOINT=1,NPOINT
          ITOTVT=(IPOINT-1)*NDOFCT+IDOFNT
          refort(itotvt,1)=refort(itotvt,1)*tsoli
         ENDDO
        ENDDO
       ENDIF
      ENDIF            ! iflagx.eq.1
C
C**** ALTERNATIVE PH-CH FRONT VELOCITY CONTROL ALGORITHM (see residut.f)
C
      IF(IFLAGX.EQ.2) THEN
       IF(TRUCHB.GT.0.0D0) THEN
C
        ITRUCB=0
        IF(IITERT.LE.1) THEN
         ITRUCB=1
        ELSE
         IF(IITERT.LT.NSOL1) ITRUCB=1
        ENDIF
C
        IF(ITRUCB.EQ.1) THEN
         tsoli=TRUCHB
         DO IDOFNT=1,NDOFNT
          DO IPOINT=1,NPOINT
           ITOTVT=(IPOINT-1)*NDOFCT+IDOFNT
           refort(itotvt,1)=refort(itotvt,1)*tsoli
          ENDDO
         ENDDO
        ENDIF
       ENDIF
      ENDIF            ! iflagx.eq.2
C
  901 FORMAT('ERROR IN OPENING VARIABLE STRATEGY INPUT FILE')
C
      RETURN
      END

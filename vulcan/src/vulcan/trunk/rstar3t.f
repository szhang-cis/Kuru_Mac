      SUBROUTINE RSTAR3T(DISPRT,DISTOT,ELPRET,ELVART,HEADST,HTLODT,
     .                   IFFIXT,REFORT,RLOADT,TLOADT,IFLAGT)
C***********************************************************************
C
C**** THIS ROUTINE REFORMS THE FOLLOWING TASKS:
C
C     NEW     RUN: WRITES THE RESULTS OF THE CURRENT TIME STEP
C     RESTART RUN: READS THE RESULTS OF THE SELECTED TIME STEP 
C                ( record pointer KPOST has been identified in 
C                  routine rschek )
C
C     - IFLAG = 1 - NEW     RUN: WRITE INTO RESTART FILE 
C     - IFLAG = 2 - RESTART RUN: READ  FROM RESTART FILE
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
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'
C
      DIMENSION HTLODT(NHLODT,NSUBFT,NFUNCT), IFFIXT(NTOTVT,2),
     .          RLOADT(NTOTVT)
      DIMENSION DISPRT(NTOTVT,3),             DISTOT(NTOTVT,3),
     .          HEADST(NPOINT,4), 
     .          REFORT(NTOTVT),               TLOADT(NTOTVT,2)
      DIMENSION ELPRET(NPREVT),               ELVART(NSTATT)
C
      CALL CPUTIMT(TIME1T)
C
      IF(IFLAGT.EQ.1)
     . OPEN(UNIT=LURSTT,FILE=CKT,STATUS='OLD',ACCESS='DIRECT',
     .      FORM='UNFORMATTED',RECL=LENRCT)
C
C**** WRITE OR READ ALL COMMON VARIABLES
C
      IF(IFLAGT.EQ.1)THEN
C
C**** WRITE COMMON BLOCK 'INTERVL_OM'
C
       WRITE(LURSTT,REC=NRECCT+1) 
C
C**** COMMON/CONSTNTA/
     .                DAMP1T,DAMP2T,DITERT,DTIMET,GRAVYT,GVECTT,
     .                STIFIT,TALFAT,TBETAT,TOLERT,XTIMET,WLUMPT,
C**** COMMON/CURRENTA/
     .                ITIMET,ISTEPT,IITERT,KDTIMT,KRESLT,KSTIFT,KUNLDT,
     .                NCHEKT,
C**** COMMON/CURRENTB/
     .                ARCLNT,FACTOT,FACPRT,PITERT,PVALUT,STICUT,
     .                STIINT,TFACTT,TTIMET,
C**** COMMON/INTERVTA/
     .                IMPOST,KALGOT,KARCLT,KDAMPT,KINTET,KOPTIT,KCONVT,
     .                LACCET,LAUTOT,LINEST,MITERT,NALGOT,NBACKT,NCDIST,
     .                NSTEPT,NOUTPT
C
C**** WRITE COMMON BLOCK 'INP_OUT_OM'
C
C**** COMMON/PLOTERTA/ & COMMON/PLOTERTB/
       WRITE(LURSTT,REC=NRECCT+2) NCOLDT,NCURVT,NPONTT,FORMAT
       WRITE(LURSTT,REC=NRECCT+3) MPLOTT
C
      ELSE
C
C**** READ COMMON BLOCK 'INTERVL_OM'
C
       READ(LURSTT,REC=NRECCT+1) 
C
C**** COMMON/CONSTNTA/
     .                DAMP1T,DAMP2T,DITERT,DTIMET,GRAVYT,GVECTT,
     .                STIFIT,TALFAT,TBETAT,TOLERT,XTIMET,WLUMPT,
C**** COMMON/CURRENTA/
     .                ITIMET,ISTEPT,IITERT,KDTIMT,KRESLT,KSTIFT,KUNLDT,
     .                NCHEKT,
C**** COMMON/CURRENTB/
     .                ARCLNT,FACTOT,FACPRT,PITERT,PVALUT,STICUT,
     .                STIINT,TFACTT,TTIMET,
C**** COMMON/INTERVTA/
     .                IMPOST,KALGOT,KARCLT,KDAMPT,KINTET,KOPTIT,KCONVT,
     .                LACCET,LAUTOT,LINEST,MITERT,NALGOT,NBACKT,NCDIST,
     .                NSTEPT,NOUTPT
C
C**** READ COMMON BLOCK 'INP_OUT_OM'
C
C**** COMMON/PLOTERTA/ & COMMON/PLOTERTA/
       READ(LURSTT,REC=NRECCT+2) NCOLDT,NCURVT,NPONTT,FORMAT
       READ(LURSTT,REC=NRECCT+3) MPLOTT
ctm
ctmC.....COMMON/I_P_CT/
ctm        READ(LURSTT,REC=NRECCT+4) IPCSTT,NTURNT,NITEMT,PNAMET,EXIMPT
ctm        READ(LURSTT,REC=NRECCT+5) TOLEDT,IPARAT
ctm
C
      ENDIF
C
C**** READ FROM PROCESS AREA AND WRITE TO CONVERGED AREA ELPRE & ELVAR
C
      IF(IFLAGT.EQ.1) THEN
       DO IELEMT=1,NELEMT
C
C**** READ ELPRE & ELVAR ( CURRENT )
C
        IF(NMEMO.EQ.1)
     .   CALL DATBAST(ELPRET,    2,    2)
        IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
     .   CALL DATBAST(ELVART,    3,    2)
C
C**** WRITE ELPRE & ELVAR
C
        IF(NMEMO.EQ.1)
     .   CALL DATRSTT(ELPRET,    2,    1)
        IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
     .   CALL DATRSTT(ELVART,    3,    1)
       ENDDO
      ENDIF
C
C**** READ OR WRITE ALL MATRICES
C
      CALL DATRSTT(DISTOT,   4,IFLAGT)
C
      CALL DATRSTT(DISPRT,   5,IFLAGT)
C
      CALL DATRSTT(HEADST,   6,IFLAGT)
C
      CALL DATRSTT(REFORT,   7,IFLAGT)
C
      CALL DATRSTT(TLOADT,   8,IFLAGT)
C
      CALL DATRSTT(HTLODT,   9,IFLAGT)
C
      CALL DATRSTT(IFFIXT,  10,IFLAGT)
C
      CALL DATRSTT(RLOADT,  11,IFLAGT)
C
C**** WRITE HEADINGS
C
      IF(IFLAGT.EQ.1) THEN
       WRITE(LURSTT,REC=     1) ITIMET,ISTEPT,TTIMET,NRECCT,KTSTET
       WRITE(LURSTT,REC=NRECCT) ITIMET,ISTEPT,TTIMET,NRECCT,KTSTET
      ENDIF
C
      IF(IFLAGT.EQ.1) CLOSE(LURSTT)
C
      CALL CPUTIMT(TIME2T)
      CPURST=CPURST+(TIME2T-TIME1T)
      RETURN
C
      END

      PROGRAM FINDDQE
C
      IMPLICIT NONE
C
      INTEGER I,J,K,NSAM,N2,NSTACK,ID,IERR,R2,RLIM2
      INTEGER SLEN2,MODE,I1,I2,L,M,IS,MRAD,EDGE,EDG2
      INTEGER JC,IN,NN,IW,XMIN,YMIN,N1,EDG3,WSAM
      INTEGER IXMIN,IXMAX,IYMIN,IYMAX,NFIT,EXCL,NI
      INTEGER NEST,IPAD,NKNOT,IPATCH,U1,IPF,IT
      INTEGER,ALLOCATABLE :: IP(:),IWRK(:)
      PARAMETER (MRAD=8,EDGE=10,EDG2=200,EDG3=200,NFIT=5)
      PARAMETER (EXCL=5,IPAD=500,NKNOT=100,IPATCH=200)
      PARAMETER (NEST=NKNOT,IPF=4)
      REAL PSIZE,THRESH,R,A,PI,N0,NH,SCAL,EQ,DQ0,IC,SINC
      REAL S,SI,SX,X(2*NFIT),E(2*NFIT),SF(4),DD,MAXS
      REAL DQ,DQI,NPS0,PFLAT,DOSE,SINCF,LOGIST
      PARAMETER (PI=3.1415926535897, PFLAT=0.01)
      REAL,ALLOCATABLE :: DATA(:),TH(:),EG(:),MA(:),DAT2(:)
      REAL,ALLOCATABLE :: BUF(:),N(:),NP(:),DLINE(:),WRK(:)
      REAL,ALLOCATABLE :: DIFF(:),MTF(:),WDAT(:),WTH(:)
      REAL,ALLOCATABLE :: SFIT(:),RS(:),IPDAT(:)
      REAL,ALLOCATABLE :: W(:),KNOTS(:),CSPLINE(:)
      DOUBLE PRECISION DMIN,DMAX,DMEAN,D1,D2,S1,S2
      DOUBLE PRECISION D12,D22,S12,S22
      DOUBLE PRECISION DRMS,NW,DMIN2,DMAX2,DMEAN2
      DOUBLE PRECISION,ALLOCATABLE :: PN(:)
      COMPLEX,ALLOCATABLE :: NQ(:),BUQ(:),WDAQ(:),IPDAQ(:)
      CHARACTER FNAME*200,FOUT*200,CFORM,TITLE*1600
      CHARACTER TOUT*200,LOUT*200
      CHARACTER LINE*200,VX*15
      LOGICAL LCOUNT
      DATA  VX/'1.05 - 21.10.13'/
C     15 chars 'X.XX - XX.XX.XX' <--<--<--<--<--<--<--<--
C*********************************************************************
C
      WRITE(*,7001) VX
7001  FORMAT(/' FindDQE - Determine DQE of a detector',
     *       /' V',A15,/,
     +       /'  Copyright 2013 Howard Hughes Medical Institute.',
     +       /'  All rights reserved.',
     +       /'  Use is subject to Janelia Farm Research Campus',
     +        ' Software Copyright 1.1',
     +       /'  license terms ( http://license.janelia.org/license',
     +        '/jfrc_copyright_1_1.html )'//)
C
      CFORM='M'
C
      WRITE(*,*)
     +  'You can input either the counts/electron'
      WRITE(*,*)
     +  '(gain conversion factor) or the total electron dose/pixel.'
      WRITE(*,*)
     +  'Enter 1 to input counts/electron, 2 for total dose/pixel'
      READ(*,*) U1
      WRITE(*,*) U1
      WRITE(*,*)  
C
      IF (U1.EQ.1) THEN
        WRITE(*,*)
     +  'Enter counts/electron (gain conversion factor)'
        READ(*,*) EQ
        WRITE(*,*) EQ
        WRITE(*,*)
      ELSEIF (U1.EQ.2) THEN
        WRITE(*,*)   
     +  'Enter total electron dose/pixel'
        READ(*,*) DOSE
        WRITE(*,*) DOSE
        WRITE(*,*)
      ELSE
        WRITE(*,*) 'Error: must input 1 or 2.'
        GOTO 9999  
      ENDIF
C
      WRITE(*,*)
     +  'Enter counts multiplier (usually 1).'
      READ(*,*) IC
      WRITE(*,*) IC
      WRITE(*,*)
C
      WRITE(*,*)    
     +  'Enter number of images'
      WRITE(*,*)  
     +  '(usually 1; 2 if also giving flat field image)'
      READ(*,*) NI
      WRITE(*,*) NI
      WRITE(*,*)
      IF ((NI.NE.1).AND.(NI.NE.2)) THEN
        WRITE(*,*) 'Error: must input 1 or 2.'
        GOTO 9999
      ENDIF
C
      WRITE(*,*)' Input image with pointer silhouette'
      READ(*,7006) FNAME
7006  FORMAT(A200)
      WRITE(*,17006) FNAME(1:SLEN2(FNAME))
17006 FORMAT(3X,A)
      CALL IOPEN(FNAME,10,CFORM,MODE,N1,N2,NSTACK,
     +           'O',PSIZE,TITLE)
C
      IF (NSTACK.GT.1) THEN
        WRITE(*,*) ' ERROR: Image must be 2D'
        GOTO 9999
      ENDIF
C
      IF (NI.NE.1) THEN
        WRITE(*,*)' 2nd (flat field) input image'
        READ(*,7006) FNAME
        WRITE(*,17006) FNAME(1:SLEN2(FNAME))
        WRITE(*,*)
      ENDIF
C
      NSAM=INT(MIN(N1,N2)/2)*2
      XMIN=INT((N1-NSAM)/2)
      YMIN=INT((N2-NSAM)/2)
      IF (N1.NE.N2) THEN
        WRITE(*,*) ' Cropping image to square size'
        WRITE(*,*) ' Xmin, Ymin, NSAM = ',XMIN,YMIN,NSAM
      ENDIF
C
      JC=NSAM/2+1
      IW=NSAM/200
      SCAL=NSAM*NSAM
C
      WRITE(*,*)
      WRITE(*,*)' Output residual image'
      READ(*,7006) FOUT
      WRITE(*,17006) FOUT(1:SLEN2(FOUT))
C     
      MODE=2
C
      ALLOCATE(DATA(NSAM*N2),TH(NSAM*N2),MA(NSAM*N2),
     +         EG(NSAM*N2),BUF(NSAM*N2*IPF*IPF),N(NSAM*N2),
     +         IP(JC),NQ(NSAM),PN(JC),NP(NSAM*NSAM/2+NSAM),
     +         BUQ(NSAM*IPF),DLINE(N1),MTF(JC),RS(JC+IPAD),
     +         W(JC+IPAD),KNOTS(NEST),CSPLINE(NEST),
     +         SFIT(JC+IPAD),WRK(JC*4+NEST*16),IWRK(NEST),
     +         STAT=IERR)
      IF (IERR.NE.0) THEN
        WRITE(*,*) ' ERROR: Memory allocation failed for arrays'
        GOTO 9999
      ENDIF
C
C  Read in image, determine min, max, mean
      WRITE(*,*)
      WRITE(*,*)' Reading in image'
C
      DMIN=1.0E30
      DMAX=-1.0E30
      DMEAN=0.0
      DO 10 I=1,NSAM
        CALL IREAD(10,DLINE,I+YMIN)
        DO 10 J=1,NSAM
          ID=J+NSAM*(I-1)
          DATA(ID)=DLINE(J+XMIN)*IC
          DMEAN=DMEAN+DATA(ID)
          IF (DATA(ID).LT.DMIN) DMIN=DATA(ID)
          IF (DATA(ID).GT.DMAX) DMAX=DATA(ID)
10    CONTINUE
      DMEAN=DMEAN/NSAM/NSAM
      CALL ICLOSE(10)
C
      WRITE(*,*)
      WRITE(*,*) 'DMIN,DMAX,DMEAN=',REAL(DMIN)/IC,
     +           REAL(DMAX)/IC,REAL(DMEAN)/IC
      FLUSH(6)
C
      IF (NI.NE.1) THEN
        CALL IOPEN(FNAME,10,CFORM,MODE,I,J,NSTACK,
     +             'O',PSIZE,TITLE)
C
        IF (NSTACK.GT.1) THEN
           WRITE(*,*) ' ERROR: Image must be 2D'
           GOTO 9999
        ENDIF
        IF ((I.NE.N1).OR.(J.NE.N2)) THEN
          WRITE(*,*) ' ERROR: Image must have same size as 1st image'
          GOTO 9999
        ENDIF
C
        ALLOCATE(DAT2(NSAM*N2),STAT=IERR)
        IF (IERR.NE.0) THEN
          WRITE(*,*) ' ERROR: Memory allocation failed for arrays'
          GOTO 9999
        ENDIF
C
        S=0.0
        DMIN2=1.0E30
        DMAX2=-1.0E30
        DMEAN2=0.0
        DO 11 I=1,NSAM
          CALL IREAD(10,DLINE,I+YMIN)
          DO 11 J=1,NSAM
            ID=J+NSAM*(I-1)
            DAT2(ID)=DLINE(J+XMIN)*IC
            DMEAN2=DMEAN2+DAT2(ID)
            IF (DAT2(ID).LT.DMIN2) DMIN2=DAT2(ID)
            IF (DAT2(ID).GT.DMAX2) DMAX2=DAT2(ID)
            S=S+(DAT2(ID)-DATA(ID))**2
11      CONTINUE
        DMEAN2=DMEAN2/NSAM/NSAM
        CALL ICLOSE(10)
C
        WRITE(*,*) 'DMIN,DMAX,DMEAN=',REAL(DMIN2)/IC,
     +             REAL(DMAX2)/IC,REAL(DMEAN2)/IC
        FLUSH(6)
C
        IF (S.EQ.0.0) THEN
          WRITE(*,*)
          WRITE(*,*) ' ERROR: 1st and 2nd images are identical'
          GOTO 9999
        ENDIF
      ENDIF
C
C  Determine image threshold
      WRITE(*,*)
      WRITE(*,*)' Determining image threshold'
C
      THRESH=DMEAN/2.0
      DO 40 J=1,5
        D1=0.0
        D2=0.0
        S1=0.0
        S2=0.0
        I1=0
        I2=0
        IF (NI.NE.1) THEN
          D12=0.0
          D22=0.0
          S12=0.0
          S22=0.0
        ENDIF
        DO 30 I=1,NSAM*NSAM
          IF (DATA(I).GE.THRESH) THEN
            I2=I2+1
            D2=D2+DATA(I)
            S2=S2+DATA(I)**2
            IF (NI.NE.1) THEN
              D22=D22+DAT2(I)
              S22=S22+DAT2(I)**2
            ENDIF
          ELSE
            I1=I1+1
            D1=D1+DATA(I)
            S1=S1+DATA(I)**2
            IF (NI.NE.1) THEN
              D12=D12+DAT2(I)
              S12=S12+DAT2(I)**2
            ENDIF
          ENDIF
30      CONTINUE
        IF (I1.NE.0) D1=D1/I1
        IF (I1.NE.0) S1=S1/I1
        IF (I2.NE.0) D2=D2/I2
        IF (I2.NE.0) S2=S2/I2
        THRESH=(D2-D1)/2.0+D1
        IF (NI.NE.1) THEN
          IF (I1.NE.0) D12=D12/I1
          IF (I1.NE.0) S12=S12/I1
          IF (I2.NE.0) D22=D22/I2
          IF (I2.NE.0) S22=S22/I2
        ENDIF
40    CONTINUE
      THRESH=(D2-D1)/2.0+D1
      S1=SQRT(S1-D1**2)
      S2=SQRT(S2-D2**2)
      IF (NI.NE.1) THEN
        S12=SQRT(S12-D12**2)
        S22=SQRT(S22-D22**2)      
      ENDIF
      WRITE(*,*)
      WRITE(*,*) 'Threshold = ',THRESH/IC
      WRITE(*,*) 'N, DMEAN, DRMS inside pointer  = ',
     +            I1,REAL(D1)/IC,REAL(S1)/IC
      WRITE(*,*) 'N, DMEAN, DRMS outside pointer = ',
     +            I2,REAL(D2)/IC,REAL(S2)/IC
      IF (U1.EQ.1) THEN
        WRITE(*,*) 'Calculated dose (e/A^2)        =  ',
     +              REAL(D2)/EQ
      ELSE
        WRITE(*,*) 'Gain conversion factor (counts/e) =  ',
     +              REAL(D2)/DOSE
                    EQ=REAL(D2)/DOSE
      ENDIF
      WRITE(*,*)
      IF (NI.NE.1) THEN
        WRITE(*,*) '2nd image:'
        WRITE(*,*) 'N, DMEAN, DRMS outside pointer = ',
     +              I2,REAL(D22)/IC,REAL(S22)/IC
        IF (U1.EQ.1) THEN
          WRITE(*,*) 'Calculated dose (e/A^2)        =  ',
     +                REAL(D22)/EQ
        ELSE
          WRITE(*,*) 'Gain conversion factor (counts/e) = ',
     +                REAL(D22)/DOSE
                      EQ=REAL(D22)/DOSE
        ENDIF
        WRITE(*,*)
        IF (ABS(D2-D22)/D2.GT.0.05) THEN
          WRITE(*,*) ' ERROR: Dose in the two images differs'
          GOTO 9999
        ENDIF
      ENDIF
      FLUSH(6)
C
C  Threshold image
      WRITE(*,*)' Thresholding image'
C
      DO 50 I=1,NSAM
        DO 50 J=1,NSAM
          ID=J+NSAM*(I-1)
          TH(ID)=D2
          IF (DATA(ID).LT.THRESH) TH(ID)=D1
50    CONTINUE
C
C  Eliminate single pixel outliers
      WRITE(*,*)' Removing single pixel outliers'
C
      IN=0
      DO 51 I=1,NSAM
        DO 51 J=1,NSAM
          ID=J+NSAM*(I-1)
          K=0
          S=0.0
          DO 52 L=-5,5
            DO 52 M=-5,5
            IF ((I+L.GE.1).AND.(I+L.LE.NSAM).AND.
     +           (J+M.GE.1).AND.(J+M.LE.NSAM).AND.
     +           (L.NE.M)) THEN
              IS=J+M+NSAM*(I+L-1)
              S=S+TH(IS)
              K=K+1
            ENDIF
52        CONTINUE
          S=(S/K-D1)/(D2-D1)
          IF (ABS(S-(TH(ID)-D1)/(D2-D1)).GT.0.7) THEN
            IF (S.GT.0.5) THEN
               TH(ID)=D2
            ELSE
               TH(ID)=D1
            ENDIF
            IN=IN+1
          ENDIF
          IF ((I.LE.EXCL).OR.(I.GT.NSAM-EXCL)) TH(ID)=D2
          IF ((J.LE.EXCL).OR.(J.GT.NSAM-EXCL)) TH(ID)=D2
51    CONTINUE
      WRITE(*,*)
      WRITE(*,*) 'Number of pixel outliers = ',IN
      IF (IN.GT.0.02*NSAM*NSAM) WRITE(*,*) 
     +  'Try a higher dose to reduce outliers to zero'
      WRITE(*,*)
C
      IXMIN=NSAM
      IYMIN=NSAM
      IXMAX=1
      IYMAX=1
      DO 240 I=1,NSAM
        DO 240 J=1,NSAM
          ID=J+NSAM*(I-1)
          IF (TH(ID).LT.(D1+D2)/2.0) THEN
            IF (IXMIN.GT.J) IXMIN=J
            IF (IXMAX.LT.J) IXMAX=J
            IF (IYMIN.GT.I) IYMIN=I
            IF (IYMAX.LT.I) IYMAX=I
          ENDIF
240   CONTINUE
C
      I=0
      IF ((IXMIN.LT.2*EDG2+1+EXCL).OR.
     +    (IXMAX.GT.NSAM-2*EDG2-EXCL)) I=I+1
      IF ((IYMIN.LT.2*EDG2+1+EXCL).OR.
     +    (IYMAX.GT.NSAM-2*EDG2-EXCL)) I=I+1
      IF (I.EQ.2) THEN
        WRITE(*,1300) '  Pointer distance from edge:',
     +    IXMIN,IXMAX,IYMIN,IYMAX
1300    FORMAT(A,4I10)
        WRITE(*,*)
     +  ' ERROR: Pointer less than',2*EDG2+EXCL,' from edge'
        GOTO 9999
      ENDIF
C
      WSAM=MIN(IXMAX-IXMIN+1,IYMAX-IYMIN+1)+4*EDG2
      IF (MOD(WSAM,2).NE.0.0) WSAM=WSAM+1
      CALL PCHECK(WSAM)
      ALLOCATE(WDAT(WSAM*WSAM),WDAQ(WSAM),WTH(WSAM*WSAM),
     +         DIFF(WSAM*WSAM),IPDAT(WSAM*WSAM*IPF*IPF),
     +         IPDAQ(WSAM*IPF),STAT=IERR)
      IF (IERR.NE.0) THEN
        WRITE(*,*) ' ERROR: Memory allocation failed for arrays'
        GOTO 9999
      ENDIF
C
C  Generate edge mask
      WRITE(*,*)' Generating smooth-edged mask for NPS calculation'
      FLUSH(6)
C
      DO 53 I=1,NSAM*NSAM
        EG(I)=0.0
        MA(I)=0.0
53    CONTINUE
C
      DD=(D2-D1)/2.0
      DO 61 I=1,NSAM
        DO 61 J=1,NSAM
          ID=J+NSAM*(I-1)
          DO 71 L=-1,1
            DO 71 M=-1,1
            IF ((I+L.GE.1).AND.(I+L.LE.NSAM).AND.
     +           (J+M.GE.1).AND.(J+M.LE.NSAM)) THEN
              IS=J+M+NSAM*(I+L-1)
              IF (TH(ID)-TH(IS).GT.DD) EG(ID)=1.0
            ENDIF
71        CONTINUE
61    CONTINUE
C
      RLIM2=(MRAD+EDG3)**2
      DO 81 I=1,NSAM
        DO 81 J=1,NSAM
          ID=J+NSAM*(I-1)
          IF (EG(ID).EQ.1.0) THEN
            DO 91 L=-MRAD-EDG3,MRAD+EDG3
              DO 91 M=-MRAD-EDG3,MRAD+EDG3
                R2=L**2+M**2
                IF (R2.LE.RLIM2) THEN
                IF ((I+L.GE.1).AND.(I+L.LE.NSAM).AND.
     +             (J+M.GE.1).AND.(J+M.LE.NSAM)) THEN
                  IS=J+M+NSAM*(I+L-1)
                  A=1.0
                  R=SQRT(REAL(R2))
                  IF (R.GT.MRAD)
     +              A=(1.0+COS((R-MRAD)/EDG3*PI))/2.0
                  IF (MA(IS).LT.A) MA(IS)=A
                ENDIF
              ENDIF
91          CONTINUE
          ENDIF
81    CONTINUE
C
      DO 100 I=1,NSAM*NSAM
        MA(I)=1.0-MA(I)
        IF (TH(I).LT.D2-DD) MA(I)=0.0
100   CONTINUE
C
C  Taper edges of image
      WRITE(*,*)' Tapering edges of image'
      FLUSH(6)
C
      DO 110 I=1,NSAM
        DO 110 J=1,NSAM
          A=1.0
          IF (I.LE.EDG2) A=A*(1.0-COS(PI*(I-1)/EDG2))/2.0
          IF (I.GT.NSAM-EDG2) A=A*(1.0-COS(PI*(NSAM-I)/EDG2))/2.0
          IF (J.LE.EDG2) A=A*(1.0-COS(PI*(J-1)/EDG2))/2.0
          IF (J.GT.NSAM-EDG2) A=A*(1.0-COS(PI*(NSAM-J)/EDG2))/2.0
          ID=J+NSAM*(I-1)
          MA(ID)=MA(ID)*A
110   CONTINUE
C
C  Apply background mask
      WRITE(*,*)' Applying mask to background'
      FLUSH(6)
C
      NW=0.0
      DMEAN=0.0
      DRMS=0.0
      IN=0
      DO 102 I=1,NSAM*NSAM
        IF (MA(I).GT.0.999) THEN
          IF (NI.EQ.1) THEN
            DMEAN=DMEAN+DATA(I)
            DRMS=DRMS+DATA(I)**2
          ELSE
            DMEAN=DMEAN+DATA(I)-DAT2(I)
            DRMS=DRMS+(DATA(I)-DAT2(I))**2
          ENDIF
          IN=IN+1
        ENDIF
        IF (NI.EQ.1) THEN
          N(I)=MA(I)*DATA(I)+(1.0-MA(I))*D2
        ELSE
          N(I)=MA(I)*(DATA(I)-DAT2(I))
        ENDIF
        NW=NW+MA(I)**2
102   CONTINUE
      DMEAN=DMEAN/IN
      DRMS=SQRT(DRMS/IN-DMEAN**2)
      IF (NI.NE.1) DRMS=DRMS/SQRT(2.0)
C  Calculate noise power spectrum
      WRITE(*,*)' Calculating noise power spectrum'
      FLUSH(6)
C
      DO 180 I=1,NSAM*NSAM
        N(I)=N(I)/SCAL/DMAX
180   CONTINUE
C
      CALL RLFT3(N,NQ,NSAM,NSAM,1,1)
C
C  THIS IS PROBABLY WHERE TO APPLY ANY FILTER
      I1=(NSAM/3)**2
      I2=(2*NSAM/3)**2
      DO 140 M=1,NSAM
        I=M-1
        IF (I.GE.JC) I=I-NSAM
        DO 150 L=1,NSAM/2
          J=L-1
          ID=(L+NSAM/2*(M-1))*2
          IS=L+(NSAM/2+1)*(M-1)
          NP(IS)=N(ID-1)**2+N(ID)**2
150     CONTINUE
        IS=JC*M
        NP(IS)=CABS(NQ(M))**2
140   CONTINUE
C
C  Calculate rotational average of power spectrum
      WRITE(*,*)' Calculating rotational average'
      FLUSH(6)
C
      NN=0
      N1=0
      N0=0.0
      NH=0.0
      I2=(NSAM/2)**2
      DO 120 K=1,JC
        PN(K)=0.0
        IP(K)=0
120   CONTINUE
      DO 130 L=1,JC
        J=L-1
        DO 130 M=1,NSAM
          I=M-1
          IF (I.GE.JC) I=I-NSAM
          R2=I**2+J**2
          IF (R2.LE.I2) THEN
            ID=L+JC*(M-1)
            K=NINT(SQRT(REAL(R2)))+1
            PN(K)=PN(K)+DBLE(NP(ID))
            IP(K)=IP(K)+1
          ENDIF
130   CONTINUE
      DO 121 K=1,JC
        IF (IP(K).NE.0) THEN
          PN(K)=PN(K)/IP(K)
        ELSE
          PN(K)=1.0
        ENDIF
        IF (NN.LT.20) THEN
          IF (IP(K).GT.50) THEN
            NN=NN+1
            N0=N0+PN(K)
          ENDIF
        ENDIF
        IF (N1.LT.20) THEN
          IF ((IP(K).GT.50).AND.(K.GT.JC-30)) THEN
            N1=N1+1
            NH=NH+PN(K)
          ENDIF
        ENDIF
121   CONTINUE
C     Power spectrum with flat start
C      PN(1:NINT(JC*PFLAT))=PN(NINT(JC*PFLAT)+1:2*NINT(JC*PFLAT))
      N0=N0/20
      NH=NH/20
      WRITE(*,*)' Fitting noise spectrum with splines'
      CALL FIT_SPLINES(JC,RS,N0,PN,W,KNOTS,NKNOT,IPAD,
     +           IPATCH,CSPLINE,SFIT,NPS0,NEST,IWRK,WRK)
      LCOUNT=.FALSE.
      IF (NH/N0.GT.0.5) THEN
        LCOUNT=.TRUE.
        WRITE(*,*) 'Counting mode detected. MTF will be'
        WRITE(*,*) 'calculated to reflect noise depression'
        WRITE(*,*) 'at low resolution due to lost counts.'
        CALL FIT_COUNTING_NPS(JC,N0,PN,SF,IT)
      ENDIF
C
C  Find MTF by fitting edge
      WRITE(*,1100) NFIT
1100  FORMAT('  Determining MTF from edge blurring using',I2,
     +       '-Gaussian model')
      FLUSH(6)
C
      XMIN=IXMIN-2*EDG2
      IF (XMIN.LE.EXCL) XMIN=EXCL+1
      YMIN=IYMIN-2*EDG2
      IF (YMIN.LE.EXCL) YMIN=EXCL+1
      DO 200 I=1,WSAM
        DO 200 J=1,WSAM
          ID=J+XMIN-1+NSAM*(I+YMIN-2)
          IS=J+WSAM*(I-1)
          WDAT(IS)=DATA(ID)
          WTH(IS)=TH(ID)
200   CONTINUE
      CALL FIT_EDGE(WSAM,EDG2,THRESH,D1,D2,LCOUNT,SF,WDAT,
     +              WDAQ,WTH,DIFF,NFIT,X,E,BUF,BUQ,
     +              IPDAT,IPDAQ,IPF,IT)
C
C  Generate MTF curve
      S=0.0
      DO 230 I=1,NFIT
        S=S+ABS(X(2*I))
230   CONTINUE
C      X(2*NFIT)=1.0-S
      DO 210 I=1,JC
        SX=0.0
        DO 220 J=1,NFIT
          SX=SX+ABS(X(2*J))*EXP(-X(2*J-1)**2*(I-1)**2/NSAM/NSAM)
220     CONTINUE
        IF (I.NE.1) THEN
          MTF(I)=SX*SINC(PI*(I-1)/NSAM)/S
        ELSE
          MTF(1)=1.0
          S=SX
        ENDIF
210   CONTINUE
      A=-1.0E30
      IF (LCOUNT) THEN
        DO 123 K=1,JC
          R=2.0*REAL(K-1)/NSAM
          IF (IT.EQ.1) THEN
            S=SINCF(SF,R)
          ELSE
            S=LOGIST(SF,R)
          ENDIF
          IF (S.GT.A) A=S
123     CONTINUE
        DO 122 K=1,JC
          R=2.0*REAL(K-1)/NSAM
          IF (IT.EQ.1) THEN
            S=SINCF(SF,R)/A
          ELSE
            S=LOGIST(SF,R)/A
          ENDIF
          MTF(K)=MTF(K)*SQRT(S)
          IF (MTF(K).GT.1.0) MTF(K)=1.0
122     CONTINUE
      ENDIF
C
      CALL IOPEN(FOUT,20,CFORM,MODE,WSAM,WSAM,1,
     +           'N',PSIZE,TITLE)
      DO 20 I=1,WSAM
        ID=1+WSAM*(I-1)
        CALL IWRITE(20,DIFF(ID),I)
20    CONTINUE
C
      CALL ICLOSE(20)
C
C  Scale noise power spectrum according to curve fit
C  Scale factor is curve fir value at zero for linear detectors
C  Scale factor is curve fit value at Nyquist for counting detectors
      N0=N0*NPS0
      NH=NH*NPS0
      DRMS=SQRT(N0)*DMAX*SCAL/SQRT(NW)
      IF (NI.NE.1) DRMS=DRMS/SQRT(2.0)
      DQ0=D2*EQ/DRMS**2*MTF(1)**2
      WRITE(*,*) 'RMS of noise adjusted for the MTF = ',
     +          REAL(DRMS)/IC
      WRITE(*,*)
C
      WRITE(*,*) 'DQE(0) = ',DQ0
      FLUSH(6)
C
      WRITE(*,*)
      WRITE(*,*) 'Frequency       MTF       NPS   NPS-Fit       DQE'
C
      MAXS=-1.0E30
      DO 250 K=1,JC
        S=SFIT(K)/NPS0
        IF (S.GT.MAXS) MAXS=S
250   CONTINUE
      DO 190 K=1,JC
        R=2.0*REAL(K-1)/NSAM
          S=SFIT(K)/NPS0
          DQ=DQ0*MTF(K)**2/S
        IF (PN(K).GT.2.0*N0) PN(K)=2.0*N0
        WRITE(*,1000) R,MTF(K),REAL(PN(K))/N0,S,DQ
1000    FORMAT(5F10.6)
190   CONTINUE
      FLUSH(6)
C
9999  CONTINUE
      END
C*********************************************************************
      SUBROUTINE BINNING(NSAM,DATA,OUT,IB)
C
      IMPLICIT NONE
C
      INTEGER NSAM,IB,I,J,ID,IS,NB,NJ,L,M,NL,NI,IB2
      REAL DATA(*),OUT(*),R,SCAL
C
      IB2=IB/2
      SCAL=1.0/IB/IB
      NB=NSAM/IB
      DO 10 J=1,NB-1
        NJ=NB*J+1
        DO 10 I=1,NB-1
          IS=I+NJ
          NI=IB*I
          R=0.0
          DO 20 L=-IB2+1,IB-IB2
            NL=NI+NSAM*(IB*J+L)+1
            DO 20 M=-IB2+1,IB-IB2
              ID=M+NL
              R=R+DATA(ID)
20        CONTINUE
          OUT(IS)=R*SCAL
10    CONTINUE
C
      RETURN
      END
C*********************************************************************
      SUBROUTINE BILIN_INTERPO(NSAM,DATA,OUT,IPF)
C
      IMPLICIT NONE
C
      INTEGER NSAM,IPF,I,J,ID,NP,NJ
      REAL DATA(*),OUT(*),X,Y,BILIN
C
      NP=IPF*NSAM
      DO 10 J=0,NP-1
        NJ=NP*J+1
        DO 10 I=0,NP-1
          ID=I+NJ
          X=REAL(I-0.5)/IPF
          Y=REAL(J-0.5)/IPF
          OUT(ID)=BILIN(NSAM,DATA,X,Y)
10    CONTINUE
C
      RETURN
      END
C*********************************************************************
      REAL FUNCTION BILIN(NSAM,DATA,X,Y)
C
      IMPLICIT NONE
C
      INTEGER NSAM,LS,LT,MS,MT,ID1,ID2,ID3,ID4,N
      REAL X,Y,X1,Y1
      REAL DATA(*),SAMP
C
c      N=NSAM-1
      LS=INT(X)
c      LS=INT(X)-1
c      IF (LS.LT.0) LS=0
      X1=X-LS
      LT=LS+1
      IF (LT.GT.N) LT=N
      MS=INT(Y)
c      MS=INT(Y)-1
c      IF (MS.LT.0) MS=0
      Y1=Y-MS
      MT=MS+1
      IF (MT.GT.N) MT=N
      ID1=LS+1+NSAM*MS
      ID2=LT+1+NSAM*MS
      ID3=LS+1+NSAM*MT
      ID4=LT+1+NSAM*MT
C
      SAMP=DATA(ID1)*(1.0-X1)*(1.0-Y1)
      SAMP=SAMP+DATA(ID2)*X1*(1.0-Y1)
      SAMP=SAMP+DATA(ID3)*(1.0-X1)*Y1
      BILIN=SAMP+DATA(ID4)*X1*Y1
C
      RETURN
      END
C*********************************************************************
      REAL FUNCTION BILIN_TEST(NSAM,DATA,X,Y)
C
      IMPLICIT NONE
C
      INTEGER NSAM,LS,LT,MS,MT
      INTEGER L,M,LL,MM,ID
      REAL X,Y,WGT
      REAL DATA(*),SAMP
C
      SAMP=0.0
      LS=INT(X)
C      IF (X.LT.0.0) LS=LS-1
      LT=LS+1
      MS=INT(Y)
C      IF (Y.LT.0.0) MS=MS-1
      MT=MS+1
C
      DO 40 L=LS,LT
        LL=L+1
C        LL=L-INT(REAL(L)/NSAM)*NSAM
C        IF (LL.LT.0) LL=LL+NSAM
        DO 40 M=MS,MT
C          MM=M-INT(REAL(M)/NSAM)*NSAM
C          IF (MM.LT.0) MM=MM+NSAM
C	  BILINEAR INTERPOLATION.....
      	  WGT=(1.0-ABS(X-L))*(1.0-ABS(Y-M))
C      	    ID=LL+1+NSAM*MM
      	    ID=LL+NSAM*M
            SAMP=SAMP+DATA(ID)*WGT
40    CONTINUE
      BILIN_TEST=SAMP
C
      RETURN
      END
C*********************************************************************
      SUBROUTINE FIT_COUNTING_NPS(JC,N0,PN,SF,IT)
C
      IMPLICIT NONE
C
      INTEGER JC,I,J,IMIN,K,IW,L,IT
      PARAMETER (IW=50)
      REAL SF(*),ESCALE,R,E(4),N0,D,SUM,MINS,S(4)
      REAL EVAL_SINCF,S2(4),EVAL_LOGIST
      DOUBLE PRECISION PN(*)
      DATA ESCALE/100.0/
C
C Sinc formula: NPS=SF1-SF2*SINC(SF3*R)^2
C Li et al., JSB 2013
C
      E(1)=0.001
      E(2)=0.001
      E(3)=0.001
      E(4)=0.001
C
      MINS=1.0E30
      IMIN=500
      DO 10 I=JC/12,JC/4
        SUM=0.0
        K=0
        DO 20 J=MAX(1,I-IW),MIN(JC,I+IW)
          SUM=SUM+PN(I)/N0
          K=K+1
20      CONTINUE
        IF (K.NE.0) SUM=SUM/K
        IF (SUM.LT.MINS) THEN
          MINS=SUM
          IMIN=I
        ENDIF
10    CONTINUE
      IMIN=1.1*IMIN
C
      MINS=1.0E30
      DO 30 I=6,14
        S(1)=REAL(I)/10.0
        DO 30 J=5,20
          S(2)=REAL(J)/100.0
          DO 30 K=3,13
            S(3)=K
            DO 30 L=1,20
              S(4)=REAL(L)/100.0
              D=EVAL_SINCF(JC,PN,N0,IMIN,S)
              IF (D.LT.MINS) THEN
                MINS=D
                SF(1)=S(1)
                SF(2)=S(2)
                SF(3)=S(3)
                SF(4)=S(4)
              ENDIF
30    CONTINUE
C
      CALL VA04A(JC,IMIN,N0,0.0D0,0.0D0,D,D,D,D,D,D,D,
     +           D,D,PN,0,SF,E,4,R,ESCALE,0,1,50,2)
      IT=1
C
C Also test logistic function fit
      S(1)=0.2
      S(2)=0.5
      S(3)=0.2
      S(4)=0.9
C
      MINS=1.0E30
      DO 31 I=1,9
        S(1)=REAL(I)/10.0
        DO 31 J=1,9
          S(2)=REAL(J)/10.0
          DO 31 K=1,9
            S(3)=REAL(K)/10.0
            DO 31 L=5,15
              S(4)=REAL(L)/10.0
              D=EVAL_LOGIST(JC,PN,N0,IMIN,S)
              IF (D.LT.MINS) THEN
                MINS=D
                S2(1)=S(1)
                S2(2)=S(2)
                S2(3)=S(3)
                S2(4)=S(4)
              ENDIF
31    CONTINUE
C
      CALL VA04A(JC,IMIN,N0,0.0D0,0.0D0,D,D,D,D,D,D,D,
     +           D,D,PN,0,S2,E,4,MINS,ESCALE,0,1,50,3)
C
      WRITE(*,*)
      IF (R.GT.MINS) THEN
        SF(1)=S2(1)
        SF(2)=S2(2)
        SF(3)=S2(3)
        SF(4)=S2(4)
        R=MINS
        IT=2
        WRITE(*,1010)(SF(I),I=1,4),R
1010    FORMAT(' Logistic function fit parameters, Chi2:',5F10.6)
      ELSE
        WRITE(*,1000)(SF(I),I=1,4),R
1000    FORMAT(' Sinc function fit parameters, Chi2:',5F10.6)
      ENDIF
C
      WRITE(*,*)
C
      RETURN
      END
C*********************************************************************
      SUBROUTINE FIT_EDGE(NSAM,EDG2,THRESH,D1,D2,LCOUNT,SF,DATA,
     +                    DATQ,TH,DIFF,NFIT,X,E,BUF,BUQ,
     +                    IPDAT,IPDAQ,IPF,IT)
C
      IMPLICIT NONE
C
      INTEGER NSAM,EDG2,NFIT,I,J,SLEN2,IPF,L,M,NP,ID,IS,IT
      INTEGER J1,J2
      REAL TH(*),DIFF(*),BUF(*),DATA(*),THRESH,F,IPDAT(*)
      REAL E(2*NFIT),ESCALE,X(2*NFIT),R,S,SUM,SCAL,SF(*),T2
      DOUBLE PRECISION D1,D2,DD,DL,SUM1,SUM2
      COMPLEX BUQ(*),DATQ(*),IPDAQ(*)
      CHARACTER*200 STRG
      LOGICAL LCOUNT
      DATA ESCALE/1000.0/
C
      NP=NSAM*IPF
C
      DD=(D2-D1)/2.0
      DL=DD/100.0
      DO 30 I=1,NSAM
        DO 30 J=1,NSAM
          ID=J+NSAM*(I-1)
          DO 40 L=-2,2
            DO 40 M=-2,2
            IF ((I+L.GE.1).AND.(I+L.LE.NSAM).AND.
     +           (J+M.GE.1).AND.(J+M.LE.NSAM)) THEN
              IS=J+M+NSAM*(I+L-1)
              IF (ABS(TH(ID)-TH(IS)).GT.DD) TH(ID)=TH(ID)+DL
            ENDIF
40        CONTINUE
30    CONTINUE
C
      DL=DL/2.0
      SUM1=0.0
      J1=0
      SUM2=0.0
      J2=0
      DO 50 I=1,NSAM
        DO 50 J=1,NSAM
          ID=J+NSAM*(I-1)
          IF ((ABS(TH(ID)-D1).GE.DL)
     +      .AND.(ABS(TH(ID)-D2).GE.DL)) THEN
            TH(ID)=1.0
          ELSEIF (ABS(TH(ID)-D1).LT.DL) THEN
            TH(ID)=0.0
            SUM1=SUM1+DATA(ID)
            J1=J1+1
          ELSE
            TH(ID)=2.0
            SUM2=SUM2+DATA(ID)
            J2=J2+1
          ENDIF
50    CONTINUE
      SUM1=SUM1/J1
      SUM2=SUM2/J2
      T2=(SUM2-SUM1)/2.0+SUM1
C
      CALL BILIN_INTERPO(NSAM,DATA,IPDAT,IPF)
C
      DO 60 I=0,NP-1
        L=INT(REAL(I)/IPF)
        DO 60 J=0,NP-1
          M=INT(REAL(J)/IPF)
          ID=J+1+NP*I
C          IS=M+1+NSAM*L
          IF (IPDAT(ID).LT.T2) THEN
            IPDAT(ID)=SUM1
          ELSE
            IPDAT(ID)=SUM2
          ENDIF
60    CONTINUE
C
      SCAL=2.0/NSAM/NSAM
      DO 20 I=1,NSAM*NSAM
        DATA(I)=(DATA(I)-T2)*SCAL
20    CONTINUE
      SCAL=2.0/NP/NP
      DO 70 I=1,NP*NP
        IPDAT(I)=(IPDAT(I)-T2)*SCAL
70    CONTINUE
C
      CALL RLFT3(DATA,DATQ,NSAM,NSAM,1,1)
      CALL RLFT3(IPDAT,IPDAQ,NP,NP,1,1)
C
      IF (LCOUNT) CALL APPLY_CNT_MTF(NSAM,DATA,DATQ,SF,IT)
      F=(REAL(EDG2)/NSAM)**2
C      CALL HIGHPASS(NSAM,DATA,DATQ,F)
      CALL RLFT3(DATA,DATQ,NSAM,NSAM,1,-1)
C
      SUM1=0.0
      J1=0
      SUM2=0.0
      J2=0
      DO 51 I=1+EDG2,NSAM-EDG2
        DO 51 J=1+EDG2,NSAM-EDG2
          ID=J+NSAM*(I-1)
          IF (TH(ID).EQ.0.0) THEN
            SUM1=SUM1+DATA(ID)
            J1=J1+1
          ELSEIF (TH(ID).EQ.2.0) THEN
            SUM2=SUM2+DATA(ID)
            J2=J2+1
          ENDIF
51    CONTINUE
      SUM1=SUM1/J1
      SUM2=SUM2/J2
      T2=(SUM2-SUM1)/2.0+SUM1
C
C      CALL HIGHPASS(NP,IPDAT,IPDAQ,F)
C
      DO 10 I=1,NFIT
        E(2*I-1)=0.01
        E(2*I)=0.01
        X(2*I-1)=10.0*(NFIT-I+1)/REAL(NFIT)
        X(2*I)=1/REAL(NFIT)
10    CONTINUE
C
      WRITE(*,*)
      WRITE(STRG,1100)('   Coeff-',I,'   Width-',I,I=1,NFIT)
1100  FORMAT(10(A9,I1))
      STRG=STRG(1:SLEN2(STRG))//'   NormRMS'
      WRITE(*,1200)STRG(1:SLEN2(STRG))
1200  FORMAT(A)
C
      CALL VA04A(NSAM,EDG2,0.0,SUM1,SUM2,DATA,
     +  DATQ,IPDAT,IPDAQ,TH,
     +  T2,DIFF,BUF,BUQ,SUM1,IPF,X,E,
     +  2*NFIT,R,ESCALE,0,1,50,1)
C
C      WRITE(*,*)
C      WRITE(*,*)'Intermediate values:'
C      WRITE(*,*)
C      WRITE(*,1200)STRG(1:SLEN2(STRG))
C      WRITE(*,1000)(ABS(X(2*I)),ABS(X(2*I-1)),I=1,NFIT),R
      DO 80 I=1,NFIT
        E(2*I-1)=0.001
        E(2*I)=0.001
80    CONTINUE
      CALL VA04A(NSAM,EDG2,0.0,SUM1,SUM2,DATA,
     +  DATQ,IPDAT,IPDAQ,TH,
     +  T2,DIFF,BUF,BUQ,SUM1,IPF,X,E,
     +  2*NFIT,R,ESCALE,0,1,50,1)
      WRITE(*,*)
      WRITE(*,*)'Final values:'
      WRITE(*,*)
      WRITE(*,1200)STRG(1:SLEN2(STRG))
      WRITE(*,1000)(ABS(X(2*I)),ABS(X(2*I-1)),I=1,NFIT),EXP(R)
1000  FORMAT(11F10.5)
      WRITE(*,*)
C
      RETURN
      END
C*********************************************************************
      SUBROUTINE CALCFX(N,X,R,NSAM,EDG2,N0,D1,D2,DATA,DATQ,
     +         TH,THQ,THM,THRESH,DIFF,BUF,BUQ,PN,IPF,IFLAG)
C
      IMPLICIT NONE
C
      INTEGER N,NSAM,EDG2,NFIT,I,NF,J,ID,IFLAG,IPF,NP,J1,J2
      REAL X(*),R,TH(*),DIFF(*),BUF(*),DATA(*),N0,RMIN
      REAL EVAL_SINCF,EVAL_LOGIST,THM(*),THRESH,SCAL
      DOUBLE PRECISION SUM,D1,D2,PN(*),SUM2
      COMPLEX BUQ(*),DATQ(*),THQ(*)
      SAVE RMIN
      DATA RMIN/999999999.0/
C
      IF (IFLAG.EQ.1) THEN
C
      NF=N/2
      NP=NSAM*IPF
      DO 80 I=1,NP*NP
        BUF(I)=TH(I)
80    CONTINUE
      DO 90 I=1,NP
        BUQ(I)=THQ(I)
90    CONTINUE
      CALL GFILTER(NP,IPF,NF,BUF,BUQ,X)
      CALL RLFT3(BUF,BUQ,NP,NP,1,-1)
      CALL BINNING(NP,BUF,BUF,IPF)
C
      SUM=0.0
      J1=0
      SUM2=0.0
      J2=0
      DO 60 I=1+EDG2,NSAM-EDG2
        DO 60 J=1+EDG2,NSAM-EDG2
          ID=J+NSAM*(I-1)
          IF (THM(ID).EQ.0.0) THEN
            SUM=SUM+BUF(ID)
            J1=J1+1
          ELSEIF (THM(ID).EQ.2.0) THEN
            SUM2=SUM2+BUF(ID)
            J2=J2+1
          ENDIF
60    CONTINUE
      SUM=SUM/J1
      SUM2=SUM2/J2
C
      SCAL=(D2-D1)/(SUM2-SUM)
      DO 70 I=1,NSAM*NSAM
        BUF(I)=(BUF(I)-SUM)*SCAL+D1-THRESH
70    CONTINUE

C
      DO 40 I=1,NSAM*NSAM
        DIFF(I)=DATA(I)-BUF(I)
40    CONTINUE
C
      SUM=0.0
      DO 50 I=1+EDG2,NSAM-EDG2
        DO 50 J=1+EDG2,NSAM-EDG2
          ID=J+NSAM*(I-1)
          SUM=SUM+DIFF(ID)**2
50    CONTINUE
      R=LOG(SQRT(SUM/(NSAM-2*EDG2)**2)/(D2-D1))
cc        WRITE(*,1000)(ABS(X(2*I)),ABS(X(2*I-1)),I=1,NF),EXP(R)
      IF (R.LT.RMIN) THEN
        WRITE(*,1000)(ABS(X(2*I)),ABS(X(2*I-1)),I=1,NF),EXP(R)
1000    FORMAT(11F10.5)
        FLUSH(6)
        RMIN=R
      ENDIF
C
      ELSEIF (IFLAG.EQ.2) THEN
C
        R=EVAL_SINCF(NSAM,PN,N0,EDG2,X)
C
      ELSE
C
        R=EVAL_LOGIST(NSAM,PN,N0,EDG2,X)
C
      ENDIF
C
      RETURN
      END
C*********************************************************************
      REAL FUNCTION EVAL_LOGIST(NSAM,PN,N0,EDG2,X)
C
      IMPLICIT NONE
C
      INTEGER NSAM,I,EDG2
      REAL N0,LOGIST,S1,F,X(*)
      DOUBLE PRECISION SUM,S,PN(*)
C
      SUM=0.0
      S1=0.0
      DO 110 I=1+EDG2,NSAM-10
C        S1=S1+I**2
        S1=S1+SQRT(REAL(NSAM-I))
        F=REAL(I-1)/(NSAM-1)
        S=LOGIST(X,F)
C        SUM=SUM+I**2*(PN(I)/N0-S)**2
        SUM=SUM+SQRT(REAL(NSAM-I))*(PN(I)/N0-S)**2
110   CONTINUE
C
      EVAL_LOGIST=SQRT(SUM/S1)
C
      RETURN
      END
C*********************************************************************
      REAL FUNCTION EVAL_SINCF(NSAM,PN,N0,EDG2,X)
C
      IMPLICIT NONE
C
      INTEGER NSAM,I,EDG2
      REAL N0,SINCF,S1,F,X(*)
      DOUBLE PRECISION SUM,S,PN(*)
C
      SUM=0.0
      S1=0.0
      DO 110 I=1+EDG2,NSAM-10
        S1=S1+SQRT(REAL(NSAM-I))
        F=REAL(I-1)/(NSAM-1)
        S=SINCF(X,F)
        SUM=SUM+SQRT(REAL(NSAM-I))*(PN(I)/N0-S)**2
110   CONTINUE
C
      EVAL_SINCF=SQRT(SUM/S1)
C
      RETURN
      END
C*********************************************************************
      REAL FUNCTION SINCF(SF,R)
C
      IMPLICIT NONE
C
      REAL SF(*),R,SINC,GAUSSIAN
C
      SINCF=(SF(1)-ABS(SF(2))*SINC(SF(3)*R))
     +      *GAUSSIAN(R,SF(4))
C
      RETURN
      END
C*********************************************************************
      REAL FUNCTION LOGIST(SF,R)
C
      IMPLICIT NONE
C
      REAL SF(*),R
C
      LOGIST=SF(1)/(1.0+EXP((SF(2)-R)/SF(3)))+SF(4)
C
      RETURN
      END
C*********************************************************************
      REAL FUNCTION GAUSSIAN(X,G)
C
      IMPLICIT NONE
C
      REAL X,G
C
      GAUSSIAN=EXP(-G**2*X**2)
C
      RETURN
      END
C*********************************************************************
      SUBROUTINE APPLY_CNT_MTF(NSAM,BUFC,BUFQ,SF,IT)
C
      IMPLICIT NONE
C
      INTEGER NSAM,L,M,LL,MM,NSAMH,JC,ID,IT
      REAL F,S,E,EL,SF(*),R,A,SINCF,LOGIST
      COMPLEX BUFC(*),BUFQ(*)
C
      NSAMH=NSAM/2
      JC=NSAMH+1
      F=(1.0/NSAMH)**2
      IF (IT.EQ.1) THEN
        A=SINCF(SF,1.0)
      ELSE
        A=LOGIST(SF,1.0)
      ENDIF
C
      DO 10 L=1,JC
        LL=L-1
        EL=LL**2
        DO 10 M=1,NSAM
          MM=M-1
          IF (MM.GE.JC) MM=MM-NSAM
          E=EL+MM**2
          R=SQRT(E*F)
          IF (IT.EQ.1) THEN
            S=SINCF(SF,R)/A
          ELSE
            S=LOGIST(SF,R)/A
          ENDIF
          IF (L.NE.JC) THEN
            ID=L+NSAMH*(M-1)
            BUFC(ID)=BUFC(ID)/SQRT(S)
          ELSE
            BUFQ(M)=BUFQ(M)/SQRT(S)
          ENDIF
10    CONTINUE
C
      RETURN
      END
C*********************************************************************
      SUBROUTINE HIGHPASS(NSAM,BUFC,BUFQ,F)
C
      IMPLICIT NONE
C
      INTEGER NSAM,L,M,LL,MM,NSAMH,JC,ID
      REAL F,W,E,EL
      COMPLEX BUFC(*),BUFQ(*)
C
      NSAMH=NSAM/2
      JC=NSAMH+1
C
      DO 10 L=1,JC
        LL=L-1
        EL=LL**2
        DO 10 M=1,NSAM
          MM=M-1
          IF (MM.GE.JC) MM=MM-NSAM
          E=EL+MM**2
          W=1.0-EXP(-E*F)
          IF (L.NE.JC) THEN
            ID=L+NSAMH*(M-1)
            BUFC(ID)=BUFC(ID)*W
          ELSE
            BUFQ(M)=BUFQ(M)*W
          ENDIF
10    CONTINUE
C
      RETURN
      END
C*********************************************************************
      SUBROUTINE SINCFILT(NSAM,BUFC,BUFQ)
C
      IMPLICIT NONE
C
      INTEGER NSAM,L,M,LL,MM,NSAMH,JC,ID
      REAL W,E,EL,PI, SINC
      PARAMETER (PI=3.1415926535897)
      COMPLEX BUFC(*),BUFQ(*)
C
      NSAMH=NSAM/2
      JC=NSAMH+1
C
      DO 10 L=1,JC
        LL=L-1
        EL=LL**2
        DO 10 M=1,NSAM
          MM=M-1
          IF (MM.GE.JC) MM=MM-NSAM
          E=EL+MM**2
          W=SINC(PI*SQRT(E)/NSAM)
C          W=1.0-EXP(-E*F)
          IF (L.NE.JC) THEN
            ID=L+NSAMH*(M-1)
            BUFC(ID)=BUFC(ID)/W
          ELSE
            BUFQ(M)=BUFQ(M)/W
          ENDIF
10    CONTINUE
C
      RETURN
      END
C*********************************************************************
      SUBROUTINE GFILTER(NSAM,IPF,NF,THC,THQ,X)
C
      IMPLICIT NONE
C
      INTEGER NSAM,L,M,LL,MM,NSAMH,JC,ID,NF,J,IPF
      REAL W,E,EL,X(*),S,SCAL,F(10),A(10)
      COMPLEX THC(*),THQ(*)
C
      NSAMH=NSAM/2
      JC=NSAMH+1
      SCAL=1.0/NSAM/NSAM*IPF*IPF
C
      S=0.0
      DO 60 J=1,NF
        S=S+ABS(X(2*J))
60    CONTINUE
      DO 30 J=1,NF
        F(J)=X(2*J-1)**2*SCAL
        A(J)=ABS(X(2*J))/S
30    CONTINUE
C
      DO 10 L=1,JC
        LL=L-1
        EL=LL**2
        DO 10 M=1,NSAM
          MM=M-1
          IF (MM.GE.JC) MM=MM-NSAM
          E=EL+MM**2
          W=0.0
          DO 20 J=1,NF
            W=W+A(J)*EXP(-E*F(J))
20        CONTINUE
          IF (L.NE.JC) THEN
            ID=L+NSAMH*(M-1)
            THC(ID)=THC(ID)*W
          ELSE
            THQ(M)=THQ(M)*W
          ENDIF
10    CONTINUE
C
      RETURN
      END
C*********************************************************************
      SUBROUTINE PCHECK(N)
C
      IMPLICIT NONE
C
      INTEGER M,N
C
99    CONTINUE
      M=N
      CALL TESTPRIME(M)
C
      IF (M.NE.1) THEN
        N=N+2
        GOTO 99
      ENDIF
C
      RETURN
      END
C*********************************************************************
      SUBROUTINE TESTPRIME(M)
C
      IMPLICIT NONE
C
      INTEGER IPRIME(10),M,I
      DATA IPRIME/2,3,5,7,11,13,17,19,23,29/
C
      DO 10 I=1,10
20      CONTINUE
        IF (MOD(M,IPRIME(I)).EQ.0) THEN
          M=M/IPRIME(I)
          GOTO 20
        ENDIF
10    CONTINUE
C
      RETURN
      END
C*********************************************************************
      REAL FUNCTION SINC(X)
C
      IMPLICIT NONE
C
      REAL X
C
      IF (X.NE.0.0) THEN
        SINC=SIN(X)/X
      ELSE
        SINC=1.0
      ENDIF
C
      RETURN   
      END
C*********************************************************************
      SUBROUTINE PVAR(NSAM,D2,N,DMEAN,DRMS,VB,EDG2)
C
      IMPLICIT NONE
C
      INTEGER NSAM,I,J,K,L,IPATCH,IP,ID,IS,EDG2
      PARAMETER (IPATCH=20)
      REAL N(*),V,VB
      DOUBLE PRECISION D2,DMEAN,DRMS,CC,S1,S2
C      DOUBLE PRECISION DM1,DM2
C
      IP=0
      DMEAN=0.0
      DRMS=0.0
      DO 20 I=1+EDG2,NSAM-IPATCH-EDG2
        DO 20 J=1+EDG2,NSAM-IPATCH-EDG2
          ID=J+NSAM*(I-1)
          IF (N(ID).NE.D2) THEN
            DMEAN=DMEAN+N(ID)
            DRMS=DRMS+N(ID)**2
            IP=IP+1
          ENDIF
20    CONTINUE
      DMEAN=DMEAN/IP
      DRMS=SQRT(DRMS/IP-DMEAN**2)
C
      VB=0.0
      DO 31 K=0,IPATCH-1
        DO 31 L=1,IPATCH-1
          IF ((K.EQ.IPATCH-1).OR.(L.EQ.IPATCH-1)) THEN
            CC=0.0
            S1=0.0
            S2=0.0
            DO 41 I=1+EDG2,NSAM-IPATCH-EDG2
              DO 41 J=1+EDG2,NSAM-IPATCH-EDG2
                ID=J+NSAM*(I-1)
                IS=J+L+NSAM*(I-1+K)
                IF ((N(ID).NE.D2).AND.(N(IS).NE.D2)) THEN
                  CC=CC+(N(ID)-DMEAN)*(N(IS)-DMEAN)
                  S1=S1+(N(ID)-DMEAN)**2
                  S2=S2+(N(IS)-DMEAN)**2
                ENDIF
41          CONTINUE
            CC=CC/SQRT(S1)/SQRT(S2)
            VB=VB+CC
          ENDIF
31    CONTINUE
      VB=VB/(2*IPATCH-1)
C
      V=0.0
      DO 30 K=0,IPATCH-1
        DO 30 L=1,IPATCH-1
          CC=0.0
          S1=0.0
          S2=0.0
          DO 40 I=1+EDG2,NSAM-IPATCH-EDG2
            DO 40 J=1+EDG2,NSAM-IPATCH-EDG2
              ID=J+NSAM*(I-1)
              IS=J+L+NSAM*(I-1+K)
              IF ((N(ID).NE.D2).AND.(N(IS).NE.D2)) THEN
                CC=CC+(N(ID)-DMEAN)*(N(IS)-DMEAN)
                S1=S1+(N(ID)-DMEAN)**2
                S2=S2+(N(IS)-DMEAN)**2
              ENDIF
40        CONTINUE
          CC=CC/SQRT(S1)/SQRT(S2)
          IF (CC.GT.3.0*VB) V=V+CC
30    CONTINUE
C
      DRMS=DRMS*SQRT(1.0+8.0*V)
C
      RETURN
      END
C*********************************************************************
      SUBROUTINE VA04A(NSAM,EDG2,N0,D1,D2,DATA,DATQ,TH,
     +                 THQ,THM,THRESH,DIFF,BUF,BUQ,PN,IPF,
     +                 X,E,N,F,ESCALE,IPRINT,ICON,MAXIT,IFLAG)
C**************************************************************************
C  STANDARD FORTRAN 66 (A VERIFIED PFORT SUBROUTINE)  - Powell minimisation
C  some changes were made to reduce diagnostic output and prevent occasional
C  crashes (Niko, 12. June 1998)
C  Calls CALCFX
C**************************************************************************
      DIMENSION W(200),X(*),E(*),XS(10)
      REAL DATA(*),TH(*),DIFF(*),BUF(*),N0,THM(*),THRESH
      DOUBLE PRECISION D1,D2,PN(*)
      COMPLEX BUQ(*),DATQ(*),THQ(*)
C**************************************************************************
C	W[N*(N+3)]
      DDMAG=0.1*ESCALE
      SCER=0.05/ESCALE
      ICNT=0
      MAXX=100*MAXIT
      DO 999 I=1,N
        XS(I)=X(I)
  999 CONTINUE
      JJ=N*N+N
      JJJ=JJ+N
      K=N+1
      NFCC=1
      IND=1
      INN=1
      DO 1 I=1,N
      DO 2 J=1,N
      W(K)=0.
      IF(I-J)4,3,4
    3 W(K)=ABS(E(I))
      W(I)=ESCALE
    4 K=K+1
    2 CONTINUE
    1 CONTINUE
      ITERC=1
      ISGRAD=2
      ICNT=ICNT+1
      IF (ICNT.GT.MAXX) GOTO 998
      CALL CALCFX(N,X,F,NSAM,EDG2,N0,D1,D2,DATA,DATQ,TH,THQ,
     +            THM,THRESH,DIFF,BUF,BUQ,PN,IPF,IFLAG)
      FKEEP=ABS(F)+ABS(F)
    5 ITONE=1
      FP=F
      SUM=0.
      IXP=JJ
      DO 6 I=1,N
      IXP=IXP+1
      W(IXP)=X(I)
    6 CONTINUE
      IDIRN=N+1
      ILINE=1
    7 DMAX=W(ILINE)
      DACC=DMAX*SCER
      DMAG=AMIN1(DDMAG,0.1*DMAX)
      DMAG=AMAX1(DMAG,20.*DACC)
      DDMAX=10.*DMAG
      GO TO (70,70,71),ITONE
   70 DL=0.
      D=DMAG
      FPREV=F
      IS=5
      FA=F
      DA=DL
    8 DD=D-DL
      DL=D
   58 K=IDIRN
      DO 9 I=1,N
      X(I)=X(I)+DD*W(K)
      K=K+1
    9 CONTINUE
      ICNT=ICNT+1
      IF (ICNT.GT.MAXX) GOTO 998
      CALL CALCFX(N,X,F,NSAM,EDG2,N0,D1,D2,DATA,DATQ,TH,THQ,
     +            THM,THRESH,DIFF,BUF,BUQ,PN,IPF,IFLAG)
      NFCC=NFCC+1
      GO TO (10,11,12,13,14,96),IS
   14 IF(F-FA)15,16,24
   16 IF (ABS(D)-DMAX) 17,17,18
   17 D=D+D
      GO TO 8
   18 CONTINUE
C      WRITE(6,19)
   19 FORMAT(5X,50HPOWELL MIN: MAXIMUM CHANGE DOES NOT ALTER FUNCTION)
C   19 FORMAT(5X,44HVA04A MAXIMUM CHANGE DOES NOT ALTER FUNCTION)
      GO TO 20
   15 FB=F
      DB=D
      GO TO 21
   24 FB=FA
      DB=DA
      FA=F
      DA=D
   21 GO TO (83,23),ISGRAD
   23 D=DB+DB-DA
      IS=1
      GO TO 8
   83 D=0.5*(DA+DB-(FA-FB)/(DA-DB))
      IS=4
      IF((DA-D)*(D-DB))25,8,8
   25 IS=1
      IF(ABS(D-DB)-DDMAX)8,8,26
   26 D=DB+SIGN(DDMAX,DB-DA)
      IS=1
      DDMAX=DDMAX+DDMAX
      DDMAG=DDMAG+DDMAG
      IF(DDMAX-DMAX)8,8,27
   27 DDMAX=DMAX
      GO TO 8
   13 IF(F-FA)28,23,23
   28 FC=FB
      DC=DB
   29 FB=F
      DB=D
      GO TO 30
   12 IF(F-FB)28,28,31
   31 FA=F
      DA=D
      GO TO 30
   11 IF(F-FB)32,10,10
   32 FA=FB
      DA=DB
      GO TO 29
   71 DL=1.
      DDMAX=5.
      FA=FP
      DA=-1.
      FB=FHOLD
      DB=0.
      D=1.
   10 FC=F
      DC=D
   30 A=(DB-DC)*(FA-FC)
      B=(DC-DA)*(FB-FC)
      IF((A+B)*(DA-DC))33,33,34
   33 FA=FB
      DA=DB
      FB=FC
      DB=DC
      GO TO 26
   34 D=0.5*(A*(DB+DC)+B*(DA+DC))/(A+B)
      DI=DB
      FI=FB
      IF(FB-FC)44,44,43
   43 DI=DC
      FI=FC
   44 GO TO (86,86,85),ITONE
   85 ITONE=2
      GO TO 45
   86 IF (ABS(D-DI)-DACC) 41,41,93
   93 IF (ABS(D-DI)-0.03*ABS(D)) 41,41,45
   45 IF ((DA-DC)*(DC-D)) 47,46,46
   46 FA=FB
      DA=DB
      FB=FC
      DB=DC
      GO TO 25
   47 IS=2
      IF ((DB-D)*(D-DC)) 48,8,8
   48 IS=3
      GO TO 8
   41 F=FI
      D=DI-DL
      DD=(DC-DB)*(DC-DA)*(DA-DB)/(A+B)
      IF (DD.LT.0.0) DD=1E-10
      DD=SQRT(DD)
      DO 49 I=1,N
      X(I)=X(I)+D*W(IDIRN)
      W(IDIRN)=DD*W(IDIRN)
      IDIRN=IDIRN+1
   49 CONTINUE
      IF (DD.EQ.0.0) DD=1E-10
      W(ILINE)=W(ILINE)/DD
      ILINE=ILINE+1
      IF(IPRINT-1)51,50,51
   50 CONTINUE
C   50 WRITE(6,52) ITERC,NFCC,F,(X(I),I=1,N)
   52 FORMAT (/1X,9HITERATION,I5,I15,16H FUNCTION VALUES,
     110X,3HF =,E21.14/(5E24.14))
      GO TO(51,53),IPRINT
   51 GO TO (55,38),ITONE
   55 IF (FPREV-F-SUM) 94,95,95
   95 SUM=FPREV-F
      JIL=ILINE
   94 IF (IDIRN-JJ) 7,7,84
   84 GO TO (92,72),IND
   92 FHOLD=F
      IS=6
      IXP=JJ
      DO 59 I=1,N
      IXP=IXP+1
      W(IXP)=X(I)-W(IXP)
   59 CONTINUE
      DD=1.
      GO TO 58
   96 GO TO (112,87),IND
  112 IF (FP-F) 37,37,91
   91 D=2.*(FP+F-2.*FHOLD)/(FP-F)**2
      IF (D*(FP-FHOLD-SUM)**2-SUM) 87,37,37
   87 J=JIL*N+1
      IF (J-JJ) 60,60,61
   60 DO 62 I=J,JJ
      K=I-N
      W(K)=W(I)
   62 CONTINUE
      DO 97 I=JIL,N
      W(I-1)=W(I)
   97 CONTINUE
   61 IDIRN=IDIRN-N
      ITONE=3
      K=IDIRN
      IXP=JJ
      AAA=0.
      DO 65 I=1,N
      IXP=IXP+1
      W(K)=W(IXP)
      IF (AAA-ABS(W(K)/E(I))) 66,67,67
   66 AAA=ABS(W(K)/E(I))
   67 K=K+1
   65 CONTINUE
      DDMAG=1.
      IF (AAA.EQ.0.0) AAA=1E-10
      W(N)=ESCALE/AAA
      ILINE=N
      GO TO 7
   37 IXP=JJ
      AAA=0.
      F=FHOLD
      DO 99 I=1,N
      IXP=IXP+1
      X(I)=X(I)-W(IXP)
      IF (AAA*ABS(E(I))-ABS(W(IXP))) 98,99,99
   98 AAA=ABS(W(IXP)/E(I))
   99 CONTINUE
      GO TO 72
   38 AAA=AAA*(1.+DI)
      GO TO (72,106),IND
   72 IF (IPRINT-2) 53,50,50
   53 GO TO (109,88),IND
  109 IF (AAA-0.1) 89,89,76
   89 GO TO (20,116),ICON
  116 IND=2
      GO TO (100,101),INN
  100 INN=2
      K=JJJ
      DO 102 I=1,N
      K=K+1
      W(K)=X(I)
      X(I)=X(I)+10.*E(I)
  102 CONTINUE
      FKEEP=F
      ICNT=ICNT+1
      IF (ICNT.GT.MAXX) GOTO 998
      CALL CALCFX(N,X,F,NSAM,EDG2,N0,D1,D2,DATA,DATQ,TH,THQ,
     +            THM,THRESH,DIFF,BUF,BUQ,PN,IPF,IFLAG)
      NFCC=NFCC+1
      DDMAG=0.
      GO TO 108
   76 IF (F-FP) 35,78,78
   78 CONTINUE
C   78 WRITE(6,80)
   80 FORMAT (5X,43HPOWELL MIN: ACCURACY LIMITED BY ERRORS IN F)
C   80 FORMAT (5X,37HVA04A ACCURACY LIMITED BY ERRORS IN F)
      GO TO 20
   88 IND=1
   35 TMP=FP-F
      IF (TMP.GT.0.0) THEN
      DDMAG=0.4*SQRT(TMP)
      ELSE
      DDMAG=0.0
      ENDIF
      ISGRAD=1
  108 ITERC=ITERC+1
      IF (ITERC-MAXIT) 5,5,81
81    CONTINUE
C   81 WRITE(6,82) MAXIT
   82 FORMAT(I5,30H ITERATIONS COMPLETED BY VA04A)
      IF (F-FKEEP) 20,20,110
  110 F=FKEEP
      DO 111 I=1,N
      JJJ=JJJ+1
      X(I)=W(JJJ)
  111 CONTINUE
      GO TO 20
  101 JIL=1
      FP=FKEEP
      IF (F-FKEEP) 105,78,104
  104 JIL=2
      FP=F
      F=FKEEP
  105 IXP=JJ
      DO 113 I=1,N
      IXP=IXP+1
      K=IXP+N
      GO TO (114,115),JIL
  114 W(IXP)=W(K)
      GO TO 113
  115 W(IXP)=X(I)
      X(I)=W(K)
  113 CONTINUE
      JIL=2
      GO TO 92
  106 IF (AAA-0.1) 20,20,107
   20 CONTINUE
      RETURN
  107 INN=1
      GO TO 35
  998 CONTINUE
      DO 997 I=1,N
        X(I)=XS(I)
  997 CONTINUE
      CALL CALCFX(N,X,F,NSAM,EDG2,N0,D1,D2,DATA,DATQ,TH,THQ,
     +            THM,THRESH,DIFF,BUF,BUQ,PN,IPF,IFLAG)
      PRINT *,'VA04A ENDLESS LOOP SAFETY CATCH: ICNT = ',ICNT
      RETURN
      END
C**************************************************************************
      SUBROUTINE FIT_SPLINES(N,R,N0,NPS,W,KNOTS,NKNOT,IPAD,
     +             IPATCH,CSPLINE,SFIT,NPS0,NEST,IWRK,WRK)
C
      IMPLICIT NONE
C
      INTEGER N,I,NEST,IER,NKNOT,IPAD,ISTART,IPATCH,LWRK
      INTEGER IWRK(*),J,IEND
C      PARAMETER (ISTART=5)
      PARAMETER (IEND=50)
      REAL R(*),W(*),S,KNOTS(*),CSPLINE(*),CH2,SFIT(*)
      REAL WRK(*),B,NPS0,N0
      DOUBLE PRECISION NPS(*)
C
      ISTART=0
      DO 70 I=1,N+IPAD
        R(I)=REAL(I-1)/(N-1)
        IF ((ISTART.EQ.0).AND.(R(I).GT.1.0/IPATCH)) ISTART=I
70    CONTINUE
      LWRK=(N+IPAD)*4+NEST*16
      DO 10 I=1,N+IPAD
C        W(I)=LOG(1.0+I)
C        W(I)=SQRT(REAL(I))
        W(I)=I
        IF (W(I).LT.20.0) W(I)=20.0
        IF (I.LE.ISTART) THEN
           S=0.0
           DO 100 J=MAX(I,2),ISTART
             IF (NPS(J)/NPS(ISTART).GT.20.0) THEN
               S=S+20.0*NPS(ISTART)/N0
             ELSE
               S=S+NPS(J)/N0
             ENDIF
100        CONTINUE
           SFIT(I)=S/(ISTART-MAX(I,2)+1)
C          SFIT(I)=(NPS(ISTART+1)+(NPS(ISTART+100)
C     +            -NPS(ISTART+1))/100.0*(I-ISTART-1))/N0
        ELSEIF ((I.LE.N-IEND).AND.(I.GT.ISTART)) THEN
          SFIT(I)=NPS(I)/N0
        ELSE
          SFIT(I)=(NPS(N-IEND)+(NPS(N-IEND)-NPS(N-IEND-200))
     +            /200.0*(I-N+IEND))/N0
        ENDIF
10    CONTINUE
C
      S=1.0
      B=REAL(N+IPAD)/(N-1)
      DO 80 I=5,NKNOT
        KNOTS(I)=B*((REAL(I-4)/(NKNOT-4))**3
     +             +REAL(I-4)/(NKNOT-4))/2.0
80    CONTINUE
      DO 81 I=1,4
        KNOTS(I)=B*REAL(I-1)/4.0*KNOTS(5)
81    CONTINUE
      CALL CURFIT(-1,N+IPAD,R,SFIT,W,0.0,B,3,S,NEST,NKNOT,
     +            KNOTS,CSPLINE,CH2,WRK,LWRK,IWRK,IER)
      IF (IER.NE.0) STOP 'ERROR: Spline fit failed'
      WRITE(*,*)
      WRITE(*,*) 'Number of knots in spline =',NKNOT
      WRITE(*,*) 'Chi2                      =',CH2
C
      CALL SPLEV(KNOTS,NKNOT,CSPLINE,3,R,SFIT,N,IER)
      IF (IER.NE.0) STOP 'ERROR: Spline evaluation failed'
C
      NPS0=SFIT(1)
      WRITE(*,*)
C
      END
C******************************************************************************
      subroutine splev(t,n,c,k,x,y,m,ier)
c  subroutine splev evaluates in a number of points x(i),i=1,2,...,m
c  a spline s(x) of degree k, given in its b-spline representation.
c
c  calling sequence:
c     call splev(t,n,c,k,x,y,m,ier)
c
c  input parameters:
c    t    : array,length n, which contains the position of the knots.
c    n    : integer, giving the total number of knots of s(x).
c    c    : array,length n, which contains the b-spline coefficients.
c    k    : integer, giving the degree of s(x).
c    x    : array,length m, which contains the points where s(x) must
c           be evaluated.
c    m    : integer, giving the number of points where s(x) must be
c           evaluated.
c
c  output parameter:
c    y    : array,length m, giving the value of s(x) at the different
c           points.
c    ier  : error flag
c      ier = 0 : normal return
c      ier =10 : invalid input data (see restrictions)
c
c  restrictions:
c    m >= 1
c    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
c
c  other subroutines required: fpbspl.
c
c  references :
c    de boor c  : on calculating with b-splines, j. approximation theory
c                 6 (1972) 50-62.
c    cox m.g.   : the numerical evaluation of b-splines, j. inst. maths
c                 applics 10 (1972) 134-149.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer n,k,m,ier
c  ..array arguments..
      real t(n),c(n),x(m),y(m)
c  ..local scalars..
      integer i,j,k1,l,ll,l1,nk1
      real arg,sp,tb,te
c  ..local array..
      real h(6)
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(m-1) 100,30,10
  10  do 20 i=2,m
        if(x(i).lt.x(i-1)) go to 100
  20  continue
  30  ier = 0
c  fetch tb and te, the boundaries of the approximation interval.
      k1 = k+1
      nk1 = n-k1
      tb = t(k1)
      te = t(nk1+1)
      l = k1
      l1 = l+1
c  main loop for the different points.
      do 80 i=1,m
c  fetch a new x-value arg.
        arg = x(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
c  search for knot interval t(l) <= arg < t(l+1)
  40    if(arg.lt.t(l1) .or. l.eq.nk1) go to 50
        l = l1
        l1 = l+1
        go to 40
c  evaluate the non-zero b-splines at arg.
  50    call fpbspl(t,n,k,arg,l,h)
c  find the value of s(x) at x=arg.
        sp = 0.
        ll = l-k1
        do 60 j=1,k1
          ll = ll+1
          sp = sp+c(ll)*h(j)
  60    continue
        y(i) = sp
  80  continue
 100  return
      end
C
      subroutine curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,
     * wrk,lwrk,iwrk,ier)
c  given the set of data points (x(i),y(i)) and the set of positive
c  numbers w(i),i=1,2,...,m,subroutine curfit determines a smooth spline
c  approximation of degree k on the interval xb <= x <= xe.
c  if iopt=-1 curfit calculates the weighted least-squares spline
c  according to a given set of knots.
c  if iopt>=0 the number of knots of the spline s(x) and the position
c  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
c  ness of s(x) is then achieved by minimalizing the discontinuity
c  jumps of the k-th derivative of s(x) at the knots t(j),j=k+2,k+3,...,
c  n-k-1. the amount of smoothness is determined by the condition that
c  f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, with s a given non-
c  negative constant, called the smoothing factor.
c  the fit s(x) is given in the b-spline representation (b-spline coef-
c  ficients c(j),j=1,2,...,n-k-1) and can be evaluated by means of
c  subroutine splev.
c
c  calling sequence:
c     call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,
c    * lwrk,iwrk,ier)
c
c  parameters:
c   iopt  : integer flag. on entry iopt must specify whether a weighted
c           least-squares spline (iopt=-1) or a smoothing spline (iopt=
c           0 or 1) must be determined. if iopt=0 the routine will start
c           with an initial set of knots t(i)=xb, t(i+k+1)=xe, i=1,2,...
c           k+1. if iopt=1 the routine will continue with the knots
c           found at the last call of the routine.
c           attention: a call with iopt=1 must always be immediately
c           preceded by another call with iopt=1 or iopt=0.
c           unchanged on exit.
c   m     : integer. on entry m must specify the number of data points.
c           m > k. unchanged on exit.
c   x     : real array of dimension at least (m). before entry, x(i)
c           must be set to the i-th value of the independent variable x,
c           for i=1,2,...,m. these values must be supplied in strictly
c           ascending order. unchanged on exit.
c   y     : real array of dimension at least (m). before entry, y(i)
c           must be set to the i-th value of the dependent variable y,
c           for i=1,2,...,m. unchanged on exit.
c   w     : real array of dimension at least (m). before entry, w(i)
c           must be set to the i-th value in the set of weights. the
c           w(i) must be strictly positive. unchanged on exit.
c           see also further comments.
c   xb,xe : real values. on entry xb and xe must specify the boundaries
c           of the approximation interval. xb<=x(1), xe>=x(m).
c           unchanged on exit.
c   k     : integer. on entry k must specify the degree of the spline.
c           1<=k<=5. it is recommended to use cubic splines (k=3).
c           the user is strongly dissuaded from choosing k even,together
c           with a small s-value. unchanged on exit.
c   s     : real.on entry (in case iopt>=0) s must specify the smoothing
c           factor. s >=0. unchanged on exit.
c           for advice on the choice of s see further comments.
c   nest  : integer. on entry nest must contain an over-estimate of the
c           total number of knots of the spline returned, to indicate
c           the storage space available to the routine. nest >=2*k+2.
c           in most practical situation nest=m/2 will be sufficient.
c           always large enough is  nest=m+k+1, the number of knots
c           needed for interpolation (s=0). unchanged on exit.
c   n     : integer.
c           unless ier =10 (in case iopt >=0), n will contain the
c           total number of knots of the spline approximation returned.
c           if the computation mode iopt=1 is used this value of n
c           should be left unchanged between subsequent calls.
c           in case iopt=-1, the value of n must be specified on entry.
c   t     : real array of dimension at least (nest).
c           on succesful exit, this array will contain the knots of the
c           spline,i.e. the position of the interior knots t(k+2),t(k+3)
c           ...,t(n-k-1) as well as the position of the additional knots
c           t(1)=t(2)=...=t(k+1)=xb and t(n-k)=...=t(n)=xe needed for
c           the b-spline representation.
c           if the computation mode iopt=1 is used, the values of t(1),
c           t(2),...,t(n) should be left unchanged between subsequent
c           calls. if the computation mode iopt=-1 is used, the values
c           t(k+2),...,t(n-k-1) must be supplied by the user, before
c           entry. see also the restrictions (ier=10).
c   c     : real array of dimension at least (nest).
c           on succesful exit, this array will contain the coefficients
c           c(1),c(2),..,c(n-k-1) in the b-spline representation of s(x)
c   fp    : real. unless ier=10, fp contains the weighted sum of
c           squared residuals of the spline approximation returned.
c   wrk   : real array of dimension at least (m*(k+1)+nest*(7+3*k)).
c           used as working space. if the computation mode iopt=1 is
c           used, the values wrk(1),...,wrk(n) should be left unchanged
c           between subsequent calls.
c   lwrk  : integer. on entry,lwrk must specify the actual dimension of
c           the array wrk as declared in the calling (sub)program.lwrk
c           must not be too small (see wrk). unchanged on exit.
c   iwrk  : integer array of dimension at least (nest).
c           used as working space. if the computation mode iopt=1 is
c           used,the values iwrk(1),...,iwrk(n) should be left unchanged
c           between subsequent calls.
c   ier   : integer. unless the routine detects an error, ier contains a
c           non-positive value on exit, i.e.
c    ier=0  : normal return. the spline returned has a residual sum of
c             squares fp such that abs(fp-s)/s <= tol with tol a relat-
c             ive tolerance set to 0.001 by the program.
c    ier=-1 : normal return. the spline returned is an interpolating
c             spline (fp=0).
c    ier=-2 : normal return. the spline returned is the weighted least-
c             squares polynomial of degree k. in this extreme case fp
c             gives the upper bound fp0 for the smoothing factor s.
c    ier=1  : error. the required storage space exceeds the available
c             storage space, as specified by the parameter nest.
c             probably causes : nest too small. if nest is already
c             large (say nest > m/2), it may also indicate that s is
c             too small
c             the approximation returned is the weighted least-squares
c             spline according to the knots t(1),t(2),...,t(n). (n=nest)
c             the parameter fp gives the corresponding weighted sum of
c             squared residuals (fp>s).
c    ier=2  : error. a theoretically impossible result was found during
c             the iteration proces for finding a smoothing spline with
c             fp = s. probably causes : s too small.
c             there is an approximation returned but the corresponding
c             weighted sum of squared residuals does not satisfy the
c             condition abs(fp-s)/s < tol.
c    ier=3  : error. the maximal number of iterations maxit (set to 20
c             by the program) allowed for finding a smoothing spline
c             with fp=s has been reached. probably causes : s too small
c             there is an approximation returned but the corresponding
c             weighted sum of squared residuals does not satisfy the
c             condition abs(fp-s)/s < tol.
c    ier=10 : error. on entry, the input data are controlled on validity
c             the following restrictions must be satisfied.
c             -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m
c             xb<=x(1)<x(2)<...<x(m)<=xe, lwrk>=(k+1)*m+nest*(7+3*k)
c             if iopt=-1: 2*k+2<=n<=min(nest,m+k+1)
c                         xb<t(k+2)<t(k+3)<...<t(n-k-1)<xe
c                       the schoenberg-whitney conditions, i.e. there
c                       must be a subset of data points xx(j) such that
c                         t(j) < xx(j) < t(j+k+1), j=1,2,...,n-k-1
c             if iopt>=0: s>=0
c                         if s=0 : nest >= m+k+1
c             if one of these conditions is found to be violated,control
c             is immediately repassed to the calling program. in that
c             case there is no approximation returned.
c
c  further comments:
c   by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the spline will be too smooth and signal will be
c   lost ; if s is too small the spline will pick up too much noise. in
c   the extreme cases the program will return an interpolating spline if
c   s=0 and the weighted least-squares polynomial of degree k if s is
c   very large. between these extremes, a properly chosen s will result
c   in a good compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the weights w(i). if these are
c   taken as 1/d(i) with d(i) an estimate of the standard deviation of
c   y(i), a good s-value should be found in the range (m-sqrt(2*m),m+
c   sqrt(2*m)). if nothing is known about the statistical error in y(i)
c   each w(i) can be set equal to one and s determined by trial and
c   error, taking account of the comments above. the best is then to
c   start with a very large value of s ( to determine the least-squares
c   polynomial and the corresponding upper bound fp0 for s) and then to
c   progressively decrease the value of s ( say by a factor 10 in the
c   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
c   approximation shows more detail) to obtain closer fits.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt=0.
c   if iopt=1 the program will continue with the set of knots found at
c   the last call of the routine. this will save a lot of computation
c   time if curfit is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   function underlying the data. but, if the computation mode iopt=1
c   is used, the knots returned may also depend on the s-values at
c   previous calls (if these were smaller). therefore, if after a number
c   of trials with different s-values and iopt=1, the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   curfit once more with the selected value for s but now with iopt=0.
c   indeed, curfit may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c
c  other subroutines required:
c    fpback,fpbspl,fpchec,fpcurf,fpdisc,fpgivs,fpknot,fprati,fprota
c
c  references:
c   dierckx p. : an algorithm for smoothing, differentiation and integ-
c                ration of experimental data using spline functions,
c                j.comp.appl.maths 1 (1975) 165-184.
c   dierckx p. : a fast algorithm for smoothing data on a rectangular
c                grid while using spline functions, siam j.numer.anal.
c                19 (1982) 1286-1304.
c   dierckx p. : an improved algorithm for curve fitting with spline
c                functions, report tw54, dept. computer science,k.u.
c                leuven, 1981.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : may 1979
c  latest update : march 1987
c
c  ..
c  ..scalar arguments..
      real xb,xe,s,fp
      integer iopt,m,k,nest,n,lwrk,ier
c  ..array arguments..
      real x(m),y(m),w(m),t(nest),c(nest),wrk(lwrk)
      integer iwrk(nest)
c  ..local scalars..
      real tol
      integer i,ia,ib,ifp,ig,iq,iz,j,k1,k2,lwest,maxit,nmin
c  ..
c  we set up the parameters tol and maxit
      maxit = 20
      tol = 0.1e-02
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(k.le.0 .or. k.gt.5) go to 50
      k1 = k+1
      k2 = k1+1
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 50
      nmin = 2*k1
      if(m.lt.k1 .or. nest.lt.nmin) go to 50
      lwest = m*k1+nest*(7+3*k)
      if(lwrk.lt.lwest) go to 50
      if(xb.gt.x(1) .or. xe.lt.x(m) .or. w(1).le.0.) go to 50
      do 10 i=2,m
         if(x(i-1).ge.x(i) .or. w(i).le.0.) go to 50
  10  continue
      if(iopt.ge.0) go to 30
      if(n.lt.nmin .or. n.gt.nest) go to 50
      j = n
      do 20 i=1,k1
         t(i) = xb
         t(j) = xe
         j = j-1
  20  continue
      call fpchec(x,m,t,n,k,ier)
      if(ier) 50,40,50
  30  if(s.lt.0.) go to 50
      if(s.eq.0. .and. nest.lt.(m+k1)) go to 50
      ier = 0
c we partition the working space and determine the spline approximation.
  40  ifp = 1
      iz = ifp+nest
      ia = iz+nest
      ib = ia+nest*k1
      ig = ib+nest*k2
      iq = ig+nest*k2
      call fpcurf(iopt,x,y,w,m,xb,xe,k,s,nest,tol,maxit,k1,k2,n,t,c,fp,
     * wrk(ifp),wrk(iz),wrk(ia),wrk(ib),wrk(ig),wrk(iq),iwrk,ier)
  50  return
      end
      subroutine fpback(a,z,n,k,c,nest)
c  subroutine fpback calculates the solution of the system of
c  equations a*c = z with a a n x n upper triangular matrix
c  of bandwidth k.
c  ..
c  ..scalar arguments..
      integer n,k,nest
c  ..array arguments..
      real a(nest,k),z(n),c(n)
c  ..local scalars..
      real store
      integer i,i1,j,k1,l,m
c  ..
      k1 = k-1
      c(n) = z(n)/a(n,1)
      i = n-1
      if(i.eq.0) go to 30
      do 20 j=2,n
        store = z(i)
        i1 = k1
        if(j.le.k1) i1 = j-1
        m = i
        do 10 l=1,i1
          m = m+1
          store = store-c(m)*a(i,l+1)
  10    continue
        c(i) = store/a(i,1)
        i = i-1
  20  continue
  30  return
      end
      subroutine fpbspl(t,n,k,x,l,h)
c  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
c  degree k at t(l) <= x < t(l+1) using the stable recurrence
c  relation of de boor and cox.
c  ..
c  ..scalar arguments..
      real x
      integer n,k,l
c  ..array arguments..
      real t(n),h(6)
c  ..local scalars..
      real f,one
      integer i,j,li,lj
c  ..local arrays..
      real hh(5)
c  ..
      one = 0.1e+01
      h(1) = one
      do 20 j=1,k
        do 10 i=1,j
          hh(i) = h(i)
  10    continue
        h(1) = 0.
        do 20 i=1,j
          li = l+i
          lj = li-j
          f = hh(i)/(t(li)-t(lj))
          h(i) = h(i)+f*(t(li)-x)
          h(i+1) = f*(x-t(lj))
  20  continue
      return
      end
      subroutine fpchec(x,m,t,n,k,ier)
c  subroutine fpchec verifies the number and the position of the knots
c  t(j),j=1,2,...,n of a spline of degree k, in relation to the number
c  and the position of the data points x(i),i=1,2,...,m. if all of the
c  following conditions are fulfilled, the error parameter ier is set
c  to zero. if one of the conditions is violated ier is set to ten.
c      1) k+1 <= n-k-1 <= m
c      2) t(1) <= t(2) <= ... <= t(k+1)
c         t(n-k) <= t(n-k+1) <= ... <= t(n)
c      3) t(k+1) < t(k+2) < ... < t(n-k)
c      4) t(k+1) <= x(i) <= t(n-k)
c      5) the conditions specified by schoenberg and whitney must hold
c         for at least one subset of data points, i.e. there must be a
c         subset of data points y(j) such that
c             t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1
c  ..
c  ..scalar arguments..
      integer m,n,k,ier
c  ..array arguments..
      real x(m),t(n)
c  ..local scalars..
      integer i,j,k1,k2,l,nk1,nk2,nk3
      real tj,tl
c  ..
      k1 = k+1
      k2 = k1+1
      nk1 = n-k1
      nk2 = nk1+1
      ier = 10
c  check condition no 1
      if(nk1.lt.k1 .or. nk1.gt.m) go to 80
c  check condition no 2
      j = n
      do 20 i=1,k
        if(t(i).gt.t(i+1)) go to 80
        if(t(j).lt.t(j-1)) go to 80
        j = j-1
  20  continue
c  check condition no 3
      do 30 i=k2,nk2
        if(t(i).le.t(i-1)) go to 80
  30  continue
c  check condition no 4
      if(x(1).lt.t(k1) .or. x(m).gt.t(nk2)) go to 80
c  check condition no 5
      if(x(1).ge.t(k2) .or. x(m).le.t(nk1)) go to 80
      i = 1
      l = k2
      nk3 = nk1-1
      if(nk3.lt.2) go to 70
      do 60 j=2,nk3
        tj = t(j)
        l = l+1
        tl = t(l)
  40    i = i+1
        if(i.ge.m) go to 80
        if(x(i).le.tj) go to 40
        if(x(i).ge.tl) go to 80
  60  continue
  70  ier = 0
  80  return
      end
      subroutine fpcurf(iopt,x,y,w,m,xb,xe,k,s,nest,tol,maxit,k1,k2,
     * n,t,c,fp,fpint,z,a,b,g,q,nrdata,ier)
c  ..
c  ..scalar arguments..
      real xb,xe,s,tol,fp
      integer iopt,m,k,nest,maxit,k1,k2,n,ier
c  ..array arguments..
      real x(m),y(m),w(m),t(nest),c(nest),fpint(nest),
     * z(nest),a(nest,k1),b(nest,k2),g(nest,k2),q(m,k1)
      integer nrdata(nest)
c  ..local scalars..
      real acc,con1,con4,con9,cos,half,fpart,fpms,fpold,fp0,f1,f2,f3,
     * one,p,pinv,piv,p1,p2,p3,rn,sin,store,term,wi,xi,yi
      integer i,ich1,ich3,it,iter,i1,i2,i3,j,k3,l,l0,
     * mk1,new,nk1,nmax,nmin,nplus,npl1,nrint,n8
c  ..local arrays..
      real h(7)
c  ..function references
      real abs,fprati
      integer max0,min0
c  ..subroutine references..
c    fpback,fpbspl,fpgivs,fpdisc,fpknot,fprota
c  ..
c  set constants
      one = 0.1e+01
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
      half = 0.5e0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 1: determination of the number of knots and their position     c
c  **************************************************************      c
c  given a set of knots we compute the least-squares spline sinf(x),   c
c  and the corresponding sum of squared residuals fp=f(p=inf).         c
c  if iopt=-1 sinf(x) is the requested approximation.                  c
c  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
c    if fp <=s we will continue with the current set of knots.         c
c    if fp > s we will increase the number of knots and compute the    c
c       corresponding least-squares spline until finally fp<=s.        c
c    the initial choice of knots depends on the value of s and iopt.   c
c    if s=0 we have spline interpolation; in that case the number of   c
c    knots equals nmax = m+k+1.                                        c
c    if s > 0 and                                                      c
c      iopt=0 we first compute the least-squares polynomial of         c
c      degree k; n = nmin = 2*k+2                                      c
c      iopt=1 we start with the set of knots found at the last         c
c      call of the routine, except for the case that s > fp0; then     c
c      we compute directly the least-squares polynomial of degree k.   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  determine nmin, the number of knots for polynomial approximation.
      nmin = 2*k1
      if(iopt.lt.0) go to 60
c  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
c  determine nmax, the number of knots for spline interpolation.
      nmax = m+k1
      if(s.gt.0.) go to 45
c  if s=0, s(x) is an interpolating spline.
c  test whether the required storage space exceeds the available one.
      n = nmax
      if(nmax.gt.nest) go to 420
c  find the position of the interior knots in case of interpolation.
  10  mk1 = m-k1
      if(mk1.eq.0) go to 60
      k3 = k/2
      i = k2
      j = k3+2
      if(k3*2.eq.k) go to 30
      do 20 l=1,mk1
        t(i) = x(j)
        i = i+1
        j = j+1
  20  continue
      go to 60
  30  do 40 l=1,mk1
        t(i) = (x(j)+x(j-1))*half
        i = i+1
        j = j+1
  40  continue
      go to 60
c  if s>0 our initial choice of knots depends on the value of iopt.
c  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
c  polynomial of degree k which is a spline without interior knots.
c  if iopt=1 and fp0>s we start computing the least squares spline
c  according to the set of knots found at the last call of the routine.
  45  if(iopt.eq.0) go to 50
      if(n.eq.nmin) go to 50
      fp0 = fpint(n)
      fpold = fpint(n-1)
      nplus = nrdata(n)
      if(fp0.gt.s) go to 60
  50  n = nmin
      fpold = 0.
      nplus = 0
      nrdata(1) = m-2
c  main loop for the different sets of knots. m is a save upper bound
c  for the number of trials.
  60  do 200 iter = 1,m
        if(n.eq.nmin) ier = -2
c  find nrint, tne number of knot intervals.
        nrint = n-nmin+1
c  find the position of the additional knots which are needed for
c  the b-spline representation of s(x).
        nk1 = n-k1
        i = n
        do 70 j=1,k1
          t(j) = xb
          t(i) = xe
          i = i-1
  70    continue
c  compute the b-spline coefficients of the least-squares spline
c  sinf(x). the observation matrix a is built up row by row and
c  reduced to upper triangular form by givens transformations.
c  at the same time fp=f(p=inf) is computed.
        fp = 0.
c  initialize the observation matrix a.
        do 80 i=1,nk1
          z(i) = 0.
          do 80 j=1,k1
            a(i,j) = 0.
  80    continue
        l = k1
        do 130 it=1,m
c  fetch the current data point x(it),y(it).
          xi = x(it)
          wi = w(it)
          yi = y(it)*wi
c  search for knot interval t(l) <= xi < t(l+1).
  85      if(xi.lt.t(l+1) .or. l.eq.nk1) go to 90
          l = l+1
          go to 85
c  evaluate the (k+1) non-zero b-splines at xi and store them in q.
  90      call fpbspl(t,n,k,xi,l,h)
          do 95 i=1,k1
            q(it,i) = h(i)
            h(i) = h(i)*wi
  95      continue
c  rotate the new row of the observation matrix into triangle.
          j = l-k1
          do 110 i=1,k1
            j = j+1
            piv = h(i)
            if(piv.eq.0.) go to 110
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(j,1),cos,sin)
c  transformations to right hand side.
            call fprota(cos,sin,yi,z(j))
            if(i.eq.k1) go to 120
            i2 = 1
            i3 = i+1
            do 100 i1 = i3,k1
              i2 = i2+1
c  transformations to left hand side.
              call fprota(cos,sin,h(i1),a(j,i2))
 100        continue
 110      continue
c  add contribution of this row to the sum of squares of residual
c  right hand sides.
 120      fp = fp+yi**2
 130    continue
        if(ier.eq.(-2)) fp0 = fp
        fpint(n) = fp0
        fpint(n-1) = fpold
        nrdata(n) = nplus
c  backward substitution to obtain the b-spline coefficients.
        call fpback(a,z,nk1,k1,c,nest)
c  test whether the approximation sinf(x) is an acceptable solution.
        if(iopt.lt.0) go to 440
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
c  if f(p=inf) < s accept the choice of knots.
        if(fpms.lt.0.) go to 250
c  if n = nmax, sinf(x) is an interpolating spline.
        if(n.eq.nmax) go to 430
c  increase the number of knots.
c  if n=nest we cannot increase the number of knots because of
c  the storage capacity limitation.
        if(n.eq.nest) go to 420
c  determine the number of knots nplus we are going to add.
        if(ier.eq.0) go to 140
        nplus = 1
        ier = 0
        go to 150
 140    npl1 = nplus*2
        rn = nplus
        if(fpold-fp.gt.acc) npl1 = rn*fpms/(fpold-fp)
        nplus = min0(nplus*2,max0(npl1,nplus/2,1))
 150    fpold = fp
c  compute the sum((w(i)*(y(i)-s(x(i))))**2) for each knot interval
c  t(j+k) <= x(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = 0.
        i = 1
        l = k2
        new = 0
        do 180 it=1,m
          if(x(it).lt.t(l) .or. l.gt.nk1) go to 160
          new = 1
          l = l+1
 160      term = 0.
          l0 = l-k2
          do 170 j=1,k1
            l0 = l0+1
            term = term+c(l0)*q(it,j)
 170      continue
          term = (w(it)*(term-y(it)))**2
          fpart = fpart+term
          if(new.eq.0) go to 180
          store = term*half
          fpint(i) = fpart-store
          i = i+1
          fpart = store
          new = 0
 180    continue
        fpint(nrint) = fpart
        do 190 l=1,nplus
c  add a new knot.
          call fpknot(x,m,t,n,fpint,nrdata,nrint,nest,1)
c  if n=nmax we locate the knots as for interpolation.
          if(n.eq.nmax) go to 10
c  test whether we cannot further increase the number of knots.
          if(n.eq.nest) go to 200
 190    continue
c  restart the computations with the new set of knots.
 200  continue
c  test whether the least-squares kth degree polynomial is a solution
c  of our approximation problem.
 250  if(ier.eq.(-2)) go to 440
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 2: determination of the smoothing spline sp(x).                c
c  ***************************************************                 c
c  we have determined the number of knots and their position.          c
c  we now compute the b-spline coefficients of the smoothing spline    c
c  sp(x). the observation matrix a is extended by the rows of matrix   c
c  b expressing that the kth derivative discontinuities of sp(x) at    c
c  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
c  ponding weights of these additional rows are set to 1/p.            c
c  iteratively we then have to determine the value of p such that      c
c  f(p)=sum((w(i)*(y(i)-sp(x(i))))**2) be = s. we already know that    c
c  the least-squares kth degree polynomial corresponds to p=0, and     c
c  that the least-squares spline corresponds to p=infinity. the        c
c  iteration process which is proposed here, makes use of rational     c
c  interpolation. since f(p) is a convex and strictly decreasing       c
c  function of p, it can be approximated by a rational function        c
c  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
c  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
c  to calculate the new value of p such that r(p)=s. convergence is    c
c  guaranteed by taking f1>0 and f3<0.                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  evaluate the discontinuity jump of the kth derivative of the
c  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
      call fpdisc(t,n,k2,b,nest)
c  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = 0.
      do 255 i=1,nk1
         p = p+a(i,1)
 255  continue
      rn = nk1
      p = rn/p
      ich1 = 0
      ich3 = 0
      n8 = n-nmin
c  iteration process to find the root of f(p) = s.
      do 360 iter=1,maxit
c  the rows of matrix b with weight 1/p are rotated into the
c  triangularised observation matrix a which is stored in g.
        pinv = one/p
        do 260 i=1,nk1
          c(i) = z(i)
          g(i,k2) = 0.
          do 260 j=1,k1
            g(i,j) = a(i,j)
 260    continue
        do 300 it=1,n8
c  the row of matrix b is rotated into triangle by givens transformation
          do 270 i=1,k2
            h(i) = b(it,i)*pinv
 270      continue
          yi = 0.
          do 290 j=it,nk1
            piv = h(1)
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,g(j,1),cos,sin)
c  transformations to right hand side.
            call fprota(cos,sin,yi,c(j))
            if(j.eq.nk1) go to 300
            i2 = k1
            if(j.gt.n8) i2 = nk1-j
            do 280 i=1,i2
c  transformations to left hand side.
              i1 = i+1
              call fprota(cos,sin,h(i1),g(j,i1))
              h(i) = h(i1)
 280        continue
            h(i2+1) = 0.
 290      continue
 300    continue
c  backward substitution to obtain the b-spline coefficients.
        call fpback(g,c,nk1,k2,c,nest)
c  computation of f(p).
        fp = 0.
        l = k2
        do 330 it=1,m
          if(x(it).lt.t(l) .or. l.gt.nk1) go to 310
          l = l+1
 310      l0 = l-k2
          term = 0.
          do 320 j=1,k1
            l0 = l0+1
            term = term+c(l0)*q(it,j)
 320      continue
          fp = fp+(w(it)*(term-y(it)))**2
 330    continue
c  test whether the approximation sp(x) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
c  test whether the maximal number of iterations is reached.
        if(iter.eq.maxit) go to 400
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 340
        if((f2-f3).gt.acc) go to 335
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p=p1*con9 + p2*con1
        go to 360
 335    if(f2.lt.0.) ich3=1
 340    if(ich1.ne.0) go to 350
        if((f1-f2).gt.acc) go to 345
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 360
        if(p.ge.p3) p = p2*con1 + p3*con9
        go to 360
 345    if(f2.gt.0.) ich1=1
c  test whether the iteration process proceeds as theoretically
c  expected.
 350    if(f2.ge.f1 .or. f2.le.f3) go to 410
c  find the new value for p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 360  continue
c  error codes and messages.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
      go to 440
 430  ier = -1
 440  return
      end
      subroutine fpdisc(t,n,k2,b,nest)
c  subroutine fpdisc calculates the discontinuity jumps of the kth
c  derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1)
c  ..scalar arguments..
      integer n,k2,nest
c  ..array arguments..
      real t(n),b(nest,k2)
c  ..local scalars..
      real an,fac,prod
      integer i,ik,j,jk,k,k1,l,lj,lk,lmk,lp,nk1,nrint
c  ..local array..
      real h(12)
c  ..
      k1 = k2-1
      k = k1-1
      nk1 = n-k1
      nrint = nk1-k
      an = nrint
      fac = an/(t(nk1+1)-t(k1))
      do 40 l=k2,nk1
        lmk = l-k1
        do 10 j=1,k1
          ik = j+k1
          lj = l+j
          lk = lj-k2
          h(j) = t(l)-t(lk)
          h(ik) = t(l)-t(lj)
  10    continue
        lp = lmk
        do 30 j=1,k2
          jk = j
          prod = h(j)
          do 20 i=1,k
            jk = jk+1
            prod = prod*h(jk)*fac
  20      continue
          lk = lp+k1
          b(lmk,j) = (t(lk)-t(lp))/prod
          lp = lp+1
  30    continue
  40  continue
      return
      end
      subroutine fpgivs(piv,ww,cos,sin)
c  subroutine fpgivs calculates the parameters of a givens
c  transformation .
c  ..
c  ..scalar arguments..
      real piv,ww,cos,sin
c  ..local scalars..
      real dd,one,store
c  ..function references..
      real abs,sqrt
c  ..
      one = 0.1e+01
      store = abs(piv)
      if(store.ge.ww) dd = store*sqrt(one+(ww/piv)**2)
      if(store.lt.ww) dd = ww*sqrt(one+(piv/ww)**2)
      cos = ww/dd
      sin = piv/dd
      ww = dd
      return
      end
      subroutine fpknot(x,m,t,n,fpint,nrdata,nrint,nest,istart)
c  subroutine fpknot locates an additional knot for a spline of degree
c  k and adjusts the corresponding parameters,i.e.
c    t     : the position of the knots.
c    n     : the number of knots.
c    nrint : the number of knotintervals.
c    fpint : the sum of squares of residual right hand sides
c            for each knot interval.
c    nrdata: the number of data points inside each knot interval.
c  istart indicates that the smallest data point at which the new knot
c  may be added is x(istart+1)
c  ..
c  ..scalar arguments..
      integer m,n,nrint,nest,istart
c  ..array arguments..
      real x(m),t(nest),fpint(nest)
      integer nrdata(nest)
c  ..local scalars..
      real an,am,fpmax
      integer ihalf,j,jbegin,jj,jk,jpoint,k,maxbeg,maxpt,
     * next,nrx,number
c  ..
      k = (n-nrint-1)/2
c  search for knot interval t(number+k) <= x <= t(number+k+1) where
c  fpint(number) is maximal on the condition that nrdata(number)
c  not equals zero.
      fpmax = 0.
      jbegin = istart
      do 20 j=1,nrint
        jpoint = nrdata(j)
        if(fpmax.ge.fpint(j) .or. jpoint.eq.0) go to 10
        fpmax = fpint(j)
        number = j
        maxpt = jpoint
        maxbeg = jbegin
  10    jbegin = jbegin+jpoint+1
  20  continue
c  let coincide the new knot t(number+k+1) with a data point x(nrx)
c  inside the old knot interval t(number+k) <= x <= t(number+k+1).
      ihalf = maxpt/2+1
      nrx = maxbeg+ihalf
      next = number+1
      if(next.gt.nrint) go to 40
c  adjust the different parameters.
      do 30 j=next,nrint
        jj = next+nrint-j
        fpint(jj+1) = fpint(jj)
        nrdata(jj+1) = nrdata(jj)
        jk = jj+k
        t(jk+1) = t(jk)
  30  continue
  40  nrdata(number) = ihalf-1
      nrdata(next) = maxpt-ihalf
      am = maxpt
      an = nrdata(number)
      fpint(number) = fpmax*an/am
      an = nrdata(next)
      fpint(next) = fpmax*an/am
      jk = next+k
      t(jk) = x(nrx)
      n = n+1
      nrint = nrint+1
      return
      end
      real function fprati(p1,f1,p2,f2,p3,f3)
c  given three points (p1,f1),(p2,f2) and (p3,f3), function fprati
c  gives the value of p such that the rational interpolating function
c  of the form r(p) = (u*p+v)/(p+w) equals zero at p.
c  ..
c  ..scalar arguments..
      real p1,f1,p2,f2,p3,f3
c  ..local scalars..
      real h1,h2,h3,p
c  ..
      if(p3.gt.0.) go to 10
c  value of p in case p3 = infinity.
      p = (p1*(f1-f3)*f2-p2*(f2-f3)*f1)/((f1-f2)*f3)
      go to 20
c  value of p in case p3 ^= infinity.
  10  h1 = f1*(f2-f3)
      h2 = f2*(f3-f1)
      h3 = f3*(f1-f2)
      p = -(p1*p2*h3+p2*p3*h1+p3*p1*h2)/(p1*h1+p2*h2+p3*h3)
c  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
  20  if(f2.lt.0.) go to 30
      p1 = p2
      f1 = f2
      go to 40
  30  p3 = p2
      f3 = f2
  40  fprati = p
      return
      end
      subroutine fprota(cos,sin,a,b)
c  subroutine fprota applies a givens rotation to a and b.
c  ..
c  ..scalar arguments..
      real cos,sin,a,b
c ..local scalars..
      real stor1,stor2
c  ..
      stor1 = a
      stor2 = b
      b = cos*stor2+sin*stor1
      a = cos*stor1-sin*stor2
      return
      end

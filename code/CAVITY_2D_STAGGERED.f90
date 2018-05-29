!==========================2D CAVITY FLOW=========================================
! PISO SCHEME FOR 2D (IN)COMPRESSIBLE FLOW
!A NEW PRESSURE-VELOCITY COUPLING PROCEDURE FOR INVICID AND VISCOUS FLOWS AT ALL SPEED
!
! U(1): PRE, U(2): XVEL; U(3): YVEL, U(4): TEMP
! STAGGERED GRID
! PRE,U,V,LOCATED DIFFERENT POSITION

MODULE PISO
parameter(MI=512,MJ=512)
REAL XL,DT,DTO,DX,DXO,DY,DYO,DZ,DZO,CFL
INTEGER NX,NXP,NXPP,NXM,NY,NYP,NYPP,NYM,NXH
REAL GAM,GAMO,GAMM,GAMMO
REAL DENI,UIN,VIN,WIN,PREI,RMACH,ENGI,RTI,REO
END MODULE
      PROGRAM SPHERE2D
	  USE PISO
	  REAL U(4,-2:MI+2,-2:MJ+2),UM(4,-2:MI+2,-2:MJ+2),UT(-2:MI+2,-2:MJ+2),UE(4,-2:MJ+2),UC(4,-2:MI+2,-2:MJ+2)
	  REAL US(4,-2:MI+2,-2:MJ+2),FB(4,-2:MI+2,-2:MJ+2),PE(-2:MI+2)
	  REAL DC(-2:MI+2,-2:MJ+2),D(-2:MI+2,-2:MJ+2),D0U(-2:MI+2,-2:MJ+2),D0V(-2:MI+2,-2:MJ+2)
	  REAL AP(-2:MI+2,-2:MJ+2),AM(-2:MI+2,-2:MJ+2)
	  REAL BP(-2:MI+2,-2:MJ+2),BM(-2:MI+2,-2:MJ+2)
	  REAL A1(-2:MI+2,-2:MJ+2),A2(-2:MI+2,-2:MJ+2)
	  REAL B1(-2:MI+2,-2:MJ+2),B2(-2:MI+2,-2:MJ+2)
	  REAL C(-2:MI+2,-2:MJ+2),PD(-2:MI+2,-2:MJ+2)
	  REAL X(-2:MI+2,-2:MJ+2),Y(-2:MI+2,-2:MJ+2),XC(-2:MI+2,-2:MJ+2),YC(-2:MI+2,-2:MJ+2)
	  PARAMETER(MAA=60,MPP=1)
	  REAL DW(-2:MI+2,-2:MJ+2)
	  REAL XE(MPP,MAA),YE(MPP,MAA),PRES(MAA)
	  REAL FCW(MPP,2),FCP(MPP,2),FCPI(MPP,MPP,2)
	  REAL SXC(MPP),SYC(MPP),SXCC(MPP),SYCC(MPP)
	  REAL PARU(MPP),PARV(MPP),PARUC(MPP),PARVC(MPP)
	  REAL FTATOLX(MPP),FTATOLY(MPP),FTATOLXC(MPP),FTATOLYC(MPP)
	  REAL DRAGX(MPP),DRAGX1(MPP),DRAGX2(MPP),DRAGY(MPP),DRAGY1(MPP),DRAGY2(MPP)
	  REAL WW(-2:MI+2,-2:MJ+2),WN(-2:MI+2,-2:MJ+2),PSI(-2:MI+2,-2:MJ+2)
      REAL:: RADIUS,GRAVND,GRAV,DENR
      REAL:: FAREA,FORCEX,FORCEY,FORCEX2,FORCEY2
	  REAL:: REON,REOP,REN,REP,DIAML,DIAMP,DIAMN,DIAMNO
      PI=4.*ATAN(1.)
      MIP=MI+1; MIH=MI/2 ; MIM=MI-1
	  MJP=MJ+1; MJH=MJ/2 ; MJM=MJ-1
	  XL=1.; DX=XL/MI; DXO=1./DX; DXO2=DXO*DXO
	  YL=1.; DY=YL/MJ; DYO=1./DY; DYO2=DYO*DYO
	  CFL=0.1 !; DT=1. ; DTO=1./DT !(CHANGE)
	  ITAT=20000000; TATOL=6000. !(ITAT=5000 OVER-----------SEE ERR)
	  IOUT=10
	  IREAD=0
      DIAMETER=10. ; RADIUS=0.5*DIAMETER
      DENI=1.; UIN=1. ; VIN=0.; PREI=1.0 !(RE=UIN/REO)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RE=100. !(RE=100, 400, 1000)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REO=UIN*XL/RE
      PR=0.72; PRO=1./PR
	  VCOX=REO*DXO*DXO; VCOY=REO*DYO*DYO
	  R1=DX/UIN
	  CDISP=0. ; CDISPP=1.1
      DT=CFL*DX/UIN
	  DTO=1/DT
	  ESLON=0.0000001; AV=0. ! AV: ARTIFICIAL VISCIUS IN FIRST VELOCITY PREDICATOR
	  DO I=0,MIP; DO J=0,MJP
	  X(I,J)=(I-0.5)*DX
	  Y(I,J)=(J-0.5)*DY
	  ENDDO; ENDDO
!========================================================================================
	  DO I=-2,MIP+1; DO J=-2,MJP+1
	  U(1,I,J)=PREI
	  U(2,I,J)=0.
	  U(3,I,J)=0.
	  U(4,I,J)=1.
	  END DO ; END DO
	  IF (IREAD.EQ.1) THEN
	  open (7,FILE='SOL0.DAT',STATUS='unknown')
	   DO I=-2,MIP+1 ; DO J=-2,MJP+1
			DO L=1,4
	        READ(7,*)U(L,I,J)
	        ENDDO
	   ENDDO; ENDDO
	   CLOSE(7)
	  END IF
	  CALL BOUNDARY(U)
	  ITER=0; TAT=0.
1000  ITER=ITER+1; TAT=TAT+DT
      UMAX=0.001; VMAX=0.001; WMAX=0.001; EMAX=0.001; DMIN=100; PMIN=100.
	  TEMPC=1.
      DO I=1,MI; DO J=1,MJ
	  TEMPU=ABS(U(2,I,J)); TEMPV=ABS(U(3,I,J))
	  EMAX=MAX(EMAX, TEMPU+TEMPC, TEMPV+TEMPC)
	  ENDDO; ENDDO
	  DT=CFL*DX/EMAX
	  DTO=1./DT
	  DO L=1,4
	  DO I=-2,MIP+1; DO J=-2,MJP+1
	  UM(L,I,J)=U(L,I,J)
	  ENDDO; ENDDO
	  ENDDO
! PREDICTOR: VELOCITY
! ----U------------------------------------------------------
	  DO J=1,MJ ; DO I=1,MIM
	  IP=I+1; IM=I-1; JP=J+1; JM=J-1
	  DUP=0.5*(U(2,IP,J)+U(2,I,J))
	  DUM=0.5*(U(2,IM,J)+U(2,I,J))
	  DUPP=0.5*(DUP+ABS(DUP)); DUPM=DUP-DUPP
	  DUMP=0.5*(DUM+ABS(DUM)); DUMM=DUM-DUMP
	  DC(I,J)=DXO*(DUPP-DUMM)+2.*VCOX
	  AP(I,J)=DXO*DUPM-VCOX
	  AM(I,J)=-DXO*DUMP-VCOX
!
	  DUP=0.5*(U(3,I,J)+U(3,IP,J))
	  DUM=0.5*(U(3,I,JM)+U(3,IP,JM))
	  DUPP=0.5*(DUP+ABS(DUP)+AV); DUPM=DUP-DUPP
	  DUMP=0.5*(DUM+ABS(DUM)+AV); DUMM=DUM-DUMP
	  DC(I,J)=DC(I,J)+DYO*(DUPP-DUMM)+2.*VCOY
	  BP(I,J)=DYO*DUPM-VCOY
	  BM(I,J)=-DYO*DUMP-VCOY
	  END DO ; END DO
!CUT
	  DO J=1,MJ;DO I=1,MIM
	  D0U(I,J)=DTO+DC(I,J) ; D(I,J)=D0U(I,J)
	  A1(I,J)=AM(I,J); A2(I,J)=AP(I,J)
	  B1(I,J)=BM(I,J); B2(I,J)=BP(I,J)
	  C(I,J)=-DXO*(U(1,I+1,J)-U(1,I,J))+DTO*U(2,I,J)
	  ENDDO; ENDDO
!	  DO J=1,MJ
!	  D(1,J)=D(1,J)-A1(1,J); A1(1,J)=0.
!	  D(MIM,J)=D(MIM,J)-A2(MIM,J); A2(MIM,J)=0.
!	  ENDDO
	  DO I=1,MIM
	  D(I,1)=D(I,1)-B1(I,1); B1(I,1)=0.
!	  D(I,MJ)=D(I,MJ)-B2(I,MJ); B2(I,MJ)=0.
	  ENDDO
      UT=0.
      DO I=0,MIP; DO J=0,MJP
	  UT(I,J)=U(2,I,J)
	  ENDDO; ENDDO
	  CALL SYUU(D,A2,A1,B2,B1,C,UT)
	  DO I=1,MI; DO J=1,MJ
	  US(2,I,J)=UT(I,J)
	  ENDDO; ENDDO
! ----V--------------------------------
      DO J=1,MJM ; DO I=1,MI
	  IP=I+1; IM=I-1; JP=J+1; JM=J-1
	  DUP=0.5*(U(2,I,J)+U(2,I,JP))
	  DUM=0.5*(U(2,IM,J)+U(2,IM,JP))
	  DUPP=0.5*(DUP+ABS(DUP)); DUPM=DUP-DUPP
	  DUMP=0.5*(DUM+ABS(DUM)); DUMM=DUM-DUMP
	  DC(I,J)=DXO*(DUPP-DUMM)+2.*VCOX
	  AP(I,J)=DXO*DUPM-VCOX
	  AM(I,J)=-DXO*DUMP-VCOX
!
	  DUP=0.5*(U(3,I,J)+U(3,I,JP))
	  DUM=0.5*(U(3,I,J)+U(3,I,JM))
	  DUPP=0.5*(DUP+ABS(DUP)+AV); DUPM=DUP-DUPP
	  DUMP=0.5*(DUM+ABS(DUM)+AV); DUMM=DUM-DUMP
	  DC(I,J)=DC(I,J)+DYO*(DUPP-DUMM)+2.*VCOY
	  BP(I,J)=DYO*DUPM-VCOY
	  BM(I,J)=-DYO*DUMP-VCOY
	  END DO ; END DO
!CUT
	  DO I=1,MI ; DO J=1,MJM
	  D0V(I,J)=DTO+DC(I,J) ; D(I,J)=D0V(I,J)
	  A1(I,J)=AM(I,J); A2(I,J)=AP(I,J)
	  B1(I,J)=BM(I,J); B2(I,J)=BP(I,J)
	  C(I,J)=-DYO*(U(1,I,J+1)-U(1,I,J))+DTO*U(3,I,J)
	  ENDDO; ENDDO
	  DO J=1,MJM
	  D(1,J)=D(1,J)-A1(1,J); A1(1,J)=0.
	  D(MI,J)=D(MI,J)-A2(MI,J); A2(MI,J)=0.
	  ENDDO
!	  DO I=1,MI
!	  D(I,1)=D(I,1)-B1(I,1); B1(I,1)=0.
!	  D(I,MJ)=D(I,MJ)-B2(I,MJ); B1(I,MJ)=0.
!	  ENDDO
      UT=0.
      DO I=0,MIP; DO J=0,MJP
	  UT(I,J)=U(3,I,J)
	  ENDDO; ENDDO
	  CALL SYUV(D,A2,A1,B2,B1,C,UT)
	  DO I=1,MI; DO J=1,MJ
	  US(3,I,J)=UT(I,J)
	  ENDDO; ENDDO
 CALL BOUNDARY(US)
! CORRECTOR: PRESSURE, VELOCITY
      DXO2=DXO*DXO; DYO2=DYO*DYO
      DO I=1,MI; DO J=1,MJ
      IM=I-1 ; JM=J-1
      D(I,J)=-2.*DXO2
	  A2(I,J)=DXO2
	  A1(I,J)=DXO2
!
	  D(I,J)=D(I,J)-2.*DYO2
	  B2(I,J)=DYO2
	  B1(I,J)=DYO2
!
      D(I,J)=D(I,J)*CDISPP
      TEMP=DXO*(US(2,I,J)-US(2,IM,J))
	  TEMP=TEMP+DYO*(US(3,I,J)-US(3,I,JM))
	  C(I,J)=TEMP
	  ENDDO; ENDDO
!MATRIX SOLVER B.C.
	  DO J=1,MJ
	  D(1,J)=D(1,J)+A1(1,J); A1(1,J)=0.
	  D(MI,J)=D(MI,J)+A2(MI,J); A2(MI,J)=0.
	  ENDDO
	  DO I=1,MI
	  D(I,1)=D(I,1)+B1(I,1); B1(I,1)=0.
	  D(I,MJ)=D(I,MJ)+B2(I,MJ); B2(I,MJ)=0.
	  ENDDO
	  DO I=0,MIP; DO J=0,MJP
	  UT(I,J)=0.
	  ENDDO; ENDDO
	  CALL SYUP(D,A2,A1,B2,B1,C,UT)
	  DO J=1,MJ
	  UT(0,J)=UT(1,J)
	  UT(MIP,J)=UT(MI,J)
	  ENDDO
	  DO I=1,MI
	  UT(I,0)=UT(I,1)
	  UT(I,MJP)=UT(I,MJ)
	  ENDDO
	  DO I=1,MI; DO J=1,MJ
	  IP=I+1; IM=I-1; JP=J+1; JM=J-1
	  U(1,I,J)=U(1,I,J)+UT(I,J)*DTO
	  U(2,I,J)=US(2,I,J)-DXO*(UT(IP,J)-UT(I,J))
	  U(3,I,J)=US(3,I,J)-DYO*(UT(I,JP)-UT(I,J))
	  ENDDO; ENDDO
	 CALL BOUNDARY(U)
!================LW  LENGTH=================================
	  ERRU=0.; ERRD=0.
	  DO I=1,MI; DO J=1,MJ
	  UTEMP=U(2,I,J)-UM(2,I,J)
	  ERRU=ERRU+UTEMP*UTEMP
	  DTEMP=U(1,I,J)-UM(1,I,J)
	  ERRD=ERRD+DTEMP*DTEMP
	  ENDDO; ENDDO
	  ERRU=SQRT(ERRU/(MI*MJ))
	  ERRD=SQRT(ERRD/(MI*MJ))
!	  WRITE(*,*)ITER,ERRU,ERRD
	  COMP=0.; COMPMAX=0.
	  DO I=2,MI-1; DO J=2,MJ-1
	  COMPIJ=0.5*DXO*(U(2,I,J)-U(2,I-1,J))+0.5*DYO*(U(3,I,J)-U(3,I,J-1))
	  COMP=COMP+COMPIJ*COMPIJ
	  COMPMAX=AMAX1(COMPMAX,COMPIJ)
	  ENDDO; ENDDO
	  ERRC=SQRT(COMP/((MI-2)*(MJ-2)))
WRITE(*,*)ITER,DT,ERRD,ERRU,ERRC
!=================================================================================
IF(ITER/1000*1000==ITER) THEN
    open (7,FILE='DEN.DAT',STATUS='unknown')
       write(7,*) 'VARIABLES = X,Y,PRE,U,V,WN,PSI'
       write(7,*) 'ZONE T="ZONE1",I=',MIP,'J=',MJP
	   DO J=0,MJ ; DO I=0,MI
	   UC(2,I,J)=0.5*(U(2,I,J)+U(2,I+1,J))
	   UC(3,I,J)=0.5*(U(3,I,J)+U(3,I,J+1))
	   WN(I,J)=0.5*(UC(3,I+1,J)-UC(3,I-1,J)-UC(2,I,J+1)+UC(2,I,J-1))
	   TEMP1=TEMP1+UC(3,I,0)
	   TEMP2=TEMP2+UC(2,I,J)
	   PSI(I,J)=TEMP2-TEMP1
	   WRITE(7,*)I*DX,J*DY,U(1,I,J),UC(2,I,J),UC(3,I,J),WN(I,J),PSI(I,J)
	   ENDDO; ENDDO
	   CLOSE(7)
	   open (7,FILE='UY.DAT',STATUS='unknown')
	   DO J=0,MJ
	   I=0.5*MI
	   WRITE(7,*)UC(2,I,J)/UIN,Y(I,J)
	   ENDDO
	   CLOSE(7)
	   open (7,FILE='VX.DAT',STATUS='unknown')
	   DO I=0,MI
	   J=0.5*MJ
	   WRITE(7,*)X(I,J),UC(3,I,J)/UIN
	   ENDDO
	   CLOSE(7)
END IF
!==============SURFACE PRESSURE
!	   IF(ERRU>0.000001) GOTO 1000
       IF(ITER<ITAT) GOTO 1000
	   open (7,FILE='SOL.DAT',STATUS='unknown')
	   DO I=-2,MIP+1; DO J=-2,MJP+1
			DO L=1,4
	        WRITE(7,*)U(L,I,J)
	        ENDDO
	   ENDDO; ENDDO
	   CLOSE(7)
100	  STOP
	  END
!-----------------------------------------------------
      SUBROUTINE BOUNDARY(W)
	  USE PISO
	  REAL W(4,-2:MI+2,-2:MJ+2)
	  MIP=MI+1; MIPP=MI+2; MJP=MJ+1; MJPP=MJ+2
      DO J=1,MJ
	  W(1,0,J)=W(1,1,J)
	  W(2,0,J)=0.
	  W(3,0,J)=-W(3,1,J)
	  ENDDO
!I=MIP
	  DO J=1,MJ
	  W(1,MIP,J)=W(1,MI,J)
	  W(2,MI,J)=0.
	  W(3,MIP,J)=-W(3,MI,J)
	  ENDDO
!J=MJP
      DO I=0,MIP
	  W(1,I,MJP)=W(1,I,MJ)
	  W(2,I,MJP)=UIN
	  W(3,I,MJ)=VIN
	  ENDDO
!J=0
	  DO I=0,MIP
	  W(1,I,0)=W(1,I,1)
	  W(2,I,0)=-W(2,I,1)
	  W(3,I,0)=0.
	  ENDDO
!I=0
	  RETURN
	  END
!-----------------------------------------------------
      SUBROUTINE SYUU(DCT,APT,AMT,BPT,BMT,CC,UT)
	  USE PISO
	  REAL DCT(-2:MI+2,-2:MJ+2),APT(-2:MI+2,-2:MJ+2),AMT(-2:MI+2,-2:MJ+2)
	  REAL BPT(-2:MI+2,-2:MJ+2),BMT(-2:MI+2,-2:MJ+2)
	  REAL CC(-2:MI+2,-2:MJ+2),UT(-2:MI+2,-2:MJ+2),UTO(-2:MI+2,-2:MJ+2)
      MIM=MI-1; MJM=MJ-1
	  ITT=0
1000  ITT=ITT+1
      DO I=1,MIM; DO J=1,MJ
	  UTO(I,J)=UT(I,J)
	  ENDDO; ENDDO
      DO J=1,MJ; DO I=1,MIM
	  IP=I+1; IM=I-1; JP=J+1; JM=J-1
	  TEMP=CC(I,J)-APT(I,J)*UT(IP,J)-AMT(I,J)*UT(IM,J)&
	               &-BPT(I,J)*UT(I,JP)-BMT(I,J)*UT(I,JM)
      UT(I,J)=TEMP/DCT(I,J)
	  ENDDO; ENDDO
	  ERR=0.
	  DO J=1,MJ; DO I=1,MIM
	  UTEMP=UT(I,J)-UTO(I,J)
	  ERR=ERR+UTEMP*UTEMP
	  ENDDO; ENDDO
	  ERR=SQRT(ERR/(MIM*MJ))
	  IF(ERR .GT. 0.000001 .OR. ITT .LT. 10) GOTO 1000
      RETURN
	  END
!=================================================================
SUBROUTINE SYUV(DCT,APT,AMT,BPT,BMT,CC,UT)
	  USE PISO
	  REAL DCT(-2:MI+2,-2:MJ+2),APT(-2:MI+2,-2:MJ+2),AMT(-2:MI+2,-2:MJ+2)
	  REAL BPT(-2:MI+2,-2:MJ+2),BMT(-2:MI+2,-2:MJ+2)
	  REAL CC(-2:MI+2,-2:MJ+2),UT(-2:MI+2,-2:MJ+2),UTO(-2:MI+2,-2:MJ+2)
      MIM=MI-1; MJM=MJ-1
	  ITT=0
1000  ITT=ITT+1
      DO I=1,MI; DO J=1,MJM
	  UTO(I,J)=UT(I,J)
	  ENDDO; ENDDO
      DO J=1,MJM; DO I=1,MI
	  IP=I+1; IM=I-1; JP=J+1; JM=J-1
	  TEMP=CC(I,J)-APT(I,J)*UT(IP,J)-AMT(I,J)*UT(IM,J)&
	               &-BPT(I,J)*UT(I,JP)-BMT(I,J)*UT(I,JM)
      UT(I,J)=TEMP/DCT(I,J)
	  ENDDO; ENDDO
	  ERR=0.
	  DO J=1,MJM; DO I=1,MI
	  UTEMP=UT(I,J)-UTO(I,J)
	  ERR=ERR+UTEMP*UTEMP
	  ENDDO; ENDDO
	  ERR=SQRT(ERR/(MI*MJM))
	  IF(ERR .GT. 0.000001 .OR. ITT .LT. 10) GOTO 1000
      RETURN
	  END
!------------------------------------------------------------
      SUBROUTINE SYUP(DCT,APT,AMT,BPT,BMT,CC,UT)
	  USE PISO
	  REAL DCT(-2:MI+2,-2:MJ+2),APT(-2:MI+2,-2:MJ+2),AMT(-2:MI+2,-2:MJ+2)
	  REAL BPT(-2:MI+2,-2:MJ+2),BMT(-2:MI+2,-2:MJ+2)
	  REAL CC(-2:MI+2,-2:MJ+2),UT(-2:MI+2,-2:MJ+2),UTO(-2:MI+2,-2:MJ+2)
	  ITT=0; IGT=MI*MJ
1000  ITT=ITT+1
      DO I=1,MI; DO J=1,MJ
	  UTO(I,J)=UT(I,J)
	  ENDDO; ENDDO
      DO I=1,MI; DO J=1,MJ
	  IP=I+1; IM=I-1; JP=J+1; JM=J-1
	  TEMP=CC(I,J)-APT(I,J)*UT(IP,J)-AMT(I,J)*UT(IM,J)&
	               &-BPT(I,J)*UT(I,JP)-BMT(I,J)*UT(I,JM)
      UT(I,J)=TEMP/DCT(I,J)
	  ENDDO; ENDDO
	  ERR=0.
	  DO J=1,MJ ; DO I=1,MI
	  UTEMP=UT(I,J)-UTO(I,J)
	  ERR=ERR+UTEMP*UTEMP
	  ENDDO; ENDDO
	  ERR=SQRT(ERR/(MI*MJ))
	  IF(ERR .GT. 0.000001 .OR. ITT .LT. 10) GOTO 1000
!	  WRITE(*,*)'GS_ITER',ITT,ERR
      RETURN
	  END
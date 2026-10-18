c-----------------------------------------------------------------------
c  Authored by IY
c  China University of Petroleum (East China), China
c-----------------------------------------------------------------------
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
	 

      INCLUDE 'ABA_PARAM.INC'
	  
!     WARNING - the aba_param.inc file declares
!        Implicit real*8(a-h,o-z)
!     This means that, by default, any variables with
!     first letter between a-h or o-z are double precision.
!     The rest are integers.
!     Note that this also means that if you type a variable
!     name incorrectly, the compiler won't catch your typo.	  

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
      INTEGER i,j
      DOUBLE PRECISION UPSTRAN(NTENS)
      DOUBLE PRECISION FI(4), DMG(4)
      DOUBLE PRECISION dft, dfc, dmt, dmc, smt, smc
      DOUBLE PRECISION df, dm, dft_eff, dfc_eff, dmt_eff, dmc_eff
      DOUBLE PRECISION rf, rm, rmt, rmc
c

      PARAMETER (ZERO = 0.0D0,ONE=1.0D0, TWO=2.0D0)
c---elastic constants
      E11 = PROPS(1)    ! YOUNGS MODULUS IN 1 DIRECTION
      E22 = PROPS(2)    ! YOUNGS MODULUS IN 2 DIRECTION
      E33 = PROPS(3)    ! YOUNGS MODULUS IN 3 DIRECTION
      PR12 = PROPS(4)   ! MAJOR POISSONS RATIO IN 1-2 PLANE
      PR13 = PROPS(5)   ! MAJOR POISSONS RATIO IN 1-3 PLANE
      PR23 = PROPS(6)   ! MAJOR POISSONS RATIO IN 2-3 PLANE
      G12 = PROPS(7)    ! SHEAR MODULUS IN 1-2 PLANE
      G13 = PROPS(8)    ! SHEAR MODULUS IN 1-3 PLANE
      G23 = PROPS(9)    ! SHEAR MODULUS IN 2-3 PLANE
c---strengths parameters    
      XT = PROPS(10)
      XC = PROPS(11)
      YT = PROPS(12)
      YC = PROPS(13)
      S12 = PROPS(14)
      S13 = PROPS(15)
      S23 = PROPS(16)
c---degradation parameters
      dft=PROPS(17)
      dfc=PROPS(18)
      dmt=PROPS(19) 
      dmc=PROPS(20)
      smt=PROPS(21)     ! recommended value is 0.9
      smc=PROPS(22)     ! recommended value is 0.5
c---
      PR21 = (E22*PR12)/E11  ! CALCULAION OF MINOR POISSONS RATIO
      PR31 = (E33*PR13)/E11
      PR32 = (E33*PR23)/E22
      DELTA = (ONE)/(ONE-PR12*PR21-PR23*PR32-PR13*PR31-TWO*PR21*PR32*PR13)
c----read damage variables from STATEV so damage history persists
c----across analysis steps; Abaqus initializes STATEV to zero at the
c----start of the analysis if no prior history is provided.
      df=MIN(MAX(STATEV(5),ZERO),ONE)  ! fiber damage variable
      dm=MIN(MAX(STATEV(6),ZERO),ONE)  ! matrix damage variable
      DMG(1)=MIN(MAX(STATEV(7),ZERO),ONE)
      DMG(2)=MIN(MAX(STATEV(8),ZERO),ONE)
      DMG(3)=MIN(MAX(STATEV(9),ZERO),ONE)
      DMG(4)=MIN(MAX(STATEV(10),ZERO),ONE)
c----initial damage variables
      dft_eff = DMG(1)*dft
      dfc_eff = DMG(2)*dfc
      dmt_eff = DMG(3)*dmt
      dmc_eff = DMG(4)*dmc
c
      df = 1D0 - (1D0-dft_eff)*(1D0-dfc_eff)
      dm = 1D0 - (1D0-dmt_eff)*(1D0-dmc_eff)
      rf = MAX(ONE-df, 1D-6)
      rm = MAX(ONE-dm, 1D-6)
      rmt = MAX(ONE-smt*dmt_eff, 1D-6)
      rmc = MAX(ONE-smc*dmc_eff, 1D-6)
c
c
c----initialization of the Jacobian matrix
      DO I=1,NTENS
          DO J=1,NTENS
              DDSDDE(I,J) = ZERO
          END DO
      END DO
C        
      DDSDDE(1,1) = (E11*DELTA*(ONE-PR23*PR32))*rf
      DDSDDE(1,2) = (E11*DELTA*(PR21+PR31*PR23))*rf*rm
      DDSDDE(1,3) = (E11*DELTA*(PR31+PR21*PR32))*rf*rm
      DDSDDE(2,1) = (E11*DELTA*(PR21+PR31*PR23))*rf*rm
      DDSDDE(2,2) = (E22*DELTA*(ONE-PR13*PR31))*rf*rm
      DDSDDE(2,3) = (E22*DELTA*(PR32+PR12*PR31))*rf*rm
      DDSDDE(3,1) = (E11*DELTA*(PR31+PR21*PR32))*rf*rm
      DDSDDE(3,2) = (E22*DELTA*(PR32+PR12*PR31))*rf*rm
      DDSDDE(3,3) = (E33*DELTA*(ONE-PR12*PR21))*rf*rm
      DDSDDE(4,4) = G12*rf*rmt*rmc
      DDSDDE(5,5) = G13*rf*rmt*rmc
      DDSDDE(6,6) = G23*rf*rmt*rmc
C    
      DO j=1,NTENS
        upstran(j) = STRAN(j) + DSTRAN(j)
      ENDDO
      DO i=1,NTENS
      	STRESS(i)=0.0D0
      	DO j=1,NTENS
      		STRESS(i)=STRESS(i)+DDSDDE(i,j)*UPSTRAN(j)
      	ENDDO
      ENDDO
c
c-----
      sig11=STRESS(1) 
      sig22=STRESS(2) 
      sig33=STRESS(3) 
      tau12=STRESS(4) 
      tau13=STRESS(5) 
      tau23=STRESS(6)
c-----
      FI(1:4)=0D0
      if(sig11.GE. 0.0) THEN
          FI(1) = (sig11/XT)**2 + (tau12/S12)**2 + (tau13/S13)**2
      ELSE
          FI(2) = (-sig11/XC)**2
      END IF
      if((sig22 + sig33) .GE. 0.0) THEN
          FI(3) = ((sig22 + sig33)/YT)**2 + (tau23**2-sig22*sig33)/(S23**2)+(tau12/S12)**2 + (tau13/S13)**2
      ELSE
          FI(4) = (ONE/YC)*((YC/(TWO*S23))**2-ONE)*(sig22+sig33)
     1      + ((sig22+sig33)**2)/((TWO*S23)**2)
     2      + (tau23**2-sig22*sig33)/(S23**2)
     3      + (tau12**2)/(S12**2)
     4      + (tau13**2)/(S13**2)
      END IF
c-----
      IF (FI(1) .GE. 1D0 .AND. DMG(1) .LT. ONE) THEN
         DMG(1)=1D0
      END IF

      IF (FI(2) .GE. 1D0 .AND. DMG(2) .LT. ONE) THEN
         DMG(2)=1D0
      END IF

      IF (FI(3) .GE. 1D0 .AND. DMG(3) .LT. ONE) THEN
         DMG(3)=1D0
      END IF

      IF (FI(4) .GE. 1D0 .AND. DMG(4) .LT. ONE) THEN
         DMG(4)=1D0
      END IF
c---update damage variables
      dft_eff = DMG(1)*dft
      dfc_eff = DMG(2)*dfc
      dmt_eff = DMG(3)*dmt
      dmc_eff = DMG(4)*dmc
      df = 1D0 - (1D0-dft_eff)*(1D0-dfc_eff)
      dm = 1D0 - (1D0-dmt_eff)*(1D0-dmc_eff)
c---keep the constitutive response from the trial state so that
c---new damage only affects the next increment, matching the older,
c---numerically smoother UMAT behavior.
c
      STATEV(1)=FI(1)
      STATEV(2)=FI(2)
      STATEV(3)=FI(3)
      STATEV(4)=FI(4)
      STATEV(5)=df
      STATEV(6)=dm
      STATEV(7)=DMG(1)
      STATEV(8)=DMG(2)
      STATEV(9)=DMG(3)
      STATEV(10)=DMG(4)
      STATEV(11)=ONE
c---
c-----
      RETURN
      END

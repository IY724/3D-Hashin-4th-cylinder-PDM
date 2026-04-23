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
      DOUBLE PRECISION FI(4),DELTAEQ(4),SIGMAEQ(4)
      DOUBLE PRECISION DELTA0(4),SIGMA0(4),UPSTRAN(NTENS)
      DOUBLE PRECISION EFFSTR(NTENS),DDSDDE0(NTENS,NTENS)
      DOUBLE PRECISION deltaf_ft, deltaf_fc, deltaf_mt, deltaf_mc
      DOUBLE PRECISION dft, dfc, dmt, dmc, df, dm, ds
      DOUBLE PRECISION dftv, dfcv, dmtv, dmcv, dfv, dmv, dsv
      DOUBLE PRECISION kappa_ft, kappa_fc, kappa_mt, kappa_mc
      DOUBLE PRECISION init_ft, init_fc, init_mt, init_mc
      DOUBLE PRECISION eps11, eps22, eps33, gam12, gam13, gam23
      DOUBLE PRECISION sig11, sig22, sig33, tau12, tau13, tau23
      DOUBLE PRECISION alpha, smt, smc, vf, vm, beta_old, beta_new
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
c---damage parameters
      G_ft=PROPS(17)
      G_fc=PROPS(18)
      G_mt=PROPS(19) 
      G_mc=PROPS(20)
c----Damage Control
      alpha=PROPS(21) ! reserved model parameter
      smt=PROPS(22)   ! recommended value is 0.9
      smc=PROPS(23)   ! recommended value is 0.5
c----viscosity coefficients for damage regularization
c----PROPS(24)=fiber viscosity, PROPS(25)=matrix viscosity
c----if only one viscosity is provided, it is used for both modes
      vf=ZERO
      vm=ZERO
      IF (NPROPS .GE. 25) THEN
          vf=MAX(PROPS(24),ZERO)
          vm=MAX(PROPS(25),ZERO)
      ELSE IF (NPROPS .GE. 24) THEN
          vf=MAX(PROPS(24),ZERO)
          vm=vf
      END IF
c----minor Poisson ratios and orthotropic scaling factor
      PR21 = (E22*PR12)/E11
      PR31 = (E33*PR13)/E11
      PR32 = (E33*PR23)/E22
      DELTA = (ONE)/(ONE-PR12*PR21-PR23*PR32-PR13*PR31
     1     -TWO*PR21*PR32*PR13)
c----define statev variables
c----damage variables
      dft=STATEV(1)  ! tensile damage variable
      dfc=STATEV(2)  ! compressive damage variable
      dmt=STATEV(3)  ! tensile damage variable in shear
      dmc=STATEV(4)  ! compressive damage variable in shear
c----equivalent strains
      kappa_ft=STATEV(5)  ! equivalent strain for tensile damage
      kappa_fc=STATEV(6)  ! equivalent strain for compressive damage
      kappa_mt=STATEV(7)  ! equivalent strain for tensile damage in shear
      kappa_mc=STATEV(8)  ! equivalent strain for compressive damage in shear
c----damage initiation flags
      init_ft=STATEV(9)   ! fiber tension initiation flag
      init_fc=STATEV(10)  ! fiber compression initiation flag
      init_mt=STATEV(11)  ! matrix tension initiation flag
      init_mc=STATEV(12)  ! matrix compression initiation flag
c----damage onset equivalent displacement
      DELTA0(1)=STATEV(13)    ! fiber tension onset displacement
      DELTA0(2)=STATEV(14)    ! fiber compression onset displacement
      DELTA0(3)=STATEV(15)    ! matrix tension onset displacement
      DELTA0(4)=STATEV(16)    ! matrix compression onset displacement
c----damage onset equivalent stress
      SIGMA0(1)=STATEV(17)    ! fiber tension onset stress
      SIGMA0(2)=STATEV(18)    ! fiber compression onset stress
      SIGMA0(3)=STATEV(19)    ! matrix tension onset stress
      SIGMA0(4)=STATEV(20)    ! matrix compression onset stress
c----viscous damage variables used in stiffness degradation
      dftv=STATEV(21)
      dfcv=STATEV(22)
      dmtv=STATEV(23)
      dmcv=STATEV(24)
c----current equivalent displacement/stress placeholders
      DELTAEQ(1)=ZERO
      DELTAEQ(2)=ZERO
      DELTAEQ(3)=ZERO
      DELTAEQ(4)=ZERO
      SIGMAEQ(1)=ZERO
      SIGMAEQ(2)=ZERO
      SIGMAEQ(3)=ZERO
      SIGMAEQ(4)=ZERO
c      
      dfv=1D0 - (1D0-dftv)*(1D0-dfcv)
      dmv=1D0 - (1D0-dmtv)*(1D0-dmcv)
      dsv = (1D0-smt*dmtv)*(1D0-smc*dmcv)
      dfv = MIN(MAX(dfv,ZERO),0.95D0)
      dmv = MIN(MAX(dmv,ZERO),0.85D0)
      dsv = MIN(MAX(dsv,0.05D0),ONE)

      CALL HASHIN_STIFFNESS(DDSDDE0,NTENS,E11,E22,E33,G12,G13,
     1     G23,PR12,PR13,PR23,PR21,PR31,PR32,DELTA,ZERO,ZERO,ONE)
      CALL HASHIN_STIFFNESS(DDSDDE,NTENS,E11,E22,E33,G12,G13,
     1     G23,PR12,PR13,PR23,PR21,PR31,PR32,DELTA,dfv,dmv,dsv)
C    
      DO I=1,NTENS
          UPSTRAN(I) = STRAN(I) + DSTRAN(I)
      END DO

      DO I=1,NTENS
          STRESS(I) = ZERO
          EFFSTR(I) = ZERO
          DO J=1,NTENS
              STRESS(I) = STRESS(I) + DDSDDE(I,J)*UPSTRAN(J)
              EFFSTR(I) = EFFSTR(I) + DDSDDE0(I,J)*UPSTRAN(J)
          END DO
      END DO
c
c-----
      eps11 = UPSTRAN(1)
      eps22 = UPSTRAN(2)
      eps33 = UPSTRAN(3)
      gam12 = UPSTRAN(4)
      gam13 = UPSTRAN(5)
      gam23 = UPSTRAN(6)

      sig11 = EFFSTR(1)
      sig22 = EFFSTR(2)
      sig33 = EFFSTR(3)
      tau12 = EFFSTR(4)
      tau13 = EFFSTR(5)
      tau23 = EFFSTR(6)      

      DELTAEQ(1) = CELENT*SQRT(MAX(eps11,ZERO)**2 + gam12**2
     1           + gam13**2)
      DELTAEQ(2) = CELENT*MAX(-eps11,ZERO)
      DELTAEQ(3) = CELENT*SQRT(MAX(eps22,ZERO)**2 + MAX(eps33,ZERO)**2
     1           + gam12**2 + gam23**2 + gam13**2)
      DELTAEQ(4) = CELENT*SQRT(MAX(-eps22,ZERO)**2 + MAX(-eps33,ZERO)**2
     1           + gam12**2 + gam23**2 + gam13**2)

      IF (DELTAEQ(1) .GT. ZERO) THEN
          SIGMAEQ(1) = CELENT*(MAX(sig11,ZERO)*MAX(eps11,ZERO)
     1              + tau12*gam12 + tau13*gam13)/DELTAEQ(1)
      END IF

      IF (DELTAEQ(2) .GT. ZERO) THEN
          SIGMAEQ(2) = CELENT*(MAX(-sig11,ZERO)*MAX(-eps11,ZERO))
     1              /DELTAEQ(2)
      END IF

      IF (DELTAEQ(3) .GT. ZERO) THEN
          SIGMAEQ(3) = CELENT*(MAX(sig22,ZERO)*MAX(eps22,ZERO)
     1              + MAX(sig33,ZERO)*MAX(eps33,ZERO)
     2              + tau12*gam12 + tau23*gam23 + tau13*gam13)
     3              /DELTAEQ(3)
      END IF

      IF (DELTAEQ(4) .GT. ZERO) THEN
          SIGMAEQ(4) = CELENT*(MAX(-sig22,ZERO)*MAX(-eps22,ZERO)
     1              + MAX(-sig33,ZERO)*MAX(-eps33,ZERO)
     2              + tau12*gam12 + tau23*gam23 + tau13*gam13)
     3              /DELTAEQ(4)
      END IF

c-----
      CALL HASHIN_CRITERION(FI,EFFSTR,NTENS,XT,XC,YT,YC,S12,S13,S23)

      IF (FI(1) .GE. ONE) init_ft = ONE
      IF (FI(2) .GE. ONE) init_fc = ONE
      IF (FI(3) .GE. ONE) init_mt = ONE
      IF (FI(4) .GE. ONE) init_mc = ONE
c
      IF (init_ft .GE. ONE .AND. DELTA0(1) .LE. ZERO) THEN
          DELTA0(1) = DELTAEQ(1)
          SIGMA0(1) = SIGMAEQ(1)
      END IF

      IF (init_fc .GE. ONE .AND. DELTA0(2) .LE. ZERO) THEN
          DELTA0(2) = DELTAEQ(2)
          SIGMA0(2) = SIGMAEQ(2)
      END IF

      IF (init_mt .GE. ONE .AND. DELTA0(3) .LE. ZERO) THEN
          DELTA0(3) = DELTAEQ(3)
          SIGMA0(3) = SIGMAEQ(3)
      END IF

      IF (init_mc .GE. ONE .AND. DELTA0(4) .LE. ZERO) THEN
          DELTA0(4) = DELTAEQ(4)
          SIGMA0(4) = SIGMAEQ(4)
      END IF
c
      IF (init_ft .GE. ONE) THEN
          kappa_ft = MAX(kappa_ft,DELTAEQ(1))
      END IF

      IF (init_fc .GE. ONE) THEN
          kappa_fc = MAX(kappa_fc,DELTAEQ(2))
      END IF

      IF (init_mt .GE. ONE) THEN
          kappa_mt = MAX(kappa_mt,DELTAEQ(3))
      END IF

      IF (init_mc .GE. ONE) THEN
          kappa_mc = MAX(kappa_mc,DELTAEQ(4))
      END IF
c
      deltaf_ft = ZERO
      deltaf_fc = ZERO
      deltaf_mt = ZERO
      deltaf_mc = ZERO      

      IF (SIGMA0(1) .GT. ZERO) THEN
          deltaf_ft = TWO*G_ft/SIGMA0(1)
      END IF

      IF (SIGMA0(2) .GT. ZERO) THEN
          deltaf_fc = TWO*G_fc/SIGMA0(2)
      END IF

      IF (SIGMA0(3) .GT. ZERO) THEN
          deltaf_mt = TWO*G_mt/SIGMA0(3)
      END IF

      IF (SIGMA0(4) .GT. ZERO) THEN
          deltaf_mc = TWO*G_mc/SIGMA0(4)
      END IF
c
      IF (init_ft .GE. ONE) THEN
          IF (deltaf_ft .GT. DELTA0(1) .AND. kappa_ft .GT. DELTA0(1))
     1    THEN
              dft = deltaf_ft*(kappa_ft-DELTA0(1))
     1              /(kappa_ft*(deltaf_ft-DELTA0(1)))
          END IF
      END IF
c
      IF (init_fc .GE. ONE) THEN
          IF (deltaf_fc .GT. DELTA0(2) .AND. kappa_fc .GT. DELTA0(2))
     1    THEN
              dfc = deltaf_fc*(kappa_fc-DELTA0(2))
     1              /(kappa_fc*(deltaf_fc-DELTA0(2)))
          END IF
      END IF
c
      IF (init_mt .GE. ONE) THEN
          IF (deltaf_mt .GT. DELTA0(3) .AND. kappa_mt .GT. DELTA0(3))
     1    THEN
              dmt = deltaf_mt*(kappa_mt-DELTA0(3))
     1              /(kappa_mt*(deltaf_mt-DELTA0(3)))
          END IF
      END IF
c
      IF (init_mc .GE. ONE) THEN
          IF (deltaf_mc .GT. DELTA0(4) .AND. kappa_mc .GT. DELTA0(4))
     1    THEN
              dmc = deltaf_mc*(kappa_mc-DELTA0(4))
     1              /(kappa_mc*(deltaf_mc-DELTA0(4)))
          END IF
      END IF     
c
      dft = MIN(MAX(dft,ZERO),0.95D0)
      dfc = MIN(MAX(dfc,ZERO),0.95D0)
      dmt = MIN(MAX(dmt,ZERO),0.85D0)
      dmc = MIN(MAX(dmc,ZERO),0.85D0)
c
      IF (vf .GT. ZERO .AND. DTIME .GT. ZERO) THEN
          beta_old = vf/(vf + DTIME)
          beta_new = DTIME/(vf + DTIME)
          dftv = beta_old*dftv + beta_new*dft
          dfcv = beta_old*dfcv + beta_new*dfc
      ELSE
          dftv = dft
          dfcv = dfc
      END IF

      IF (vm .GT. ZERO .AND. DTIME .GT. ZERO) THEN
          beta_old = vm/(vm + DTIME)
          beta_new = DTIME/(vm + DTIME)
          dmtv = beta_old*dmtv + beta_new*dmt
          dmcv = beta_old*dmcv + beta_new*dmc
      ELSE
          dmtv = dmt
          dmcv = dmc
      END IF

      dftv = MIN(MAX(dftv,ZERO),0.95D0)
      dfcv = MIN(MAX(dfcv,ZERO),0.95D0)
      dmtv = MIN(MAX(dmtv,ZERO),0.85D0)
      dmcv = MIN(MAX(dmcv,ZERO),0.85D0)

c
      STATEV(1) = dft
      STATEV(2) = dfc
      STATEV(3) = dmt
      STATEV(4) = dmc
      STATEV(5) = kappa_ft
      STATEV(6) = kappa_fc
      STATEV(7) = kappa_mt
      STATEV(8) = kappa_mc
      STATEV(9)  = init_ft
      STATEV(10) = init_fc
      STATEV(11) = init_mt
      STATEV(12) = init_mc
      STATEV(13) = DELTA0(1)
      STATEV(14) = DELTA0(2)
      STATEV(15) = DELTA0(3)
      STATEV(16) = DELTA0(4)
      STATEV(17) = SIGMA0(1)
      STATEV(18) = SIGMA0(2)
      STATEV(19) = SIGMA0(3)
      STATEV(20) = SIGMA0(4)
      STATEV(21) = dftv
      STATEV(22) = dfcv
      STATEV(23) = dmtv
      STATEV(24) = dmcv
c---

c-----
      RETURN
      END
c-----Hashin stiffness subroutine
      SUBROUTINE HASHIN_STIFFNESS(DDSDDE,NTENS,E11,E22,E33,
     1 G12,G13,G23,PR12,PR13,PR23,PR21,PR31,PR32,DELTA,df,dm,ds)

      INCLUDE 'ABA_PARAM.INC'

      DIMENSION DDSDDE(NTENS,NTENS)
      INTEGER I,J
      DOUBLE PRECISION E11,E22,E33,G12,G13,G23
      DOUBLE PRECISION PR12,PR13,PR23,PR21,PR31,PR32,DELTA,df,dm,ds

      PARAMETER (ZERO = 0.0D0,ONE=1.0D0, TWO=2.0D0)

c----initialize the constitutive matrix
      DO I=1,NTENS
          DO J=1,NTENS
              DDSDDE(I,J) = ZERO
          END DO
      END DO

      DDSDDE(1,1) = (E11*DELTA*(ONE-PR23*PR32))*(ONE-df)
      DDSDDE(1,2) = (E11*DELTA*(PR21+PR31*PR23))*(ONE-df)*(ONE-dm)
      DDSDDE(1,3) = (E11*DELTA*(PR31+PR21*PR32))*(ONE-df)*(ONE-dm)
      DDSDDE(2,1) = (E11*DELTA*(PR21+PR31*PR23))*(ONE-df)*(ONE-dm)
      DDSDDE(2,2) = (E22*DELTA*(ONE-PR13*PR31))*(ONE-dm)*(ONE-df)
      DDSDDE(2,3) = (E22*DELTA*(PR32+PR12*PR31))*(ONE-dm)*(ONE-df)
      DDSDDE(3,1) = (E11*DELTA*(PR31+PR21*PR32))*(ONE-df)*(ONE-dm)
      DDSDDE(3,2) = (E22*DELTA*(PR32+PR12*PR31))*(ONE-dm)*(ONE-df)
      DDSDDE(3,3) = (E33*DELTA*(ONE-PR12*PR21))*(ONE-dm)*(ONE-df)
      DDSDDE(4,4) = G12*ds*(ONE-df)
      DDSDDE(5,5) = G13*ds*(ONE-df)
      DDSDDE(6,6) = G23*ds*(ONE-df)

      RETURN
      END
c
c----Hashin failure criterion subroutine
      SUBROUTINE HASHIN_CRITERION(FI,STRESS,NTENS,XT,XC,YT,YC,S12,
     1 S13,S23)

      INCLUDE 'ABA_PARAM.INC'

      DIMENSION STRESS(NTENS)
      DOUBLE PRECISION FI(4)
      DOUBLE PRECISION XT,XC,YT,YC,S12,S13,S23
      DOUBLE PRECISION sig11,sig22,sig33,tau12,tau13,tau23

      PARAMETER (ZERO = 0.0D0,ONE=1.0D0, TWO=2.0D0)
c
      sig11=STRESS(1) 
      sig22=STRESS(2) 
      sig33=STRESS(3) 
      tau12=STRESS(4) 
      tau13=STRESS(5) 
      tau23=STRESS(6)
c
      FI(1) = ZERO
      FI(2) = ZERO
      FI(3) = ZERO
      FI(4) = ZERO

      IF (sig11 .GE. ZERO) THEN
          FI(1) = (sig11/XT)**2 + (tau12/S12)**2 + (tau13/S13)**2
      ELSE
          FI(2) = (-sig11/XC)**2
      END IF

      IF ((sig22 + sig33) .GE. ZERO) THEN
          FI(3) = ((sig22 + sig33)/YT)**2
     1          + (tau23**2-sig22*sig33)/(S23**2)
     2          + (tau12/S12)**2 + (tau13/S13)**2
      ELSE
          FI(4) = (ONE/YC)*((YC/(TWO*S23))**2-ONE)*(sig22+sig33)
     1          + ((sig22+sig33)**2)/((TWO*S23)**2)
     2          + (tau23**2-sig22*sig33)/(S23**2)
     3          + (tau12**2)/(S12**2)
     4          + (tau13**2)/(S13**2)
      END IF

      RETURN
      END

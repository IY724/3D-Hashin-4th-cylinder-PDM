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

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS), STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS), DDSDDT(NTENS), DRPLDE(NTENS),
     2 STRAN(NTENS), DSTRAN(NTENS), TIME(2), PREDEF(1), DPRED(1),
     3 PROPS(NPROPS), COORDS(3), DROT(3,3), DFGRD0(3,3), DFGRD1(3,3),
     4 JSTEP(4)

      INTEGER I, J
      DOUBLE PRECISION FI(4)
      DOUBLE PRECISION E11, E22, E33, PR12, PR13, PR23, G12, G13, G23
      DOUBLE PRECISION XT, XC, YT, YC, S12, S13, S23
      DOUBLE PRECISION PR21, PR31, PR32, DELTA
      DOUBLE PRECISION SIG11, SIG22, SIG33, TAU12, TAU13, TAU23
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)

C     PROPS(1:9): orthotropic elastic constants
      E11  = PROPS(1)
      E22  = PROPS(2)
      E33  = PROPS(3)
      PR12 = PROPS(4)
      PR13 = PROPS(5)
      PR23 = PROPS(6)
      G12  = PROPS(7)
      G13  = PROPS(8)
      G23  = PROPS(9)

C     PROPS(10:16): Hashin strengths
      XT  = PROPS(10)
      XC  = PROPS(11)
      YT  = PROPS(12)
      YC  = PROPS(13)
      S12 = PROPS(14)
      S13 = PROPS(15)
      S23 = PROPS(16)

C     Build elastic Jacobian
      PR21 = E22*PR12/E11
      PR31 = E33*PR13/E11
      PR32 = E33*PR23/E22
      DELTA = ONE/(ONE-PR12*PR21-PR23*PR32-PR13*PR31
     1        -TWO*PR21*PR32*PR13)

      DO I = 1, NTENS
         DO J = 1, NTENS
            DDSDDE(I,J) = ZERO
         END DO
      END DO

      DDSDDE(1,1) = E11*DELTA*(ONE-PR23*PR32)
      DDSDDE(1,2) = E11*DELTA*(PR21+PR31*PR23)
      DDSDDE(1,3) = E11*DELTA*(PR31+PR21*PR32)
      DDSDDE(2,1) = DDSDDE(1,2)
      DDSDDE(2,2) = E22*DELTA*(ONE-PR13*PR31)
      DDSDDE(2,3) = E22*DELTA*(PR32+PR12*PR31)
      DDSDDE(3,1) = DDSDDE(1,3)
      DDSDDE(3,2) = DDSDDE(2,3)
      DDSDDE(3,3) = E33*DELTA*(ONE-PR12*PR21)
      DDSDDE(4,4) = G12
      DDSDDE(5,5) = G13
      DDSDDE(6,6) = G23

C     Stress update: sigma_{n+1} = sigma_n + C : dEps
      DO I = 1, NTENS
         DO J = 1, NTENS
            STRESS(J) = STRESS(J) + DDSDDE(J,I)*DSTRAN(I)
         END DO
      END DO

      SIG11 = STRESS(1)
      SIG22 = STRESS(2)
      SIG33 = STRESS(3)
      TAU12 = STRESS(4)
      TAU13 = STRESS(5)
      TAU23 = STRESS(6)

C     Hashin failure indices
      FI(1) = ZERO
      FI(2) = ZERO
      FI(3) = ZERO
      FI(4) = ZERO

C     Fiber mode
      IF (SIG11 .GE. ZERO) THEN
         FI(1) = (SIG11/XT)**2 + (TAU12/S12)**2 + (TAU13/S13)**2
      ELSE
         FI(2) = (-SIG11/XC)**2
      END IF

C     Matrix mode
      IF ((SIG22+SIG33) .GE. ZERO) THEN
         FI(3) = ((SIG22+SIG33)/YT)**2
     1         + (TAU23**2-SIG22*SIG33)/(S23**2)
     2         + (TAU12/S12)**2 + (TAU13/S13)**2
      ELSE
         FI(4) = (ONE/YC)*((YC/(TWO*S23))**2-ONE)*(SIG22+SIG33)
     1         + ((SIG22+SIG33)**2)/((TWO*S23)**2)
     2         + (TAU23**2-SIG22*SIG33)/(S23**2)
     3         + (TAU12**2)/(S12**2)
     4         + (TAU13**2)/(S13**2)
      END IF

C     Output Hashin indices to STATEV(9:12)
      STATEV(9)  = FI(1)
      STATEV(10) = FI(2)
      STATEV(11) = FI(3)
      STATEV(12) = FI(4)

      RETURN
      END

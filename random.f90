
MODULE random
IMPLICIT NONE
REAL, PRIVATE      :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0,   &
                      vsmall = TINY(1.0), vlarge = HUGE(1.0)
PRIVATE            :: integral
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)


CONTAINS


FUNCTION random_binomial2(n, pp, first) RESULT(ival)
!**********************************************************************
!     Translated to Fortran 90 by Alan Miller from:
!                              RANLIB
!
!     Library of Fortran Routines for Random Number Generation
!
!                      Compiled and Written by:
!
!                           Barry W. Brown
!                            James Lovato
!
!               Department of Biomathematics, Box 237
!               The University of Texas, M.D. Anderson Cancer Center
!               1515 Holcombe Boulevard
!               Houston, TX      77030
!
! This work was supported by grant CA-16672 from the National Cancer Institute.

!                    GENerate BINomial random deviate

!                              Function

!     Generates a single random deviate from a binomial
!     distribution whose number of trials is N and whose
!     probability of an event in each trial is P.

!                              Arguments

!     N  --> The number of trials in the binomial distribution
!            from which a random deviate is to be generated.
!                              INTEGER N

!     P  --> The probability of an event in each trial of the
!            binomial distribution from which a random deviate
!            is to be generated.
!                              REAL P

!     FIRST --> Set FIRST = .TRUE. for the first call to perform initialization
!               the set FIRST = .FALSE. for further calls using the same pair
!               of parameter values (N, P).
!                              LOGICAL FIRST

!     random_binomial2 <-- A random deviate yielding the number of events
!                from N independent trials, each of which has
!                a probability of event P.
!                              INTEGER random_binomial

!                              Method

!     This is algorithm BTPE from:

!         Kachitvichyanukul, V. and Schmeiser, B. W.
!         Binomial Random Variate Generation.
!         Communications of the ACM, 31, 2 (February, 1988) 216.

!**********************************************************************

!*****DETERMINE APPROPRIATE ALGORITHM AND WHETHER SETUP IS NECESSARY

!     ..
!     .. Scalar Arguments ..
REAL, INTENT(IN)    :: pp
INTEGER, INTENT(IN) :: n
LOGICAL, INTENT(IN) :: first
INTEGER             :: ival
!     ..
!     .. Local Scalars ..
REAL            :: alv, amaxp, f, f1, f2, u, v, w, w2, x, x1, x2, ynorm, z, z2
REAL, PARAMETER :: zero = 0.0, half = 0.5, one = 1.0
INTEGER         :: i, ix, ix1, k, mp
INTEGER, SAVE   :: m
REAL, SAVE      :: p, q, xnp, ffm, fm, xnpq, p1, xm, xl, xr, c, al, xll,  &
                   xlr, p2, p3, p4, qn, r, g

!     ..
!     .. Executable Statements ..

!*****SETUP, PERFORM ONLY WHEN PARAMETERS CHANGE

IF (first) THEN
  p = MIN(pp, one-pp)
  q = one - p
  xnp = n * p
END IF

IF (xnp > 30.) THEN
  IF (first) THEN
    ffm = xnp + p
    m = ffm
    fm = m
    xnpq = xnp * q
    p1 = INT(2.195*SQRT(xnpq) - 4.6*q) + half
    xm = fm + half
    xl = xm - p1
    xr = xm + p1
    c = 0.134 + 20.5 / (15.3 + fm)
    al = (ffm-xl) / (ffm - xl*p)
    xll = al * (one + half*al)
    al = (xr - ffm) / (xr*q)
    xlr = al * (one + half*al)
    p2 = p1 * (one + c + c)
    p3 = p2 + c / xll
    p4 = p3 + c / xlr
  END IF

!*****GENERATE VARIATE, Binomial mean at least 30.

  20 CALL RANDOM_NUMBER(u)
  u = u * p4
  CALL RANDOM_NUMBER(v)

!     TRIANGULAR REGION

  IF (u <= p1) THEN
    ix = xm - p1 * v + u
    GO TO 110
  END IF

!     PARALLELOGRAM REGION

  IF (u <= p2) THEN
    x = xl + (u-p1) / c
    v = v * c + one - ABS(xm-x) / p1
    IF (v > one .OR. v <= zero) GO TO 20
    ix = x
  ELSE

!     LEFT TAIL

    IF (u <= p3) THEN
      ix = xl + LOG(v) / xll
      IF (ix < 0) GO TO 20
      v = v * (u-p2) * xll
    ELSE

!     RIGHT TAIL

      ix = xr - LOG(v) / xlr
      IF (ix > n) GO TO 20
      v = v * (u-p3) * xlr
    END IF
  END IF

!*****DETERMINE APPROPRIATE WAY TO PERFORM ACCEPT/REJECT TEST

  k = ABS(ix-m)
  IF (k <= 20 .OR. k >= xnpq/2-1) THEN

!     EXPLICIT EVALUATION

    f = one
    r = p / q
    g = (n+1) * r
    IF (m < ix) THEN
      mp = m + 1
      DO i = mp, ix
        f = f * (g/i-r)
      END DO

    ELSE IF (m > ix) THEN
      ix1 = ix + 1
      DO i = ix1, m
        f = f / (g/i-r)
      END DO
    END IF

    IF (v > f) THEN
      GO TO 20
    ELSE
      GO TO 110
    END IF
  END IF

!     SQUEEZING USING UPPER AND LOWER BOUNDS ON LOG(F(X))

  amaxp = (k/xnpq) * ((k*(k/3. + .625) + .1666666666666)/xnpq + half)
  ynorm = -k * k / (2.*xnpq)
  alv = LOG(v)
  IF (alv<ynorm - amaxp) GO TO 110
  IF (alv>ynorm + amaxp) GO TO 20

!     STIRLING'S (actually de Moivre's) FORMULA TO MACHINE ACCURACY FOR
!     THE FINAL ACCEPTANCE/REJECTION TEST

  x1 = ix + 1
  f1 = fm + one
  z = n + 1 - fm
  w = n - ix + one
  z2 = z * z
  x2 = x1 * x1
  f2 = f1 * f1
  w2 = w * w
  IF (alv - (xm*LOG(f1/x1) + (n-m+half)*LOG(z/w) + (ix-m)*LOG(w*p/(x1*q)) +    &
      (13860.-(462.-(132.-(99.-140./f2)/f2)/f2)/f2)/f1/166320. +               &
      (13860.-(462.-(132.-(99.-140./z2)/z2)/z2)/z2)/z/166320. +                &
      (13860.-(462.-(132.-(99.-140./x2)/x2)/x2)/x2)/x1/166320. +               &
      (13860.-(462.-(132.-(99.-140./w2)/w2)/w2)/w2)/w/166320.) > zero) THEN
    GO TO 20
  ELSE
    GO TO 110
  END IF

ELSE
!     INVERSE CDF LOGIC FOR MEAN LESS THAN 30
  IF (first) THEN
    qn = q ** n
    r = p / q
    g = r * (n+1)
  END IF

  90 ix = 0
  f = qn
  CALL RANDOM_NUMBER(u)
  100 IF (u >= f) THEN
    IF (ix > 110) GO TO 90
    u = u - f
    ix = ix + 1
    f = f * (g/ix - r)
    GO TO 100
  END IF
END IF

110 IF (pp > half) ix = n - ix
ival = ix
RETURN

END FUNCTION random_binomial2


END MODULE random

! Copyright 2013-2014, Schmidt


! Alpha*X*Y
! (D)ouble precision (A)lpha (X) (T)imes (Y)
! Modified from DAXPY, originally: jack dongarra, linpack, 3/11/78.
      SUBROUTINE DAXTY(N, DA, DX, INCX, DY, INCY)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER INCX,INCY,N
      DOUBLE PRECISION DA, DX(N), DY(N)
      ! Local
      INTEGER I, IX, IY, M
      ! Intrinsic
      INTRINSIC MOD
      
      
      IF (N.LE.0) RETURN
      IF (DA.EQ.0.0D0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
        DO I = M, N, 1
          DY(I) = DA*DY(I)*DX(I)
        END DO
      ELSE
        IX = 1
        IY = 1
        IF (INCX.LT.0) IX = (-N+1)*INCX + 1
        IF (INCY.LT.0) IY = (-N+1)*INCY + 1
        DO I = 1, N
          DY(IY) = DY(IY)*DA*DX(IX)
          IX = IX + INCX
          IY = IY + INCY
        END DO
      END IF
      
      RETURN
      END



! Modification of DNRM2 for more general Minkowski 
! Original code is:
C*This version written on 25-October-1982.
!     Modified on 14-October-1993 to inline the call to DLASSQ.
!     Sven Hammarling, Nag Ltd.
! Modified code by Drew Schmidt, 2014
      DOUBLE PRECISION FUNCTION DNRM3(N, X, INCX, P)
      IMPLICIT NONE
!     .. Scalar Arguments ..
      INTEGER INCX,N, P
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION X(*)
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION ABSXI,NORM,SCALE,SSQ
      INTEGER IX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
!     ..
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE IF (N.EQ.1) THEN
          NORM = ABS(X(1))
      ELSE
          SCALE = ZERO
          SSQ = ONE
!        The following loop is equivalent to this call to the LAPACK
!        auxiliary routine:
!        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (X(IX).NE.ZERO) THEN
                  ABSXI = ABS(X(IX))
                  IF (SCALE.LT.ABSXI) THEN
                      SSQ = ONE + SSQ* (SCALE/ABSXI)**P
                      SCALE = ABSXI
                  ELSE
                      SSQ = SSQ + (ABSXI/SCALE)**P
                  END IF
              END IF
   10     CONTINUE
          NORM = SCALE*((SSQ)**DBLE(1.0d0/DBLE(P)))
      END IF
      dnrm3 = NORM
      RETURN
!     End of dnrm3.
      END






! Location of element having min absolute value
      INTEGER FUNCTION IDAMIN(N,DX,INCX)
      IMPLICIT NONE
      INTEGER INCX,N
      DOUBLE PRECISION DX(*)
!
! Original code is:
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
! Modified code by Drew Schmidt, 2013
      DOUBLE PRECISION DMAX
      INTEGER I,IX
      INTRINSIC DABS
!     ..
      IDAMIN = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IDAMIN = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) THEN
!        code for increment equal to 1
         DMAX = DABS(DX(1))
         DO I = 2,N
            IF (DABS(DX(I)).LT.DMAX) THEN
               IDAMIN = I
               DMAX = DABS(DX(I))
            END IF
         END DO
      ELSE
!        code for increment not equal to 1
         IX = 1
         DMAX = DABS(DX(1))
         IX = IX + INCX
         DO I = 2,N
            IF (DABS(DX(IX)).LT.DMAX) THEN
               IDAMIN = I
               DMAX = DABS(DX(IX))
            END IF
            IX = IX + INCX
         END DO
      END IF
      RETURN
      END



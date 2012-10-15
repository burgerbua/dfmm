! Copyright (c) 2012, Matthias Messner
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
! 
!  * Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!  * Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
! THE POSSIBILITY OF SUCH DAMAGE.

      SUBROUTINE ZHAD(TYPE,N,ZA,ZX,INCX,ZY,INCY)
*     .. Scalar Arguments ..
      DOUBLE COMPLEX ZA
      INTEGER INCX,INCY,N
      CHARACTER TYPE
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX ZX(*),ZY(*)
*     ..
*
*  Purpose
*  =======
*
*     constant times a vector times|divided a vector.
*     matthias messner, 3/16/2011
*

*     .. Local Scalars ..
      INTEGER I,IX,IY
*     ..
*     .. External Functions ..
      DOUBLE PRECISION DCABS1
      EXTERNAL DCABS1
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
      IF (N.LE.0) RETURN
      IF (DCABS1(ZA).EQ.0.0d0) RETURN
*     ..

*     
*     multiply
*     
      IF (LSAME(TYPE,'M')) THEN
         IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*     
*     code for unequal increments or equal increments
*     not equal to 1
*     
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO 10 I = 1,N
            ZY(IY) = ZY(IY) * ZA*ZX(IX)
            IX = IX + INCX
            IY = IY + INCY
 10      CONTINUE
         RETURN
*     
*     code for both increments equal to 1
*     
 20      DO 30 I = 1,N
            ZY(I) = ZY(I) * ZA*ZX(I)
 30      CONTINUE
         RETURN
*     
*     divide
*     
      ELSE IF (LSAME(TYPE,'D')) THEN
         IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 50
*     
*     code for unequal increments or equal increments
*     not equal to 1
*     
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO 40 I = 1,N
            ZY(IY) = ZY(IY) / (ZA*ZX(IX))
            IX = IX + INCX
            IY = IY + INCY
 40      CONTINUE
         RETURN
*     
*     code for both increments equal to 1
*     
 50      DO 60 I = 1,N
            ZY(I) = ZY(I) / (ZA*ZX(I))
 60      CONTINUE
         RETURN
      ELSE
         CALL XERBLA('ZHAD ',1)
         RETURN
      END IF
      END
      
*     ..

! MacKinnon response-surface coefficient source: statsmodels 0.14.6.
! Copyright (C) 2006, Jonathan E. Taylor
! All rights reserved.
!
! Copyright (c) 2006-2008 Scipy Developers.
! All rights reserved.
!
! Copyright (c) 2009-2018 statsmodels Developers.
! All rights reserved.
!
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
!   a. Redistributions of source code must retain the above copyright notice,
!      this list of conditions and the following disclaimer.
!   b. Redistributions in binary form must reproduce the above copyright
!      notice, this list of conditions and the following disclaimer in the
!      documentation and/or other materials provided with the distribution.
!   c. Neither the name of statsmodels nor the names of its contributors
!      may be used to endorse or promote products derived from this software
!      without specific prior written permission.
!
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL STATSMODELS OR CONTRIBUTORS BE LIABLE FOR
! ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
! OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
! DAMAGE.

module ADFStationarity
   use SMConstants, only: RP
   use, intrinsic                                     :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   public :: ADFTest, ADFPValue

!
!  ========
   contains
!  ========
!

!
!///////////////////////////////////////////////////////////////////////////
!
!           ADF STATIONARITY PROCEDURES
!
!///////////////////////////////////////////////////////////////////////////
!
      real(kind=RP) function ADFPValue( statistic ) result(p)
!
!        *******************************************************************
!        Evaluate the MacKinnon asymptotic p-value for the ADF statistic.
!        *******************************************************************
!
         implicit none
         real(kind=RP), intent(in)                          :: statistic
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)                                      :: polynomial

         ! MacKinnon (1994), regression=c, N=1, as in statsmodels 0.14.6.
         if ( statistic > 2.74_RP ) then
            p = 1
         else if ( statistic < -18.83_RP ) then
            p = 0
         else
            if ( statistic <= -1.61_RP ) then
               polynomial = 2.1659_RP + statistic * (1.4412_RP + statistic * .038269_RP)
            else
               polynomial = 1.7339_RP + statistic * (.93202_RP + statistic * (-.12745_RP - statistic * .010368_RP))
            end if
            p = .5_RP * erfc(-polynomial / sqrt(2.0_RP))
         end if
      end function ADFPValue

      subroutine ADFTest( x, lag, alpha, statistic, critical, valid )
!
!        *******************************************************************
!        Fit the augmented regression and compute its stationarity statistic.
!        *******************************************************************
!
         implicit none
         real(kind=RP), intent(in)                          :: x(:), alpha
         integer, intent(in)                                :: lag
         real(kind=RP), intent(out)                         :: statistic, critical
         logical, intent(out)                               :: valid
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP), allocatable                         :: z(:), a(:,:), y(:), v(:), coef(:), inverse_row(:)
         real(kind=RP)                                      :: sd, norm, diagonal, factor, denom, sse, variance, scale
         integer                                            :: n, m, k, i, j, l, index

         ! MacKinnon (2010), regression=c, N=1. Coefficients transcribed from
         ! statsmodels 0.14.6 tsa/adfvalues.py; copyright and BSD terms retained above.
         real(kind=RP), parameter                           :: cv(4, 3) = reshape([ &
            -3.43035_RP, -6.5393_RP, -16.786_RP, -79.433_RP, &
            -2.86154_RP, -2.8903_RP, -4.234_RP, -40.040_RP, &
            -2.56677_RP, -1.5384_RP, -2.809_RP, 0.0_RP], [4, 3])

         valid = .false.
         statistic = 0
         critical = 0
         n = size(x)
         m = n - lag - 1
         k = lag + 2
         if ( lag < 0 .or. m <= k .or. .not. all(ieee_is_finite(x)) ) return
         index = 0
         if ( abs(alpha - .01_RP) < 1.0e-12_RP ) index = 1
         if ( abs(alpha - .05_RP) < 1.0e-12_RP ) index = 2
         if ( abs(alpha - .10_RP) < 1.0e-12_RP ) index = 3
         if ( index == 0 ) return
         critical = cv(1, index) + cv(2, index) / m + cv(3, index) / real(m, RP)**2 + cv(4, index) / real(m, RP)**3
         allocate ( z(n), a(m, k), y(m), v(m), coef(k), inverse_row(k) )
         z = x - sum(x) / n
         sd = sqrt(sum(z * z) / n)
         if ( sd <= tiny(sd) ) return
         z = z / sd
         do i = 1, m
            j = i + lag + 1
            y(i) = z(j) - z(j - 1)
            a(i, 1) = z(j - 1)
            a(i, 2) = 1
            do l = 1, lag
               a(i, l + 2) = z(j - l) - z(j - l - 1)
            end do
         end do
         ! Column-scaled Householder QR, never normal equations X'X.
         do j = 1, k
            norm = sqrt(sum(a(:, j)**2))
            if ( norm <= tiny(norm) ) return
            a(:, j) = a(:, j) / norm
         end do
         do j = 1, k
            norm = sqrt(sum(a(j:m, j)**2))
            ! Columns have unit norm. Use a precision-relative rank floor;
            ! a fixed 1e-10 rejects resolvable, smooth native 200-lag records.
            if ( norm <= epsilon(1.0_RP)*real(k,RP) ) return ! Inconclusive rank.
            diagonal = -sign(norm, a(j, j))
            v(j:m) = a(j:m, j)
            v(j) = v(j) - diagonal
            denom = sum(v(j:m)**2)
            if ( denom <= tiny(denom) ) return
            do l = j + 1, k
               factor = 2 * dot_product(v(j:m), a(j:m, l)) / denom
               a(j:m, l) = a(j:m, l) - factor * v(j:m)
            end do
            factor = 2 * dot_product(v(j:m), y(j:m)) / denom
            y(j:m) = y(j:m) - factor * v(j:m)
            a(j, j) = diagonal
            a(j + 1:m, j) = 0
         end do
         do j = k, 1, -1
            coef(j) = (y(j) - dot_product(a(j, j + 1:k), coef(j + 1:k))) / a(j, j)
         end do
         ! Solve R^T w=e1: ||w||^2 is (R^-1 R^-T)_(1,1).
         inverse_row = 0
         do j = 1, k
            scale = 0
            if ( j == 1 ) scale = 1
            inverse_row(j) = (scale - dot_product(a(1:j - 1, j), inverse_row(1:j - 1))) / a(j, j)
         end do
         sse = sum(y(k + 1:m)**2)
         variance = sse / (m - k)
         denom = sqrt(variance * sum(inverse_row**2))
         if ( denom <= tiny(denom) ) return
         statistic = coef(1) / denom
         valid = ieee_is_finite(statistic)
      end subroutine ADFTest
end module ADFStationarity

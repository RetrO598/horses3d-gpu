module FlowTransientDetector
   use SMConstants, only: RP
   use ADFStationarity, only: ADFTest, ADFPValue
   implicit none
   private
   public :: FlowDetector_t
   ! Growing-prefix ADF with an additional relative fluctuation requirement.
!
!  ************************************
!  Transient detection class definition
!  ************************************
!
   type FlowDetector_t
      integer                                            :: lag = 200
      integer                                            :: minimum_samples = 1000
      real(kind=RP)                                      :: alpha = 0.01_RP
      real(kind=RP)                                      :: relative_std_tolerance = .05_RP

      contains
         procedure   :: Check        => DetectorCheck
   end type FlowDetector_t

!
!  ========
   contains
!  ========
!

!
!///////////////////////////////////////////////////////////////////////////
!
!           TRANSIENT DETECTION PROCEDURES
!
!///////////////////////////////////////////////////////////////////////////
!
      subroutine DetectorCheck( self, x, pass, statistic, critical, p_value, reason, relative_std )
!
!        *******************************************************************
!        Check stationarity and the relative fluctuation requirement.
!        *******************************************************************
!
         implicit none
         class(FlowDetector_t), intent(in)                  :: self
         real(kind=RP), intent(in)                          :: x(:)
         logical, intent(out)                               :: pass
         real(kind=RP), intent(out)                         :: statistic, critical, p_value
         ! Reasons: 0 pass, 4 ADF not rejected, 5 undefined, 6 relative std, 7 insufficient samples.
         integer, intent(out)                              :: reason
         real(kind=RP), intent(out), optional               :: relative_std
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)                                      :: mean, sd, ratio
         logical                                            :: valid

         pass = .false.
         statistic = 0
         critical = 0
         p_value = -1
         reason = 7
         if (present(relative_std)) relative_std = -1
         ! This floor applies even to callers outside the monitor manager.
        if (size(x) < self%minimum_samples) return
         ratio = huge(1.0_RP)
         mean = sum(x) / max(1, size(x))
         sd = sqrt(sum((x - mean)**2) / max(1, size(x)))
         if ( abs(mean) > tiny(mean) ) ratio = sd / abs(mean)
         if ( present(relative_std) ) relative_std = ratio
         call ADFTest( x, self % lag, self % alpha, statistic, critical, valid )
         pass = .false.
         p_value = -1
         reason = 5
         if ( .not. valid ) return
         p_value = ADFPValue(statistic)
         reason = 4
         if ( p_value >= self % alpha ) return
         reason = 6
         pass = ratio < self % relative_std_tolerance
         if ( pass ) reason = 0
      end subroutine DetectorCheck
end module FlowTransientDetector

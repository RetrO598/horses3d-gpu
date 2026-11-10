module MassFlowAveraging
   use SMConstants, only: RP

   implicit none
   private
   public :: FlowAverage_t

!
!  *******************************
!  Time averaging class definition
!  *******************************
!
   type FlowAverage_t
      logical                                            :: active = .false.
      real(kind=RP)                                      :: start = 0.0_RP
      real(kind=RP)                                      :: previous_time = 0.0_RP
      real(kind=RP)                                      :: duration = 0.0_RP
      real(kind=RP)                                      :: mean = 0.0_RP
      real(kind=RP)                                      :: m2 = 0.0_RP
      real(kind=RP)                                      :: previous = 0.0_RP

      contains
         procedure   :: Begin        => AverageBegin
         procedure   :: Add          => AverageAdd
   end type FlowAverage_t

!
!  ========
   contains
!  ========
!

!
!///////////////////////////////////////////////////////////////////////////
!
!           TIME AVERAGING PROCEDURES
!
!///////////////////////////////////////////////////////////////////////////
!
      subroutine AverageBegin( self, t, x )
!
!        *******************************************************************
!        Start a new time-averaging interval at the current observation.
!        *******************************************************************
!
         implicit none
         class(FlowAverage_t), intent(inout)                :: self
         real(kind=RP), intent(in)                          :: t, x

         self % active = .true.
         self % start = t
         self % previous_time = t
         self % duration = 0.0_RP
         self % mean = 0.0_RP
         self % m2 = 0.0_RP
         self % previous = x
      end subroutine AverageBegin

      subroutine AverageAdd( self, t, x )
!
!        *******************************************************************
!        Accumulate the piecewise-linear time mean and second central moment.
!        *******************************************************************
!
         implicit none
         class(FlowAverage_t), intent(inout)                :: self
         real(kind=RP), intent(in)                          :: t, x
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)                                      :: h, new_duration, delta, center, local_m2

         if ( .not. self % active ) return
         h = t - self % previous_time
         if ( h <= 0.0_RP ) return

         center = (x + self % previous) / 2
         local_m2 = h * (x - self % previous)**2 / 12
         new_duration = self % duration + h
         delta = center - self % mean
         self % m2 = self % m2 + local_m2 + delta**2 * self % duration * h / new_duration
         self % mean = self % mean + delta * h / new_duration
         self % duration = new_duration
         self % previous_time = t
         self % previous = x
      end subroutine AverageAdd

end module MassFlowAveraging

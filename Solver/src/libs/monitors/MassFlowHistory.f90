module MassFlowHistory
   use SMConstants, only: RP
   implicit none
   private
   public :: FlowHistory_t
!
!  **********************************
!  Mass-flow history class definition
!  **********************************
!
   type FlowHistory_t
      integer                                            :: capacity = 0
      integer                                            :: count = 0
      real(kind=RP)                                      :: last_time = 0
      real(kind=RP), allocatable                         :: values(:)
      real(kind=RP), allocatable                         :: times(:)

      contains
         procedure          :: Init         => HistoryInit
         procedure          :: Reset        => HistoryReset
         procedure          :: Push         => HistoryPush
         procedure, private :: Reserve      => HistoryReserve
         procedure          :: Window       => HistoryWindow
   end type FlowHistory_t

!
!  ========
   contains
!  ========
!

!
!///////////////////////////////////////////////////////////////////////////
!
!           MASS-FLOW HISTORY PROCEDURES
!
!///////////////////////////////////////////////////////////////////////////
!
      subroutine HistoryInit( self, capacity )
!
!        *******************************************************************
!        Allocate storage for actual complete-step observations.
!        *******************************************************************
!
         implicit none
         class(FlowHistory_t)                               :: self
         integer, intent(in)                                :: capacity

         self % capacity = capacity
         if ( allocated(self % values) ) deallocate ( self % values, self % times )
         allocate ( self % values(capacity), self % times(capacity) )
         self % values = 0
         self % times = 0
         call self % Reset( )
      end subroutine HistoryInit

      subroutine HistoryReset( self )
!
!        *******************************************************************
!        Clear the native observation prefix.
!        *******************************************************************
!
         implicit none
         class(FlowHistory_t)                               :: self

         self % count = 0
         self % last_time = 0
      end subroutine HistoryReset

      subroutine HistoryReserve( self, required )
!
!        *******************************************************************
!        Grow storage while preserving chronological sample order.
!        *******************************************************************
!
         implicit none
         class(FlowHistory_t), intent(inout)                :: self
         integer, intent(in)                                :: required
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP), allocatable                         :: values(:), times(:)
         integer                                            :: capacity

         if ( required <= self % capacity ) return
         capacity = max(required, 2 * self % capacity)
         allocate ( values(capacity), times(capacity) )
         values = 0
         times = 0
         values(1:self%count) = self%values(1:self%count)
         times(1:self%count) = self%times(1:self%count)
         call move_alloc( values, self % values )
         call move_alloc( times, self % times )
         self % capacity = capacity
      end subroutine HistoryReserve

      subroutine HistoryPush( self, t, x )
!
!        Retain one actual observation per complete CFD step. There is no
!        time grid, thinning, interpolation, or synthesis of extra samples.
!
         class(FlowHistory_t), intent(inout) :: self
         real(kind=RP), intent(in) :: t, x

         call self%Reserve(self%count+1)
         self%count = self%count+1
         self%values(self%count) = x
         self%times(self%count) = t
         self%last_time = t
      end subroutine HistoryPush

      subroutine HistoryWindow( self, x )
!
!        *******************************************************************
!        Copy the requested history window in chronological order.
!        *******************************************************************
!
         implicit none
         class(FlowHistory_t), intent(in)                   :: self
         real(kind=RP), intent(out)                         :: x(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                                            :: n

         n = size(x)
         x = self%values(self%count-n+1:self%count)
      end subroutine HistoryWindow

end module MassFlowHistory

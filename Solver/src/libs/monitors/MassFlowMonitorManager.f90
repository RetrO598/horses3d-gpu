#include "Includes.h"
module MassFlowMonitorManager
   use SMConstants, only: RP, LINE_LENGTH
   use FTValueDictionaryClass, only: FTValueDictionary
   use HexMeshClass, only: HexMesh
   use ParamfileRegions, only: readValueInRegion
   use MonitorsClass, only: Monitor_t
   use MPI_Process_Info, only: MPI_Process
   use FileReaders, only: controlFileName
   use Utilities, only: toLower
   use, intrinsic                                       :: iso_fortran_env, only: int64
   use MassFlowHistory, only: FlowHistory_t
   use FlowTransientDetector, only: FlowDetector_t
   use MassFlowAveraging, only: FlowAverage_t
   use, intrinsic                                       :: ieee_arithmetic, only: ieee_is_finite
#ifdef _HAS_MPI_
   use mpi
#endif
   implicit none
   private
   public :: MassFlowManager_t, MassFlowStatistics_t, GetSolutionExtensionBounds

   ! Standalone monitoring: WAITING -> AVERAGING.
   ! Controller-managed monitoring: CONTROL_WAITING <-> DETECTED -> BLOCK_AVERAGING.
   integer, parameter                                   :: WAITING = 0, DETECTED = 1, AVERAGING = 2
   integer, parameter                                   :: CONTROL_WAITING = 3, BLOCK_AVERAGING = 4
   integer, parameter                                   :: OUTPUT_RECORD_LENGTH = 512

!
!  Per-boundary sampling, stationarity detection and averaging state
!  ---------------------------------------------------------------
!
   ! One completed averaging interval; both internal storage and public results
   ! use this representation so that its fields cannot drift out of alignment.
   type FlowBlock_t
      real(kind=RP)                                      :: mean = 0, std = 0, start = 0, finish = 0
      integer                                            :: samples = 0
   end type FlowBlock_t

   type MassFlowStatistics_t
      logical   :: stable = .false.
      integer                                            :: block_count = 0
      type(FlowBlock_t)                                  :: blocks(2)
   end type MassFlowStatistics_t

   type BoundaryMassFlow_t
      logical                                            :: enabled = .false.
      logical                                            :: diagnostic = .false.
      logical                                            :: seen = .false.

      integer                                            :: zone = 0
      integer                                            :: block_count = 0, block_samples = 0
      integer                                            :: check_iteration = -1
      real(kind=RP)                                      :: block_duration = 0
      type(FlowBlock_t)                                  :: blocks(2)
      integer                                            :: monitor = 0
      integer                                            :: hop = 1000
      integer                                            :: state = WAITING
      integer                                            :: last_check = 0
      integer                                            :: last_step = -1
      integer                                            :: reason = -1
      integer                                            :: check_samples = 0

      real(kind=RP)                                      :: diagnostics(4) = [0.0_RP, 0.0_RP, -1.0_RP, -1.0_RP]
      real(kind=RP)                                      :: check_time = -1

      real(kind=RP)                                      :: time_scale = 1
      real(kind=RP)                                      :: flux_scale = 1
      real(kind=RP)                                      :: last_time = 0
      real(kind=RP)                                      :: last_x = 0
      real(kind=RP)                                      :: identity(11) = 0

      character(len=:), allocatable                      :: mesh_name
      character(len=:), allocatable                      :: boundary_name

      type(FlowHistory_t)                                :: history
      type(FlowDetector_t)                               :: detector
      type(FlowAverage_t)                                :: average

      logical                                            :: output_open = .false.
      integer                                            :: output_unit = 0
      integer                                            :: output_count = 0
      integer                                            :: last_output_step = -1
      integer                                            :: output_interval = 100
      integer                                            :: flush_interval = 100
      character(len=:), allocatable                      :: output_name
      character(len=OUTPUT_RECORD_LENGTH), allocatable   :: output_buffer(:)

      contains
         procedure   :: Construct    => FlowConstruct
         procedure   :: Observe      => FlowObserve
         procedure   :: Save         => FlowSave
   end type BoundaryMassFlow_t

!
!  Manager for the enabled boundary monitors
!  -----------------------------------------
!
   type MassFlowManager_t
      type(BoundaryMassFlow_t), allocatable              :: boundaries(:)

      contains
         procedure   :: Construct    => BoundariesConstruct
         procedure   :: Observe      => BoundariesObserve
         procedure   :: Save         => BoundariesSave
         procedure   :: WriteToFile  => BoundariesWriteToFile
         procedure   :: Close        => BoundariesClose
         procedure   :: Bind          => BoundariesBind
         procedure   :: GetStatistics => BoundariesGetStatistics
         procedure   :: BeginAverage  => BoundariesBeginAverage
         procedure   :: Reset         => BoundariesReset
   end type MassFlowManager_t

!
!  ========
   contains
!  ========
!

!
!///////////////////////////////////////////////////////////////////////////
!
!           MASS-FLOW MANAGER PROCEDURES
!
!///////////////////////////////////////////////////////////////////////////
!
      subroutine BoundariesConstruct( self, cfg, monitors, mesh, time_scale, flux_scale, iteration, time )
!
!        *******************************************************************
!        Read boundary settings and bind the active mass-flow monitors.
!        *******************************************************************
!
         implicit none
         class(MassFlowManager_t), intent(inout)            :: self
         type(FTValueDictionary), intent(in)                :: cfg
         class(Monitor_t), intent(in)                       :: monitors
         class(HexMesh), intent(in)                         :: mesh
         real(kind=RP), intent(in)                          :: time_scale, flux_scale, time
         integer, intent(in)                                :: iteration
!
!        ---------------
!        Local variables
!        ---------------
!
         type(FTValueDictionary)                            :: local
         character(len=LINE_LENGTH)                         :: region, value, name
         character(len=32), parameter                       :: options(*) = [character(len=32) :: &
            "history capacity", "hop samples", "adf lags", "significance", &
            "restart policy", "pseudo time diagnostic", "relative std tolerance", "minimum samples"]
         integer                                            :: z, j, k, id, count
         logical                                            :: enabled

!
!        Read enabled boundaries and bind their surface monitors
!        -------------------------------------------------------
!
         allocate ( self % boundaries(size(mesh % zones)) )
         count = 0
         do z = 1, size(mesh % zones)
            name = trim(mesh % zones(z) % name)
            region = "#define boundary " // trim(name)
            call readValueInRegion( trim(controlFileName), "mass flow statistics", value, trim(region), "#end" )
            if ( len_trim(value) == 0 ) cycle ! Default is disabled, independently for every boundary.
            read ( value, *, iostat = j ) enabled
            if ( j /= 0 ) call Fatal( "Invalid mass flow statistics switch on boundary " // trim(name) )
            if ( .not. enabled ) cycle
            id = 0
            do k = 1, monitors % no_of_surfaceMonitors
               if ( .not. monitors % surfaceMonitors(k) % active .or. monitors % surfaceMonitors(k) % IBM ) cycle
               if ( monitors % surfaceMonitors(k) % marker /= z ) cycle
               if ( trim(monitors % surfaceMonitors(k) % variable) /= "mass-flow" ) cycle
               if ( id /= 0 ) call Fatal( "Multiple active mass-flow monitors for boundary " // trim(name) )
               id = k
            end do
            if ( id == 0 ) &
               call Fatal( "Enabled boundary requires an active mass-flow surface monitor: " // trim(name) )
            call local % init( )
            do j = 1, size(options)
               call readValueInRegion( trim(controlFileName), "mass flow statistics " // trim(options(j)), &
                  value, trim(region), "#end" )
               if ( len_trim(value) > 0 ) &
                  call local % addValueForKey( trim(value), "mass flow statistics " // trim(options(j)) )
            end do
            count = count + 1
            call self % boundaries(count) % Construct( cfg, local, monitors, id, name, &
                                                      time_scale, flux_scale, iteration, time )
            call local % destruct( )
         end do
!
!        Display the boundary statistics configuration
!        ---------------------------------------------
!
         if ( count == 0 .or. .not. MPI_Process % isRoot ) return
         print *, "Boundary mass-flow statistics (0=off, 1=on; monitor ID):"
         do z = 1, size(mesh % zones)
            id = 0
            do k = 1, count
               if ( monitors % surfaceMonitors(self % boundaries(k) % monitor) % marker == z ) id = self % boundaries(k) % monitor
            end do
            print *, trim(mesh % zones(z) % name), merge(1, 0, id > 0), id
         end do
      end subroutine BoundariesConstruct

      subroutine BoundariesObserve( self, monitors, iteration, time )
!
!        *******************************************************************
!        Process the current observation for each enabled boundary.
!        *******************************************************************
!
         implicit none
         class(MassFlowManager_t), intent(inout)            :: self
         class(Monitor_t), intent(in)                       :: monitors
         integer, intent(in)                                :: iteration
         real(kind=RP), intent(in)                          :: time
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                                            :: k

         if ( .not. allocated(self % boundaries) ) return
         do k = 1, size(self % boundaries)
            call self % boundaries(k) % Observe( monitors, iteration, time )
         end do
      end subroutine BoundariesObserve

      integer function BoundaryIndex(self, zone) result(id)
         class(MassFlowManager_t), intent(in) :: self
         integer, intent(in) :: zone
         integer :: k
         id = 0
         if (.not. allocated(self%boundaries)) return
         do k = 1, size(self%boundaries)
            if (.not. self%boundaries(k)%enabled) cycle
            if (self%boundaries(k)%zone /= zone) cycle
            id = k
            return
         end do
      end function BoundaryIndex

      subroutine BoundariesBind(self, zones, ids)
         class(MassFlowManager_t), intent(in) :: self
         integer, intent(in) :: zones(:)
         integer, intent(out) :: ids(size(zones))
         integer :: j
         ! Resolve and validate once; managed boundaries remain fixed thereafter.
         do j = 1, size(zones)
            ids(j) = BoundaryIndex(self, zones(j))
            if (ids(j) == 0) call Fatal("Every observed boundary needs enabled mass-flow statistics")
         end do
      end subroutine BoundariesBind

      subroutine BoundariesGetStatistics(self, ids, statistics)
         class(MassFlowManager_t), intent(in) :: self
         integer, intent(in) :: ids(:)
         type(MassFlowStatistics_t), intent(out) :: statistics(size(ids))
         integer :: j
         do j = 1, size(ids)
            associate(b => self%boundaries(ids(j)))
               statistics(j)%stable = (b%state == DETECTED .or. b%state == BLOCK_AVERAGING .or. b%state == AVERAGING) .and. b%reason == 0
               statistics(j)%block_count = b%block_count
               statistics(j)%blocks = b%blocks
            end associate
         end do
      end subroutine BoundariesGetStatistics

      subroutine BoundariesBeginAverage(self, ids, duration)
         class(MassFlowManager_t), intent(inout) :: self
         integer, intent(in) :: ids(:)
         real(kind=RP), intent(in) :: duration
         integer :: j
         if (.not. MPI_Process%isRoot) return
         ! The controller starts the entire bound group from one complete step.
         do j = 1, size(ids)
            associate(b => self%boundaries(ids(j)))
               call WriteAverage(b)
               b%state = BLOCK_AVERAGING
               b%block_duration = duration
               b%block_count = 0
               b%block_samples = 1
               b%blocks = FlowBlock_t()
               b%last_output_step = -1
               call b%average%Begin(b%last_time, b%last_x)
            end associate
         end do
      end subroutine BoundariesBeginAverage

      subroutine BoundariesReset(self, ids, iteration, time)
         class(MassFlowManager_t), intent(inout) :: self
         integer, intent(in) :: ids(:), iteration
         real(kind=RP), intent(in) :: time
         integer :: j
         if (.not. MPI_Process%isRoot) return
         do j = 1, size(ids)
            associate(b => self%boundaries(ids(j)))
               call WriteAverage(b)
               call b%history%Reset()
               b%average = FlowAverage_t()
               b%state = CONTROL_WAITING
               b%check_iteration = -1
               b%last_check = 0
               b%reason = -1
               b%diagnostics = [0.0_RP, 0.0_RP, -1.0_RP, -1.0_RP]
               b%check_time = -1
               b%check_samples = 0
               b%seen = .false.
               b%last_step = iteration
               b%last_x = 0
               b%last_time = time*b%time_scale
               b%last_output_step = -1
               b%block_duration = 0
               b%block_count = 0
               b%block_samples = 0
               b%blocks = FlowBlock_t()
            end associate
         end do
      end subroutine BoundariesReset

      subroutine BoundariesSave( self, u, iteration, time )
!
!        *******************************************************************
!        Write boundary records to the shared solution extension.
!        *******************************************************************
!
         implicit none
         class(MassFlowManager_t), intent(inout)            :: self
         integer, intent(in)                                :: u, iteration
         real(kind=RP), intent(in)                          :: time
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                                            :: k, io, number

         number = 0
         if (allocated(self%boundaries)) number = count(self%boundaries%enabled)
         call self % WriteToFile( force = .true. )
         write ( u, iostat=io ) number
         if (io /= 0) call Fatal("Cannot write embedded boundary count")
         if (.not. allocated(self%boundaries)) return
         do k = 1, size(self % boundaries)
            if ( self % boundaries(k) % enabled ) call self % boundaries(k) % Save( u, iteration, time )
         end do
      end subroutine BoundariesSave

      subroutine BoundariesWriteToFile( self, force )
!
!        *******************************************************************
!        Flush buffered records and optionally record the final averages.
!        *******************************************************************
!
         implicit none
         class(MassFlowManager_t), intent(inout)            :: self
         logical, intent(in), optional                      :: force
!
!        ---------------
!        Local variables
!        ---------------
!
         logical                                            :: force_output
         integer                                            :: k

         if ( .not. MPI_Process % isRoot .or. .not. allocated(self % boundaries) ) return
         force_output = .false.
         if ( present(force) ) force_output = force
         do k = 1, size(self % boundaries)
            if ( .not. self % boundaries(k) % enabled ) cycle
            if ( force_output ) call WriteAverage( self % boundaries(k) )
            call FlowWriteToFile( self % boundaries(k), force_output .or. &
                                 modulo(self % boundaries(k) % last_step, self % boundaries(k) % flush_interval) == 0 )
         end do
      end subroutine BoundariesWriteToFile

      subroutine BoundariesClose( self )
!
!        *******************************************************************
!        Flush the remaining records and close the statistics files.
!        *******************************************************************
!
         implicit none
         class(MassFlowManager_t), intent(inout)            :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                                            :: k, io

         if ( .not. MPI_Process % isRoot .or. .not. allocated(self % boundaries) ) return
         call self % WriteToFile( force = .true. )
         do k = 1, size(self % boundaries)
            associate ( boundary => self % boundaries(k) )
               if ( .not. boundary % output_open ) cycle
               close ( boundary % output_unit, iostat = io )
               if ( io /= 0 ) call Fatal( "Cannot close mass-flow statistics: " // boundary % output_name )
               boundary % output_open = .false.
               deallocate ( boundary % output_buffer )
            end associate
         end do
      end subroutine BoundariesClose

!
!///////////////////////////////////////////////////////////////////////////
!
!           BOUNDARY OBSERVATION PROCEDURES
!
!///////////////////////////////////////////////////////////////////////////
!
      subroutine FlowConstruct( self, cfg, options, monitors, monitor, boundary_name, time_scale, flux_scale, iteration, time )
!
!        *******************************************************************
!        Initialize one boundary and restore its statistical state if requested.
!        *******************************************************************
!
         implicit none
         class(BoundaryMassFlow_t), intent(inout)           :: self
         type(FTValueDictionary), intent(in)                :: cfg, options

         integer, intent(in)                                :: monitor
         character(len=*), intent(in)                       :: boundary_name
         class(Monitor_t), intent(in)                       :: monitors
         real(kind=RP), intent(in)                          :: time_scale, flux_scale, time
         integer, intent(in)                                :: iteration
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                                            :: capacity, output_interval, flush_interval
         character(len=LINE_LENGTH)                         :: policy, restart
         character(len=:), allocatable                      :: output_name

!
!        Read output settings and validate the integration mode
!        ------------------------------------------------------
!
         self % enabled = .true.
         output_interval = cfg % getValueOrDefault("output interval", 100)
         flush_interval = cfg % getValueOrDefault("monitors flush interval", 100)
         if ( output_interval < 1 .or. flush_interval < 1 ) &
            call Fatal( "Mass-flow output intervals must be positive" )
         self % diagnostic = options % getValueOrDefault("mass flow statistics pseudo time diagnostic", .false.)
         policy = cfg % stringValueForKey("time integration", LINE_LENGTH)
         call toLower( policy )
         if ( trim(policy) /= "explicit" ) call Fatal( "Only explicit integration is supported" )
         policy = cfg % stringValueForKey("simulation type", LINE_LENGTH)
         call toLower( policy )
         if ( trim(policy) == "steady-state" .and. .not. self % diagnostic ) &
            call Fatal( "Steady-state sampling requires explicit pseudo time diagnostic mode" )

         ! Read boundary options, using the same dictionary defaults as other source modules.
         self % monitor = monitor
         self % zone = monitors % surfaceMonitors(monitor) % marker
         self % hop = options % getValueOrDefault("mass flow statistics hop samples", 1000)
         capacity = options % getValueOrDefault("mass flow statistics history capacity", max(4096, self % hop))
         self % detector % lag = options % getValueOrDefault("mass flow statistics adf lags", self%detector%lag)
         self % detector % minimum_samples = options % getValueOrDefault("mass flow statistics minimum samples", &
                                                                      max(1000,2*self%detector%lag+4))
         self % detector % alpha = options % getValueOrDefault("mass flow statistics significance", self % detector % alpha)
         self % detector % relative_std_tolerance = &
            options % getValueOrDefault("mass flow statistics relative std tolerance", self % detector % relative_std_tolerance)

         if ( capacity < 1 .or. self % hop < 1 .or. self % detector % lag < 0 ) &
            call Fatal( "Invalid ADF history/lag/sample configuration" )
         if (time_scale <= 0 .or. flux_scale <= 0) call Fatal("Reference scales must be positive")
         if (self%detector%minimum_samples < 1000 .or. self%detector%minimum_samples <= 2*self%detector%lag+3) &
            call Fatal("Minimum samples must be >= 1000 and provide residual degrees of freedom for the ADF lag")
         if ( all(abs(self % detector % alpha - [0.01_RP, 0.05_RP, 0.10_RP]) > 1.0e-12_RP) ) &
            call Fatal( "ADF significance must be 0.01, 0.05 or 0.10" )
         if ( .not. ieee_is_finite(self % detector % relative_std_tolerance) .or. self % detector % relative_std_tolerance <= 0 ) &
            call Fatal( "Relative standard deviation tolerance must be finite and positive" )

!
!        Store reference scales and checkpoint configuration
!        ---------------------------------------------------
!
         self % time_scale = time_scale
         self % flux_scale = flux_scale
         self % boundary_name = trim(boundary_name)
         self % mesh_name = trim(cfg % stringValueForKey("mesh file name", LINE_LENGTH))
         self % identity = [real(self % monitor, RP), real(capacity, RP), real(self % hop, RP), real(self % detector % lag, RP), &
            real(self%detector%minimum_samples,RP), time_scale, flux_scale, self % detector % alpha, &
            merge(1.0_RP, 0.0_RP, self % diagnostic), &
            real(monitors % surfaceMonitors(self % monitor) % marker, RP), &
            self % detector % relative_std_tolerance]
         if ( .not. all(ieee_is_finite(self % identity)) ) call Fatal( "Nonfinite ADF configuration" )

         if ( .not. MPI_Process % isRoot ) return
         call self % history % Init( capacity )

!
!        Restore statistical state when requested
!        ----------------------------------------
!
         policy = "reset"
         if ( options % containsKey("mass flow statistics restart policy") ) &
            policy = options % stringValueForKey("mass flow statistics restart policy", LINE_LENGTH)
         call toLower( policy )
         if ( trim(policy) /= "reset" .and. trim(policy) /= "continue" ) &
            call Fatal( "Restart policy must be reset or continue" )
         if ( trim(policy) == "continue" .and. cfg % logicalValueForKey("restart") ) then
            restart = cfg % stringValueForKey("restart file name", LINE_LENGTH)
            call FlowRestore( self, trim(restart), iteration, time * time_scale )
         end if

!
!        Open the statistics output
!        --------------------------
!
         output_name = trim(monitors % solution_file) // "." // self % boundary_name // ".massflow"
         call FlowOpenOutput( self, output_name, iteration, time * time_scale, policy, output_interval, flush_interval )

         if ( self % diagnostic ) print *, "Mass-flow time axis is PSEUDO TIME, not physical time."
      end subroutine FlowConstruct

      subroutine FlowObserve( self, monitors, iteration, time )
!
!        *******************************************************************
!        Process one observation, then advance detection or time averaging.
!        *******************************************************************
!
         implicit none
         class(BoundaryMassFlow_t), intent(inout)           :: self
         class(Monitor_t), intent(in)                       :: monitors
         integer, intent(in)                                :: iteration
         real(kind=RP), intent(in)                          :: time
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)                                      :: t, x, statistic, critical, p_value, relative_std
         real(kind=RP), allocatable                         :: prefix(:)
         logical                                            :: pass
         integer                                            :: reason, b

         if ( .not. self % enabled .or. .not. MPI_Process % isRoot ) return
!
!        Read and dimensionalize the current observation
!        -----------------------------------------------
!
         b = monitors % bufferLine
         t = time * self % time_scale
         x = monitors % surfaceMonitors(self % monitor) % values(b) * self % flux_scale
         if ( .not. ieee_is_finite(x) ) call Fatal( "Nonfinite monitor sample" )
!
!        Check repeated observations and time ordering
!        ---------------------------------------------
!
         if ( self % seen ) then
            if ( iteration == self % last_step ) return
            if ( t <= self % last_time .or. iteration <= self % last_step ) &
               call Fatal( "Nonmonotone flow history" )
         end if
!
!        Retain the current observation
!        ------------------------------
!
         self % last_step = iteration
         self % last_time = t
         self % last_x = x
         self % seen = .true.
!
!        Advance detection or time averaging
!        -----------------------------------
!
         if ( self % state == AVERAGING .or. self%state == BLOCK_AVERAGING ) then
            ! Both criteria passed on the same prefix. Continue averaging thereafter.
            call self % average % Add( t, x )
            if (self%state == BLOCK_AVERAGING) then
               self%block_samples = self%block_samples + 1
               if (self%average%duration >= self%block_duration) then
                  self%blocks(1) = self%blocks(2)
                  self%blocks(2) = FlowBlock_t(self%average%mean, &
                     sqrt(max(0.0_RP,self%average%m2/self%average%duration)), &
                     self%average%start,t,self%block_samples)
                  self%block_count = self%block_count + 1
                  call WriteAverage(self)
                  call self%average%Begin(t, x)
                  self%block_samples = 1
               end if
            end if
         else
            ! Controlled boundaries remain in detection (including rechecks of
            ! a passed boundary) until the controller releases the entire group.
            ! Uncontrolled monitors retain the original start-on-first-pass rule.
            call self % history % Push( t, x )
            if (self%history%count >= self%detector%minimum_samples .and. &
                self%history%count-self%last_check >= self%hop) then
               self % last_check = self % history % count
               allocate ( prefix(self % history % count) )
               call self % history % Window( prefix )
               call self % detector % Check( prefix, pass, statistic, critical, p_value, reason, relative_std )
               self%check_iteration = iteration
               self % reason = reason
               self % diagnostics = [statistic, critical, p_value, relative_std]
               self % check_time = t
               self % check_samples = self % history % count
               call WriteCheck( self )
               if (self%state == CONTROL_WAITING .or. self%state == DETECTED) then
                  if (pass) then
                     self%state = DETECTED
                  else
                     self%state = CONTROL_WAITING
                  end if
               else if ( pass ) then
                  self % state = AVERAGING
                  call self % average % Begin( t, x )
                  call WriteEvent( self, "average-start" )
               end if
            end if
         end if
!
!        Record averages at the configured output interval
!        -------------------------------------------------
!
         if ( modulo(iteration, self % output_interval) == 0 ) call WriteAverage( self )
      end subroutine FlowObserve

!
!///////////////////////////////////////////////////////////////////////////
!
!           STATISTICS OUTPUT PROCEDURES
!
!///////////////////////////////////////////////////////////////////////////
!
      subroutine FlowOpenOutput( self, name, iteration, time, policy, output_interval, flush_interval )
!
!        *******************************************************************
!        Open the statistics output and write the run metadata.
!        *******************************************************************
!
         implicit none
         class(BoundaryMassFlow_t), intent(inout)           :: self
         character(len=*), intent(in)                       :: name, policy
         integer, intent(in)                                :: iteration, output_interval, flush_interval
         real(kind=RP), intent(in)                          :: time
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                                            :: io
         character(len=OUTPUT_RECORD_LENGTH)                :: record

         self % output_name = name
         self % output_interval = output_interval
         self % flush_interval = flush_interval
         allocate ( self % output_buffer(flush_interval) )
         open ( newunit = self % output_unit, file = self % output_name, status = "unknown", position = "append", &
            action = "write", iostat = io )
         if ( io /= 0 ) call Fatal( "Cannot create mass-flow statistics: " // self % output_name )
         self % output_open = .true.
         call QueueOutput( self, "# Boundary: " // self % boundary_name )
         if ( self % diagnostic ) then
            call QueueOutput( self, "# Time basis: pseudo-time; units: s, kg/s" )
         else
            call QueueOutput( self, "# Time basis: physical time; units: s, kg/s" )
         end if
         write ( record, '(A,I0,1X,ES24.16,1X,A)' ) "# Run start: ", iteration, time, trim(policy )
         call QueueOutput( self, record )
         write ( record, '(A,ES24.16,A,ES24.16)' ) &
            "# Criteria: p < ", self % detector % alpha, "; relative std < ", self % detector % relative_std_tolerance
         call QueueOutput( self, record )
         write ( record, '(A,I0,A,I0,A,I0)' ) &
            "# Native complete-step sampling: hop=",self%hop,"; lags=",self%detector%lag, &
            "; minimum samples=",self%detector%minimum_samples
         call QueueOutput( self, record )
         call QueueOutput( self, "# check: iteration time samples adf critical p_value relative_std reason" )
         call QueueOutput( self, "# average-start: iteration time" )
         call QueueOutput( self, "# average: iteration time duration mean std rms" )
         if ( self % state == AVERAGING .or. self%state == BLOCK_AVERAGING ) then
            write ( record, '(A,ES24.16)' ) "# Restored average start: ", self % average % start
            call QueueOutput( self, record )
         end if
         call QueueOutput( self, "# reason: 0=pass, 4=ADF-not-rejected, 5=undefined-regression, 6=relative-std-too-large, 7=insufficient-samples" )
         call FlowWriteToFile( self, .true. )
      end subroutine FlowOpenOutput

      subroutine QueueOutput( self, record )
!
!        *******************************************************************
!        Buffer a formatted record and write it when the buffer is full.
!        *******************************************************************
!
         implicit none
         class(BoundaryMassFlow_t), intent(inout)           :: self
         character(len=*), intent(in)                       :: record

         self % output_count = self % output_count + 1
         self % output_buffer(self % output_count) = record
         call FlowWriteToFile( self, .false. )
      end subroutine QueueOutput

      subroutine FlowWriteToFile( self, force )
!
!        *******************************************************************
!        Write pending records and flush the output stream.
!        *******************************************************************
!
         implicit none
         class(BoundaryMassFlow_t), intent(inout)           :: self
         logical, intent(in)                                :: force
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                                            :: i, io

         if ( .not. self % output_open ) return
         if ( .not. force .and. self % output_count < size(self % output_buffer) ) return
         do i = 1, self % output_count
            write ( self % output_unit, '(A)', iostat = io ) trim(self % output_buffer(i) )
            if ( io /= 0 ) call Fatal( "Cannot write mass-flow statistics: " // self % output_name )
         end do
         self % output_count = 0
         flush ( self % output_unit, iostat = io )
         if ( io /= 0 ) call Fatal( "Cannot flush mass-flow statistics: " // self % output_name )
      end subroutine FlowWriteToFile

      subroutine WriteCheck( self )
!
!        *******************************************************************
!        Record the latest ADF and relative standard deviation results.
!        *******************************************************************
!
         implicit none
         class(BoundaryMassFlow_t), intent(inout)           :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=OUTPUT_RECORD_LENGTH)                :: record

         write ( record, '(A,1X,I0,1X,ES24.16,1X,I0,4(1X,ES24.16),1X,I0)' ) &
            "check", self % last_step, self % check_time, self % check_samples, &
            self % diagnostics, self % reason
         call QueueOutput( self, record )
      end subroutine WriteCheck

      subroutine WriteEvent( self, event )
!
!        *******************************************************************
!        Record the start of standalone time averaging.
!        *******************************************************************
!
         implicit none
         class(BoundaryMassFlow_t), intent(inout)           :: self
         character(len=*), intent(in)                       :: event
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=OUTPUT_RECORD_LENGTH)                :: record

         write ( record, '(A,1X,I0,1X,ES24.16)' ) event, self % last_step, self % last_time
         call QueueOutput( self, record )
      end subroutine WriteEvent

      subroutine WriteAverage( self )
!
!        *******************************************************************
!        Record the current time-weighted mean, standard deviation and RMS.
!        *******************************************************************
!
         implicit none
         class(BoundaryMassFlow_t), intent(inout)           :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)                                      :: std, rms
         character(len=OUTPUT_RECORD_LENGTH)                :: record

         if ( .not. self % output_open .or. self % average % duration <= 0.0_RP ) return
         if ( self % last_output_step == self % last_step ) return
         std = sqrt(max(0.0_RP, self % average % m2 / self % average % duration))
         rms = sqrt(std**2 + self % average % mean**2)
         write ( record, '(A,1X,I0,5(1X,ES24.16))' ) "average", self % last_step, self % last_time, &
            self % average % duration, self % average % mean, std, rms
         call QueueOutput( self, record )
         self % last_output_step = self % last_step
      end subroutine WriteAverage

!
!///////////////////////////////////////////////////////////////////////////
!
!           CHECKPOINT PROCEDURES
!
!///////////////////////////////////////////////////////////////////////////
!
      subroutine FlowRestore( self, restart, iteration, time )
!
!        *******************************************************************
!        Restore the boundary history and averaging state from a checkpoint.
!        *******************************************************************
!
         implicit none
         class(BoundaryMassFlow_t), intent(inout)           :: self
         character(len=*), intent(in)                       :: restart
         integer, intent(in)                                :: iteration
         real(kind=RP), intent(in)                          :: time
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                                            :: io, u, stored_step, stored_capacity
         integer(int64)                                     :: record_end, record_pos
         real(kind=RP)                                      :: old_identity(11), saved_time
         character(len=LINE_LENGTH)                         :: mesh

!
!        Open the checkpoint and locate the boundary record
!        --------------------------------------------------
!
         open ( newunit = u, file = trim(restart), &
                status = "old", access = "stream", &
                form = "unformatted", action = "read", iostat = io )
         if ( io /= 0 ) call Fatal( "Cannot open restart solution: " // trim(restart) )
         call LocateBoundary( u, self % boundary_name, record_end )
!
!        Read and validate the checkpoint configuration
!        ----------------------------------------------
!
         read ( u, iostat = io ) &
            old_identity, &
            mesh, &
            stored_step, &
            saved_time, &
            stored_capacity
         if ( io /= 0 ) call Fatal( "Cannot read embedded ADF header" )
         if ( any(abs(old_identity - self % identity) > &
                  1.0e-12_RP * max(abs(self % identity), tiny(1.0_RP))) ) &
            call Fatal( "Flow checkpoint configuration mismatch" )
         if ( trim(mesh) /= self % mesh_name .or. &
              stored_step /= iteration .or. &
              abs(saved_time - time) > &
                 1.0e-10_RP * max(abs(saved_time), self % time_scale) ) &
            call Fatal( "Flow checkpoint mesh/iteration/time mismatch" )
!
!        Validate storage requirements before allocation
!        -----------------------------------------------
!
         inquire ( unit = u, pos = record_pos )
         if ( stored_capacity < 1 .or. &
              int(stored_capacity, int64) > &
                 (record_end - record_pos) / (2 * storage_size(1.0_RP) / 8) ) &
            call Fatal( "Invalid embedded ADF history capacity" )
!
!        Restore the history, detection and averaging state
!        --------------------------------------------------
!
         call self % history % Init( stored_capacity )
         read ( u, iostat = io ) &
            self % history % count, &
            self % history % last_time, &
            self % history % values, &
            self % history % times, &
            self % state, &
            self % last_check, &
            self % last_step, &
            self % last_time, &
            self % last_x, &
            self % seen, &
            self % average, &
            self % reason, &
            self % diagnostics, &
            self % check_time, &
            self % check_samples, &
            self % check_iteration, &
            self % block_duration, self % block_count, self % block_samples, self % blocks
         inquire ( unit = u, pos = record_pos )
         close ( u )
!
!        Check the record boundary and restored state
!        --------------------------------------------
!
         if ( io /= 0 .or. record_pos /= record_end ) &
            call Fatal( "Truncated or mismatched embedded ADF record" )
         if ( self % history % count < 0 .or. &
              self % history % count > stored_capacity .or. &
              (.not. any(self%state == [WAITING,DETECTED,AVERAGING,CONTROL_WAITING,BLOCK_AVERAGING])) ) &
            call Fatal( "Invalid embedded ADF state" )
      end subroutine FlowRestore

      subroutine FlowSave( self, u, iteration, time )
!
!        *******************************************************************
!        Write one indexed boundary checkpoint record.
!        *******************************************************************
!
         implicit none
         class(BoundaryMassFlow_t), intent(inout)           :: self
         integer, intent(in)                                :: u, iteration
         real(kind=RP), intent(in)                          :: time
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                                            :: io
         integer(int64)                                     :: next, pos, bytes, index_bytes
         character(len=LINE_LENGTH)                         :: mesh, name

         mesh = self % mesh_name
         name = self % boundary_name
         inquire ( iolength = bytes ) &
            self % identity, &
            mesh, &
            iteration, &
            time * self % time_scale, &
            self % history % capacity, &
            self % history % count, &
            self % history % last_time, &
            self % history % values, &
            self % history % times, &
            self % state, &
            self % last_check, &
            self % last_step, &
            self % last_time, &
            self % last_x, &
            self % seen, &
            self % average, &
            self % reason, &
            self % diagnostics, &
            self % check_time, &
            self % check_samples, &
            self % check_iteration, &
            self % block_duration, self % block_count, self % block_samples, self % blocks
         inquire ( iolength = index_bytes ) next, name
         inquire ( unit = u, pos = pos )
         next = pos + index_bytes + bytes
         write ( u, iostat = io ) &
            next, &
            name, &
            self % identity, &
            mesh, &
            iteration, &
            time * self % time_scale, &
            self % history % capacity, &
            self % history % count, &
            self % history % last_time, &
            self % history % values, &
            self % history % times, &
            self % state, &
            self % last_check, &
            self % last_step, &
            self % last_time, &
            self % last_x, &
            self % seen, &
            self % average, &
            self % reason, &
            self % diagnostics, &
            self % check_time, &
            self % check_samples, &
            self % check_iteration, &
            self % block_duration, self % block_count, self % block_samples, self % blocks
         if ( io /= 0 ) call Fatal( "Cannot save embedded boundary ADF statistics" )
      end subroutine FlowSave

      subroutine GetSolutionExtensionBounds(u, start, finish, status)
         integer, intent(in) :: u
         integer(int64), intent(out) :: start, finish
         integer, intent(out) :: status
         integer(int64) :: bytes, footer_bytes, integer_bytes
         integer :: io, terminator

         status = 1
         inquire (iolength=footer_bytes) start, terminator
         inquire (iolength=integer_bytes) terminator
         inquire (unit=u, size=bytes)
         if (bytes < footer_bytes + integer_bytes) return
         finish = bytes - footer_bytes + 1
         read (u, pos=finish, iostat=io) start, terminator
         if (io /= 0) return
         if (terminator /= 99) return
         if (start <= integer_bytes .or. start >= finish) return
         read (u, pos=start - integer_bytes, iostat=io) terminator
         if (io /= 0) return
         if (terminator /= 99) return
         status = 0
      end subroutine GetSolutionExtensionBounds

      subroutine LocateBoundary( u, name, record_end )
!
!        *******************************************************************
!        Locate a boundary record using the checkpoint footer and index.
!        *******************************************************************
!
         implicit none
         integer, intent(in)                                :: u
         character(len=*), intent(in)                       :: name
         integer(int64), intent(out)                        :: record_end
!
!        ---------------
!        Local variables
!        ---------------
!
         integer(int64)                                     :: start, finish, next, pos, index_bytes
         integer                                            :: io, count, k
         character(len=LINE_LENGTH)                         :: stored_name

         call GetSolutionExtensionBounds(u, start, finish, io)
         if (io /= 0) call Fatal("Invalid solution statistics extension")
         read ( u, pos = start, iostat = io ) count
         if ( io /= 0 ) call Fatal( "Cannot read embedded statistics header" )
         inquire (iolength=index_bytes) next, stored_name
         if (count < 1 .or. int(count, int64) > (finish-start)/index_bytes) &
            call Fatal( "Invalid embedded boundary count" )
         do k = 1, count
            inquire ( unit = u, pos = pos )
            if ( pos >= finish ) call Fatal( "Truncated embedded boundary index" )
            read ( u, iostat = io ) next, stored_name
            if ( io /= 0 ) call Fatal( "Cannot read embedded boundary index" )
            inquire ( unit = u, pos = pos )
            if ( next <= pos .or. next > finish ) &
               call Fatal( "Invalid embedded boundary record size" )
            if ( trim(stored_name) == trim(name) ) then
               record_end = next
               return
            end if
            ! POS on an empty I/O list seeks to the next indexed record.
            read ( u, pos = next, iostat = io )
            if ( io /= 0 ) call Fatal( "Cannot seek embedded boundary record" )
         end do
         call Fatal( "No embedded statistics for boundary " // trim(name) // "; use reset" )
      end subroutine LocateBoundary

!
!///////////////////////////////////////////////////////////////////////////
!
!           ERROR HANDLING
!
!///////////////////////////////////////////////////////////////////////////
!
      subroutine Fatal( message )
!
!        *******************************************************************
!        Report a statistics error and terminate the calculation.
!        *******************************************************************
!
         implicit none
         character(len=*), intent(in)                       :: message
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                                            :: ierr
#endif
         print *, "MASS FLOW ANALYSIS: ", message
#ifdef _HAS_MPI_
         if ( MPI_Process % DoMPIAction ) call MPI_Abort( MPI_COMM_WORLD, 99, ierr )
#endif
         error stop 99
      end subroutine Fatal

end module MassFlowMonitorManager

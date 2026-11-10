#include "Includes.h"
! Fixed-pressure CFD response -> signed mass residual -> bounded pressure update.
! Statistics and boundary physics remain owned by their existing modules.
module MassFlowBalanceController
   use SMConstants, only: RP, LINE_LENGTH
   use FTValueDictionaryClass, only: FTValueDictionary
   use ParamfileRegions, only: readValueInRegion
   use FileReaders, only: controlFileName
   use HexMeshClass, only: HexMesh
   use BoundaryConditions, only: BCs
   use OutflowBCClass, only: OutflowBC_t
   use MassFlowMonitorManager, only: MassFlowManager_t, MassFlowStatistics_t, GetSolutionExtensionBounds
   use MPI_Process_Info, only: MPI_Process
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int64
#ifdef _HAS_MPI_
   use mpi
#endif
   implicit none
   private
   public :: MassFlowBalanceController_t

   integer, parameter :: OFF = 0, WAIT_STABLE = 1, COLLECT_MEAN = 2, CONFIRM = 3, BALANCED = 4, ENDED = 5
   integer, parameter :: NO_ACTION = 0, START_MEAN = 1, UPDATE_PRESSURE = 2
   character(len=16), parameter :: RECORD_TAG = "MASSFLOW-CONTROL"
   integer, parameter :: RECORD_VERSION = 1
   character(len=*), parameter :: REGION = "#define mass flow balance"

   type MassFlowBalanceController_t
      logical :: enabled = .false.
      integer :: state = OFF, epoch = 0, evaluations = 0, phase_start = 0
      integer :: minimum_steps = 50
      integer :: log_unit = 0, last_block = 0
      integer, allocatable :: controlled(:), monitored(:)
      real(kind=RP), allocatable :: coefficients(:), pressures(:), residuals(:)
      real(kind=RP) :: pressure = 0, target = 0, flow_scale = 0, flow_floor = 1e-12_RP
      real(kind=RP) :: pmin = 0, pmax = 0, probe = 0, max_step = 0, relaxation = 0.5_RP
      real(kind=RP) :: rtol = 1e-4_RP, atol = 0, mean_fraction = 0.25_RP, block_duration = 0
      character(len=64) :: reason = "disabled", method = "initial"
      character(len=32) :: constraint = "target-mass-flow"
      character(len=LINE_LENGTH), allocatable :: observed_names(:), controlled_names(:)
      character(len=64) :: last_limit_reason = ""
      type(MassFlowStatistics_t), allocatable :: statistics(:)
      character(len=:), allocatable :: log_name
   contains
      procedure :: Construct => BalanceConstruct
      procedure :: AfterStep => BalanceAfterStep
      procedure :: Finish => BalanceFinish
      procedure :: Save => BalanceSave
      procedure :: Close => BalanceClose
   end type MassFlowBalanceController_t


contains

   subroutine BalanceConstruct(self, cfg, mesh, massFlow, iteration, time, solution_name)
      class(MassFlowBalanceController_t), intent(inout) :: self
      type(FTValueDictionary), intent(in) :: cfg
      class(HexMesh), intent(in) :: mesh
      type(MassFlowManager_t), intent(inout) :: massFlow
      integer, intent(in) :: iteration
      real(kind=RP), intent(in) :: time
      character(len=*), intent(in) :: solution_name
      type(FTValueDictionary) :: options
      character(len=LINE_LENGTH) :: value, constraint, names, restart_policy
      character(len=40), parameter :: keys(*) = [character(len=40) :: &
                                   "algorithm", "constraint", "observed boundaries", "controlled boundaries", "flow coefficients", &
                                                 "target mass flow (kg/s)", "reference mass flow (kg/s)", "flow floor (kg/s)", &
                                                 "minimum pressure (pa)", "maximum pressure (pa)", "probe pressure step (pa)", &
                                                 "maximum pressure step (pa)", "relaxation", "relative mass tolerance", &
                                            "absolute mass tolerance (kg/s)", "mean error fraction", "average block duration (s)", &
                                                 "minimum settle steps"]
      integer, allocatable :: observed(:)
      integer :: j, io, bc_id
      logical :: has_inlet, has_outlet
      real(kind=RP) :: initial_pressure, area, global_area
#ifdef _HAS_MPI_
      integer :: ierr
#endif
      call readValueInRegion(trim(controlFileName), "enabled", value, REGION, "#end")
      if (len_trim(value) == 0) return
      read (value, *, iostat=io) self%enabled
      if (io /= 0) call Fatal("Invalid enabled switch")
      if (.not. self%enabled) return
      call options%init()
      do j = 1, size(keys)
         call readValueInRegion(trim(controlFileName), trim(keys(j)), value, REGION, "#end")
         if (len_trim(value) > 0) call options%addValueForKey(trim(value), trim(keys(j)))
      end do
      if (options%containsKey("algorithm")) then
         value = options%stringValueForKey("algorithm", LINE_LENGTH)
         if (trim(value) /= "safeguarded-secant") call Fatal("Unsupported pressure update algorithm")
      end if
      do j = 1, 6
         select case (j)
         case (1); names = "controlled boundaries"
         case (2); names = "minimum pressure (pa)"
         case (3); names = "maximum pressure (pa)"
         case (4); names = "probe pressure step (pa)"
         case (5); names = "maximum pressure step (pa)"
         case (6); names = "average block duration (s)"
         end select
         if (.not. options%containsKey(trim(names))) call Fatal("Missing "//trim(names))
      end do
      names = options%stringValueForKey("controlled boundaries", LINE_LENGTH)
      call ResolveZones(mesh, names, self%controlled)
      constraint = "target-mass-flow"
      if (options%containsKey("constraint")) constraint = options%stringValueForKey("constraint", LINE_LENGTH)
      if (.not. options%containsKey("observed boundaries")) &
         call Fatal("Observed boundaries must explicitly include the inlet and outlet")
      names = options%stringValueForKey("observed boundaries", LINE_LENGTH)
      call ResolveZones(mesh, names, observed)
      allocate (self%coefficients(size(observed)), self%monitored(size(observed)))
      call massFlow%Bind(observed, self%monitored)
      allocate (self%observed_names(size(observed)), self%controlled_names(size(self%controlled)))
      do j = 1, size(observed)
         self%observed_names(j) = mesh%zones(observed(j))%name
      end do
      do j = 1, size(self%controlled)
         self%controlled_names(j) = mesh%zones(self%controlled(j))%name
      end do
      self%constraint = constraint
      ! Observation gates include both inlet and outlet. By default only
      ! controlled outlets contribute to a prescribed-flow residual.
      self%coefficients = 0
      has_inlet = .false.
      has_outlet = .false.
      do j = 1, size(observed)
         bc_id = BoundaryId(mesh, observed(j))
         select case (trim(BCs(bc_id)%bc%BCType))
         case ("inflow", "totalinflow")
            has_inlet = .true.
         case ("outflow")
            has_outlet = .true.
         end select
         if (any(self%controlled == observed(j))) self%coefficients(j) = 1
      end do
      if (.not. has_inlet .or. .not. has_outlet) call Fatal("Observed boundaries must contain an inlet and an outlet")
      do j = 1, size(self%controlled)
         if (.not. any(observed == self%controlled(j))) &
            call Fatal("Every controlled outlet must participate in detection")
      end do
      if (options%containsKey("flow coefficients")) then
         value = options%stringValueForKey("flow coefficients", LINE_LENGTH)
         do j = 1, len_trim(value)
            if (value(j:j) == '[' .or. value(j:j) == ']') value(j:j) = ' '
         end do
         read (value, *, iostat=io) self%coefficients
         if (io /= 0) call Fatal("Invalid flow coefficients")
      end if
      select case (trim(constraint))
      case ("target-mass-flow")
         if (.not. options%containsKey("target mass flow (kg/s)")) call Fatal("Missing target mass flow (kg/s)")
         self%target = options%doublePrecisionValueForKey("target mass flow (kg/s)")
         if (self%target <= 0) call Fatal("Target flow must be positive; zero-flow control is not supported")
      case ("equal-mass-flow")
         if (.not. options%containsKey("flow coefficients")) &
            call Fatal("Equal-flow control requires explicit signed coefficients")
         self%target = 0
      case default
         call Fatal("Unsupported constraint: "//trim(constraint))
      end select
      self%pmin = options%doublePrecisionValueForKey("minimum pressure (pa)")
      self%pmax = options%doublePrecisionValueForKey("maximum pressure (pa)")
      self%probe = options%doublePrecisionValueForKey("probe pressure step (pa)")
      self%max_step = options%doublePrecisionValueForKey("maximum pressure step (pa)")
      self%block_duration = options%doublePrecisionValueForKey("average block duration (s)")
      self%relaxation = options%getValueOrDefault("relaxation", 0.5_RP)
      self%rtol = options%getValueOrDefault("relative mass tolerance", 1e-4_RP)
      self%atol = options%getValueOrDefault("absolute mass tolerance (kg/s)", 0.0_RP)
      self%flow_floor = options%getValueOrDefault("flow floor (kg/s)", 1e-12_RP)
      self%flow_scale = options%getValueOrDefault("reference mass flow (kg/s)", abs(self%target))
      self%mean_fraction = options%getValueOrDefault("mean error fraction", 0.25_RP)
      self%minimum_steps = options%getValueOrDefault("minimum settle steps", 50)
      if (.not. all(ieee_is_finite([self%pmin, self%pmax, self%probe, self%max_step, self%target, self%rtol, &
             self%atol, self%flow_scale, self%flow_floor, self%mean_fraction, self%block_duration, self%relaxation])) .or. &
          .not. all(ieee_is_finite(self%coefficients))) call Fatal("Nonfinite configuration")
      if (self%pmin <= 0 .or. self%pmax <= self%pmin .or. self%probe <= 0 .or. self%max_step <= 0 .or. &
          self%block_duration <= 0 .or. self%relaxation <= 0 .or. self%relaxation > 1 .or. self%rtol < 0 .or. &
          self%atol < 0 .or. self%rtol + self%atol <= 0 .or. self%flow_scale < 0 .or. self%flow_floor <= 0 .or. &
          self%mean_fraction <= 0 .or. self%mean_fraction >= 1 .or. self%minimum_steps < 1 .or. &
          all(self%coefficients == 0)) call Fatal("Invalid pressure, tolerance or iteration limits")
      do j = 1, size(self%controlled)
         bc_id = BoundaryId(mesh, self%controlled(j))
         select type (bc => BCs(bc_id)%bc)
         type is (OutflowBC_t)
            initial_pressure = bc%GetPressurePa()
         class default
            call Fatal("Controlled boundary must be outflow")
         end select
         if (.not. ieee_is_finite(initial_pressure)) call Fatal("Nonfinite initial pressure")
         if (j == 1) self%pressure = initial_pressure
         if (abs(initial_pressure - self%pressure) > 1e-10_RP*self%pressure) &
            call Fatal("Controlled patches must share one initial pressure")
      end do
      do j = 1, size(observed)
         area = mesh%zones(observed(j))%no_of_faces
         global_area = area
#ifdef _HAS_MPI_
         if (MPI_Process%doMPIAction) call MPI_Allreduce(area, global_area, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
         if (global_area <= 0) call Fatal("Observed boundary has no global faces")
      end do
      ! Store BC IDs after name resolution; monitor zone IDs remain separate.
      do j = 1, size(self%controlled)
         self%controlled(j) = BoundaryId(mesh, self%controlled(j))
      end do
      allocate (self%pressures(16), self%residuals(16))
      if (MPI_Process%isRoot) allocate (self%statistics(size(observed)))
      self%pressures = 0
      self%residuals = 0
      self%state = WAIT_STABLE
      self%phase_start = iteration
      self%reason = "waiting-for-stability"
      restart_policy = "resume"
      call readValueInRegion(trim(controlFileName), "restart policy", value, REGION, "#end")
      if (len_trim(value) > 0) restart_policy = trim(value)
      if (trim(restart_policy) /= "resume" .and. trim(restart_policy) /= "reset") call Fatal("Invalid restart policy")
      self%log_name = trim(solution_name)//".massflow-control"
      if (cfg%logicalValueForKey("restart")) then
         if (trim(restart_policy) == "resume") then
            value = cfg%stringValueForKey("restart file name", LINE_LENGTH)
            call Restore(self, trim(value))
         end if
      end if
      if (self%pressure < self%pmin .or. self%pressure > self%pmax) call Fatal("Pressure outside configured limits")
      if (cfg%logicalValueForKey("restart")) call ApplyPressure(self)
      call massFlow%Reset(self%monitored, iteration, time)
      if (MPI_Process%isRoot) then
         open (newunit=self%log_unit, file=self%log_name, status="unknown", position="append", action="write", iostat=io)
         if (io /= 0) call Fatal("Cannot open control log")
         write (self%log_unit, '(A)') &
            "# iteration epoch evaluation pressure(Pa) R(kg/s) tolerance drift std state method", &
            "# control-start / mean-start: iteration epoch pressure(Pa)", &
            "# pressure-update: iteration epoch old-pressure(Pa) new-pressure(Pa) method", &
            "# pressure-restricted: iteration epoch pressure(Pa) reason"
         write (self%log_unit, '(A,2(1X,I0),1X,ES24.16)') "control-start", iteration, self%epoch, self%pressure
         flush (self%log_unit)
         print *, "Mass-flow balance enabled; pressure (Pa): ", self%pressure
      end if
      call options%destruct()
   end subroutine BalanceConstruct

   subroutine ResolveZones(mesh, names, zones)
      class(HexMesh), intent(in) :: mesh
      character(len=*), intent(in) :: names
      integer, allocatable, intent(out) :: zones(:)
      character(len=LINE_LENGTH) :: text, name
      integer :: j, count, start, stop, z
      text = trim(names)
      do j = 1, len_trim(text)
         if (text(j:j) == '[' .or. text(j:j) == ']') text(j:j) = ' '
      end do
      count = 1
      do j = 1, len_trim(text)
         if (text(j:j) == ',') count = count + 1
      end do
      allocate (zones(count))
      start = 1
      do j = 1, count
         stop = index(text(start:), ',')
         if (stop == 0) stop = len_trim(text) - start + 2
         name = trim(adjustl(text(start:start + stop - 2)))
         zones(j) = 0
         do z = 1, size(mesh%zones)
            if (trim(mesh%zones(z)%name) == trim(name)) zones(j) = z
         end do
         if (zones(j) == 0) call Fatal("Unknown boundary: "//trim(name))
         if (any(zones(:j - 1) == zones(j))) call Fatal("Duplicate boundary: "//trim(name))
         start = start + stop
      end do
   end subroutine ResolveZones

   integer function BoundaryId(mesh, zone) result(id)
      class(HexMesh), intent(in) :: mesh
      integer, intent(in) :: zone
      integer :: j
      id = 0
      do j = 1, size(BCs)
         if (trim(BCs(j)%bc%bname) == trim(mesh%zones(zone)%name)) id = j
      end do
      if (id == 0) call Fatal("Boundary has no BC object")
   end function BoundaryId

   subroutine BalanceAfterStep(self, massFlow, iteration, time)
      class(MassFlowBalanceController_t), intent(inout) :: self
      type(MassFlowManager_t), intent(inout) :: massFlow
      integer, intent(in) :: iteration
      real(kind=RP), intent(in) :: time
      integer :: action
      logical :: changed, confirming
      real(kind=RP) :: error, tolerance, drift, fluctuation, throughput
#ifdef _HAS_MPI_
      integer :: ierr, packet(3)
#endif
      if (.not. self%enabled) return
      if (self%state == BALANCED .or. self%state == ENDED) return
      action = NO_ACTION
      if (MPI_Process%isRoot) then
         associate(statistics => self%statistics)
            call massFlow%GetStatistics(self%monitored, statistics)
            if (self%state == WAIT_STABLE) then
               if (all(statistics%stable) .and. iteration - self%phase_start >= self%minimum_steps) then
                  self%state = COLLECT_MEAN
                  action = START_MEAN
               end if
            else if (all(statistics%block_count >= 2) .and. statistics(1)%block_count > self%last_block) then
               self%last_block = statistics(1)%block_count
               call EvaluateResidual(self, statistics, error, tolerance, drift, fluctuation, throughput)
               if (.not. all(ieee_is_finite([error, tolerance, drift, fluctuation, throughput]))) &
                  call Fatal("Nonfinite flow response")
               if (.not. ResponseIsUsable(self, statistics, tolerance, drift, fluctuation)) then
                  self%reason = "waiting-for-mean-quality"
               else if (throughput <= self%flow_floor) then
                  self%state = COLLECT_MEAN
                  call RestrictPressure(self, iteration, "no-throughflow")
               else
                  confirming = self%state == CONFIRM
                  call RecordResponse(self, error)
                  if (abs(error) <= tolerance) then
                     if (confirming) then
                        self%state = BALANCED
                        self%reason = "balanced"
                     else
                        self%state = CONFIRM
                        self%reason = "confirming"
                        action = START_MEAN
                     end if
                     call LogResponse(self, iteration, error, tolerance, drift, fluctuation)
                  else
                     ! A failed confirmation is already a usable pressure response.
                     self%state = COLLECT_MEAN
                     call LogResponse(self, iteration, error, tolerance, drift, fluctuation)
                     call ProposePressure(self, tolerance, iteration, changed)
                     if (changed) then
                        self%epoch = self%epoch + 1
                        self%state = WAIT_STABLE
                        self%phase_start = iteration
                        action = UPDATE_PRESSURE
                     end if
                  end if
               end if
            end if
         end associate
      end if
#ifdef _HAS_MPI_
      if (MPI_Process%doMPIAction) then
         packet = [action, self%state, self%epoch]
         call MPI_Bcast(packet, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         action = packet(1)
         self%state = packet(2)
         self%epoch = packet(3)
         call MPI_Bcast(self%pressure, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(self%reason, len(self%reason), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      end if
#endif
      if (action /= NO_ACTION) self%last_block = 0
      select case (action)
      case (START_MEAN)
         call massFlow%BeginAverage(self%monitored, self%block_duration)
         if (MPI_Process%isRoot) then
            write (self%log_unit, '(A,2(1X,I0),1X,ES24.16)') "mean-start", iteration, self%epoch, self%pressure
            flush (self%log_unit)
         end if
      case (UPDATE_PRESSURE)
         call ApplyPressure(self)
         call massFlow%Reset(self%monitored, iteration, time)
      end select
   end subroutine BalanceAfterStep

   subroutine EnsureCapacity(self, needed)
      class(MassFlowBalanceController_t), intent(inout) :: self
      integer, intent(in) :: needed
      real(kind=RP), allocatable :: pressures(:), residuals(:)
      integer :: capacity
      if (needed <= size(self%pressures)) return
      capacity = max(needed, 2*size(self%pressures))
      allocate (pressures(capacity), residuals(capacity))
      pressures = 0; residuals = 0
      pressures(:self%evaluations) = self%pressures(:self%evaluations)
      residuals(:self%evaluations) = self%residuals(:self%evaluations)
      call move_alloc(pressures, self%pressures)
      call move_alloc(residuals, self%residuals)
   end subroutine EnsureCapacity

   subroutine RecordResponse(self, error)
      class(MassFlowBalanceController_t), intent(inout) :: self
      real(kind=RP), intent(in) :: error
      logical :: repeat_pressure
      repeat_pressure = .false.
      if (self%evaluations > 0) repeat_pressure = self%pressures(self%evaluations) == self%pressure
      if (.not. repeat_pressure) then
         call EnsureCapacity(self, self%evaluations + 1)
         self%evaluations = self%evaluations + 1
      end if
      self%pressures(self%evaluations) = self%pressure
      self%residuals(self%evaluations) = error
   end subroutine RecordResponse

   subroutine RestrictPressure(self, iteration, reason)
      class(MassFlowBalanceController_t), intent(inout) :: self
      integer, intent(in) :: iteration
      character(len=*), intent(in) :: reason
      self%reason = reason
      if (self%last_limit_reason == reason) return
      self%last_limit_reason = reason
      write (self%log_unit, '(A,2(1X,I0),1X,ES24.16,1X,A)') &
         "pressure-restricted", iteration, self%epoch, self%pressure, trim(reason)
      flush (self%log_unit)
   end subroutine RestrictPressure

   subroutine EvaluateResidual(self, statistics, error, tolerance, drift, fluctuation, throughput)
      class(MassFlowBalanceController_t), intent(inout) :: self
      type(MassFlowStatistics_t), intent(in) :: statistics(:)
      real(kind=RP), intent(out) :: error, tolerance, drift, fluctuation, throughput
      integer :: j
      real(kind=RP) :: b, weight, quality_weight
      b = 0; drift = 0; fluctuation = 0; throughput = 0
      do j = 1, size(statistics)
         weight = abs(self%coefficients(j))
         b = b + self%coefficients(j)*statistics(j)%blocks(2)%mean
         ! Zero residual coefficients do not exempt an inlet from mean quality
         ! checks. Every observed boundary participates in both gates.
         quality_weight = max(weight, 1.0_RP)
         drift = drift + quality_weight*abs(statistics(j)%blocks(2)%mean - statistics(j)%blocks(1)%mean)
         fluctuation = fluctuation + quality_weight*maxval(statistics(j)%blocks%std)
         throughput = throughput + weight*abs(statistics(j)%blocks(2)%mean)
      end do
      if (self%flow_scale <= 0) self%flow_scale = max(throughput/sum(abs(self%coefficients)), self%flow_floor)
      error = b - self%target
      tolerance = self%atol + self%rtol*max(self%flow_scale, self%flow_floor)
   end subroutine EvaluateResidual

   logical function ResponseIsUsable(self, statistics, tolerance, drift, fluctuation) result(usable)
      class(MassFlowBalanceController_t), intent(in) :: self
      type(MassFlowStatistics_t), intent(in) :: statistics(:)
      real(kind=RP), intent(in) :: tolerance, drift, fluctuation
      integer :: j
      usable = .false.
      if (drift > self%mean_fraction*tolerance .or. fluctuation > self%mean_fraction*tolerance) return
      do j = 1, size(statistics)
         if (any(statistics(j)%blocks%samples < 3)) return
      end do
      usable = .true.
   end function ResponseIsUsable

   subroutine ProposePressure(self, tolerance, iteration, changed)
      class(MassFlowBalanceController_t), intent(inout) :: self
      real(kind=RP), intent(in) :: tolerance
      logical, intent(out) :: changed
      integer, intent(in) :: iteration
      integer :: n, i, j, ia, ib
      real(kind=RP) :: next, old, dp, dr, width, best_width, lo, hi, center
      changed = .false.
      n = self%evaluations
      old = self%pressure
      ia = 0; ib = 0; best_width = huge(1.0_RP)
      do i = 1, n - 1
         do j = i + 1, n
            if ((self%residuals(i) > 0) .eqv. (self%residuals(j) > 0)) cycle
            width = abs(self%pressures(i) - self%pressures(j))
            if (width >= best_width) cycle
            ia = i; ib = j; best_width = width
         end do
      end do
      self%method = "probe"
      if (ia > 0) then
         lo = min(self%pressures(ia), self%pressures(ib))
         hi = max(self%pressures(ia), self%pressures(ib))
         next = self%pressures(ib) - self%residuals(ib)*(self%pressures(ib) - self%pressures(ia))/ &
                (self%residuals(ib) - self%residuals(ia))
         self%method = "bracket-secant"
         if (next <= lo + 0.1_RP*(hi - lo) .or. next >= hi - 0.1_RP*(hi - lo)) then
            next = lo + 0.5_RP*(hi - lo)
            self%method = "bisection"
         end if
         ! Limit toward the candidate without leaving the retained bracket.
         next = max(lo, min(hi, max(old - self%max_step, min(old + self%max_step, next))))
      else
         j = 0
         do i = n - 1, 1, -1
            dp = old - self%pressures(i)
            dr = self%residuals(n) - self%residuals(i)
            if (abs(dp) <= 64*epsilon(old)*max(abs(old), 1.0_RP)) cycle
            if (abs(dr) <= 2*self%mean_fraction*tolerance) cycle
            j = i
            exit
         end do
         if (j > 0) then
            next = old - self%relaxation*self%residuals(n)*(old - self%pressures(j))/ &
                   (self%residuals(n) - self%residuals(j))
            self%method = "relaxed-secant"
         else
            ! Bounded expanding exploration, alternating around the first
            ! pressure if increasing pressure has no resolvable response.
            center = self%pressures(1)
            dp = min(self%max_step, self%probe*2.0_RP**min(n - 1, 20))
            if (n == 1) then
               next = old + dp
               if (next > self%pmax) next = old - dp
            else if (mod(n, 2) == 0) then
               next = center - dp
            else
               next = center + dp
            end if
         end if
         next = max(self%pmin, min(self%pmax, max(old - self%max_step, min(old + self%max_step, next))))
      end if
      if (.not. ieee_is_finite(next)) then
         call Fatal("Nonfinite pressure proposal")
      else if (abs(next - old) <= 64*epsilon(old)*max(abs(old), 1.0_RP)) then
         call RestrictPressure(self, iteration, "pressure-limit-or-insensitive")
      else
         ! Keep observing at the current pressure when exploration is restricted.
         if (any(abs(self%pressures(:n) - next) <= 64*epsilon(old)*max(abs(old), 1.0_RP))) then
            call RestrictPressure(self, iteration, "no-new-pressure-point")
            return
         end if
         changed = .true.
         self%last_limit_reason = ""
         self%pressure = next
         self%reason = "waiting-for-stability"
         write (self%log_unit, '(A,2(1X,I0),2(1X,ES24.16),1X,A)') &
            "pressure-update", iteration, self%epoch + 1, old, next, trim(self%method)
         flush (self%log_unit)
      end if
   end subroutine ProposePressure

   subroutine ApplyPressure(self)
      class(MassFlowBalanceController_t), intent(in) :: self
      integer :: j, status
      do j = 1, size(self%controlled)
         select type (bc => BCs(self%controlled(j))%bc)
         type is (OutflowBC_t)
            call bc%SetPressurePa(self%pressure, status)
            if (status /= 0) call Fatal("Pressure setter rejected proposal")
         end select
      end do
   end subroutine ApplyPressure

   subroutine LogResponse(self, iteration, error, tolerance, drift, fluctuation)
      class(MassFlowBalanceController_t), intent(in) :: self
      integer, intent(in) :: iteration
      real(kind=RP), intent(in) :: error, tolerance, drift, fluctuation
      write (self%log_unit, '(3(I0,1X),5(ES24.16,1X),I0,1X,A)') iteration, self%epoch, self%evaluations, &
         self%pressures(self%evaluations), error, tolerance, drift, fluctuation, self%state, trim(self%method)
      flush (self%log_unit)
   end subroutine LogResponse

   subroutine BalanceFinish(self)
      class(MassFlowBalanceController_t), intent(inout) :: self
      if (.not. self%enabled) return
      if (self%state == BALANCED) return
      self%state = ENDED
      self%reason = "calculation-ended-unbalanced"
   end subroutine BalanceFinish

   subroutine BalanceClose(self)
      class(MassFlowBalanceController_t), intent(inout) :: self
      if (.not. self%enabled) return
      if (MPI_Process%isRoot) then
         write (self%log_unit, '(A,1X,A,1X,ES24.16)') "final", trim(self%reason), self%pressure
         close (self%log_unit)
         print *, "MASS FLOW BALANCE: ", trim(self%reason), "; pressure (Pa) = ", self%pressure
      end if
   end subroutine BalanceClose

   integer(int64) function NamesBytes(names) result(bytes)
      character(len=*), intent(in) :: names(:)
      integer(int64) :: item_bytes
      integer :: j
      inquire (iolength=bytes) size(names)
      do j = 1, size(names)
         inquire (iolength=item_bytes) len_trim(names(j)), trim(names(j))
         bytes = bytes + item_bytes
      end do
   end function NamesBytes

   subroutine WriteNames(u, names)
      integer, intent(in) :: u
      character(len=*), intent(in) :: names(:)
      integer :: j, io
      write (u, iostat=io) size(names)
      if (io /= 0) call Fatal("Cannot save control boundary count")
      do j = 1, size(names)
         write (u, iostat=io) len_trim(names(j)), trim(names(j))
         if (io /= 0) call Fatal("Cannot save control boundary name")
      end do
   end subroutine WriteNames

   subroutine CheckNames(u, names, record_end)
      integer, intent(in) :: u
      character(len=*), intent(in) :: names(:)
      integer(int64), intent(in) :: record_end
      character(len=LINE_LENGTH) :: name
      integer(int64) :: position
      integer :: n, length, j, io
      read (u, iostat=io) n
      if (io /= 0) call Fatal("Cannot read control boundary count")
      if (n /= size(names)) call Fatal("Control boundary count mismatch")
      do j = 1, n
         read (u, iostat=io) length
         if (io /= 0) call Fatal("Cannot read control boundary name length")
         inquire (unit=u, pos=position)
         if (length < 1 .or. length > len(name) .or. int(length, int64) > record_end-position) &
            call Fatal("Invalid control boundary name length")
         read (u, iostat=io) name(:length)
         if (io /= 0) call Fatal("Cannot read control boundary name")
         if (name(:length) /= trim(names(j))) call Fatal("Control boundary name/order mismatch")
      end do
   end subroutine CheckNames

   subroutine BalanceSave(self, u)
      class(MassFlowBalanceController_t), intent(in) :: self
      integer, intent(in) :: u
      integer :: io
      integer(int64) :: position, bytes, record_end
      if (.not. self%enabled) return
      inquire (unit=u, pos=position)
      inquire (iolength=bytes) RECORD_TAG, RECORD_VERSION, record_end, self%constraint, self%target, &
         self%coefficients, self%pressure, self%flow_scale, self%epoch, self%evaluations, &
         self%pressures(:self%evaluations), self%residuals(:self%evaluations)
      record_end = position + bytes + NamesBytes(self%observed_names) + NamesBytes(self%controlled_names)
      write (u, iostat=io) RECORD_TAG, RECORD_VERSION, record_end, self%constraint, self%target
      if (io /= 0) call Fatal("Cannot save embedded control header")
      call WriteNames(u, self%observed_names)
      call WriteNames(u, self%controlled_names)
      write (u, iostat=io) self%coefficients, self%pressure, self%flow_scale, self%epoch, self%evaluations, &
         self%pressures(:self%evaluations), self%residuals(:self%evaluations)
      if (io /= 0) call Fatal("Cannot save embedded control history")
   end subroutine BalanceSave

   subroutine Restore(self, solution)
      class(MassFlowBalanceController_t), intent(inout) :: self
      character(len=*), intent(in) :: solution
      character(len=16) :: tag
      character(len=32) :: constraint
      character(len=LINE_LENGTH) :: boundary_name
      integer :: u, io, version, count, j
      integer(int64) :: start, finish, position, next, index_bytes, record_end, history_bytes
      real(kind=RP) :: target, coefficients(size(self%coefficients))
#ifdef _HAS_MPI_
      integer :: ierr
#endif
      if (MPI_Process%isRoot) then
         open (newunit=u, file=trim(solution), status="old", access="stream", form="unformatted", &
               action="read", iostat=io)
         if (io /= 0) call Fatal("Cannot open control restart solution")
         call GetSolutionExtensionBounds(u, start, finish, io)
         if (io /= 0) call Fatal("Invalid solution extension")
         read (u, pos=start, iostat=io) count
         if (io /= 0) call Fatal("Cannot read embedded boundary count")
         inquire (iolength=index_bytes) next, boundary_name
         if (count < 0 .or. int(count, int64) > (finish-start)/index_bytes) &
            call Fatal("Invalid embedded boundary count")
         do j = 1, count
            read (u, iostat=io) next, boundary_name
            if (io /= 0) call Fatal("Cannot read embedded boundary index")
            inquire (unit=u, pos=position)
            if (next <= position .or. next > finish) call Fatal("Invalid embedded boundary record size")
            read (u, pos=next, iostat=io)
            if (io /= 0) call Fatal("Cannot seek embedded boundary record")
         end do
         inquire (unit=u, pos=position)
         if (position >= finish) call Fatal("Missing embedded control record; use restart policy = reset")
         read (u, iostat=io) tag, version, record_end
         if (io /= 0) call Fatal("Cannot read embedded control header")
         if (tag /= RECORD_TAG .or. version /= RECORD_VERSION .or. record_end /= finish) &
            call Fatal("Invalid embedded control record")
         read (u, iostat=io) constraint, target
         if (io /= 0) call Fatal("Cannot read embedded control configuration")
         if (.not. ieee_is_finite(target)) call Fatal("Nonfinite control target")
         if (constraint /= self%constraint .or. target /= self%target) call Fatal("Control constraint/target mismatch")
         call CheckNames(u, self%observed_names, record_end)
         call CheckNames(u, self%controlled_names, record_end)
         read (u, iostat=io) coefficients, self%pressure, self%flow_scale, self%epoch, self%evaluations
         if (io /= 0) call Fatal("Cannot read embedded control state")
         if (.not. all(ieee_is_finite(coefficients))) call Fatal("Nonfinite control coefficients")
         if (any(coefficients /= self%coefficients)) call Fatal("Control coefficients mismatch")
         if (self%evaluations < 0 .or. self%epoch < 0) call Fatal("Invalid controller history")
         inquire (unit=u, pos=position)
         history_bytes = 2_int64*int(self%evaluations, int64)*(storage_size(self%pressure)/8)
         if (record_end-position /= history_bytes) call Fatal("Invalid embedded control history length")
         deallocate (self%pressures, self%residuals)
         allocate (self%pressures(max(16, self%evaluations)), self%residuals(max(16, self%evaluations)))
         self%pressures = 0; self%residuals = 0
         read (u, iostat=io) self%pressures(:self%evaluations), self%residuals(:self%evaluations)
         if (io /= 0) call Fatal("Cannot read embedded control history")
         close (u)
         if (.not. all(ieee_is_finite([self%pressure, self%flow_scale, &
               self%pressures(:self%evaluations), self%residuals(:self%evaluations)]))) &
            call Fatal("Nonfinite controller history")
         if (self%flow_scale < 0) call Fatal("Invalid restored flow scale")
         if (any(self%pressures(:self%evaluations) < self%pmin) .or. &
             any(self%pressures(:self%evaluations) > self%pmax)) call Fatal("History pressure outside limits")
         self%reason = "restart-redetect"
      end if
#ifdef _HAS_MPI_
      if (MPI_Process%doMPIAction) then
         call MPI_Bcast(self%pressure, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(self%epoch, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      end if
#endif
      ! Retain completed responses; detection and partial means restart at restored p.
   end subroutine Restore

   subroutine Fatal(message)
      character(len=*), intent(in) :: message
#ifdef _HAS_MPI_
      integer :: ierr
#endif
      print *, "MASS FLOW BALANCE: ", message
#ifdef _HAS_MPI_
      if (MPI_Process%doMPIAction) call MPI_Abort(MPI_COMM_WORLD, 99, ierr)
#endif
      error stop 99
   end subroutine Fatal
end module MassFlowBalanceController

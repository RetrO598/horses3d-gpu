#include "Includes.h"
module TotalInflowBCClass
#if defined(NAVIERSTOKES) && !defined(CAHNHILLIARD)
   use SMConstants
   use PhysicsStorage
   use FileReaders,            only: controlFileName
   use FileReadingUtilities,   only: GetKeyword, GetValueAsString, PreprocessInputLine, CheckIfBoundaryNameIsContained
   use FTValueDictionaryClass, only: FTValueDictionary
   use Utilities,             only: toLower
   use GenericBoundaryConditionClass
   use FluidData
   use HexMeshClass
   use ZoneClass
   implicit none
!
!  *****************************
!  Default everything to private
!  *****************************
!
   private
!
!  ******************
!  Public definitions
!  ******************
!
   public TotalInflowBC_t, ComputeTotalInletState
!
!  ****************
!  Class definition
!  ****************
!
   type, extends(GenericBC_t) :: TotalInflowBC_t
      real(kind=RP)              :: p0
      real(kind=RP)              :: T0
      real(kind=RP)              :: direction(NDIM)
#if defined(SPALARTALMARAS)
      real(kind=RP)              :: eddy_theta
#endif
      contains
         procedure         :: Destruct          => TotalInflowBC_Destruct
         procedure         :: Describe          => TotalInflowBC_Describe
         procedure         :: FlowState         => TotalInflowBC_FlowState
         procedure         :: FlowNeumann       => TotalInflowBC_FlowNeumann
         procedure         :: CreateDeviceData  => TotalInflowBC_CreateDeviceData
         procedure         :: ExitDeviceData    => TotalInflowBC_ExitDeviceData
   end type TotalInflowBC_t
!
!  *******************************************************************
!  Traditionally, constructors are exported with the name of the class
!  *******************************************************************
!
   interface TotalInflowBC_t
      module procedure ConstructTotalInflowBC
   end interface TotalInflowBC_t
!
!  ========
   contains
!  ========
!
!/////////////////////////////////////////////////////////
!
!        Class constructor
!        -----------------
!
!/////////////////////////////////////////////////////////
!
      function ConstructTotalInflowBC(bname)
!
!        ********************************************************************
!        Definition of the total-condition inlet in the control file:
!              #define boundary bname
!                 type              = totalinflow
!                 total pressure    = #value        (Pa, absolute)
!                 total temperature = #value        (K)
!                 AoAPhi            = #value        (degrees)
!                 AoATheta          = #value        (degrees)
!                 Turbulence parameter theta = #value   (SA, same scaling as Inflow)
!              #end
!        ********************************************************************
!
         implicit none
         type(TotalInflowBC_t)        :: ConstructTotalInflowBC
         character(len=*), intent(in) :: bname
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                    :: fid, io
         character(len=LINE_LENGTH) :: currentLine
         character(len=LINE_LENGTH) :: keyword, keyval
         logical                    :: inside
         real(kind=RP)              :: AoAPhi, AoATheta
         type(FTValueDictionary)    :: bcdict

         if ( .not. flowIsNavierStokes ) then
            print*, "TotalInflow is only supported for compressible Navier-Stokes (laminar or SA)."
            errorMessage(STD_OUT)
            error stop 99
         end if

         open(newunit = fid, file = trim(controlFileName), status = "old", action = "read")

         call bcdict % InitWithSize(16)

         ConstructTotalInflowBC % bname  = bname
         ConstructTotalInflowBC % BCType = "totalinflow"
         call toLower(ConstructTotalInflowBC % bname)
!
!        Navigate until the "#define boundary bname" sentinel is found
!        -------------------------------------------------------------
         inside = .false.
         do
            read(fid, '(A)', iostat=io) currentLine

            if ( io .ne. 0 ) exit

            call PreprocessInputLine(currentLine)
            call toLower(currentLine)

            if ( index(trim(currentLine),"#define boundary") .ne. 0 ) then
               inside = CheckIfBoundaryNameIsContained(trim(currentLine), trim(ConstructTotalInflowBC % bname))
            end if
!
!           Get all keywords inside the zone
!           --------------------------------
            if ( inside ) then
               if ( trim(currentLine) .eq. "#end" ) exit

               keyword = adjustl(GetKeyword(currentLine))
               keyval  = adjustl(GetValueAsString(currentLine))
               call bcdict % AddValueForKey(keyval, trim(keyword))
            end if
         end do
!
!        Analyze the gathered data
!        -------------------------
#if defined(SPALARTALMARAS)
         call GetValueWithDefault(bcdict, "Turbulence parameter theta", refValues % mu, ConstructTotalInflowBC % eddy_theta)
         ConstructTotalInflowBC % eddy_theta = ConstructTotalInflowBC % eddy_theta / refValues % mu
#endif
         if ( .not. bcdict % ContainsKey("total pressure") .or. &
              .not. bcdict % ContainsKey("total temperature") ) then
            print*, "Missing total pressure or total temperature for boundary ", trim(bname)
            errorMessage(STD_OUT)
            error stop 99
         end if

         ConstructTotalInflowBC % p0 = bcdict % DoublePrecisionValueForKey("total pressure")
         ConstructTotalInflowBC % T0 = bcdict % DoublePrecisionValueForKey("total temperature")
         call GetValueWithDefault(bcdict, "aoaphi"  , refValues % AoAPhi  , AoAPhi  )
         call GetValueWithDefault(bcdict, "aoatheta", refValues % AoATheta, AoATheta)

         ConstructTotalInflowBC % p0 = ConstructTotalInflowBC % p0 / refValues % p
         ConstructTotalInflowBC % T0 = ConstructTotalInflowBC % T0 / refValues % T
         AoAPhi   = AoAPhi * PI / 180.0_RP
         AoATheta = AoATheta * PI / 180.0_RP

         ConstructTotalInflowBC % direction(IX) = cos(AoATheta) * cos(AoAPhi)
         ConstructTotalInflowBC % direction(IY) = sin(AoATheta) * cos(AoAPhi)
         ConstructTotalInflowBC % direction(IZ) = sin(AoAPhi)
         ConstructTotalInflowBC % constructed = .true.

         call bcdict % Destruct
         close(fid)

      end function ConstructTotalInflowBC
!
!/////////////////////////////////////////////////////////
!
!        Boundary condition description
!        ------------------------------
!
!/////////////////////////////////////////////////////////
!
      subroutine TotalInflowBC_Describe(self)
         implicit none
         class(TotalInflowBC_t), intent(in) :: self

         write(STD_OUT,'(30X,A,A28,A)')    "->", " Boundary condition type: ", "TotalInflow"
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", " Total pressure (Pa): ", self % p0 * refValues % p
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", " Total temperature (K): ", self % T0 * refValues % T
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", " AoaPhi: ", asin(self % direction(IZ)) * 180.0_RP / PI
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", " AoaTheta: ", &
                                          atan2(self % direction(IY), self % direction(IX)) * 180.0_RP / PI

#if defined(SPALARTALMARAS)
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", " SA boundary theta: ", 3.0_RP * self % eddy_theta
#endif

      end subroutine TotalInflowBC_Describe
!
!/////////////////////////////////////////////////////////
!
!        Class destructor and device data
!        --------------------------------
!
!/////////////////////////////////////////////////////////
!
      subroutine TotalInflowBC_Destruct(self)
         implicit none
         class(TotalInflowBC_t) :: self

         if ( .not. self % constructed ) return
         call self % ExitDeviceData()
         self % constructed = .false.

      end subroutine TotalInflowBC_Destruct

      subroutine TotalInflowBC_CreateDeviceData(self)
         implicit none
         class(TotalInflowBC_t), intent(in) :: self

         !$acc enter data copyin(self)

      end subroutine TotalInflowBC_CreateDeviceData

      subroutine TotalInflowBC_ExitDeviceData(self)
         implicit none
         class(TotalInflowBC_t), intent(in) :: self

         !$acc wait(1)
         !$acc exit data delete(self)

      end subroutine TotalInflowBC_ExitDeviceData
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for compressible Navier--Stokes equations
!        -----------------------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
      subroutine TotalInflowBC_FlowState(self, mesh, zone)
         implicit none
         class(TotalInflowBC_t), intent(in) :: self
         type(HexMesh),       intent(inout) :: mesh
         type(Zone_t),           intent(in) :: zone
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: Q(NCONS), Qext(NCONS), nHat(NDIM)
         integer       :: i, j, fID, zonefID

         !$acc parallel loop gang present(mesh, self, zone) private(fID) async(1)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(2) private(Q, Qext, nHat)
            do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1)

               Q    = mesh % faces(fID) % storage(1) % Q(:,i,j)
               nHat = mesh % faces(fID) % geom % normal(:,i,j)

               call ComputeTotalInletState(Q, nHat, self % direction, self % p0, self % T0, &
                                          thermodynamics % gamma, dimensionless % gammaM2, Qext)

#if defined(SPALARTALMARAS)
               Qext(IRHOTHETA) = Qext(IRHO) * 3.0_RP * self % eddy_theta
#endif
               mesh % faces(fID) % storage(2) % Q(:,i,j) = Qext
            end do ; end do
         end do
         !$acc end parallel loop

      end subroutine TotalInflowBC_FlowState

      pure subroutine ComputeTotalInletState(Q, nHat, direction, p0, T0, gamma, gammaM2, Qext)
         !$acc routine seq
         implicit none
         real(kind=RP), intent(in)  :: Q(NCONS), nHat(NDIM), direction(NDIM)
         real(kind=RP), intent(in)  :: p0, T0, gamma, gammaM2
         real(kind=RP), intent(out) :: Qext(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: vel(NDIM), vel2, alpha, rPlus
         real(kind=RP) :: p, a2, a02, gm1
         real(kind=RP) :: aa, bb, cc, dd, Mach2
         real(kind=RP) :: Tb, pb, rho, speed
!
!        Extrapolate the outgoing acoustic invariant from the interior
!        -------------------------------------------------------------
         gm1    = gamma - 1.0_RP
         vel    = Q(IRHOU:IRHOW) / Q(IRHO)
         vel2   = dot_product(vel, vel)
         p      = gm1 * (Q(IRHOE) - 0.5_RP * Q(IRHO) * vel2)
         a2     = gamma * p / Q(IRHO)
         rPlus  = dot_product(vel, nHat) + 2.0_RP * sqrt(a2) / gm1
!
!        Total sound speed and prescribed inlet direction
         a02   = gamma * T0 / gammaM2
         alpha = dot_product(nHat, direction)
!
!        Solve the quadratic equation for the inlet speed
!        ------------------------------------------------
         aa = 1.0_RP + 0.5_RP * gm1 * alpha*alpha
         bb = -gm1 * alpha * rPlus
         cc = 0.5_RP * gm1 * rPlus*rPlus - 2.0_RP * a02 / gm1
         dd = sqrt(max(0.0_RP, bb*bb - 4.0_RP * aa * cc))

         speed = max(0.0_RP, (-bb + dd) / (2.0_RP * aa))
         vel2  = speed*speed
         a2    = a02 - 0.5_RP * gm1 * vel2
!
!        Apply the Mach-number limit and update the sound speed
!        ----------------------------------------------------------
         Mach2 = min(1.0_RP, vel2 / a2)
         vel2  = Mach2 * a2
         speed = sqrt(vel2)
         a2    = a02 - 0.5_RP * gm1 * vel2
!
!        Static boundary values and conservative exterior state
!        ------------------------------------------------------
!        The solver uses T = gammaM2*p/rho.
         Tb    = gammaM2 * a2 / gamma
         pb    = p0 * (Tb / T0)**(gamma / gm1)
         rho   = gammaM2 * pb / Tb

#if defined(SPALARTALMARAS)
!        The caller replaces this extrapolated SA value with the prescribed inlet value.
         Qext(IRHOTHETA) = rho * Q(IRHOTHETA) / Q(IRHO)
#endif
         Qext(IRHO)        = rho
         Qext(IRHOU:IRHOW) = rho * speed * direction
         Qext(IRHOE)       = pb / gm1 + 0.5_RP * rho * vel2

      end subroutine ComputeTotalInletState

      subroutine TotalInflowBC_FlowNeumann(self, mesh, zone)
!
!        *******************************************
!        Cancel out the viscous flux at the inlet
!        *******************************************
!
         implicit none
         class(TotalInflowBC_t), intent(in) :: self
         type(HexMesh),       intent(inout) :: mesh
         type(Zone_t),           intent(in) :: zone
!
!        ---------------
!        Local variables
!        ---------------
!
         integer :: i, j, fID, zonefID

         !$acc parallel loop gang present(mesh, self, zone) private(fID) async(1)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(2)
            do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1)
               mesh % faces(fID) % storage(2) % FStar(:,i,j) = 0.0_RP
            end do ; end do
         end do
         !$acc end parallel loop

      end subroutine TotalInflowBC_FlowNeumann
#else
   implicit none
   private
#endif
end module TotalInflowBCClass

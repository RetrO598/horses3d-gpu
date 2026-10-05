#include "Includes.h"
module FarfieldBCClass
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
   use Physics,                only: ViscousFlux_STATE
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
   public FarfieldBC_t, ComputeFarfieldState
!
!  ****************
!  Class definition
!  ****************
!
   type, extends(GenericBC_t) :: FarfieldBC_t
      real(kind=RP)              :: p
      real(kind=RP)              :: rho
      real(kind=RP)              :: v
      real(kind=RP)              :: AoAPhi
      real(kind=RP)              :: AoATheta
      real(kind=RP)              :: QInf(NCONS)
#if defined(SPALARTALMARAS)
      real(kind=RP)              :: eddy_theta
#endif
      contains
         procedure         :: Destruct          => FarfieldBC_Destruct
         procedure         :: Describe          => FarfieldBC_Describe
         procedure         :: FlowState         => FarfieldBC_FlowState
         procedure         :: FlowGradVars      => FarfieldBC_FlowGradVars
         procedure         :: FlowNeumann       => FarfieldBC_FlowNeumann
         procedure         :: CreateDeviceData  => FarfieldBC_CreateDeviceData
         procedure         :: ExitDeviceData    => FarfieldBC_ExitDeviceData
   end type FarfieldBC_t
!
!  *******************************************************************
!  Traditionally, constructors are exported with the name of the class
!  *******************************************************************
!
   interface FarfieldBC_t
      module procedure ConstructFarfieldBC
   end interface FarfieldBC_t
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
      function ConstructFarfieldBC(bname)
!
!        ********************************************************************
!        Definition of the characteristic farfield in the control file:
!              #define boundary bname
!                 type              = farfield
!                 pressure          = #value        (Pa, absolute static pressure)
!                 density           = #value        (kg/m^3)
!                 mach              = #value
!                 AoAPhi            = #value        (degrees)
!                 AoATheta          = #value        (degrees)
!                 Turbulence parameter theta = #value   (SA, same scaling as Inflow)
!              #end
!        Omitted values use the global reference freestream.
!        This boundary supports stationary meshes and laminar or SA Navier-Stokes.
!        ********************************************************************
!
         implicit none
         type(FarfieldBC_t)        :: ConstructFarfieldBC
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
         real(kind=RP)              :: velocity(NDIM)
         type(FTValueDictionary)    :: bcdict

         if ( .not. flowIsNavierStokes ) then
            print*, "Farfield is only supported for compressible Navier-Stokes (laminar or SA)."
            errorMessage(STD_OUT)
            error stop 99
         end if

         open(newunit = fid, file = trim(controlFileName), status = "old", action = "read")

         call bcdict % InitWithSize(16)

         ConstructFarfieldBC % bname  = bname
         ConstructFarfieldBC % BCType = "farfield"
         call toLower(ConstructFarfieldBC % bname)
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
               inside = CheckIfBoundaryNameIsContained(trim(currentLine), trim(ConstructFarfieldBC % bname))
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
         call GetValueWithDefault(bcdict, "Turbulence parameter theta", refValues % mu, ConstructFarfieldBC % eddy_theta)
         ConstructFarfieldBC % eddy_theta = ConstructFarfieldBC % eddy_theta / refValues % mu
#endif
         call GetValueWithDefault(bcdict, "pressure", refValues % p / dimensionless % gammaM2, ConstructFarfieldBC % p)
         call GetValueWithDefault(bcdict, "density" , refValues % rho                        , ConstructFarfieldBC % rho)
         call GetValueWithDefault(bcdict, "mach"    , dimensionless % Mach                   , ConstructFarfieldBC % v)
         call GetValueWithDefault(bcdict, "aoaphi"  , refValues % AoAPhi                     , ConstructFarfieldBC % AoAPhi)
         call GetValueWithDefault(bcdict, "aoatheta", refValues % AoATheta                   , ConstructFarfieldBC % AoATheta)

         associate ( bc => ConstructFarfieldBC )
         bc % p        = bc % p / refValues % p
         bc % rho      = bc % rho / refValues % rho
         bc % v        = bc % v * sqrt(thermodynamics % gamma * bc % p / bc % rho)
         bc % AoAPhi   = bc % AoAPhi * PI / 180.0_RP
         bc % AoATheta = bc % AoATheta * PI / 180.0_RP

         velocity(IX) = bc % v * cos(bc % AoATheta) * cos(bc % AoAPhi)
         velocity(IY) = bc % v * sin(bc % AoATheta) * cos(bc % AoAPhi)
         velocity(IZ) = bc % v * sin(bc % AoAPhi)
         bc % QInf(IRHO)        = bc % rho
         bc % QInf(IRHOU:IRHOW) = bc % rho * velocity
         bc % QInf(IRHOE)       = bc % p / thermodynamics % gammaMinus1 + 0.5_RP * bc % rho * bc % v**2
#if defined(SPALARTALMARAS)
         bc % QInf(IRHOTHETA) = bc % rho * 3.0_RP * bc % eddy_theta
#endif
         end associate
         ConstructFarfieldBC % constructed = .true.

         call bcdict % Destruct
         close(fid)

      end function ConstructFarfieldBC
!
!/////////////////////////////////////////////////////////
!
!        Boundary condition description
!        ------------------------------
!
!/////////////////////////////////////////////////////////
!
      subroutine FarfieldBC_Describe(self)
         implicit none
         class(FarfieldBC_t), intent(in) :: self

         write(STD_OUT,'(30X,A,A28,A)')    "->", " Boundary condition type: ", "Farfield"
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", " Pressure (Pa): ", self % p * refValues % p
         write(STD_OUT,'(30X,A,A28,F10.4)') "->", " Density (kg/m^3): ", self % rho * refValues % rho
         write(STD_OUT,'(30X,A,A28,F10.4)') "->", " Mach number: ", self % v / sqrt(thermodynamics % gamma * self % p / self % rho)
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", " AoaPhi: ", self % AoAPhi * 180.0_RP / PI
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", " AoaTheta: ", self % AoATheta * 180.0_RP / PI

#if defined(SPALARTALMARAS)
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", " SA boundary theta: ", 3.0_RP * self % eddy_theta
#endif

      end subroutine FarfieldBC_Describe
!
!/////////////////////////////////////////////////////////
!
!        Class destructor and device data
!        --------------------------------
!
!/////////////////////////////////////////////////////////
!
      subroutine FarfieldBC_Destruct(self)
         implicit none
         class(FarfieldBC_t) :: self

         if ( .not. self % constructed ) return
         call self % ExitDeviceData()
         self % constructed = .false.

      end subroutine FarfieldBC_Destruct

      subroutine FarfieldBC_CreateDeviceData(self)
         implicit none
         class(FarfieldBC_t), intent(in) :: self

         !$acc enter data copyin(self)

      end subroutine FarfieldBC_CreateDeviceData

      subroutine FarfieldBC_ExitDeviceData(self)
         implicit none
         class(FarfieldBC_t), intent(in) :: self

         !$acc wait(1)
         !$acc exit data delete(self)

      end subroutine FarfieldBC_ExitDeviceData
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for compressible Navier--Stokes equations
!        -----------------------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
      subroutine FarfieldBC_FlowState(self, mesh, zone)
         implicit none
         class(FarfieldBC_t), intent(in) :: self
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

               call ComputeFarfieldState(Q, self % QInf, nHat, thermodynamics % gamma, Qext)

               mesh % faces(fID) % storage(2) % Q(:,i,j) = Qext
            end do ; end do
         end do
         !$acc end parallel loop

      end subroutine FarfieldBC_FlowState

      pure subroutine ComputeFarfieldState(Q, QInf, nHat, gamma, Qext)
         !$acc routine seq
         implicit none
         real(kind=RP), intent(in)  :: Q(NCONS), QInf(NCONS), nHat(NDIM), gamma
         real(kind=RP), intent(out) :: Qext(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: vel(NDIM), velInf(NDIM), velExt(NDIM)
         real(kind=RP) :: p, pInf, a, aInf, un, unInf, rPlus, rMinus
         real(kind=RP) :: gm1, entropy, rho, unExt, aExt, pExt

         gm1    = gamma - 1.0_RP
         vel    = Q(IRHOU:IRHOW) / Q(IRHO)
         velInf = QInf(IRHOU:IRHOW) / QInf(IRHO)
         p      = gm1 * (Q(IRHOE) - 0.5_RP * Q(IRHO) * dot_product(vel, vel))
         pInf   = gm1 * (QInf(IRHOE) - 0.5_RP * QInf(IRHO) * dot_product(velInf, velInf))
         a      = sqrt(gamma * p / Q(IRHO))
         aInf   = sqrt(gamma * pInf / QInf(IRHO))
         un     = dot_product(vel, nHat)
         unInf  = dot_product(velInf, nHat)
!
!        Select acoustic invariants using the freestream normal velocity
!        ---------------------------------------------------------------
!        The mesh is stationary and nHat points out of the domain.
         if ( unInf > -aInf ) then
            rPlus = un + 2.0_RP * a / gm1
         else
            rPlus = unInf + 2.0_RP * aInf / gm1
         end if

         if ( unInf > aInf ) then
            rMinus = un - 2.0_RP * a / gm1
         else
            rMinus = unInf - 2.0_RP * aInf / gm1
         end if

         unExt = 0.5_RP * (rPlus + rMinus)
         aExt  = 0.25_RP * gm1 * (rPlus - rMinus)
!
!        Entropy and tangential velocity come from the upstream state
!        ------------------------------------------------------------
         if ( unInf > 0.0_RP ) then
            velExt  = vel + (unExt - un) * nHat
            entropy = Q(IRHO)**gamma / p
         else
            velExt  = velInf + (unExt - unInf) * nHat
            entropy = QInf(IRHO)**gamma / pInf
         end if

         rho  = (entropy * aExt*aExt / gamma)**(1.0_RP / gm1)
         pExt = rho * aExt*aExt / gamma
#if defined(SPALARTALMARAS)
         if ( unInf > 0.0_RP ) then
            Qext(IRHOTHETA) = rho * Q(IRHOTHETA) / Q(IRHO)
         else
            Qext(IRHOTHETA) = rho * QInf(IRHOTHETA) / QInf(IRHO)
         end if
#endif
         Qext(IRHO)        = rho
         Qext(IRHOU:IRHOW) = rho * velExt
         Qext(IRHOE)       = pExt / gm1 + 0.5_RP * rho * dot_product(velExt, velExt)

      end subroutine ComputeFarfieldState

      subroutine FarfieldBC_FlowGradVars(self, mesh, zone)
         implicit none
         class(FarfieldBC_t), intent(in) :: self
         type(HexMesh),    intent(inout) :: mesh
         type(Zone_t),        intent(in) :: zone
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: i, j, d, fID, zonefID
         real(kind=RP) :: jump(NGRAD)
!
!        Average interior and exterior states, as in GenericBC_FlowGradVars
!        -----------------------------------------------------------------
         !$acc parallel loop gang present(mesh, zone) private(fID) async(1)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(2) private(jump)
            do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1)
               jump = 0.5_RP * (mesh % faces(fID) % storage(2) % Q(:,i,j) - &
                               mesh % faces(fID) % storage(1) % Q(:,i,j))
               !$acc loop seq
               do d = 1, NDIM
                  mesh % faces(fID) % storage(1) % unStar(:,d,i,j) = jump * &
                         mesh % faces(fID) % geom % normal(d,i,j) * mesh % faces(fID) % geom % jacobian(i,j)
               end do
            end do ; end do
         end do
         !$acc end parallel loop

      end subroutine FarfieldBC_FlowGradVars

      subroutine FarfieldBC_FlowNeumann(self, mesh, zone)
         implicit none
         class(FarfieldBC_t), intent(in) :: self
         type(HexMesh),    intent(inout) :: mesh
         type(Zone_t),        intent(in) :: zone
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: i, j, fID, zonefID
         real(kind=RP) :: flux(NCONS,NDIM), Q(NCONS), Qext(NCONS), velocityCorrection(NDIM)
!
!        Extrapolate interior primitive gradients and transport coefficients
!        -------------------------------------------------------------------
!        Stress and heat conduction use interior gradients. Viscous work uses
!        the mean boundary velocity. This adapts gradient extrapolation to DG;
!        the finite-volume reflected-point correction is not used.
         !$acc parallel loop gang present(mesh, zone) private(fID) async(1)
         do zonefID = 1, zone % no_of_faces
            fID = zone % faces(zonefID)
            !$acc loop vector collapse(2) private(flux, Q, Qext, velocityCorrection)
            do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1)
               Q    = mesh % faces(fID) % storage(1) % Q(:,i,j)
               Qext = mesh % faces(fID) % storage(2) % Q(:,i,j)
               call ViscousFlux_STATE(NCONS, NGRAD, Q, &
                                     mesh % faces(fID) % storage(1) % U_x(:,i,j), &
                                     mesh % faces(fID) % storage(1) % U_y(:,i,j), &
                                     mesh % faces(fID) % storage(1) % U_z(:,i,j), &
                                     mesh % faces(fID) % storage(1) % mu_NS(1,i,j), 0.0_RP, &
                                     mesh % faces(fID) % storage(1) % mu_NS(2,i,j), flux)

#if defined(SPALARTALMARAS)
!              Cancel the SA diffusive boundary flux, as in InflowBC.
               flux(IRHOTHETA,:) = 0.0_RP
#endif
               velocityCorrection = 0.5_RP * (Qext(IRHOU:IRHOW) / Qext(IRHO) - Q(IRHOU:IRHOW) / Q(IRHO))
               flux(IRHOE,:) = flux(IRHOE,:) + matmul(velocityCorrection, flux(IRHOU:IRHOW,:))
               mesh % faces(fID) % storage(2) % FStar(:,i,j) = matmul(flux, mesh % faces(fID) % geom % normal(:,i,j))
            end do ; end do
         end do
         !$acc end parallel loop

      end subroutine FarfieldBC_FlowNeumann
#else
   implicit none
   private
#endif
end module FarfieldBCClass

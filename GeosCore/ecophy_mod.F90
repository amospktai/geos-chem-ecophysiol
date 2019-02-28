!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: ecophy_mod.F90
!
! !DESCRIPTION: Module ECOPHY\_MOD contains variables and routines for the
!  GEOS-Chem ecophysiology scheme. It modifies the bulk canopy stomatal 
!  resistance RIX in the dry deposition scheme.
!\\
!\\
! !INTERFACE: 
!
      MODULE ECOPHY_MOD
! 
! !USES:
! 
      USE CMN_SIZE_MOD                          ! Size parameters
      USE ERROR_MOD                             ! Error handling routines
      USE PhysConstants                         ! Physical constants
      USE PRECISION_MOD                         ! For GEOS-Chem Precision (fp)
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: DO_ECOPHY
      PUBLIC :: INIT_ECOPHY
      PUBLIC :: CLEANUP_ECOPHY
!
! !REMARKS:
!  References:
!  ============================================================================
!EOP
!------------------------------------------------------------------------------
!BOC
!
!DEFINED PARAMETERS:
!
      !---------------------------------------------------------------------------------------
      ! NUMPFT        : Total number of PFTs                              []
      ! IS_C3_PLANT   : IS_C3_PLANT = 1 for C3 plants, else 0             []
      ! V_CMAX25      : V_CMAX at 25 deg C                                [mol CO2 m^-2 s^-1]
      ! ALPHA         : Quantum efficiency of photosynthesis              [mol CO2 mol^-1 PAR]
      ! T_UPP         : PFT-specific parameter for V_CMAX                 [deg C]
      ! T_LOW         : PFT-specific parameter for V_CMAX                 [deg C]
      ! F_DARKRESP    : Dark respiration coefficient                      []
      ! D_STAR        : PFT-specific parameter for closure eq.            [kg H2O / kg air]
      ! f0            : PFT-specific parameter for closure eq.            []
      ! G_LEAF_MIN    : PFT-specific min leaf conductance for closure eq. [m s^-1]
      ! K_EXTINCT     : Light extinction coefficient                      []
      ! PARAM_A       : PFT-specific parameter for ozone damage scheme    [m^2 s nmol^-1]
      ! FLUXO3_CRIT   : PFT-specific threshold for ozone uptake           [nmol m^-2 s^-1]
      ! THETA_WILT    : Volumetric soil moisture at wilting point         [m^3 water / m^3 soil]
      ! THETA_CRIT    : Critical value of volumetric soil moisture        [m^3 water / m^3 soil]
      !                 at which photosynthesis is not limited by it
      ! THRESHOLD     : Threshold of relative error                       []
      !---------------------------------------------------------------------------------------
      INTEGER,  PARAMETER   :: NUMPFT                 = 5 ! Switch to call CMN_SIZE_MOD.F later
      INTEGER,  PARAMETER   :: IS_C3_PLANT   (NUMPFT) = (/ 1,1,1,0,1 /)
      REAL,     PARAMETER   :: ALPHA         (NUMPFT) = (/ 0.08, 0.08, 0.12, 0.06, 0.08 /)
      REAL,     PARAMETER   :: V_CMAX25      (NUMPFT) = (/ 0.046, 0.033, 0.073, 0.060, 0.060 /) * &
                                                        (/ 0.0008, 0.0008, 0.0008, 0.0004, 0.0008 /)
      REAL,     PARAMETER   :: T_UPP         (NUMPFT) = (/ 36.0, 26.0, 36.0, 45.0, 36.0 /)
      REAL,     PARAMETER   :: T_LOW         (NUMPFT) = (/ 0.0, -10.0, 0.0,  13.0, 0.0  /)
      REAL,     PARAMETER   :: F_DARKRESP    (NUMPFT) = (/ .015, .015, .015, .025, .015 /)
      REAL,     PARAMETER   :: D_STAR        (NUMPFT) = (/ 0.09, 0.06, 0.1,  .075, 0.1  /)
      REAL,     PARAMETER   :: f0            (NUMPFT) = (/ .875, .875, 0.9,  0.8,  0.9  /)
      REAL,     PARAMETER   :: G_LEAF_MIN    (NUMPFT) = 1.0e-6
      REAL,     PARAMETER   :: K_EXTINCT     (NUMPFT) = 0.5
      REAL,     PARAMETER   :: PARAM_A       (NUMPFT) = (/ 0.04, 0.02, 0.25, 0.13, 0.03 /)
      REAL,     PARAMETER   :: FLUXO3_CRIT   (NUMPFT) = (/ 1.6,  1.6,  5.0,  5.0,  1.6  /)
      REAL,     PARAMETER   :: THRESHOLD              = 1.0e-3
      ! Constants
      ! REAL,     PARAMETER   :: RSTARG = 8.31446        ! Switch to call physconstant.F later (in Headers) 
      REAL,     PARAMETER   :: CO2_O2_RATIO = 1.6
!
! PRIVATE TYPES:
!
      !========================================================================
      ! MODULE VARIABLES:
      ! THETA_SATU    : Soil moisture at saturation point                 [mm]
      ! THETA_CRIT    : Soil moisture at critical point                   [mm]
      ! THETA_WILT    : Soil moisture at wilting point                    [mm]
      ! LO3_DAMAGE    : Logical switch for ozone damage scheme            []
      ! NUMPFT        : Total number of PFTs                              []
      !========================================================================
      REAL,     ALLOCATABLE :: THETA_SATU    ( :,: )    
      REAL,     ALLOCATABLE :: THETA_CRIT    ( :,: )    
      REAL,     ALLOCATABLE :: THETA_WILT    ( :,: )  
      INTEGER,  ALLOCATABLE :: IPFT          ( :   )
      LOGICAL               :: LECOPHY
      LOGICAL               :: LO3_DAMAGE
      ! INTEGER               :: NUMPFT

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_ecophy
!
! !DESCRIPTION: Subroutine DO\_ECOPHY is the interface between dry 
!  deposition module and the ecophysiology module. It computes the
!  bulk canopy stomatal resistance r_s according to meterological inputs, 
!  plant functional types, and soil types.
!\\
!\\
! !INTERFACE:
!      
      SUBROUTINE DO_ECOPHY ( am_I_Root, Input_Opt,  State_Met, &
                             State_Chm, State_Diag, RC,        & 
                             I, J,      LDT, RS                )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE Species_Mod,        ONLY : Species
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE State_Diag_Mod,     ONLY : DgnState
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
      INTEGER,        INTENT(IN)    :: I           ! longitude index
      INTEGER,        INTENT(IN)    :: J           ! latitude index
      INTEGER,        INTENT(IN)    :: PFT         ! PFT index
      INTEGER,        INTENT(IN)    :: LDT         ! Land type index
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
      REAL,           INTENT(OUT)   :: RS          ! Bulk canopy stomatal resistance
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      !  Note: Read subroutine DO_PHOTOSYNTHESIS for descriptions. 
      REAL       :: TEMPK         
      REAL       :: SPHU          
      REAL       :: RA            
      REAL       :: PAR_ABSORBED  
      REAL       :: PRESSURE      
      REAL       :: CO2           
      REAL       :: O2            
      REAL       :: O3            
      REAL       :: LAI           
      LOGICAL    :: LO3_DAMAGE    
      REAL       :: SOIL_WETNESS  
      ! INTEGER    :: PFT     
      REAL       :: G_CAN_OUT     
      REAL       :: G_LEAF_OUT    
      REAL       :: CO2_IN        
      REAL       :: A_CAN_OUT    
      REAL       :: A_NET_OUT     
      REAL       :: RESP_CAN_OUT  
      REAL       :: RESP_OUT      
      REAL       :: FLUXO3_CAN    
      REAL       :: FLUXO3        
      REAL       :: FACTOR_O3     
      REAL       :: BETA       
      ! Arrays

      ! Pointers
      ! REAL(fp), POINTER :: G_CANOPY       ( :,:,: )
      ! REAL(fp), POINTER :: A_CANOPY       ( :,:,: )
      ! REAL(fp), POINTER :: R_CANOPY       ( :,:,: )
      ! REAL(fp), POINTER :: A_NET_CANOPY   ( :,:,: )
      ! REAL(fp), POINTER :: FLXO3_CANOPY   ( :,:,: )
      ! REAL(fp), POINTER :: BETA_O3        ( :,:,: )
      ! REAL(fp), POINTER :: BETA_SM        ( :,:,: )   

      ! For ESMF, need to assign these from Input_Opt (copied from drydep_mod)
      LOGICAL       :: LPRT

      !=================================================================
      ! DO_ECOPHY begins here!
      !=================================================================
      ! Assume success 
      RC = GC_SUCCESS

      ! Initialize
      SpcInfo => NULL()
      ErrMsg  = ''
      ThisLoc = &
      ' -> at Do_ECOPHY (in module GeosCore/ecophysiology.F90)' 

      ! Point to columns of derived-type object fields
      ! G_CANOPY     => State_Chm%G_CAN     
      ! A_CANOPY     => State_Chm%A_CAN     
      ! R_CANOPY     => State_Chm%RESP     
      ! A_NET_CANOPY => State_Chm%A_NET 
      ! FLXO3_CANOPY => State_Chm%FLXO3 
      ! BETA_O3      => State_Chm%BETA_O3      
      ! BETA_SM      => State_Chm%BETA_SM   

      ! get inputs for the module
      CALL GET_ECOPHY_INPUTS( State_Met,    State_Chm, Input_Opt,&
                              TEMPK,        SPHU,                &
                              PAR_ABSORBED, PRESSURE,  CO2,      &
                              O2,           LAI,       O3,       &
                              LO3_DAMAGE,   SOIL_WETNESS         &
                              )  

      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in call to "GET_ECOPHY_INPUTS!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF 

      ! Get plant functional type index
      PFT = IPFT(LDT)

      ! simulate plant processes
      CALL DO_PHOTOSYNTHESIS( TEMPK,        SPHU,       RA,           &
                              PAR_ABSORBED, PRESSURE,   CO2,          &
                              O2,           LAI,        O3,           &
                              LO3_DAMAGE,   SOIL_WETNESS,             &
                              G_CAN_OUT,    A_CAN_OUT,  RESP_CAN_OUT, &
                              G_LEAF_OUT,   CO2_IN,     A_NET_OUT,    &
                              RESP_OUT,     FLUXO3_CAN, FLUXO3,       &
                              FACTOR_O3,    BETA,       PFT           &
                              )

      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in call to "GET_ECOPHY_INPUTS!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

! #if defined( NC_DIAG )
!             ! send to diagnostics outputs
!             State_Diag%G_CAN(I,J) = G_CAN_OUT
! #endif

      ! Write a subroutine to get PFT-weighted values

      ! Output RS to dry deposition module
      RS = 1.0 / G_CAN_OUT

      !### Debug
      IF ( LPRT .and. am_I_Root ) THEN
         CALL DEBUG_MSG( '### DO_DRYDEP: after dry dep' )
      ENDIF

      ! Nullify pointers
      NULLIFY( G_CANOPY     )
      NULLIFY( A_CANOPY     )
      NULLIFY( R_CANOPY     )
      NULLIFY( A_NET_CANOPY )
      NULLIFY( FLXO3_CANOPY )
      NULLIFY( BETA_O3      )
      NULLIFY( BETA_SM      )

      END SUBROUTINE DO_ECOPHY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_photosynthesis
!
! !DESCRIPTION: Subroutine DO\_PHOTOSYNTHESIS is the main driver of this 
!  module. 
!\\
!\\
! !INTERFACE:
!      
      SUBROUTINE DO_PHOTOSYNTHESIS( TEMPK,        SPHU,       RA,           &
                                    PAR_ABSORBED, PRESSURE,   CO2,          &
                                    O2,           LAI,        O3,           &
                                    LO3_DAMAGE,   SOIL_WETNESS,             &
                                    G_CAN_OUT,    A_CAN_OUT,  RESP_CAN_OUT, &
                                    G_LEAF_OUT,   CO2_IN,     A_NET_OUT,    &
                                    RESP_OUT,     FLUXO3_CAN, FLUXO3,       &
                                    FACTOR_O3,    BETA,       PFT           &
                                    )
! Main driver of the photosynthesis-stomatal conductance model
!
!INPUT PARAMETERS:
!
!
!INPUT PARAMETERS:
!
      !---------------------------------------------------------------------------------------
      ! TEMPK         : Leaf temperature in Kelvin                        [K]
      ! SPHU          : Specific humidity in canopy layer                 [kg H2O / kg air]
      ! RA            : Aerodynamic and boundary layer resistance         [s m^-1]
      ! PAR_ABSORBED  : Absorbed PAR                                      [W m^-2]
      ! PRESSURE      : Atmospheric Pressure in canopy layer              [Pa]
      ! CO2           : Ambient CO2 mole fraction                         [mol/mol air]
      ! O2            : Ambient O2 mole fraction                          [mol/mol air]
      ! O3            : Ozone mole fraction in canopy layer               [mol/mol air]
      ! LAI           : Leaf area index for the PFT                       [m^2 m^-2]
      ! LO3_DAMAGE    : Logical switch for ozone damage scheme            []
      ! SOIL_WETNESS  : Fraction of moisture in soil pores                []
      ! PFT           : Index for PFT                                     []
      !---------------------------------------------------------------------------------------
      REAL,     INTENT(IN)  :: TEMPK         
      REAL,     INTENT(IN)  :: SPHU          
      REAL,     INTENT(IN)  :: RA            
      REAL,     INTENT(IN)  :: PAR_ABSORBED  
      REAL,     INTENT(IN)  :: PRESSURE      
      REAL,     INTENT(IN)  :: CO2           
      REAL,     INTENT(IN)  :: O2            
      REAL,     INTENT(IN)  :: O3            
      REAL,     INTENT(IN)  :: LAI           
      LOGICAL,  INTENT(IN)  :: LO3_DAMAGE    
      REAL,     INTENT(IN)  :: SOIL_WETNESS  
      INTEGER,  INTENT(IN)  :: PFT           
!
!OUTPUT PARAMETERS:
!
      !---------------------------------------------------------------------------------------
      ! G_CAN_OUT     : Canopy conductance for H2O (output)               [m s^-1]
      ! G_LEAF_OUT    : Leaf level stomatal conductance for H2O (output)  [m s^-1]
      ! CO2_IN        : Leaf internal partial pressure of CO2             [Pa]
      ! A_CAN_OUT       : Canopy net photosynthetic rate (output)           [mol CO2 m^-2 s^-1]
      ! A_NET_OUT       : Leaf level net photosynthetic rate (output)       [mol CO2 m^-2 s^-1]
      ! RESP_CAN_OUT  : Canopy dark respiration (output)                  [mol CO2 m^-2 s^-1]
      ! RESP_OUT      : Leaf level dark respiration (output)              [mol CO2 m^-2 s^-1]
      ! FLUXO3_CAN    : Canopy ozone uptake                               [nmol m^-2 s^-1]
      ! FLUXO3        : Stomatal ozone uptake                             [nmol m^-2 s^-1]
      ! FACTOR_O3     : Ozone damage factor                               []
      ! BETA          : Soil moisture stress factor                       []
      !---------------------------------------------------------------------------------------
      REAL,     INTENT(OUT) :: G_CAN_OUT     
      REAL,     INTENT(OUT) :: G_LEAF_OUT    
      REAL,     INTENT(OUT) :: CO2_IN        
      REAL,     INTENT(OUT) :: A_CAN_OUT    
      REAL,     INTENT(OUT) :: A_NET_OUT     
      REAL,     INTENT(OUT) :: RESP_CAN_OUT  
      REAL,     INTENT(OUT) :: RESP_OUT      
      REAL,     INTENT(OUT) :: FLUXO3_CAN    
      REAL,     INTENT(OUT) :: FLUXO3        
      REAL,     INTENT(OUT) :: FACTOR_O3     
      REAL,     INTENT(OUT) :: BETA          
!
!LOCAL VARIABLES:
!
      !---------------------------------------------------------------------------------------
      ! TEMPC         : Leaf temperature in degree Celsius                [deg C]
      ! SPHU_SAT      : SPHU at saturation in canopy layer                [kg H2O / kg air]
      ! DEFICIT_Q     : Specific humidity deficit at leaf surface         [kg H2O / kg air]
      ! APAR          : Absorbed PAR                                      [mol m^-2 s^-1]
      ! CO2_AMBIENT   : Ambient CO2 partial pressure                      [Pa]
      ! O3_CONC       : Ozone concentration in canopy layer               [nmol m^-3]
      ! BIGLEAFSCALE  : Scaling with Big-leaf approach                    []
      ! G_CAN         : Canopy conductance for H2O                        [m s^-1]
      ! G_LEAF        : Leaf level stomatal conductance for H2O           [m s^-1]
      ! G_LEAF_PREV   : G_LEAF in previous iteration                      [m s^-1]
      ! CO2_IN_PREV   : CO2_IN in previous iteration                      [Pa]
      ! A_NET         : Leaf level net photosynthetic rate                [mol CO2 m^-2 s^-1]
      ! RESP          : Leaf level dark respiration                       [mol CO2 m^-2 s^-1]
      ! A_NET_PREV    : A_NET in previous iteration                       [mol CO2 m^-2 s^-1]
      ! V_CMAX        : Max Rubisco carboxylation rate                    [mol CO2 m^-2 s^-1]
      ! CO2_GAMMA     : CO2 Compensation point                            [Pa]
      ! RATE_LIGHT    : Light-limited rate                                [mol CO2 m^-2 s^-1]
      ! RATE_PRODUCT  : Product-limited rate                              [mol CO2 m^-2 s^-1]
      ! RATE_RUBISCO  : Rubisco-limited rate                              [mol CO2 m^-2 s^-1]
      ! A_GROSS       : Leaf level gross photosynthesis                   [mol CO2 m^-2 s^-1]
      ! TAU           : Rubisco specificity for CO2 to O2                 []
      ! DENOM         : Denominator for calculating V_CMAX                []
      ! ITER          : No. of iterations in calculating conductance      []
      ! ERR1          : Relative change for G_LEAF between iterations     []
      ! ERR2          : Relative change for CO2_IN between iterations     []
      ! ERR3          : Relative change for A_NET between iterations      []
      ! DELTA         : Maximum of the 3 relative changes                 []
      !---------------------------------------------------------------------------------------
      REAL                  :: TEMPC         
      REAL                  :: SPHU_SAT      
      REAL                  :: DEFICIT_Q     
      REAL                  :: APAR          
      REAL                  :: CO2_AMBIENT   
      REAL                  :: O3_CONC       
      REAL                  :: BIGLEAFSCALE  
      REAL                  :: G_CAN         
      REAL                  :: G_LEAF        
      REAL                  :: G_LEAF_PREV   
      REAL                  :: CO2_IN_PREV   
      REAL                  :: A_NET         
      REAL                  :: RESP          
      REAL                  :: A_NET_PREV    
      REAL                  :: V_CMAX        
      REAL                  :: CO2_GAMMA     
      REAL                  :: RATE_LIGHT    
      REAL                  :: RATE_PRODUCT  
      REAL                  :: RATE_RUBISCO  
      REAL                  :: A_GROSS       
      REAL                  :: TAU           
      REAL                  :: DENOM         
      INTEGER               :: ITER          
      REAL                  :: ERR1          
      REAL                  :: ERR2          
      REAL                  :: ERR3          
      REAL                  :: DELTA         




      TEMPC             = TEMPK - 273.15
      ! Calculate V_CMAX and respiration which depends on V_CMAX only
      DENOM             = ( 1.0 + EXP( 0.3*( TEMPC - T_UPP(PFT) ) ) ) &
                        * ( 1.0 + EXP( 0.3*( T_LOW(PFT) - TEMPC ) ) )
      V_CMAX            = V_CMAX25(PFT) * FACTOR_Q10( 2.0, TEMPC ) / DENOM
      RESP              = F_DARKRESP(PFT) * V_CMAX
      ! Calculate CO2 compensation point
      TAU               = 2600.0 * FACTOR_Q10( 0.57, TEMPC )                    ! CO2/O2 Specificity Ratio
      CO2_GAMMA         = IS_C3_PLANT(PFT) / ( 2.0 * TAU ) &
                        * PRESSURE * O2 
      ! Calculate canopy scaling factor 
      BIGLEAFSCALE      = ( 1 - EXP( -K_EXTINCT(PFT) * LAI ) ) / K_EXTINCT(PFT)
      ! Convert unit of absorbed PAR to mol photon m^-2 s^-1
      APAR              = 4.6e-6 * PAR_ABSORBED
      ! Calculate CO2 partial pressure in ambient air
      CO2_AMBIENT       = PRESSURE * CO2
      ! Calculate O3 molar concentration in canopy layer
      O3_CONC           = O3 * PRESSURE / RSTARG / TEMPK * 1e9
    
      ! To modify net photosynthesis rate by soil moisture stress later
      ! Not needed to be inside the loop
      CALL MOIST_STRESS( SOIL_WETNESS, THETA_SATU, THETA_CRIT, &
                         THETA_WILT,   BETA                    )
      
      ! Iterate to find a self-consistent set of photosynthesis,
      ! stomatal conductance and leaf internal CO2 concentration 
      ! Initial guess: G_LEAF = 0 and other initializations
      ITER              = 1 
      G_LEAF            = 0.0   
      G_CAN             = 0.0
      CO2_IN_PREV       = 0.0
      A_NET_PREV        = 0.0
      G_LEAF_PREV       = 0.0
      ERR1              = 1.0
      ERR2              = 1.0
      ERR3              = 1.0
      DELTA             = 1.0
      DO WHILE ( DELTA >= THRESHOLD .AND. ITER <= 100 )
        ! Step 1: Closure condition by Jacobs (1994)
        SPHU_SAT          = 0.622 * E_SAT( TEMPC ) / PRESSURE
        G_CAN             = G_LEAF * BIGLEAFSCALE
        DEFICIT_Q         = ( SPHU_SAT - SPHU ) / ( 1 + RA * G_CAN )     
        CO2_IN            = CO2_GAMMA + f0(PFT)*( 1 - DEFICIT_Q / D_STAR(PFT) ) &
                          * ( CO2_AMBIENT - CO2_GAMMA ) 
        IF ( BETA == 0.0 .OR. DEFICIT_Q >= D_STAR(PFT) .OR. PAR_ABSORBED == 0.0 ) THEN 
          ! Close stomata if the above conditions are satisfied
          A_NET           = - RESP * BETA 
          G_LEAF          = G_LEAF_MIN(PFT)
          PRINT *, "Stomata is closed."
        ELSE 
          ! Step 2: Photosynthesis model
          CALL PHOTOSYNTHESIS_LIMITS( CO2_IN,       CO2_GAMMA,    &
                                      O2, APAR,     PRESSURE,     &
                                      TEMPC,        V_CMAX,       & 
                                      PFT,         RATE_LIGHT,    &
                                      RATE_PRODUCT, RATE_RUBISCO  )                                
          CALL SOLVE_COLIMIT( RATE_LIGHT,   RATE_PRODUCT, &
                              RATE_RUBISCO, A_GROSS      )
          PRINT *, "Photosynthesis calculated."
          A_NET    = ( A_GROSS - RESP ) * BETA
          ! Step 3: Diffusive CO2 flux thru open stomata
          CALL LEAF_CONDUCTANCE( A_NET, CO2_AMBIENT, CO2_IN,  &
                                 TEMPK, G_LEAF                )
          ! Close stomata if net photosynthesis <= 0 or 
          ! stomatal conductance is too small
          IF ( A_NET <= 0 .OR. G_LEAF <= G_LEAF_MIN(PFT) ) THEN
            A_NET     = - RESP * BETA 
            G_LEAF    = G_LEAF_MIN(PFT)
          END IF
          ! Apply ozone damage scheme by Sitch et al. (2007)
          IF ( LO3_DAMAGE ) THEN
            CALL OZONE_DAMAGE ( O3_CONC,  RA,          &    
                                G_LEAF,   PFT,         &
                                FLUXO3,   FACTOR_O3    )
            PRINT *, "Ozone damage calculated."
            A_NET_OUT     = FACTOR_O3 * A_NET
            RESP_OUT      = FACTOR_O3 * RESP    
            G_LEAF_OUT    = FACTOR_O3 * G_LEAF
              ! Close stomata if net photosynthesis <= 0 or 
              ! stomatal conductance is too small
              IF ( A_NET_OUT <= 0 .OR. G_LEAF_OUT <= G_LEAF_MIN(PFT) ) THEN
                A_NET_OUT   = - RESP * BETA 
                G_LEAF_OUT  = G_LEAF_MIN(PFT)
              END IF
          ELSE 
            A_NET_OUT     = A_NET
            RESP_OUT      = RESP
            G_LEAF_OUT    = G_LEAF
          END IF  ! O3 damage
        END IF    ! Open or closed stomata
        IF ( ITER >= 2 ) THEN 
        ! calculate error from step 2 onwards
          ERR1  = REL_ERR( G_LEAF, G_LEAF_PREV )
          ERR2  = REL_ERR( CO2_IN, CO2_IN_PREV )
          ERR3  = REL_ERR( A_NET,  A_NET_PREV  )
          DELTA = ABS( MAX( ERR1, ERR2, ERR3 ) )
        END IF
        CO2_IN_PREV  = CO2_IN
        A_NET_PREV   = A_NET
        G_LEAF_PREV  = G_LEAF        
        PRINT *, "iteration = ", ITER
        PRINT *, "Delta = ", DELTA
        ITER = ITER + 1
      END DO  ! Do while loop
      
      ! Canopy scaling
      A_CAN_OUT     = BIGLEAFSCALE * A_NET_OUT
      G_CAN_OUT     = BIGLEAFSCALE * G_LEAF_OUT
      RESP_CAN_OUT  = BIGLEAFSCALE * RESP_OUT
      FLUXO3_CAN    = BIGLEAFSCALE * FLUXO3 
      
      ! Deal with diagnoses

      ! Print for debugs
#if defined (DEBUG) 
        PRINT *, "TEMPC = ", TEMPC
        PRINT *, "SPHU_SAT = ", SPHU_SAT     
        PRINT *, "DEFICIT_Q = ", DEFICIT_Q    
        PRINT *, "G_LEAF_PREV = ", G_LEAF_PREV   
        PRINT *, "CO2_IN_PREV = ",  CO2_IN_PREV   
        PRINT *, "A_NET_PREV = ", A_NET_PREV    
        PRINT *, "V_CMAX = ", V_CMAX        
        PRINT *, "CO2_GAMMA = ", CO2_GAMMA     
        PRINT *, "RATE_LIGHT = ", RATE_LIGHT    
        PRINT *, "RATE_PRODUCT = ", RATE_PRODUCT  
        PRINT *, "RATE_RUBISCO = ", RATE_RUBISCO  
        PRINT *, "A_GROSS = ", A_GROSS       
        PRINT *, "TAU = ", TAU           
        PRINT *, "DENOM = ", DENOM         
        PRINT *, "ITER = ", ITER          
        PRINT *, "ERR1 = ", ERR1          
        PRINT *, "ERR2 = ", ERR2          
        PRINT *, "ERR3 = ", ERR3          
        PRINT *, "DELTA = ", DELTA 
#endif
      END SUBROUTINE DO_PHOTOSYNTHESIS
      
      SUBROUTINE PHOTOSYNTHESIS_LIMITS( CO2_IN,       CO2_GAMMA,    &
                                        O2, APAR,     PRESSURE,     &
                                        TEMPC,        V_CMAX,       & 
                                        PFT,          RATE_LIGHT,   &
                                        RATE_PRODUCT, RATE_RUBISCO  )
! Determine potential leaf-level photosynthesis, according to C3 and C4
! photosynthesis model from Collatz et al. (1991) and Collatz et al. (1992)
!
! !INPUT PARAMETERS:
!
      !---------------------------------------------------------------------------------------
      ! CO2_IN        : Leaf internal partial pressure of CO2 [Pa]
      ! CO2_GAMMA     : CO2 Compensation point [Pa]
      ! O2            : Ambient O2 mole fraction                          [mol/mol air]
      ! APAR          : Absorbed photosynthetically active radiation (PAR) [mol photon m^-2 s^-1]
      ! PRESSURE      : Surface air pressure [Pa]
      ! TEMPC         : Temperature [deg C]
      ! V_CMAX        : Maximum rate of carboxylation of Rubisco [mol CO2 m^-2 s^-1]
      ! PFT           : Index for PFT                                     []
      !---------------------------------------------------------------------------------------
      REAL,     INTENT(IN)  :: CO2_IN
      REAL,     INTENT(IN)  :: CO2_GAMMA
      REAL,     INTENT(IN)  :: O2
      REAL,     INTENT(IN)  :: APAR       
      REAL,     INTENT(IN)  :: PRESSURE  
      REAL,     INTENT(IN)  :: TEMPC    
      REAL,     INTENT(IN)  :: V_CMAX
      INTEGER,  INTENT(IN)  :: PFT
!
! !OUTPUT PARAMETERS:
!
      !---------------------------------------------------------------------------------------
      ! RATE_LIGHT         : Light-limited rate [mol CO2 m^-2 s^-1]
      ! RATE_PRODUCT       : Product-limited rate [mol CO2 m^-2 s^-1]
      ! RATE_RUBISCO       : Rubisco-limited rate [mol CO2 m^-2 s^-1]
      !---------------------------------------------------------------------------------------
      REAL,     INTENT(OUT) :: RATE_LIGHT
      REAL,     INTENT(OUT) :: RATE_PRODUCT
      REAL,     INTENT(OUT) :: RATE_RUBISCO
      ! Success or failure flag
!      INTEGER,  INTENT(OUT) :: RC
!
! !LOCAL VARIABLES: 
! 
      ! Michaelis-Menten parameters for CO2 and O2 respectively
      REAL                  :: K_C 
      REAL                  :: K_O 
      
      !=================================================================
      ! PHOTOSYNTHESIS_LIMITS begins here!
      !=================================================================
      ! Assume success
!      RC                    = GC_SUCCESS
      K_C                   = 30.0 * FACTOR_Q10( 2.1, TEMPC )
      K_O                   = 30000.0 * FACTOR_Q10( 1.2, TEMPC )
      ! For C4 plants
      IF ( IS_C3_PLANT(PFT) == 0 ) THEN
        RATE_RUBISCO        = V_CMAX
        RATE_LIGHT          = ALPHA(PFT) * APAR
        RATE_PRODUCT        = 20000.0 * V_CMAX * CO2_IN / PRESSURE
      ELSE    ! For C3 plants
        RATE_RUBISCO        = V_CMAX * ( CO2_IN - CO2_GAMMA )  &
                            / ( CO2_IN + K_C * ( 1 + PRESSURE * O2 / K_O ) )
        RATE_LIGHT          = ALPHA(PFT) * APAR               &
                            * ( CO2_IN -     CO2_GAMMA )       &
                            / ( CO2_IN + 2 * CO2_GAMMA )
        RATE_RUBISCO        = MAX( RATE_RUBISCO, 0.0 )
        RATE_LIGHT          = MAX( RATE_LIGHT, 0.0 )         
        RATE_PRODUCT        = 0.5 * V_CMAX
      END IF
      END SUBROUTINE PHOTOSYNTHESIS_LIMITS
      
      SUBROUTINE SOLVE_COLIMIT( RATE_LIGHT,   RATE_PRODUCT, &
                                RATE_RUBISCO, A_GROSS       )
      !---------------------------------------------------------------------------------------
      ! RATE_LIGHT         : Light-limited rate [mol CO2 m^-2 s^-1]
      ! RATE_PRODUCT       : Product-limited rate [mol CO2 m^-2 s^-1]
      ! RATE_RUBISCO       : Rubisco-limited rate [mol CO2 m^-2 s^-1]
      ! A_GROSS            : Gross rate of photosynthesis [mol CO2 m^-2 s^-1]
      !---------------------------------------------------------------------------------------
      REAL,     INTENT(IN)  :: RATE_LIGHT
      REAL,     INTENT(IN)  :: RATE_PRODUCT
      REAL,     INTENT(IN)  :: RATE_RUBISCO
      REAL,     INTENT(OUT) :: A_GROSS 
      
      ! Parameters
      REAL,     PARAMETER   :: BETA1 = 0.83
      REAL,     PARAMETER   :: BETA2 = 0.93
      ! Local parameter
      REAL                  :: TEMP
      REAL                  :: B
      REAL                  :: C
      ! 1st quadratic
      B = -( RATE_RUBISCO + RATE_LIGHT ) / BETA1
      C = RATE_RUBISCO * RATE_LIGHT / BETA1
      ! Note that C > 0, SQRT( B^2 - 4*C ) < ABS(B)
      ! Take smaller root
      TEMP = 0.5 * ( - B - SQRT( B * B - 4 * C ) ) 
      
      ! 2nd quadratic 
      B = - ( TEMP + RATE_PRODUCT ) / BETA2
      C = TEMP * RATE_PRODUCT / BETA2
      ! Note that C > 0, SQRT( B^2 - 4*C ) < ABS(B)
      ! Take smaller root
      A_GROSS  = 0.5 * ( - B - SQRT( B * B - 4 * C ) ) 
      END SUBROUTINE SOLVE_COLIMIT
      
      SUBROUTINE MOIST_STRESS( SOIL_WETNESS, THETA_SATU, THETA_CRIT, &
                               THETA_WILT,   BETA                    )
! Calculate the moisture stress factor BETA

      !---------------------------------------------------------------------------------------
      ! SOIL_WETNESS    : Volumetric mean moisture concentration in root zone divided by porosity 
      ! THETA_SATU      : Volumetric soil moisture at saturation (= porosity)
      ! THETA_CRIT      : Volumetric soil moisture at critical point (above which
      !                   plants are no longer stressed by soil moisture)
      ! THETA_WILT      : Volumetric soil moisture at wilting point (below which
      !                   photosynthesis is stopped by limited soil moisture)
      ! BETA            : Moisture stress factor
      !---------------------------------------------------------------------------------------
      REAL,     INTENT(IN)  :: SOIL_WETNESS 
      REAL,     INTENT(IN)  :: THETA_SATU 
      REAL,     INTENT(IN)  :: THETA_CRIT 
      REAL,     INTENT(IN)  :: THETA_WILT   
      REAL,     INTENT(OUT) :: BETA
      
      BETA =  ( SOIL_WETNESS - THETA_WILT ) / ( THETA_CRIT - THETA_WILT ) 
      BETA = MIN( MAX( 0.0, BETA ), 1.0 )
      END SUBROUTINE MOIST_STRESS
      
      SUBROUTINE LEAF_CONDUCTANCE( A_NET, CO2_AMBIENT, CO2_IN,  &
                                   TEMPK, G_LEAF                )
! Calculate the leaf conductance based on net photosynthesis and CO2 partial pressure deficit
      !---------------------------------------------------------------------------------------
      ! A_NET             : Net photosynthesis                      [mol CO2 m^-2 s^-1]
      ! CO2_AMBIENT       : Ambient CO2 partial pressure            [Pa] 
      ! CO2_IN            : Leaf internal CO2 partial pressure      [Pa]
      ! TEMPK             : Leaf temperature                        [K]
      ! G_LEAF            : Leaf conductance for H2O                [m s^-1]
      ! RSTARG            : Universal Gas Constant                  [J K^-1 mol^-1]
      ! CO2_O2_RATIO      : Ratio of leaf resistance for CO2 to H2O 
      !---------------------------------------------------------------------------------------
      REAL,     INTENT(IN)      :: A_NET
      REAL,     INTENT(IN)      :: CO2_AMBIENT
      REAL,     INTENT(IN)      :: CO2_IN
      REAL,     INTENT(IN)      :: TEMPK
      REAL,     INTENT(OUT)     :: G_LEAF
      G_LEAF = CO2_O2_RATIO * RSTARG * TEMPK  &
             * A_NET / ( CO2_AMBIENT - CO2_IN )
      END SUBROUTINE LEAF_CONDUCTANCE
      
      SUBROUTINE OZONE_DAMAGE ( O3_CONC, RA,          &                          ! Should be Ra + Rb? Need to check
                                G_LEAF,  PFT,         &
                                FLUXO3,  FACTOR_O3    )
! Calculate ozone damage factor based on Sitch et al. (2008)
      !---------------------------------------------------------------------------------------
      ! O3_CONC         : Ozone concentration in canopy layer                     [nmol m^-3]
      ! RA              : Aerodynamic and boundary resistance                     [s m^-1]
      ! G_LEAF          : Leaf conductance for H2O in the absence of O3 effects   [m s^-1]
      ! PFT             : Index for PFT                                           []
      ! FLUXO3          : Leaf uptake of O3                                       [nmol m^-2 s^-1]
      ! FACTOR_O3       : Ozone damage factor                                     []
      !---------------------------------------------------------------------------------------        
      REAL,     INTENT(IN)    :: O3_CONC         
      REAL,     INTENT(IN)    :: RA              
      REAL,     INTENT(IN)    :: G_LEAF
      INTEGER,  INTENT(IN)    :: PFT
      REAL,     INTENT(OUT)   :: FLUXO3
      REAL,     INTENT(OUT)   :: FACTOR_O3
      ! Local variables
      REAL :: B       ! Second coefficient in quadratic equation
      REAL :: C       ! Constant in quadratic equation
      REAL :: TEMP1
      REAL :: TEMP2
      REAL :: F
      TEMP1       = 1 + PARAM_A(PFT) * FLUXO3_CRIT(PFT)
      TEMP2       = 1.67 / G_LEAF        
      ! Calculate coefficients for quadratic equation F^2 + B*F + C = 0
      IF ( ABS(RA) < EPSILON(1.0) ) THEN
!        RA        = 0.0  
        FACTOR_O3 = TEMP1 / ( 1 + PARAM_A(PFT) * O3_CONC / TEMP2 )
      ELSE
        B         = TEMP2 / RA - TEMP1 + PARAM_A(PFT) * O3_CONC / RA 
        C         = - TEMP1 * TEMP2 / RA
        ! Note that C < 0, SQRT( B^2 - 4*C ) > ABS(B)
        ! Take positive root
        F         = 0.5 * ( SQRT( B * B - 4 * C ) - B ) 
        FACTOR_O3 = MIN( MAX( F, 0.0 ), 1.0 )
      END IF
      FLUXO3      = O3_CONC / ( RA + TEMP2 / FACTOR_O3 )      ! MAYBE NOT NEEDED?
      END SUBROUTINE OZONE_DAMAGE
      
      FUNCTION FACTOR_Q10( Q10, TEMPC ) RESULT( FACTOR )
! Calculate the standard Q10 temperature dependence
      REAL, INTENT(IN)    :: Q10
      REAL, INTENT(IN)    :: TEMPC              ! Temperature [deg C]
      REAL                :: FACTOR
      FACTOR = Q10**( 0.1 * ( TEMPC - 25 ) )
      END FUNCTION FACTOR_Q10
      
      FUNCTION REL_ERR( ITEM, ITEM_PREV ) RESULT( ERROR )
! Calculate the relative error of a quantity between two iterations
      REAL, INTENT(IN)    :: ITEM
      REAL, INTENT(IN)    :: ITEM_PREV
      REAL                :: ERROR
      ERROR = ( ITEM - ITEM_PREV ) / ITEM_PREV 
      END FUNCTION REL_ERR      
      
      FUNCTION E_SAT( TEMPC ) RESULT ( Esat )
! Calculate the saturation vapour pressure using the empirical formula by Lowe and Ficke (1974)
      REAL, INTENT(IN)  :: TEMPC
      REAL              :: Esat
      ! Local parameters
      REAL, PARAMETER   :: a0 = 6.107799961
      REAL, PARAMETER   :: a1 = 4.436518521e-1
      REAL, PARAMETER   :: a2 = 1.428945805e-2
      REAL, PARAMETER   :: a3 = 2.650648471e-4
      REAL, PARAMETER   :: a4 = 3.031240396e-6
      REAL, PARAMETER   :: a5 = 2.034080948e-8
      REAL, PARAMETER   :: a6 = 6.136820929e-11
      Esat = 100.0 * ( a0 + TEMPC * ( a1 + TEMPC * ( a2 + TEMPC   &
                   * ( a3 + TEMPC * ( a4 + TEMPC * ( a5 + TEMPC * a6 ) ) ) ) ) )
      END FUNCTION E_SAT
      
!      SUBROUTINE PFT_PARAMS( NUMPFT,      K_EXTINCT,   ALPHA,       &
!                             OMEGA,       F_DARKRESP,  K_EXTINCT,   &
!                             R_GROWTH,    LEAF_N_TOP,  T_LOW,       &
!                             T_UPP,       FLUXO3_CRIT, FACTOR_O3,   &
!                             G_LEAF_MIN,  IS_C4_PLANT, f0,          &
!                             )
!! Define all PFT-specific parameters
!! Input parameters
!      INTEGER,  INTENT(IN)    :: NUMPFT
!! Output parameters 
!      REAL,     INTENT(OUT)   :: K_EXTINCT   (NUMPFT) ! PAR extinction coefficient [m^-1]
!      REAL,     INTENT(OUT)   :: ALPHA       (NUMPFT) ! Quantum efficiency [mol CO2 m^-2 s^-1]
!      REAL,     INTENT(OUT)   :: OMEGA       (NUMPFT) ! Leaf scattering coefficient for PAR
!      REAL,     INTENT(OUT)   :: F_DARKRESP  (NUMPFT) ! Dark respiration coefficient
!      REAL,     INTENT(OUT)   :: R_GROWTH    (NUMPFT) ! Growth respiration coefficient
!      REAL,     INTENT(OUT)   :: LEAF_N_TOP  (NUMPFT) ! Leaf N concentration at canopy top
!      REAL,     INTENT(OUT)   :: T_LOW       (NUMPFT) ! Lower temperature parameter [deg C]
!      REAL,     INTENT(OUT)   :: T_UPP       (NUMPFT) ! Upper temperature parameter [deg C]
!      REAL,     INTENT(OUT)   :: FLUXO3_CRIT (NUMPFT) ! Threshold ozone flux [nmol m^-2 s^-1]
!      REAL,     INTENT(OUT)   :: FACTOR_O3   (NUMPFT) ! Ozone factor [mmol^-1 m^-2]
!      REAL,     INTENT(OUT)   :: G_LEAF_MIN  (NUMPFT) ! Minimum leaf conductance for H2O [m s^-1]
!      INTEGER,  INTENT(OUT)   :: IS_C4_PLANT (NUMPFT) ! IS_C4_PLANT = 1 if the PFT is C4 plant, 
!                                                      !             = 0 if it is not
!      REAL,     INTENT(OUT)   :: f0          (NUMPFT) ! Calibration parameter, from Cox et al. (1998)
!      REAL,     INTENT(OUT)   :: D_STAR      (NUMPFT) ! Calibration parameter, from Cox et al. (1998)
!      
!      END SUBROUTINE PFT_PARAMS

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ecophy_inputs
!
! !DESCRIPTION: Subroutine GET\_ECOPHY_INPUTS get time-dependent inputs 
!  from Met and Chem State objects and Input Options.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_ECOPHY_INPUTS( State_Met,    State_Chm,           &
                                    TEMPK,        SPHU,                &
                                    PAR_ABSORBED, PRESSURE,  CO2,      &
                                    O2,           LAI,       O3,       &
                                    LO3_DAMAGE,   SOIL_WETNESS         &
                                    )
!
! !USES:
!
      USE State_Chm_Mod,      ONLY : ChmState, Ind_
      USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
      !---------------------------------------------------------------------------------------
      ! State_Met     : Meteorology State Object
      ! State_Chm     : Chemistry State Object
      !---------------------------------------------------------------------------------------
      Type(MetState), INTENT(IN)  :: State_Met
      Type(ChmState), INTENT(IN)  :: State_Chm
!
! !OUTPUT PARAMETERS:
!
      !---------------------------------------------------------------------------------------
      ! TEMPK         : Temperature in Kelvin                             [K]
      ! SPHU          : Specific humidity in canopy layer                 [kg H2O / kg air]
      ! RA            : Aerodynamic and boundary layer resistance         [s m^-1]
      ! PAR_ABSORBED  : Absorbed PAR                                      [W m^-2]
      ! PRESSURE      : Atmospheric Pressure in canopy layer              [Pa]
      ! CO2           : Ambient CO2 mole fraction                         [kg / kg dry]
      ! O2            : Ambient O2 mole fraction                          [kg / kg dry]
      ! O3            : Ozone mole fraction in canopy layer               [kg / kg dry]
      ! LAI           : Leaf area index for the PFT                       [m^2 m^-2]
      ! LO3_DAMAGE    : Logical switch for ozone damage scheme            []
      ! SOIL_WETNESS  : Fraction of moisture in soil pores                []
      !---------------------------------------------------------------------------------------
      REAL,     INTENT(OUT) :: TEMPK         
      REAL,     INTENT(OUT) :: SPHU          
      REAL,     INTENT(OUT) :: RA            
      REAL,     INTENT(OUT) :: PAR_ABSORBED  
      REAL,     INTENT(OUT) :: PRESSURE      
      REAL,     INTENT(OUT) :: CO2           
      REAL,     INTENT(OUT) :: O2            
      REAL,     INTENT(OUT) :: O3            
      REAL,     INTENT(OUT) :: LAI           
      LOGICAL,  INTENT(OUT) :: LO3_DAMAGE    
      REAL,     INTENT(OUT) :: SOIL_WETNESS  
!
! !LOCAL VARIABLES:
!
      REAL    :: PARDR, PARDF, PAR
      INTEGER :: I, J
      INTEGER :: id_CO2, id_O2, id_O3

      ! Find tracer indices with function the Ind_() function
      id_CO2    = IND_( 'CO2' )
      id_O2     = IND_( 'O2'  )
      id_O3     = IND_( 'O3'  )

      ! ! Loop over surface grid boxes
      ! DO J = 1, JJPAR
      ! DO I = 1, IIPAR
      ! Surface temperature [K]
      TEMPK         = State_Met%TS( I,J )
      ! Specific humidity [kg/kg]  
      SPHU          = State_Met%SPHU( I,J,1 ) / 1000
      ! Photosynthetically active radiation absorbed [W m^-2] 
      PARDR         = State_Met%PARDR( I,J ) 
      PARDF         = State_Met%PARDF( I,J )
      PAR           = PARDR + PARDF
      ! Pressure [Pa]
      PRESSURE      = State_Met%SLP( I,J )
      ! CO2 mole fraction [mol/mol]
      CO2           = State_Chm%Species( I,J,1,id_CO2 ) * AIRMW &
                    \ State_Chm%SpcData( id_CO2 )%Info%MW_g
      ! O2 mole fraction [mol/mol]
      O2            = State_Chm%Species( I,J,1,id_O2  ) * AIRMW &
                    \ State_Chm%SpcData( id_O2  )%Info%MW_g
      ! O3 mole fraction [mol/mol]
      O3            = State_Chm%Species( I,J,1,id_O3  ) * AIRMW &
                    \ State_Chm%SpcData( id_O3  )%Info%MW_g
      ! LAI [m^2 m^-2]
      LAI           = State_Met%XLAI( I,J,: )
      ! Root zone soil wetness
      SOIL_WETNESS  = State_Met%GWETROOT(I,J)
       ! END DO
       ! END DO
      END SUBROUTINE GET_ECOPHY_INPUTS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_ecophy
!
! !DESCRIPTION: Subroutine INIT\_ECOPHY initializes certain variables for the
!  GEOS-CHEM ecophysiology subroutines. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_ECOPHY( am_I_Root, Input_Opt, 
                              State_Chm, State_Diag, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod,  ONLY : OptInput
      USE Species_Mod,    ONLY : Species
      USE State_Chm_Mod,  ONLY : ChmState
      USE State_Chm_Mod,  ONLY : Ind_
      USE State_Diag_Mod, ONLY : DgnState
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL                :: LECOPHY
      LOGICAL                :: LO3_DAMAGE
      INTEGER                :: NUMPFT
      INTEGER                :: N

      ! Strings
      CHARACTER(LEN=255)     :: Msg, ErrMsg, ThisLoc

      !=================================================================
      ! INIT_ECOPHY begins here!
      !=================================================================
      
      ! Initialize
      RC        = GC_SUCCESS
      LECOPHY   = Input_Opt%LECOPHY
      LO3_DAMAGE= Input_Opt%LO3_DAMAGE
      NUMPFT    = 5
      ErrMsg    = ''
      ThisLoc   = ' -> at Init_Ecophy (in module GeosCore/ecophysiology.F90)'
      
      !===================================================================
      ! Arrays that hold inputs for the main driver do_ecophy
      ! Only allocate these if ecophysiology is activated
      !===================================================================
      ALLOCATE( THETA_SATU( IIPAR, JJPAR ), STAT=RC )
      CALL GC_CheckVar( 'ecophy_mod:THETA_SATU', 0, RC)
      IF ( RC /= GC_SUCCESS ) RETURN 
      THETA_SATU( :,: ) = 0e+0_f8

      ALLOCATE( THETA_CRIT( IIPAR, JJPAR ), STAT=RC )
      CALL GC_CheckVar( 'ecophy_mod:THETA_CRIT', 0, RC)
      IF ( RC /= GC_SUCCESS ) RETURN 
      THETA_CRIT( :,: ) = 0e+0_f8

      ALLOCATE( THETA_WILT( IIPAR, JJPAR ), STAT=RC )
      CALL GC_CheckVar( 'ecophy_mod:THETA_WILT', 0, RC)
      IF ( RC /= GC_SUCCESS ) RETURN 
      THETA_WILT( :,: ) = 0e+0_f8

      ALLOCATE( IPFT( NUMPFT ), STAT=RC )
      CALL GC_CheckVar( 'ecophy_mod:IPFT', 0, RC)
      IF ( RC /= GC_SUCCESS ) RETURN 
      IPFT( : ) = 0

      ! Get soil map
      ! Get PFT mapping?

!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CLEANUP\_ECOPHY
!
! !DESCRIPTION: Subroutine CLEANUP_ECOPHY deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_ECOPHY()
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC

      !=================================================================
      ! CLEANUP_ECOPHY begins here!
      !=================================================================
      IF ( ALLOCATED ( THETA_SATU ) ) DEALLOCATE ( THETA_SATU )
      IF ( ALLOCATED ( THETA_CRIT ) ) DEALLOCATE ( THETA_CRIT )
      IF ( ALLOCATED ( THETA_WILT ) ) DEALLOCATE ( THETA_WILT )
      IF ( ALLOCATED ( IPFT       ) ) DEALLOCATE ( IPFT       )

      ! Return to calling program
      END SUBROUTINE CLEANUP_ECOPHY
!EOC
      END MODULE ECOPHY_MOD


      

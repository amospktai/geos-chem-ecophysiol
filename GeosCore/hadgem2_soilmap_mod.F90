!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hadgem2_soilmap_mod.F90
!
! !DESCRIPTION: Module HADGEM2\_SOILMAP\_MOD reads the HadGEM2 soil map and
!  computes the THETA_WILT and THETA_CRIT State\_Met arrays. 
!\\
!\\
! !INTERFACE: 
!
MODULE HadGEM2_SoilMap_Mod
!
! !USES:
!
!  USE CMN_SIZE_MOD                      ! Size parameters
  USE ERROR_MOD                         ! Error checking routines
  USE GC_GRID_MOD                       ! Horizontal grid definition
  USE MAPPING_MOD                       ! Mapping weights & areas
!  USE PhysConstants                     ! Physical constants
  USE PRECISION_MOD                     ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
!#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
  PUBLIC  :: Compute_Soilmap_GCHP
!#else
  PUBLIC  :: Init_Soilmap
  PUBLIC  :: Compute_Soilmap
  PUBLIC  :: Cleanup_SoilMap
!#endif
!                                                                                              
!  The variables are defined as follows:
!      State_Met%THETA_WILT(I,J)    : Volumetric soil moisture content at wilting point 
!      State_Met%THETA_CRIT(I,J)    : Volumetric soil moisture content at critical point
!                                                  
!  NOTES: 
!  (1) THETA_WILT and THETA_CRIT are used by ecophysiology modules
!
! !REVISION HISTORY: 
!  09 Feb 2019 - J. Lam - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Scalars
  INTEGER              :: I_SOIL       ! # of lons (0.5 x 0.5)
  INTEGER              :: J_SOIL       ! # of lats (0.5 x 0.5)
  REAL(fp)             :: D_LON         ! Delta longitude, Olson grid [degrees]
  REAL(fp)             :: D_LAT         ! Delta latitude,  Olson grid [degrees]

  ! Arrays
  REAL*4,  ALLOCATABLE :: lon  (:    )  ! Lon centers, Olson grid [degrees]
  REAL*4,  ALLOCATABLE :: lat  (  :  )  ! Lat centers, Olson grid [degrees]
  INTEGER, ALLOCATABLE :: OLSON(:,:,:)  ! Olson land types
  REAL*4,  ALLOCATABLE :: A_CM2(:,:,:)  ! Surface areas [cm2]

CONTAINS
!EOC
!#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: compute_soilmap_gchp
!
! !DESCRIPTION: Subroutine COMPUTE\_SOILMAP\_GCHP computes the 
!  GEOS-Chem State\_Met variables THETA_WILT and THETA_CRIT.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Soilmap_GCHP( am_I_Root, State_Met, RC )
!
! !USES:
!
    USE State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN)    :: am_I_Root    ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState),  INTENT(INOUT) :: State_Met    ! Meteorology State object
    INTEGER,         INTENT(INOUT) :: RC
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  27 Sep 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER        :: I, J, T
    INTEGER        :: typeCounter, maxFracInd(1), sumIUSE

    !======================================================================
    ! Initialize
    !======================================================================

    
    ! Loop over all grid cells to set State_Met variables
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Initialize fraction land for this grid cell
       State_Met%FRCLND(I,J) = 1.e+0_fp ! Initialized as all land

       ! Initialize local variables
       typeCounter = 0       ! Tally of number of types found in the cell
       maxFracInd  = 0       ! type index with greatest coverage
       sumIUSE     = 0       ! total coverage across all types [mil]

       ! Loop over all landmap types to set IREG, ILAND, and IUSE
       DO T = 1, NSURFTYPE

          ! If this type has non-zero coverage in this grid box, update vars
          IF ( State_Met%LandTypeFrac(I,J,T) > 0.e+0_fp ) THEN

             ! Increment number of types in this cell
             typeCounter = typeCounter + 1

             ! Set IREG to number of types
             State_Met%IREG(I,J) = typeCounter

             ! Store type index in ILAND array for this grid cell.
             ! Use 0-based index for compatibility with legacy drydep code.
             State_Met%ILAND(I,J,typeCounter) = T-1
             
             ! Store fractional coverage in IUSE array for this grid cell.
             ! Units are [mil] for compatibility with legacy drydep code.
             State_Met%IUSE(I,J,typeCounter) = State_Met%LandTypeFrac(I,J,T) &
                                               * 1000

             ! If this type is water, set fraction land
             IF ( T .eq. 1 ) THEN
                State_Met%FRCLND(I,J) = 1.e+0_fp                          &
                                        - State_Met%LandTypeFrac(I,J,T)
             ENDIF

          ENDIF
       ENDDO

       ! Get IUSE type index with maximum coverage [mil]
       maxFracInd  = MAXLOC(State_Met%IUSE(I,J,1:State_Met%IREG(I,J)))

       ! Force IUSE to sum to 1000 by updating max value if necessary
       sumIUSE =  SUM(State_Met%IUSE(I,J,1:State_Met%IREG(I,J)))
       IF ( sumIUSE /= 1000 ) THEN
          State_Met%IUSE(I,J,maxFracInd) = State_Met%IUSE(I,J,maxFracInd) &
                                           + ( 1000 - sumIUSE )

       ENDIF

    ENDDO
    ENDDO

  END SUBROUTINE Compute_Soilmap_GCHP
!EOC
!#else
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_soilmap
!
! !DESCRIPTION: Subroutine INIT\_SOILMAP reads HadGEM2 soil map 
! information from disk (in netCDF format).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_SoilMap( am_I_Root, Input_Opt, RC )
!
! !USES:
!

    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE m_netcdf_io_open
    USE m_netcdf_io_read
    USE m_netcdf_io_readattr
    USE m_netcdf_io_close
    
    IMPLICIT NONE
    
#   include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!  
! !REVISION HISTORY:
!  09 Feb 2019 - J. Lam - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    !======================================================================
    ! Variable declarations
    !======================================================================
    
    ! Scalars
    INTEGER            :: I, J               ! Loop indices
    INTEGER            :: fId                ! netCDF file ID
    INTEGER            :: as                 ! Allocation status
    
    ! Character strings
    CHARACTER(LEN=255) :: nc_dir             ! netCDF directory name
    CHARACTER(LEN=255) :: nc_file            ! netCDF file name
    CHARACTER(LEN=255) :: nc_path            ! netCDF path name
    CHARACTER(LEN=255) :: v_name             ! netCDF variable name 
    CHARACTER(LEN=255) :: a_name             ! netCDF attribute name
    CHARACTER(LEN=255) :: a_val              ! netCDF attribute value
     
    ! Arrays for netCDF start and count values
    INTEGER            :: st1d(1), ct1d(1)   ! For 1D arrays    
    INTEGER            :: st2d(2), ct2d(2)   ! For 2D arrays 
     
    !======================================================================
    ! Initialize variables
    !======================================================================

    I_SOIL  = 192                                      ! # lons (0.625x0.9375)
    J_SOIL  = 145                                      ! # lats (0.625x0.9375)
    D_LON   = 0.9375e+0_fp                             ! Delta lon [degrees]
    D_LAT   = 0.625e+0_fp                              ! Delta lat [degrees]
    nc_file = 'HadGEM2ES_Soil_Ancil.nc'                ! Input file name

    ! Allocate arrays
    ALLOCATE( lon       ( I_SOIL         ), STAT=as )  
    ALLOCATE( lat       ( J_SOIL         ), STAT=as )
    ALLOCATE( THETA_WILT( I_SOIL, J_SOIL ), STAT=as ) 
    ALLOCATE( THETA_CRIT( I_SOIL, J_SOIL ), STAT=as )

    !======================================================================
    ! Open and read data from the netCDF file
    !======================================================================

    ! Construct file path from directory & file name
    nc_dir  = TRIM( Input_Opt%CHEM_INPUTS_DIR ) // 'Soil_Map/'
    nc_path = TRIM( nc_dir )                    // TRIM( nc_file )

    ! Open file for read
    CALL Ncop_Rd( fId, TRIM(nc_path) )
     
    ! Echo info to stdout
    IF ( am_I_Root ) THEN
       WRITE( 6, 100 ) REPEAT( '%', 79 )
       WRITE( 6, 110 ) TRIM(nc_file)
       WRITE( 6, 120 ) TRIM(nc_dir)
    ENDIF

    !----------------------------------------
    ! VARIABLE: lon
    !----------------------------------------
     
    ! Variable name
    v_name = "lon"
    
    ! Read lon from file
    st1d   = (/ 1      /)
    ct1d   = (/ I_SOIL /)
    CALL NcRd( lon, fId, TRIM(v_name), st1d, ct1d )
 
    ! Read the lon:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
    
    ! Echo info to stdout
    IF ( am_I_Root ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)     
    ENDIF

    !----------------------------------------
    ! VARIABLE: lat
    !----------------------------------------
    
    ! Variable name
    v_name = "lat"
    
    ! Read lat from file
    st1d   = (/ 1      /)
    ct1d   = (/ J_SOIL /)
    CALL NcRd( lat, fId, TRIM(v_name), st1d, ct1d )
     
    ! Read the lat:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
    
    ! Echo info to stdout
    IF ( am_I_Root ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val) 
    ENDIF

    !----------------------------------------
    ! VARIABLE: THETA_WILT
    !----------------------------------------

    ! Variable name
    v_name = "THETA_WILT"

    ! Read THETA_WILT from file
    st2d   = (/ 1, 1 /)
    ct2d   = (/ I_SOIL, J_SOIL /)
    CALL NcRd( THETA_WILT, fId, TRIM(v_name), st2d, ct2d )

    ! Read the THETA_WILT:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
    
    ! Echo info to stdout
    IF ( am_I_Root ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val) 
    ENDIF

    !----------------------------------------
    ! VARIABLE: THETA_CRIT
    !----------------------------------------

    ! Variable name
    v_name = "THETA_CRIT"

    ! Read THETA_CRIT from file
    st2d   = (/ 1, 1 /)
    ct2d   = (/ I_SOIL, J_SOIL /)
    CALL NcRd( THETA_CRIT, fId, TRIM(v_name), st2d, ct2d )

    ! Read the THETA_CRIT:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
    
    ! Echo info to stdout
    IF ( am_I_Root ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val) 
    ENDIF
    
    !=================================================================
    ! Cleanup and quit
    !=================================================================
    
    ! Close netCDF file
    CALL NcCl( fId )
    
    ! Echo info to stdout
    IF ( am_I_Root ) THEN
       WRITE( 6, 140 )
       WRITE( 6, 100 ) REPEAT( '%', 79 )
    ENDIF

    ! FORMAT statements
100 FORMAT( a                                              )
110 FORMAT( '%% Opening file  : ',         a               )
120 FORMAT( '%%  in directory : ',         a, / , '%%'     )
130 FORMAT( '%% Successfully read ',       a, ' [', a, ']' )
140 FORMAT( '%% Successfully closed file!'                 )

  END SUBROUTINE Init_SoilMap
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: compute_soilmap
!
! !DESCRIPTION: Subroutine COMPUTE\_SOILMAP computes the GEOS-Chem
!  arrays IREG, ILAND, IUSE (and corresponding 1-D arrays IJREG, IJLAND, 
!  IJUSE) on-the-fly from the Olson Land map file.  This routine, which is
!  intended to facilitate the Grid-Independent GEOS-Chem, replaces
!  the old rdland.F, which read from pre-computed "vegtype.global" files.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_SoilMap( am_I_Root, mapping, State_Met )
!
! !USES:
!
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN)    :: am_I_Root    ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MapWeight), POINTER       :: mapping(:,:) ! "fine" -> "coarse" mapping
    TYPE(MetState),  INTENT(INOUT) :: State_Met    ! Meteorology State object
!
! !REMARKS:
!  This routine supplies arrays that are required for legacy code routines:
!  (1) IREG,  ILAND,  IUSE are used by the Soil NOx routines
!  (2) IJREG, IJLAND, IJUSE are used by the dry deposition routines
! 
! !REVISION HISTORY: 
!  13 Mar 2012 - R. Yantosca - Initial version
!  19 Mar 2012 - R. Yantosca - Reorder ILAND, IUSE, IJLAND, IJUSE to be
!                              consistent w/ the leaf area indices
!  19 Mar 2012 - R. Yantosca - Compute the FRCLND array (from CMN_DEP_mod.F)
!  21 Mar 2012 - R. Yantosca - Now use REAL*4 for computation, to reduce
!                              roundoff errors at high-resolution
!  22 Mar 2012 - R. Yantosca - Now get surface area directly from variable
!                              A_CM2 (read from disk) instead of computing it
!  02 Apr 2012 - R. Yantosca - Now pass MAP (mapping weight object) via the
!                              arg list, to save the mapping info for later
!  09 Apr 2012 - R. Yantosca - Remove IJLOOP variable
!  09 Apr 2012 - R. Yantosca - Now do not compute IJREG, IJLAND, IJUSE; these
!                              are replaced by IREG, ILAND, IUSE arrays
!  17 Apr 2012 - R. Yantosca - Rename "map" object to "mapping" to avoid name
!                              confusion with an F90 intrinsic function
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  29 Nov 2012 - R. Yantosca - Added am_I_Root argument
!  12 Dec 2012 - R. Yantosca - Now get IREG, ILAND, IUSE from State_Met
!  20 Mar 2014 - R. Yantosca - Add shunts in lat & lon to reduce wall time
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL :: isGlobal
    INTEGER :: I,         J,         II,       III
    INTEGER :: JJ,        T,         N,        type
    INTEGER :: uniqOlson, sumIuse,   C,        IG
    REAL*4  :: xedge_w,   xedge_e,   yedge_s,  yedge_n
    REAL*4  :: xedgeC_w,  xedgeC_e,  yedgeC_s, yedgeC_n
    REAL*4  :: dxdy,      dxdy4,     mapWt,    area
    REAL*4  :: sumArea
    
    ! Generic arrays
    INTEGER :: maxIuse(1)
    
    ! Arrays on the Olson land map NATIVE GRID
    INTEGER :: indLon  (I_SOIL                     ) ! Index array for lons
    INTEGER :: shiftLon(I_SOIL                     ) ! Shifted indLon array
    REAL*4  :: lonedge (I_SOIL+1                   ) ! Lon edges   [degrees]
    REAL*4  :: latedge (          J_SOIL+1         ) ! Lat edges   [degrees]
    
    ! Arrays on the GEOS-CHEM GRID                 
    INTEGER :: ctOlson (IIPAR, JJPAR, 0:NSURFTYPE-1) ! Count of land types/box
    REAL*4  :: frOlson (IIPAR, JJPAR, 0:NSURFTYPE-1) ! Frac of land types/box
    INTEGER :: ordOlson(IIPAR, JJPAR, 0:NSURFTYPE-1) ! Order of land types

    ! Pointers
    INTEGER,  POINTER :: IREG(:,:)
    INTEGER,  POINTER :: ILAND(:,:,:)
    INTEGER,  POINTER :: IUSE(:,:,:)
    REAL(fp), POINTER :: FRCLND(:,:)
!
! !DEFINED PARAMETERS:
!
! The following parameters are used to skip over Olson NATIVE GRID boxes
! that are too far away from the GEOS-CHEM GRID BOX.  This can speed up
! the Olson computation by a factor of 100 or more!
!
#if defined( GRID05x0625 ) || defined( GRID025x03125 )
    REAL(fp), PARAMETER :: latThresh = 1e+0_fp   ! Lat threshold, nested grid
    REAL(fp), PARAMETER :: lonThresh = 1e+0_fp   ! Lon threshold, nested grid
#else
    REAL(fp), PARAMETER :: latThresh = 5e+0_fp   ! Lat threshold, global
    REAL(fp), PARAMETER :: lonThresh = 6e+0_fp   ! Lon threshold, global
#endif

    !======================================================================
    ! NATIVE GRID parameters (i.e. 0.5 x 0.5 "GENERIC")
    !======================================================================

    ! Be lazy, construct lon edges from lon centers
    DO I = 1, I_SOIL
       lonedge(I)      = DBLE( lon(I) ) - ( D_LON * 0.5e+0_fp )
       indLon(I)       = I
    ENDDO
    lonedge(I_SOIL+1) = lonedge(I_SOIL) + D_LON
    
    ! Be lazy, construct lat edges from lat centers
    DO J = 1, J_SOIL
       latedge(J)      = DBLE( lat(J) ) - ( D_LAT * 0.5e+0_fp )
    ENDDO
    latedge(J_SOIL+1) = latedge(J_SOIL) + D_LAT
    
    ! Shift longitudes by 2 degrees to the west for date-line handling
    shiftLon           = CSHIFT( indLon, -20 )

    !======================================================================
    ! Initialize variables outside of the main loop 
    !======================================================================

    ! Initialize pointers
    IREG    => State_Met%IREG
    ILAND   => State_Met%ILAND
    IUSE    => State_Met%IUSE
    FRCLND  => State_Met%FRCLND

    IREG     = 0
    ILAND    = 0
    IUSE     = 0
    FRCLND   = 1000e0
    ctOlson  = 0
    frOlson  = 0e0
    ordOlson = -999
    isGlobal = ( .not. ITS_A_NESTED_GRID() )

    !======================================================================
    ! Loop over all GEOS-CHEM GRID BOXES and initialize variables
    !======================================================================
    !$OMP PARALLEL DO                                                  &
    !$OMP DEFAULT( SHARED )                                            &
    !$OMP PRIVATE( I,        J,         xedgeC_w, yedgeC_s, xedgeC_e ) &
    !$OMP PRIVATE( yedgeC_n, dxdy4,     sumArea,  JJ,       III      ) &
    !$OMP PRIVATE( dxdy,     mapWt,     II,       xedge_w,  yedge_s  ) &
    !$OMP PRIVATE( xedge_e,  yedge_n,   area,     type,     maxIuse  ) &
    !$OMP PRIVATE( sumIUse,  uniqOlson, C,        IG                 ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Global lon index (needed for when running in ESMF)
       IG = I + I_LO - 1

       ! Edges of this GEOS-CHEM GRID box
       xedgeC_w  = GET_XEDGE( I,   J,   1 )          ! W edge
       yedgeC_s  = GET_YEDGE( I,   J,   1 )          ! S edge
       xedgeC_e  = GET_XEDGE( I+1, J,   1 )          ! E edge
       yedgeC_n  = GET_YEDGE( I,   J+1, 1 )          ! N edge
       
       ! "Area" of the GEOS-CHEM GRID box in degrees (DLON * DLAT)
       dxdy4     = ( xedgeC_e - xedgeC_w ) * ( yedgeC_n - yedgeC_s )
     
       ! Zero the summing array
       sumArea   = 0e0

       ! Reset counter of olson land types found per box
       uniqOlson = 0e0

       ! Counter for mapping object
       C         = 0

       !===================================================================
       ! Find each NATIVE GRID BOX that fits into the GEOS-CHEM GRID BOX.  
       ! Keep track of the land types and coverage fractions.
       !===================================================================

       ! Loop over latitudes on the NATIVE GRID
       DO JJ  = 1, J_SOIL

          ! Latitude edges of this NATIVE GRID box
          yedge_s    = latedge(JJ  )                ! S edge
          yedge_n    = latedge(JJ+1)                ! N edge

          !%%%%%% LATITUDE SHUNT TO REDUCE WALL TIME (bmy, 3/20/14) %%%%%%%%%%
          !%%%
          !%%% Skip further computations unless we are within LATTHRESH 
          !%%% degrees of the western edge of box (I,J).  This prevents 
          !%%% excess computations and subroutine calls.
          !%%%
          IF ( ABS( yedge_s - yedgeC_s ) > latThresh ) CYCLE
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          ! Loop over longitudes on the NATIVE GRID
          DO III = 1, I_SOIL

             ! Initialize
             dxdy       = 0e0
             mapWt      = 0e0
             
             ! Find the NATIVE GRID longitude index for use below.  Account for 
             ! the first GEOS-CHEM GRID box, which straddles the date line.
             IF ( isGlobal .and.  IG == 1 ) THEN
                II      = shiftLon(III)
             ELSE
                II      = indLon(III)
             ENDIF
             
             ! Edges of this NATIVE GRID box
             xedge_w    = lonedge(II  )                ! W edge
             xedge_e    = lonedge(II+1)                ! E edge

             ! Because the first GEOS-CHEM GRID BOX straddles the date line,
             ! we have to adjust the W and E edges of the NATIVE GRID BOX to
             ! be in monotonically increasing order.  This will prevent
             ! erronous results from being returned by GET_MAP_WT below.
             IF ( isGlobal .and. IG == 1 .and. II >= shiftLon(1) )  THEN
                xedge_w = xedge_w - 360e0
                xedge_e = xedge_e - 360e0
             ENDIF
         
             !%%%%%% LONGITUDE SHUNT TO REDUCE WALL TIME (bmy, 3/20/14) %%%%%%
             !%%%
             !%%% Skip further computations unless we are within LONTHRESH 
             !%%% degrees of the western edge of box (I,J).  This prevents 
             !%%% excess computations and subroutine calls.
             !%%%
             IF ( ABS( xedge_w - xedgeC_w ) > lonThresh ) CYCLE
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             ! "Area" of the NATIVE GRID BOX in degrees (DLON * DLAT)
             dxdy       = ( xedge_e - xedge_w ) * ( yedge_n - yedge_s )

             ! Get the mapping weight (i.e. The fraction of the NATIVE 
             ! GRID BOX that lies w/in the GEOS-CHEM GRID BOX)
             CALL GET_MAP_WT( xedge_w, xedge_e, xedgeC_w, xedgeC_e,  &
                              yedge_s, yedge_n, yedgeC_s, yedgeC_n,  &
                              mapWt                                 )

             ! Skip unless part (or all) of the NATIVE GRID BOX
             ! actually fits into the GEOS-CHEM GRID BOX
             IF ( mapWt <= 0e0 .or. mapWt > 1e0 ) CYCLE

             ! Area of the NATIVE GRID BOX that lies
             ! within the GEOS-CHEM GRID BOX
             area              = A_CM2(II,JJ,1) * mapWt
           
             ! Keep a total of the area
             sumArea           = sumArea + area

             ! Olson land map type on the NATIVE GRID
             type              = OLSON(II,JJ,1)
           
             ! Increment count of Olson types
             ctOlson(I,J,type) = ctOlson(I,J,type) + 1

             ! Add area covered by this olson type
             frOlson(I,J,type) = frOlson(I,J,type) + area

             ! Preserve ordering for backwards-compatibility w/ LAI data
             IF ( ordOlson(I,J,type) < 0 ) THEN 

                ! Counter of land types we have encountered for the first time
                uniqOlson          = uniqOlson + 1

                ! Record the order in which this land type was first encountered
                ordOlson(I,J,type) = uniqOlson

             ENDIF
             
             ! Save mapping information for later use in modis_lai_mod.F90
             ! in order to prepare the State_Met%XLAI array for use with the 
             ! legacy dry-deposition and soil NOx emissions codes.
             C                     = C + 1
             mapping(I,J)%count    = C
             mapping(I,J)%II(C)    = II
             mapping(I,J)%JJ(C)    = JJ
             mapping(I,J)%olson(C) = type
             mapping(I,J)%area(C)  = area
             mapping(I,J)%sumarea  = sumarea

          ENDDO
       ENDDO

       !===================================================================
       ! Construct GEOS-Chem type output arrays from the binning that we 
       ! just have completed.  Preserve the ordering from "vegtype.global"
       ! for backwards compatibility w/ existing code.
       !===================================================================
     
       ! Land type index for ILAND & IUSE
       maxIUse = 0

       ! Loop over all land types
       DO T = 0, NSURFTYPE-1

          ! Save the ordering of Olson land types for later use 
          ! by routines in the module modis_lai_mod.F90
          mapping(I,J)%ordOlson(T) = ordOlson(I,J,T)

          ! Normalize the land type coverage 
          frOlson(I,J,T)                =  &
               INT( ( ( frOlson(I,J,T) / sumArea ) * 1e3 ) + 0.5e0 )
 
          ! If land type T is represented in this box ...
          IF ( ctOlson(I,J,T) > 0 .and. ordOlson(I,J,T) > 0 ) THEN 
 
             ! Increment the count of Olson types in the box 
             IREG(I,J)                  = IREG(I,J) + 1
             
             ! Save land type into ILAND
             ILAND(I,J,ordOlson(I,J,T)) = T
             
             ! Save the fraction (in mils) of this land type
             IUSE(I,J,ordOlson(I,J,T))  = frOlson(I,J,T)

          ENDIF
       ENDDO

       ! Land type with the largest coverage in the GEOS-CHEM GRID BOX
       maxIuse = MAXLOC( IUSE( I, J, 1:IREG(I,J) ) )

       ! Sum of all land types in the GEOS-CHEM GRID BOX (should be 1000)
       sumIUse = SUM   ( IUSE( I, J, 1:IREG(I,J) ) )

       ! Make sure everything adds up to 1000.  If not, then adjust
       ! the land type w/ the largest coverage accordingly.
       ! This follows the algorithm from "regridh_lai.pro".
       IF ( sumIUse /= 1000 ) THEN
          IUSE(I,J,maxIUse) = IUSE(I,J,maxIUse) &
                            + ( 1000 - sumIUse )
       ENDIF
      
       ! Loop over land types in the GEOS-CHEM GRID BOX
       DO T = 1, IREG(I,J)

          ! If the current Olson land type is water (type 0),
          ! subtract the coverage fraction (IUSE) from FRCLND.
          IF ( ILAND(I,J,T) == 0 ) THEN
             FRCLND(I,J) = FRCLND(I,J)  - IUSE(I,J,T)
          ENDIF
       ENDDO

       ! Normalize FRCLND into the range of 0-1
       ! NOTE: Use REAL*4 for backwards compatibility w/ existing code!
       FRCLND(I,J) = FRCLND(I,J) / 1000e0

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO
  
!### Save code here for debugging
!###    do j = 1, jjpar
!###    do i = 1, iipar
!###       write( 800, '(2i5, f13.6)' ) i, j, frclnd(i,j)
!###       
!###       write( 810, '(25i4)'       ) i, j, ireg(i,j),                 &
!###                                    ( iland(i,j,t), t=1,ireg(i,j) ), &
!###                                    ( iuse (i,j,t), t=1,ireg(i,j) )
!###    enddo
!###    enddo
!###   
!###    ! ### DEBUG OUTPUT
!###    C = mapping(23,34)%count
!###    print*, '### count   : ', C
!###    print*, '### II      : ', mapping(23,34)%II(1:C)
!###    print*, '### JJ      : ', mapping(23,34)%JJ(1:C)
!###    print*, '### area    : ', mapping(23,34)%area(1:C)
!###    print*, '### sumarea : ', mapping(23,34)%sumarea

  END SUBROUTINE Compute_SoilMap
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_soilmap
!
! !DESCRIPTION: Subroutine CLEANUP\_SOILMAP deallocates all allocated
!  global module variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_SoilMap( am_I_Root )
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN) :: am_I_Root   ! Are we on the root CPU?
!
! !REVISION HISTORY:'
!  09 Feb 2019 - J. Lam - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( ALLOCATED( lon        ) ) DEALLOCATE( lon        )
    IF ( ALLOCATED( lat        ) ) DEALLOCATE( lat        )
    IF ( ALLOCATED( THETA_WILT ) ) DEALLOCATE( THETA_WILT )
    IF ( ALLOCATED( THETA_CRIT ) ) DEALLOCATE( THETA_CRIT )

  END SUBROUTINE Cleanup_SoilMap
!EOC
!#endif
END MODULE HadGEM2_SoilMap_Mod

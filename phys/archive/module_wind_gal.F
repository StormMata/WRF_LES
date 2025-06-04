!WRF:MODEL_LAYER:PHYSICS
MODULE module_wind_gal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This is the main module for generalized actuator line (GAL) scheme. !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================================================================
! nampelist.input options:
!   &physics
!   wind_wtp_opt(max_domains)           = 1,        ! 0 = deactivates GAL/GAD/GADrs/SAD models
!                                                     1 = activates GAL model
!                                                     2 = activates GAD model
!                                                     3 = activates GADrs model
!                                                     4 = activates SAD model
!   wind_wtp_spec_rho                   = -1.0,     ! <=0 air density calculated
!                                                         by the model
!                                                     >0  user-specified air density
!   wind_wtp_spec_turbdir               = -1.0,     ! <0  direction the wind turbine
!                                                         is facing (calculated by 
!                                                         the model)
!                                                     >=0 user-specified direction
!   wind_wtp_bladesmear_opt             = 0,        ! 0 = Gaussian smearing function
!                                                     1 = Polynomial smearing function
!                                                     2 = Generalized Gaussian smearing
!                                                         function
!                                                     3 = Chord-thickness Gaussian
!                                                         smearing function
!   wind_wtp_nacellesmear_opt           = 0,        ! 0 = Gaussian smearing function
!                                                     1 = Polynomial smearing function
!   wind_wtp_towersmear_opt             = 0,        ! 0 = Gaussian smearing function
!                                                     1 = Polynomial smearing function
!   wind_wtp_aerodynloss_opt            = 3,        ! 0 = no aerodynamic loss model activated
!                                                     1 = only hub loss model from Prandtl
!                                                     2 = only tip loss model from Prandtl
!                                                     3 = both hub and tip loss models from Prandtl
!                                                     4 = tip loss model from Shen et al. (2015)
!   wind_wtp_correction3D_opt           = 0,        ! 0 = no three-dimensional correction
!                                                         for airfoils
!                                                     1 = Snel et al. (1993) correction for cl
!                                                         - inviscid cl based on thin airfoil theory
!                                                     2 = Lindenburg (2004) correction for cl
!                                                         - inviscid cl based on thin airfoil theory
!                                                     3 = Du & Selig (1998) correction for cl and cd
!                                                         - inviscid cl based on thin airfoil theory
!                                                     4 = Du & Selig (1998) correction for cl and cd
!                                                         - inviscid cl based on Branlard (2011) approximation
!                                                     5 = Chaviaropoulos and Hansen (2000) correction for cl and cd
!                                                         - inviscid cl based on thin airfoil theory
!   wind_wtp_printdiag_opt              = 0,        ! 0 = no diagnostics printed
!                                                     1 = only necessary information
!                                                     2 = full diagnostics printed on
!                                                         processor(s) where turbine(s)
!                                                         is/are located
!   wind_wtp_yaw_opt                    = .true.,   ! turbine yaw mechanism is turned on
!   wind_wtp_hubpar_opt                 = .false.,  ! no hub parameterization 
!   wind_wtp_towpar_opt                 = .false.,  ! no tower parameterization
!   wind_wtp_errorcheck_opt             = .false.,  ! no fatal error check for GAD/GAL
!   wind_wtp_speedcontrol_opt           = .false.,  ! no speed controller
!   wind_wtp_pitchcontrol_opt           = .false.,  ! no pitch controller
!   wind_wtp_nTurbines                  = 1,        ! number of wind turbines
!   wind_wtp_nSections                  = 3,        ! number of sections
!   wind_wtp_nElements                  = 64,       ! number of blade-elements
!   wind_wtp_nSeries                    = 192,      ! length of storage matrix,
!                                                      wind_wtp_nTurbines &
!                                                     *wind_wtp_nSections &
!                                                     *wind_wtp_nElements
!   wind_wtp_nTowpoints                 = 40,       ! number of tower actuator points
!   wind_wtp_nSeries_tow                = 40,       ! length of storage matrix,
!                                                      wind_wtp_nTurbines &
!                                                     *wind_wtp_nTowpoints
!   wind_wtp_nNacpoints                 = 10,       ! number of nacelle actuator points
!   wind_wtp_nSeries_nac                = 10,       ! length of storage matrix,
!                                                      wind_wtp_nTurbines &
!                                                     *wind_wtp_nNacpoints
!   wind_wtp_nAvesamples                = 100,      ! number of samples for running average
!   wind_wtp_nSeries_tme                = 100,      ! length of storage matrix,
!                                                      wind_wtp_nTurbines &
!                                                     *wind_wtp_nAvesamples
!   windturbines_latlon                 = 0,        ! 0 = turbine location in SI units
!                                                         (in meters), used for ideal cases
!                                                     1 = turbine location in lat and lon
!                                                         format, used for real cases
!   &dynamics
!   wind_wtp_tsout_opt                  = .false.,  ! time series output for wind turbine
!                                                     variables is turned off
!   wind_wtp_ts_buf_size                = 1,        ! buffer size for time series output
!                                                     e.g., 1=every timestep, 100=every
!                                                     100 timesteps
!   wind_wtp_pln_tsout_opt              = .false.,  ! Switch to turn on time series output
!                                                     for three velocity components at the
!                                                     specified plane location
!   wind_wtp_pln_ts_buf_size            = 1,        ! buffer size for time series output
!                                                     e.g., 1=every timestep, 100=every
!                                                     100 timesteps
!   wind_wtp_pln_nTurbines              = 1,        ! total number of wind turbines for which
!                                                     time series of velocity components to be
!                                                     calculated
!   wind_wtp_pln_turbine_ids            = 1,        ! wind turbine id for which time series
!                                                     of velocity components to be calculated
!   wind_wtp_pln_nz                     = 20,       ! number of plane-element ponts in the
!                                                     vertical direction
!   wind_wtp_pln_nxy                    = 20,       ! number of plane-element points in the
!                                                     horizontal space
!   wind_wtp_pln_nSeries_tme            = 400,    ! length of storage matrix,
!                                                      wind_wtp_pln_nTurbines &
!                                                     *wind_wtp_pln_nz &
!                                                     *wind_wtp_pln_nxy
!   wind_wtp_pln_width                  = 300.0,    ! plane width in meter
!   wind_wtp_pln_height                 = 300.0,    ! plane height in meter
!   wind_wtp_pln_dist                   = 100.0,    ! distance between reference plane and
!                                                     corresponding wind turbine
!=======================================================================
! Note: If you want to operate the wind turbine at a constant rotational
!       speed or blade pitch angle, set the second or third column in
!       the windTurbineData/operationalData.tbl to the same value for
!       all wind speeds.
!=======================================================================
! BROADCASTING
!=======================================================================
! IF nprocX .GT. 3 and commSize .GT. 10 THEN
!   IF sender .NE. root ---> broadcast to 10 processors including root.
!   IF sender .EQ. root ---> broadcast to 9 processors without root.
! IF commSize .LE. 8 THEN
!                       ---> broadcast to MPI_COMM_WORLD.
! ENDIF
!=======================================================================
!     EXAMPLE   : A new communicator is created for a single turbine using
!                 DMPAR decomposition and 16 processors.
!      o        : Inflow wind location
!      T        : Stand-alone turbine location
!      P0-15    : Available processors
!      (1)-(9)  : Grouped processors within the new communicator
!      P0       : Root processor (also involved in the communication)  
!      comm_new : Includes ---> P0 + (P5, P6, P7, P9, P10, P11, P13, P14, P15)
!      #####
!      #   #    : Area needed for the GAD/GAL computations including turbine yaw
!      #####
!      +------------+------------+------------+------------+
!      |P12         |P13         |         P14|         P15|
!      |            |            |            |            |
!      |            |            |            |            |
!      |            |   ###########################        |
!      |            |   #(7)     |(8)         |(9)#        |
!      +------------+---#--------+------------+---#--------+ 
!      |P8          |P9 #    o   |P10         |   #     P11|
!      |            |   #      * |      /     |   #        |
!      |            |   #        *     /      |   #        |
!      |            |   #        | *  /       |   #        |
!      |            |   #(4)     |   T     (5)|(6)#        |
!      +------------+---#--------+--/---------+---#--------+
!      |P4          |P5 #        | /        P6|   #      P7|
!      |            |   #        |/           |   #        |
!      |            |   #        |            |   #        |
!      |            |   #        |            |   #        |
!      |            |   #(1)     |(2)         |(3)#        |
!      +------------+---###########################--------+
!      |P0          |P1          |          P2|          P3|
!      |            |            |            |            |
!      |            |            |            |            |
!      |            |            |            |            |
!      |            |            |            |            |
!      +------------+------------+------------+------------+
!      |<---------------------nprocX--------------------->|
!      REMARKS:
!       (1) This procedure is repeated for each turbine.
!       (2) If turbine is located at P0, "newComm" includes only
!           processors P1, P4, P5 and broadcasts to 3 processors.
!       (3) If turbine is located within one of the corner processors 
!           except P0 (i.e., P3, P12, or P15), "newComm" includes
!           contiguous processors + P0 and broadcasts to 4 processors.
!       (4) If turbine is located within one of the "near-wall" processors 
!           except from corners (i.e., P1, P2, P4, P7, P8, P11, P13, P14),
!           "newComm" includes contiguous processors + P0 and broadcasts
!           to 6 processors. For instance, if a turbine is located at P4,
!           "newComm" will include processors P0, P1, P4, P5, P8, and P9.
!=======================================================================
#ifdef DM_PARALLEL
  USE module_dm
#endif
  USE module_llxy
  USE module_domain, ONLY : domain
  USE module_configure, ONLY : model_config_rec, grid_config_rec_type
  USE module_model_constants, ONLY : degrad
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.d0)
  DOUBLE PRECISION, PARAMETER :: piconst = 4.0_dp*ATAN(1.0_dp)
  DOUBLE PRECISION, PARAMETER :: rpmradsec = (2.0_dp*piconst)/60.0_dp
!
  INTEGER, PARAMETER :: smoothing = CEILING( SQRT( LOG(1.0/0.001) ) )
#ifdef DM_PARALLEL
  INTEGER, PARAMETER :: root = 0 ! Specifies the root processor
#endif
  INTEGER, PARAMETER :: maxVals1 = 1E2 ! Maximum allowed vector size
  INTEGER, PARAMETER :: maxVals2 = 1E3 ! Maximum allowed vector size
!
  REAL :: specRho, specDir, angle
  REAL :: thetacOut, omegaOut, aeroPowOut, aeroThrOut
  REAL :: planeWidth, planeHeight, planeDistance
!
  INTEGER :: bladeSmearOpt, nacelleSmearOpt, towerSmearOpt
  INTEGER :: aerodynLossOpt, correction3DOpt, printDiagOpt
  INTEGER :: nTurbineTypes, nTurbines, nElements, nSections
  INTEGER :: nTowPoints, nNacPoints, nAveSamples, wtID
  INTEGER :: ils, ile, jls, jle, kls, kle
!
  INTEGER :: nPlnWT, nPlnWTcounter, nHorizPoints, nVertPoints
#ifdef DM_PARALLEL
  INTEGER :: ierr, rank, commSize, worldGroup
  INTEGER :: procsW, procsE, procsN, procsS
  INTEGER :: procsNW, procsNE, procsSW, procsSE
#endif
  LOGICAL :: yawOpt, hubParOpt, towParOpt, errorCheckOpt, tsoutOpt
  LOGICAL :: restartController, restartController1
!
  LOGICAL :: tsoutPlnOpt
!
  REAL, DIMENSION(:), ALLOCATABLE :: solidity, projectionRadius, nacelleLength, contribution
  REAL, DIMENSION(:), ALLOCATABLE :: epsilonMin, epsilonMin2, polA, polB, epsilonFactor
  REAL, DIMENSION(:), ALLOCATABLE :: tiltAngle, preconeAngle, hubDiameter, dzConst
  REAL, DIMENSION(:), ALLOCATABLE :: towerDiameter, rotorDiameter, mechEfficiency, overhang
  REAL, DIMENSION(:), ALLOCATABLE :: rTip, rHub, rTow, dr, towerDs, nacelleDs, undSling
  REAL, DIMENSION(:), ALLOCATABLE :: inflowLocation, yawRate, thresholdYaw, cutIn, cutOut
  REAL, DIMENSION(:), ALLOCATABLE :: initGenTorque, ratedRotSpeed
  REAL, DIMENSION(:), ALLOCATABLE :: GBEfficiency, GBRatio, ratedGenPow
  REAL, DIMENSION(:), ALLOCATABLE :: genIner, hubIner, bladeIner, genEfficiency 
  REAL, DIMENSION(:), ALLOCATABLE :: rateLimitGenTorque, cutInGenSpeed, region2StartGenSpeed
  REAL, DIMENSION(:), ALLOCATABLE :: region2EndGenSpeed, cutInGenTorque, ratedGenTorque
  REAL, DIMENSION(:), ALLOCATABLE :: maxGenTorque, KGen, speedFilterCornerFrequency, pitchRate
  REAL, DIMENSION(:), ALLOCATABLE :: rateLimitBladePitch, pitchMin, pitchMax, pitchK
  REAL, DIMENSION(:), ALLOCATABLE :: pitchControlKP, pitchControlKI, pitchControlKD
  REAL, DIMENSION(:), ALLOCATABLE :: nacDragCoeff, towDragCoeff
  REAL, DIMENSION(:), ALLOCATABLE :: nacelleEpsilon, nacelleCepsilon
  REAL, DIMENSION(:), ALLOCATABLE :: towerEpsilon, towerCepsilon
  REAL, DIMENSION(:), ALLOCATABLE :: uvShaftDir, lat, lon
  REAL, DIMENSION(:), ALLOCATABLE :: xMinTile, xMaxTile, yMinTile, yMaxTile
  REAL, DIMENSION(:), ALLOCATABLE :: xMinTrb, xMaxTrb, yMinTrb, yMaxTrb
  REAL, DIMENSION(:), ALLOCATABLE :: hubHeight, hubHeightASL, hubHeightASLLocal
  REAL, DIMENSION(:), ALLOCATABLE :: towerToRotorApex, undSlingRotor
  REAL, DIMENSION(:), ALLOCATABLE :: windDirFstAveLocal, uFstAveLocal
  REAL, DIMENSION(:), ALLOCATABLE :: vFstAveLocal, wFstAveLocal
  REAL, DIMENSION(:), ALLOCATABLE :: rotSpeedLocal, rotSpeedF, driveTrainIner
  REAL, DIMENSION(:), ALLOCATABLE :: speedErrorInit, trqGenInit, GK
  REAL, DIMENSION(:), ALLOCATABLE :: rotSpeedInit, pitchInit
  REAL, DIMENSION(:), ALLOCATABLE :: speedErrorLocal, intSpeedErrorLocal
  REAL, DIMENSION(:), ALLOCATABLE :: speedErrorLast, derivSpeedError
  REAL, DIMENSION(:), ALLOCATABLE :: pitchP, pitchI, pitchD
  REAL, DIMENSION(:), ALLOCATABLE :: pitchCommanded, pitchRateCommanded
  REAL, DIMENSION(:), ALLOCATABLE :: deltaAzimuthLocal, psiOldLocal
  REAL, DIMENSION(:), ALLOCATABLE :: azimuthLocal, capPhiLocal, delta, beta
  REAL, DIMENSION(:), ALLOCATABLE :: rhoFstAveLocal, trbDirLocal
  REAL, DIMENSION(:), ALLOCATABLE :: wsRotorLocal, v0FstLocal
  REAL, DIMENSION(:), ALLOCATABLE :: eYawLocal, deltaCapPhiLocal
  REAL, DIMENSION(:), ALLOCATABLE :: tsrSavedLocal, trqGenCommandedOldLocal
  REAL, DIMENSION(:), ALLOCATABLE :: thrSavedLocal, trqSavedLocal
  REAL, DIMENSION(:), ALLOCATABLE :: powSavedLocal, mechPowSavedLocal
  REAL, DIMENSION(:), ALLOCATABLE :: powCoefSavedLocal, thrCoefSavedLocal
  REAL, DIMENSION(:), ALLOCATABLE :: trqCoefSavedLocal, genPowSavedLocal
  REAL, DIMENSION(:), ALLOCATABLE :: genSpeedOldLocal, genSpeedFOldLocal
  REAL, DIMENSION(:), ALLOCATABLE :: pitchSavedLocal, pitchOld, pitchLocal
  REAL, DIMENSION(:), ALLOCATABLE :: trqGenCommandedOldOld, genSpeedFOldOld
  REAL, DIMENSION(:), ALLOCATABLE :: Lx, Ly, Lz, deltaYawLocal, deltaYaw
  REAL, DIMENSION(:), ALLOCATABLE :: psiOld, azimuthOld, deltaAzimuthOld
!
  INTEGER, DIMENSION(:), ALLOCATABLE :: specTrb
  INTEGER, DIMENSION(:), ALLOCATABLE :: nBlade, kTrbBot, kTrbHub, kTrbTop
  INTEGER, DIMENSION(:), ALLOCATABLE :: kHubBot, kHubTop, kTowBot, kTowTop
  INTEGER, DIMENSION(:), ALLOCATABLE :: xLocRot, xLocHub, xLocTow
  INTEGER, DIMENSION(:), ALLOCATABLE :: yLocRot, yLocHub, yLocTow
  INTEGER, DIMENSION(:), ALLOCATABLE :: operDataVectLen, trbDataVectLen
  INTEGER, DIMENSION(:), ALLOCATABLE :: iTow, jTow, iRot, jRot
  INTEGER, DIMENSION(:), ALLOCATABLE :: iFstInit, jFstInit, kFst
  INTEGER, DIMENSION(:), ALLOCATABLE :: iFstLocal, jFstLocal
  INTEGER, DIMENSION(:), ALLOCATABLE :: turbineTypeID, turbineTypeVect
#ifdef DM_PARALLEL
  INTEGER, DIMENSION(:), ALLOCATABLE :: sender, procID, procIDFst, procIDFstInit
  INTEGER, DIMENSION(:), ALLOCATABLE :: proc, procFst, ProcFstInit
  INTEGER, DIMENSION(:), ALLOCATABLE :: procInit, procIDInit, newGroup, newComm
#endif
! 
  INTEGER, DIMENSION(:, :), ALLOCATABLE :: foilType, foilPropVectLen
!
  LOGICAL, DIMENSION(:), ALLOCATABLE :: speedCtrlOpt, pitchCtrlOpt
  LOGICAL, DIMENSION(:), ALLOCATABLE :: speedControlFileExists, pitchControlFileExists
  LOGICAL, DIMENSION(:), ALLOCATABLE :: doYaw, doYawLocal
  LOGICAL, DIMENSION(:), ALLOCATABLE :: initialYaw, initialYawLocal
  LOGICAL, DIMENSION(:), ALLOCATABLE :: reOperate, initialYawDone
!
  REAL, DIMENSION(:, :), ALLOCATABLE :: c, t, thetaS
  REAL, DIMENSION(:, :), ALLOCATABLE :: v1Vect, thetapVect, omegaVect
  REAL, DIMENSION(:, :), ALLOCATABLE :: aeroPowVect, aeroThrVect
  REAL, DIMENSION(:, :), ALLOCATABLE :: bladeEpsilon, bladeCepsilon
  REAL, DIMENSION(:, :), ALLOCATABLE :: bladeRadius, rVect, chordVect
  REAL, DIMENSION(:, :), ALLOCATABLE :: thicknessVect, thetaGVect
  REAL, DIMENSION(:, :), ALLOCATABLE :: bladeRoot, rotorApexInit, posRotorApex
  REAL, DIMENSION(:, :), ALLOCATABLE :: rotorApexLocal, towerBaseLocation
  REAL, DIMENSION(:, :), ALLOCATABLE :: uvShaft, uvTower, towerShaftIntersect
  REAL, DIMENSION(:, :), ALLOCATABLE :: rhoFstTSLocal, uFstTSLocal
  REAL, DIMENSION(:, :), ALLOCATABLE :: vFstTSLocal, wFstTSLocal, windDirFstTSLocal
!
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: uPlaneSavedLocal
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: vPlaneSavedLocal
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: wPlaneSavedLocal
!
#ifdef DM_PARALLEL
  INTEGER, DIMENSION(:, :), ALLOCATABLE :: procs, procsComm
#endif
!
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: alphaVect, cdVect, clVect, cmVect
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: fnSavedLocal, ftSavedLocal, FSavedLocal
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: phiSavedLocal, alphaSavedLocal
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: clSavedLocal, cdSavedLocal, cmSavedLocal
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: LSavedLocal, DSavedLocal, v1SavedLocal 
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: vTSavedLocal, vRelSavedLocal
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: vThetaSavedLocal, vZSavedLocal
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: uSavedLocal, vSavedLocal, wSavedLocal
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: towerPoints, towerPointsLocal
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: nacellePoints, nacellePointsLocal
!
  REAL, DIMENSION(:, :, :, :), ALLOCATABLE :: planePoints, planePointsLocal
!
  REAL, DIMENSION(:, :, :, :), ALLOCATABLE :: bladePoints, bladePointsLocal
CONTAINS
  SUBROUTINE calc_forces_gal(                                  &
        &  dx, dy, rdzw, dt, muu, muv, mut, ht, ru_tendf       &
        &, rv_tendf, rw_tendf, u_2, v_2, w_2                   &
        &, ids, ide, jds, jde, kds, kde                        &
        &, ims, ime, jms, jme, kms, kme                        &
        &, its, ite, jts, jte, kts, kte                        &
        &, itimestep, restart, iFst, jFst                      &
        &, windDirFstAve, rhoFstAve, v0Fst, rhoFstTS, uFstTS   &
        &, vFstTS, wFstTS, windDirFstTS, wsRotor, trbDirSaved  &
        &, deltaAzimuth, capPhi, eYaw, rho, tsrSaved, thrSaved &
        &, trqSaved, powSaved, mechPowSaved, genPowSaved       &
        &, thrCoeffSaved, trqCoeffSaved, powCoeffSaved         &
        &, FnSaved, FtSaved, FSaved, phiSaved, alphaSaved      &
        &, clSaved, cdSaved, cmSaved, LSaved, DSaved, v1Saved  &
        &, vThetaSaved, vTSaved, vRelSaved, vZSaved, uSaved    &
        &, vSaved, wSaved, rotSpeed, psi, azimuth              &
        &, trqGenCommandedOld, genSpeedOld, genSpeedFOld       &
        &, pitch, speedError, intSpeedError, rotorApexX        &
        &, rotorApexY, rotorApexZ, bladePointsX, bladePointsY  &
        &, bladePointsZ, towerPointsX, towerPointsY            &
        &, towerPointsZ, nacellePointsX, nacellePointsY        &
        &, nacellePointsZ, wtp_history_flag, wtp_restart_flag  &
        &, readFromRSTFile, uPlaneSaved, vPlaneSaved           &
        &, wPlaneSaved, planePointsX, planePointsY             &
        &, planePointsZ                                        &
        &)
  IMPLICIT NONE
#if defined(DM_PARALLEL) && !defined(STUBMPI)
  INCLUDE "mpif.h"
#endif
    REAL, INTENT(IN) :: dx, dy, dt
!
    INTEGER, INTENT(IN) :: ids, ide, jds, jde, kds, kde
    INTEGER, INTENT(IN) :: its, ite, jts, jte, kts, kte
    INTEGER, INTENT(IN) :: ims, ime, jms, jme, kms, kme
    INTEGER, INTENT(IN) :: itimestep
!
    LOGICAL, INTENT(IN) :: restart, wtp_history_flag, wtp_restart_flag
!
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: muu, muv, mut, ht
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: rdzw, rho, u_2, v_2, w_2
!
    INTEGER, INTENT(INOUT) :: readFromRSTFile
!
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT) :: ru_tendf, rv_tendf, rw_tendf
!
    REAL, DIMENSION(1:nTurbines), INTENT(OUT) :: windDirFstAve, rhoFstAve, v0Fst, capPhi
    REAL, DIMENSION(1:nTurbines), INTENT(OUT) :: trbDirSaved, thrSaved, wsRotor
    REAL, DIMENSION(1:nTurbines), INTENT(OUT) :: trqSaved, powSaved, mechPowSaved
    REAL, DIMENSION(1:nTurbines), INTENT(OUT) :: genPowSaved, thrCoeffSaved
    REAL, DIMENSION(1:nTurbines), INTENT(OUT) :: trqCoeffSaved, powCoeffSaved
    REAL, DIMENSION(1:nTurbines), INTENT(INOUT) :: deltaAzimuth, tsrSaved
    REAL, DIMENSION(1:nTurbines), INTENT(INOUT) :: genSpeedOld, genSpeedFOld
    REAL, DIMENSION(1:nTurbines), INTENT(INOUT) :: speedError, intSpeedError
    REAL, DIMENSION(1:nTurbines), INTENT(INOUT) :: pitch, psi, azimuth
    REAL, DIMENSION(1:nTurbines), INTENT(INOUT) :: eYaw, rotSpeed, trqGenCommandedOld
    REAL, DIMENSION(1:nTurbines), INTENT(INOUT) :: rotorApexX, rotorApexY, rotorApexZ
!
    INTEGER, DIMENSION(1:nTurbines), INTENT(INOUT) :: iFst, jFst
!
    REAL, DIMENSION(1:nTurbines, 1:nTowPoints), INTENT(OUT) :: towerPointsX
    REAL, DIMENSION(1:nTurbines, 1:nTowPoints), INTENT(OUT) :: towerPointsY
    REAL, DIMENSION(1:nTurbines, 1:nTowPoints), INTENT(OUT) :: towerPointsZ
!
    REAL, DIMENSION(1:nTurbines, 1:nNacPoints), INTENT(OUT) :: nacellePointsX
    REAL, DIMENSION(1:nTurbines, 1:nNacPoints), INTENT(OUT) :: nacellePointsY
    REAL, DIMENSION(1:nTurbines, 1:nNacPoints), INTENT(OUT) :: nacellePointsZ
!
    REAL, DIMENSION(1:nTurbines, 1:nAveSamples), INTENT(INOUT) :: rhoFstTS, windDirFstTS
    REAL, DIMENSION(1:nTurbines, 1:nAveSamples), INTENT(INOUT) :: uFstTS, vFstTS, wFstTS
!
    REAL, DIMENSION(1:nTurbines, 1:nSections, 1:nElements), INTENT(OUT) :: FnSaved, FtSaved, FSaved
    REAL, DIMENSION(1:nTurbines, 1:nSections, 1:nElements), INTENT(OUT) :: phiSaved, alphaSaved
    REAL, DIMENSION(1:nTurbines, 1:nSections, 1:nElements), INTENT(OUT) :: clSaved, cdSaved, cmSaved
    REAL, DIMENSION(1:nTurbines, 1:nSections, 1:nElements), INTENT(OUT) :: LSaved, DSaved, v1Saved
    REAL, DIMENSION(1:nTurbines, 1:nSections, 1:nElements), INTENT(OUT) :: vThetaSaved, vTSaved
    REAL, DIMENSION(1:nTurbines, 1:nSections, 1:nElements), INTENT(OUT) :: vRelSaved, vZSaved
    REAL, DIMENSION(1:nTurbines, 1:nSections, 1:nElements), INTENT(OUT) :: uSaved, vSaved, wSaved
    REAL, DIMENSION(1:nTurbines, 1:nSections, 1:nElements), INTENT(OUT) :: bladePointsX
    REAL, DIMENSION(1:nTurbines, 1:nSections, 1:nElements), INTENT(OUT) :: bladePointsY
    REAL, DIMENSION(1:nTurbines, 1:nSections, 1:nElements), INTENT(OUT) :: bladePointsZ
!
    REAL, DIMENSION(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints), INTENT(OUT) :: uPlaneSaved
    REAL, DIMENSION(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints), INTENT(OUT) :: vPlaneSaved
    REAL, DIMENSION(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints), INTENT(OUT) :: wPlaneSaved
    REAL, DIMENSION(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints), INTENT(OUT) :: planePointsX
    REAL, DIMENSION(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints), INTENT(OUT) :: planePointsY
    REAL, DIMENSION(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints), INTENT(OUT) :: planePointsZ
!
    CHARACTER*256  :: outstring
! ... Local
    REAL, PARAMETER :: thresholdNOTYawing = 0.001 ! [deg]
    REAL, PARAMETER :: thresholdCapPhi = 0.100 ! [deg]
    REAL, PARAMETER :: alphaMin = -5.000 ! alphaMin is the angle of attack where linear portion of lift curve slope begins [deg]
    REAL, PARAMETER :: alphaMax = 5.000 ! alphaMax is the angle of attack where linear portion of lift curve slope ends [deg]
    REAL, PARAMETER :: alphaThreshold = 30.000 ! maximum angle of attack to apply full correction [deg]
!
    REAL :: nacelleFrontalArea, clOut, cdOut, cmOut
    REAL :: projectionDistanceU, projectionDistanceV, projectionDistanceW
    REAL :: etaEpsilonU3D, etaEpsilonV3D, etaEpsilonW3D
    REAL :: dist, projectionLimit, coeffU3D, coeffV3D, coeffW3D
    REAL :: lowPassCoeff, region2StartGenTorque
    REAL :: region2EndGenTorque, trqSlope, ratedGenSpeed
    REAL :: region4EndGenTorque
    REAL :: fnSum, ftSum, FSum, phiSum, alphaSum
    REAL :: clSum, cdSum, cmSum, LSum, DSum
    REAL :: v1Sum, vTSum, vRelSum, vThetaSum, vZSum
    REAL :: d1, d2, d3, clMin, clMax
    REAL :: alphaZero, cdZero, clInv, deltaCd, deltaCl, slope, adj
    REAL :: rhoFstSum, uFstSum
    REAL :: vFstSum, wFstSum, windDirFstSum
    REAL :: sinFstSum, cosFstSum
!
    INTEGER :: i, j, k, n, e, m, kt, nt
    INTEGER :: itf, jtf, ktf, ii, jj, kk, ilev, jlev, klev
    INTEGER :: indexStart, indexEnd
!
    LOGICAL :: pitch_sw
!
    INTEGER, DIMENSION(1:nTurbines) :: iRotLocal, jRotLocal, iSmear, jSmear
!
    INTEGER, DIMENSION(1:nTurbineTypes) :: kEnd
!
    ! Rotation matrices: 
    REAL, DIMENSION(1:3) :: velocityFstXYZ, rotVelocityFstXYZ
    REAL, DIMENSION(1:3) :: velocityXYZ, rotVelocityXYZ
    REAL, DIMENSION(1:3) :: posXYZ, rotPosXYZ
    REAL, DIMENSION(1:3) :: forceXYZ, rotForceXYZ, vectZero
    REAL, DIMENSION(1:3) :: disVector, dir1, dir2, dir3
!
    REAL, DIMENSION(1:nTurbineTypes) :: zMaxTrb
!
    REAL, DIMENSION(1:nTurbines) :: rhoFst, windDirFst
    REAL, DIMENSION(1:nTurbines) :: uAtRotorLocal, vAtRotorLocal, wAtRotorLocal
    REAL, DIMENSION(1:nTurbines) :: uAtRotor, vAtRotor, wAtRotor
    REAL, DIMENSION(1:nTurbines) :: uFst, vFst, wFst
    REAL, DIMENSION(1:nTurbines) :: uFstLocal, vFstLocal, wFstLocal, rhoFstLocal
    REAL, DIMENSION(1:nTurbines) :: inflowLocation_x, inflowLocation_y, inflowLocation_z
    REAL, DIMENSION(1:nTurbines) :: genSpeed, genSpeedF, dGenSpeed
    REAL, DIMENSION(1:nTurbines) :: trqGenCommanded, trqGenRateCommanded
    REAL, DIMENSION(1:nTurbines) :: thrAcc, trqAcc, deltaCapPhi
    REAL, DIMENSION(1:nTurbines) :: powGenerator, yawdir
!
    REAL, DIMENSION(1:nNacPoints) :: vNacX, vNacY, vNacZ, vNac, capPhiAtNac
    REAL, DIMENSION(1:nNacPoints) :: fnNacelle, fxNacelle, fyNacelle, fzNacelle
!
    REAL, DIMENSION(1:nTowPoints) :: vTowX, vTowY, vTow, capPhiAtTow
    REAL, DIMENSION(1:nTowPoints) :: fnTower, fxTower, fyTower, fzTower
!
    REAL, DIMENSION(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints) :: uPlaneLocal
    REAL, DIMENSION(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints) :: vPlaneLocal
    REAL, DIMENSION(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints) :: wPlaneLocal
!
    REAL, DIMENSION(1:nTurbines, 1:nElements) :: fnSectAveLocal, ftSectAveLocal
    REAL, DIMENSION(1:nTurbines, 1:nElements) :: FSectAveLocal, phiSectAveLocal
    REAL, DIMENSION(1:nTurbines, 1:nElements) :: alphaSectAveLocal, clSectAveLocal
    REAL, DIMENSION(1:nTurbines, 1:nElements) :: cdSectAveLocal, cmSectAveLocal
    REAL, DIMENSION(1:nTurbines, 1:nElements) :: LSectAveLocal, DSectAveLocal
    REAL, DIMENSION(1:nTurbines, 1:nElements) :: v1SectAveLocal, vTSectAveLocal
    REAL, DIMENSION(1:nTurbines, 1:nElements) :: vRelSectAveLocal, vThetaSectAveLocal
    REAL, DIMENSION(1:nTurbines, 1:nElements) :: vZSectAveLocal
!
    REAL, DIMENSION(1:nSections, 1:nElements) :: phi, theta, alphaIn, F, Fcd, Fcl
    REAL, DIMENSION(1:nSections, 1:nElements) :: cl, cd, cm, Cn, Ct, L, D, fn, ft
    REAL, DIMENSION(1:nSections, 1:nElements) :: Posx, Posy, Posz, fx, fy, fz
    REAL, DIMENSION(1:nSections, 1:nElements) :: thrust, torque
!
    REAL, DIMENSION(1:nTurbines, 1:nNacPoints) :: uAtNacLocal, vAtNacLocal, wAtNacLocal
    REAL, DIMENSION(1:nTurbines, 1:nNacPoints) :: uAtNac, vAtNac, wAtNac
!
    REAL, DIMENSION(1:nTurbines, 1:nTowPoints) :: uAtTowLocal, vAtTowLocal, wAtTowLocal
    REAL, DIMENSION(1:nTurbines, 1:nTowPoints) :: uAtTow, vAtTow, wAtTow
!
    REAL, DIMENSION(1:nTurbines, 1:nSections, 1:nElements) :: uLESLocal, vLESLocal, wLESLocal
    REAL, DIMENSION(1:nTurbines, 1:nSections, 1:nElements) :: uLES, vLES, wLES
    REAL, DIMENSION(1:nTurbines, 1:nSections, 1:nElements) :: v1, vT, vRel, vTheta, vZ
    REAL, DIMENSION(1:nTurbines, 1:nSections, 1:nElements) :: uLocal, vLocal, wLocal
!
    REAL, DIMENSION(1:nTurbines, 1:nSections, 1:nElements, 1:3, 1:3) :: bladeAlignedVectors
!
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme) :: fxDist, fyDist, fzDist
!
    REAL, DIMENSION(its-1:ite+1, kts:kte, jts-1:jte+1) :: w_z ! z-coordinate of w-nodes (m) ! Grid point (node)
    REAL, DIMENSION(its-1:ite+1, kts:kte, jts-1:jte+1) :: c_z ! z-coordinate of c-nodes (m) ! Like z_at_w but at cell centers (at half level)
    REAL, DIMENSION(its-1:ite+1, kts:kte, jts-1:jte+1) :: dz  ! Thickness of each grid cell (m)
!
    LOGICAL, EXTERNAL :: wrf_dm_on_monitor
!=======================================================================
    itf = MIN0(ite, ide-1)
    jtf = MIN0(jte, jde-1)
    ktf = MIN0(kte, kde-1)
    IF( printDiagOpt .EQ. 1 )THEN
      IF ( wrf_dm_on_monitor() ) THEN
        WRITE(outstring, *)''
        CALL wrf_message ( TRIM(outstring) )
      ENDIF
    ENDIF
    DO j = jts-1, jtf+1   
      DO i = its-1, itf+1
        w_z(i, kts, j) = 0.0      ! meter above ground, add hgt to get height ASL
      ENDDO
    ENDDO
    DO k = kts+1, ktf+1           ! compute 1D column of z-values to determine kEnd for the tile
      w_z(its, k, jts) = w_z(its, k-1, jts) + 1/rdzw(its, k-1, jts)
    ENDDO
    DO j = jts, jtf               ! compute the rest of w_z up to kEnd
      DO k = kts+1, ktf+1
        DO i = its, itf
          w_z(i, k, j) = w_z(i, k-1, j) + 1/rdzw(i, k-1, j)
        ENDDO
      ENDDO
    ENDDO
    DO j = jts, jtf               ! compute the i_start and i_end ghost nodes
      DO k = kts+1, ktf+1
        w_z(its-1 ,k, j) = w_z(its, k, j)
        w_z(itf+1, k, j) = w_z(itf, k, j)
      ENDDO
    ENDDO
    DO k = kts+1, ktf+1           ! compute the j_start and j_end ghost nodes
      DO i = its, itf  
        w_z(i, k, jts-1) = w_z(i, k, jts)
        w_z(i, k, jtf+1) = w_z(i, k, jtf)
      ENDDO
    ENDDO
    DO k = kts+1, ktf+1           ! compute get the corner ghost nodes
      w_z(its-1, k, jts-1) = 0.5*(w_z(its-1, k, jts) + w_z(its, k, jts-1) )
      w_z(its-1, k, jtf+1) = 0.5*(w_z(its-1, k, jtf) + w_z(its, k, jtf+1) )
      w_z(itf+1, k, jts-1) = 0.5*(w_z(itf+1, k, jts) + w_z(itf, k, jts-1) ) 
      w_z(itf+1, k, jtf+1) = 0.5*(w_z(itf+1, k, jtf) + w_z(itf, k, jtf+1) ) 
    ENDDO
    DO i = its-1, itf+1
      DO k = kts, ktf
        DO j = jts-1, jtf+1
          c_z(i, k, j) = 0.5*( w_z(i, k, j) + w_z(i, k+1, j) )
          dz (i, k, j) = w_z(i, k+1, j) - w_z(i, k, j)
        ENDDO
      ENDDO
    ENDDO 
    nPlnWTcounter = 1
    DO kt = 1, nTurbines
      wtID = FINDLOC(turbineTypeVect, VALUE=turbineTypeID(kt), DIM=1)
      zMaxTrb(wtID) = hubHeight(wtID) + (rotorDiameter(wtID)/2) + projectionRadius(wtID)
      DO k = kts+1, ktf+1           ! kEnd = 1 index greater to allow for nonuniform heights within the tile
        IF ( ( w_z(its, k-1, jts) .LT. zMaxTrb(wtID) ) .AND. ( w_z(its, k, jts) .GE. zMaxTrb(wtID) ) ) kEnd(wtID) = k+1
      ENDDO 
      ! Calculation of the boundaries:
      xMinTile(kt) = (its-1)*dx
      yMinTile(kt) = (jts-1)*dy
      xMaxTile(kt) = itf*dx
      yMaxTile(kt) = jtf*dy
#ifdef DM_PARALLEL
      IF( newComm(kt) .NE. MPI_COMM_NULL ) THEN
        IF( readFromRSTFile .EQ. 0 ) THEN
          rhoFstTSLocal(kt, :) = 0.0
          uFstTSLocal(kt, :) = 0.0
          vFstTSLocal(kt, :) = 0.0
          wFstTSLocal(kt, :) = 0.0
          windDirFstTSLocal(kt, :) = 0.0
          eYawLocal(kt) = 0.0
        ENDIF
        IF( restartController .AND. readFromRSTFile .EQ. 1 ) THEN
          iFstLocal(kt) = iFst(kt)
          jFstLocal(kt) = jFst(kt)
          rhoFstTSLocal(kt, :) = rhoFstTS(kt, :)
          uFstTSLocal(kt, :) = uFstTS(kt, :)
          vFstTSLocal(kt, :) = vFstTS(kt, :)
          wFstTSLocal(kt, :) = wFstTS(kt, :)
          windDirFstTSLocal(kt, :) = windDirFstTS(kt, :)
          eYawLocal(kt) = eYaw(kt)
          rotSpeedLocal(kt) = rotSpeed(kt)
          trqGenCommandedOldOld(kt) = trqGenCommandedOld(kt)
          genSpeedOldLocal(kt) = genSpeedOld(kt)
          genSpeedFOldOld(kt) = genSpeedFOld(kt)
          pitchLocal(kt) = pitch(kt)
          pitchOld(kt) = pitch(kt)
          psiOld(kt) = psi(kt)
          tsrSavedLocal(kt) = tsrSaved(kt)
          deltaAzimuthOld(kt) = deltaAzimuth(kt)
          azimuthOld(kt) = azimuth(kt)
          initialYawDone(kt) = .TRUE.
        ENDIF
        uAtRotor(kt) = 0.0
        vAtRotor(kt) = 0.0
        wAtRotor(kt) = 0.0
        uAtRotorLocal(kt) = 0.0
        vAtRotorLocal(kt) = 0.0
        wAtRotorLocal(kt) = 0.0
        ! Interpolation of the u-component:
        ii = INT( (posRotorApex(kt, 1)/dx) )
        jj = INT( (posRotorApex(kt, 2) + 0.5*dy)/dy )
        IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
          kk = MINLOC( ABS(c_z(ii, kts:kEnd(wtID), jj) + 0.5*dzConst(wtID) - posRotorApex(kt, 3) ), 1)
          IF( ( c_z(ii, kk, jj) + 0.5*dzConst(wtID) ) .GT. posRotorApex(kt, 3) ) kk = kk - 1
          CALL interpolateVelocity(                                           &
              &  ii, jj, kk, dx, dy, posRotorApex(kt, 1), posRotorApex(kt, 2) &
              &, posRotorApex(kt, 3), u_2(ii:ii+1, kk:kk+1, jj:jj+1)          &
              &, dz(ii, kk, jj), c_z(ii, kk, jj), 'u', uAtRotorLocal(kt)      &
              &)
        ELSE
          uAtRotorLocal(kt) = -999999.999
        ENDIF
        ! Interpolation of the v-component:
        ii = INT( (posRotorApex(kt, 1) + 0.5*dx)/dx )
        jj = INT( (posRotorApex(kt, 2)/dy) )
        IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
          kk = MINLOC( ABS(c_z(ii, kts:kEnd(wtID), jj) + 0.5*dzConst(wtID) - posRotorApex(kt, 3) ), 1)
          IF( ( c_z(ii, kk, jj) + 0.5*dzConst(wtID) ) .GT. posRotorApex(kt, 3) ) kk = kk - 1
          CALL interpolateVelocity(                                           &
              &  ii, jj, kk, dx, dy, posRotorApex(kt, 1), posRotorApex(kt, 2) &
              &, posRotorApex(kt, 3), v_2(ii:ii+1, kk:kk+1, jj:jj+1)          &
              &, dz(ii, kk, jj), c_z(ii, kk, jj), 'v', vAtRotorLocal(kt)      &
              &)
        ELSE
          vAtRotorLocal(kt) = -999999.999
        ENDIF
        ! Interpolation of the w-component:
        ii = INT( (posRotorApex(kt, 1) + 0.5*dx)/dx )
        jj = INT( (posRotorApex(kt, 2) + 0.5*dy)/dy )
        IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
          kk = MINLOC( ABS(w_z(ii, kts:kEnd(wtID), jj) - posRotorApex(kt, 3) ), 1)
          IF( w_z(ii, kk, jj) .GT. posRotorApex(kt, 3) ) kk = kk - 1
          CALL interpolateVelocity(                             &
              &  ii, jj, kk, dx, dy, posRotorApex(kt, 1), posRotorApex(kt, 2) &
              &, posRotorApex(kt, 3), w_2(ii:ii+1, kk:kk+1, jj:jj+1)          &
              &, dz(ii, kk, jj), w_z(ii, kk, jj), 'w', wAtRotorLocal(kt)      &
              &)
        ELSE
          wAtRotorLocal(kt) = -999999.999
        ENDIF
        CALL MPI_ALLREDUCE( uAtRotorLocal(kt), uAtRotor(kt), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
        CALL MPI_ALLREDUCE( vAtRotorLocal(kt), vAtRotor(kt), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
        CALL MPI_ALLREDUCE( wAtRotorLocal(kt), wAtRotor(kt), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
        wsRotorLocal(kt) = SQRT( uAtRotor(kt)**2 + vAtRotor(kt)**2 + wAtRotor(kt)**2 )
        IF( readFromRSTFile .EQ. 0 ) THEN
          trbDirLocal(kt) = ATAN2( uAtRotor(kt), vAtRotor(kt) )/degrad + 180.0
          IF( specDir .GE. 0.0 ) trbDirLocal(kt) = specDir
          CALL compassToStandard( trbDirLocal(kt) )
        ELSE
          trbDirLocal(kt) = trbDirSaved(kt)
        ENDIF
        IF( doYaw(kt) .OR. ( restartController1 .AND. .NOT. yawOpt ) ) THEN
          trbDirLocal(kt) = trbDirSaved(kt)
        ENDIF
        IF( ( its .LE. iRot(kt) .AND. iRot(kt) .LE. itf ) .AND. ( jts .LE. jRot(kt) .AND. jRot(kt) .LE. jtf ) ) THEN
          hubHeightASLLocal(kt) = ht( iRot(kt), jRot(kt) ) + towerShaftIntersect(kt, 3)
        ELSE
          hubHeightASLLocal(kt) = -999999.999
        ENDIF
        CALL MPI_ALLREDUCE( hubHeightASLLocal(kt), hubHeightASL(kt), 1, MPI_REAL , MPI_MAX, newComm(kt), ierr )
        uFst(kt) = 0.0
        vFst(kt) = 0.0
        wFst(kt) = 0.0
        rhoFst(kt) = 0.0
        uFstLocal(kt) = 0.0
        vFstLocal(kt) = 0.0
        wFstLocal(kt) = 0.0
        rhoFstLocal(kt) = 0.0
        inflowLocation_z(kt) = hubHeightASL(kt)
        ! Interpolation of the u-component:
        inflowLocation_x(kt) = posRotorApex(kt, 1) - inflowLocation(wtID)*COS(trbDirLocal(kt)*degrad)
        inflowLocation_y(kt) = posRotorApex(kt, 2) - inflowLocation(wtID)*SIN(trbDirLocal(kt)*degrad) + 0.5*dy
        ii = INT(inflowLocation_x(kt)/dx)
        jj = INT(inflowLocation_y(kt)/dy)
        ! Set horizontal boundary index in x-direction to i_domain_start if ii is less than or equal to ids
        IF( ii .LE. 0 ) THEN
          ii = ids
          inflowLocation_x(kt) = (ids-1)*dx
        ENDIF
        ! Set horizontal boundary index in x-direction to i_domain_end-1 if ii is greater than ide
        IF( ii .GE. ide ) THEN
          ii = ide-1
          inflowLocation_x(kt) = (ide-1)*dx
        ENDIF
        ! Set horizontal boundary index in y-direction to j_domain_start if jj is less than or equal to jds
        IF( jj .LE. 0 ) THEN
          jj = jds
          inflowLocation_y(kt) = (jds-1)*dy
        ENDIF
        ! Set horizontal boundary index in y-direction to j_domain_end-1 if jj is greater than jde
        IF( jj .GE. jde ) THEN
          jj = jde-1
          inflowLocation_y(kt) = (jde-1)*dy
        ENDIF
        IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
          kk = MINLOC( ABS( c_z(ii, kts:kEnd(wtID), jj) + 0.5*dzConst(wtID) + ht(ii, jj) - inflowLocation_z(kt) ), 1 )
          IF( ( c_z(ii, kk, jj) + 0.5*dzConst(wtID) + ht(ii, jj) ) .GT. inflowLocation_z(kt) ) kk = kk - 1
          CALL interpolateVelocity(                                             &
              &  ii, jj, kk, dx, dy, inflowLocation_x(kt), inflowLocation_y(kt) &
              &, inflowLocation_z(kt), u_2(ii:ii+1, kk:kk+1, jj:jj+1)           &
              &, dz(ii, kk, jj), c_z(ii, kk, jj), 'u', uFstLocal(kt)            &
              &)
        ELSE
          uFstLocal(kt) = -999999.999
        ENDIF
        ! Interpolation of the v-component:
        inflowLocation_x(kt) = posRotorApex(kt, 1) - inflowLocation(wtID)*COS(trbDirLocal(kt)*degrad) + 0.5*dx
        inflowLocation_y(kt) = posRotorApex(kt, 2) - inflowLocation(wtID)*SIN(trbDirLocal(kt)*degrad)
        ii = INT(inflowLocation_x(kt)/dx)
        jj = INT(inflowLocation_y(kt)/dy)
        ! Set horizontal boundary index in x-direction to i_domain_start if ii is less than or equal to ids
        IF( ii .LE. 0 ) THEN
          ii = ids
          inflowLocation_x(kt) = (ids-1)*dx
        ENDIF
        ! Set horizontal boundary index in x-direction to i_domain_end-1 if ii is greater than ide
        IF( ii .GE. ide ) THEN
          ii = ide-1
          inflowLocation_x(kt) = (ide-1)*dx
        ENDIF
        ! Set horizontal boundary index in y-direction to j_domain_start if jj is less than or equal to jds
        IF( jj .LE. 0 ) THEN
          jj = jds
          inflowLocation_y(kt) = (jds-1)*dy
        ENDIF
        ! Set horizontal boundary index in y-direction to j_domain_end-1 if jj is greater than jde
        IF( jj .GE. jde ) THEN
          jj = jde-1
          inflowLocation_y(kt) = (jde-1)*dy
        ENDIF
        IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
          kk = MINLOC( ABS( c_z(ii, kts:kEnd(wtID), jj) + 0.5*dzConst(wtID) + ht(ii, jj) - inflowLocation_z(kt) ), 1 )
          IF( ( c_z(ii, kk, jj) + 0.5*dzConst(wtID) + ht(ii, jj) ) .GT. inflowLocation_z(kt) ) kk = kk - 1
          CALL interpolateVelocity(                                             &
              &  ii, jj, kk, dx, dy, inflowLocation_x(kt), inflowLocation_y(kt) &
              &, inflowLocation_z(kt), v_2(ii:ii+1, kk:kk+1, jj:jj+1)           &
              &, dz(ii, kk, jj), c_z(ii, kk, jj), 'v', vFstLocal(kt)            &
              &)
        ELSE
          vFstLocal(kt) = -999999.999
        ENDIF
        ! Interpolation of the w-component:
        inflowLocation_x(kt) = posRotorApex(kt, 1) - inflowLocation(wtID)*COS(trbDirLocal(kt)*degrad) + 0.5*dx
        inflowLocation_y(kt) = posRotorApex(kt, 2) - inflowLocation(wtID)*SIN(trbDirLocal(kt)*degrad) + 0.5*dy
        ii = INT(inflowLocation_x(kt)/dx)
        jj = INT(inflowLocation_y(kt)/dy)
        ! Set horizontal boundary index in x-direction to i_domain_start if ii is less than or equal to ids
        IF( ii .LE. 0 ) THEN
          ii = ids
          inflowLocation_x(kt) = (ids-1)*dx
        ENDIF
        ! Set horizontal boundary index in x-direction to i_domain_end-1 if ii is greater than ide
        IF( ii .GE. ide ) THEN
          ii = ide-1
          inflowLocation_x(kt) = (ide-1)*dx
        ENDIF
        ! Set horizontal boundary index in y-direction to j_domain_start if jj is less than or equal to jds
        IF( jj .LE. 0 ) THEN
          jj = jds
          inflowLocation_y(kt) = (jds-1)*dy
        ENDIF
        ! Set horizontal boundary index in y-direction to j_domain_end-1 if jj is greater than jde
        IF( jj .GE. jde ) THEN
          jj = jde-1
          inflowLocation_y(kt) = (jde-1)*dy
        ENDIF
        IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
          kk = MINLOC( ABS( ( w_z(ii, kts:kEnd(wtID), jj) + ht(ii, jj) ) - inflowLocation_z(kt) ), 1)
          IF( ( w_z(ii, kk, jj) + ht(ii, jj) ) .GT. inflowLocation_z(kt) ) kk = kk - 1
          CALL interpolateVelocity(                                             &
              &  ii, jj, kk, dx, dy, inflowLocation_x(kt), inflowLocation_y(kt) &
              &, inflowLocation_z(kt), w_2(ii:ii+1, kk:kk+1, jj:jj+1)           &
              &, dz(ii, kk, jj), w_z(ii, kk, jj), 'w', wFstLocal(kt)            &
              &)
        ELSE
          wFstLocal(kt) = -999999.999
        ENDIF
        ! Interpolation of density:
        IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
          kk = MINLOC( ABS( c_z(ii, kts:kEnd(wtID), jj) + 0.5*dzConst(wtID) + ht(ii, jj) - inflowLocation_z(kt) ), 1 )
          IF( ( c_z(ii, kk, jj) + 0.5*dzConst(wtID) + ht(ii, jj) ) .GT. inflowLocation_z(kt) ) kk = kk - 1
          CALL interpolateScalar(                                               &
              &  ii, jj, kk, dx, dy, inflowLocation_x(kt), inflowLocation_y(kt) &
              &, inflowLocation_z(kt), rho(ii:ii+1, kk:kk+1, jj:jj+1)           &
              &, dz(ii, kk, jj), c_z(ii, kk, jj), rhoFstLocal(kt)               &
              &)
        ELSE
          rhoFstLocal(kt) = -999999.999
        ENDIF
        CALL MPI_ALLREDUCE( uFstLocal(kt), uFst(kt), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
        CALL MPI_ALLREDUCE( vFstLocal(kt), vFst(kt), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
        CALL MPI_ALLREDUCE( wFstLocal(kt), wFst(kt), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
        CALL MPI_ALLREDUCE( rhoFstLocal(kt), rhoFst(kt), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
        iFstLocal(kt) = ii
        jFstLocal(kt) = jj
        IF( ( its .LE. iFstLocal(kt) .AND. iFstLocal(kt) .LE. itf ) .AND. ( jts .LE. jFstLocal(kt) .AND. jFstLocal(kt) .LE. jtf ) ) THEN
          CALL wrf_get_myproc( procFst(kt) )
        ELSE
          procFst(kt) = -1
        ENDIF
        CALL MPI_ALLREDUCE( procFst(kt), procIDFst(kt), 1, MPI_INT , MPI_MAX, newComm(kt), ierr )
        windDirFst(kt) = ATAN2( uFst(kt), vFst(kt) )/degrad + 180.0
        ! Shift rho,u,v, and w for each turbine back by one unit to
        ! make room for the current values: 
        DO n = nAveSamples-1, 1, -1
          rhoFstTSLocal(kt, n+1) = rhoFstTSLocal(kt, n)
          uFstTSLocal(kt, n+1) = uFstTSLocal(kt, n)
          vFstTSLocal(kt, n+1) = vFstTSLocal(kt, n)
          wFstTSLocal(kt, n+1) = wFstTSLocal(kt, n)
          windDirFstTSLocal(kt, n+1) = windDirFstTSLocal(kt, n)
        ENDDO
        ! Sum all values to check whether this is the
        ! first timestep or not:
        rhoFstSum = 0.0
        uFstSum = 0.0
        vFstSum = 0.0
        wFstSum = 0.0
        windDirFstSum = 0.0
        DO n = 1, nAveSamples
          rhoFstSum = rhoFstSum + rhoFstTSLocal(kt, n)
          uFstSum = uFstSum + uFstTSLocal(kt, n)
          vFstSum = vFstSum + vFstTSLocal(kt, n)
          wFstSum = wFstSum + wFstTSLocal(kt, n)
          windDirFstSum = windDirFstSum + windDirFstTSLocal(kt, n)
        ENDDO
        IF( rhoFstSum .EQ. 0.0 ) rhoFstTSLocal(kt, :) = (/(rhoFst(kt), n=1, nAveSamples)/)
        IF( uFstSum .EQ. 0.0 ) uFstTSLocal(kt, :) = (/(uFst(kt), n=1, nAveSamples)/)
        IF( vFstSum .EQ. 0.0 ) vFstTSLocal(kt, :) = (/(vFst(kt), n=1, nAveSamples)/)
        IF( wFstSum .EQ. 0.0 ) wFstTSLocal(kt, :) = (/(wFst(kt), n=1, nAveSamples)/)
        IF( windDirFstSum .EQ. 0.0 ) windDirFstTSLocal(kt, :) = (/(windDirFst(kt), n=1, nAveSamples)/)
        ! Assign the first elements to current values:
        rhoFstTSLocal(kt, 1) = rhoFst(kt)
        uFstTSLocal(kt, 1) = uFst(kt)
        vFstTSLocal(kt, 1) = vFst(kt)
        wFstTSLocal(kt, 1) = wFst(kt)
        windDirFstTSLocal(kt, 1) = windDirFst(kt)
        ! Running time-average:
        rhoFstSum = 0.0
        uFstSum = 0.0
        vFstSum = 0.0
        wFstSum = 0.0
        windDirFstSum = 0.0
        sinFstSum = 0.0
        cosFstSum = 0.0          
        DO n = 1, nAveSamples
          rhoFstSum = rhoFstSum + rhoFstTSLocal(kt, n)
          uFstSum = uFstSum + uFstTSLocal(kt, n)
          vFstSum = vFstSum + vFstTSLocal(kt, n)
          wFstSum = wFstSum + wFstTSLocal(kt, n)
          windDirFstSum = windDirFstSum + windDirFstTSLocal(kt, n)
          sinFstSum = sinFstSum + SIN( windDirFstTSLocal(kt, n)*degrad )
          cosFstSum = cosFstSum + COS( windDirFstTSLocal(kt, n)*degrad )
        ENDDO
        rhoFstAveLocal(kt) = rhoFstSum/REAL(nAveSamples)
        uFstAveLocal(kt) = uFstSum/REAL(nAveSamples)
        vFstAveLocal(kt) = vFstSum/REAL(nAveSamples)
        wFstAveLocal(kt) = wFstSum/REAL(nAveSamples)
        windDirFstAveLocal(kt) = MOD(ATAN2(sinFstSum, cosFstSum)/degrad + 360.0, 360.0)
        CALL compassToStandard( windDirFstAveLocal(kt) )
        capPhiLocal(kt) = trbDirLocal(kt)
        ! Find local velocity matrix:
        velocityFstXYZ(1) = uFstAveLocal(kt)
        velocityFstXYZ(2) = vFstAveLocal(kt)
        velocityFstXYZ(3) = wFstAveLocal(kt)
        CALL rotGlobalToLocal(-beta(wtID), 0.0, -delta(wtID), capPhiLocal(kt), velocityFstXYZ, rotVelocityFstXYZ)
        ! Unperturbed hub-height upstream velocity:
        v0FstLocal(kt) = rotVelocityFstXYZ(1)
        IF( ( its .LE. iRot(kt) .AND. iRot(kt) .LE. itf ) .AND. ( jts .LE. jRot(kt) .AND. jRot(kt) .LE. jtf ) ) THEN
          IF( doYaw(kt) ) THEN
            eYawLocal(kt) = eYaw(kt)
            trbDirLocal(kt) = trbDirSaved(kt)
            rotorApexLocal(kt, 1) = rotorApexX(kt)
            rotorApexLocal(kt, 2) = rotorApexY(kt)
            rotorApexLocal(kt, 3) = rotorApexZ(kt)
          ENDIF
          deltaCapPhiLocal(kt) = trbDirLocal(kt) - windDirFstAveLocal(kt)
          IF( deltaCapPhiLocal(kt) .LT. -180.0 ) THEN 
            deltaCapPhiLocal(kt) = trbDirLocal(kt) + ( 360.0 - windDirFstAveLocal(kt) )
          ELSEIF (deltaCapPhiLocal(kt) .GE. 180.0 ) THEN 
            deltaCapPhiLocal(kt) = trbDirLocal(kt) - ( 360.0 + windDirFstAveLocal(kt) )
          ENDIF
          ! Turbine yawing:
          IF( ABS( deltaCapPhiLocal(kt) ) .GT. 90.0 ) THEN
            initialYawLocal(kt) = .TRUE.
            eYawLocal(kt) = thresholdYaw(wtID) + thresholdNOTYawing
            IF( ABS( deltaCapPhiLocal(kt) ) .LT. thresholdCapPhi ) eYawLocal(kt) = 0.0
            IF( eYawLocal(kt) .NE. 0.0 ) THEN
              ! Check in which direction to yaw:
              yawdir(kt) = -SIGN( 1.0, deltaCapPhiLocal(kt) )
              ! Start turbine yaw:
              deltaYawLocal(kt) = yawdir(kt)*yawRate(wtID)*dt
              trbDirLocal(kt) = trbDirLocal(kt) + deltaYawLocal(kt)
              ! Keep angle between 0 and 360 degrees:
              IF( trbDirLocal(kt) .GE. 360.0 ) trbDirLocal(kt) = trbDirLocal(kt) - 360.0
              IF( trbDirLocal(kt) .LT.   0.0 ) trbDirLocal(kt) = trbDirLocal(kt) + 360.0
              doYawLocal(kt) = .TRUE.
            ENDIF
            IF( doYawLocal(kt) ) THEN
              ! Update rotor apex position:
              rotorApexLocal(kt, :) = rotatePoint( rotorApexLocal(kt, :), towerShaftIntersect(kt, :), uvTower(kt, :), deltaYawLocal(kt) )
              iRotLocal(kt) = INT( (rotorApexLocal(kt, 1) + 0.5*dx)/dx)
              jRotLocal(kt) = INT( (rotorApexLocal(kt, 2) + 0.5*dy)/dy)
            ENDIF
          ELSE
            initialYawLocal(kt) = .FALSE.
          ENDIF
          ! Turbine not yawing:
          IF( ABS( eYawLocal(kt) ) .LE. thresholdYaw(wtID) ) THEN
            IF( ABS( deltaCapPhiLocal(kt) ) .GE. thresholdNOTYawing ) THEN
              eYawLocal(kt) = eYawLocal(kt) + ( -SIGN( deltaCapPhiLocal(kt)          &
                            , deltaCapPhiLocal(kt) ) / ABS( deltaCapPhiLocal(kt) ) ) &
                            * (deltaCapPhiLocal(kt)**2)*dt
            ENDIF
          ! Turbine yawing:
          ELSE
            IF( yawOpt ) THEN
              IF( .NOT. initialYawLocal(kt) ) THEN
                IF( ABS( deltaCapPhiLocal(kt) ) .LT. thresholdCapPhi ) eYawLocal(kt) = 0.0
                IF( eYawLocal(kt) .NE. 0.0 ) THEN
                  ! Check in which direction to yaw:
                  yawdir(kt) = -SIGN( 1.0, deltaCapPhiLocal(kt) )
                  ! Start turbine yaw:
                  deltaYawLocal(kt) = yawdir(kt)*yawRate(wtID)*dt
                  trbDirLocal(kt) = trbDirLocal(kt) + deltaYawLocal(kt)
                  ! Keep angle between 0 and 360 degrees:
                  IF( trbDirLocal(kt) .GE. 360.0 ) trbDirLocal(kt) = trbDirLocal(kt) - 360.0
                  IF( trbDirLocal(kt) .LT.   0.0 ) trbDirLocal(kt) = trbDirLocal(kt) + 360.0
                  doYawLocal(kt) = .TRUE.
                ENDIF
                IF( doYawLocal(kt) ) THEN
                  ! Update rotor apex position:
                  rotorApexLocal(kt, :) = rotatePoint( rotorApexLocal(kt, :), towerShaftIntersect(kt, :), uvTower(kt, :), deltaYawLocal(kt) )
                  iRotLocal(kt) = INT( (rotorApexLocal(kt, 1) + 0.5*dx)/dx)
                  jRotLocal(kt) = INT( (rotorApexLocal(kt, 2) + 0.5*dy)/dy)
                ENDIF
              ENDIF
            ENDIF
          ENDIF
          capPhiLocal(kt) = trbDirLocal(kt)
          CALL wrf_get_myproc( proc(kt) )
        ELSE
          initialYawLocal(kt) = .FALSE.
          doYawLocal(kt) = .FALSE.
          iRotLocal(kt) = -999999
          jRotLocal(kt) = -999999
          rotorApexLocal(kt, :) = -999999.999
          capPhiLocal(kt) = -999999.999
          deltaCapPhiLocal(kt) = -999999.999
          deltaYawLocal(kt) = -999999.999
          eYawLocal(kt) = -999999.999
          trbDirLocal(kt) = -999999.999
        ENDIF
        CALL MPI_ALLREDUCE( doYawLocal(kt), doYaw(kt), 1, MPI_LOGICAL, MPI_LOR, newComm(kt), ierr )
        CALL MPI_ALLREDUCE( initialYawLocal(kt), initialYaw(kt), 1, MPI_LOGICAL, MPI_LOR, newComm(kt), ierr )
        CALL MPI_ALLREDUCE( deltaCapPhiLocal(kt), deltaCapPhi(kt), 1, MPI_LOGICAL, MPI_LOR, newComm(kt), ierr )
        IF( doYaw(kt) ) THEN
          CALL MPI_ALLREDUCE( iRotLocal(kt), iRot(kt), 1, MPI_INT , MPI_MAX, newComm(kt), ierr )
          CALL MPI_ALLREDUCE( jRotLocal(kt), jRot(kt), 1, MPI_INT , MPI_MAX, newComm(kt), ierr )
          CALL MPI_ALLREDUCE( deltaYawLocal(kt), deltaYaw(kt), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
          CALL MPI_ALLREDUCE( rotorApexLocal(kt, 1), rotorApexX(kt), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
          CALL MPI_ALLREDUCE( rotorApexLocal(kt, 2), rotorApexY(kt), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
          CALL MPI_ALLREDUCE( rotorApexLocal(kt, 3), rotorApexZ(kt), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
          doYawLocal(kt) = .FALSE.
        ENDIF
        CALL MPI_ALLREDUCE( proc(kt), procID(kt), 1, MPI_INT , MPI_MAX, newComm(kt), ierr )
        CALL MPI_ALLREDUCE( trbDirLocal(kt), trbDirSaved(kt), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
        CALL MPI_ALLREDUCE( capPhiLocal(kt), capPhi(kt), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
        CALL MPI_ALLREDUCE( eYawLocal(kt), eYaw(kt), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
        IF( tsoutPlnOpt .AND. ANY(specTrb .EQ. kt) ) THEN
          CALL calc_wind_plane_tsout_vel(                    &
              &  ims, ime, jms, jme, kms, kme                &
              &, its, itf, jts, jtf, kts, ktf                &
              &, nPlnWT, nPlnWTcounter                       &
              &, nVertPoints, nHorizPoints, dx, dy           &
              &, dzConst(wtID)                               &
              &, dz(its:itf, kts:ktf, jts:jtf)               &
              &, doYaw(kt), planePoints                      &
              &, towerShaftIntersect(kt, :)                  &
              &, uvTower(kt, :), deltaYaw(kt)                &
              &, c_z(its:itf, kts:ktf, jts:jtf)              &
              &, w_z(its:itf, kts:ktf, jts:jtf)              &
              &, u_2(ims:ime, kms:kme, jms:jme)              &
              &, v_2(ims:ime, kms:kme, jms:jme)              &
              &, w_2(ims:ime, kms:kme, jms:jme), newComm(kt) &
              &, ierr, uPlaneLocal, vPlaneLocal, wPlaneLocal &
              &)
          uPlaneSavedLocal(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints) = uPlaneLocal(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints)
          vPlaneSavedLocal(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints) = vPlaneLocal(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints)
          wPlaneSavedLocal(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints) = wPlaneLocal(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints)
          planePointsLocal(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints, 1:3) = planePoints(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints, 1:3)
        ENDIF
      ENDIF
      proc(kt) = -1
      procFst(kt) = -1
#endif
      IF( specRho .GT. 0.0 ) rhoFstAveLocal(kt) = specRho
#ifdef DM_PARALLEL
      IF( newComm(kt) .NE. MPI_COMM_NULL ) THEN
        IF( .NOT. initialYaw(kt) ) THEN
          IF( v0FstLocal(kt) .GE. cutIn(wtID) .AND. v0FstLocal(kt) .LE. cutOut(wtID) ) THEN
            ! Loop over blade radial sections:
            DO e = 1, nElements
              psi(kt) = psiOld(kt)
              ! Loop over blade sectors:
              DO n = 1, nSections
                uLES(kt, n, e) = 0.0
                vLES(kt, n, e) = 0.0
                wLES(kt, n, e) = 0.0
                uLESLocal(kt, n, e) = 0.0
                vLESLocal(kt, n, e) = 0.0
                wLESLocal(kt, n, e) = 0.0
                IF( doYaw(kt) ) THEN
                  ! Define the vector along the shaft pointing in the direction of the wind
                  posRotorApex(kt, 1) = rotorApexX(kt)
                  posRotorApex(kt, 2) = rotorApexY(kt)
                  posRotorApex(kt, 3) = rotorApexZ(kt)
                  uvShaft(kt, :) = vectorAdd( posRotorApex(kt, :), -towerShaftIntersect(kt, :) )
                  uvShaft(kt, :) = vectorDivide( uvShaft(kt, :), vectorMag( uvShaft(kt, :) ) )
                  uvShaft(kt, :) = vectorMultiply( uvShaft(kt, :), uvShaftDir(wtID) )
                  ! Update blade element position based on yaw motion:
                  bladePoints(kt, n, e, :) = rotatePoint( bladePoints(kt, n, e, :), towerShaftIntersect(kt, :), uvTower(kt, :), deltaYaw(kt) )
                ENDIF
                ! Rotate turbine blades, blade by blade, point by point, and update blade element position:
                bladePoints(kt, n, e, :) = rotatePoint( bladePoints(kt, n, e, :), posRotorApex(kt, :), uvShaft(kt, :), deltaAzimuthOld(kt) )
                ! Define a vector in z':
                bladeAlignedVectors(kt, n, e, 3, :) = vectorAdd( bladePoints(kt, n, e, :), posRotorApex(kt, :) )
                bladeAlignedVectors(kt, n, e, 3, :) = vectorDivide( bladeAlignedVectors(kt, n, e, 3, :), &
                                                      vectorMag( bladeAlignedVectors(kt, n, e, 3, :) ) )
                ! Define a vector in y'
                bladeAlignedVectors(kt, n, e, 2, :) = crossProduct( bladeAlignedVectors(kt, n, e, 3, :), uvShaft(kt, :) )
                bladeAlignedVectors(kt, n, e, 2, :) = vectorDivide( bladeAlignedVectors(kt, n, e, 2, :), &
                                                      vectorMag( bladeAlignedVectors(kt, n, e, 2, :) ) )
                ! Define a vector in x'
                bladeAlignedVectors(kt, n, e, 1, :) = crossProduct( bladeAlignedVectors(kt, n, e, 2, :), &
                                                      bladeAlignedVectors(kt, n, e, 3, :) )
                bladeAlignedVectors(kt, n, e, 1, :) = vectorDivide( bladeAlignedVectors(kt, n, e, 1, :), &
                                                      vectorMag( bladeAlignedVectors(kt, n, e, 1, :) ) )
                ! Interpolation of the u-component:
                ii = INT( (bladePoints(kt, n, e, 1)/dx ) )
                jj = INT( (bladePoints(kt, n, e, 2) + 0.5*dy)/dy )
                IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
                  kk = MINLOC( ABS(c_z(ii, kts:kEnd(wtID), jj) + 0.5*dzConst(wtID) - bladePoints(kt, n, e, 3) ), 1)
                  IF( ( c_z(ii, kk, jj) + 0.5*dzConst(wtID) ) .GT. bladePoints(kt, n, e, 3) ) kk = kk - 1
                  CALL interpolateVelocity(                                 &
                      &  ii, jj, kk, dx, dy, bladePoints(kt, n, e, 1)       &
                      &, bladePoints(kt, n, e, 2), bladePoints(kt, n, e, 3) &
                      &, u_2(ii:ii+1, kk:kk+1, jj:jj+1), dz(ii, kk, jj)     &
                      &, c_z(ii, kk, jj), 'u', uLESLocal(kt, n, e)          &
                      &)
                ELSE
                  uLESLocal(kt, n, e) = -999999.999
                ENDIF
                ! Interpolation of the v-component:
                ii = INT( (bladePoints(kt, n, e, 1) + 0.5*dx)/dx )
                jj = INT( (bladePoints(kt, n, e, 2)/dy ) )
                IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
                  kk = MINLOC( ABS(c_z(ii, kts:kEnd(wtID), jj) + 0.5*dzConst(wtID) - bladePoints(kt, n, e, 3) ), 1)
                  IF( ( c_z(ii, kk, jj) + 0.5*dzConst(wtID) ) .GT. bladePoints(kt, n, e, 3) ) kk = kk - 1
                  CALL interpolateVelocity(                                 &
                      &  ii, jj, kk, dx, dy, bladePoints(kt, n, e, 1)       &
                      &, bladePoints(kt, n, e, 2), bladePoints(kt, n, e, 3) &
                      &, v_2(ii:ii+1, kk:kk+1, jj:jj+1), dz(ii, kk, jj)     &
                      &, c_z(ii, kk, jj), 'v', vLESLocal(kt, n, e)          &
                      &)
                ELSE
                  vLESLocal(kt, n, e) = -999999.999
                ENDIF
                ! Interpolation of the w-component:
                ii = INT( (bladePoints(kt, n, e, 1) + 0.5*dx)/dx )
                jj = INT( (bladePoints(kt, n, e, 2) + 0.5*dy)/dy )
                IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
                  kk = MINLOC( ABS(w_z(ii, kts:kEnd(wtID), jj) - bladePoints(kt, n, e, 3) ), 1)
                  IF( w_z(ii, kk, jj) .GT. bladePoints(kt, n, e, 3) ) kk = kk - 1
                  CALL interpolateVelocity(                                 &
                      &  ii, jj, kk, dx, dy, bladePoints(kt, n, e, 1)       &
                      &, bladePoints(kt, n, e, 2), bladePoints(kt, n, e, 3) &
                      &, w_2(ii:ii+1, kk:kk+1, jj:jj+1), dz(ii, kk, jj)     &
                      &, w_z(ii, kk, jj), 'w', wLESLocal(kt, n, e)          &
                      &)
                ELSE
                  wLESLocal(kt, n, e) = -999999.999
                ENDIF
                CALL MPI_ALLREDUCE( uLESLocal(kt, n, e), uLES(kt, n, e), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
                CALL MPI_ALLREDUCE( vLESLocal(kt, n, e), vLES(kt, n, e), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
                CALL MPI_ALLREDUCE( wLESLocal(kt, n, e), wLES(kt, n, e), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
                ! Find local velocity matrix:
                velocityXYZ(1) = uLES(kt, n, e)
                velocityXYZ(2) = vLES(kt, n, e)
                velocityXYZ(3) = wLES(kt, n, e)
                CALL rotGlobalToLocal(-beta(wtID), psi(kt), -delta(wtID), capPhi(kt), velocityXYZ, rotVelocityXYZ)
                ! Store axial, tangential, and relative velocities for each turbine:
                ! Velocity normal to rotor - from LES
                v1(kt, n, e) = rotVelocityXYZ(1)
                ! Tangential velocity without rotation - from LES
                vTheta(kt, n, e) = rotVelocityXYZ(2)
                ! Velocity perpendicular to rotor - from LES
                vZ(kt, n, e) = rotVelocityXYZ(3)
                ! Longitudinal velocity at the blade element location:
                uLocal(kt, n, e) = uLES(kt, n, e)
                ! Lateral velocity at the blade element location:
                vLocal(kt, n, e) = vLES(kt, n, e)
                ! Vertical velocity at the blade element location:
                wLocal(kt, n, e) = wLES(kt, n, e)
                psi(kt) = psi(kt) + angle
                ! Keep angle between 0 and 360 degrees:
                IF( psi(kt) .GE. 360.0 ) psi(kt) = psi(kt) - 360.0
                IF( psi(kt) .LT.   0.0 ) psi(kt) = psi(kt) + 360.0
              ENDDO
            ENDDO
            IF( .NOT. ( ( xMaxTile(kt) .LT. xMinTrb(kt) ) .OR. ( xMinTile(kt) .GE. xMaxTrb(kt) ) .OR. &
                        ( yMaxTile(kt) .LT. yMinTrb(kt) ) .OR. ( yMinTile(kt) .GE. yMaxTrb(kt) ) ) ) THEN
              ! Calculate rotational speed and pitch angle from look-up table:
              IF( .NOT. speedCtrlOpt(wtID) .AND. .NOT. pitchCtrlOpt(wtID) ) THEN
                CALL interpolate(operDataVectLen(wtID), v0FstLocal(kt), v1Vect(wtID, :), thetapVect(wtID, :), thetacOut)
                CALL interpolate(operDataVectLen(wtID), v0FstLocal(kt), v1Vect(wtID, :), omegaVect(wtID, :), omegaOut)
!                CALL interpolate(operDataVectLen(wtID), v0FstLocal(kt), v1Vect(wtID, :), aeroPowVect(wtID, :), aeroPowOut)
!                CALL interpolate(operDataVectLen(wtID), v0FstLocal(kt), v1Vect(wtID, :), aeroThrVect(wtID, :), aeroThrOut)
                rotSpeedLocal(kt) = omegaOut*rpmradsec
                tsrSavedLocal(kt) = rotSpeedLocal(kt)*rTip(wtID)/v0FstLocal(kt)
                pitchOld(kt) = thetacOut
                pitchLocal(kt) = thetacOut
              ELSEIF( speedCtrlOpt(wtID) .OR. pitchCtrlOpt(wtID) ) THEN
                ! Use speed controller without pitch controller, and calculate pitch angle from look-up table:
                IF( speedCtrlOpt(wtID) .AND. .NOT. pitchCtrlOpt(wtID) ) THEN
                  IF( readFromRSTFile .EQ. 0 .OR. .NOT. initialYawDone(kt) .OR. reOperate(kt) ) THEN
                    rotSpeedLocal(kt) = rotSpeedInit(wtID)*rpmradsec
                    tsrSavedLocal(kt) = rotSpeedLocal(kt)*rTip(wtID)/v0FstLocal(kt)
                  ENDIF
                  CALL interpolate(operDataVectLen(wtID), v0FstLocal(kt), v1Vect(wtID, :), thetapVect(wtID, :), thetacOut)
                  pitchOld(kt) = thetacOut
                  pitchLocal(kt) = pitchOld(kt)
                ! Use speed and pitch controllers together:
                ELSEIF( speedCtrlOpt(wtID) .AND. pitchCtrlOpt(wtID) ) THEN
                  IF( readFromRSTFile .EQ. 0 .OR. .NOT. initialYawDone(kt) .OR. reOperate(kt) ) THEN
                    rotSpeedLocal(kt) = rotSpeedInit(wtID)*rpmradsec
                    tsrSavedLocal(kt) = rotSpeedLocal(kt)*rTip(wtID)/v0FstLocal(kt)
                    pitchOld(kt) = pitchInit(wtID)
                    pitchLocal(kt) = pitchOld(kt)
                  ENDIF
                ENDIF
              ENDIF
              thrAcc(kt) = 0.0
              trqAcc(kt) = 0.0
              pitch_sw = .FALSE.
  pitch_loop: DO
                IF( pitchCtrlOpt(wtID) ) THEN
                  IF( .NOT. pitchControlFileExists(wtID) ) THEN
                    IF( pitch_sw ) THEN
                      thrAcc(kt) = 0.0
                      trqAcc(kt) = 0.0
                      pitchLocal(kt) = pitchLocal(kt) + 0.250 ! add a small pitch gain
                    ELSE
!                      ! Limit minimum pitch angle as a function of wind speed using operationalData.tbl
!                      CALL interpolate(operDataVectLen(wtID), v0FstLocal(kt), v1Vect(wtID, :), thetapVect(wtID, :), thetacOut)
!                      pitchLocal(kt) = MAX( pitchOld(kt) - pitchRate(wtID)*dt, thetacOut )
                      pitchLocal(kt) = MAX( pitchOld(kt) - pitchRate(wtID)*dt, MINVAL(thetapVect(wtID, :), DIM=1) )
                    ENDIF
                  ENDIF
                ENDIF
                ! Loop over blade radial sections:
                DO e = 1, nElements
                  psi(kt) = psiOld(kt)
                  ! Loop over blade sectors:
                  DO n = 1, nSections
                    ! Calculate tangential velocity including rotational effects and blade precone:
                    vT(kt, n, e) = rotSpeedLocal(kt)*bladeRadius(wtID, e)*COS(preconeAngle(wtID)*degrad) - vTheta(kt, n, e)
                    ! Calculate relative velocity:
                    vRel(kt, n, e) = SQRT( v1(kt, n, e)**2 + vT(kt, n, e)**2 )
                    ! Calculate flow angle:
                    phi(n, e) = ATAN2( v1(kt, n, e), vT(kt, n, e) )/degrad ! [deg]
                    ! Calculate local pitch angle:
                    theta(n, e) = pitchLocal(kt) + thetaS(wtID, e) ! theta is the local pitch angle which is the
                                                                   ! summation of blade pitch angle (control) and
                                                                   ! structural twist angle of the blade section.
                    ! Calculate local angle of attack (AoA):
                    alphaIn(n, e) = phi(n, e) - theta(n, e)  ! alphaIn is the local AoA and phi is the flow
                                                             ! angle which is the angle between the plane of
                                                             ! rotation and the relative velocity, vRel.
                    ! Calculate root/tip losses:
                    CALL calcAeroLosses(                                                    &
                        &  bladeRadius(wtID, e), phi(n, e), tsrSavedLocal(kt), nBlade(wtID) &
                        &, rTip(wtID), rHub(wtID), aerodynLossOpt, F(n, e)                  &
                        &)
                    ! Calculate airfoil coefficients:
                    CALL interpolate(                                                                         &
                        &  foilPropVectLen( wtID, foilType(wtID, e) ), alphaIn(n, e)                          &
                        &, alphaVect( wtID, foilType(wtID, e), 1:foilPropVectLen( wtID, foilType(wtID, e) ) ) &
                        &, clVect( wtID, foilType(wtID, e), 1:foilPropVectLen( wtID, foilType(wtID, e) ) )    &
                        &, clOut                                                                              )
                    CALL interpolate(                                                                         &
                        &  foilPropVectLen( wtID, foilType(wtID, e) ), alphaIn(n, e)                          &
                        &, alphaVect( wtID, foilType(wtID, e), 1:foilPropVectLen( wtID, foilType(wtID, e) ) ) &
                        &, cdVect( wtID, foilType(wtID, e), 1:foilPropVectLen( wtID, foilType(wtID, e) ) )    &
                        &, cdOut                                                                              )
                    CALL interpolate(                                                                         &
                        &  foilPropVectLen( wtID, foilType(wtID, e) ), alphaIn(n, e)                          &
                        &, alphaVect( wtID, foilType(wtID, e), 1:foilPropVectLen( wtID, foilType(wtID, e) ) ) &
                        &, cmVect( wtID, foilType(wtID, e), 1:foilPropVectLen( wtID, foilType(wtID, e) ) )    &
                        &, cmOut                                                                              )
                    ! Applies 3-D corrections for rotating sections from the 2-D data:
                    IF( correction3DOpt .GT. 0 ) THEN
                      ! Find the angle of attack corresponding to zero lift:
                      indexStart = MINLOC( ABS( alphaVect(wtID, foilType(wtID, e), :) - alphaMin ), DIM=1 )
                      indexEnd   = MINLOC( ABS( alphaVect(wtID, foilType(wtID, e), :) - alphaMax ), DIM=1 )
                      CALL interpolate(                                              &
                          &  (indexEnd-indexStart+1), 0.0                            &
                          &, clVect(wtID, foilType(wtID, e), indexStart:indexEnd)    &
                          &, alphaVect(wtID, foilType(wtID, e), indexStart:indexEnd) &
                          &, alphaZero                                               &
                          &)
                      ! Find the drag coefficient corresponding to zero angle of attack:
                      CALL interpolate(                                                                         &
                          &  foilPropVectLen( wtID, foilType(wtID, e) ), 0.0                                    &
                          &, alphaVect( wtID, foilType(wtID, e), 1:foilPropVectLen( wtID, foilType(wtID, e) ) ) &
                          &, cdVect( wtID, foilType(wtID, e), 1:foilPropVectLen( wtID, foilType(wtID, e) ) )    &
                          &, cdZero                                                                             &
                          &)
                      IF( correction3DOpt .EQ. 4 ) THEN
                        ! Find the slope for linear region in cl-AoA curve for the 2D data:
                        clMin = clVect( wtID, foilType(wtID, e), indexStart )
                        clMax = clVect( wtID, foilType(wtID, e), indexEnd )
                        IF( clOut .NE. 0.0 ) THEN
                          slope = (alphaMax - alphaMin)/(clMax - clMin)
                        ELSE
                          slope = 0.0
                        ENDIF
                        ! Inviscid lift coefficient from MSc thesis of E. Branlard (2011)
                        clInv = slope*SIN( (alphaIn(n, e) - alphaZero - 0.01 )*degrad )
                        ! Comes from AirfoilPrep spreadsheet:
                        IF( alphaIn(n, e) .LE. alphaThreshold ) THEN
                          adj = 1.0
                        ELSE
                          adj = ( ( (piconst/2) - (alphaIn(n, e)*degrad) )/( (piconst/2)/(alphaThreshold*degrad) ) )**2
                        ENDIF
                      ELSE
                        slope = 2*piconst
                        ! Inviscid lift coefficient:
                        clInv = slope*(alphaIn(n, e) - alphaZero)*degrad
                        adj = 1.0
                      ENDIF
                      ! Find three-dimensional correction factors for cl and cd:
                      CALL correction3D(                                     &
                          &  bladeRadius(wtID, e)                            &
                          &, c(wtID, e), rotSpeedLocal(kt), v0FstLocal(kt)   &
                          &, tsrSavedLocal(kt), vRel(kt, n, e), vT(kt, n, e) &
                          &, theta(n, e), alphaZero, slope, correction3DOpt  &
                          &, rTip(wtID), Fcl(n, e), Fcd(n, e)                &
                          &)
                      ! Apply three-dimensional correction to cl:
                      deltaCl = clInv - clOut
                      clOut = clOut + Fcl(n, e)*deltaCl*adj
                      ! Apply three-dimensional correction to cd:
                      deltaCd = cdOut - cdZero
                      cdOut = cdOut + Fcd(n, e)*deltaCd
                    ENDIF
                    ! Apply root/tip losses:
                    cl(n, e) = F(n, e)*clOut
                    cd(n, e) = F(n, e)*cdOut
                    cm(n, e) = F(n, e)*cmOut
                    ! Calculate normal and tangential coefficients:
                    Cn(n, e) = cl(n, e)*COS(phi(n, e)*degrad) + cd(n, e)*SIN(phi(n, e)*degrad)
                    Ct(n, e) = cl(n, e)*SIN(phi(n, e)*degrad) - cd(n, e)*COS(phi(n, e)*degrad)
                    ! Calculate lift and drag forces per length:
                    L(n, e) = 0.5*rhoFstAveLocal(kt)*c(wtID, e)*cl(n, e)*(vRel(kt, n, e)**2)
                    D(n, e) = 0.5*rhoFstAveLocal(kt)*c(wtID, e)*cd(n, e)*(vRel(kt, n, e)**2)
                    ! Calculate normal and tangential forces per length:
                    fn(n, e) = L(n, e)*COS(phi(n, e)*degrad) + D(n, e)*SIN(phi(n, e)*degrad)
                    ft(n, e) = L(n, e)*SIN(phi(n, e)*degrad) - D(n, e)*COS(phi(n, e)*degrad)
                    ! Find global force matrix:
                    forceXYZ(1) =  fn(n, e)
                    forceXYZ(2) = -ft(n, e)
                    forceXYZ(3) =  0.0
                    CALL rotLocalToGlobal(capPhi(kt), -delta(wtID), psi(kt), -beta(wtID), forceXYZ, rotForceXYZ)
                    ! Assign forces to x,y,z components:
                    fx(n, e) = rotForceXYZ(1)
                    fy(n, e) = rotForceXYZ(2)
                    fz(n, e) = rotForceXYZ(3)
                    ! Calculate thrust and torque at actuator points:
                    thrust(n, e) = fn(n, e)*dr(wtID)
                    torque(n, e) = ft(n, e)*bladeRadius(wtID, e)*dr(wtID)
                    ! Store thrust and torque for each turbine:
                    thrAcc(kt) = thrAcc(kt) + thrust(n, e)
                    trqAcc(kt) = trqAcc(kt) + torque(n, e)
                    psi(kt) = psi(kt) + angle
                    ! Keep angle between 0 and 360 degrees:
                    IF( psi(kt) .GE. 360.0 ) psi(kt) = psi(kt) - 360.0
                    IF( psi(kt) .LT.   0.0 ) psi(kt) = psi(kt) + 360.0
                  ENDDO
                ENDDO
                ! Store calculated variables for each turbine. To be written into wrfout file:
                fnSavedLocal(kt, 1:nSections, 1:nElements) = fn
                ftSavedLocal(kt, 1:nSections, 1:nElements) = ft
                FSavedLocal(kt, 1:nSections, 1:nElements) = F
                phiSavedLocal(kt, 1:nSections, 1:nElements) = phi
                alphaSavedLocal(kt, 1:nSections, 1:nElements) = alphaIn
                clSavedLocal(kt, 1:nSections, 1:nElements) = cl
                cdSavedLocal(kt, 1:nSections, 1:nElements) = cd
                cmSavedLocal(kt, 1:nSections, 1:nElements) = cm
                LSavedLocal(kt, 1:nSections, 1:nElements) = L
                DSavedLocal(kt, 1:nSections, 1:nElements) = D
                v1SavedLocal(kt, 1:nSections, 1:nElements) = v1(kt, 1:nSections, 1:nElements)
                vTSavedLocal(kt, 1:nSections, 1:nElements) = vT(kt, 1:nSections, 1:nElements)
                vRelSavedLocal(kt, 1:nSections, 1:nElements) = vRel(kt, 1:nSections, 1:nElements)
                vThetaSavedLocal(kt, 1:nSections, 1:nElements) = vTheta(kt, 1:nSections, 1:nElements)
                vZSavedLocal(kt, 1:nSections, 1:nElements) = vZ(kt, 1:nSections, 1:nElements)
                uSavedLocal(kt, 1:nSections, 1:nElements) = uLocal(kt, 1:nSections, 1:nElements)
                vSavedLocal(kt, 1:nSections, 1:nElements) = vLocal(kt, 1:nSections, 1:nElements)
                wSavedLocal(kt, 1:nSections, 1:nElements) = wLocal(kt, 1:nSections, 1:nElements)
                bladePointsLocal(kt, 1:nSections, 1:nElements, 1:3) = bladePoints(kt, 1:nSections, 1:nElements, 1:3)
                ! Compute thrust and torque and their non-dimensional coefficients
                ! and store them for each turbine. To be written into wrfout file:
                thrSavedLocal(kt) = thrAcc(kt)
                trqSavedLocal(kt) = trqAcc(kt)*COS(preconeAngle(wtID)*degrad)
                thrCoefSavedLocal(kt) = thrSavedLocal(kt)/(0.5*rhoFstAveLocal(kt)*(v0FstLocal(kt)**2)*piconst*rTip(wtID)**2)
                trqCoefSavedLocal(kt) = trqSavedLocal(kt)/(0.5*rhoFstAveLocal(kt)*(v0FstLocal(kt)**2)*piconst*rTip(wtID)**3)
                IF( speedCtrlOpt(wtID) ) THEN
                  ratedGenSpeed = ratedRotSpeed(wtID)*GBRatio(wtID)
                  IF( readFromRSTFile .EQ. 0 .OR. .NOT. initialYawDone(kt) .OR. reOperate(kt) ) THEN
                    rotSpeedF(kt) = rotSpeedLocal(kt)
                    ! The generator speed is given by the product of gear GBRatio and rotor speed:
                    genSpeed(kt) = (rotSpeedLocal(kt)/rpmradsec)*GBRatio(wtID)
                    genSpeedFOldOld(kt) = genSpeed(kt)
                    genSpeedF(kt) = genSpeed(kt)
                    IF( genSpeedF(kt) .LT. ratedGenSpeed ) THEN
                      trqGenCommandedOldOld(kt) = trqGenInit(wtID)
                    ELSE
                      trqGenCommandedOldOld(kt) = ratedGenTorque(wtID)
                    ENDIF
                    reOperate(kt) = .FALSE.
                  ELSE
                    ! Calculation of the current generator speed for rotor speed control:
                    ! The acceleration of the rotor speed is calculated from the force balance
                    ! of the accelerating torque and the torque of the rotating rotor and generator:
                    rotSpeedLocal(kt) = rotSpeedLocal(kt) + ( dt/driveTrainIner(wtID) ) &
                                      * ( trqSavedLocal(kt)*GBEfficiency(wtID) &
                                      - GBRatio(wtID)*trqGenCommandedOldOld(kt) )
                    ! Calculate new rotor speed and confine with rated rotor speed:
!                    rotSpeedLocal(kt) = MIN( MAX(rotSpeedLocal(kt), 0.0), ratedRotSpeed(wtID)*rpmradsec )
                    ! The generator speed is given by the product of gear GBRatio and rotor speed:
                    genSpeed(kt) = (rotSpeedLocal(kt)/rpmradsec)*GBRatio(wtID)
                    lowPassCoeff = EXP( -dt*( 2*piconst*speedFilterCornerFrequency(wtID) ) )
                    genSpeedF(kt) = (1.0 - lowPassCoeff)*genSpeed(kt) + lowPassCoeff*genSpeedFOldOld(kt)
                    rotSpeedF(kt) = (genSpeedF(kt)*rpmradsec)/GBRatio(wtID)
                  ENDIF
                  ! Region 1: Generator torque is set to zero to accelerate the rotor:
                  IF( genSpeedF(kt) .LT. cutInGenSpeed(wtID) ) THEN
                    trqGenCommanded(kt) = cutInGenTorque(wtID)
                  ! Region 1.5: Generator torque is increasing linearly with rotor speed:
                  ELSEIF( genSpeedF(kt) .GE. cutInGenSpeed(wtID) .AND. genSpeedF(kt) .LT. region2StartGenSpeed(wtID) ) THEN
                    dGenSpeed(kt) = genSpeedF(kt) - cutInGenSpeed(wtID)
                    region2StartGenTorque = KGen(wtID)*region2StartGenSpeed(wtID)**2
                    trqSlope = ( region2StartGenTorque - cutInGenTorque(wtID) )/( region2StartGenSpeed(wtID) - cutInGenSpeed(wtID) )
                    trqGenCommanded(kt) = cutInGenTorque(wtID) + trqSlope*dGenSpeed(kt)
                  ! Region 2: Generator torque is increased by the square of the generator speed to keep the TSR optimal:
                  ELSEIF( genSpeedF(kt) .GE. region2StartGenSpeed(wtID) .AND. genSpeedF(kt) .LT. region2EndGenSpeed(wtID) ) THEN
                    trqGenCommanded(kt) = KGen(wtID)*genSpeedF(kt)**2
                  ! Region 2.5: Slipage region between 2 and 3:
                  ELSEIF( genSpeedF(kt) .GE. region2EndGenSpeed(wtID) .AND. genSpeedF(kt) .LT. ratedGenSpeed ) THEN
                    dGenSpeed(kt) = genSpeedF(kt) - region2EndGenSpeed(wtID)
                    region2EndGenTorque = KGen(wtID)*region2EndGenSpeed(wtID)**2
                    trqSlope = (ratedGenTorque(wtID) - region2EndGenTorque)/( ratedGenSpeed - region2EndGenSpeed(wtID) )
                    trqGenCommanded(kt) = region2EndGenTorque + trqSlope*dGenSpeed(kt)
                  ! Region 3: Generator torque is antiproportional to the rotor speed to keep the power constant:
                  ELSEIF( genSpeedF(kt) .GE. ratedGenSpeed ) THEN
!                    trqGenCommanded(kt) = ratedGenPow(wtID)/(genSpeedF(kt)*rpmradsec)
!                    trqGenCommanded(kt) = ratedGenTorque(wtID)
                    trqGenCommanded(kt) = (ratedGenTorque(wtID)*ratedGenSpeed)/genSpeedF(kt)
                  ENDIF
                  ! Calculate torque rate and confine with a max:
                  trqGenRateCommanded(kt) = ( trqGenCommanded(kt) - trqGenCommandedOldOld(kt) )/dt
                  trqGenRateCommanded(kt) = MIN( MAX(trqGenRateCommanded(kt), -1.0*rateLimitGenTorque(wtID) ), rateLimitGenTorque(wtID) )
                  ! Calculate new gen torque and confine with max torque:
                  trqGenCommanded(kt) = trqGenCommandedOldOld(kt) + trqGenRateCommanded(kt)*dt
                  trqGenCommanded(kt) = MIN( trqGenCommanded(kt), maxGenTorque(wtID) )
                  ! Compute the generator power:
!                  powGenerator(kt) = trqGenCommanded(kt)*rotSpeedLocal(kt)*GBRatio(wtID)*genEfficiency(wtID)
                  powGenerator(kt) = trqGenCommanded(kt)*(genSpeed(kt)*rpmradsec)*genEfficiency(wtID)
                  ! Store tip speed ratio for next timestep:
                  tsrSavedLocal(kt) = rotSpeedLocal(kt)*rTip(wtID)/v0FstLocal(kt)
                  ! Store filtered generator speed for next timestep:
                  genSpeedFOldOld(kt) = genSpeedF(kt)
                  ! Store filtered generator speed for next timestep:
                  trqGenCommandedOldOld(kt) = trqGenCommanded(kt)
                  ! Store filtered generator speed to write into wrfrst:
                  genSpeedFOldLocal(kt) = genSpeedF(kt)
                  ! Store commanded generator torque to write into wrfrst:
                  trqGenCommandedOldLocal(kt) = trqGenCommanded(kt)
                  ! Store generator speed to write into wrfrst:
                  genSpeedOldLocal(kt) = genSpeed(kt)
                  ! Store generator power to write into wrfout:
                  genPowSavedLocal(kt) = powGenerator(kt)
                ENDIF
                IF( pitchCtrlOpt(wtID) ) THEN
                  IF( pitchControlFileExists(wtID) ) THEN
                    IF( readFromRSTFile .EQ. 0 .OR. .NOT. initialYawDone(kt) .OR. reOperate(kt) ) THEN
                      speedError(kt) = speedErrorInit(wtID)
                      ! Set integrated speed error such that if there is no speed error
                      ! on the first real control step (see the else section below), there
                      ! will not be a change in pitch:
                      GK(kt) = 1/( 1.0 + ( pitchLocal(kt)/pitchK(wtID) ) )
                      intSpeedError(kt) = (pitchLocal(kt)*degrad)/( GK(kt)*pitchControlKI(wtID) )
                      IF( .NOT. speedCtrlOpt(wtID) ) THEN
                        reOperate(kt) = .FALSE.
                      ENDIF
                    ELSEIF( readFromRSTFile .EQ. 1 .OR. initialYawDone(kt) ) THEN
                      ! Set the gain scheduling variable:
                      GK(kt) = 1/( 1.0 + ( pitchLocal(kt)/pitchK(wtID) ) )
                      ! Store the old value of speed error:
                      speedErrorLast(kt) = speedError(kt)
                      ! Compute the low speed shaft speed error:
                      speedError(kt) = rotSpeedF(kt) - (ratedRotSpeed(wtID)*rpmradsec)
                      ! Numerically integrate the speed error over time:
                      intSpeedError(kt) = intSpeedError(kt) + speedError(kt)*dt
                      ! Numerically take the deriviative of speed error w.r.t time:
                      derivSpeedError(kt) = ( speedError(kt) - speedErrorLast(kt) )/dt
                      ! Saturate the integrated speed error based on pitch saturation:
                      intSpeedError(kt) = MIN( MAX( intSpeedError(kt), (pitchMin(wtID)*degrad)/(GK(kt)*pitchControlKI(wtID) ) ) &
                                        , (pitchMax(wtID)*degrad)/(GK(kt)*pitchControlKI(wtID) ) )
                      ! Compute the pitch components from the proportional, integral, and derivative parts and sum them:
                      pitchP(kt) = GK(kt)*pitchControlKP(wtID)*speedError(kt)
                      pitchI(kt) = GK(kt)*pitchControlKI(wtID)*intSpeedError(kt)
                      pitchD(kt) = GK(kt)*pitchControlKD(wtID)*derivSpeedError(kt)
                      pitchCommanded(kt) = pitchP(kt) + pitchI(kt) + pitchD(kt)
                      ! Saturate the pitch based on the pitch limits of the pitch actuator:
                      pitchCommanded(kt) = MIN( MAX( pitchCommanded(kt), (pitchMin(wtID)*degrad) ), (pitchMax(wtID)*degrad) )
                      ! Limit the change in blade pitch:
                      pitchRateCommanded(kt) = ( pitchCommanded(kt) - (pitchLocal(kt)*degrad) )/dt
                      pitchRateCommanded(kt) = MIN( MAX(pitchRateCommanded(kt), -1.0*(rateLimitBladePitch(wtID)*degrad) ) &
                                             , (rateLimitBladePitch(wtID)*degrad) )
                      pitchCommanded(kt) = (pitchLocal(kt)*degrad) + pitchRateCommanded(kt)*dt
                      pitchLocal(kt) = pitchCommanded(kt)/degrad
                    ENDIF
                    ! Store the old value of speed error to write into wrfrst:
                    speedErrorLocal(kt) = speedErrorLast(kt)
                    ! Store the old value of integral speed error to write into wrfrst:
                    intSpeedErrorLocal(kt) = intSpeedError(kt)
                  ELSE
                    IF( speedCtrlOpt(wtID) ) THEN
                      ! If the current generator speed is above rated, the pitch is not saturated and the
                      ! change from the last time step is within the maximum pitch rate, then the pitch calculation
                      ! is repeated with a pitch gain:
                      IF( ( genSpeed(kt) .GE. ratedGenSpeed ) .AND. ( pitchLocal(kt) .LT. MAXVAL(thetapVect(wtID, :), DIM=1) + 5.0 ) .AND. &
                        ( pitchLocal(kt) .LT. (pitchOld(kt) + pitchRate(wtID)*dt) ) ) THEN
                        pitch_sw = .TRUE.
                        ! Go back to beginning of pit_loop:
                        CYCLE pitch_loop
                      ENDIF
                      pitchCommanded(kt) = pitchLocal(kt)*degrad
                      pitchOld(kt) = pitchLocal(kt)
                      pitch_sw = .FALSE.
                    ENDIF
                  ENDIF
                ENDIF
                EXIT pitch_loop
              ENDDO pitch_loop
              deltaAzimuthOld(kt) = rotSpeedLocal(kt)/degrad*dt
              psiOld(kt) = psiOld(kt) + deltaAzimuthOld(kt)
              azimuthOld(kt) = azimuthOld(kt) - deltaAzimuthOld(kt)
              ! Keep angle between 0 and 360 degrees:
              IF( azimuthOld(kt) .GE. 360.0 ) azimuthOld(kt) = azimuthOld(kt) - 360.0
              IF( azimuthOld(kt) .LT.   0.0 ) azimuthOld(kt) = azimuthOld(kt) + 360.0
              IF( psiOld(kt) .GE. 360.0 ) psiOld(kt) = psiOld(kt) - 360.0
              IF( psiOld(kt) .LT.   0.0 ) psiOld(kt) = psiOld(kt) + 360.0
              ! To be written into wrfout file:
              pitchSavedLocal(kt) = pitchLocal(kt)
              ! Compute power and its non-dimensional coefficient and store them for each turbine:
              ! To be written into wrfout file:
              powSavedLocal(kt) = trqSavedLocal(kt)*rotSpeedLocal(kt)
              mechPowSavedLocal(kt) = powSavedLocal(kt)*mechEfficiency(wtID)
              powCoefSavedLocal(kt) = powSavedLocal(kt)/(0.5*rhoFstAveLocal(kt)*(v0FstLocal(kt)**3)*piconst*rTip(wtID)**2)
              ! To be written into wrfout file:
              psiOldLocal(kt) = psiOld(kt)
              azimuthLocal(kt) = azimuthOld(kt)
              deltaAzimuthLocal(kt) = deltaAzimuthOld(kt)
              posXYZ(1) = 0.0
              posXYZ(2) = (rotorDiameter(wtID)/2)
              posXYZ(3) = 0.0
              CALL rotLocalToGlobal(capPhi(kt), 0.0, 0.0, 0.0, posXYZ, rotPosXYZ)
              ! Calculation of the boundaries:
              iSmear(kt) = CEILING( ( ABS( rotPosXYZ(1) ) + projectionRadius(wtID) )/dx )
              jSmear(kt) = CEILING( ( ABS( rotPosXYZ(2) ) + projectionRadius(wtID) )/dy )
              ils = MAX( its, iRot(kt) - iSmear(kt) )
              ile = MIN( ite, iRot(kt) + iSmear(kt) )
              jls = MAX( jts, jRot(kt) - jSmear(kt) )
              jle = MIN( jte, jRot(kt) + jSmear(kt) )
              kls = kTrbBot(wtID) - smoothing
              ! Set lower boundary index to 1 if kls is less than MIN(e_vert)
              IF( kls .LE. kts ) kls = kts
              ! Set upper boundary index to maximum if kle is greater than MAX(e_vert)
              kle = kTrbTop(wtID) + smoothing
              IF( kle .GE. ktf ) kle = ktf
              ! Calculate fx, fy and fz forces which are projected onto WRF's [x, y, z] coordinates:
              ! Loop over j turbine blade levels:
              DO jlev = jls, jle
                ! Loop over k turbine blade levels:
                DO klev = kls, kle
                  ! Loop over i turbine blade levels:
                  DO ilev = ils, ile
                    fxDist(ilev, klev, jlev) = 0.0
                    fyDist(ilev, klev, jlev) = 0.0
                    fzDist(ilev, klev, jlev) = 0.0
                    ! Loop over blade sectors:
                    DO n = 1, nSections
                      ! Loop over blade radial sections:
                      DO e = 1, nElements
                        projectionDistanceU = SQRT( ( (ilev      )*dx - bladePoints(kt, n, e, 1) )**2 &
                                            +       ( (jlev - 0.5)*dy - bladePoints(kt, n, e, 2) )**2 &
                                            +       ( (klev - 0.5)*dzConst(wtID)                      &
                                            -          bladePoints(kt, n, e, 3) )**2                  )
                        projectionDistanceV = SQRT( ( (ilev - 0.5)*dx - bladePoints(kt, n, e, 1) )**2 &
                                            +       ( (jlev      )*dy - bladePoints(kt, n, e, 2) )**2 &
                                            +       ( (klev - 0.5)*dzConst(wtID)                      &
                                            -          bladePoints(kt, n, e, 3) )**2                  )
                        projectionDistanceW = SQRT( ( (ilev - 0.5)*dx - bladePoints(kt, n, e, 1) )**2 &
                                            +       ( (jlev - 0.5)*dy - bladePoints(kt, n, e, 2) )**2 &
                                            +       ( (klev      )*dzConst(wtID)                      &
                                            -          bladePoints(kt, n, e, 3) )**2                  )
                        IF(bladeSmearOpt .EQ. 0) THEN ! 3D Gaussian smearing
                          ! x-direction:
                          IF( projectionDistanceU .LE. projectionRadius(wtID) ) THEN
                            etaEpsilonU3D = EXP( -( projectionDistanceU/bladeEpsilon(wtID, 1) )**2 )/( (bladeEpsilon(wtID, 1)**3)*(piconst**1.5) )
                            coeffU3D = etaEpsilonU3D*dr(wtID)*solidity(wtID)
                            fxDist(ilev, klev, jlev) = fxDist(ilev, klev, jlev) + fx(n, e)*coeffU3D
                          ENDIF
                          ! y-direction:
                          IF( projectionDistanceV .LE. projectionRadius(wtID) ) THEN
                            etaEpsilonV3D = EXP( -( projectionDistanceV/bladeEpsilon(wtID, 1) )**2 )/( (bladeEpsilon(wtID, 1)**3)*(piconst**1.5) )
                            coeffV3D = etaEpsilonV3D*dr(wtID)*solidity(wtID)
                            fyDist(ilev, klev, jlev) = fyDist(ilev, klev, jlev) + fy(n, e)*coeffV3D
                          ENDIF
                          ! z-direction:
                          IF( projectionDistanceW .LE. projectionRadius(wtID) ) THEN
                            etaEpsilonW3D = EXP( -( projectionDistanceW/bladeEpsilon(wtID, 1) )**2 )/( (bladeEpsilon(wtID, 1)**3)*(piconst**1.5) )
                            coeffW3D = etaEpsilonW3D*dr(wtID)*solidity(wtID)
                            fzDist(ilev, klev, jlev) = fzDist(ilev, klev, jlev) + fz(n, e)*coeffW3D
                          ENDIF
                        ELSEIF(bladeSmearOpt .EQ. 1) THEN ! 3D polynomial smearing
                          ! x-direction:
                          IF( (projectionDistanceU**2) .LE. epsilonMin2(wtID) ) THEN
                            etaEpsilonU3D = ( ( polA(wtID)*(projectionDistanceU**2) - polB(wtID) )*(projectionDistanceU**2) + 1.0 )*epsilonFactor(wtID)
                            coeffU3D = etaEpsilonU3D*dr(wtID)*solidity(wtID)
                            fxDist(ilev, klev, jlev) = fxDist(ilev, klev, jlev) + fx(n, e)*coeffU3D
                          ENDIF
                          ! y-direction:
                          IF( (projectionDistanceV**2) .LE. epsilonMin2(wtID) ) THEN
                            etaEpsilonV3D = ( ( polA(wtID)*(projectionDistanceV**2) - polB(wtID) )*(projectionDistanceV**2) + 1.0 )*epsilonFactor(wtID)
                            coeffV3D = etaEpsilonV3D*dr(wtID)*solidity(wtID)
                            fyDist(ilev, klev, jlev) = fyDist(ilev, klev, jlev) + fy(n, e)*coeffV3D
                          ENDIF
                          ! z-direction:
                          IF( (projectionDistanceW**2) .LE. epsilonMin2(wtID) ) THEN
                            etaEpsilonW3D = ( ( polA(wtID)*(projectionDistanceW**2) - polB(wtID) )*(projectionDistanceW**2) + 1.0 )*epsilonFactor(wtID)
                            coeffW3D = etaEpsilonW3D*dr(wtID)*solidity(wtID)
                            fzDist(ilev, klev, jlev) = fzDist(ilev, klev, jlev) + fz(n, e)*coeffW3D
                          ENDIF
                        ELSEIF(bladeSmearOpt .EQ. 2) THEN
                          vectZero(1) = 0.0
                          vectZero(2) = 0.0
                          vectZero(3) = 0.0
                          dir1 = bladeAlignedVectors(kt, n, e, 2, :)
                          dir2 = bladeAlignedVectors(kt, n, e, 1, :)
                          dir3 = bladeAlignedVectors(kt, n, e, 3, :)
                          dir1 = rotatePoint(dir1, vectZero, dir3, -( pitchLocal(kt) + thetaS(wtID, e) ) )
                          dir2 = rotatePoint(dir2, vectZero, dir3, -( pitchLocal(kt) + thetaS(wtID, e) ) )
                          ! x-direction:
                          IF( projectionDistanceU .LE. projectionRadius(wtID) ) THEN
                            disVector(1) = (ilev      )*dx - bladePoints(kt, n, e, 1)
                            disVector(2) = (jlev - 0.5)*dy - bladePoints(kt, n, e, 2)
                            disVector(3) = (klev - 0.5)*dzConst(wtID) - bladePoints(kt, n, e, 3)
                            d1 = DOT_PRODUCT(disVector, dir1)
                            d2 = DOT_PRODUCT(disVector, dir2)
                            d3 = DOT_PRODUCT(disVector, dir3)
                            etaEpsilonU3D = EXP( -( d1/bladeEpsilon(wtID, 1) )**2  &
                                          - ( d2/bladeEpsilon(wtID, 2) )**2        &
                                          - ( d3/bladeEpsilon(wtID, 3) )**2 )      &
                                          / ( ( bladeEpsilon(wtID, 1)              &
                                          * bladeEpsilon(wtID, 2)                  &
                                          * bladeEpsilon(wtID, 3) )*(piconst**1.5) )
                            coeffU3D = etaEpsilonU3D*dr(wtID)*solidity(wtID)
                            fxDist(ilev, klev, jlev) = fxDist(ilev, klev, jlev) + fx(n, e)*coeffU3D
                          ENDIF
                          ! y-direction:
                          IF( projectionDistanceV .LE. projectionRadius(wtID) ) THEN
                            disVector(1) = (ilev - 0.5)*dx - bladePoints(kt, n, e, 1)
                            disVector(2) = (jlev      )*dy - bladePoints(kt, n, e, 2)
                            disVector(3) = (klev - 0.5)*dzConst(wtID) - bladePoints(kt, n, e, 3)
                            d1 = DOT_PRODUCT(disVector, dir1)
                            d2 = DOT_PRODUCT(disVector, dir2)
                            d3 = DOT_PRODUCT(disVector, dir3)
                            etaEpsilonV3D = EXP( -( d1/bladeEpsilon(wtID, 1) )**2  &
                                          - ( d2/bladeEpsilon(wtID, 2) )**2        &
                                          - ( d3/bladeEpsilon(wtID, 3) )**2 )      &
                                          / ( ( bladeEpsilon(wtID, 1)              &
                                          * bladeEpsilon(wtID, 2)                  &
                                          * bladeEpsilon(wtID, 3) )*(piconst**1.5) )
                            coeffV3D = etaEpsilonV3D*dr(wtID)*solidity(wtID)
                            fyDist(ilev, klev, jlev) = fyDist(ilev, klev, jlev) + fy(n, e)*coeffV3D
                          ENDIF
                          ! z-direction:
                          IF( projectionDistanceW .LE. projectionRadius(wtID) ) THEN
                            disVector(1) = (ilev - 0.5)*dx - bladePoints(kt, n, e, 1)
                            disVector(2) = (jlev - 0.5)*dy - bladePoints(kt, n, e, 2)
                            disVector(3) = (klev      )*dzConst(wtID) - bladePoints(kt, n, e, 3)
                            d1 = DOT_PRODUCT(disVector, dir1)
                            d2 = DOT_PRODUCT(disVector, dir2)
                            d3 = DOT_PRODUCT(disVector, dir3)
                            etaEpsilonW3D = EXP( -( d1/bladeEpsilon(wtID, 1) )**2  &
                                          - ( d2/bladeEpsilon(wtID, 2) )**2        &
                                          - ( d3/bladeEpsilon(wtID, 3) )**2 )      &
                                          / ( ( bladeEpsilon(wtID, 1)              &
                                          * bladeEpsilon(wtID, 2)                  &
                                          * bladeEpsilon(wtID, 3) )*(piconst**1.5) )
                            coeffW3D = etaEpsilonW3D*dr(wtID)*solidity(wtID)
                            fzDist(ilev, klev, jlev) = fzDist(ilev, klev, jlev) + fz(n, e)*coeffW3D
                          ENDIF
                        ELSEIF(bladeSmearOpt .EQ. 3) THEN
                          vectZero(1) = 0.0
                          vectZero(2) = 0.0
                          vectZero(3) = 0.0
                          dir1 = bladeAlignedVectors(kt, n, e, 2, :)
                          dir2 = bladeAlignedVectors(kt, n, e, 1, :)
                          dir3 = bladeAlignedVectors(kt, n, e, 3, :)
                          dir1 = rotatePoint(dir1, vectZero, dir3, -( pitchLocal(kt) + thetaS(wtID, e) ) )
                          dir2 = rotatePoint(dir2, vectZero, dir3, -( pitchLocal(kt) + thetaS(wtID, e) ) )
                          bladeEpsilon(wtID, 1) = bladeCepsilon(wtID, 1)*c(wtID, e)
                          bladeEpsilon(wtID, 2) = bladeCepsilon(wtID, 2)*c(wtID, e)*t(wtID, e)
                          bladeEpsilon(wtID, 3) = bladeCepsilon(wtID, 3)*dr(wtID)
                          ! x-direction:
                          IF( projectionDistanceU .LE. ( MAXVAL( c(wtID, :) )/MAXVAL( bladeEpsilon( wtID, : ) ) )*projectionRadius(wtID) ) THEN
                            disVector(1) = (ilev      )*dx - bladePoints(kt, n, e, 1)
                            disVector(2) = (jlev - 0.5)*dy - bladePoints(kt, n, e, 2)
                            disVector(3) = (klev - 0.5)*dzConst(wtID) - bladePoints(kt, n, e, 3)
                            d1 = DOT_PRODUCT(disVector, dir1)
                            d2 = DOT_PRODUCT(disVector, dir2)
                            d3 = DOT_PRODUCT(disVector, dir3)
                            etaEpsilonU3D = EXP( -( d1/bladeEpsilon(wtID, 1) )**2  &
                                          - ( d2/bladeEpsilon(wtID, 2) )**2        &
                                          - ( d3/bladeEpsilon(wtID, 3) )**2 )      &
                                          / ( ( bladeEpsilon(wtID, 1)              &
                                          * bladeEpsilon(wtID, 2)                  &
                                          * bladeEpsilon(wtID, 3) )*(piconst**1.5) )
                            coeffU3D  = etaEpsilonU3D*dr(wtID)*solidity(wtID)
                            fxDist(ilev, klev, jlev) = fxDist(ilev, klev, jlev) + fx(n, e)*coeffU3D
                          ENDIF
                          ! y-direction:
                          IF( projectionDistanceV .LE. ( MAXVAL( c(wtID, :) )/MAXVAL( bladeEpsilon( wtID, : ) ) )*projectionRadius(wtID) ) THEN
                            disVector(1) = (ilev - 0.5)*dx - bladePoints(kt, n, e, 1)
                            disVector(2) = (jlev      )*dy - bladePoints(kt, n, e, 2)
                            disVector(3) = (klev - 0.5)*dzConst(wtID) - bladePoints(kt, n, e, 3)
                            d1 = DOT_PRODUCT(disVector, dir1)
                            d2 = DOT_PRODUCT(disVector, dir2)
                            d3 = DOT_PRODUCT(disVector, dir3)
                            etaEpsilonV3D = EXP( -( d1/bladeEpsilon(wtID, 1) )**2  &
                                          - ( d2/bladeEpsilon(wtID, 2) )**2        &
                                          - ( d3/bladeEpsilon(wtID, 3) )**2 )      &
                                          / ( ( bladeEpsilon(wtID, 1)              &
                                          * bladeEpsilon(wtID, 2)                  &
                                          * bladeEpsilon(wtID, 3) )*(piconst**1.5) )
                            coeffV3D = etaEpsilonV3D*dr(wtID)*solidity(wtID)
                            fyDist(ilev, klev, jlev) = fyDist(ilev, klev, jlev) + fy(n, e)*coeffV3D
                          ENDIF
                          ! z-direction:
                          IF( projectionDistanceW .LE. ( MAXVAL( c(wtID, :) )/MAXVAL( bladeEpsilon( wtID, : ) ) )*projectionRadius(wtID) ) THEN
                            disVector(1) = (ilev - 0.5)*dx - bladePoints(kt, n, e, 1)
                            disVector(2) = (jlev - 0.5)*dy - bladePoints(kt, n, e, 2)
                            disVector(3) = (klev      )*dzConst(wtID) - bladePoints(kt, n, e, 3)
                            d1 = DOT_PRODUCT(disVector, dir1)
                            d2 = DOT_PRODUCT(disVector, dir2)
                            d3 = DOT_PRODUCT(disVector, dir3)
                            etaEpsilonW3D = EXP( -( d1/bladeEpsilon(wtID, 1) )**2  &
                                          - ( d2/bladeEpsilon(wtID, 2) )**2        &
                                          - ( d3/bladeEpsilon(wtID, 3) )**2 )      &
                                          / ( ( bladeEpsilon(wtID, 1)              &
                                          * bladeEpsilon(wtID, 2)                  &
                                          * bladeEpsilon(wtID, 3) )*(piconst**1.5) )
                            coeffW3D = etaEpsilonW3D*dr(wtID)*solidity(wtID)
                            fzDist(ilev, klev, jlev) = fzDist(ilev, klev, jlev) + fz(n, e)*coeffW3D
                          ENDIF
                        ENDIF
                      ENDDO
                    ENDDO
                    ! Apply the x-component of the force to the u-tendency of the flow:
                    ru_tendf(ilev, klev, jlev) = ru_tendf(ilev, klev, jlev) - muu(ilev, jlev)*fxDist(ilev, klev, jlev)/rhoFstAveLocal(kt)
                    ! Apply the y-component of the force to the v-tendency of the flow:
                    rv_tendf(ilev, klev, jlev) = rv_tendf(ilev, klev, jlev) - muv(ilev, jlev)*fyDist(ilev, klev, jlev)/rhoFstAveLocal(kt)
                    ! Apply the z-component of the force to the w-tendency of the flow:
                    rw_tendf(ilev, klev, jlev) = rw_tendf(ilev, klev, jlev) - mut(ilev, jlev)*fzDist(ilev, klev, jlev)/rhoFstAveLocal(kt)
                  ENDDO
                ENDDO
              ENDDO
!=======================================================================
              IF( printDiagOpt .EQ. 2 ) THEN
                DO e = 1, nElements
                  fnSum = 0.0
                  ftSum = 0.0
                  FSum = 0.0
                  phiSum = 0.0
                  alphaSum = 0.0
                  clSum = 0.0
                  cdSum = 0.0
                  cmSum = 0.0
                  LSum = 0.0
                  DSum = 0.0
                  v1Sum = 0.0
                  vtSum = 0.0
                  vRelSum = 0.0
                  vThetaSum = 0.0
                  vZSum = 0.0
                  DO n = 1, nSections
                    fnSum = fnSum + fnSavedLocal(kt, n, e)
                    ftSum = ftSum + ftSavedLocal(kt, n, e)
                    FSum = FSum + FSavedLocal(kt, n, e)
                    phiSum = phiSum + phiSavedLocal(kt, n, e)
                    alphaSum = alphaSum + alphaSavedLocal(kt, n, e)
                    clSum = clSum + clSavedLocal(kt, n, e)
                    cdSum = cdSum + cdSavedLocal(kt, n, e)
                    cmSum = cmSum + cmSavedLocal(kt, n, e)
                    LSum = LSum + LSavedLocal(kt, n, e)
                    DSum = DSum + DSavedLocal(kt, n, e)
                    v1Sum = v1Sum + v1SavedLocal(kt, n, e)
                    vtSum = vtSum + vTSavedLocal(kt, n, e)
                    vRelSum = vRelSum + vRelSavedLocal(kt, n, e)
                    vThetaSum = vThetaSum + vThetaSavedLocal(kt, n, e)
                    vZSum = vZSum + vZSavedLocal(kt, n, e)
                  ENDDO
                  fnSectAveLocal(kt, e) = fnSum/nSections
                  ftSectAveLocal(kt, e) = ftSum/nSections
                  FSectAveLocal(kt, e) = FSum/nSections
                  phiSectAveLocal(kt, e) = phiSum/nSections
                  alphaSectAveLocal(kt, e) = alphaSum/nSections
                  clSectAveLocal(kt, e) = clSum/nSections
                  cdSectAveLocal(kt, e) = cdSum/nSections
                  cmSectAveLocal(kt, e) = cmSum/nSections
                  LSectAveLocal(kt, e) = LSum/nSections
                  DSectAveLocal(kt, e) = DSum/nSections
                  v1SectAveLocal(kt, e) = v1Sum/nSections
                  vTSectAveLocal(kt, e) = vtSum/nSections
                  vRelSectAveLocal(kt, e) = vRelSum/nSections
                  vThetaSectAveLocal(kt, e) = vThetaSum/nSections
                  vZSectAveLocal(kt, e) = vZSum/nSections
                ENDDO
              ENDIF
!=======================================================================
              IF( printDiagOpt .EQ. 2 ) THEN
                WRITE(outstring, *)''
                CALL wrf_message ( TRIM(outstring) )
                WRITE(outstring, "(A, I10, A, A, I10, A)")                    &
                                 ' Wind turbine #   = ', kt, ' [ - ],',       &
                                 ' ITIMESTEP      = ', itimestep, '  [ - ]'
                CALL wrf_message ( TRIM(outstring) )
                WRITE(outstring, "(A, I10, A, A, I10, A, A, I7, A)")          &
                                 ' i                = ', iRot(kt), ' [ - ],', &
                                 ' j              = ', jRot(kt), '  [ - ],',  &
                                 ' procID    = ', procID(kt), '  [ - ]'
                CALL wrf_message ( TRIM(outstring) )
                WRITE(outstring, "(A, I10, A, A, I10, A, A, I7, A)")               &
                                 ' iFst             = ', iFstLocal(kt), ' [ - ],', &
                                 ' jFst           = ', jFstLocal(kt), '  [ - ],',  &
                                 ' procIDFst = ', procIDFst(kt), '  [ - ]'
                CALL wrf_message ( TRIM(outstring) )
                WRITE(outstring, "(A, A, A, A, A, A, A, A, A, A, A, A, A, A, A)")' El. # [ - ]', ' r/R [ - ]', ' v1 [m/s]',       &
                                 ' vTheta [m/s]', ' vT [m/s]', ' Vr [m/s]', ' vZ [m/s]', ' phi [deg]', ' AoA [deg]', ' cl [ - ]', &
                                 ' cd [ - ]', '    L [N/m]', '    D [N/m]', '   ft [N/m]', '   fn [N/m]'
                CALL wrf_message ( TRIM(outstring) )
                WRITE(outstring, "(A, A, A, A, A, A, A, A, A, A, A, A, A, A, A)")' -----------', ' ---------', ' --------', ' ------------', &
                                 ' --------', ' --------', ' --------', ' ---------', ' ---------', ' --------', ' --------', '    -------', &
                                 '    -------', '   --------', '   --------'
                CALL wrf_message ( TRIM(outstring) )
                DO e = 1, nElements
                  WRITE(outstring, "(I12, F10.3, F9.3, F13.3, F9.3, F9.3, F9.3, F10.3, F10.3, F9.3, F9.3, F11.3, F11.3, F11.3, F11.3)") e, &
                                   bladeRadius(wtID, e)/rTip(wtID), v1SectAveLocal(kt, e), vThetaSectAveLocal(kt, e),                      &
                                   vTSectAveLocal(kt, e), vRelSectAveLocal(kt, e), vZSectAveLocal(kt, e), phiSectAveLocal(kt, e),          &
                                   alphaSectAveLocal(kt, e), clSectAveLocal(kt, e), cdSectAveLocal(kt, e), LSectAveLocal(kt, e),           &
                                   DSectAveLocal(kt, e), ftSectAveLocal(kt, e), fnSectAveLocal(kt, e)
                  CALL wrf_message ( TRIM(outstring) )
                ENDDO
                WRITE(outstring, "(A, F10.3, A, A, F9.3, A, A, F8.3, A)")                             &
                                 ' TSR                 = ', tsrSavedLocal(kt), ' [ - ],',             &
                                 ' Total_thrst    = ', (thrSavedLocal(kt)/1000.0), '   [ kN],',       &
                                 ' CT          = ', thrCoefSavedLocal(kt), ' [ - ]'
                CALL wrf_message ( TRIM(outstring) )
                WRITE(outstring, "(A, F10.3, A, A, F9.3, A, A, F8.3, A)")                             &
                                 ' Rot_speed           = ', rotSpeedLocal(kt)/rpmradsec, ' [RPM],',   &
                                 ' Aero_power     = ', (powSavedLocal(kt)/1000.0), '   [ kW],',       &
                                 ' CP          = ', powCoefSavedLocal(kt), ' [ - ]'
                CALL wrf_message ( TRIM(outstring) )
                WRITE(outstring, "(A, F10.3, A, A, F9.3, A, A, F8.3, A)")                             &
                                 ' Wind_dir            = ', windDirFstAveLocal(kt),                   &
                                 ' [deg],', ' delta_capPhi   = ', deltaCapPhi(kt),                    &
                                 '   [deg],', ' V0_LES      = ', v0FstLocal(kt), ' [m/s]'
                CALL wrf_message ( TRIM(outstring) )
                WRITE(outstring, "(A, F10.3, A, A, F9.3, A, A, F8.3, A)")                             &
                                 ' Turbine_dir         = ', trbDirSaved(kt), ' [deg],',               &
                                 ' wsRotor        = ', wsRotorLocal(kt), '   [m/s],',                 &
                                 ' CapPhi      = ', capPhi(kt), ' [deg]'
                CALL wrf_message ( TRIM(outstring) )
                WRITE(outstring, "(A, F10.3, A, A, F9.3, A, A, F8.3, A)")                             &
                                 ' Yaw_error           = ', eYaw(kt), ' [ - ],',                      &
                                 ' Rho            = ', rhoFstAveLocal(kt), ' [kg/m3],',               &
                                 ' Azimuth     = ', azimuthLocal(kt), ' [deg]'
                CALL wrf_message ( TRIM(outstring) )
                IF( speedCtrlOpt(wtID) ) THEN
                  WRITE(outstring, "(A, F10.3, A, A, F9.3, A, A, F8.3, A)")                           &
                                   ' Gen_power           = ', genPowSavedLocal(kt)/1000.0, ' [ kW],', &
                                   ' Genspeed       = ', genSpeedOldLocal(kt), '   [RPM],',           &
                                   ' Pitch       = ', pitchLocal(kt), ' [deg]'
                  CALL wrf_message ( TRIM(outstring) )
                ELSE
                  WRITE(outstring, "(A, F10.3, A)")                                                   &
                                   ' Pitch               = ', pitchLocal(kt), ' [deg]'
                  CALL wrf_message ( TRIM(outstring) )
                ENDIF
              ENDIF
            ELSE
              rotSpeedLocal(kt) = -999999.999
              pitchSavedLocal(kt) = -999999.999
              tsrSavedLocal(kt) = -999999.999
              IF( speedCtrlOpt(wtID) ) THEN
                genSpeedOldLocal(kt) = -999999.999
              ENDIF
            ENDIF
!=======================================================================
            ! Hub Parameterization
            IF( hubParOpt ) THEN
              projectionLimit = rHub(wtID) + projectionRadius(wtID)
              nacelleFrontalArea = piconst*rHub(wtID)**2
              dist = 0.0
              ! Loop over nacelle points:
              DO m = 1, nNacPoints
                uAtNac(kt, m) = 0.0
                vAtNac(kt, m) = 0.0
                wAtNac(kt, m) = 0.0
                uAtNacLocal(kt, m) = 0.0
                vAtNacLocal(kt, m) = 0.0
                wAtNacLocal(kt, m) = 0.0
                IF( doYaw(kt) ) THEN
                  ! Define the vector along the shaft pointing in the direction of the wind
                  posRotorApex(kt, 1) = rotorApexX(kt)
                  posRotorApex(kt, 2) = rotorApexY(kt)
                  posRotorApex(kt, 3) = rotorApexZ(kt)
                  uvShaft(kt, :) = vectorAdd( posRotorApex(kt, :), -towerShaftIntersect(kt, :) )
                  uvShaft(kt, :) = vectorDivide( uvShaft(kt, :), vectorMag( uvShaft (kt, :) ) )
                  uvShaft(kt, :) = vectorMultiply( uvShaft(kt, :), uvShaftDir(wtID) )
                  ! Update blade element position:
                  nacellePoints(kt, m, :) = rotatePoint( nacellePoints(kt, m, :), towerShaftIntersect(kt, :) &
                                          ,              uvTower(kt, :), deltaYaw(kt) )
                ENDIF
                ! Interpolation of the u-component:
                ii = INT( (nacellePoints(kt, m, 1)/dx ) )
                jj = INT( (nacellePoints(kt, m, 2) + 0.5*dy)/dy )
                IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. &
                    ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
                  kk = MINLOC( ABS(c_z(ii, kts:kHubTop(wtID), jj) + 0.5*dzConst(wtID) &
                     - nacellePoints(kt, m, 3) ), 1)
                  IF( ( c_z(ii, kk, jj) + 0.5*dzConst(wtID) ) .GT. nacellePoints(kt, m, 3) ) kk = kk - 1
                  CALL interpolateVelocity(                               &
                      &  ii, jj, kk, dx, dy, nacellePoints(kt, m, 1)      &
                      &, nacellePoints(kt, m, 2), nacellePoints(kt, m, 3) &
                      &, u_2(ii:ii+1, kk:kk+1, jj:jj+1), dz(ii, kk, jj)   &
                      &, c_z(ii, kk, jj), 'u', uAtNacLocal(kt, m)         &
                      &)
                ELSE
                  uAtNacLocal(kt, m) = -999999.999
                ENDIF
                ! Interpolation of the v-component:
                ii = INT( (nacellePoints(kt, m, 1) + 0.5*dx)/dx )
                jj = INT( (nacellePoints(kt, m, 2)/dy ) )
                IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. &
                    ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
                  kk = MINLOC( ABS(c_z(ii, kts:kHubTop(wtID), jj) + 0.5*dzConst(wtID) &
                     - nacellePoints(kt, m, 3) ), 1)
                  IF( ( c_z(ii, kk, jj) + 0.5*dzConst(wtID) ) .GT. nacellePoints(kt, m, 3) ) kk = kk - 1
                  CALL interpolateVelocity(                               &
                      &  ii, jj, kk, dx, dy, nacellePoints(kt, m, 1)      &
                      &, nacellePoints(kt, m, 2), nacellePoints(kt, m, 3) &
                      &, v_2(ii:ii+1, kk:kk+1, jj:jj+1), dz(ii, kk, jj)   &
                      &, c_z(ii, kk, jj), 'v', vAtNacLocal(kt, m)         &
                      &)
                ELSE
                  vAtNacLocal(kt, m) = -999999.999
                ENDIF
                ! Interpolation of the w-component:
                ii = INT( (nacellePoints(kt, m, 1) + 0.5*dx)/dx )
                jj = INT( (nacellePoints(kt, m, 2) + 0.5*dy)/dy )
                IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. &
                    ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
                  kk = MINLOC( ABS(w_z(ii, kts:kHubTop(wtID), jj) - nacellePoints(kt, m, 3) ), 1)
                  IF( w_z(ii, kk, jj) .GT. nacellePoints(kt, m, 3) ) kk = kk - 1
                  IF( kk .EQ. 0 ) kk = kts
                  CALL interpolateVelocity(                               &
                      &  ii, jj, kk, dx, dy, nacellePoints(kt, m, 1)      &
                      &, nacellePoints(kt, m, 2), nacellePoints(kt, m, 3) &
                      &, w_2(ii:ii+1, kk:kk+1, jj:jj+1), dz(ii, kk, jj)   &
                      &, w_z(ii, kk, jj), 'w', wAtNacLocal(kt, m)         &
                      &)
                ELSE
                  wAtNacLocal(kt, m) = -999999.999
                ENDIF
                CALL MPI_ALLREDUCE( uAtNacLocal(kt, m), uAtNac(kt, m), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
                CALL MPI_ALLREDUCE( vAtNacLocal(kt, m), vAtNac(kt, m), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
                CALL MPI_ALLREDUCE( wAtNacLocal(kt, m), wAtNac(kt, m), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
                ! Find local wind direction for nacelle-elements:
                capPhiAtNac(m) = ATAN2( uAtNac(kt, m), vAtNac(kt, m) )/degrad + 180.0
                CALL compassToStandard( capPhiAtNac(m) )
                ! Find local velocity matrix:
                velocityXYZ(1) = uAtNac(kt, m)
                velocityXYZ(2) = vAtNac(kt, m)
                velocityXYZ(3) = wAtNac(kt, m)
                CALL rotGlobalToLocal(0.0, 0.0, -delta(wtID), capPhiAtNac(m), velocityXYZ, rotVelocityXYZ)
                ! Velocity perpendicular to the nacelle:
                vNacX(m) = rotVelocityXYZ(1)
                ! Velocity in y-direction tangential to the nacelle:
                vNacY(m) = rotVelocityXYZ(2)
                ! Velocity in z-direction tangential to the nacelle:
                vNacZ(m) = rotVelocityXYZ(3)
                ! Find the local velocity magnitude:
                vNac(m) = SQRT(vNacX(m)**2 + vNacY(m)**2 + vNacZ(m)**2)
                ! Calculate drag force generated by the nacelle:
                fnNacelle(m) = contribution(wtID)*0.5*rhoFstAveLocal(kt) &
                             * nacDragCoeff(wtID)*(vNac(m)**2)*nacelleFrontalArea
                ! Find global force matrix:
                forceXYZ(1) = fnNacelle(m)
                forceXYZ(2) = 0.0
                forceXYZ(3) = 0.0
                CALL rotLocalToGlobal(capPhiAtNac(m), -delta(wtID), 0.0, 0.0, forceXYZ, rotForceXYZ)
                ! Assign forces to x,y,z components:
                fxNacelle(m) = rotForceXYZ(1)
                fyNacelle(m) = rotForceXYZ(2)
                fzNacelle(m) = rotForceXYZ(3)
              ENDDO
              posXYZ(1)  = 0.0
              posXYZ(2)  = (hubDiameter(wtID)/2)
              posXYZ(3)  = 0.0
              CALL rotLocalToGlobal(capPhi(kt), 0.0, 0.0, 0.0, posXYZ, rotPosXYZ)
              ! Calculation of the boundaries:
              iSmear(kt) = CEILING( ( ABS( rotPosXYZ(1) ) + projectionRadius(wtID) )/dx )
              jSmear(kt) = CEILING( ( ABS( rotPosXYZ(2) ) + projectionRadius(wtID) )/dy )
              ils = MAX( its, iTow(kt) - iSmear(kt) )
              ile = MIN( ite, iTow(kt) + iSmear(kt) )
              jls = MAX( jts, jTow(kt) - jSmear(kt) )
              jle = MIN( jte, jTow(kt) + jSmear(kt) )
              kls = kHubBot(wtID) - smoothing
              kle = kHubTop(wtID) + smoothing
              ! loop over j turbine blade levels
              DO jlev = jls, jle
                ! loop over k turbine blade levels
                DO klev = kls, kle
                  ! loop over i turbine blade levels
                  DO ilev = ils, ile
                    fxDist(ilev, klev, jlev) = 0.0
                    fyDist(ilev, klev, jlev) = 0.0
                    fzDist(ilev, klev, jlev) = 0.0
                    ! Loop over nacelle points:
                    DO m = 1, nNacPoints
                      projectionDistanceU = SQRT( ( (ilev      )*dx - nacellePoints(kt, m, 1) )**2 &
                                          +       ( (jlev - 0.5)*dy - nacellePoints(kt, m, 2) )**2 &
                                          +       ( (klev - 0.5)*dzConst(wtID)                     &
                                          - nacellePoints(kt, m, 3) )**2 )
                      projectionDistanceV = SQRT( ( (ilev - 0.5)*dx - nacellePoints(kt, m, 1) )**2 &
                                          +       ( (jlev      )*dy - nacellePoints(kt, m, 2) )**2 &
                                          +       ( (klev - 0.5)*dzConst(wtID)                     &
                                          - nacellePoints(kt, m, 3) )**2 )
                      projectionDistanceW = SQRT( ( (ilev - 0.5)*dx - nacellePoints(kt, m, 1) )**2 &
                                          +       ( (jlev - 0.5)*dy - nacellePoints(kt, m, 2) )**2 &
                                          +       ( (klev      )*dzConst(wtID)                     &
                                          - nacellePoints(kt, m, 3) )**2 )
                      IF(nacelleSmearOpt .EQ. 0) THEN ! 3D Gaussian smearing
                        ! x-direction:
                        IF(projectionDistanceU .LE. projectionLimit) THEN
                          etaEpsilonU3D = EXP( -( projectionDistanceU/nacelleEpsilon(wtID) )**2 )/( (nacelleEpsilon(wtID)**3)*(piconst**1.5) )
                          coeffU3D = etaEpsilonU3D
                          fxDist(ilev, klev, jlev) = fxDist(ilev, klev, jlev) + fxNacelle(m)*coeffU3D
                        ENDIF
                        ! y-direction:
                        IF(projectionDistanceV .LE. projectionLimit) THEN
                          etaEpsilonV3D = EXP( -( projectionDistanceV/nacelleEpsilon(wtID) )**2 )/( (nacelleEpsilon(wtID)**3)*(piconst**1.5) )
                          coeffV3D = etaEpsilonV3D
                          fyDist(ilev, klev, jlev) = fyDist(ilev, klev, jlev) + fyNacelle(m)*coeffV3D
                        ENDIF
                        ! z-direction:
                        IF(projectionDistanceW .LE. projectionLimit) THEN
                          etaEpsilonW3D = EXP( -( projectionDistanceW/nacelleEpsilon(wtID) )**2 )/( (nacelleEpsilon(wtID)**3)*(piconst**1.5) )
                          coeffW3D = etaEpsilonW3D
                          fzDist(ilev, klev, jlev) = fzDist(ilev, klev, jlev) + fzNacelle(m)*coeffW3D
                        ENDIF
                      ELSEIF(nacelleSmearOpt .EQ. 1) THEN ! 3D polynomial smearing
                        ! x-direction:
                        IF( (projectionDistanceU**2) .LE. epsilonMin2(wtID) ) THEN
                          etaEpsilonU3D = ( ( polA(wtID) *(projectionDistanceU**2) - polB(wtID) )*(projectionDistanceU**2) + 1.0 )*epsilonFactor(wtID)
                          coeffU3D = etaEpsilonU3D
                          fxDist(ilev, klev, jlev) = fxDist(ilev, klev, jlev) + fxNacelle(m)*coeffU3D
                        ENDIF
                        ! y-direction:
                        IF( (projectionDistanceV**2) .LE. epsilonMin2(wtID) ) THEN
                          etaEpsilonV3D = ( ( polA(wtID)*(projectionDistanceV**2) - polB(wtID) )*(projectionDistanceV**2) + 1.0 )*epsilonFactor(wtID)
                          coeffV3D = etaEpsilonV3D
                          fyDist(ilev, klev, jlev) = fyDist(ilev, klev, jlev) + fyNacelle(m)*coeffV3D
                        ENDIF
                        ! z-direction:
                        IF( (projectionDistanceW**2) .LE. epsilonMin2(wtID) ) THEN
                          etaEpsilonW3D = ( ( polA(wtID)*(projectionDistanceW**2) - polB(wtID) )*(projectionDistanceW**2) + 1.0 )*epsilonFactor(wtID)
                          coeffW3D = etaEpsilonW3D
                          fzDist(ilev, klev, jlev) = fzDist(ilev, klev, jlev) + fzNacelle(m)*coeffW3D
                        ENDIF
                      ENDIF
                    ENDDO
                    ! Apply the x,y,z-components of the nacelle forces to the u,v,w-tendencies of the flow:
                    ru_tendf(ilev, klev, jlev) = ru_tendf(ilev, klev, jlev) - muu(ilev, jlev)*fxDist(ilev, klev, jlev)/rhoFstAveLocal(kt)
                    rv_tendf(ilev, klev, jlev) = rv_tendf(ilev, klev, jlev) - muv(ilev, jlev)*fyDist(ilev, klev, jlev)/rhoFstAveLocal(kt)
                    rw_tendf(ilev, klev, jlev) = rw_tendf(ilev, klev, jlev) - mut(ilev, jlev)*fzDist(ilev, klev, jlev)/rhoFstAveLocal(kt)
                  ENDDO
                ENDDO
              ENDDO
              nacellePointsLocal(kt, 1:nNacPoints, 1:3) = nacellePoints(kt, 1:nNacPoints, 1:3)
            ENDIF
!=======================================================================
            ! Tower Parameterization
            IF( towParOpt ) THEN
              projectionLimit = rTow(wtID) + projectionRadius(wtID)
              dist = 0.0
              ! Loop over tower points:
              DO m = 1, nTowPoints
                uAtTow(kt, m) = 0.0
                vAtTow(kt, m) = 0.0
                wAtTow(kt, m) = 0.0
                uAtTowLocal(kt, m) = 0.0
                vAtTowLocal(kt, m) = 0.0
                wAtTowLocal(kt, m) = 0.0
                ! Interpolation of the u-component:
                ii = INT( (towerPoints(kt, m, 1)/dx ) )
                jj = INT( (towerPoints(kt, m, 2) + 0.5*dy)/dy )
                IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
                  kk = MINLOC( ABS(c_z(ii, kts:kTowTop(wtID), jj) + 0.5*dzConst(wtID) - towerPoints(kt, m, 3) ), 1)
                  IF( ( c_z(ii, kk, jj) + 0.5*dzConst(wtID) ) .GT. towerPoints(kt, m, 3) ) kk = kk - 1
                  IF( kk .EQ. 0 ) kk = kts
                  CALL interpolateVelocity(                             &
                      &  ii, jj, kk, dx, dy, towerPoints(kt, m, 1)      &
                      &, towerPoints(kt, m, 2), towerPoints(kt, m, 3)   &
                      &, u_2(ii:ii+1, kk:kk+1, jj:jj+1), dz(ii, kk, jj) &
                      &, c_z(ii, kk, jj), 'u', uAtTowLocal(kt, m)       &
                      &)
                ELSE
                  uAtTowLocal(kt, m) = -999999.999
                ENDIF
                ! Interpolation of the v-component:
                ii = INT( (towerPoints(kt, m, 1) + 0.5*dx)/dx )
                jj = INT( (towerPoints(kt, m, 2)/dy ) )
                IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
                  kk = MINLOC( ABS(c_z(ii, kts:kTowTop(wtID), jj) + 0.5*dzConst(wtID) - towerPoints(kt, m, 3) ), 1)
                  IF( ( c_z(ii, kk, jj) + 0.5*dzConst(wtID) ) .GT. towerPoints(kt, m, 3) ) kk = kk - 1
                  IF( kk .EQ. 0 ) kk = kts
                  CALL interpolateVelocity(                             &
                      &  ii, jj, kk, dx, dy, towerPoints(kt, m, 1)      &
                      &, towerPoints(kt, m, 2), towerPoints(kt, m, 3)   &
                      &, v_2(ii:ii+1, kk:kk+1, jj:jj+1), dz(ii, kk, jj) &
                      &, c_z(ii, kk, jj), 'v', vAtTowLocal(kt, m)       &
                      &)
                ELSE
                  vAtTowLocal(kt, m) = -999999.999
                ENDIF
                ! Interpolation of the w-component:
                ii = INT( (towerPoints(kt, m, 1) + 0.5*dx)/dx )
                jj = INT( (towerPoints(kt, m, 2) + 0.5*dy)/dy )
                IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
                  kk = MINLOC( ABS(w_z(ii, kts:kTowTop(wtID), jj) - towerPoints(kt, m, 3) ), 1)
                  IF( w_z(ii, kk, jj) .GT. towerPoints(kt, m, 3) ) kk = kk - 1
                  IF( kk .EQ. 0 ) kk = kts
                  CALL interpolateVelocity(                             &
                      &  ii, jj, kk, dx, dy, towerPoints(kt, m, 1)      &
                      &, towerPoints(kt, m, 2), towerPoints(kt, m, 3)   &
                      &, w_2(ii:ii+1, kk:kk+1, jj:jj+1), dz(ii, kk, jj) &
                      &, w_z(ii, kk, jj), 'w', wAtTowLocal(kt, m)       &
                      &)
                ELSE
                  wAtTowLocal(kt, m) = -999999.999
                ENDIF
                CALL MPI_ALLREDUCE( uAtTowLocal(kt, m), uAtTow(kt, m), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
                CALL MPI_ALLREDUCE( vAtTowLocal(kt, m), vAtTow(kt, m), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
                CALL MPI_ALLREDUCE( wAtTowLocal(kt, m), wAtTow(kt, m), 1, MPI_REAL, MPI_MAX, newComm(kt), ierr )
                ! Find local wind direction for tower-elements:
                capPhiAtTow(m) = ATAN2( uAtTow(kt, m), vAtTow(kt, m) )/degrad + 180.0
                CALL compassToStandard( capPhiAtTow(m) )
                ! Find local velocity matrix:
                velocityXYZ (1) = uAtTow(kt, m)
                velocityXYZ (2) = vAtTow(kt, m)
                velocityXYZ (3) = wAtTow(kt, m)
                CALL rotGlobalToLocal(0.0, 0.0, 0.0, capPhiAtTow(m), velocityXYZ, rotVelocityXYZ)
                ! Velocity perpendicular to the tower:
                vTowX(m) = rotVelocityXYZ(1)
                ! Velocity tangential to the tower:
                vTowY(m) = rotVelocityXYZ(2)
                ! Local velocity magnitude composed of only the horizontal part of the flow:
                vTow(m) = SQRT(vTowX(m)**2 + vTowY(m)**2)
                ! Calculate thrust force generated by the tower:
                fnTower(m) = 0.5*rhoFstAveLocal(kt)*towDragCoeff(wtID)*(vTow(m)**2) &
                           * towerDiameter(wtID)*towerDs(wtID)
                ! Find global force matrix:
                forceXYZ(1) = fnTower(m)
                forceXYZ(2) = 0.0
                forceXYZ(3) = 0.0
                CALL rotLocalToGlobal(capPhiAtTow(m), 0.0, 0.0, 0.0, forceXYZ, rotForceXYZ)
                ! Assign forces to x,y,z components:
                fxTower(m) = rotForceXYZ(1)
                fyTower(m) = rotForceXYZ(2)
                fzTower(m) = rotForceXYZ(3)
              ENDDO
              posXYZ (1) = 0.0
              posXYZ (2) = (towerDiameter(wtID)/2)
              posXYZ (3) = 0.0
              CALL rotLocalToGlobal(capPhi(kt), 0.0, 0.0, 0.0, posXYZ, rotPosXYZ)
              ! Calculation of the boundaries:
              iSmear(kt) = CEILING( ( ABS( rotPosXYZ(1) ) + projectionRadius(wtID) )/dx )
              jSmear(kt) = CEILING( ( ABS( rotPosXYZ(2) ) + projectionRadius(wtID) )/dy )
              ils = MAX( its, iTow(kt) - iSmear(kt) )
              ile = MIN( ite, iTow(kt) + iSmear(kt) )
              jls = MAX( jts, jTow(kt) - jSmear(kt) )
              jle = MIN( jte, jTow(kt) + jSmear(kt) )
              kls  = kTowBot(wtID) - smoothing
              ! Set lower boundary index to 1 if kls is less than MIN(e_vert)
              IF( kls .LE. kts ) kls = kts
              kle = kTowTop(wtID) + smoothing
              ! loop over j turbine blade levels
              DO jlev = jls, jle
                ! loop over k turbine blade levels
                DO klev = kls, kle
                  ! loop over i turbine blade levels
                  DO ilev = ils, ile
                    fxDist(ilev, klev, jlev) = 0.0
                    fyDist(ilev, klev, jlev) = 0.0
                    fzDist(ilev, klev, jlev) = 0.0
                    ! Loop over tower points:
                    DO m = 1, nTowPoints
                      projectionDistanceU = SQRT( ( (ilev      )*dx - towerPoints(kt, m, 1) )**2 &
                                          +       ( (jlev - 0.5)*dy - towerPoints(kt, m, 2) )**2 &
                                          +       ( (klev - 0.5)*dzConst(wtID)                   &
                                          - towerPoints(kt, m, 3) )**2 )
                      projectionDistanceV = SQRT( ( (ilev - 0.5)*dx - towerPoints(kt, m, 1) )**2 &
                                          +       ( (jlev      )*dy - towerPoints(kt, m, 2) )**2 &
                                          +       ( (klev - 0.5)*dzConst(wtID)                   &
                                          - towerPoints(kt, m, 3) )**2 )
                      projectionDistanceW = SQRT( ( (ilev - 0.5)*dx - towerPoints(kt, m, 1) )**2 &
                                          +       ( (jlev - 0.5)*dy - towerPoints(kt, m, 2) )**2 &
                                          +       ( (klev      )*dzConst(wtID)                   &
                                          - towerPoints(kt, m, 3) )**2 )
                      IF(towerSmearOpt .EQ. 0) THEN ! 3D Gaussian smearing
                        ! x-direction:
                        IF(projectionDistanceU .LE. projectionLimit) THEN
                          etaEpsilonU3D = EXP( -( projectionDistanceU/towerEpsilon(wtID) )**2 )/( (towerEpsilon(wtID)**3)*(piconst**1.5) )
                          coeffU3D = etaEpsilonU3D
                          fxDist(ilev, klev, jlev) = fxDist(ilev, klev, jlev) + fxTower(m)*coeffU3D
                        ENDIF
                        ! y-direction:
                        IF(projectionDistanceV .LE. projectionLimit) THEN
                          etaEpsilonV3D = EXP( -( projectionDistanceV/towerEpsilon(wtID) )**2 )/( (towerEpsilon(wtID)**3)*(piconst**1.5) )
                          coeffV3D = etaEpsilonV3D
                          fyDist(ilev, klev, jlev) = fyDist(ilev, klev, jlev) + fyTower(m)*coeffV3D
                        ENDIF
                        ! z-direction:
                        IF(projectionDistanceW .LE. projectionLimit) THEN
                          etaEpsilonW3D = EXP( -( projectionDistanceW/towerEpsilon(wtID) )**2 )/( (towerEpsilon(wtID)**3)*(piconst**1.5) )
                          coeffW3D = etaEpsilonW3D
                          fzDist(ilev, klev, jlev) = fzDist(ilev, klev, jlev) + fzTower(m)*coeffW3D
                        ENDIF
                      ELSEIF(towerSmearOpt .EQ. 1) THEN ! 3D polynomial smearing
                        ! x-direction:
                        IF( (projectionDistanceU**2) .LE. epsilonMin2(wtID) ) THEN
                          etaEpsilonU3D = ( ( polA(wtID)*(projectionDistanceU**2) - polB(wtID) )*(projectionDistanceU**2) + 1.0 )*epsilonFactor(wtID)
                          coeffU3D = etaEpsilonU3D
                          fxDist(ilev, klev, jlev) = fxDist(ilev, klev, jlev) + fxTower(m)*coeffU3D
                        ENDIF
                        ! y-direction:
                        IF( (projectionDistanceV**2) .LE. epsilonMin2(wtID) ) THEN
                          etaEpsilonV3D = ( ( polA(wtID)*(projectionDistanceV**2) - polB(wtID) )*(projectionDistanceV**2) + 1.0 )*epsilonFactor(wtID)
                          coeffV3D = etaEpsilonV3D
                          fyDist(ilev, klev, jlev) = fyDist(ilev, klev, jlev) + fyTower(m)*coeffV3D
                        ENDIF
                        ! z-direction:
                        IF( (projectionDistanceW**2) .LE. epsilonMin2(wtID) ) THEN
                          etaEpsilonW3D = ( ( polA(wtID)*(projectionDistanceW**2) - polB(wtID) )*(projectionDistanceW**2) + 1.0 )*epsilonFactor(wtID)
                          coeffW3D = etaEpsilonW3D
                          fzDist(ilev, klev, jlev) = fzDist(ilev, klev, jlev) + fzTower(m)*coeffW3D
                        ENDIF
                      ENDIF
                    ENDDO
                    ! Apply the x,y,z-components of the tower forces to the u,v,w-tendencies of the flow:
                    ru_tendf(ilev, klev, jlev) = ru_tendf(ilev, klev, jlev) - muu(ilev, jlev)*fxDist(ilev, klev, jlev)/rhoFstAveLocal(kt)
                    rv_tendf(ilev, klev, jlev) = rv_tendf(ilev, klev, jlev) - muv(ilev, jlev)*fyDist(ilev, klev, jlev)/rhoFstAveLocal(kt)
                    rw_tendf(ilev, klev, jlev) = rw_tendf(ilev, klev, jlev) - mut(ilev, jlev)*fzDist(ilev, klev, jlev)/rhoFstAveLocal(kt)
                  ENDDO
                ENDDO
              ENDDO
              towerPointsLocal(kt, 1:nTowPoints, 1:3) = towerPoints(kt, 1:nTowPoints, 1:3)
            ENDIF
          ELSE
            IF( v0FstLocal(kt) .LT. cutIn(wtID) ) THEN
              WRITE(outstring, "(A, I4, A, F8.3, A)") &
                               'WTP-WARNING ***** Wind turbine # ', kt, ' stopped due to &
                                very low instantaneous wind velocity component normal to &
                                the rotor (v0 = ', v0FstLocal(kt), ' [m/s] < cut-in speed). *****'
              CALL wrf_message ( TRIM(outstring) )
            ELSEIF( v0FstLocal(kt) .GT. cutOut(wtID) ) THEN
              WRITE(outstring, "(A, I4, A, F8.3, A)") &
                               'WTP-WARNING ***** Wind turbine # ', kt, ' stopped due to &
                                very low instantaneous wind velocity component normal to &
                                the rotor (v0 = ', v0FstLocal(kt), ' [m/s] > cut-out speed). *****'
              CALL wrf_message ( TRIM(outstring) )
            ENDIF
            ! If freestream wind speed is below cut-in speed or above cut-out speed of wind turbine
            ! set the variables to -999999.999 for proper reduction to all processors.
            tsrSavedLocal(kt) = -999999.999
            thrSavedLocal(kt) = -999999.999
            trqSavedLocal(kt) = -999999.999
            powSavedLocal(kt) = -999999.999
            mechPowSavedLocal(kt) = -999999.999
            thrCoefSavedLocal(kt) = -999999.999
            trqCoefSavedLocal(kt) = -999999.999
            powCoefSavedLocal(kt) = -999999.999
            rotSpeedLocal(kt) = -999999.999
            pitchSavedLocal(kt) = -999999.999
            fnSavedLocal(kt, :, :) = -999999.999
            ftSavedLocal(kt, :, :) = -999999.999
            FSavedLocal(kt, :, :) = -999999.999
            phiSavedLocal(kt, :, :) = -999999.999
            alphaSavedLocal(kt, :, :) = -999999.999
            clSavedLocal(kt, :, :) = -999999.999
            cdSavedLocal(kt, :, :) = -999999.999
            cmSavedLocal(kt, :, :) = -999999.999
            LSavedLocal(kt, :, :) = -999999.999
            DSavedLocal(kt, :, :) = -999999.999
            v1SavedLocal(kt, :, :) = -999999.999
            vThetaSavedLocal(kt, :, :) = -999999.999
            vTSavedLocal(kt, :, :) = -999999.999
            vRelSavedLocal(kt, :, :) = -999999.999
            vZSavedLocal(kt, :, :) = -999999.999
            uSavedLocal(kt, :, :) = -999999.999
            vSavedLocal(kt, :, :) = -999999.999
            wSavedLocal(kt, :, :) = -999999.999
            IF( speedCtrlOpt(wtID) ) THEN
              trqGenCommandedOldLocal(kt) = -999999.999
              genPowSavedLocal(kt) = -999999.999
              genSpeedOldLocal(kt) = -999999.999
              genSpeedFOldLocal(kt) = -999999.999
            ENDIF
            IF( pitchCtrlOpt(wtID) ) THEN
              speedErrorLocal(kt) = -999999.999
              intspeedErrorLocal(kt) = -999999.999
            ENDIF
            reOperate(kt) = .TRUE.
          ENDIF
          initialYawDone(kt) = .TRUE.
        ELSE
          IF( printDiagOpt .EQ. 1 .OR. printDiagOpt .EQ. 2 ) THEN
            WRITE(outstring, "(A, I4, A, F8.3, A, A, F8.3, A, A, F8.3, A, A, F10.3, A)") &
                             'Wind turbine # ', kt, ' is yawing... No force is &
                              calculated before MOD(ABS(deltaCapPhi), 360.0) &
                              < 90.0 [deg] ---> Wind_dir = ', windDirFstAveLocal(kt),  &
                             ' [deg],', ' trbDirSaved = ', trbDirSaved(kt),            &
                             ' [deg],', ' deltaCapPhi = ', deltaCapPhi(kt), ' [deg],', &
                             ' eYaw = ', eYaw(kt), ' [-]'
            CALL wrf_message ( TRIM(outstring) )
          ENDIF
        ENDIF
      ELSE
        ! These variables are set to a very small number for proper reduction to all processors
        ! when writing into wrfout or wrfrst files.
        iFstLocal(kt) = -999999
        jFstLocal(kt) = -999999
        rhoFstAveLocal(kt) = -999999.999
        capPhiLocal(kt) = -999999.999
        deltaCapPhiLocal(kt) = -999999.999
        eYawLocal(kt) = -999999.999
        uFstAveLocal(kt) = -999999.999
        vFstAveLocal(kt) = -999999.999
        wFstAveLocal(kt) = -999999.999
        windDirFstAveLocal(kt) = -999999.999
        tsrSavedLocal(kt) = -999999.999
        thrSavedLocal(kt) = -999999.999
        trqSavedLocal(kt) = -999999.999
        powSavedLocal(kt) = -999999.999
        mechPowSavedLocal(kt) = -999999.999
        thrCoefSavedLocal(kt) = -999999.999
        trqCoefSavedLocal(kt) = -999999.999
        powCoefSavedLocal(kt) = -999999.999
        wsRotorLocal(kt) = -999999.999
        trbDirLocal(kt) = -999999.999
        rotSpeedLocal(kt) = -999999.999
        v0FstLocal(kt) = -999999.999
        psiOldLocal(kt) = -999999.999
        azimuthLocal(kt) = -999999.999
        deltaAzimuthLocal(kt) = -999999.999
        pitchSavedLocal(kt) = -999999.999
        rhoFstTSLocal(kt, :) = -999999.999
        uFstTSLocal(kt, :) = -999999.999
        vFstTSLocal(kt, :) = -999999.999
        wFstTSLocal(kt, :) = -999999.999
        windDirFstTSLocal(kt, :) = -999999.999
        rotorApexLocal(kt, :) = -999999.999
        fnSavedLocal(kt, :, :) = -999999.999
        ftSavedLocal(kt, :, :) = -999999.999
        FSavedLocal(kt, :, :) = -999999.999
        phiSavedLocal(kt, :, :) = -999999.999
        alphaSavedLocal(kt, :, :) = -999999.999
        clSavedLocal(kt, :, :) = -999999.999
        cdSavedLocal(kt, :, :) = -999999.999
        cmSavedLocal(kt, :, :) = -999999.999
        LSavedLocal(kt, :, :) = -999999.999
        DSavedLocal(kt, :, :) = -999999.999
        v1SavedLocal(kt, :, :) = -999999.999
        vThetaSavedLocal(kt, :, :) = -999999.999
        vTSavedLocal(kt, :, :) = -999999.999
        vRelSavedLocal(kt, :, :) = -999999.999
        vZSavedLocal(kt, :, :) = -999999.999
        uSavedLocal(kt, :, :) = -999999.999
        vSavedLocal(kt, :, :) = -999999.999
        wSavedLocal(kt, :, :) = -999999.999
        bladePointsLocal(kt, :, :, :) = -999999.999
        IF( tsoutPlnOpt .AND. ANY(specTrb .EQ. kt) ) THEN
          nt = nPlnWTcounter
          uPlaneSavedLocal(nt, :, :) = -999999.999
          vPlaneSavedLocal(nt, :, :) = -999999.999
          wPlaneSavedLocal(nt, :, :) = -999999.999
          planePointsLocal(nt, :, :, :) = -999999.999
        ENDIF
        IF( towParOpt ) THEN
          towerPointsLocal(kt, :, :) = -999999.999
        ENDIF
        IF( hubParOpt ) THEN
          nacellePointsLocal(kt, :, :) = -999999.999
        ENDIF
        IF( speedCtrlOpt(wtID) ) THEN
          trqGenCommandedOldLocal(kt) = -999999.999
          genPowSavedLocal(kt) = -999999.999
          genSpeedOldLocal(kt) = -999999.999
          genSpeedFOldLocal(kt) = -999999.999
        ENDIF
        IF( pitchCtrlOpt(wtID) ) THEN
          speedErrorLocal(kt) = -999999.999
          intspeedErrorLocal(kt) = -999999.999
        ENDIF
      ENDIF
      IF( tsoutPlnOpt .AND. ANY(specTrb .EQ. kt) ) THEN
        nt = nPlnWTcounter
        CALL MPI_REDUCE( uPlaneSavedLocal(nt, :, :), uPlaneSaved(nt, :, :), nHorizPoints*nVertPoints, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( vPlaneSavedLocal(nt, :, :), vPlaneSaved(nt, :, :), nHorizPoints*nVertPoints, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( wPlaneSavedLocal(nt, :, :), wPlaneSaved(nt, :, :), nHorizPoints*nVertPoints, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( planePointsLocal(nt, :, :, 1), planePointsX(nt, :, :), nHorizPoints*nVertPoints, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( planePointsLocal(nt, :, :, 2), planePointsY(nt, :, :), nHorizPoints*nVertPoints, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( planePointsLocal(nt, :, :, 3), planePointsZ(nt, :, :), nHorizPoints*nVertPoints, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        nPlnWTcounter = nPlnWTcounter + 1
      ENDIF
      IF( tsoutOpt ) THEN
        IF( specRho .LE. 0.0 ) THEN
          CALL MPI_REDUCE( rhoFstAveLocal(kt), rhoFstAve(kt), 1, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        ENDIF
        CALL MPI_REDUCE( capPhiLocal(kt), capPhi(kt), 1, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( windDirFstAveLocal(kt), windDirFstAve(kt), 1, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( tsrSavedLocal(kt), tsrSaved(kt), 1, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( thrSavedLocal(kt), thrSaved(kt), 1, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( trqSavedLocal(kt), trqSaved(kt), 1, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( powSavedLocal(kt), powSaved(kt), 1, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( mechPowSavedLocal(kt), mechPowSaved(kt), 1, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( thrCoefSavedLocal(kt), thrCoeffSaved(kt), 1, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( trqCoefSavedLocal(kt), trqCoeffSaved(kt), 1, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( powCoefSavedLocal(kt), powCoeffSaved(kt), 1, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( wsRotorLocal(kt), wsRotor(kt), 1, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( trbDirLocal(kt), trbDirSaved(kt), 1, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( rotSpeedLocal(kt), rotSpeed(kt), 1, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( v0FstLocal(kt), v0Fst(kt), 1, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( pitchSavedLocal(kt), pitch(kt), 1, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( rotorApexLocal(kt, 1), rotorApexX(kt), 1, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( rotorApexLocal(kt, 2), rotorApexY(kt), 1, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( rotorApexLocal(kt, 3), rotorApexZ(kt), 1, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        IF( speedCtrlOpt(wtID) ) THEN
          CALL MPI_REDUCE( genPowSavedLocal(kt), genPowSaved(kt), 1, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
          CALL MPI_REDUCE( trqGenCommandedOldLocal(kt), trqGenCommandedOld(kt), 1, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        ENDIF
        CALL MPI_REDUCE( fnSavedLocal(kt, :, :), FnSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( ftSavedLocal(kt, :, :), FtSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( FSavedLocal(kt, :, :), FSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( phiSavedLocal(kt, :, :), phiSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( alphaSavedLocal(kt, :, :), alphaSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( clSavedLocal(kt, :, :), clSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( cdSavedLocal(kt, :, :), cdSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( cmSavedLocal(kt, :, :), cmSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( LSavedLocal(kt, :, :), LSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( DSavedLocal(kt, :, :), DSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( v1SavedLocal(kt, :, :), v1Saved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( vThetaSavedLocal(kt, :, :), vThetaSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( vTSavedLocal(kt, :, :), vTSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( vRelSavedLocal(kt, :, :), vRelSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( vZSavedLocal(kt, :, :), vZSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( uSavedLocal(kt, :, :), uSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( vSavedLocal(kt, :, :), vSaved (kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( wSavedLocal(kt, :, :), wSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( bladePointsLocal(kt, :, :, 1), bladePointsX(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( bladePointsLocal(kt, :, :, 2), bladePointsY(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( bladePointsLocal(kt, :, :, 3), bladePointsZ(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, root, MPI_COMM_WORLD, ierr )
      ENDIF
      IF( wtp_history_flag ) THEN
        CALL MPI_ALLREDUCE( iFstLocal(kt), iFst(kt), 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( jFstLocal(kt), jFst(kt), 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD, ierr )
        IF( specRho .LE. 0.0 ) THEN
          CALL MPI_ALLREDUCE( rhoFstAveLocal(kt), rhoFstAve(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        ENDIF
        CALL MPI_ALLREDUCE( capPhiLocal(kt), capPhi(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( windDirFstAveLocal(kt), windDirFstAve(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( tsrSavedLocal(kt), tsrSaved(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( thrSavedLocal(kt), thrSaved(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( trqSavedLocal(kt), trqSaved(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( powSavedLocal(kt), powSaved(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( mechPowSavedLocal(kt), mechPowSaved(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( thrCoefSavedLocal(kt), thrCoeffSaved(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( trqCoefSavedLocal(kt), trqCoeffSaved(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( powCoefSavedLocal(kt), powCoeffSaved (kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( wsRotorLocal(kt), wsRotor(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( trbDirLocal(kt), trbDirSaved(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( rotSpeedLocal(kt), rotSpeed(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( v0FstLocal(kt), v0Fst(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( pitchSavedLocal(kt), pitch(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( rotorApexLocal(kt, 1), rotorApexX(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( rotorApexLocal(kt, 2), rotorApexY(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( rotorApexLocal(kt, 3), rotorApexZ(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        IF( speedCtrlOpt(wtID) ) THEN
          CALL MPI_ALLREDUCE( genPowSavedLocal(kt), genPowSaved(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        ENDIF
        CALL MPI_ALLREDUCE( fnSavedLocal(kt, :, :), FnSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( ftSavedLocal(kt, :, :), FtSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( FSavedLocal(kt, :, :), FSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( phiSavedLocal(kt, :, :), phiSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( alphaSavedLocal(kt, :, :), alphaSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( clSavedLocal(kt, :, :), clSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( cdSavedLocal(kt, :, :), cdSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( cmSavedLocal(kt, :, :), cmSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( LSavedLocal(kt, :, :), LSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( DSavedLocal(kt, :, :), DSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( v1SavedLocal(kt, :, :), v1Saved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( vThetaSavedLocal(kt, :, :), vThetaSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( vTSavedLocal(kt, :, :), vTSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( vRelSavedLocal(kt, :, :), vRelSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( vZSavedLocal(kt, :, :), vZSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( uSavedLocal(kt, :, :), uSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( vSavedLocal(kt, :, :), vSaved (kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( wSavedLocal(kt, :, :), wSaved(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( bladePointsLocal(kt, :, :, 1), bladePointsX(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( bladePointsLocal(kt, :, :, 2), bladePointsY(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( bladePointsLocal(kt, :, :, 3), bladePointsZ(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        IF( hubParOpt ) THEN
          CALL MPI_ALLREDUCE( nacellePointsLocal(kt, :, 1), nacellePointsX(kt, :), nNacPoints, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
          CALL MPI_ALLREDUCE( nacellePointsLocal(kt, :, 2), nacellePointsY(kt, :), nNacPoints, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
          CALL MPI_ALLREDUCE( nacellePointsLocal(kt, :, 3), nacellePointsZ(kt, :), nNacPoints, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        ENDIF
        IF( towParOpt ) THEN
          CALL MPI_ALLREDUCE( towerPointsLocal(kt, :, 1), towerPointsX(kt, :), nTowPoints, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
          CALL MPI_ALLREDUCE( towerPointsLocal(kt, :, 2), towerPointsY(kt, :), nTowPoints, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
          CALL MPI_ALLREDUCE( towerPointsLocal(kt, :, 3), towerPointsZ(kt, :), nTowPoints, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        ENDIF
      ENDIF
      IF( wtp_restart_flag ) THEN           
        IF( specRho .LE. 0.0 ) THEN
          CALL MPI_ALLREDUCE( rhoFstAveLocal(kt), rhoFstAve(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        ENDIF
        CALL MPI_ALLREDUCE( capPhiLocal(kt), capPhi(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( trbDirLocal(kt), trbDirSaved(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( eYawLocal(kt), eYaw(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( tsrSavedLocal(kt), tsrSaved(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( rotSpeedLocal(kt), rotSpeed(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( azimuthLocal(kt), azimuth(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( deltaAzimuthLocal(kt), deltaAzimuth(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( psiOldLocal(kt), psi(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( pitchSavedLocal(kt), pitch(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( rotorApexLocal(kt, 1), rotorApexX(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( rotorApexLocal(kt, 2), rotorApexY(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( rotorApexLocal(kt, 3), rotorApexZ(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( rhoFstTSLocal(kt, :), rhoFstTS(kt, :), nAveSamples, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( uFstTSLocal(kt, :), uFstTS(kt, :), nAveSamples, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( vFstTSLocal(kt, :), vFstTS(kt, :), nAveSamples, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( wFstTSLocal(kt, :), wFstTS(kt, :), nAveSamples, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( windDirFstTSLocal(kt, :), windDirFstTS(kt, :), nAveSamples, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        IF( speedCtrlOpt(wtID) ) THEN
          CALL MPI_ALLREDUCE( trqGenCommandedOldLocal(kt), trqGenCommandedOld(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
          CALL MPI_ALLREDUCE( genSpeedOldLocal(kt), genSpeedOld(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
          CALL MPI_ALLREDUCE( genSpeedFOldLocal(kt), genSpeedFOld(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        ENDIF
        IF( pitchCtrlOpt(wtID) ) THEN
          CALL MPI_ALLREDUCE( speedErrorLocal(kt), speedError(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
          CALL MPI_ALLREDUCE( intSpeedErrorLocal(kt), intSpeedError(kt), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        ENDIF
        CALL MPI_ALLREDUCE( bladePointsLocal(kt, :, :, 1), bladePointsX(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( bladePointsLocal(kt, :, :, 2), bladePointsY(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( bladePointsLocal(kt, :, :, 3), bladePointsZ(kt, :, :), nElements*nSections, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        IF( hubParOpt ) THEN
          CALL MPI_ALLREDUCE( nacellePointsLocal(kt, :, 1), nacellePointsX(kt, :), nNacPoints, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
          CALL MPI_ALLREDUCE( nacellePointsLocal(kt, :, 2), nacellePointsY(kt, :), nNacPoints, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
          CALL MPI_ALLREDUCE( nacellePointsLocal(kt, :, 3), nacellePointsZ(kt, :), nNacPoints, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        ENDIF
        IF( towParOpt ) THEN
          CALL MPI_ALLREDUCE( towerPointsLocal(kt, :, 1), towerPointsX(kt, :), nTowPoints, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
          CALL MPI_ALLREDUCE( towerPointsLocal(kt, :, 2), towerPointsY(kt, :), nTowPoints, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
          CALL MPI_ALLREDUCE( towerPointsLocal(kt, :, 3), towerPointsZ(kt, :), nTowPoints, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        ENDIF
      ENDIF
#endif
      IF( printDiagOpt .EQ. 1 )THEN
        IF( wrf_dm_on_monitor() ) THEN
          IF( v0FstLocal(kt) .GE. cutIn(wtID) .AND. v0FstLocal(kt) .LE. cutOut(wtID) ) THEN
            WRITE(outstring, "(A, I4, A, A, F8.3, A, A, F7.3, A, A, F8.3, A, A, F7.3, A,A, F10.3, A)") &
                             ' Wind turbine # ', kt, ',',                                              &
                             ' Wind_dir = ', windDirFstAveLocal(kt), ' [deg],',                        &
                             ' v0 = ', v0FstLocal(kt), ' [m/s],',                                      &
                             ' Turbine_dir = ', trbDirSaved(kt), ' [deg],',                            &
                             ' wsRotor  = ', wsRotorLocal(kt), ' [m/s],',                              &
                             ' Yaw_error = ', eYaw(kt), ' [ - ]'
            CALL wrf_message ( TRIM(outstring) )
          ELSEIF( v0FstLocal(kt) .LT. cutIn(wtID) ) THEN
            WRITE(outstring, "(A, I4, A, F8.3, A, A)")                                 &
                             'WTP-WARNING ***** Wind turbine # ', kt, ' stopped due to &
                              very low instantaneous wind velocity component normal to &
                              the rotor (v0 = ', v0FstLocal(kt), ' [m/s]', &
                              ' < cut-in speed). *****'
            CALL wrf_message ( TRIM(outstring) )
          ELSEIF( v0FstLocal(kt) .GT. cutOut(wtID) ) THEN
            WRITE(outstring, "(A, I4, A, F8.3, A, A)")                                 &
                             'WTP-WARNING ***** Wind turbine # ', kt, ' stopped due to &
                              very high instantaneous wind velocity component normal to &
                              the rotor (v0 = ', v0FstLocal(kt), ' [m/s]', &
                              ' > cut-out speed). *****'
            CALL wrf_message ( TRIM(outstring) )
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    restartController = .FALSE.
    readFromRSTFile = 1
  END SUBROUTINE calc_forces_gal
!=======================================================================
! This subroutine interpolates scalars to desired i, j, k levels.
!=======================================================================
  SUBROUTINE interpolateScalar(        &
        &  i, j, k, dx, dy, Lx, Ly, Lz &
        &, val_in, dz, z, val_out      &
        &)
  IMPLICIT NONE
    INTEGER, INTENT(IN)  :: i, j, k
    REAL, DIMENSION(i:i+1, k:k+1, j:j+1), INTENT(IN) :: val_in
    REAL, INTENT(IN) :: dx, dy, dz, Lx, Ly, Lz, z
    REAL, INTENT(OUT) :: val_out
! Local
    REAL :: aa, bb, cc, dd, val_lower, val_upper
    aa = ( (i + 0.5)*dx - Lx )*( (j + 0.5)*dy - Ly )
    bb = ( Lx - (i - 0.5)*dx )*( (j + 0.5)*dy - Ly )
    cc = ( (i + 0.5)*dx - Lx )*( Ly - (j - 0.5)*dy )
    dd = ( Lx - (i - 0.5)*dx )*( Ly - (j - 0.5)*dy )
    val_lower = ( aa*val_in(i, k  , j  ) + bb*val_in(i+1, k  , j  ) &
              +   cc*val_in(i, k  , j+1) + dd*val_in(i+1, k  , j+1) )/(dx*dy)
    val_upper = ( aa*val_in(i, k+1, j  ) + bb*val_in(i+1, k+1, j  ) &
              +   cc*val_in(i, k+1, j+1) + dd*val_in(i+1, k+1, j+1) )/(dx*dy)
    val_out = val_lower + ( (Lz - z)/dz )*(val_upper - val_lower)
  END SUBROUTINE interpolateScalar
!=======================================================================
! This subroutine interpolates velocity components u, v, w to desired
! i, j, k levels.
!=======================================================================
  SUBROUTINE interpolateVelocity(      &
        &  i, j, k, dx, dy, Lx, Ly, Lz &
        &, vel_in, dz, z, interp_field &
        &, vel_out                     &
        &)
  IMPLICIT NONE
    INTEGER, INTENT(IN)  :: i, j, k
    CHARACTER, INTENT(IN) :: interp_field
    REAL, DIMENSION(i:i+1, k:k+1, j:j+1), INTENT(IN) :: vel_in
    REAL, INTENT(IN) :: dx, dy, dz, Lx, Ly, Lz, z
    REAL, INTENT(OUT) :: vel_out
! Local
    REAL :: aa, bb, cc, dd, vel_lower, vel_upper
    IF( interp_field .EQ. 'u' ) THEN
      aa = ( i*dx - Lx )*( (j + 0.5)*dy - Ly )
      bb = ( Lx - (i-1)*dx )*( (j + 0.5)*dy - Ly )
      cc = ( i*dx - Lx )*( Ly - (j - 0.5)*dy )
      dd = ( Lx - (i-1)*dx )*( Ly - (j - 0.5)*dy )
    ELSEIF( interp_field .EQ. 'v' ) THEN
      aa = ( (i + 0.5)*dx - Lx )*( j*dy - Ly )
      bb = ( Lx - (i - 0.5)*dx )*( j*dy - Ly )
      cc = ( (i + 0.5)*dx - Lx )*( Ly - (j-1)*dy )
      dd = ( Lx - (i - 0.5)*dx )*( Ly - (j-1)*dy )
    ELSEIF( interp_field .EQ. 'w' ) THEN
      aa = ( (i + 0.5)*dx - Lx )*( (j + 0.5)*dy - Ly )
      bb = ( Lx - (i - 0.5)*dx )*( (j + 0.5)*dy - Ly )
      cc = ( (i + 0.5)*dx - Lx )*( Ly - (j - 0.5)*dy )
      dd = ( Lx - (i - 0.5)*dx )*( Ly - (j - 0.5)*dy )
    ENDIF
    vel_lower = ( aa*vel_in(i, k  , j  ) + bb*vel_in(i+1, k  , j  ) &
              +   cc*vel_in(i, k  , j+1) + dd*vel_in(i+1, k  , j+1) )/(dx*dy)
    vel_upper = ( aa*vel_in(i, k+1, j  ) + bb*vel_in(i+1, k+1, j  ) &
              +   cc*vel_in(i, k+1, j+1) + dd*vel_in(i+1, k+1, j+1) )/(dx*dy)
    vel_out = vel_lower + ( (Lz - z)/dz )*(vel_upper - vel_lower)
  END SUBROUTINE interpolateVelocity
!=======================================================================
  FUNCTION vectorAdd(a, b)
! This function adds 2 vectors ( arrays real, dimension(3) )
!=======================================================================
  IMPLICIT NONE
    REAL, DIMENSION(3), INTENT(IN) :: a,b
    REAL, DIMENSION(3) :: vectorAdd
    vectorAdd(1) = a(1) + b(1)
    vectorAdd(2) = a(2) + b(2)
    vectorAdd(3) = a(3) + b(3)
    RETURN
  END FUNCTION vectorAdd
!=======================================================================
  FUNCTION vectorDivide(a, b)
! This function divides one vector ( array real, dimension(3) by a
! number )
!=======================================================================
  IMPLICIT NONE
    REAL, DIMENSION(3), INTENT(IN) :: a
    REAL, INTENT(IN) :: b
    REAL, DIMENSION(3) :: vectorDivide
    vectorDivide(1) = a(1)/b
    vectorDivide(2) = a(2)/b
    vectorDivide(3) = a(3)/b
    RETURN
  END FUNCTION vectorDivide
!=======================================================================
  FUNCTION vectorMultiply(a, b)
! This function multiplies one vector (array real, dimension(3) by
! a real number)
!=======================================================================
  IMPLICIT NONE
    REAL, DIMENSION(3), INTENT(IN) :: a
    REAL, INTENT(IN) :: b
    REAL, DIMENSION(3) :: vectorMultiply
    vectorMultiply(1) = a(1)*b
    vectorMultiply(2) = a(2)*b
    vectorMultiply(3) = a(3)*b
    RETURN
  END FUNCTION vectorMultiply
!=======================================================================
  FUNCTION vectorMag(a)
! This function calculates the magnitude of a vector
!=======================================================================
  IMPLICIT NONE
    REAL, DIMENSION(3), INTENT(IN) :: a
    REAL :: vectorMag
    vectorMag = ABS( SQRT( a(1)**2 + a(2)**2 + a(3)**2 ) )
    RETURN
  END FUNCTION vectorMag
!=======================================================================
  FUNCTION rotatePoint(point_in, rotationPoint, axis, angle)
! This function performs rotation of a point with respect to an axis or
! rotation and a certain angle
!=======================================================================
  IMPLICIT NONE
    REAL, DIMENSION(3), INTENT(IN) :: point_in
    REAL, DIMENSION(3), INTENT(IN) :: rotationPoint
    REAL, DIMENSION(3), INTENT(IN) :: axis
    REAL, INTENT(IN) :: angle
    REAL, DIMENSION(3, 3) :: RM ! Rotation Matrix tensor
    REAL, DIMENSION(3) :: rotatePoint, point
    point = point_in
    RM(1, 1) = axis(1) * axis(1) + ( 1.0 - axis(1) * axis(1) )           * COS(angle*degrad)
    RM(1, 2) = axis(1) * axis(2) * ( 1.0 - COS(angle*degrad) ) - axis(3) * SIN(angle*degrad)
    RM(1, 3) = axis(1) * axis(3) * ( 1.0 - COS(angle*degrad) ) + axis(2) * SIN(angle*degrad)
    RM(2, 1) = axis(1) * axis(2) * ( 1.0 - COS(angle*degrad) ) + axis(3) * SIN(angle*degrad)
    RM(2, 2) = axis(2) * axis(2) + ( 1.0 - axis(2) * axis(2) )           * COS(angle*degrad)
    RM(2, 3) = axis(2) * axis(3) * ( 1.0 - COS(angle*degrad) ) - axis(1) * SIN(angle*degrad)
    RM(3, 1) = axis(1) * axis(3) * ( 1.0 - COS(angle*degrad) ) - axis(2) * SIN(angle*degrad)
    RM(3, 2) = axis(2) * axis(3) * ( 1.0 - COS(angle*degrad) ) + axis(1) * SIN(angle*degrad)
    RM(3, 3) = axis(3) * axis(3) + ( 1.0 - axis(3) * axis(3) )           * COS(angle*degrad)
    ! Rotation matrices make a rotation about the origin, so need to subtract
    ! rotation point off the point to be rotated
    point = vectorAdd(point, -rotationPoint)
    ! Perform rotation (multiplication matrix and vector)
    point = matrixVector(RM, point)
    ! Return the rotated point to its new location relative to the rotation point
    rotatePoint = point + rotationPoint
    RETURN
  END FUNCTION rotatePoint
!=======================================================================
  FUNCTION matrixVector(RM, point)
! This function multiplies a matrix and a vector
!=======================================================================
  IMPLICIT NONE
    REAL, DIMENSION(3, 3), INTENT(IN) :: RM ! Matrix
    REAL, DIMENSION(3), INTENT(IN) :: point ! vector point
    REAL, DIMENSION(3) :: matrixVector
    ! Perform rotation
    matrixVector(1) = RM(1, 1)*point(1) + RM(1, 2)*point(2) + RM(1, 3)*point(3)
    matrixVector(2) = RM(2, 1)*point(1) + RM(2, 2)*point(2) + RM(2, 3)*point(3)
    matrixVector(3) = RM(3, 1)*point(1) + RM(3, 2)*point(2) + RM(3, 3)*point(3)
    RETURN
  END FUNCTION matrixVector
!=======================================================================
  FUNCTION crossProduct(u, v)
! This function calculates the cross product of 2 vectors
!=======================================================================
  IMPLICIT NONE
    REAL, INTENT(IN) :: u(3), v(3)
    REAL :: crossProduct(3)
    crossProduct(1) = u(2)*v(3) - u(3)*v(2)
    crossProduct(2) = u(3)*v(1) - u(1)*v(3)
    crossProduct(3) = u(1)*v(2) - u(2)*v(1)
    RETURN
  END FUNCTION crossProduct
!=======================================================================
  SUBROUTINE compassToStandard(dir)
! This subroutine converts wind/turbine direction from compass
! directions to the standard convention of 0 degrees on the + x axis
! with positive degrees in the counter-clockwise direction.
!=======================================================================
  IMPLICIT NONE
    REAL, INTENT(INOUT) :: dir
    dir = dir + 180.0
    IF(dir .GE. 360.0) THEN
      dir = dir - 360.0
    ENDIF
    dir = 90.0 - dir
    IF(dir < 0.0) THEN
      dir = dir + 360.0
    ENDIF
  END SUBROUTINE compassToStandard
!=======================================================================
! This subroutine is used to perform local to global rotation.
!=======================================================================
  SUBROUTINE rotLocalToGlobal(                     &
        &  phi, delta, psi, beta, local_variableIN &
        &, global_variableOUT                      &
        &)
  IMPLICIT NONE
    REAL, INTENT(IN) :: phi, psi, delta, beta
    REAL, DIMENSION(3), INTENT(IN)  :: local_variableIN
    REAL, DIMENSION(3), INTENT(OUT) :: global_variableOUT
! local
    REAL, DIMENSION(3, 3) :: yaw, tilt, azimuth, precone
    REAL, DIMENSION(3, 3) :: rotYawTilt, rotYawTiltAzimuth
    REAL, DIMENSION(3, 3) :: rotYawTiltAzimuthPrecone
    ! Positive rotation around z-axis
    yaw(1, 1:3) = (/ COS(phi*degrad), - SIN(phi*degrad), 0.0 /)
    yaw(2, 1:3) = (/ SIN(phi*degrad),   COS(phi*degrad), 0.0 /)
    yaw(3, 1:3) = (/             0.0,               0.0, 1.0 /)
    ! Positive rotation around y-axis
    tilt(1, 1:3) = (/   COS(delta*degrad), 0.0, SIN(delta*degrad) /)
    tilt(2, 1:3) = (/                 0.0, 1.0,               0.0 /)
    tilt(3, 1:3) = (/ - SIN(delta*degrad), 0.0, COS(delta*degrad) /)
    ! Positive rotation around x-axis
    azimuth(1, 1:3) = (/ 1.0,             0.0,               0.0 /)
    azimuth(2, 1:3) = (/ 0.0, COS(psi*degrad), - SIN(psi*degrad) /)
    azimuth(3, 1:3) = (/ 0.0, SIN(psi*degrad),   COS(psi*degrad) /)
    ! Negative rotation around y-axis
    precone(1, 1:3) = (/ COS(beta*degrad), 0.0, - SIN(beta*degrad) /)
    precone(2, 1:3) = (/              0.0, 1.0,                0.0 /)
    precone(3, 1:3) = (/ SIN(beta*degrad), 0.0,   COS(beta*degrad) /)
    rotYawTilt = MATMUL( yaw, tilt )
    rotYawTiltAzimuth = MATMUL( rotYawTilt, azimuth )
    rotYawTiltAzimuthPrecone = MATMUL( rotYawTiltAzimuth, precone )
    global_variableOUT = MATMUL( rotYawTiltAzimuthPrecone &
                              , local_variableIN )
END SUBROUTINE rotLocalToGlobal
!=======================================================================
! This subroutine is used to perform global to local rotation.
!=======================================================================
  SUBROUTINE rotGlobalToLocal(                      &
        &  beta, psi, delta, phi, global_variableIN &
        &, local_variableOUT                        &
        &)
  IMPLICIT NONE
    REAL, INTENT(IN) :: phi, psi, delta, beta
    REAL, DIMENSION(3), INTENT(IN) :: global_variableIN
    REAL, DIMENSION(3), INTENT(OUT) :: local_variableOUT
! local
    REAL, DIMENSION(3, 3) :: transposeYaw, transposeTilt
    REAL, DIMENSION(3, 3) :: transposeAzimuth, transposePrecone
    REAL, DIMENSION(3, 3) :: rotPreconeAzimuth, rotPreconeAzimuthTilt
    REAL, DIMENSION(3, 3) :: rotPreconeAzimuthTiltYaw
    ! Negative rotation around z-axis
    transposeYaw(1, 1:3) = (/   COS(phi*degrad), SIN(phi*degrad), 0.0 /)
    transposeYaw(2, 1:3) = (/ - SIN(phi*degrad), COS(phi*degrad), 0.0 /)
    transposeYaw(3, 1:3) = (/               0.0,             0.0, 1.0 /)
    ! Negative rotation around y-axis
    transposeTilt(1, 1:3) = (/ COS(delta*degrad), 0.0, - SIN(delta*degrad) /)
    transposeTilt(2, 1:3) = (/               0.0, 1.0,                 0.0 /)
    transposeTilt(3, 1:3) = (/ SIN(delta*degrad), 0.0,   COS(delta*degrad) /)
    ! Negative rotation around x-axis
    transposeAzimuth(1, 1:3) = (/ 1.0,               0.0,             0.0 /)
    transposeAzimuth(2, 1:3) = (/ 0.0,   COS(psi*degrad), SIN(psi*degrad) /)
    transposeAzimuth(3, 1:3) = (/ 0.0, - SIN(psi*degrad), COS(psi*degrad) /)
    ! Positive rotation around y-axis
    transposePrecone(1, 1:3) = (/   COS(beta*degrad), 0.0, SIN(beta*degrad) /)
    transposePrecone(2, 1:3) = (/                0.0, 1.0,              0.0 /)
    transposePrecone(3, 1:3) = (/ - SIN(beta*degrad), 0.0, COS(beta*degrad) /)
    rotPreconeAzimuth = MATMUL( transposePrecone, transposeAzimuth  )
    rotPreconeAzimuthTilt = MATMUL( rotPreconeAzimuth, transposeTilt )
    rotPreconeAzimuthTiltYaw = MATMUL( rotPreconeAzimuthTilt, transposeYaw )
    local_variableOUT = MATMUL( rotPreconeAzimuthTiltYaw &
                              , global_variableIN )
  END SUBROUTINE rotGlobalToLocal
!=======================================================================
! This subroutine calculates aerodynamic losses of the turbine blades.
!=======================================================================
  SUBROUTINE calcAeroLosses(             &
        &  tSectRad, phi, TSR, nBlade    &
        &, rTip, rHub, aerodynLossOpt, F &
        &)
  IMPLICIT NONE
    INTEGER, INTENT(IN) :: aerodynLossOpt, nBlade
    REAL, INTENT(IN) :: tSectRad, phi, TSR, rTip, rHub
    REAL, INTENT(OUT) :: F
! Local
    REAL :: f_tip, f_hub, Ftip, Fhub
    REAL :: c1, c2, c3, g
    f_tip = nBlade*(rTip - tSectRad)/( 2*tSectRad*ABS( SIN(phi*degrad) ) )
    f_hub = nBlade*(tSectRad - rHub)/( 2*tSectRad*ABS( SIN(phi*degrad) ) )
    IF( aerodynLossOpt .GE. 0 .AND. aerodynLossOpt .LE. 3 ) THEN ! Glauert root/tip correction
      Ftip = (2/piconst)*ACOS( EXP(-f_tip) )
      Fhub = (2/piconst)*ACOS( EXP(-f_hub) )
    ELSE ! Shen et al. (2005) tip correction
      c1   =  0.125
      c2   = 21.000
      c3   =  0.100
      g    = EXP( -c1*( (nBlade*TSR) - c2 ) ) + c3
      Ftip = (2/piconst)*ACOS( EXP(-g*f_tip) )
    ENDIF
    IF( aerodynLossOpt .EQ. 0 ) THEN
      F = 1.0
    ELSEIF( aerodynLossOpt .EQ. 1 ) THEN
      F = Fhub
    ELSEIF( aerodynLossOpt .EQ. 2 .OR. aerodynLossOpt .EQ. 4 ) THEN
      F = Ftip
    ELSEIF( aerodynLossOpt .EQ. 3 ) THEN
      F = Ftip*Fhub
    ENDIF
  END SUBROUTINE calcAeroLosses
!=======================================================================
! This subroutine applies three-dimensional corrections to the airfoils.
!=======================================================================
  SUBROUTINE correction3D( tSectRad, chord, omega, v0, TSR  &
        &, Vr, vT, theta, alphaZero, slope, correction3DOpt &
        &, rTip, Fcl, Fcd                                   )
  IMPLICIT NONE
    INTEGER, INTENT(IN) :: correction3DOpt
    REAL, INTENT(IN) :: tSectRad, chord, theta, alphaZero
    REAL, INTENT(IN) :: omega, v0, TSR, slope, Vr, vT, rTip
    REAL, INTENT(OUT) :: Fcd, Fcl
! Local
    REAL :: a, b, d, h, n, lambda
    lambda = (omega*rTip)/SQRT(v0**2 + vT**2)
    IF( alphaZero .EQ. 0.0 ) THEN
      Fcl = 0.0
      Fcd = 0.0
    ELSE
      IF( correction3DOpt .EQ. 1 ) THEN ! Three-dimensional correction for lift from Snel et al. (1993)
        Fcl = 3.0*( (chord/tSectRad)**2 )
        Fcd = 0.0
      ELSEIF( correction3DOpt .EQ. 2 ) THEN ! Three-dimensional correction for lift from Lindenburg (2004)
        Fcl = 3.1*( (omega*tSectRad/Vr) )**2*( (chord/tSectRad)**2 )
        Fcd = 0.0
      ELSEIF( correction3DOpt .EQ. 3 .OR. correction3DOpt .EQ. 4 ) THEN ! Three-dimensional corrections for lift and drag from Du and Selig (1998)
        a = 1.0
        b = 1.0
        d = 1.0
        Fcl = (1/slope)*( ( 1.6*(chord/tSectRad)/0.1267 ) &
            * ( a - (chord/tSectRad)**( (d*rTip)/(lambda*tSectRad) ) ) &
            / ( b + (chord/tSectRad)**( (d*rTip)/(lambda*tSectRad) ) ) - 1 )
        Fcd = - (1/slope)*( ( 1.6*(chord/tSectRad)/0.1267 ) &
            * ( a - (chord/tSectRad)**( (d*rTip)/(2*lambda*tSectRad) ) ) &
            / ( b + (chord/tSectRad)**( (d*rTip)/(2*lambda*tSectRad) ) ) - 1 )
      ELSEIF( correction3DOpt .EQ. 5 ) THEN ! Three-dimensional corrections for lift and drag from Chaviaropoulos and Hansen (2000)
        a = 2.2
        h = 1
        n = 4
        Fcl = a*( (chord/tSectRad)**h )*( COS(theta*degrad)**n )
        Fcd = Fcl
      ENDIF
    ENDIF
  END SUBROUTINE correction3D
!=======================================================================
! This subroutine interpolates a real variable w.r.t to the input vector.
!=======================================================================
  SUBROUTINE interpolate(VectLen, xNew, xOld, yOld, yNew)
  IMPLICIT NONE
    INTEGER, INTENT(IN) :: VectLen
    REAL, INTENT(IN) :: xNew
    REAL, DIMENSION(VectLen), INTENT(IN) :: xOld, yOld
    REAL, INTENT(OUT) :: yNew
! local
    INTEGER :: index, indexP, indexM
    index = MINLOC(ABS(xNew - xOld), DIM=1)
    IF( xNew .LT. xOld(index) ) THEN
      IF( index .EQ. 1 ) THEN
        indexP = index + 1
        indexM = index
      ELSE
        indexP = index
        indexM = index - 1
      ENDIF
      yNew = yOld(indexM) + ( ( yOld(indexP) - yOld(indexM) ) &
           / ( xOld(indexP) - xOld(indexM) ) )*( xNew - xOld(indexM) )
    ELSEIF( xNew .GT. xOld(index) ) THEN
      IF( index .EQ. SIZE(xOld) ) THEN
        indexP = SIZE(xOld)
        indexM = indexP - 1
      ELSE
        indexP = index + 1
        indexM = index
      ENDIF
      yNew = yOld(indexM) + ( ( yOld(indexP) - yOld(indexM) ) &
           / ( xOld(indexP) - xOld(indexM) ) )*( xNew - xOld(indexM) )
    ELSEIF( xNew .EQ. xOld(index) ) THEN
      yNew = yOld(index)
    ELSE
      yNew = 0.0
    ENDIF
END SUBROUTINE interpolate
!=======================================================================
! This subroutine is used to initialize Generalized Actuator Line module.
!=======================================================================
  SUBROUTINE init_module_wind_gal(            &
        &  grid, config_flags, restart        &
        &, xlong_u, xlat_v, ht, dx, dy        &
        &, ids, ide, jds, jde, kds, kde       &
        &, ims, ime, jms, jme, kms, kme       &
        &, its, ite, jts, jte, kts, kte       &
#ifdef DM_PARALLEL
        &, nprocX                             &
#endif
        &, dz8w, u_2, v_2, capPhi             &
        &, rotorApexX, rotorApexY, rotorApexZ &
        &, bladePointsX, bladePointsY         &
        &, bladePointsZ                       &
        &, towerPointsX, towerPointsY         &
        &, towerPointsZ                       &
        &, nacellePointsX, nacellePointsY     &
        &, nacellePointsZ                     &
        &, planePointsX, planePointsY         &
        &, planePointsZ                       &
        &)
!
  IMPLICIT NONE
#if defined(DM_PARALLEL) && !defined(STUBMPI)
  INCLUDE "mpif.h"
#endif
    TYPE (domain) :: grid
    TYPE (grid_config_rec_type) :: config_flags
    TYPE (PROJ_INFO) :: ts_proj
!
    REAL, INTENT(IN) :: dx, dy
!
    INTEGER, INTENT(IN) :: its, ite, jts, jte, kts, kte
    INTEGER, INTENT(IN) :: ims, ime, jms, jme, kms, kme
    INTEGER, INTENT(IN) :: ids, ide, jds, jde, kds, kde
!
#ifdef DM_PARALLEL
    INTEGER, INTENT(IN) :: nprocX
#endif
!
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: xlong_u, xlat_v, ht
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: dz8w, u_2, v_2
!
    REAL, DIMENSION(1:config_flags%wind_wtp_nTurbines), INTENT(IN) :: rotorApexX
    REAL, DIMENSION(1:config_flags%wind_wtp_nTurbines), INTENT(IN) :: rotorApexY
    REAL, DIMENSION(1:config_flags%wind_wtp_nTurbines), INTENT(IN) :: rotorApexZ
    REAL, DIMENSION(1:config_flags%wind_wtp_nTurbines), INTENT(IN) :: capPhi
!
    REAL, DIMENSION(1:config_flags%wind_wtp_nTurbines, 1:config_flags%wind_wtp_nSections &
                  , 1:config_flags%wind_wtp_nElements), INTENT(IN) :: bladePointsX
    REAL, DIMENSION(1:config_flags%wind_wtp_nTurbines, 1:config_flags%wind_wtp_nSections &
                  , 1:config_flags%wind_wtp_nElements), INTENT(IN) :: bladePointsY
    REAL, DIMENSION(1:config_flags%wind_wtp_nTurbines, 1:config_flags%wind_wtp_nSections &
                  , 1:config_flags%wind_wtp_nElements), INTENT(IN) :: bladePointsZ
!
    REAL, DIMENSION(1:config_flags%wind_wtp_nTurbines, 1:config_flags%wind_wtp_nTowpoints) &
                  , INTENT(IN) :: towerPointsX
    REAL, DIMENSION(1:config_flags%wind_wtp_nTurbines, 1:config_flags%wind_wtp_nTowpoints) &
                  , INTENT(IN) :: towerPointsY
    REAL, DIMENSION(1:config_flags%wind_wtp_nTurbines, 1:config_flags%wind_wtp_nTowpoints) &
                  , INTENT(IN) :: towerPointsZ
!
    REAL, DIMENSION(1:config_flags%wind_wtp_nTurbines, 1:config_flags%wind_wtp_nNacpoints) &
                  , INTENT(IN) :: nacellePointsX
    REAL, DIMENSION(1:config_flags%wind_wtp_nTurbines, 1:config_flags%wind_wtp_nNacpoints) &
                  , INTENT(IN) :: nacellePointsY
    REAL, DIMENSION(1:config_flags%wind_wtp_nTurbines, 1:config_flags%wind_wtp_nNacpoints) &
                  , INTENT(IN) :: nacellePointsZ
!
    REAL, DIMENSION(1:config_flags%wind_wtp_pln_nTurbines, &
                    1:config_flags%wind_wtp_pln_nz,        &
                    1:config_flags%wind_wtp_pln_nxy        ), INTENT(IN) :: planePointsX
    REAL, DIMENSION(1:config_flags%wind_wtp_pln_nTurbines, &
                    1:config_flags%wind_wtp_pln_nz,        &
                    1:config_flags%wind_wtp_pln_nxy        ), INTENT(IN) :: planePointsY
    REAL, DIMENSION(1:config_flags%wind_wtp_pln_nTurbines, &
                    1:config_flags%wind_wtp_pln_nz,        &
                    1:config_flags%wind_wtp_pln_nxy        ), INTENT(IN) :: planePointsZ
!
    LOGICAL, INTENT(IN) :: restart
! Local
    REAL :: dist, ts_rx, ts_ry, known_lat, known_lon
!
    real :: towerToPlane
!
    INTEGER :: kt, nVal0, nVal1, nVal2, nVal3, nVal4, nVal5
    INTEGER :: i, j, e, k, l, m, n, nt, p, o, nSize
    INTEGER :: itf, jtf, ktf
!
    CHARACTER*256 :: input0, input1, input2, input3, input4, input5
    CHARACTER*256 :: message_windturbine, outstring
!
    LOGICAL, EXTERNAL :: wrf_dm_on_monitor
!
    REAL, DIMENSION(:), ALLOCATABLE :: windDir, capPhi_init, zeta
    REAL, DIMENSION(:), ALLOCATABLE :: zMaxTrb, bladeRad, bladeChord
    REAL, DIMENSION(:), ALLOCATABLE :: bladeThickness, bladeTwist
    REAL, DIMENSION(:), ALLOCATABLE :: trbWindSpeed, trbPitch, trbRotSpeed, trbPow, trbThr
    REAL, DIMENSION(:), ALLOCATABLE :: bladeUpper, bladeLower, hubLower, hubUpper
#ifdef DM_PARALLEL
    INTEGER, DIMENSION(:), ALLOCATABLE :: procsCommTmp, procsTmp
    INTEGER, DIMENSION(:), ALLOCATABLE :: arrayIN, arrayOUT
#endif
    INTEGER, DIMENSION(:), ALLOCATABLE ::  kEnd, turbineType, airfoilType, nFoilTypes
!
    CHARACTER*256, DIMENSION(:), ALLOCATABLE :: turbineName, airfoilName
    CHARACTER*256, DIMENSION(:), ALLOCATABLE :: turbineNameVect
    CHARACTER*256, DIMENSION(:), ALLOCATABLE :: foilNameOUT, foilNameReduced
    CHARACTER*256, DIMENSION(:), ALLOCATABLE :: foilNameVect
!
    REAL, DIMENSION(:, :), ALLOCATABLE :: foilTypeREAL, foilTypeVectREAL
    REAL, DIMENSION(:, :), ALLOCATABLE :: bladeAoA, bladeLift
    REAL, DIMENSION(:, :), ALLOCATABLE :: bladeDrag, bladeMoment
!
    INTEGER, DIMENSION(:, :), ALLOCATABLE ::  foilTypeVect
!
    REAL, DIMENSION(its-1:ite+1, kts:kte, jts-1:jte+1) :: w_z ! z-coordinate of w-nodes (m) ! Grid point (node)
    REAL, DIMENSION(its-1:ite+1, kts:kte, jts-1:jte+1) :: dz  ! Thickness of each grid cell (m)
!=======================================================================
    IF( restart ) THEN
      restartController = .TRUE.
      restartController1 = .TRUE.
    ENDIF
    ! namelist.input options for the WTP model:
    bladeSmearOpt = config_flags%wind_wtp_bladesmear_opt
    nacelleSmearOpt = config_flags%wind_wtp_nacellesmear_opt
    towerSmearOpt = config_flags%wind_wtp_towersmear_opt 
    yawOpt = config_flags%wind_wtp_yaw_opt
    aerodynLossOpt = config_flags%wind_wtp_aerodynloss_opt
    correction3DOpt = config_flags%wind_wtp_correction3D_opt
    hubParOpt = config_flags%wind_wtp_hubpar_opt
    towParOpt = config_flags%wind_wtp_towpar_opt
    printDiagOpt = config_flags%wind_wtp_printdiag_opt
    errorCheckOpt = config_flags%wind_wtp_errorcheck_opt
    nTurbines = config_flags%wind_wtp_nTurbines
    nSections = config_flags%wind_wtp_nSections
    nElements = config_flags%wind_wtp_nElements
    nTowPoints = config_flags%wind_wtp_nTowpoints
    nNacPoints = config_flags%wind_wtp_nNacpoints
    nAveSamples = config_flags%wind_wtp_nAvesamples
    tsoutOpt = config_flags%wind_wtp_tsout_opt
    specRho = config_flags%wind_wtp_spec_rho
    specDir = config_flags%wind_wtp_spec_turbdir
!
    tsoutPlnOpt = config_flags%wind_wtp_pln_tsout_opt
    nPlnWT = config_flags%wind_wtp_pln_nTurbines
    nHorizPoints = config_flags%wind_wtp_pln_nxy
    nVertPoints = config_flags%wind_wtp_pln_nz
    planeWidth = config_flags%wind_wtp_pln_width
    planeHeight = config_flags%wind_wtp_pln_height
    planeDistance = config_flags%wind_wtp_pln_dist
    IF( wrf_dm_on_monitor() ) THEN
! ... Open the file with the location of the wind turbines ...
      IF( config_flags%windturbines_latlon .EQ. 0 ) THEN
        OPEN(70, FILE='windturbines-ij.dat', FORM='FORMATTED', STATUS='OLD')
      ELSE
        OPEN(70, FILE='windturbines-latlon.dat', FORM='FORMATTED', STATUS='OLD')
      ENDIF
 100  CONTINUE
      REWIND (70)
    ENDIF
! ... When doing spin-up runs in real configuration, real.exe sometimes initializes the finest domain
! ... with GAD/GAL twice. For that reason, first, deallocate all allocatable arrays to work around the
! ... allocated array issue.
    IF( ALLOCATED(specTrb) ) DEALLOCATE(specTrb)
    IF( ALLOCATED(Lx) ) DEALLOCATE(Lx)
    IF( ALLOCATED(Ly) ) DEALLOCATE(Ly)
    IF( ALLOCATED(Lz) ) DEALLOCATE(Lz)
    IF( ALLOCATED(turbineTypeID) ) DEALLOCATE(turbineTypeID)
    IF( ALLOCATED(iTow) ) DEALLOCATE(iTow)
    IF( ALLOCATED(jTow) ) DEALLOCATE(jTow)
    IF( ALLOCATED(lat) ) DEALLOCATE(lat)
    IF( ALLOCATED(lon) ) DEALLOCATE(lon)
    IF( ALLOCATED(hubHeight) ) DEALLOCATE(hubHeight)
    IF( ALLOCATED(hubHeightASL) ) DEALLOCATE(hubHeightASL)
    IF( ALLOCATED(hubHeightASLLocal) ) DEALLOCATE(hubHeightASLLocal)
    IF( ALLOCATED(turbineTypeVect) ) DEALLOCATE(turbineTypeVect)
    IF( ALLOCATED(turbineNameVect) ) DEALLOCATE(turbineNameVect)
    IF( ALLOCATED(speedCtrlOpt) ) DEALLOCATE(speedCtrlOpt)
    IF( ALLOCATED(pitchCtrlOpt) ) DEALLOCATE(pitchCtrlOpt)
    IF( ALLOCATED(cutIn) ) DEALLOCATE(cutIn)
    IF( ALLOCATED(cutOut) ) DEALLOCATE(cutOut)
    IF( ALLOCATED(operDataVectLen) ) DEALLOCATE(operDataVectLen)
    IF( ALLOCATED(nBlade) ) DEALLOCATE(nBlade)
    IF( ALLOCATED(hubHeight) ) DEALLOCATE(hubHeight)
    IF( ALLOCATED(towerDiameter) ) DEALLOCATE(towerDiameter)
    IF( ALLOCATED(rotorDiameter) ) DEALLOCATE(rotorDiameter)
    IF( ALLOCATED(hubDiameter) ) DEALLOCATE(hubDiameter)
    IF( ALLOCATED(nacelleLength) ) DEALLOCATE(nacelleLength)
    IF( ALLOCATED(nacDragCoeff) ) DEALLOCATE(nacDragCoeff)
    IF( ALLOCATED(towDragCoeff) ) DEALLOCATE(towDragCoeff)
    IF( ALLOCATED(overhang) ) DEALLOCATE(overhang)
    IF( ALLOCATED(undSling) ) DEALLOCATE(undSling)
    IF( ALLOCATED(mechEfficiency) ) DEALLOCATE(mechEfficiency)
    IF( ALLOCATED(tiltAngle) ) DEALLOCATE(tiltAngle)
    IF( ALLOCATED(preconeAngle) ) DEALLOCATE(preconeAngle)
    IF( ALLOCATED(nacelleCepsilon) ) DEALLOCATE(nacelleCepsilon)
    IF( ALLOCATED(towerCepsilon) ) DEALLOCATE(towerCepsilon)
    IF( ALLOCATED(inflowLocation) ) DEALLOCATE(inflowLocation)
    IF( ALLOCATED(yawRate) ) DEALLOCATE(yawRate)
    IF( ALLOCATED(thresholdYaw) ) DEALLOCATE(thresholdYaw)
    IF( ALLOCATED(trbDataVectLen) ) DEALLOCATE(trbDataVectLen)
    IF( ALLOCATED(rHub) ) DEALLOCATE(rHub)
    IF( ALLOCATED(rTip) ) DEALLOCATE(rTip)
    IF( ALLOCATED(rTow) ) DEALLOCATE(rTow)
    IF( ALLOCATED(foilTypeVect) ) DEALLOCATE(foilTypeVect)
    IF( ALLOCATED(foilNameVect) ) DEALLOCATE(foilNameVect)
    IF( ALLOCATED(nFoilTypes) ) DEALLOCATE(nFoilTypes)
    IF( ALLOCATED(foilPropVectLen) ) DEALLOCATE(foilPropVectLen)
    IF( ALLOCATED(ratedRotSpeed) ) DEALLOCATE(ratedRotSpeed)
    IF( ALLOCATED(GBEfficiency) ) DEALLOCATE(GBEfficiency)
    IF( ALLOCATED(GBRatio) ) DEALLOCATE(GBRatio)
    IF( ALLOCATED(ratedGenPow) ) DEALLOCATE(ratedGenPow)
    IF( ALLOCATED(genIner) ) DEALLOCATE(genIner)
    IF( ALLOCATED(hubIner) ) DEALLOCATE(hubIner)
    IF( ALLOCATED(bladeIner) ) DEALLOCATE(bladeIner)
    IF( ALLOCATED(genEfficiency) ) DEALLOCATE(genEfficiency)
    IF( ALLOCATED(rateLimitGenTorque) ) DEALLOCATE(rateLimitGenTorque)
    IF( ALLOCATED(cutInGenSpeed) ) DEALLOCATE(cutInGenSpeed)
    IF( ALLOCATED(region2StartGenSpeed) ) DEALLOCATE(region2StartGenSpeed)
    IF( ALLOCATED(region2EndGenSpeed) ) DEALLOCATE(region2EndGenSpeed)
    IF( ALLOCATED(initGenTorque) ) DEALLOCATE(initGenTorque)
    IF( ALLOCATED(cutInGenTorque) ) DEALLOCATE(cutInGenTorque)
    IF( ALLOCATED(ratedGenTorque) ) DEALLOCATE(ratedGenTorque)
    IF( ALLOCATED(maxGenTorque) ) DEALLOCATE(maxGenTorque)
    IF( ALLOCATED(KGen) ) DEALLOCATE(KGen)
    IF( ALLOCATED(speedFilterCornerFrequency) ) DEALLOCATE(speedFilterCornerFrequency)
    IF( ALLOCATED(pitchRate) ) DEALLOCATE(pitchRate)
    IF( ALLOCATED(rateLimitBladePitch) ) DEALLOCATE(rateLimitBladePitch)
    IF( ALLOCATED(pitchMin) ) DEALLOCATE(pitchMin)
    IF( ALLOCATED(pitchMax) ) DEALLOCATE(pitchMax)
    IF( ALLOCATED(pitchK) ) DEALLOCATE(pitchK)
    IF( ALLOCATED(pitchControlKP) ) DEALLOCATE(pitchControlKP)
    IF( ALLOCATED(pitchControlKI) ) DEALLOCATE(pitchControlKI)
    IF( ALLOCATED(pitchControlKD) ) DEALLOCATE(pitchControlKD)
    IF( ALLOCATED(speedControlFileExists) ) DEALLOCATE(speedControlFileExists)
    IF( ALLOCATED(pitchControlFileExists) ) DEALLOCATE(pitchControlFileExists)
    IF( ALLOCATED(v1Vect) ) DEALLOCATE(v1Vect)
    IF( ALLOCATED(thetapVect) ) DEALLOCATE(thetapVect)
    IF( ALLOCATED(omegaVect) ) DEALLOCATE(omegaVect)
    IF( ALLOCATED(aeroPowVect) ) DEALLOCATE(aeroPowVect)
    IF( ALLOCATED(aeroThrVect) ) DEALLOCATE(aeroThrVect)
    IF( ALLOCATED(bladeCepsilon) ) DEALLOCATE(bladeCepsilon)
    IF( ALLOCATED(bladeEpsilon) ) DEALLOCATE(bladeEpsilon)
    IF( ALLOCATED(rVect) ) DEALLOCATE(rVect)
    IF( ALLOCATED(chordVect) ) DEALLOCATE(chordVect)
    IF( ALLOCATED(thicknessVect) ) DEALLOCATE(thicknessVect)
    IF( ALLOCATED(thetaGVect) ) DEALLOCATE(thetaGVect)
    IF( ALLOCATED(alphaVect) ) DEALLOCATE(alphaVect)
    IF( ALLOCATED(clVect) ) DEALLOCATE(clVect)
    IF( ALLOCATED(cdVect) ) DEALLOCATE(cdVect)
    IF( ALLOCATED(cmVect) ) DEALLOCATE(cmVect)
    IF( ALLOCATED(foilPropVectLen) ) DEALLOCATE(foilPropVectLen)
    IF( ALLOCATED(bladeRadius) ) DEALLOCATE(bladeRadius)
    IF( ALLOCATED(c) ) DEALLOCATE(c)
    IF( ALLOCATED(t) ) DEALLOCATE(t)
    IF( ALLOCATED(thetaS) ) DEALLOCATE(thetaS)
    IF( ALLOCATED(foilType) ) DEALLOCATE(foilType)
    IF( ALLOCATED(foilTypeREAL) ) DEALLOCATE(foilTypeREAL)
    IF( ALLOCATED(foilTypeVectREAL) ) DEALLOCATE(foilTypeVectREAL)
    IF( ALLOCATED(doYaw) ) DEALLOCATE(doYaw)
    IF( ALLOCATED(doYawLocal) ) DEALLOCATE(doYawLocal)
    IF( ALLOCATED(reOperate) ) DEALLOCATE(reOperate)
    IF( ALLOCATED(initialYaw) ) DEALLOCATE(initialYaw)
    IF( ALLOCATED(initialYawLocal) ) DEALLOCATE(initialYawLocal)
    IF( ALLOCATED(initialYawDone) ) DEALLOCATE(initialYawDone)
    IF( ALLOCATED(hubLower) ) DEALLOCATE(hubLower)
    IF( ALLOCATED(hubUpper) ) DEALLOCATE(hubUpper)
    IF( ALLOCATED(bladeUpper) ) DEALLOCATE(bladeUpper)
    IF( ALLOCATED(bladeLower) ) DEALLOCATE(bladeLower)
    IF( ALLOCATED(speedErrorInit) ) DEALLOCATE(speedErrorInit)
    IF( ALLOCATED(trqGenInit) ) DEALLOCATE(trqGenInit)
    IF( ALLOCATED(rotSpeedInit) ) DEALLOCATE(rotSpeedInit)
    IF( ALLOCATED(pitchInit) ) DEALLOCATE(pitchInit)
    IF( ALLOCATED(speedErrorLocal) ) DEALLOCATE(speedErrorLocal)
    IF( ALLOCATED(intSpeedErrorLocal) ) DEALLOCATE(intSpeedErrorLocal)
    IF( ALLOCATED(GK) ) DEALLOCATE(GK )
    IF( ALLOCATED(speedErrorLast) ) DEALLOCATE(speedErrorLast)
    IF( ALLOCATED(derivSpeedError) ) DEALLOCATE(derivSpeedError)
    IF( ALLOCATED(pitchP) ) DEALLOCATE(pitchP)
    IF( ALLOCATED(pitchI) ) DEALLOCATE(pitchI)
    IF( ALLOCATED(pitchD) ) DEALLOCATE(pitchD)
    IF( ALLOCATED(pitchCommanded) ) DEALLOCATE(pitchCommanded)
    IF( ALLOCATED(pitchRateCommanded) ) DEALLOCATE(pitchRateCommanded)
    IF( ALLOCATED(deltaAzimuthOld) ) DEALLOCATE(deltaAzimuthOld)
    IF( ALLOCATED(deltaAzimuthLocal) ) DEALLOCATE(deltaAzimuthLocal)
    IF( ALLOCATED(psiOld) ) DEALLOCATE(psiOld)
    IF( ALLOCATED(psiOldLocal) ) DEALLOCATE(psiOldLocal)
    IF( ALLOCATED(azimuthOld) ) DEALLOCATE(azimuthOld)
    IF( ALLOCATED(azimuthLocal) ) DEALLOCATE(azimuthLocal)
    IF( ALLOCATED(rotSpeedLocal) ) DEALLOCATE(rotSpeedLocal)
    IF( ALLOCATED(rotSpeedF) ) DEALLOCATE(rotSpeedF)
    IF( ALLOCATED(driveTrainIner) ) DEALLOCATE(driveTrainIner)
    IF( ALLOCATED(windDirFstAveLocal) ) DEALLOCATE(windDirFstAveLocal)
    IF( ALLOCATED(uFstAveLocal) ) DEALLOCATE(uFstAveLocal)
    IF( ALLOCATED(vFstAveLocal) ) DEALLOCATE(vFstAveLocal)
    IF( ALLOCATED(wFstAveLocal) ) DEALLOCATE(wFstAveLocal)
    IF( ALLOCATED(rhoFstTSLocal) ) DEALLOCATE(rhoFstTSLocal)
    IF( ALLOCATED(uFstTSLocal) ) DEALLOCATE(uFstTSLocal)
    IF( ALLOCATED(vFstTSLocal) ) DEALLOCATE(vFstTSLocal)
    IF( ALLOCATED(wFstTSLocal) ) DEALLOCATE(wFstTSLocal)
    IF( ALLOCATED(windDirFstTSLocal) ) DEALLOCATE(windDirFstTSLocal)
#ifdef DM_PARALLEL
    IF( ALLOCATED(sender) ) DEALLOCATE(sender)
    IF( ALLOCATED(procID) ) DEALLOCATE(procID)
    IF( ALLOCATED(procInit) ) DEALLOCATE(procInit)
    IF( ALLOCATED(procIDInit) ) DEALLOCATE(procIDInit)
    IF( ALLOCATED(procIDFstInit) ) DEALLOCATE(procIDFstInit)
    IF( ALLOCATED(ProcFstInit) ) DEALLOCATE(ProcFstInit)
    IF( ALLOCATED(procIDFst) ) DEALLOCATE(procIDFst)
    IF( ALLOCATED(proc) ) DEALLOCATE(proc)
    IF( ALLOCATED(procFst) ) DEALLOCATE(procFst)
    IF( ALLOCATED(newGroup) ) DEALLOCATE(newGroup)
    IF( ALLOCATED(newComm) ) DEALLOCATE(newComm)
#endif
    IF( ALLOCATED(iRot) ) DEALLOCATE(iRot)
    IF( ALLOCATED(jRot) ) DEALLOCATE(jRot)
    IF( ALLOCATED(iFstLocal) ) DEALLOCATE(iFstLocal)
    IF( ALLOCATED(jFstLocal) ) DEALLOCATE(jFstLocal)
    IF( ALLOCATED(iFstInit) ) DEALLOCATE(iFstInit)
    IF( ALLOCATED(jFstInit) ) DEALLOCATE(jFstInit)
    IF( ALLOCATED(kFst) ) DEALLOCATE(kFst)
    IF( ALLOCATED(dr) ) DEALLOCATE(dr)
    IF( ALLOCATED(delta) ) DEALLOCATE(delta)
    IF( ALLOCATED(windDir) ) DEALLOCATE(windDir)
    IF( ALLOCATED(uvShaftDir) ) DEALLOCATE(uvShaftDir)
    IF( ALLOCATED(towerToRotorApex) ) DEALLOCATE(towerToRotorApex)
    IF( ALLOCATED(undSlingRotor) ) DEALLOCATE(undSlingRotor)
    IF( ALLOCATED(xMinTile) ) DEALLOCATE(xMinTile)
    IF( ALLOCATED(xMaxTile) ) DEALLOCATE(xMaxTile)
    IF( ALLOCATED(yMinTile) ) DEALLOCATE(yMinTile)
    IF( ALLOCATED(yMaxTile) ) DEALLOCATE(yMaxTile)
    IF( ALLOCATED(xMinTrb) ) DEALLOCATE(xMinTrb)
    IF( ALLOCATED(xMaxTrb) ) DEALLOCATE(xMaxTrb)
    IF( ALLOCATED(yMinTrb) ) DEALLOCATE(yMinTrb)
    IF( ALLOCATED(yMaxTrb) ) DEALLOCATE(yMaxTrb)
    IF( ALLOCATED(deltaYawLocal) ) DEALLOCATE(deltaYawLocal)
    IF( ALLOCATED(deltaYaw) ) DEALLOCATE(deltaYaw)
    IF( ALLOCATED(rotorApexInit) ) DEALLOCATE(rotorApexInit)
    IF( ALLOCATED(rotorApexLocal) ) DEALLOCATE(rotorApexLocal)
    IF( ALLOCATED(posRotorApex) ) DEALLOCATE(posRotorApex)
    IF( ALLOCATED(towerBaseLocation) ) DEALLOCATE(towerBaseLocation)
    IF( ALLOCATED(uvShaft) ) DEALLOCATE(uvShaft)
    IF( ALLOCATED(uvTower) ) DEALLOCATE(uvTower)
    IF( ALLOCATED(towerShaftIntersect) ) DEALLOCATE(towerShaftIntersect)
    IF( ALLOCATED(bladeRoot) ) DEALLOCATE(bladeRoot)
    IF( ALLOCATED(capPhi_init) ) DEALLOCATE(capPhi_init)
!
    IF( ALLOCATED(planePoints) ) DEALLOCATE(planePoints)
    IF( ALLOCATED(planePointsLocal) ) DEALLOCATE(planePointsLocal)
!
    IF( ALLOCATED(bladePoints) ) DEALLOCATE(bladePoints)
    IF( ALLOCATED(bladePointsLocal) ) DEALLOCATE(bladePointsLocal)
    IF( ALLOCATED(towerPoints) ) DEALLOCATE(towerPoints)
    IF( ALLOCATED(towerPointsLocal) ) DEALLOCATE(towerPointsLocal)
    IF( ALLOCATED(nacellePoints) ) DEALLOCATE(nacellePoints)
    IF( ALLOCATED(nacellePointsLocal) ) DEALLOCATE(nacellePointsLocal)
    IF( ALLOCATED(rhoFstAveLocal) ) DEALLOCATE(rhoFstAveLocal)
    IF( ALLOCATED(trbDirLocal) ) DEALLOCATE(trbDirLocal)
    IF( ALLOCATED(beta) ) DEALLOCATE(beta)
    IF( ALLOCATED(zeta) ) DEALLOCATE(zeta)
    IF( ALLOCATED(wsRotorLocal) ) DEALLOCATE(wsRotorLocal)
    IF( ALLOCATED(v0FstLocal) ) DEALLOCATE(v0FstLocal)
    IF( ALLOCATED(capPhiLocal) ) DEALLOCATE(capPhiLocal)
    IF( ALLOCATED(eYawLocal) ) DEALLOCATE(eYawLocal)
    IF( ALLOCATED(deltaCapPhiLocal) ) DEALLOCATE(deltaCapPhiLocal)
    IF( ALLOCATED(tsrSavedLocal) ) DEALLOCATE(tsrSavedLocal)
    IF( ALLOCATED(thrSavedLocal) ) DEALLOCATE(thrSavedLocal)
    IF( ALLOCATED(trqSavedLocal) ) DEALLOCATE(trqSavedLocal)
    IF( ALLOCATED(powSavedLocal) ) DEALLOCATE(powSavedLocal)
    IF( ALLOCATED(mechPowSavedLocal) ) DEALLOCATE(mechPowSavedLocal)
    IF( ALLOCATED(genPowSavedLocal) ) DEALLOCATE(genPowSavedLocal)
    IF( ALLOCATED(thrCoefSavedLocal) ) DEALLOCATE(thrCoefSavedLocal)
    IF( ALLOCATED(trqCoefSavedLocal) ) DEALLOCATE(trqCoefSavedLocal)
    IF( ALLOCATED(powCoefSavedLocal) ) DEALLOCATE(powCoefSavedLocal)
    IF( ALLOCATED(trqGenCommandedOldOld) ) DEALLOCATE(trqGenCommandedOldOld)
    IF( ALLOCATED(trqGenCommandedOldLocal) ) DEALLOCATE(trqGenCommandedOldLocal)
    IF( ALLOCATED(genSpeedOldLocal) ) DEALLOCATE(genSpeedOldLocal)
    IF( ALLOCATED(genSpeedFOldOld) ) DEALLOCATE(genSpeedFOldOld)
    IF( ALLOCATED(genSpeedFOldLocal) ) DEALLOCATE(genSpeedFOldLocal)
    IF( ALLOCATED(pitchOld) ) DEALLOCATE(pitchOld)
    IF( ALLOCATED(pitchLocal) ) DEALLOCATE(pitchLocal)
    IF( ALLOCATED(pitchSavedLocal) ) DEALLOCATE(pitchSavedLocal)
    IF( ALLOCATED(fnSavedLocal) ) DEALLOCATE(fnSavedLocal)
    IF( ALLOCATED(ftSavedLocal) ) DEALLOCATE(ftSavedLocal)
    IF( ALLOCATED(FSavedLocal) ) DEALLOCATE(FSavedLocal)
    IF( ALLOCATED(phiSavedLocal) ) DEALLOCATE(phiSavedLocal)
    IF( ALLOCATED(alphaSavedLocal) ) DEALLOCATE(alphaSavedLocal)
    IF( ALLOCATED(clSavedLocal) ) DEALLOCATE(clSavedLocal)
    IF( ALLOCATED(cdSavedLocal) ) DEALLOCATE(cdSavedLocal)
    IF( ALLOCATED(cmSavedLocal) ) DEALLOCATE(cmSavedLocal)
    IF( ALLOCATED(LSavedLocal) ) DEALLOCATE(LSavedLocal)
    IF( ALLOCATED(DSavedLocal) ) DEALLOCATE(DSavedLocal)
    IF( ALLOCATED(v1SavedLocal) ) DEALLOCATE(v1SavedLocal)
    IF( ALLOCATED(vTSavedLocal) ) DEALLOCATE(vTSavedLocal)
    IF( ALLOCATED(vRelSavedLocal) ) DEALLOCATE(vRelSavedLocal)
    IF( ALLOCATED(vThetaSavedLocal) ) DEALLOCATE(vThetaSavedLocal)
    IF( ALLOCATED(vZSavedLocal) ) DEALLOCATE(vZSavedLocal)
    IF( ALLOCATED(uSavedLocal) ) DEALLOCATE(uSavedLocal)
    IF( ALLOCATED(vSavedLocal) ) DEALLOCATE(vSavedLocal)
    IF( ALLOCATED(wSavedLocal) ) DEALLOCATE(wSavedLocal)
    IF( ALLOCATED(solidity) ) DEALLOCATE(solidity)
    IF( ALLOCATED(nacelleEpsilon) ) DEALLOCATE(nacelleEpsilon)
    IF( ALLOCATED(towerEpsilon) ) DEALLOCATE(towerEpsilon)
    IF( ALLOCATED(projectionRadius) ) DEALLOCATE(projectionRadius)
    IF( ALLOCATED(epsilonMin) ) DEALLOCATE(epsilonMin)
    IF( ALLOCATED(epsilonMin2) ) DEALLOCATE(epsilonMin2)
    IF( ALLOCATED(polA) ) DEALLOCATE(polA)
    IF( ALLOCATED(polB) ) DEALLOCATE(polB)
    IF( ALLOCATED(epsilonFactor) ) DEALLOCATE(epsilonFactor)
    IF( ALLOCATED(xLocRot) ) DEALLOCATE(xLocRot)
    IF( ALLOCATED(yLocRot) ) DEALLOCATE(yLocRot)
    IF( ALLOCATED(xLocHub) ) DEALLOCATE(xLocHub)
    IF( ALLOCATED(yLocHub) ) DEALLOCATE(yLocHub)
    IF( ALLOCATED(xLocTow) ) DEALLOCATE(xLocTow)
    IF( ALLOCATED(yLocTow) ) DEALLOCATE(yLocTow)
    IF( ALLOCATED(kEnd) ) DEALLOCATE(kEnd)
    IF( ALLOCATED(towerDs) ) DEALLOCATE(towerDs)
    IF( ALLOCATED(nacelleDs) ) DEALLOCATE(nacelleDs)
    IF( ALLOCATED(contribution) ) DEALLOCATE(contribution)
    IF( ALLOCATED(zmaxtrb) ) DEALLOCATE(zmaxtrb)
    IF( ALLOCATED(kTrbBot) ) DEALLOCATE(kTrbBot)
    IF( ALLOCATED(kTrbHub) ) DEALLOCATE(kTrbHub)
    IF( ALLOCATED(kTrbTop) ) DEALLOCATE(kTrbTop)
    IF( ALLOCATED(kHubBot) ) DEALLOCATE(kHubBot)
    IF( ALLOCATED(kHubTop) ) DEALLOCATE(kHubTop)
    IF( ALLOCATED(kTowBot) ) DEALLOCATE(kTowBot)
    IF( ALLOCATED(kTowTop) ) DEALLOCATE(kTowTop)
    IF( ALLOCATED(dzconst) ) DEALLOCATE(dzconst)
!
    IF( ALLOCATED(uPlaneSavedLocal) ) DEALLOCATE(uPlaneSavedLocal)
    IF( ALLOCATED(vPlaneSavedLocal) ) DEALLOCATE(vPlaneSavedLocal)
    IF( ALLOCATED(wPlaneSavedLocal) ) DEALLOCATE(wPlaneSavedLocal)
! ... Initialize the configuration of the wind turbine(s) ...
    ALLOCATE( Lx(1:nTurbines) )
    ALLOCATE( Ly(1:nTurbines) )
    ALLOCATE( Lz(1:nTurbines) )
    ALLOCATE( turbineTypeID(1:nTurbines) )
    ALLOCATE( iTow(1:nTurbines) )
    ALLOCATE( jTow(1:nTurbines) )
    ALLOCATE( lat(1:nTurbines) )
    ALLOCATE( lon(1:nTurbines) )
    ALLOCATE( hubHeightASL(1:nTurbines) )
    ALLOCATE( hubHeightASLLocal(1:nTurbines) )
    IF( wrf_dm_on_monitor() ) THEN
      DO k = 1, nTurbines
        IF( config_flags%windturbines_latlon .EQ. 0 ) THEN
          READ(70, *) Lx(k), Ly(k), Lz(k), turbineTypeID(k)
          iTow(k) = INT( (Lx(k) + 0.5*dx)/dx )
          jTow(k) = INT( (Ly(k) + 0.5*dy)/dy )
        ELSE
          READ(70, *) lat(k), lon(k), Lz(k), turbineTypeID(k)
          CALL map_init(ts_proj)
          known_lat = xlat_v(its, jts)
          known_lon = xlong_u(its, jts)
          ! Mercator
          IF(config_flags%map_proj == PROJ_MERC) THEN
             CALL map_set(PROJ_MERC, ts_proj,               &
                          truelat1 = config_flags%truelat1, &
                          lat1     = known_lat,             &
                          lon1     = known_lon,             &
                          knowni   = REAL(its),             &
                          knownj   = REAL(jts),             &
                          dx       = config_flags%dx)
          ! Lambert conformal
          ELSEIF(config_flags%map_proj == PROJ_LC) THEN
             CALL map_set(PROJ_LC, ts_proj,                  &
                          truelat1 = config_flags%truelat1,  &
                          truelat2 = config_flags%truelat2,  &
                          stdlon   = config_flags%stand_lon, &
                          lat1     = known_lat,              &
                          lon1     = known_lon,              &
                          knowni   = REAL(its),              &
                          knownj   = REAL(jts),              &
                          dx       = config_flags%dx)
!         ! Polar stereographic
          ELSEIF(config_flags%map_proj == PROJ_PS) THEN
             CALL map_set(PROJ_PS, ts_proj,                  &
                          truelat1 = config_flags%truelat1,  &
                          stdlon   = config_flags%stand_lon, &
                          lat1     = known_lat,              &
                          lon1     = known_lon,              &
                          knowni   = REAL(its),              &
                          knownj   = REAL(jts),              &
                          dx       = config_flags%dx)
#if (EM_CORE == 1)
          ! Cassini (global ARW)
          ELSEIF(config_flags%map_proj == PROJ_CASSINI) THEN
             CALL map_set(PROJ_CASSINI, ts_proj,                          &
                          latinc   = grid%dy*360.0/(2*EARTH_RADIUS_M*PI), &
                          loninc   = grid%dx*360.0/(2*EARTH_RADIUS_M*PI), &
                          lat1     = known_lat,                           &
                          lon1     = known_lon,                           &
                          lat0     = config_flags%pole_lat,               &
                          lon0     = config_flags%pole_lon,               &
                          knowni   = 1.,                                  &
                          knownj   = 1.,                                  &
                          stdlon   = config_flags%stand_lon)
#endif
          ! Rotated latitude-longitude
          ELSEIF (config_flags%map_proj == PROJ_ROTLL) THEN
             CALL map_set(PROJ_ROTLL, ts_proj,                      &
! I have no idea how this should work for NMM nested domains
                          ixdim    = grid%e_we-1,                   &
                          jydim    = grid%e_sn-1,                   &
                          phi      = REAL(grid%e_sn-2)*grid%dy/2,   &
                          lambda   = REAL(grid%e_we-2)*grid%dx,     &
                          lat1     = config_flags%cen_lat,          &
                          lon1     = config_flags%cen_lon,          &
                          latinc   = grid%dy,                       &
                          loninc   = grid%dx,                       &
                          stagger  = HH)
          ENDIF
          CALL latlon_to_ij(ts_proj, lat(k), lon(k), ts_rx, ts_ry)
          Lx(k) = ts_rx*dx
          Ly(k) = ts_ry*dy
          iTow(k) = INT( (Lx(k) + 0.5*dx)/dx )
          jTow(k) = INT( (Ly(k) + 0.5*dy)/dy )
        ENDIF
      ENDDO
      CLOSE(70)
    ENDIF
#ifdef DM_PARALLEL
    CALL wrf_dm_bcast_integer( turbineTypeID, nTurbines )
#endif
    CALL removeDupsINT(turbineTypeID, turbineTypeVect, nTurbineTypes)
    input0 = "windTurbines/windTurbineTypes.dat"
    OPEN(FILE=input0, UNIT=14, FORM='FORMATTED', STATUS='OLD')
    nVal0 = MAXVAL(turbineTypeID)
    ALLOCATE( turbineType(1:nVal0), turbineName(1:nVal0) )
    turbineType(:) = -1
    DO n = 1, nVal0
      READ(14, *, ERR=109) turbineType(n), turbineName(n)
    ENDDO
    CLOSE(14)
    ALLOCATE( turbineNameVect(1:nVal0) )
    turbineNameVect = turbineName(1:nVal0)
    DEALLOCATE( turbineType, turbineName )
 109  CONTINUE  
    ALLOCATE( speedCtrlOpt(1:nTurbineTypes) )
    ALLOCATE( pitchCtrlOpt(1:nTurbineTypes) )
    ALLOCATE( cutIn(1:nTurbineTypes) )
    ALLOCATE( cutOut(1:nTurbineTypes) )
    ALLOCATE( operDataVectLen(1:nTurbineTypes) )
    ALLOCATE( v1Vect(1:nTurbineTypes, 1:maxVals2) )
    ALLOCATE( thetapVect(1:nTurbineTypes, 1:maxVals2) )
    ALLOCATE( omegaVect(1:nTurbineTypes, 1:maxVals2) )
    ALLOCATE( aeroPowVect(1:nTurbineTypes, 1:maxVals2) )
    ALLOCATE( aeroThrVect(1:nTurbineTypes, 1:maxVals2) )
    ALLOCATE( bladeCepsilon(1:nTurbineTypes, 1:3) )
    ALLOCATE( bladeEpsilon(1:nTurbineTypes, 1:3) )
    ALLOCATE( nBlade(1:nTurbineTypes) )
    ALLOCATE( hubHeight(1:nTurbineTypes) )
    ALLOCATE( towerDiameter(1:nTurbineTypes) )
    ALLOCATE( rotorDiameter(1:nTurbineTypes) )
    ALLOCATE( hubDiameter(1:nTurbineTypes) )
    ALLOCATE( nacelleLength(1:nTurbineTypes) )
    ALLOCATE( nacDragCoeff(1:nTurbineTypes) )
    ALLOCATE( towDragCoeff(1:nTurbineTypes) )
    ALLOCATE( overhang(1:nTurbineTypes) )
    ALLOCATE( undSling(1:nTurbineTypes) )
    ALLOCATE( mechEfficiency(1:nTurbineTypes) )
    ALLOCATE( tiltAngle(1:nTurbineTypes) )
    ALLOCATE( preconeAngle(1:nTurbineTypes) )
    ALLOCATE( nacelleCepsilon(1:nTurbineTypes) )
    ALLOCATE( towerCepsilon(1:nTurbineTypes) )
    ALLOCATE( inflowLocation(1:nTurbineTypes) )
    ALLOCATE( yawRate(1:nTurbineTypes) )
    ALLOCATE( thresholdYaw(1:nTurbineTypes) )
    ALLOCATE( rVect(1:nTurbineTypes, 1:maxVals2) )
    ALLOCATE( chordVect(1:nTurbineTypes, 1:maxVals2) )
    ALLOCATE( thicknessVect(1:nTurbineTypes, 1:maxVals2) )
    ALLOCATE( thetaGVect(1:nTurbineTypes, 1:maxVals2) )
    ALLOCATE( trbDataVectLen(1:nTurbineTypes) )
    ALLOCATE( rHub(1:nTurbineTypes) )
    ALLOCATE( rTip(1:nTurbineTypes) )
    ALLOCATE( rTow(1:nTurbineTypes) )
    ALLOCATE( foilTypeVect(1:nTurbineTypes, 1:maxVals2) )
    ALLOCATE( foilNameVect(1:maxVals2) )
    ALLOCATE( nFoilTypes(1:nTurbineTypes) )
    ALLOCATE( alphaVect(1:nTurbineTypes, 1:maxVals1, 1:maxVals2) )
    ALLOCATE( clVect(1:nTurbineTypes, 1:maxVals1, 1:maxVals2) )
    ALLOCATE( cdVect(1:nTurbineTypes, 1:maxVals1, 1:maxVals2) )
    ALLOCATE( cmVect(1:nTurbineTypes, 1:maxVals1, 1:maxVals2) )
    ALLOCATE( foilPropVectLen(1:nTurbineTypes, 1:maxVals1) )
    ALLOCATE( ratedRotSpeed(1:nTurbineTypes) )
    ALLOCATE( GBEfficiency(1:nTurbineTypes) )
    ALLOCATE( GBRatio(1:nTurbineTypes) )
    ALLOCATE( ratedGenPow(1:nTurbineTypes) )
    ALLOCATE( genIner(1:nTurbineTypes) )
    ALLOCATE( hubIner(1:nTurbineTypes) )
    ALLOCATE( bladeIner(1:nTurbineTypes) )
    ALLOCATE( genEfficiency(1:nTurbineTypes) )
    ALLOCATE( rateLimitGenTorque(1:nTurbineTypes) )
    ALLOCATE( cutInGenSpeed(1:nTurbineTypes) )
    ALLOCATE( region2StartGenSpeed(1:nTurbineTypes) )
    ALLOCATE( region2EndGenSpeed(1:nTurbineTypes) )
    ALLOCATE( initGenTorque(1:nTurbineTypes) )
    ALLOCATE( cutInGenTorque(1:nTurbineTypes) )
    ALLOCATE( ratedGenTorque(1:nTurbineTypes) )
    ALLOCATE( maxGenTorque(1:nTurbineTypes) )
    ALLOCATE( KGen(1:nTurbineTypes) )
    ALLOCATE( speedFilterCornerFrequency(1:nTurbineTypes) )
    ALLOCATE( pitchRate(1:nTurbineTypes) )
    ALLOCATE( rateLimitBladePitch(1:nTurbineTypes) )
    ALLOCATE( pitchMin(1:nTurbineTypes) )
    ALLOCATE( pitchMax(1:nTurbineTypes) )
    ALLOCATE( pitchK(1:nTurbineTypes) )
    ALLOCATE( pitchControlKP(1:nTurbineTypes) )
    ALLOCATE( pitchControlKI(1:nTurbineTypes) )
    ALLOCATE( pitchControlKD(1:nTurbineTypes) )
    ALLOCATE( speedControlFileExists(1:nTurbineTypes) )
    ALLOCATE( pitchControlFileExists(1:nTurbineTypes) )
    ALLOCATE( solidity(1:nTurbineTypes) )
    ALLOCATE( nacelleEpsilon(1:nTurbineTypes) )
    ALLOCATE( towerEpsilon(1:nTurbineTypes) )
    ALLOCATE( projectionRadius(1:nTurbineTypes) )
    ALLOCATE( epsilonMin(1:nTurbineTypes) )
    ALLOCATE( epsilonMin2(1:nTurbineTypes) )
    ALLOCATE( polA(1:nTurbineTypes) )
    ALLOCATE( polB(1:nTurbineTypes) )
    ALLOCATE( epsilonFactor(1:nTurbineTypes) )
    ALLOCATE( xLocRot(1:nTurbineTypes) )
    ALLOCATE( yLocRot(1:nTurbineTypes) )
    ALLOCATE( xLocHub(1:nTurbineTypes) )
    ALLOCATE( yLocHub(1:nTurbineTypes) )
    ALLOCATE( xLocTow(1:nTurbineTypes) )
    ALLOCATE( yLocTow(1:nTurbineTypes) )
    ALLOCATE( kEnd(1:nTurbineTypes) )
!=======================================================================
! Initialize allocated variables:
!=======================================================================
    cutIn(:) = 0.0
    cutOut(:) = 0.0
    operDataVectLen(:) = 0.0
    v1Vect(:, :) = 0.0
    thetapVect(:, :) = 0.0
    omegaVect(:, :) = 0.0
    aeroPowVect(:, :) = 0.0
    aeroThrVect(:, :) = 0.0
    bladeCepsilon(:, :) = 0.0
    bladeEpsilon(:, :) = 0.0
    nBlade(:) = 0.0
    hubHeight(:) = 0.0
    towerDiameter(:) = 0.0
    rotorDiameter(:) = 0.0
    hubDiameter(:) = 0.0
    nacelleLength(:) = 0.0
    nacDragCoeff(:) = 0.0
    towDragCoeff(:) = 0.0
    overhang(:) = 0.0
    undSling(:) = 0.0
    mechEfficiency(:) = 0.0
    tiltAngle(:) = 0.0
    preconeAngle(:) = 0.0
    nacelleCepsilon(:) = 0.0
    towerCepsilon(:) = 0.0
    inflowLocation(:) = 0.0
    yawRate(:) = 0.0
    thresholdYaw(:) = 0.0
    rVect(:, :) = 0.0
    chordVect(:, :) = 0.0
    thicknessVect(:, :) = 0.0
    thetaGVect(:, :) = 0.0
    trbDataVectLen(:) = 0.0
    rHub(:) = 0.0
    rTip(:) = 0.0
    rTow(:) = 0.0
    foilTypeVect(:, :) = -1
    nFoilTypes(:) = -1
    alphaVect(:, :, :) = 0.0
    clVect(:, :, :) = 0.0
    cdVect(:, :, :) = 0.0
    cmVect(:, :, :) = 0.0
    foilPropVectLen(:, :) = 0
    ratedRotSpeed(:) = 0.0
    GBEfficiency(:) = 0.0
    GBRatio(:) = 0.0
    ratedGenPow(:) = 0.0
    genIner(:) = 0.0
    hubIner(:) = 0.0
    bladeIner(:) = 0.0
    genEfficiency(:) = 0.0
    rateLimitGenTorque(:) = 0.0
    cutInGenSpeed(:) = 0.0
    region2StartGenSpeed(:) = 0.0
    region2EndGenSpeed(:) = 0.0
    initGenTorque(:) = 0.0
    cutInGenTorque(:) = 0.0
    ratedGenTorque(:) = 0.0
    maxGenTorque(:) = 0.0
    KGen(:) = 0.0
    speedFilterCornerFrequency(:) = 0.0
    pitchRate(:) = 0.0
    rateLimitBladePitch(:) = 0.0
    pitchMin(:) = 0.0
    pitchMax(:) = 0.0
    pitchK(:) = 0.0
    pitchControlKP(:) = 0.0
    pitchControlKI(:) = 0.0
    pitchControlKD(:) = 0.0
    speedControlFileExists(:) = .FALSE.
    pitchControlFileExists(:) = .FALSE.
    speedCtrlOpt(:) = config_flags%wind_wtp_speedcontrol_opt
    pitchCtrlOpt(:) = config_flags%wind_wtp_pitchcontrol_opt
    solidity(:) = 0.0
    nacelleEpsilon(:) = 0.0
    towerEpsilon(:) = 0.0
    projectionRadius(:) = 0.0
    epsilonMin(:) = 0.0
    epsilonMin2(:) = 0.0
    polA(:) = 0.0
    polB(:) = 0.0
    epsilonFactor(:) = 0.0
    xLocRot(:) = 0.0
    yLocRot(:) = 0.0
    xLocHub(:) = 0.0
    yLocHub(:) = 0.0
    xLocTow(:) = 0.0
    yLocTow(:) = 0.0
    kEnd(:) = 0.0
    DO n = 1, nTurbineTypes
      input1 = "windTurbines/"//TRIM( turbineNameVect( turbineTypeVect(n) ) )//"/operationalData.tbl"
      OPEN(FILE=input1, UNIT=15, FORM='FORMATTED', STATUS='OLD')
      READ(15, *, ERR=110) nVal1
      ALLOCATE( trbWindSpeed(1:nVal1) )
      ALLOCATE( trbPitch(1:nVal1) )
      ALLOCATE( trbRotSpeed(1:nVal1) )
      ALLOCATE( trbPow(1:nVal1) )
      ALLOCATE( trbThr(1:nVal1) )
      trbWindSpeed(:) = 0.0
      trbPitch(:) = 0.0
      trbRotSpeed(:) = 0.0
      trbPow(:) = 0.0
      trbThr(:) = 0.0
      DO l = 1, nVal1
        READ(15, *, ERR=110) trbWindSpeed(l), trbPitch(l), trbRotSpeed(l), &
                             trbPow(l), trbThr(l)
      ENDDO
      CLOSE(15)
      operDataVectLen(n) = nVal1
      cutIn(n) = trbWindSpeed(1)
      cutOut(n) = trbWindSpeed( operDataVectLen(n) )
      v1Vect(n, 1:nVal1) = trbWindSpeed( 1:operDataVectLen(n) )
      thetapVect(n, 1:nVal1) = trbPitch( 1:operDataVectLen(n) )
      omegaVect(n, 1:nVal1) = trbRotSpeed( 1:operDataVectLen(n) )
      aeroPowVect(n, 1:nVal1) = trbPow( 1:operDataVectLen(n) )
      aeroThrVect(n, 1:nVal1) = trbThr(1 :operDataVectLen(n) )
      DEALLOCATE( trbWindSpeed, trbPitch, trbRotSpeed, trbPow, trbThr )
 110  CONTINUE
      input2 = "windTurbines/"//TRIM( turbineNameVect( turbineTypeVect(n) ) )//"/turbineProperties.tbl"
      OPEN(FILE=input2, UNIT=16, FORM='FORMATTED', STATUS='OLD')
      READ(16, *, ERR=111) nVal2
      READ(16, *, ERR=111) nVal3
      READ(16, *, ERR=111) nBlade(n)
      READ(16, *, ERR=111) hubHeight(n)
      READ(16, *, ERR=111) towerDiameter(n)
      READ(16, *, ERR=111) rotorDiameter(n)
      READ(16, *, ERR=111) hubDiameter(n)
      READ(16, *, ERR=111) nacelleLength(n)
      READ(16, *, ERR=111) nacDragCoeff(n)
      READ(16, *, ERR=111) towDragCoeff(n)
      READ(16, *, ERR=111) overhang(n)
      READ(16, *, ERR=111) undSling(n)
      READ(16, *, ERR=111) mechEfficiency(n)
      READ(16, *, ERR=111) tiltAngle(n)
      READ(16, *, ERR=111) preconeAngle(n)
      READ(16, *, ERR=111) bladeCepsilon(n, 1)
      READ(16, *, ERR=111) bladeCepsilon(n, 2)
      READ(16, *, ERR=111) bladeCepsilon(n, 3)
      READ(16, *, ERR=111) nacelleCepsilon(n)
      READ(16, *, ERR=111) towerCepsilon(n)
      READ(16, *, ERR=111) inflowLocation(n)
      READ(16, *, ERR=111) yawRate(n)
      READ(16, *, ERR=111) thresholdYaw(n)
      ALLOCATE( bladeRad(1:nVal2) )
      ALLOCATE( bladeChord(1:nVal2) )
      ALLOCATE( bladeThickness(1:nVal2) )
      ALLOCATE( bladeTwist(1:nVal2) )
      ALLOCATE( airfoilName(1:nVal2), airfoilType(1:nVal2) )
      bladeRad(:) = 0.0
      bladeChord(:) = 0.0
      bladeThickness(:) = 0.0
      bladeTwist(:) = 0.0
      airfoilType(:) = -1
      DO m = 1, nVal2
        READ(16, *, ERR=111) bladeRad(m), bladeChord(m), &
                             bladeThickness(m), bladeTwist(m), &
                             airfoilName(m), airfoilType(m)
      ENDDO
      CLOSE(16)
      trbDataVectLen(n) = nVal2
      rHub(n) = hubDiameter(n)/2
      rTip(n) = rotorDiameter(n)/2
      rTow(n) = towerDiameter(n)/2
      rVect(n, 1:nVal2) = bladeRad( 1:trbDataVectLen(n) )
      chordVect(n, 1:nVal2) = bladeChord( 1:trbDataVectLen(n) )
      thicknessVect(n, 1:nVal2) = bladeThickness( 1:trbDataVectLen(n) )
      thetaGVect(n, 1:nVal2) = bladeTwist( 1:trbDataVectLen(n) )
      foilTypeVect(n, 1:nVal2) = airfoilType(1:nVal2)
      foilNameVect(1:nVal2) = airfoilName(1:nVal2)
      nFoilTypes(n) = nVal3
      DEALLOCATE( bladeRad, bladeChord, bladeThickness, bladeTwist )
      DEALLOCATE( airfoilType, airfoilName )
 111  CONTINUE
      CALL removeDupsSTR(foilNameVect, foilNameOUT, nSize)
      foilNameReduced = foilNameOUT(1:nSize)
      ALLOCATE( bladeAoA(1:nFoilTypes(n), 1:maxVals2) )
      ALLOCATE( bladeLift(1:nFoilTypes(n), 1:maxVals2) )
      ALLOCATE( bladeDrag (1:nFoilTypes(n), 1:maxVals2) )
      ALLOCATE( bladeMoment(1:nFoilTypes(n), 1:maxVals2) )
      bladeAoA(:, :) = 0.0
      bladeLift(:, :) = 0.0
      bladeDrag(:, :) = 0.0
      bladeMoment(:, :) = 0.0
      DO o = 1, nFoilTypes(n)
        OPEN(UNIT=17, FILE="windTurbines/"//TRIM( turbineNameVect( turbineTypeVect(n) ) )//"/"//TRIM(foilNameReduced(o))//".dat", FORM='FORMATTED', STATUS='OLD')
        READ(17, *, ERR=112) nVal4
        DO p = 1, nVal4
          READ(17, *, ERR=112) bladeAoA(o, p), bladeLift(o, p), bladeDrag(o, p), bladeMoment(o, p)
        ENDDO
        CLOSE(17)
        alphaVect(n, o, 1:nVal4) = bladeAoA(o, 1:nVal4)
        clVect(n, o, 1:nVal4) = bladeLift(o, 1:nVal4)
        cdVect(n, o, 1:nVal4) = bladeDrag(o, 1:nVal4)
        cmVect(n, o, 1:nVal4) = bladeMoment(o, 1:nVal4)
        foilPropVectLen(n, o) = nVal4
      ENDDO
      DEALLOCATE( bladeAoA, bladeLift, bladeDrag, bladeMoment )
 112  CONTINUE
      CALL fileCheck( "windTurbines/"//TRIM( turbineNameVect( turbineTypeVect(n) ) )//"/speedController.tbl", speedControlFileExists(n) )
      IF( speedCtrlOpt(n) .AND. speedControlFileExists(n) ) THEN
        input3 = "windTurbines/"//TRIM( turbineNameVect( turbineTypeVect(n) ) )//"/speedController.tbl"
        OPEN(FILE=input3, UNIT=18, FORM='FORMATTED', STATUS='OLD')
        READ(18, *, ERR=113) ratedRotSpeed(n)
        READ(18, *, ERR=113) GBEfficiency(n)
        READ(18, *, ERR=113) GBRatio(n)
        READ(18, *, ERR=113) ratedGenPow(n)
        READ(18, *, ERR=113) genIner(n)
        READ(18, *, ERR=113) hubIner(n)
        READ(18, *, ERR=113) bladeIner(n)
        READ(18, *, ERR=113) genEfficiency(n)
        READ(18, *, ERR=113) rateLimitGenTorque(n)
        READ(18, *, ERR=113) cutInGenSpeed(n)
        READ(18, *, ERR=113) region2StartGenSpeed(n)
        READ(18, *, ERR=113) region2EndGenSpeed(n)
        READ(18, *, ERR=113) initGenTorque(n)
        READ(18, *, ERR=113) cutInGenTorque(n)
        READ(18, *, ERR=113) ratedGenTorque(n)
        READ(18, *, ERR=113) maxGenTorque(n)
        READ(18, *, ERR=113) KGen(n)
        READ(18, *, ERR=113) speedFilterCornerFrequency(n)
        READ(18, *, ERR=113) pitchRate(n)
        CLOSE(18)
 113    CONTINUE
      ENDIF
      IF( .NOT. speedControlFileExists(n) ) speedCtrlOpt(n) = .FALSE.
      IF( .NOT. speedCtrlOpt(n) .AND. pitchCtrlOpt(n) ) pitchCtrlOpt(n) = .FALSE.
      CALL fileCheck( "windTurbines/"//TRIM( turbineNameVect( turbineTypeVect(n) ) )//"/pitchController.tbl", pitchControlFileExists(n) )
      IF( pitchCtrlOpt(n) .AND. pitchControlFileExists(n) ) THEN
        input4 = "windTurbines/"//TRIM( turbineNameVect( turbineTypeVect(n) ) )//"/pitchController.tbl"
        OPEN(FILE=input4, UNIT=19, FORM='FORMATTED', STATUS='OLD')
        READ(19, *, ERR=114) rateLimitBladePitch(n)
        READ(19, *, ERR=114) pitchMin(n)
        READ(19, *, ERR=114) pitchMax(n)
        READ(19, *, ERR=114) pitchK(n)
        READ(19, *, ERR=114) pitchControlKP(n)
        READ(19, *, ERR=114) pitchControlKI(n)
        READ(19, *, ERR=114) pitchControlKD(n)
        CLOSE(19)
 114    CONTINUE
      ENDIF
    ENDDO
!=======================================================================
#ifdef DM_PARALLEL
      CALL wrf_dm_bcast_real( Ly, nTurbines )
      CALL wrf_dm_bcast_real( Lx, nTurbines )
      CALL wrf_dm_bcast_real( Lz, nTurbines )
      CALL wrf_dm_bcast_integer( iTow, nTurbines )
      CALL wrf_dm_bcast_integer( jTow, nTurbines )
    DO n = 1, nTurbineTypes
      CALL wrf_dm_bcast_real( hubHeight(n), 1 )
      CALL wrf_dm_bcast_real( towerDiameter(n), 1 )
      CALL wrf_dm_bcast_real( rotorDiameter(n), 1 )
      CALL wrf_dm_bcast_real( hubDiameter(n), 1 )
      CALL wrf_dm_bcast_real( nacelleLength(n), 1 )
      CALL wrf_dm_bcast_real( nacDragCoeff(n), 1 )
      CALL wrf_dm_bcast_real( towDragCoeff(n), 1 )
      CALL wrf_dm_bcast_real( overhang(n), 1 )
      CALL wrf_dm_bcast_real( undSling(n), 1 )
      CALL wrf_dm_bcast_real( mechEfficiency(n), 1 )
      CALL wrf_dm_bcast_real( tiltAngle(n), 1 )
      CALL wrf_dm_bcast_real( preconeAngle(n), 1 )
      CALL wrf_dm_bcast_real( bladeCepsilon(n, :), 1 )
      CALL wrf_dm_bcast_real( nacelleCepsilon(n), 1 )
      CALL wrf_dm_bcast_real( towerCepsilon(n), 1 )
      CALL wrf_dm_bcast_real( inflowLocation(n), 1 )
      CALL wrf_dm_bcast_real( yawRate(n), 1 )
      CALL wrf_dm_bcast_real( thresholdYaw(n), 1 )
      CALL wrf_dm_bcast_real( cutIn(n), 1 )
      CALL wrf_dm_bcast_real( cutOut(n), 1 )
      CALL wrf_dm_bcast_integer( nBlade(n), 1 )
      IF( speedCtrlOpt(n) .AND. speedControlFileExists(n) ) THEN
        CALL wrf_dm_bcast_real( ratedRotSpeed(n), 1 )
        CALL wrf_dm_bcast_real( GBEfficiency(n), 1 )
        CALL wrf_dm_bcast_real( GBRatio(n), 1 )
        CALL wrf_dm_bcast_real( ratedGenPow(n), 1 )
        CALL wrf_dm_bcast_real( genIner(n), 1 )
        CALL wrf_dm_bcast_real( hubIner(n), 1 )
        CALL wrf_dm_bcast_real( bladeIner(n), 1 )
        CALL wrf_dm_bcast_real( genEfficiency(n), 1 )
        CALL wrf_dm_bcast_real( rateLimitGenTorque(n), 1 )
        CALL wrf_dm_bcast_real( cutInGenSpeed(n), 1 )
        CALL wrf_dm_bcast_real( region2StartGenSpeed(n), 1 )
        CALL wrf_dm_bcast_real( region2EndGenSpeed(n), 1 )
        CALL wrf_dm_bcast_real( initGenTorque(n), 1 )
        CALL wrf_dm_bcast_real( cutInGenTorque(n), 1 )
        CALL wrf_dm_bcast_real( ratedGenTorque(n), 1 )
        CALL wrf_dm_bcast_real( maxGenTorque(n), 1 )
        CALL wrf_dm_bcast_real( KGen(n), 1 )
        CALL wrf_dm_bcast_real( speedFilterCornerFrequency(n), 1 )
        CALL wrf_dm_bcast_real( pitchRate(n), 1 )
      ENDIF
      IF( pitchCtrlOpt(n) .AND. pitchControlFileExists(n) ) THEN
        CALL wrf_dm_bcast_real( rateLimitBladePitch(n), 1 )
        CALL wrf_dm_bcast_real( pitchMin(n), 1 )
        CALL wrf_dm_bcast_real( pitchMax(n), 1 )
        CALL wrf_dm_bcast_real( pitchK(n), 1 )
        CALL wrf_dm_bcast_real( pitchControlKP(n), 1 )
        CALL wrf_dm_bcast_real( pitchControlKI(n), 1 )
        CALL wrf_dm_bcast_real( pitchControlKD(n), 1 )
      ENDIF
    ENDDO
#endif
    ! Calculate section angle:
    angle = 360.0/nSections
    DO n = 1, nTurbineTypes
      ! Calculate solidity:
      solidity(n) = REAL( nBlade(n) )/REAL( nSections )
      bladeEpsilon(n, 1) = bladeCepsilon(n, 1)*dx
      bladeEpsilon(n, 2) = bladeCepsilon(n, 2)*dx
      bladeEpsilon(n, 3) = bladeCepsilon(n, 3)*dx
      nacelleEpsilon(n) = nacelleCepsilon(n)*dx
      towerEpsilon(n) = towerCepsilon(n)*dx
      ! First compute the radius of the force projection (to the radius
      ! where the projection is only 0.001 its maximum value - this seems
      ! recover 99.9% of the total forces when integrated).    
      projectionRadius(n) = SQRT( LOG(1.0/0.001) )*MAXVAL( bladeEpsilon(n, :) )
      !-- The zero point (epsilonMin) of the polynomial function must be the following if the integral of
      !-- the polynomial function (for values < epsilonMin) shall be equal to the integral of the Gaussian
      !-- function used before:
      epsilonMin(n) = ( (105.0/32.0)**(1.0/3.0) )*( piconst**(1.0/6.0) )*MINVAL( bladeEpsilon(n, :) )
      !-- Square of epsilonMin:
      epsilonMin2(n) = epsilonMin(n)**2
      !-- Parameters in the polynomial function:
      polA(n) = 1.0/(epsilonMin(n)**4)
      polB(n) = 2.0/(epsilonMin(n)**2)
      !-- Normalization factor which is the inverse of the integral of the smearing function:
      epsilonFactor(n) = 105.0/( 32.0*piconst*(epsilonMin(n)**3) )
      ! Check this later for better interpolation.
      xLocRot(n) = CEILING(rTip(n)/dx)
      yLocRot(n) = CEILING(rTip(n)/dy)
      xLocHub(n) = CEILING(rHub(n)/dx)
      yLocHub(n) = CEILING(rHub(n)/dy)
      xLocTow(n) = CEILING(rTow(n)/dx)
      yLocTow(n) = CEILING(rTow(n)/dy)
    ENDDO
!
    ALLOCATE( specTrb(1:nPlnWT) )
!
    ALLOCATE( bladeRadius(1:nTurbineTypes, 1:nElements) )
    ALLOCATE( c(1:nTurbineTypes, 1:nElements) )
    ALLOCATE( t(1:nTurbineTypes, 1:nElements) )
    ALLOCATE( thetaS(1:nTurbineTypes, 1:nElements) )
    ALLOCATE( foilType(1:nTurbineTypes, 1:nElements) )
    ALLOCATE( foilTypeREAL(1:nTurbineTypes, 1:nElements) )
    ALLOCATE( foilTypeVectREAL(1:nTurbineTypes, 1:nElements) )
!
    ALLOCATE( dr(1:nTurbineTypes) )
    ALLOCATE( towerDs(1:nTurbineTypes) )
    ALLOCATE( nacelleDs(1:nTurbineTypes) )
    ALLOCATE( contribution(1:nTurbineTypes) )
    ALLOCATE( zmaxtrb(1:nTurbineTypes) )
    ALLOCATE( kTrbBot(1:nTurbineTypes) )
    ALLOCATE( kTrbHub(1:nTurbineTypes) )
    ALLOCATE( kTrbTop(1:nTurbineTypes) )
    ALLOCATE( kHubBot(1:nTurbineTypes) )
    ALLOCATE( kHubTop(1:nTurbineTypes) )
    ALLOCATE( kTowBot(1:nTurbineTypes) )
    ALLOCATE( kTowTop(1:nTurbineTypes) )
    ALLOCATE( dzconst(1:nTurbineTypes) )
    ALLOCATE( hubLower(1:nTurbineTypes) )
    ALLOCATE( hubUpper(1:nTurbineTypes) )
    ALLOCATE( bladeUpper(1:nTurbineTypes) )
    ALLOCATE( bladeLower(1:nTurbineTypes) )
    ALLOCATE( speedErrorInit(1:nTurbineTypes) )
    ALLOCATE( trqGenInit(1:nTurbineTypes) )
    ALLOCATE( rotSpeedInit(1:nTurbineTypes) )
    ALLOCATE( pitchInit(1:nTurbineTypes) )
    ALLOCATE( driveTrainIner(1:nTurbineTypes) )
    ALLOCATE( towerToRotorApex(1:nTurbineTypes) )
    ALLOCATE( undSlingRotor(1:nTurbineTypes) )
    ALLOCATE( delta(1:nTurbineTypes) )
    ALLOCATE( beta(1:nTurbineTypes) )
    ALLOCATE( zeta(1:nTurbineTypes) )
    ALLOCATE( uvShaftDir(1:nTurbineTypes) )
!
    ALLOCATE( doYaw(1:nTurbines) )
    ALLOCATE( doYawLocal(1:nTurbines) )
    ALLOCATE( reOperate(1:nTurbines) )
    ALLOCATE( initialYaw(1:nTurbines) )
    ALLOCATE( initialYawLocal(1:nTurbines) )
    ALLOCATE( initialYawDone(1:nTurbines) )
    ALLOCATE( GK(1:nTurbines) )
    ALLOCATE( speedErrorLast(1:nTurbines) )
    ALLOCATE( speedErrorLocal(1:nTurbines) )
    ALLOCATE( intspeedErrorLocal(1:nTurbines) )
    ALLOCATE( derivSpeedError(1:nTurbines) )
    ALLOCATE( pitchP(1:nTurbines) )
    ALLOCATE( pitchI(1:nTurbines) )
    ALLOCATE( pitchD(1:nTurbines) )
    ALLOCATE( pitchCommanded(1:nTurbines) )
    ALLOCATE( pitchRateCommanded(1:nTurbines) )
    ALLOCATE( rotSpeedLocal(1:nTurbines) )
    ALLOCATE( rotSpeedF(1:nTurbines) )
    ALLOCATE( windDirFstAveLocal(1:nTurbines) )
    ALLOCATE( uFstAveLocal(1:nTurbines) )
    ALLOCATE( vFstAveLocal(1:nTurbines) )
    ALLOCATE( wFstAveLocal(1:nTurbines) )
    ALLOCATE( deltaAzimuthOld(1:nTurbines) )
    ALLOCATE( deltaAzimuthLocal(1:nTurbines) )
    ALLOCATE( psiOld(1:nTurbines) )
    ALLOCATE( psiOldLocal(1:nTurbines) )
    ALLOCATE( azimuthOld(1:nTurbines) )
    ALLOCATE( azimuthLocal(1:nTurbines) )
!
    ALLOCATE( rhoFstTSLocal(1:nTurbines, 1:nAveSamples) )
    ALLOCATE( uFstTSLocal(1:nTurbines, 1:nAveSamples) )
    ALLOCATE( vFstTSLocal(1:nTurbines, 1:nAveSamples) )
    ALLOCATE( wFstTSLocal(1:nTurbines, 1:nAveSamples) )
    ALLOCATE( windDirFstTSLocal(1:nTurbines, 1:nAveSamples) )
!
#ifdef DM_PARALLEL
    ALLOCATE( sender(1:nTurbines) )
    ALLOCATE( procID(1:nTurbines) )
    ALLOCATE( procInit(1:nTurbines) )
    ALLOCATE( procIDInit(1:nTurbines) )
    ALLOCATE( procIDFstInit(1:nTurbines) )
    ALLOCATE( ProcFstInit(1:nTurbines) )
    ALLOCATE( procIDFst(1:nTurbines) )
    ALLOCATE( proc(1:nTurbines) )
    ALLOCATE( procFst(1:nTurbines) )
    ALLOCATE( newGroup(1:nTurbines) )
    ALLOCATE( newComm(1:nTurbines) )
#endif
    ALLOCATE( iRot(1:nTurbines) )
    ALLOCATE( jRot(1:nTurbines) )
    ALLOCATE( iFstLocal(1:nTurbines) )
    ALLOCATE( jFstLocal(1:nTurbines) )
    ALLOCATE( iFstInit(1:nTurbines) )
    ALLOCATE( jFstInit(1:nTurbines) )
    ALLOCATE( kFst(1:nTurbines) )
    ALLOCATE( windDir(1:nTurbines) )
    ALLOCATE( xMinTile(1:nTurbines) )
    ALLOCATE( xMaxTile(1:nTurbines) )
    ALLOCATE( yMinTile(1:nTurbines) )
    ALLOCATE( yMaxTile(1:nTurbines) )
    ALLOCATE( xMinTrb(1:nTurbines) )
    ALLOCATE( xMaxTrb(1:nTurbines) )
    ALLOCATE( yMinTrb(1:nTurbines) )
    ALLOCATE( yMaxTrb(1:nTurbines) )
    ALLOCATE( deltaYawLocal(1:nTurbines) )
    ALLOCATE( deltaYaw(1:nTurbines) )
    ALLOCATE( rotorApexInit(1:nTurbines, 1:3) )
    ALLOCATE( rotorApexLocal(1:nTurbines, 1:3) )
    ALLOCATE( posRotorApex(1:nTurbines, 1:3) )
    ALLOCATE( towerBaseLocation(1:nTurbines, 1:3) )
    ALLOCATE( uvShaft(1:nTurbines, 1:3) )
    ALLOCATE( uvTower(1:nTurbines, 1:3) )
    ALLOCATE( towerShaftIntersect(1:nTurbines, 1:3) )
    ALLOCATE( bladeRoot(1:nTurbines, 1:3) )
    ALLOCATE( capPhi_init(1:nTurbines) )
    ALLOCATE( bladePoints(1:nTurbines, 1:nSections, 1:nElements, 1:3) )
    ALLOCATE( bladePointsLocal(1:nTurbines, 1:nSections, 1:nElements, 1:3) )
    ALLOCATE( towerPoints(1:nTurbines, 1:nTowPoints, 1:3) )
    ALLOCATE( towerPointsLocal(1:nTurbines, 1:nTowPoints, 1:3) )
    ALLOCATE( nacellePoints(1:nTurbines, 1:nNacPoints, 1:3) )
    ALLOCATE( nacellePointsLocal(1:nTurbines, 1:nNacPoints, 1:3) )
!
    ALLOCATE( planePoints(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints, 1:3) )
    ALLOCATE( planePointsLocal(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints, 1:3) )
!
    ALLOCATE( rhoFstAveLocal(1:nTurbines) )
    ALLOCATE( trbDirLocal(1:nTurbines) )
    ALLOCATE( wsRotorLocal(1:nTurbines) )
    ALLOCATE( v0FstLocal(1:nTurbines) )
    ALLOCATE( capPhiLocal(1:nTurbines) )
    ALLOCATE( eYawLocal(1:nTurbines) )
    ALLOCATE( deltaCapPhiLocal(1:nTurbines) )
    ALLOCATE( tsrSavedLocal(1:nTurbines) )
    ALLOCATE( thrSavedLocal(1:nTurbines) )
    ALLOCATE( trqSavedLocal(1:nTurbines) )
    ALLOCATE( powSavedLocal(1:nTurbines) )
    ALLOCATE( mechPowSavedLocal(1:nTurbines) )
    ALLOCATE( genPowSavedLocal(1:nTurbines) )
    ALLOCATE( thrCoefSavedLocal(1:nTurbines) )
    ALLOCATE( trqCoefSavedLocal(1:nTurbines) )
    ALLOCATE( powCoefSavedLocal(1:nTurbines) )
    ALLOCATE( trqGenCommandedOldOld(1:nTurbines) )
    ALLOCATE( trqGenCommandedOldLocal(1:nTurbines) )
    ALLOCATE( genSpeedOldLocal(1:nTurbines) )
    ALLOCATE( genSpeedFOldOld(1:nTurbines) )
    ALLOCATE( genSpeedFOldLocal(1:nTurbines) )
    ALLOCATE( pitchOld(1:nTurbines) )
    ALLOCATE( pitchLocal(1:nTurbines) )
    ALLOCATE( pitchSavedLocal(1:nTurbines) )
!
    ALLOCATE( fnSavedLocal(1:nTurbines, 1:nSections, 1:nElements) )
    ALLOCATE( ftSavedLocal(1:nTurbines, 1:nSections, 1:nElements) )
    ALLOCATE( FSavedLocal(1:nTurbines, 1:nSections, 1:nElements) )
    ALLOCATE( phiSavedLocal(1:nTurbines, 1:nSections, 1:nElements) )
    ALLOCATE( alphaSavedLocal(1:nTurbines, 1:nSections, 1:nElements) )
    ALLOCATE( clSavedLocal(1:nTurbines, 1:nSections, 1:nElements) )
    ALLOCATE( cdSavedLocal(1:nTurbines, 1:nSections, 1:nElements) )
    ALLOCATE( cmSavedLocal(1:nTurbines, 1:nSections, 1:nElements) )
    ALLOCATE( LSavedLocal(1:nTurbines, 1:nSections, 1:nElements) )
    ALLOCATE( DSavedLocal(1:nTurbines, 1:nSections, 1:nElements) )
    ALLOCATE( v1SavedLocal(1:nTurbines, 1:nSections, 1:nElements) )
    ALLOCATE( vTSavedLocal(1:nTurbines, 1:nSections, 1:nElements) )
    ALLOCATE( vRelSavedLocal(1:nTurbines, 1:nSections, 1:nElements) )
    ALLOCATE( vThetaSavedLocal(1:nTurbines, 1:nSections, 1:nElements) )
    ALLOCATE( vZSavedLocal(1:nTurbines, 1:nSections, 1:nElements) )
    ALLOCATE( uSavedLocal(1:nTurbines, 1:nSections, 1:nElements) )
    ALLOCATE( vSavedLocal(1:nTurbines, 1:nSections, 1:nElements) )
    ALLOCATE( wSavedLocal(1:nTurbines, 1:nSections, 1:nElements) )
!
    ALLOCATE( uPlaneSavedLocal(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints) )
    ALLOCATE( vPlaneSavedLocal(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints) )
    ALLOCATE( wPlaneSavedLocal(1:nPlnWT, 1:nVertPoints, 1:nHorizPoints) )
!=======================================================================
! Initialize allocated variables:
!=======================================================================
    doYaw(:) = .FALSE.
    doYawLocal(:) = .FALSE.
    initialYaw(:) = .FALSE.
    initialYawLocal(:) = .FALSE.
    initialYawDone(:) = .FALSE.
    reOperate(:) = .FALSE.
#ifdef DM_PARALLEL
    sender(:) = 0
    procID(:) = 0
    procIDFst(:) = 0
    procIDFstInit(:) = 0
    procInit(:) = 0
    procIDInit(:) = 0
    proc(:) = 0
    procFst(:) = 0
    ProcFstInit(:) = 0
    newGroup(:) = 0
    newComm(:) = 0
#endif
    iRot(:) = 0
    jRot(:) = 0
    iFstLocal(:) = 0
    jFstLocal(:) = 0
    iFstInit(:) = 0
    jFstInit(:) = 0
    kFst(:) = 0
    dr(:) = 0.0
    towerDs(:) = 0.0
    nacelleDs(:) = 0.0
    contribution(:) = 0.0
    zmaxtrb(:) = 0.0
    kTrbBot(:) = 0.0
    kTrbHub(:) = 0.0
    kTrbTop(:) = 0.0
    kHubBot(:) = 0.0
    kHubTop(:) = 0.0
    kTowBot(:) = 0.0
    kTowTop(:) = 0.0
    dzconst(:) = 0.0
    speedErrorInit(:) = 0.0
    trqGenInit(:) = 0.0
    rotSpeedInit(:) = 0.0
    pitchInit(:) = 0.0
    speedErrorLast(:) = 0.0
    derivSpeedError(:) = 0.0
    pitchP(:) = 0.0
    pitchI(:) = 0.0
    pitchD(:) = 0.0
    pitchCommanded(:) = 0.0
    pitchRateCommanded(:) = 0.0
    trqGenCommandedOldOld(:) = 0.0
    genSpeedFOldOld(:) = 0.0
    pitchOld(:) = 0.0
    azimuthOld(:) = 90.0 ! Assuming the first actuator line is aligned with z-axis (standard axis notation).
                         ! Only used for diagnostics
    psiOld(:) = 0.0 ! Assuming the first actuator line is aligned with z-axis
    deltaAzimuthOld(:) = 0.0 
    rotSpeedF(:) = 0.0
    driveTrainIner(:) = 0.0
    windDir(:) = 0.0
    delta(:) = 0.0
    capPhi_init(:) = 0.0
    deltaYaw(:) = 0.0
    deltaYawLocal(:) = 0.0
    towerToRotorApex(:) = 0.0
    undSlingRotor(:) = 0.0
    xMinTile(:) = 0.0
    xMaxTile(:) = 0.0
    yMinTile(:) = 0.0
    yMaxTile(:) = 0.0
    xMinTrb(:) = 0.0
    xMaxTrb(:) = 0.0
    yMinTrb(:) = 0.0
    yMaxTrb(:) = 0.0
    beta(:) = 0.0
    zeta(:) = 0.0
    wsRotorLocal(:) = 0.0
    hubLower(:) = 0.0
    hubUpper(:) = 0.0
    bladeLower(:) = 0.0
    bladeUpper(:) = 0.0
    uvShaftDir(:) = 0.0
    bladeRoot(:, :) = 0.0
    rotorApexInit(:, :) = 0.0
    posRotorApex(:, :) = 0.0
    towerBaseLocation(:, :) = 0.0
    uvShaft(:, :) = 0.0
    uvTower(:, :) = 0.0
    towerShaftIntersect(:, :) = 0.0
    foilType(:, :) = -1
    bladeRadius(:, :) = 0.0
    c(:, :) = 0.0
    t(:, :) = 0.0
    thetaS(:, :) = 0.0
    foilTypeREAL(:, :) = -1.0
    foilTypeVectREAL(:, :) = -1.0
    towerPoints(:, :, :) = 0.0
    nacellePoints(:, :, :) = 0.0
    bladePoints(:, :, :, :) = 0.0
!
    planePoints(:, :, :, :) = 0.0
!
    specTrb(:) = -999999
!
    hubHeightASL(:) = -999999.999
    hubHeightASLLocal(:) = -999999.999
    ! History & Restart
    rhoFstAveLocal(:) = -999999.999
    capPhiLocal(:) = -999999.999
    deltaCapPhiLocal(:) = -999999.999
    eYawLocal(:) = -999999.999
    uFstAveLocal(:) = -999999.999
    vFstAveLocal(:) = -999999.999
    wFstAveLocal(:) = -999999.999
    windDirFstAveLocal(:) = -999999.999
    tsrSavedLocal(:) = -999999.999
    thrSavedLocal(:) = -999999.999
    trqSavedLocal(:) = -999999.999
    powSavedLocal(:) = -999999.999
    mechPowSavedLocal(:) = -999999.999
    genPowSavedLocal(:) = -999999.999
    thrCoefSavedLocal(:) = -999999.999
    trqCoefSavedLocal(:) = -999999.999
    powCoefSavedLocal(:) = -999999.999
    trbDirLocal(:) = -999999.999
    rotSpeedLocal(:) = -999999.999
    v0FstLocal(:) = -999999.999
    pitchLocal(:) = -999999.999
    pitchSavedLocal(:) = -999999.999
    trqGenCommandedOldLocal(:) = -999999.999
    genSpeedOldLocal(:) = -999999.999
    genSpeedFOldLocal(:) = -999999.999
    speedErrorLocal(:) = -999999.999
    intspeedErrorLocal(:) = -999999.999
    rhoFstTSLocal(:, :) = -999999.999
    psiOldLocal(:) = -999999.999
    deltaAzimuthLocal(:) = -999999.999
    azimuthLocal(:) = -999999.999
    rhoFstTSLocal(:, :) = -999999.999
    uFstTSLocal(:, :) = -999999.999
    vFstTSLocal(:, :) = -999999.999
    wFstTSLocal(:, :) = -999999.999
    windDirFstTSLocal(:, :) = -999999.999
    rotorApexLocal(:, :) = -999999.999
    fnSavedLocal(:, :, :) = -999999.999
    ftSavedLocal(:, :, :) = -999999.999
    FSavedLocal(:, :, :) = -999999.999
    phiSavedLocal(:, :, :) = -999999.999
    alphaSavedLocal(:, :, :) = -999999.999
    clSavedLocal(:, :, :) = -999999.999
    cdSavedLocal(:, :, :) = -999999.999
    cmSavedLocal(:, :, :) = -999999.999
    LSavedLocal(:, :, :) = -999999.999
    DSavedLocal(:, :, :) = -999999.999
    v1SavedLocal(:, :, :) = -999999.999
    vThetaSavedLocal(:, :, :) = -999999.999
    vTSavedLocal(:, :, :) = -999999.999
    vRelSavedLocal(:, :, :) = -999999.999
    vZSavedLocal(:, :, :) = -999999.999
    uSavedLocal(:, :, :) = -999999.999
    vSavedLocal(:, :, :) = -999999.999
    wSavedLocal(:, :, :) = -999999.999
    towerPointsLocal(:, :, :) = -999999.999
    nacellePointsLocal(:, :, :) = -999999.999
    bladePointsLocal(:, :, :, :) = -999999.999
!
    uPlaneSavedLocal(:, :, :) = -999999.999
    vPlaneSavedLocal(:, :, :) = -999999.999
    wPlaneSavedLocal(:, :, :) = -999999.999
    planePointsLocal(:, :, :, :) = -999999.999
!=======================================================================
    foilTypeVectREAL = foilTypeVect
    ! Uniform point distribution:
    DO n = 1, nTurbineTypes
      dr(n) = ( rTip(n)-rHub(n) )/nElements
      dist = 0.0
      WRITE(0,*) TRIM( turbineNameVect( turbineTypeVect(n) ) )
      WRITE(0,*) 'blade-element # ',' r/R     [-]',' chord     [m]',' thickness     [%]', ' twist     [deg]', '      airfoil #'
      DO e = 1, nElements
        dist = dist + 0.5*dr(n)
        bladeRadius(n, e) = rHub(n) + dist
        dist = dist + 0.5*dr(n)
        CALL interpolate( trbDataVectLen(n), bladeRadius(n, e), rVect(n, :), chordVect(n, :), c(n, e) )
        CALL interpolate( trbDataVectLen(n), bladeRadius(n, e), rVect(n, :), thicknessVect(n, :), t(n, e) )
        CALL interpolate( trbDataVectLen(n), bladeRadius(n, e), rVect(n, :), thetaGVect(n, :), thetaS(n, e) )
        CALL interpolate( trbDataVectLen(n), bladeRadius(n, e), rVect(n, :), foilTypeVectREAL(n, :), foilTypeREAL(n, e) )
        foilType(n, e) = MIN( MAX( 1, NINT( foilTypeREAL(n, e) ) ), nFoilTypes(n) )
        WRITE(0,*) e, bladeRadius(n, e)/(rotorDiameter(n)/2), c(n, e), t(n, e), thetaS(n, e), foilType(n, e)
      ENDDO
    ENDDO
#ifdef DM_PARALLEL
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, commSize, ierr)
#endif
    itf = MIN0(ite, ide-1)
    jtf = MIN0(jte, jde-1)
    ktf = MIN0(kte, kde-1)
    DO j = jts-1, jtf+1   
      DO i = its-1, itf+1
         w_z(i, kts, j) = 0.0     ! meter above ground, add hgt to get height ASL
      ENDDO
    ENDDO
    DO k = kts+1, ktf+1           ! compute 1D column of z-values to determine kEnd for the tile
      w_z(its, k, jts) = w_z(its, k-1, jts) + dz8w(its, k-1, jts)
    ENDDO
    DO j = jts, jtf               ! compute the rest of w_z up to kEnd
      DO k = kts+1, ktf+1
        DO i = its, itf
          w_z(i, k, j) = w_z(i, k-1, j) + dz8w(i, k-1, j)
        ENDDO
      ENDDO
    ENDDO
    DO j = jts, jtf               ! compute the i_start and i_end ghost nodes
      DO k = kts+1, ktf+1
        w_z(its-1 ,k, j) = w_z(its, k, j)
        w_z(itf+1, k, j) = w_z(itf, k, j)
      ENDDO
    ENDDO
    DO k = kts+1, ktf+1           ! compute the j_start and j_end ghost nodes
      DO i = its, itf  
        w_z(i, k, jts-1) = w_z(i, k, jts)
        w_z(i, k, jtf+1) = w_z(i, k, jtf)
      ENDDO
    ENDDO
    DO k = kts+1, ktf+1           ! compute get the corner ghost nodes
      w_z(its-1, k, jts-1) = 0.5*(w_z(its-1, k, jts) + w_z(its, k, jts-1) )
      w_z(its-1, k, jtf+1) = 0.5*(w_z(its-1, k, jtf) + w_z(its, k, jtf+1) )
      w_z(itf+1, k, jts-1) = 0.5*(w_z(itf+1, k, jts) + w_z(itf, k, jts-1) ) 
      w_z(itf+1, k, jtf+1) = 0.5*(w_z(itf+1, k, jtf) + w_z(itf, k, jtf+1) ) 
    ENDDO
    DO i = its-1, itf+1
      DO k = kts, ktf
        DO j = jts-1, jtf+1
          dz(i, k, j) = w_z(i, k+1, j) - w_z(i, k, j)
        ENDDO
      ENDDO
    ENDDO
    specTrb = model_config_rec%wind_wtp_pln_turbine_ids(1:nPlnWT)
    nPlnWTcounter = 1
    DO kt = 1, nTurbines
      wtID = FINDLOC(turbineTypeVect, VALUE=turbineTypeID(kt), DIM=1)
      zMaxTrb(wtID) = hubHeight(wtID) + (rotorDiameter(wtID)/2) + smoothing*dx
      DO k = kts+1, ktf+1           ! kEnd = 1 index greater to allow for nonuniform heights within the tile
        IF( (w_z(its, k-1, jts) .LT. zMaxTrb(wtID) ) .AND. (w_z(its, k, jts) .GE. zMaxTrb(wtID) ) ) kEnd(wtID) = k+1
      ENDDO 
      reOperate(kt) = .FALSE.
      IF( pitchCtrlOpt(wtID) ) THEN
        initialYawDone(kt) = .FALSE.
        pitchInit(wtID) = 0.0
      ENDIF
      IF( speedCtrlOpt(wtID) ) THEN
        initialYawDone(kt) = .FALSE.
        speedErrorInit(wtID) = 0.0
        rotSpeedInit(wtID) = 0.0
        trqGenInit(wtID) = initGenTorque(wtID)
        driveTrainIner(wtID) = (nBlade(wtID)*bladeIner(wtID) + hubIner(wtID) + GBRatio(wtID)*GBRatio(wtID)*genIner(wtID) )
      ENDIF
      bladeLower(wtID) = hubHeight(wtID) - rTip(wtID) ! Height of lower blade tip above ground (m)
      bladeUpper(wtID) = hubHeight(wtID) + rTip(wtID) ! Height of upper blade tip above ground (m)
      hubLower(wtID) = hubHeight(wtID) - rHub(wtID) ! Height of lower hub tip above ground (m)
      hubUpper(wtID) = hubHeight(wtID) + rHub(wtID) ! Height of upper hub tip above ground (m)
      xMinTrb(kt) = Lx(kt) - (rotorDiameter(wtID)/2) - smoothing*dx
      yMinTrb(kt) = Ly(kt) - (rotorDiameter(wtID)/2) - smoothing*dy
      xMaxTrb(kt) = Lx(kt) + (rotorDiameter(wtID)/2) + smoothing*dx
      yMaxTrb(kt) = Ly(kt) + (rotorDiameter(wtID)/2) + smoothing*dy 
      towerBaseLocation(kt, 1) = Lx(kt)
      towerBaseLocation(kt, 2) = Ly(kt)
      towerBaseLocation(kt, 3) = Lz(kt)
      towerShaftIntersect(kt, 1) = towerBaseLocation(kt, 1)
      towerShaftIntersect(kt, 2) = towerBaseLocation(kt, 2)
      towerShaftIntersect(kt, 3) = towerBaseLocation(kt, 3) + hubHeight(wtID)
      ! Define vector aligned with the tower pointing from the ground to the nacelle
      uvTower(kt, :) = vectorAdd( towerShaftIntersect(kt, :), -towerBaseLocation(kt, :) )
      uvTower(kt, :) = vectorDivide( uvTower(kt, :), vectorMag( uvTower(kt, :) ) )
      IF( ( its .LE. iTow(kt) .AND. iTow(kt) .LE. itf ) .AND. ( jts .LE. jTow(kt) .AND. jTow(kt) .LE. jtf ) ) THEN
        kTrbHub(wtID) = MINLOC( ABS( w_z(iTow(kt), kts:ktf, jTow(kt) ) - hubHeight(wtID) ), 1 )
        kHubBot(wtID) = MINLOC( ABS( w_z(iTow(kt), kts:kEnd(wtID), jTow(kt) ) - hubLower(wtID) ), 1 )
        kHubTop(wtID) = MINLOC( ABS( w_z(iTow(kt), kts:kEnd(wtID), jTow(kt) ) - hubUpper(wtID) ), 1 )
        kTowBot(wtID) = MINLOC( ABS( w_z(iTow(kt), kts:kEnd(wtID), jTow(kt) ) - towerBaseLocation(kt, 3) ), 1 )
        kTowTop(wtID) = MINLOC( ABS( w_z(iTow(kt), kts:kEnd(wtID), jTow(kt) ) - hubLower(wtID) ), 1 )
        IF( w_z( iTow(kt), kHubBot(wtID), jTow(kt) ) .GT. hubLower(wtID) ) kHubBot(wtID) = kHubBot(wtID) - 1
        IF( w_z( iTow(kt), kHubTop(wtID), jTow(kt) ) .GT. hubUpper(wtID) ) kHubTop(wtID) = kHubTop(wtID) - 1
        IF( w_z( iTow(kt), kTowBot(wtID), jTow(kt) ) .GT. towerBaseLocation(kt, 3) ) kTowBot(wtID) = kTowBot(wtID) - 1
        IF( w_z( iTow(kt), kTowTop(wtID), jTow(kt) ) .GT. hubLower(wtID) ) kTowTop(wtID) = kTowTop(wtID) - 1
        IF( w_z( iTow(kt), kTrbHub(wtID), jTow(kt) ) .GT. hubHeight(wtID) ) kTrbHub(wtID) = kTrbHub(wtID) - 1
        IF( kTowBot(wtID) .EQ. 0 ) kTowBot(wtID) = kts
        windDir(kt) = ATAN2( u_2( iTow(kt), kTrbHub(wtID), jTow(kt) ), v_2( iTow(kt), kTrbHub(wtID), jTow(kt) ) )/degrad + 180.0
#ifdef DM_PARALLEL
        CALL wrf_get_myproc( procInit(kt) )
#endif
      ENDIF
#ifdef DM_PARALLEL
      procIDInit(kt) = wrf_dm_max_int( procInit(kt) )
      CALL MPI_BCAST(windDir(kt), 1, MPI_REAL, procIDInit(kt), MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(kHubBot(wtID), 1, MPI_INT, procIDInit(kt), MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(kHubTop(wtID), 1, MPI_INT, procIDInit(kt), MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(kTowBot(wtID), 1, MPI_INT, procIDInit(kt), MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(kTowTop(wtID), 1, MPI_INT, procIDInit(kt), MPI_COMM_WORLD, ierr)
#endif
      IF( specDir .GE. 0.0 ) windDir(kt) = specDir
      CALL compassToStandard( windDir(kt) )
      uvShaftDir(wtID) = -1.0
      delta(wtID) = uvShaftDir(wtID)*tiltAngle(wtID)
      beta(wtID) = uvShaftDir(wtID)*preconeAngle(wtID)
      zeta(wtID) = tiltAngle(wtID)-preconeAngle(wtID)
      towerDs(wtID) = ( hubHeight(wtID)-rHub(wtID) )/nTowPoints
      nacelleDs(wtID) = nacelleLength(wtID)/nNacPoints
      contribution(wtID) = nacelleDs(wtID)/nacelleLength(wtID)
      IF( .NOT. restart .OR. (rotorApexX(kt) .EQ. 0.0 .AND. rotorApexY(kt) .EQ. 0.0 .AND. rotorApexZ(kt) .EQ. 0.0) ) THEN
        capPhi_init(kt) = windDir(kt)
        towerToRotorApex(wtID) = uvShaftDir(wtID)*overhang(wtID)
        undSlingRotor(wtID) = uvShaftDir(wtID)*undSling(wtID)
        rotorApexInit(kt, 1) = towerBaseLocation(kt, 1) + ( towerToRotorApex(wtID) + undSlingRotor(wtID) )*COS(delta(wtID)*degrad)
        rotorApexInit(kt, 2) = towerBaseLocation(kt, 2)
        rotorApexInit(kt, 3) = towerBaseLocation(kt, 3) + hubHeight(wtID) + ( towerToRotorApex(wtID) + undSlingRotor(wtID) )*SIN(delta(wtID)*degrad)
        bladeRoot(kt, 1) = rotorApexInit(kt, 1) + rHub(wtID)*SIN(zeta(wtID)*degrad)
        bladeRoot(kt, 2) = rotorApexInit(kt, 2)
        bladeRoot(kt, 3) = rotorApexInit(kt, 3) + rHub(wtID)*COS(zeta(wtID)*degrad)
        ! Define the vector along the shaft pointing in the direction of the wind
        uvShaft(kt, :) = vectorAdd( rotorApexInit(kt, :), -towerShaftIntersect(kt, :) )
        uvShaft(kt, :) = vectorDivide( uvShaft(kt, :), vectorMag( uvShaft(kt, :) ) )
        uvShaft(kt, :) = vectorMultiply( uvShaft(kt, :), uvShaftDir(wtID) )
        ! Loop over blade sectors:
        DO n = 1, nSections
          dist = 0.0
          ! Loop over blade radial sections:
          DO e = 1, nElements
            dist = dist + 0.5*dr(wtID)
            bladePoints(kt, n, e, 1) = bladeRoot(kt, 1) + dist*SIN(zeta(wtID)*degrad)
            bladePoints(kt, n, e, 2) = bladeRoot(kt, 2)
            bladePoints(kt, n, e, 3) = bladeRoot(kt, 3) + dist*COS(zeta(wtID)*degrad)
            dist = dist + 0.5*dr(wtID)
            ! If there are more than one blade create the points of other blades by
            ! rotating the points of the first blade
            IF( n .GT. 1) THEN
              bladePoints(kt, n, e, :) = rotatePoint( bladePoints(kt, n, e, :), rotorApexInit(kt, :), uvShaft(kt, :), (360.0/nSections)*(n-1) )
            ENDIF
            bladePoints(kt, n, e, :) = rotatePoint( bladePoints(kt, n, e, :), towerShaftIntersect(kt, :), uvTower(kt, :), capPhi_init(kt) )
          ENDDO
        ENDDO
        IF( towParOpt ) THEN
          dist = 0.0
          ! Loop over tower points:
          DO m = 1, nTowPoints
            dist = dist + 0.5*towerDs(wtID)
            towerPoints(kt, m, 1) = towerBaseLocation(kt, 1)
            towerPoints(kt, m, 2) = towerBaseLocation(kt, 2)
            towerPoints(kt, m, 3) = towerBaseLocation(kt, 3) + dist
            dist = dist + 0.5*towerDs(wtID)
          ENDDO
        ENDIF
        IF( hubParOpt ) THEN
          dist = 0.0
          ! Loop over nacelle points:
          DO m = 1, nNacPoints
            dist = dist + 0.5*nacelleDs(wtID)
            nacellePoints(kt, m, 1) = rotorApexInit(kt, 1) + dist*COS(delta(wtID)*degrad)
            nacellePoints(kt, m, 2) = rotorApexInit(kt, 2)
            nacellePoints(kt, m, 3) = rotorApexInit(kt, 3) + dist*SIN(delta(wtID)*degrad)
            dist = dist + 0.5*nacelleDs(wtID)
            nacellePoints(kt, m, :) = rotatePoint( nacellePoints(kt, m, :), towerShaftIntersect(kt, :), uvTower(kt, :), capPhi_init(kt) )
          ENDDO
        ENDIF
        ! Update rotor apex position:
        rotorApexInit(kt, :) = rotatePoint( rotorApexInit(kt, :), towerShaftIntersect(kt, :), uvTower(kt, :), capPhi_init(kt) )
        ! Update rotor shaft vector:
        uvShaft(kt, :) = vectorAdd( rotorApexInit(kt, :), -towerShaftIntersect(kt, :) )
        uvShaft(kt, :) = vectorDivide( uvShaft(kt, :), vectorMag( uvShaft (kt, :) ) )
        uvShaft(kt, :) = vectorMultiply( uvShaft(kt, :), uvShaftDir(wtID) )
        iRot(kt) = INT( (rotorApexInit(kt, 1) + 0.5*dx)/dx)
        jRot(kt) = INT( (rotorApexInit(kt, 2) + 0.5*dy)/dy)
        rotorApexLocal(kt, :) = rotorApexInit (kt, :)
        posRotorApex(kt, :) = rotorApexLocal(kt, :)
        IF( tsoutPlnOpt .AND. ANY(specTrb .EQ. kt) ) THEN
          towerToPlane = uvShaftDir(wtID)*planeDistance
          CALL init_wind_plane_tsout_locs(                  &
              &  nPlnWT, nPlnWTcounter, nVertPoints         &
              &, nHorizPoints, dx, dy, planeWidth           &
              &, planeHeight, towerToPlane                  &
              &, towerShaftIntersect(kt, :), uvTower(kt, :) &
              &, capPhi_init(kt), Lx(kt), Ly(kt), Lz(kt)    &
              &, planePoints                                &
              &)
          nPlnWTcounter = nPlnWTcounter + 1
        ENDIF
      ELSE
        rotorApexInit(kt, 1) = rotorApexX(kt)
        rotorApexInit(kt, 2) = rotorApexY(kt)
        rotorApexInit(kt, 3) = rotorApexZ(kt)
        iRot(kt) = INT( (rotorApexX(kt) + 0.5*dx)/dx)
        jRot(kt) = INT( (rotorApexY(kt) + 0.5*dy)/dy)
        rotorApexLocal(kt, :) = rotorApexInit(kt, :)
        posRotorApex(kt, 1) = rotorApexX(kt)
        posRotorApex(kt, 2) = rotorApexY(kt)
        posRotorApex(kt, 3) = rotorApexZ(kt)
        ! Define the vector along the shaft pointing in the direction of the wind
        uvShaft(kt, :) = vectorAdd( rotorApexLocal(kt, :), -towerShaftIntersect(kt, :) )
        uvShaft(kt, :) = vectorDivide( uvShaft(kt, :), vectorMag( uvShaft(kt, :) ) )
        uvShaft(kt, :) = vectorMultiply( uvShaft(kt, :), uvShaftDir(wtID) )
        ! Loop over blade radial sections:
        DO e = 1, nElements
          ! Loop over blade sectors:
          DO n = 1, nSections
            bladePoints(kt, n, e, 1) = bladePointsX(kt, n, e)
            bladePoints(kt, n, e, 2) = bladePointsY(kt, n, e)
            bladePoints(kt, n, e, 3) = bladePointsZ(kt, n, e)
          ENDDO
        ENDDO
        IF( towParOpt ) THEN
          ! Loop over tower points:
          DO m = 1, nTowPoints
            towerPoints(kt, m, 1) = towerPointsX(kt, m)
            towerPoints(kt, m, 2) = towerPointsY(kt, m)
            towerPoints(kt, m, 3) = towerPointsZ(kt, m)
          ENDDO
        ENDIF
        IF( hubParOpt ) THEN
          ! Loop over nacelle points:
          DO m = 1, nNacPoints
            nacellePoints(kt, m, 1) = nacellePointsX(kt, m)
            nacellePoints(kt, m, 2) = nacellePointsY(kt, m)
            nacellePoints(kt, m, 3) = nacellePointsZ(kt, m)
          ENDDO
        ENDIF
        IF( tsoutPlnOpt .AND. ANY(specTrb .EQ. kt) ) THEN
          IF( (planePointsX(nPlnWTcounter, nVertPoints, nHorizPoints) .EQ. 0.0) .AND. &
              (planePointsY(nPlnWTcounter, nVertPoints, nHorizPoints) .EQ. 0.0) .AND. &
              (planePointsZ(nPlnWTcounter, nVertPoints, nHorizPoints) .EQ. 0.0) ) THEN
            capPhi_init(kt) = capPhi(kt)
            towerToPlane = uvShaftDir(wtID)*planeDistance
            CALL init_wind_plane_tsout_locs(                  &
                &  nPlnWT, nPlnWTcounter, nVertPoints         &
                &, nHorizPoints, dx, dy, planeWidth           &
                &, planeHeight, towerToPlane                  &
                &, towerShaftIntersect(kt, :), uvTower(kt, :) &
                &, capPhi_init(kt), Lx(kt), Ly(kt), Lz(kt)    &
                &, planePoints                                &
                &)
            nPlnWTcounter = nPlnWTcounter + 1
          ELSE
            nt = nPlnWTcounter
            DO m = 1, nHorizPoints
              DO n = 1, nVertPoints
                planePoints(nt, n, m, 1) = planePointsX(nt, n, m)
                planePoints(nt, n, m, 2) = planePointsY(nt, n, m)
                planePoints(nt, n, m, 3) = planePointsZ(nt, n, m)
              ENDDO
            ENDDO
            nPlnWTcounter = nPlnWTcounter + 1
          ENDIF
        ENDIF
      ENDIF
      IF( ( its .LE. iRot(kt) .AND. iRot(kt) .LE. itf ) .AND. &
          ( jts .LE. jRot(kt) .AND. jRot(kt) .LE. jtf ) ) THEN
        kTrbBot(wtID) = MINLOC( ABS( w_z(iRot(kt), kts:kEnd(wtID), jRot(kt) ) - bladeLower(wtID) ), 1 )
        kTrbHub(wtID) = MINLOC( ABS( w_z(iRot(kt), kts:kEnd(wtID), jRot(kt) ) - hubHeight(wtID) ), 1 )
        kTrbTop(wtID) = MINLOC( ABS( w_z(iRot(kt), kts:kEnd(wtID), jRot(kt) ) - bladeUpper(wtID) ), 1 )
        IF( w_z( iRot(kt), kTrbBot(wtID), jRot(kt) ) .GT. bladeLower(wtID) ) kTrbBot(wtID) = kTrbBot(wtID) - 1
        IF( w_z( iRot(kt), kTrbHub(wtID), jRot(kt) ) .GT. hubHeight(wtID) ) kTrbHub(wtID) = kTrbHub(wtID) - 1
        IF( w_z( iRot(kt), kTrbTop(wtID), jRot(kt) ) .GT. bladeUpper(wtID) ) kTrbTop(wtID) = kTrbTop(wtID) - 1
        dzConst(wtID) = dz( iRot(kt), kTrbHub(wtID), jRot(kt) )
        IF( errorCheckOpt ) THEN
          ! Check division by zero for running averages:
          IF( nAveSamples .LE. 0 ) THEN
            CALL wrf_error_fatal('WTP-ERROR ***** wind_wtp_nAvesamples and wind_wtp_nSeries_tme &
                                  are not set in namelist.input properly. *****')
          ENDIF
          ! Check indices of the turbine dimensions:
          IF( ( xMinTrb(kt) .LT. (ids-1)*dx .OR. xMaxTrb(kt) .GT. ide*dx ) .OR. &
              ( yMinTrb(kt) .LT. (jds-1)*dy .OR. yMaxTrb(kt) .GT. jde*dy ) ) THEN
            CALL wrf_error_fatal('WTP-ERROR ***** Turbine dimensions are outside of the domain. &
                                  Relocate wind turbine(s) properly. *****')
          ENDIF
          ! Check non-positive indices in vertical direction:
          IF( kHubBot(wtID) .LE. 0 .OR. kTowBot(wtID) .LE. 0 .OR. &
              kTowTop(wtID) .LE. 0 .OR. kTrbBot(wtID) .LE. 0 .OR. &
              kTrbHub(wtID) .LE. 0 .OR. kHubTop(wtID) .LE. 0 .OR. &
              kTrbTop(wtID) .LE. 0 ) THEN
            CALL wrf_error_fatal('WTP-ERROR ***** Non-positive "k" in vertical direction is found. &
                                  Check WRF-LES setup prior to the WTP simulation. *****')
          ENDIF
          DO k = kTrbBot(wtID), kTrbTop(wtID)
            ! Check model grid spacing:
            IF( NINT( rotorDiameter(wtID)/dx) .LT. 2 .OR. NINT(rotorDiameter(wtID)/dy ) .LT. 2 .OR. &
                NINT( rotorDiameter(wtID)/dzConst(wtID) ) .LT. 2 ) THEN
              WRITE(outstring, *)''
              CALL wrf_message ( TRIM(outstring) )
              WRITE(outstring, "(A, I3, A, I3, A, I3, A)") &
                               'WTP-WARNING ***** kTrbBot = ', kTrbBot(wtID), &
                               ', kTrbHub = ', kTrbHub(wtID), ', kTrbTop = ',  &
                               kTrbTop(wtID), ' *****'
              CALL wrf_message ( TRIM(outstring) )
              WRITE(outstring, "(A, F7.3, A, A, F7.3, A, A, F7.3, A, A, F7.3, A, I3, A)")        &
                               'WTP-WARNING ***** rotorDiameter = ', rotorDiameter(wtID), ' m,', &
                               ' dx = ', dx, ' m,', ' dy = ', dy, ' m,', ' dz = ', dzConst(wtID), &
                               ' m at k = ', k, ' *****'
              CALL wrf_message ( TRIM(outstring) )
              WRITE(outstring, "(A, I3, A, I3, A, I3, A)") &
                               'WTP-WARNING ***** (rotordiamater/dx) = ', NINT(rotorDiameter(wtID)/dx), &
                               ' (rotorDiameter/dy) = ', NINT(rotorDiameter(wtID)/dy),                  &
                               ' and (rotorDiameter/dz) = ',                                            &
                               NINT( rotorDiameter(wtID)/dzConst(wtID) ),                               &
                               ' must be equal to or greater than 3. *****'
              CALL wrf_message ( TRIM(outstring) )
              WRITE(outstring, *)''
              CALL wrf_message ( TRIM(outstring) )
              CALL wrf_error_fatal('WTP-ERROR ***** Increase spatial resolution for correct &
                                    modeling. dx, dy, dz values should be as close to each &
                                    other as possible. *****')
            ENDIF
            ! Check model grid aspect ratios:
            IF( ABS( dx-dzConst(wtID) ) .GT. 5.0 .OR. ABS( dy-dzConst(wtID) ) .GT. 5.0 ) THEN
              IF( dx .GT. dzConst(wtID) .OR. dy .GT. dzConst(wtID) ) THEN
                CALL wrf_error_fatal('WTP-ERROR ***** Increase horizontal resolution or decrease &
                                      vertical resolution for correct modeling. *****' )
              ELSEIF( dx .LT. dzConst(wtID) .OR. dy .LT. dzConst(wtID) ) THEN
                CALL wrf_error_fatal('WTP-ERROR ***** Increase vertical resolution or decrease &
                                      horizontal resolution for correct modeling. *****' )
              ENDIF
            ENDIF
          ENDDO
        ENDIF
#ifdef DM_PARALLEL
        CALL wrf_get_myproc( proc(kt) )
#endif
      ENDIF
      iFstInit(kt) = INT( ( rotorApexInit(kt, 1) - inflowLocation(wtID)*COS(capPhi_init(kt)*degrad) + 0.5*dx )/dx )
      ! Set horizontal boundary index in x-direction to i_tile_start if iFstInit(kt)
      ! is less than its
      IF( iFstInit(kt) .LE. 0 ) iFstInit(kt) = 1
      ! Set horizontal boundary index in x-direction to i_domain_end-1 if iFstInit(kt)
      ! is greater than ide
      IF( iFstInit(kt) .GE. ide ) iFstInit(kt) = ide-1
      jFstInit(kt) = INT( ( rotorApexInit(kt, 2) - inflowLocation(wtID)*SIN(capPhi_init(kt)*degrad) + 0.5*dy )/dy )
      ! Set horizontal boundary index in y-direction to 1 if jFstInit(kt)
      ! is less than 0
      IF( jFstInit(kt) .LE. 0 ) jFstInit(kt) = 1
      ! Set horizontal boundary index in y-direction to j_domain_end-1 if jFstInit(kt)
      ! is greater than jde
      IF( jFstInit(kt) .GE. jde ) jFstInit(kt) = jde-1
      IF( ( its .LE. iFstInit(kt) .AND. iFstInit(kt) .LE. itf ) .AND. &
          ( jts .LE. jFstInit(kt) .AND. jFstInit(kt) .LE. jtf ) ) THEN
#ifdef DM_PARALLEL
        CALL wrf_get_myproc( ProcFstInit(kt) )
#endif
      ENDIF
#ifdef DM_PARALLEL
      procID(kt) = wrf_dm_max_int( proc(kt) )
      procIDFstInit(kt) = wrf_dm_max_int( ProcFstInit(kt) )
      CALL MPI_BCAST(kTrbBot(wtID), 1, MPI_INT, procID(kt), MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(kTrbHub(wtID), 1, MPI_INT, procID(kt), MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(kTrbTop(wtID), 1, MPI_INT, procID(kt), MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(dzConst(wtID), 1, MPI_REAL, procID(kt), MPI_COMM_WORLD, ierr)
#endif
#ifdef DM_PARALLEL
      IF( nprocX .GE. 3 .AND. commSize .GT. 10 ) THEN
        procsW = procID(kt) - 1
        procsE = procID(kt) + 1
        procsN = procID(kt) + nprocX
        procsS = procID(kt) - nprocX
        procsNW = procsN - 1
        procsNE = procsN + 1
        procsSW = procsS - 1
        procsSE = procsS + 1
        IF( MOD(procID(kt), nprocX ) .EQ. 0 ) THEN
          procsW = -1
          procsNW = -1
          procsSW = -1
        ELSEIF( MOD(procID(kt), nprocX ) .EQ. nprocX-1 ) THEN
          procsE = -1
          procsNE = -1
          procsSE = -1
        ENDIF
        IF( procID(kt) .EQ. root .OR. procID(kt) .EQ. root+1 ) THEN
          ALLOCATE( procs(1:nTurbines, 1:9) )
          procs(kt, :) = (/procsW, procsE, procsN, procsS, procsNW, procsNE, procsSW, procsSE, procID(kt)/)
        ELSE
          ALLOCATE( procs(1:nTurbines, 1:10) )
          procs(kt, :) = (/procsW, procsE, procsN, procsS, procsNW, procsNE, procsSW, procsSE, root, procID(kt)/)
        ENDIF
        IF( ANY( procs(kt, :) .LT. root ) .OR. ANY( procs(kt, :) .GE. commSize ) ) THEN
          procsTmp = PACK( procs(kt, :), MASK = ( procs(kt, :) .GE. root .AND. procs(kt, :) .LT. commSize ) )
          ALLOCATE( procsCommTmp( 1:SIZE(procsTmp, DIM=1) ) )
          arrayIN = procsTmp
          CALL removeDupsINT(arrayIN, arrayOUT, nSize)
          procsCommTmp = arrayOUT(1:nSize)
          ALLOCATE( procsComm( 1:nTurbines, 1:nSize ) )
          procsComm(kt, 1:nSize) = procsCommTmp(1:nSize)
        ELSE
          ALLOCATE( procsCommTmp( 1:SIZE(procs, DIM=2) ) )
          arrayIN = procs(kt, :)
          CALL removeDupsINT(arrayIN, arrayOUT, nSize)
          procsCommTmp = arrayOUT(1:nSize)
          ALLOCATE( procsComm( 1:nTurbines, 1:nSize ) )
          procsComm(kt, 1:nSize) = procsCommTmp(1:nSize)
        ENDIF
        IF( errorCheckOpt ) THEN
          IF( FINDLOC(procsComm(kt, :), VALUE=procIDFstInit(kt), DIM=1) .EQ. 0 ) THEN
            CALL wrf_error_fatal('WTP-ERROR ***** Inflow point cannot be placed outside &
                                  of processor boundaries within communicators. &
                                  Decrease inflow distance in &
                                  windTurbines/WTNAME/turbineProperties.tbl or increase &
                                  domain size. *****')
          ENDIF
        ENDIF
        ! Get the group or processes of the default communicator
        CALL MPI_COMM_GROUP(MPI_COMM_WORLD, worldGroup, ierr)
        ! Keep only the processes of "procsComm" in the new group.
        CALL MPI_GROUP_INCL(worldGroup, SIZE( procsComm(kt, 1:nSize) ), procsComm(kt, 1:nSize), newGroup(kt), ierr)
        ! Create the new communicator from that group of processes.
        CALL MPI_COMM_CREATE(MPI_COMM_WORLD, newGroup(kt), newComm(kt), ierr)
        sender(kt) = FINDLOC(procsComm(kt, 1:nSize), VALUE=procID(kt), DIM=1) - 1
        DEALLOCATE(procs, procsCommTmp, procsComm)
      ELSE
        sender(kt) = procID(kt)
        newComm(kt) = MPI_COMM_WORLD
      ENDIF
      IF ( wrf_dm_on_monitor() ) THEN
        IF ( config_flags%windturbines_latlon .EQ. 0 ) THEN
          WRITE(message_windturbine, "(A, I5, A, F9.3, A, A, F9.3, A, A, F9.3, A, A, F9.3, A)")             &
                                     'Wind Turbine #', kt, ': LXrotor = ', rotorApexInit(kt, 1), ' m',      &
                                     ', LYrotor = ',  rotorApexInit(kt, 2), ' m',                           &
                                     ', LXtower = ', towerBaseLocation(kt, 1), ' m', ', LYtower = ',        &
                                     towerBaseLocation(kt, 2), ' m'
          CALL wrf_message(message_windturbine)
          WRITE(message_windturbine, "(A, I9, A, I9, A, I9, A, I9, A, I4)")                                  &
                                     '                     Irotor = ', iRot(kt), '  ,  Jrotor = ', jRot(kt), &
                                     '  ,  Itower = ', iTow(kt), '  ,  Jtower = ', jTow(kt),                 &
                                     '  , ProcID = ', procID(kt)
          CALL wrf_message(message_windturbine)
        ELSE
          WRITE(message_windturbine, "(A, I5, A, F13.8, A, F13.8, A, I4, A, I4, A, I4)")                            &
                                     'Wind Turbine #', kt ,': Lat = ', lat(kt), ', Lon = ', lon(kt), ', Irotor = ', &
                                      iRot(kt), ', Jrotor = ', jRot(kt), ', ProcID = ', procID(kt)
          CALL wrf_message(message_windturbine)
        ENDIF
      ENDIF
#endif
    ENDDO
!=======================================================================
#ifdef DM_PARALLEL
    CALL wrf_dm_bcast_real(hubHeight, nTurbineTypes)
    CALL wrf_dm_bcast_real(uvShaftDir, nTurbineTypes)
    CALL wrf_dm_bcast_real(delta, nTurbineTypes)
    CALL wrf_dm_bcast_real(beta, nTurbineTypes)
    CALL wrf_dm_bcast_real(xMinTrb, nTurbines)
    CALL wrf_dm_bcast_real(xMaxTrb, nTurbines)
    CALL wrf_dm_bcast_real(yMinTrb, nTurbines)
    CALL wrf_dm_bcast_real(yMaxTrb, nTurbines)
    CALL wrf_dm_bcast_real(uvTower, 3*nTurbines)
    CALL wrf_dm_bcast_real(uvShaft, 3*nTurbines)
    CALL wrf_dm_bcast_real(towerBaseLocation, 3*nTurbines)
    CALL wrf_dm_bcast_real(posRotorApex, 3*nTurbines)
    CALL wrf_dm_bcast_real(rotorApexInit, 3*nTurbines)
    CALL wrf_dm_bcast_real(rotorApexLocal, 3*nTurbines)
    CALL wrf_dm_bcast_real(bladePoints, 3*nTurbines*nSections*nElements)
    CALL wrf_dm_bcast_real(nacellePoints, 3*nTurbines*nNacPoints)
    CALL wrf_dm_bcast_real(towerPoints, 3*nTurbines*nTowPoints)
    CALL wrf_dm_bcast_real(planePoints, 3*nPlnWT*nVertPoints*nHorizPoints)
    CALL wrf_dm_bcast_integer(procID, nTurbines)
    CALL wrf_dm_bcast_integer(iRot, nTurbines)
    CALL wrf_dm_bcast_integer(jRot, nTurbines)
#endif
  END SUBROUTINE init_module_wind_gal
!=======================================================================
  SUBROUTINE init_wind_plane_tsout_locs(          &
        &  nPlnWT, nPlnWTcounter, nz, nxy, dx, dy &
        &, planeWidth, planeHeight, dist_plane    &
        &, towerShaftIntersect, uvTower           &
        &, capPhi, Lx, Ly, Lz, planePoints        &
        &)
!
  IMPLICIT NONE
!
    INTEGER, INTENT(IN) :: nPlnWT, nPlnWTcounter, nz, nxy
    REAL, INTENT(IN) :: planeWidth, planeHeight, dist_plane
    REAL, INTENT(IN) :: dx, dy, capPhi, Lx, Ly, Lz
    REAL, DIMENSION(1:3), INTENT(IN) :: towerShaftIntersect, uvTower
    REAL, DIMENSION(1:nPlnWT, 1:nz, 1:nxy, 1:3), INTENT(INOUT) :: planePoints
! local
    INTEGER :: k, n, m
    REAL :: dist_z, dist_xy
    REAL :: dxy, dz
!
    dxy = planeWidth/nxy
    dz = planeHeight/nz
    dist_xy = 0.0
    k = nPlnWTcounter
    DO m = 1, nxy
      dist_xy = dist_xy + 0.5*dxy
      dist_z = 0.0
      DO n = 1, nz
        dist_z = dist_z + 0.5*dz
        planePoints(k, n, m, 1) = Lx + dist_plane
        planePoints(k, n, m, 2) = Ly - (planeWidth/2) + dist_xy
        planePoints(k, n, m, 3) = Lz + dist_z
        dist_z = dist_z + 0.5*dz
        planePoints(k, n, m, :) = rotatePoint( planePoints(k, n, m, :), towerShaftIntersect, uvTower, capPhi )
      ENDDO
      dist_xy = dist_xy + 0.5*dxy
    ENDDO
    IF( errorCheckOpt ) THEN
    ! spec_bdy_width = 5
      IF( (Lx-dist_plane) .LE. (5*dx) .OR. ( Ly-(planeWidth/2) ) .LE. (5*dy) ) THEN
        CALL wrf_error_fatal('WTP-ERROR ***** Plane must not be placed &
                              very close to domain boundaries. *****')
      ENDIF
    ENDIF
  END SUBROUTINE init_wind_plane_tsout_locs
!=======================================================================
  SUBROUTINE calc_wind_plane_tsout_vel(    &
        &  ims, ime, jms, jme, kms, kme    &
        &, its, itf, jts, jtf, kts, ktf    &
        &, nPlnWT, nPlnWTcounter, nz, nxy  &
        &, dx, dy, dz_cst, dz, wt_do_yaw   &
        &, plane_points                    &
        &, dist_wt_base_to_hub             &
        &, wt_unit_vector, wt_delta_yaw    &
        &, c_z, w_z, u_2, v_2, w_2         &
        &, mpi_new_comm, ierr              &
        &, u_plane_l, v_plane_l, w_plane_l &
        &)
!
  IMPLICIT NONE
#if defined(DM_PARALLEL) && !defined(STUBMPI)
  INCLUDE "mpif.h"
#endif
    LOGICAL, INTENT(IN) :: wt_do_yaw
    INTEGER, INTENT(IN) :: ims, ime, jms, jme, kms, kme
    INTEGER, INTENT(IN) :: nPlnWT, nz, nxy, its, itf, jts, jtf, kts, ktf
    INTEGER, INTENT(IN) :: nPlnWTcounter, mpi_new_comm, ierr
    REAL, INTENT(IN) :: dx, dy, dz_cst, wt_delta_yaw
    REAL, DIMENSION(1:3), INTENT(IN) :: dist_wt_base_to_hub, wt_unit_vector
    REAL, DIMENSION(1:nPlnWT, 1:nz, 1:nxy, 1:3), INTENT(INOUT) :: plane_points
    REAL, DIMENSION(its:itf, kts:ktf, jts:jtf), INTENT(IN) :: dz, c_z, w_z
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: u_2, v_2, w_2
    REAL, DIMENSION(1:nPlnWT, 1:nz, 1:nxy), INTENT(OUT) :: u_plane_l, v_plane_l, w_plane_l
! local
    INTEGER :: k, m, n, ii, jj, kk
    REAL, DIMENSION(1:nz, 1:nxy) :: u_plane_les, v_plane_les, w_plane_les
    REAL, DIMENSION(1:nz, 1:nxy) :: u_plane_les_l, v_plane_les_l, w_plane_les_l
!
    k = nPlnWTcounter
    ! Loop over plane element lateral positions:
    DO m = 1, nxy
      ! Loop over plane element vertical positions:
      DO n = 1, nz
        u_plane_les(n, m) = 0.0
        v_plane_les(n, m) = 0.0
        w_plane_les(n, m) = 0.0
        u_plane_les_l(n, m) = 0.0
        v_plane_les_l(n, m) = 0.0
        w_plane_les_l(n, m) = 0.0
        IF( wt_do_yaw ) THEN
          ! Update plane element position based on yaw motion:
          plane_points(k, n, m, :) = rotatePoint(plane_points(k, n, m, :), dist_wt_base_to_hub, wt_unit_vector, wt_delta_yaw)
        ENDIF
        ! Interpolation of the u-component:
        ii = INT( (plane_points(k, n, m, 1)/dx ) )
        jj = INT( (plane_points(k, n, m, 2) + 0.5*dy)/dy )
        IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. &
            ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
          kk = MINLOC( ABS(c_z(ii, :, jj) + 0.5*dz_cst - plane_points(k, n, m, 3) ), 1)
          IF( ( c_z(ii, kk, jj) + 0.5*dz_cst ) .GT. plane_points(k, n, m, 3) ) kk = kk - 1
          IF( kk .EQ. 0 ) kk = kts
          CALL interpolateVelocity(                                 &
              &  ii, jj, kk, dx, dy, plane_points(k, n, m, 1)       &
              &, plane_points(k, n, m, 2), plane_points(k, n, m, 3) &
              &, u_2(ii:ii+1, kk:kk+1, jj:jj+1), dz(ii, kk, jj)     &
              &, c_z(ii, kk, jj), 'u', u_plane_les_l(n, m)          &
              &)
        ELSE
          u_plane_les_l(n, m) = -999999.999
        ENDIF
        ! Interpolation of the v-component:
        ii = INT( (plane_points(k, n, m, 1) + 0.5*dx)/dx )
        jj = INT( (plane_points(k, n, m, 2)/dy ) )
        IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. &
            ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
          kk = MINLOC( ABS(c_z(ii, kts:ktf, jj) + 0.5*dz_cst - plane_points(k, n, m, 3) ), 1)
          IF( ( c_z(ii, kk, jj) + 0.5*dz_cst ) .GT. plane_points(k, n, m, 3) ) kk = kk - 1
          IF( kk .EQ. 0 ) kk = kts
          CALL interpolateVelocity(                                 &
              &  ii, jj, kk, dx, dy, plane_points(k, n, m, 1)       &
              &, plane_points(k, n, m, 2), plane_points(k, n, m, 3) &
              &, v_2(ii:ii+1, kk:kk+1, jj:jj+1), dz(ii, kk, jj)     &
              &, c_z(ii, kk, jj), 'v', v_plane_les_l(n, m)          &
              &)
        ELSE
          v_plane_les_l(n, m) = -999999.999
        ENDIF
        ! Interpolation of the w-component:
        ii = INT( (plane_points(k, n, m, 1) + 0.5*dx)/dx )
        jj = INT( (plane_points(k, n, m, 2) + 0.5*dy)/dy )
        IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. &
            ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
          kk = MINLOC( ABS(w_z(ii, kts:ktf, jj) - plane_points(k, n, m, 3) ), 1)
          IF( w_z(ii, kk, jj) .GT. plane_points(k, n, m, 3) ) kk = kk - 1
          IF( kk .EQ. 0 ) kk = kts
          CALL interpolateVelocity(                                 &
              &  ii, jj, kk, dx, dy, plane_points(k, n, m, 1)       &
              &, plane_points(k, n, m, 2), plane_points(k, n, m, 3) &
              &, w_2(ii:ii+1, kk:kk+1, jj:jj+1), dz(ii, kk, jj)     &
              &, w_z(ii, kk, jj), 'w', w_plane_les_l(n, m)          &
              &)
        ELSE
          w_plane_les_l(n, m) = -999999.999
        ENDIF
      ENDDO
    ENDDO
    CALL MPI_ALLREDUCE( u_plane_les_l, u_plane_les, nxy*nz, MPI_REAL, MPI_MAX, mpi_new_comm, ierr )
    CALL MPI_ALLREDUCE( v_plane_les_l, v_plane_les, nxy*nz, MPI_REAL, MPI_MAX, mpi_new_comm, ierr )
    CALL MPI_ALLREDUCE( w_plane_les_l, w_plane_les, nxy*nz, MPI_REAL, MPI_MAX, mpi_new_comm, ierr )
    ! Longitudinal velocity at the plane:
    u_plane_l(k, :, :) = u_plane_les
    ! Lateral velocity at the plane:
    v_plane_l(k, :, :) = v_plane_les
    ! Vertical velocity at the plane:
    w_plane_l(k, :, :) = w_plane_les
  END SUBROUTINE calc_wind_plane_tsout_vel
!=======================================================================
  SUBROUTINE removeDupsINT(arrayIN, arrayOUT, nSize)
  IMPLICIT NONE
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN)  :: arrayIN
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: arrayOUT
    INTEGER, INTENT(OUT) :: nSize
! local
    INTEGER :: i, j, k
    INTEGER, DIMENSION( SIZE(arrayIN, DIM=1) ) :: res
    k = 1
    res(1) = arrayIN(1)
    outer: DO i = 2, SIZE(arrayIN, DIM=1)
      DO j = 1, k
        IF( res(j) .EQ. arrayIN(i) ) THEN
          ! Found a match so start looking again
          CYCLE outer
        ENDIF
      ENDDO
      ! No match found so add it to the output
      k = k + 1
      res(k) = arrayIN(i)
    ENDDO outer
    nSize = k
    arrayOUT = res(1:nSize)
  END SUBROUTINE removeDupsINT
!=======================================================================
  SUBROUTINE removeDupsSTR(arrayIN, arrayOUT, nSize)
  IMPLICIT NONE
    CHARACTER*256, DIMENSION(:), ALLOCATABLE, INTENT(IN)  :: arrayIN
    CHARACTER*256, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: arrayOUT
    INTEGER, INTENT(OUT) :: nSize
! local
    INTEGER :: i, j, k
    CHARACTER*256, DIMENSION( SIZE(arrayIN, DIM=1) ) :: res
    k = 1
    res(1) = arrayIN(1)
    outer: DO i = 2, SIZE(arrayIN, DIM=1)
      DO j = 1, k
        IF( res(j) .EQ. arrayIN(i) ) THEN
          ! Found a match so start looking again
          CYCLE outer
        ENDIF
      ENDDO
      ! No match found so add it to the output
      k = k + 1
      res(k) = arrayIN(i)
    ENDDO outer
    nSize = k
    arrayOUT = res(1:nSize)
  END SUBROUTINE removeDupsSTR
!=======================================================================
  SUBROUTINE fileCheck(filename, exists)
  IMPLICIT NONE
    CHARACTER(LEN = *), INTENT(IN) :: filename
    LOGICAL, INTENT(OUT) :: exists
    INQUIRE(FILE="./"//TRIM(filename), EXIST=exists)
  END SUBROUTINE fileCheck
!=======================================================================
END MODULE module_wind_gal

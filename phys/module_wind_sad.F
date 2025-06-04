!WRF:MODEL_LAYER:PHYSICS
MODULE module_wind_sad
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This is the main module for simplified actuator disk (SAD) scheme.  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================================================================
! nampelist.input options:
!   &physics
!   wind_wtp_opt(max_domains)           = 4,        ! 0 = deactivates GAL/GAD/GADrs/SAD models
!                                                     1 = activates GAL model
!                                                     2 = activates GAD model
!                                                     3 = activates GADrs model
!                                                     4 = activates SAD model
!   wind_wtp_spec_rho                   = -1.0,     ! <=0 air density calculated
!                                                         by the model
!                                                     >0  user-specified air density
!                                                     >=0 user-specified direction
!   wind_wtp_printdiag_opt              = 0,        ! 0 = no diagnostics printed
!                                                     1 = only necessary information
!   wind_wtp_yaw_opt                    = .true.,   ! turbine yaw mechanism is turned on
!   wind_wtp_nTurbines                  = 10,       ! number of wind turbines
!   wind_wtp_tke_fact                   = 0.25,     ! turbulence kinetic energy factor
!                                                     constant
!   wind_wtp_dist_cst                   = 1.0,      ! Constant to calculate distance between
!                                                     wind turbine and freestream location.
!                                                     dist=wind_wtp_dist_cst*rotor_diameter
!   wind_wtp_yaw_rate                   = 0.5,      ! Wind turbine yaw rate [deg/s]
!   wind_wtp_yaw_threshold              = 1000.0,   ! Yaw error threshold beyond which wind
!                                                     turbine starts to yaw [deg^2/s]
!   wind_wtp_aif_opt                    = 1,        ! Switch between axial induction factor
!                                                     calculation methods
!                                                     1 = direct calculation axial induction
!                                                         factor from hub height freestream
!                                                         wind velocity and hub height wind
!                                                         velocity at rotor (recommended)
!                                                     2 = freestream wind velocity is calculated
!                                                         from tabulated thrust coefficient and
!                                                         hub height horizontal wind velocity
!                                                         at the rotor disk.
!   wind_wtp_nAvesamples                = 100,      ! number of samples for running average
!   wind_wtp_nSeries_tme                = 1000,     ! length of storage matrix,
!                                                      wind_wtp_nTurbines &
!                                                     *wind_wtp_nAvesamples
!   windfarm_ij                         = 0,        ! 1 = turbine location in SI units
!                                                         (in meters), used for ideal cases
!                                                     0 = turbine location in lat and lon
!                                                         format, used for real cases
!   &dynamics
!   wind_wtp_sad_tsout_opt              = .false.,  ! time series output for wind turbine
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
!   wind_wtp_pln_nTurbines              = 3,        ! total number of wind turbines for which
!                                                     time series of velocity components to be
!                                                     calculated
!   wind_wtp_pln_turbine_ids            = 1, 5, 9,  ! wind turbine id for which time series
!                                                     of velocity components to be calculated
!   wind_wtp_pln_nz                     = 20,       ! number of plane-element ponts in the
!                                                     vertical direction
!   wind_wtp_pln_nxy                    = 20,       ! number of plane-element points in the
!                                                     horizontal space
!   wind_wtp_pln_nSeries_tme            = 1200,     ! length of storage matrix,
!                                                      wind_wtp_pln_nTurbines &
!                                                     *wind_wtp_pln_nz &
!                                                     *wind_wtp_pln_nxy
!   wind_wtp_pln_width                  = 300.0,    ! plane width in meter
!   wind_wtp_pln_height                 = 300.0,    ! plane height in meter
!   wind_wtp_pln_dist                   = 100.0,    ! distance between reference plane and
!                                                     corresponding wind turbine
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
!      P0       : mpi_proc_root processor (also involved in the communication)  
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
!      |<---------------------n_proc_x--------------------->|
!      REMARKS:
!       (1) This procedure is repeated for each turbine.
!       (2) If turbine is located at P0, "mpi_new_comm" includes only
!           processors P1, P4, P5 and broadcasts to 3 processors.
!       (3) If turbine is located within one of the corner processors 
!           except P0 (i.e., P3, P12, or P15), "mpi_new_comm" includes
!           contiguous processors + P0 and broadcasts to 4 processors.
!       (4) If turbine is located within one of the "near-wall" processors 
!           except from corners (i.e., P1, P2, P4, P7, P8, P11, P13, P14),
!           "mpi_new_comm" includes contiguous processors + P0 and broadcasts
!           to 6 processors. For instance, if a turbine is located at P4,
!           "mpi_new_comm" will include processors P0, P1, P4, P5, P8, and P9.
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
!
#ifdef DM_PARALLEL
  INTEGER, PARAMETER :: mpi_proc_root = 0 ! Specifies the root processor
#endif
  INTEGER, PARAMETER :: n_max_val = 1E2 ! Maximum allowed vector size
!
  LOGICAL :: opt_yaw, opt_tsout, opt_tsout_plane
  LOGICAL :: opt_error_check, rst_controller, rst_controller1
!
  INTEGER :: opt_aif, opt_print_diag, n_wt, n_wt_types, n_ave_samples
  INTEGER :: ils, ile, jls, jle, kls, kle
  INTEGER :: n_pln_wt, n_pln_wt_counter, n_horz_points, n_vert_points
  INTEGER :: sad_max_points, row_max_points
!
#ifdef DM_PARALLEL
  INTEGER :: ierr, mpi_size_comm, mpi_world_group
  INTEGER :: procs_W, procs_E, procs_N, procs_S
  INTEGER :: procs_NW, procs_NE, procs_SW, procs_SE
#endif
!
!  REAL :: max_radius, min_dx, min_dz
  REAL :: plane_width, plane_height, plane_distance, spec_rho
  REAL :: wt_ck_cst
!
  INTEGER, DIMENSION(:), ALLOCATABLE :: wt_type_id, wt_i, wt_j
  INTEGER, DIMENSION(:), ALLOCATABLE :: wt_tab_n_lines, spec_wt
  INTEGER, DIMENSION(:), ALLOCATABLE :: wt_k_bot, wt_k_hub, wt_k_top
  INTEGER, DIMENSION(:), ALLOCATABLE :: sad_n_points, sad_n_points_l
  INTEGER, DIMENSION(:), ALLOCATABLE :: fst_i_l, fst_j_l
!
  INTEGER, DIMENSION(:, :), ALLOCATABLE :: sad_i_l, sad_j_l, sad_k_l
  INTEGER, DIMENSION(:, :), ALLOCATABLE :: sad_i, sad_j, sad_k
!
#ifdef DM_PARALLEL
  INTEGER, DIMENSION(:), ALLOCATABLE :: mpi_sender, mpi_proc, mpi_proc_id, mpi_proc_fst_init
  INTEGER, DIMENSION(:), ALLOCATABLE :: mpi_proc_id_fst_init, mpi_proc_fst, mpi_proc_id_fst
  INTEGER, DIMENSION(:), ALLOCATABLE :: mpi_new_group, mpi_new_comm
#endif
!
  REAL, DIMENSION(:), ALLOCATABLE :: wt_aif_l, wt_ct_l, wt_cp_l, wt_ck_l, wt_power_l
  REAL, DIMENSION(:), ALLOCATABLE :: wt_yaw_angle_l, fst_uinf_l
  REAL, DIMENSION(:), ALLOCATABLE :: wt_vel_hub, wt_ck, dz_cst
  REAL, DIMENSION(:), ALLOCATABLE :: wt_x, wt_y, wt_z, wt_lat, wt_lon
  REAL, DIMENSION(:), ALLOCATABLE :: wt_delta_yaw_angle_l
  REAL, DIMENSION(:), ALLOCATABLE :: wt_yaw_rate, wt_yaw_error_l
  REAL, DIMENSION(:), ALLOCATABLE :: wt_yaw_error_threshold, wt_cutin, wt_cutout
  REAL, DIMENSION(:), ALLOCATABLE :: wt_hub_height, wt_mod_hub_height
  REAL, DIMENSION(:), ALLOCATABLE :: wt_rotor_diameter, wt_std_ct, wt_rated_power
  REAL, DIMENSION(:), ALLOCATABLE :: wt_norm_radius, wt_area, wt_dir_l
  REAL, DIMENSION(:), ALLOCATABLE :: fst_dist_cst, fst_wdir_ave_l, fst_rho_ave_l
!
  LOGICAL, DIMENSION(:), ALLOCATABLE :: wt_do_yaw, wt_do_yaw_l
!
  REAL, DIMENSION(:, :), ALLOCATABLE :: wt_unit_vector, dist_wt_base_to_hub
  REAL, DIMENSION(:, :), ALLOCATABLE :: sad_area_fact_l
  REAL, DIMENSION(:, :), ALLOCATABLE :: wt_tab_wsph, wt_tab_power, wt_tab_ct
  REAL, DIMENSION(:, :), ALLOCATABLE :: wt_tab_wsph_1, wt_tab_aif
  REAL, DIMENSION(:, :), ALLOCATABLE :: fst_wdir_ts_l, fst_u_ts_l
  REAL, DIMENSION(:, :), ALLOCATABLE :: fst_v_ts_l, fst_rho_ts_l
!
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: u_plane_l, v_plane_l, w_plane_l
!
  REAL, DIMENSION(:, :, :, :), ALLOCATABLE :: plane_points
CONTAINS
  SUBROUTINE calc_power_ct_sad(                     &
        &  dx, dy, rdzw, dt, muu, muv, mut, ht      &
        &, ru_tendf, rv_tendf, tke_tend             &
        &, u_2, v_2, w_2                            &
        &, ids, ide, jds, jde, kds, kde             &
        &, ims, ime, jms, jme, kms, kme             &
        &, its, ite, jts, jte, kts, kte             &
        &, fst_i, fst_j, fst_wdir_ave, fst_rho_ave  &
        &, fst_uinf, fst_rho_ts, fst_u_ts, fst_v_ts &
        &, fst_wdir_ts, wt_dir, wt_yaw_error, rho   &
        &, wt_power, wt_cp, wt_ct, wt_aif           &
        &, wtp_history_flag, wtp_restart_flag       &
        &, read_from_rst_file, u_plane, v_plane     &
        &, w_plane, plane_points_x, plane_points_y  &
        &, plane_points_z                           &
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
!
    LOGICAL, INTENT(IN) :: wtp_history_flag, wtp_restart_flag
!
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: muu, muv, mut, ht
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: rdzw, rho, u_2, v_2, w_2
!
    INTEGER, INTENT(INOUT) :: read_from_rst_file
!
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT) :: ru_tendf, rv_tendf, tke_tend
!
    REAL, DIMENSION(1:n_wt), INTENT(OUT) :: fst_wdir_ave, fst_rho_ave, fst_uinf
    REAL, DIMENSION(1:n_wt), INTENT(OUT) :: wt_dir, wt_power, wt_cp, wt_ct, wt_aif
    REAL, DIMENSION(1:n_wt), INTENT(INOUT) :: wt_yaw_error
!
    INTEGER, DIMENSION(1:n_wt), INTENT(INOUT) :: fst_i, fst_j
!
    REAL, DIMENSION(1:n_wt, 1:n_ave_samples), INTENT(INOUT) :: fst_rho_ts, fst_u_ts
    REAL, DIMENSION(1:n_wt, 1:n_ave_samples), INTENT(INOUT) :: fst_v_ts, fst_wdir_ts
!
    REAL, DIMENSION(1:n_pln_wt, 1:n_vert_points, 1:n_horz_points), INTENT(OUT) :: u_plane
    REAL, DIMENSION(1:n_pln_wt, 1:n_vert_points, 1:n_horz_points), INTENT(OUT) :: v_plane
    REAL, DIMENSION(1:n_pln_wt, 1:n_vert_points, 1:n_horz_points), INTENT(OUT) :: w_plane
    REAL, DIMENSION(1:n_pln_wt, 1:n_vert_points, 1:n_horz_points), INTENT(OUT) :: plane_points_x
    REAL, DIMENSION(1:n_pln_wt, 1:n_vert_points, 1:n_horz_points), INTENT(OUT) :: plane_points_y
    REAL, DIMENSION(1:n_pln_wt, 1:n_vert_points, 1:n_horz_points), INTENT(OUT) :: plane_points_z
!
    CHARACTER*256  :: outstring
! ... Local
    REAL, PARAMETER :: wt_not_yaw_threshold = 0.001 ! [deg]
    REAL, PARAMETER :: wt_yaw_angle_threshold = 0.100 ! [deg]
!
    REAL :: fst_rho_sum, fst_u_sum
    REAL :: fst_v_sum, wFstSum, fst_wdir_sum
    REAL :: fst_sin_sum, fst_cos_sum
    REAL :: grid_ijk_wsph, wt_aif_tmp
!
    INTEGER :: i, j, k, m, n
    INTEGER :: itf, jtf, ktf, ii, jj, kk
    INTEGER :: seg_k_bot, seg_k_top
!
    INTEGER, DIMENSION(1:n_wt) :: fst_k
!
    ! Rotation matrices: 
    REAL, DIMENSION(1:3) :: fst_uv, fst_uv_rot
!
    REAL, DIMENSION(1:n_wt) :: fst_rho, fst_wdir, fst_u, fst_v
    REAL, DIMENSION(1:n_wt) :: fst_rho_l, fst_wdir_l, fst_u_l, fst_v_l
    REAL, DIMENSION(1:n_wt) :: fst_x, fst_y, fst_z
    REAL, DIMENSION(1:n_wt) :: wt_yaw_dir, fst_u_ave_l, fst_v_ave_l
    REAL, DIMENSION(1:n_wt) :: wt_u_hub_l, wt_v_hub_l
    REAL, DIMENSION(1:n_wt) :: wt_u_hub, wt_v_hub, wt_wdir_hub
    REAL, DIMENSION(1:n_wt) :: wt_delta_yaw, wt_delta_yaw_l
    REAL, DIMENSION(1:n_wt) :: fst_hub_height_asl, fst_hub_height_asl_l
!
    REAL, DIMENSION(kts:kte-1) :: wt_Faz, wt_Rz
!
    REAL, DIMENSION(1:n_wt, 1:sad_max_points) :: sad_area_fact
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
    n_pln_wt_counter = 1 
    DO n = 1, n_wt
#ifdef DM_PARALLEL
      IF( mpi_new_comm(n) .NE. MPI_COMM_NULL ) THEN
        IF( read_from_rst_file .EQ. 0 ) THEN
          fst_rho_ts_l(n, :) = 0.0
          fst_u_ts_l(n, :) = 0.0
          fst_v_ts_l(n, :) = 0.0
          fst_wdir_ts_l(n, :) = 0.0
          wt_yaw_error_l(n) = 0.0
        ENDIF
        IF( rst_controller .AND. read_from_rst_file .EQ. 1 ) THEN
          fst_i_l(n) = fst_i(n)
          fst_j_l(n) = fst_j(n)
          fst_rho_ts_l(n, :) = fst_rho_ts(n, :)
          fst_u_ts_l(n, :) = fst_u_ts(n, :)
          fst_v_ts_l(n, :) = fst_v_ts(n, :)
          fst_wdir_ts_l(n, :) = fst_wdir_ts(n, :)
          wt_yaw_error_l(n) = wt_yaw_error(n)
        ENDIF      
        wt_u_hub(n) = 0.0
        wt_v_hub(n) = 0.0
        wt_u_hub_l(n) = 0.0
        wt_v_hub_l(n) = 0.0
        ! Interpolation of the u-component:
        ii = INT( (wt_x(n)/dx ) )
        jj = INT( (wt_y(n) + 0.5*dy)/dy )
        IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
          kk = MINLOC( ABS(c_z(ii, kts:ktf, jj) + 0.5*dz_cst(n) - wt_mod_hub_height( wt_type_id(n) ) ), 1)
          IF( ( c_z(ii, kk, jj) + 0.5*dz_cst(n) ) .GT. wt_mod_hub_height( wt_type_id(n) ) ) kk = kk - 1
          CALL interp_velocity(                       &
              &  ii, jj, kk, dx, dy, wt_x(n), wt_y(n) &
              &, wt_mod_hub_height( wt_type_id(n) )   &
              &, u_2(ii:ii+1, kk:kk+1, jj:jj+1)       &
              &, dz(ii, kk, jj), c_z(ii, kk, jj), 'u' &
              &, wt_u_hub_l(n)                        &
              &)
        ELSE
          wt_u_hub_l(n) = -999999.999
        ENDIF
        ! Interpolation of the v-component:
        ii = INT( (wt_x(n) + 0.5*dx)/dx )
        jj = INT( (wt_y(n)/dy ) )
        IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
          kk = MINLOC( ABS(c_z(ii, kts:ktf, jj) + 0.5*dz_cst(n) - wt_mod_hub_height( wt_type_id(n) ) ), 1)
          IF( ( c_z(ii, kk, jj) + 0.5*dz_cst(n) ) .GT. wt_mod_hub_height( wt_type_id(n) ) ) kk = kk - 1
          CALL interp_velocity(                       &
              &  ii, jj, kk, dx, dy, wt_x(n), wt_y(n) &
              &, wt_mod_hub_height( wt_type_id(n) )   &
              &, v_2(ii:ii+1, kk:kk+1, jj:jj+1)       &
              &, dz(ii, kk, jj), c_z(ii, kk, jj), 'v' &
              &, wt_v_hub_l(n)                        &
              &)
        ELSE
          wt_v_hub_l(n) = -999999.999
        ENDIF
        CALL MPI_ALLREDUCE( wt_u_hub_l(n), wt_u_hub(n), 1, MPI_REAL, MPI_MAX, mpi_new_comm(n), ierr )
        CALL MPI_ALLREDUCE( wt_v_hub_l(n), wt_v_hub(n), 1, MPI_REAL, MPI_MAX, mpi_new_comm(n), ierr )
!        wt_vel_hub(n) = SQRT( wt_u_hub(n)**2 + wt_v_hub(n)**2 )
        wt_wdir_hub(n) = ATAN2( wt_u_hub(n), wt_v_hub(n) )/degrad + 180.0
        CALL compass_to_standard( wt_wdir_hub(n) )
        ! Find local velocity matrix:
        fst_uv(1) = wt_u_hub(n)
        fst_uv(2) = wt_v_hub(n)
        fst_uv(3) = 0.0
        CALL rot_global_to_local(wt_wdir_hub(n), fst_uv, fst_uv_rot)
        ! Unperturbed hub-height freestream velocity:
        wt_vel_hub(n) = fst_uv_rot(1)
        IF( read_from_rst_file .EQ. 0 ) THEN
          wt_dir_l(n) = ATAN2( wt_u_hub(n), wt_v_hub(n) )/degrad + 180.0
          CALL compass_to_standard( wt_dir_l(n) )
        ELSE
          wt_dir_l(n) = wt_dir(n)
        ENDIF
        IF( wt_do_yaw(n) .OR. ( rst_controller1 .AND. .NOT. opt_yaw ) ) THEN
          wt_dir_l(n) = wt_dir(n)
        ENDIF
        IF( ( its .LE. wt_i(n) .AND. wt_i(n) .LE. itf ) .AND. ( jts .LE. wt_j(n) .AND. wt_j(n) .LE. jtf ) ) THEN
          fst_hub_height_asl_l(n) = ht( wt_i(n), wt_j(n) ) + wt_mod_hub_height( wt_type_id(n) )
        ELSE
          fst_hub_height_asl_l(n) = -999999.999
        ENDIF
        CALL MPI_ALLREDUCE( fst_hub_height_asl_l(n), fst_hub_height_asl(n), 1, MPI_REAL , MPI_MAX, mpi_new_comm(n), ierr )
        fst_u(n) = 0.0
        fst_v(n) = 0.0
        fst_rho(n) = 0.0
        fst_u_l(n) = 0.0
        fst_v_l(n) = 0.0
        fst_rho_l(n) = 0.0
        fst_z(n) = fst_hub_height_asl(n)
        ! Interpolation of the u-component:
        fst_x(n) = wt_x(n) - fst_dist_cst(n)*wt_rotor_diameter( wt_type_id(n) )*COS(wt_dir_l(n)*degrad)
        fst_y(n) = wt_y(n) - fst_dist_cst(n)*wt_rotor_diameter( wt_type_id(n) )*SIN(wt_dir_l(n)*degrad) + 0.5*dy
        ii = INT(fst_x(n)/dx)
        jj = INT(fst_y(n)/dy)
        ! Set horizontal boundary index in x-direction to i_domain_start if ii is less than or equal to ids
        IF( ii .LE. 0 ) THEN
          ii = ids
          fst_x(n) = (ids-1)*dx
        ENDIF
        ! Set horizontal boundary index in x-direction to i_domain_end-1 if ii is greater than ide
        IF( ii .GE. ide ) THEN
          ii = ide-1
          fst_x(n) = (ide-1)*dx
        ENDIF
        ! Set horizontal boundary index in y-direction to j_domain_start if jj is less than or equal to jds
        IF( jj .LE. 0 ) THEN
          jj = jds
          fst_y(n) = (jds-1)*dy
        ENDIF
        ! Set horizontal boundary index in y-direction to j_domain_end-1 if jj is greater than jde
        IF( jj .GE. jde ) THEN
          jj = jde-1
          fst_y(n) = (jde-1)*dy
        ENDIF
        IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
          kk = MINLOC( ABS( c_z(ii, kts:ktf, jj) + 0.5*dz_cst(n) + ht(ii, jj) - fst_z(n) ), 1 )
          IF( ( c_z(ii, kk, jj) + 0.5*dz_cst(n) + ht(ii, jj) ) .GT. fst_z(n) ) kk = kk - 1
          CALL interp_velocity(                                   &
              &  ii, jj, kk, dx, dy, fst_x(n), fst_y(n)           &
              &, fst_z(n), u_2(ii:ii+1, kk:kk+1, jj:jj+1)         &
              &, dz(ii, kk, jj), c_z(ii, kk, jj), 'u', fst_u_l(n) &
              &)
        ELSE
          fst_u_l(n) = -999999.999
        ENDIF
        ! Interpolation of the v-component:
        fst_x(n) = wt_x(n) - fst_dist_cst(n)*wt_rotor_diameter( wt_type_id(n) )*COS(wt_dir_l(n)*degrad) + 0.5*dx
        fst_y(n) = wt_y(n) - fst_dist_cst(n)*wt_rotor_diameter( wt_type_id(n) )*SIN(wt_dir_l(n)*degrad)
        ii = INT(fst_x(n)/dx)
        jj = INT(fst_y(n)/dy)
        ! Set horizontal boundary index in x-direction to i_domain_start if ii is less than or equal to ids
        IF( ii .LE. 0 ) THEN
          ii = ids
          fst_x(n) = (ids-1)*dx
        ENDIF
        ! Set horizontal boundary index in x-direction to i_domain_end-1 if ii is greater than ide
        IF( ii .GE. ide ) THEN
          ii = ide-1
          fst_x(n) = (ide-1)*dx
        ENDIF
        ! Set horizontal boundary index in y-direction to j_domain_start if jj is less than or equal to jds
        IF( jj .LE. 0 ) THEN
          jj = jds
          fst_y(n) = (jds-1)*dy
        ENDIF
        ! Set horizontal boundary index in y-direction to j_domain_end-1 if jj is greater than jde
        IF( jj .GE. jde ) THEN
          jj = jde-1
          fst_y(n) = (jde-1)*dy
        ENDIF
        IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
          kk = MINLOC( ABS( c_z(ii, kts:ktf, jj) + 0.5*dz_cst(n) + ht(ii, jj) - fst_z(n) ), 1 )
          IF( ( c_z(ii, kk, jj) + 0.5*dz_cst(n) + ht(ii, jj) ) .GT. fst_z(n) ) kk = kk - 1
          CALL interp_velocity(                                   &
              &  ii, jj, kk, dx, dy, fst_x(n), fst_y(n)           &
              &, fst_z(n), v_2(ii:ii+1, kk:kk+1, jj:jj+1)         &
              &, dz(ii, kk, jj), c_z(ii, kk, jj), 'v', fst_v_l(n) &
              &)
        ELSE
          fst_v_l(n) = -999999.999
        ENDIF
        ! Interpolation of density:
        fst_x(n) = wt_x(n) - fst_dist_cst(n)*wt_rotor_diameter( wt_type_id(n) )*COS(wt_dir_l(n)*degrad) + 0.5*dx
        fst_y(n) = wt_y(n) - fst_dist_cst(n)*wt_rotor_diameter( wt_type_id(n) )*SIN(wt_dir_l(n)*degrad) + 0.5*dy
        ii = INT(fst_x(n)/dx)
        jj = INT(fst_y(n)/dy)
        ! Set horizontal boundary index in x-direction to i_domain_start if ii is less than or equal to ids
        IF( ii .LE. 0 ) THEN
          ii = ids
          fst_x(n) = (ids-1)*dx
        ENDIF
        ! Set horizontal boundary index in x-direction to i_domain_end-1 if ii is greater than ide
        IF( ii .GE. ide ) THEN
          ii = ide-1
          fst_x(n) = (ide-1)*dx
        ENDIF
        ! Set horizontal boundary index in y-direction to j_domain_start if jj is less than or equal to jds
        IF( jj .LE. 0 ) THEN
          jj = jds
          fst_y(n) = (jds-1)*dy
        ENDIF
        ! Set horizontal boundary index in y-direction to j_domain_end-1 if jj is greater than jde
        IF( jj .GE. jde ) THEN
          jj = jde-1
          fst_y(n) = (jde-1)*dy
        ENDIF
        IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
          kk = MINLOC( ABS( c_z(ii, kts:ktf, jj) + 0.5*dz_cst(n) + ht(ii, jj) - fst_z(n) ), 1 )
          IF( ( c_z(ii, kk, jj) + 0.5*dz_cst(n) + ht(ii, jj) ) .GT. fst_z(n) ) kk = kk - 1
          CALL interp_scalar(                                  &
              &  ii, jj, kk, dx, dy, fst_x(n), fst_y(n)        &
              &, fst_z(n), rho(ii:ii+1, kk:kk+1, jj:jj+1)      &
              &, dz(ii, kk, jj), c_z(ii, kk, jj), fst_rho_l(n) &
              &)
        ELSE
          fst_rho_l(n) = -999999.999
        ENDIF
        CALL MPI_ALLREDUCE( fst_u_l(n), fst_u(n), 1, MPI_REAL, MPI_MAX, mpi_new_comm(n), ierr )
        CALL MPI_ALLREDUCE( fst_v_l(n), fst_v(n), 1, MPI_REAL, MPI_MAX, mpi_new_comm(n), ierr )
        CALL MPI_ALLREDUCE( fst_rho_l(n), fst_rho(n), 1, MPI_REAL, MPI_MAX, mpi_new_comm(n), ierr )
        fst_i_l(n) = ii
        fst_j_l(n) = jj
        fst_wdir(n) = ATAN2( fst_u(n), fst_v(n) )/degrad + 180.0
        ! Shift rho,u,v, and w for each turbine back by one unit to
        ! make room for the current values: 
        DO m = n_ave_samples-1, 1, -1
          fst_rho_ts_l(n, m+1) = fst_rho_ts_l(n, m)
          fst_u_ts_l(n, m+1) = fst_u_ts_l(n, m)
          fst_v_ts_l(n, m+1) = fst_v_ts_l(n, m)
          fst_wdir_ts_l(n, m+1) = fst_wdir_ts_l(n, m)
        ENDDO
        ! Sum all values to check whether this is the
        ! first timestep or not:
        fst_rho_sum = 0.0
        fst_u_sum = 0.0
        fst_v_sum = 0.0
        fst_wdir_sum = 0.0
        DO m = 1, n_ave_samples
          fst_rho_sum = fst_rho_sum + fst_rho_ts_l(n, m)
          fst_u_sum = fst_u_sum + fst_u_ts_l(n, m)
          fst_v_sum = fst_v_sum + fst_v_ts_l(n, m)
          fst_wdir_sum = fst_wdir_sum + fst_wdir_ts_l(n, m)
        ENDDO
        IF( fst_rho_sum .EQ. 0.0 ) fst_rho_ts_l(n, :) = (/(fst_rho(n), m=1, n_ave_samples)/)
        IF( fst_u_sum .EQ. 0.0 ) fst_u_ts_l(n, :) = (/(fst_u(n), m=1, n_ave_samples)/)
        IF( fst_v_sum .EQ. 0.0 ) fst_v_ts_l(n, :) = (/(fst_v(n), m=1, n_ave_samples)/)
        IF( fst_wdir_sum .EQ. 0.0 ) fst_wdir_ts_l(n, :) = (/(fst_wdir(n), m=1, n_ave_samples)/)
        ! Assign the first elements to current values:
        fst_rho_ts_l(n, 1) = fst_rho(n)
        fst_u_ts_l(n, 1) = fst_u(n)
        fst_v_ts_l(n, 1) = fst_v(n)
        fst_wdir_ts_l(n, 1) = fst_wdir(n)
        ! Running time-average:
        fst_rho_sum = 0.0
        fst_u_sum = 0.0
        fst_v_sum = 0.0
        fst_wdir_sum = 0.0
        fst_sin_sum = 0.0
        fst_cos_sum = 0.0          
        DO m = 1, n_ave_samples
          fst_rho_sum = fst_rho_sum + fst_rho_ts_l(n, m)
          fst_u_sum = fst_u_sum + fst_u_ts_l(n, m)
          fst_v_sum = fst_v_sum + fst_v_ts_l(n, m)
          fst_wdir_sum = fst_wdir_sum + fst_wdir_ts_l(n, m)
          fst_sin_sum = fst_sin_sum + SIN( fst_wdir_ts_l(n, m)*degrad )
          fst_cos_sum = fst_cos_sum + COS( fst_wdir_ts_l(n, m)*degrad )
        ENDDO
        fst_rho_ave_l(n) = fst_rho_sum/REAL(n_ave_samples)
        fst_u_ave_l(n) = fst_u_sum/REAL(n_ave_samples)
        fst_v_ave_l(n) = fst_v_sum/REAL(n_ave_samples)
        fst_wdir_ave_l(n) = MOD(ATAN2(fst_sin_sum, fst_cos_sum)/degrad + 360.0, 360.0)
        CALL compass_to_standard( fst_wdir_ave_l(n) )
        wt_yaw_angle_l(n) = wt_dir_l(n)
        ! Find local velocity matrix:
        fst_uv(1) = fst_u_ave_l(n)
        fst_uv(2) = fst_v_ave_l(n)
        fst_uv(3) = 0.0
        CALL rot_global_to_local(wt_yaw_angle_l(n), fst_uv, fst_uv_rot)
        ! Unperturbed hub-height freestream velocity:
        fst_uinf_l(n) = fst_uv_rot(1)
        IF( ( its .LE. wt_i(n) .AND. wt_i(n) .LE. itf ) .AND. ( jts .LE. wt_j(n) .AND. wt_j(n) .LE. jtf ) ) THEN
          IF( wt_do_yaw(n) ) THEN
            wt_yaw_error_l(n) = wt_yaw_error(n)
            wt_dir_l(n) = wt_dir(n)
          ENDIF
          wt_delta_yaw_angle_l(n) = wt_dir_l(n) - fst_wdir_ave_l(n)
          IF( wt_delta_yaw_angle_l(n) .LT. -180.0 ) THEN 
            wt_delta_yaw_angle_l(n) = wt_dir_l(n) + ( 360.0 - fst_wdir_ave_l(n) )
          ELSEIF(wt_delta_yaw_angle_l(n) .GE. 180.0 ) THEN 
            wt_delta_yaw_angle_l(n) = wt_dir_l(n) - ( 360.0 + fst_wdir_ave_l(n) )
          ENDIF
          CALL calc_vert_area_fact(                   &
              &  wt_mod_hub_height( wt_type_id(n) )   &
              &, wt_rotor_diameter( wt_type_id(n) )/2 &
              &, w_z( wt_i(n), kts:ktf, wt_j(n) )     &
              &, dx, dy, kts, ktf, seg_k_bot          &
              &, seg_k_top, wt_Faz, wt_Rz             &
              &)
          CALL sad_yaw(                               &
              &  seg_k_bot, seg_k_top, wt_i(n)        &
              &, wt_j(n), wt_Faz(seg_k_bot:seg_k_top) &
              &, wt_Rz(seg_k_bot:seg_k_top)           &
              &, wt_dir_l(n), sad_n_points_l(n)       &
              &, sad_i_l(n, :), sad_j_l(n, :)         &
              &, sad_k_l(n, :), sad_area_fact_l(n, :) &
              &)
          ! Turbine not yawing:
          IF( ABS( wt_yaw_error_l(n) ) .LE. wt_yaw_error_threshold( wt_type_id(n) ) ) THEN
            IF( ABS( wt_delta_yaw_angle_l(n) ) .GE. wt_not_yaw_threshold ) THEN
              wt_yaw_error_l(n) = wt_yaw_error_l(n) + ( -SIGN( wt_delta_yaw_angle_l(n)         &
                                , wt_delta_yaw_angle_l(n) ) / ABS( wt_delta_yaw_angle_l(n) ) ) &
                                * wt_delta_yaw_angle_l(n)**2*dt
            ENDIF
          ! Turbine yawing:
          ELSE
            IF( opt_yaw ) THEN
              IF( ABS( wt_delta_yaw_angle_l(n) ) .LT. wt_yaw_angle_threshold ) wt_yaw_error_l(n) = 0.0
              IF( wt_yaw_error_l(n) .NE. 0.0 ) THEN
                ! Check in which direction to yaw:
                wt_yaw_dir(n) = -SIGN( 1.0, wt_delta_yaw_angle_l(n) )
                ! Start turbine yaw:
                wt_delta_yaw_l(n) = wt_yaw_dir(n)*wt_yaw_rate( wt_type_id(n) )*dt
                wt_dir_l(n) = wt_dir_l(n) + wt_delta_yaw_l(n)
                ! Keep angle between 0 and 360 degrees:
                IF( wt_dir_l(n) .GE. 360.0 ) wt_dir_l(n) = wt_dir_l(n) - 360.0
                IF( wt_dir_l(n) .LT.   0.0 ) wt_dir_l(n) = wt_dir_l(n) + 360.0
                CALL sad_yaw(                               &
                    &  seg_k_bot, seg_k_top, wt_i(n)        &
                    &, wt_j(n), wt_Faz(seg_k_bot:seg_k_top) &
                    &, wt_Rz(seg_k_bot:seg_k_top)           &
                    &, wt_dir_l(n), sad_n_points_l(n)       &
                    &, sad_i_l(n, :), sad_j_l(n, :)         &
                    &, sad_k_l(n, :), sad_area_fact_l(n, :) &
                    &)
                wt_do_yaw_l(n) = .TRUE.
              ENDIF
            ENDIF
          ENDIF
        ELSE
          wt_do_yaw_l(n) = .FALSE.
          wt_delta_yaw_l(n) = -999999.999
          wt_yaw_error_l(n) = -999999.999
          wt_dir_l(n) = -999999.999
          sad_n_points_l(n) = -999999
          sad_i_l(n, :) = -999999
          sad_j_l(n, :) = -999999
          sad_k_l(n, :) = -999999
          sad_area_fact_l(n, :) = -999999.999
        ENDIF
        CALL MPI_ALLREDUCE( wt_do_yaw_l(n), wt_do_yaw(n), 1, MPI_LOGICAL, MPI_LOR, mpi_new_comm(n), ierr )
        IF( wt_do_yaw(n) ) THEN
          CALL MPI_ALLREDUCE( wt_delta_yaw_l(n), wt_delta_yaw(n), 1, MPI_REAL, MPI_MAX, mpi_new_comm(n), ierr )
          wt_do_yaw_l(n) = .FALSE.
        ENDIF
        IF( spec_rho .LE. 0.0 ) fst_rho_ave_l(n) = spec_rho
        CALL MPI_ALLREDUCE( wt_dir_l(n), wt_dir(n), 1, MPI_REAL, MPI_MAX, mpi_new_comm(n), ierr )
        CALL MPI_ALLREDUCE( wt_yaw_error_l(n), wt_yaw_error(n), 1, MPI_REAL, MPI_MAX, mpi_new_comm(n), ierr )
        CALL MPI_ALLREDUCE( sad_i_l(n, :), sad_i(n, :), sad_max_points, MPI_INT, MPI_MAX, mpi_new_comm(n), ierr )
        CALL MPI_ALLREDUCE( sad_j_l(n, :), sad_j(n, :), sad_max_points, MPI_INT, MPI_MAX, mpi_new_comm(n), ierr )
        CALL MPI_ALLREDUCE( sad_k_l(n, :), sad_k(n, :), sad_max_points, MPI_INT, MPI_MAX, mpi_new_comm(n), ierr )
        CALL MPI_ALLREDUCE( sad_n_points_l(n), sad_n_points(n), 1, MPI_INT, MPI_MAX, mpi_new_comm(n), ierr )
        CALL MPI_ALLREDUCE( sad_area_fact_l(n, :), sad_area_fact(n, :), sad_max_points, MPI_REAL, MPI_MAX, mpi_new_comm(n), ierr )
        ! wt_delta_yaw(n) is set to 0.0 initially. If wind turbine is not yawing, wt_do_yaw(n) = .FALSE., initial
        ! value of wt_delta_yaw(n) is preserved.
        IF( opt_tsout_plane .AND. ANY(spec_wt .EQ. n) ) THEN
          CALL calc_wind_plane_tsout_vel(        &
              &  ims, ime, jms, jme, kms, kme    &
              &, its, itf, jts, jtf, kts, ktf    &
              &, n_pln_wt, n_pln_wt_counter      &
              &, n_vert_points, n_horz_points    &
              &, dx, dy, dz_cst( wt_type_id(n) ) &
              &, dz(its:itf, kts:ktf, jts:jtf)   &
              &, wt_do_yaw(n), plane_points      &
              &, dist_wt_base_to_hub(n, :)       &
              &, wt_unit_vector(n, :)            &
              &, wt_delta_yaw(n)                 &
              &, c_z(its:itf, kts:ktf, jts:jtf)  &
              &, w_z(its:itf, kts:ktf, jts:jtf)  &
              &, u_2(ims:ime, kms:kme, jms:jme)  &
              &, v_2(ims:ime, kms:kme, jms:jme)  &
              &, w_2(ims:ime, kms:kme, jms:jme)  &
              &, mpi_new_comm(n), ierr           &
              &, u_plane_l, v_plane_l, w_plane_l &
              &)
        ENDIF
        IF( fst_uinf_l(n) .GE. wt_cutin( wt_type_id(n) ) .AND. fst_uinf_l(n) .LE. wt_cutout( wt_type_id(n) ) ) THEN
          ils = wt_i(n) - wt_norm_radius( wt_type_id(n) )
          ile = wt_i(n) + wt_norm_radius( wt_type_id(n) )
          jls = wt_j(n) - wt_norm_radius( wt_type_id(n) )
          jle = wt_j(n) + wt_norm_radius( wt_type_id(n) )
          ! If wind turbine is within the current tile:
          IF(ils .LE. ite .AND. ile .GE. its .AND. jls .LE. jte .AND. jle .GE. jts) THEN
            wt_aif_l(n) = calc_wt_aif(wt_type_id(n), wt_vel_hub(n) )
            IF( opt_aif .EQ. 1 ) THEN
              wt_aif_tmp = ( fst_uinf_l(n)-wt_vel_hub(n) )/fst_uinf_l(n)
              IF( wt_aif_tmp .LE. wt_aif_l(n) ) wt_aif_l(n) = wt_aif_tmp
              wt_ct_l(n) = calc_wt_ct(wt_type_id(n), fst_uinf_l(n), 1)
            ELSEIF( opt_aif .EQ. 2 ) THEN
              ! freestream wind velocity is calculated from tabulated thrust coefficient and hub height
              ! horizontal wind velocity at the rotor disk. This means freestream wind velocity obtained
              ! from LES flow field is overwritten (it does not seem to be accurate for me; thus, I
              ! recommend using direct calculation of axial induction factor from LES flow field,
              ! opt_aif = 1).
              fst_uinf_l(n) = wt_vel_hub(n)/( 1-wt_aif_l(n) )
              wt_ct_l(n) = calc_wt_ct(wt_type_id(n), wt_vel_hub(n), 2)
            ENDIF
            wt_power_l(n) = calc_wt_power(wt_type_id(n), fst_uinf_l(n) )*1000.0 ! convert to Watt
            wt_cp_l(n) = wt_power_l(n)/(0.5*fst_rho_ave_l(n)*wt_area( wt_type_id(n) )*fst_uinf_l(n)**3)
            wt_ck_l(n) = wt_ck_cst*( wt_ct_l(n)-wt_cp_l(n) )
            IF(wt_ck_l(n) .LT. 0.0) wt_ck_l(n) = 0.0
            ! Now update tendency for u,v, and tke:
            DO m = 1, sad_n_points(n)
              ii = sad_i(n, m)
              kk = sad_k(n, m)
              jj = sad_j(n, m)
              ! Update tendencies within the current tile only:
              IF(ii .LE. ite .AND. ii .GE. its .AND. jj .LE. jte .AND. jj .GE. jts) THEN
                grid_ijk_wsph = SQRT( u_2(ii, kk, jj)**2 + v_2(ii, kk, jj)**2 )
                ru_tendf(ii, kk, jj) = ru_tendf(ii, kk, jj) - 0.5*muu(ii, jj)*wt_ct_l(n)*grid_ijk_wsph*u_2(ii, kk, jj)*sad_area_fact(n, m)/( 1-wt_aif_l(n) )**2 
                rv_tendf(ii, kk, jj) = rv_tendf(ii, kk, jj) - 0.5*muv(ii, jj)*wt_ct_l(n)*grid_ijk_wsph*v_2(ii, kk, jj)*sad_area_fact(n, m)/( 1-wt_aif_l(n) )**2
                tke_tend(ii, kk, jj) = tke_tend(ii, kk, jj) + 0.5*mut(ii, jj)*wt_ck_l(n)*(grid_ijk_wsph**3)*sad_area_fact(n, m)/(1-wt_aif_l(n) )**3
              ENDIF
            ENDDO
          ELSE
            wt_aif_l(n) = -999999.999
            wt_ct_l(n) = -999999.999
            wt_cp_l(n) = -999999.999
            wt_ck_l(n) = -999999.999
            wt_power_l(n) = -999999.999
          ENDIF
          CALL MPI_ALLREDUCE( wt_aif_l(n), wt_aif(n), 1 ,MPI_REAL, MPI_MAX, mpi_new_comm(n), ierr )
          CALL MPI_ALLREDUCE( wt_ct_l(n), wt_ct(n), 1 ,MPI_REAL, MPI_MAX, mpi_new_comm(n), ierr )
          CALL MPI_ALLREDUCE( wt_cp_l(n), wt_cp(n), 1 ,MPI_REAL, MPI_MAX, mpi_new_comm(n), ierr )
          CALL MPI_ALLREDUCE( wt_power_l(n), wt_power(n), 1, MPI_REAL, MPI_MAX, mpi_new_comm(n), ierr )
        ELSE
          IF( fst_uinf_l(n) .LT. wt_cutin( wt_type_id(n) ) ) THEN
            WRITE(outstring, "(A, I4, A, F8.3, A)") &
                             'WTP-WARNING ***** Wind turbine # ', n, ' stopped due to &
                              very low instantaneous wind velocity component normal to &
                              the rotor (fst_uinf = ', fst_uinf_l(n), ' [m/s] < cut-in speed). *****'
            CALL wrf_message ( TRIM(outstring) )
          ELSEIF( fst_uinf_l(n) .GT. wt_cutout( wt_type_id(n) ) ) THEN
            WRITE(outstring, "(A, I4, A, F8.3, A)") &
                             'WTP-WARNING ***** Wind turbine # ', n, ' stopped due to &
                              very high instantaneous wind velocity component normal to &
                              the rotor (fst_uinf = ', fst_uinf_l(n), ' [m/s] > cut-out speed). *****'
            CALL wrf_message ( TRIM(outstring) )
          ENDIF
          ! If freestream wind speed is below cut-in speed or above cut-out speed of wind turbine
          ! set the variables to -999999.999 for proper reduction to all processors.
          wt_ct_l(n) = -999999.999
          wt_cp_l(n) = -999999.999
          wt_power_l(n) = -999999.999
        ENDIF
      ELSE
        ! These variables are set to a very small number for proper reduction to all processors
        ! when writing into wrfout or wrfrst files.
        fst_i_l(n) = -999999
        fst_j_l(n) = -999999
        wt_yaw_error_l(n) = -999999.999
        fst_rho_ave_l(n) = -999999.999
        fst_wdir_ave_l(n) = -999999.999
        wt_aif_l(n) = -999999.999
        wt_ct_l(n) = -999999.999
        wt_cp_l(n) = -999999.999
        wt_power_l(n) = -999999.999
        wt_dir_l(n) = -999999.999
        fst_uinf_l(n) = -999999.999
        fst_rho_ts_l(n, :) = -999999.999
        fst_u_ts_l(n, :) = -999999.999
        fst_v_ts_l(n, :) = -999999.999
        fst_wdir_ts_l(n, :) = -999999.999
        IF( opt_tsout_plane .AND. ANY(spec_wt .EQ. n) ) THEN
          k = n_pln_wt_counter
          u_plane_l(k, :, :) = -999999.999
          v_plane_l(k, :, :) = -999999.999
          w_plane_l(k, :, :) = -999999.999
          plane_points(k, :, :, :) = -999999.999
        ENDIF
      ENDIF
      IF( opt_tsout_plane .AND. ANY(spec_wt .EQ. n) ) THEN
        k = n_pln_wt_counter
        CALL MPI_REDUCE( u_plane_l(k, :, :), u_plane(k, :, :), n_horz_points*n_vert_points, MPI_REAL, MPI_MAX, mpi_proc_root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( v_plane_l(k, :, :), v_plane(k, :, :), n_horz_points*n_vert_points, MPI_REAL, MPI_MAX, mpi_proc_root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( w_plane_l(k, :, :), w_plane(k, :, :), n_horz_points*n_vert_points, MPI_REAL, MPI_MAX, mpi_proc_root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( plane_points(k, :, :, 1), plane_points_x(k, :, :), n_horz_points*n_vert_points, MPI_REAL, MPI_MAX, mpi_proc_root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( plane_points(k, :, :, 2), plane_points_y(k, :, :), n_horz_points*n_vert_points, MPI_REAL, MPI_MAX, mpi_proc_root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( plane_points(k, :, :, 3), plane_points_z(k, :, :), n_horz_points*n_vert_points, MPI_REAL, MPI_MAX, mpi_proc_root, MPI_COMM_WORLD, ierr )
        n_pln_wt_counter = n_pln_wt_counter + 1
      ENDIF
      IF( opt_tsout ) THEN
        CALL MPI_REDUCE( fst_rho_ave_l(n), fst_rho_ave(n), 1, MPI_REAL, MPI_MAX, mpi_proc_root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( fst_wdir_ave_l(n), fst_wdir_ave(n), 1, MPI_REAL, MPI_MAX, mpi_proc_root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( wt_power_l(n), wt_power(n), 1, MPI_REAL, MPI_MAX, mpi_proc_root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( wt_aif_l(n), wt_aif(n), 1, MPI_REAL, MPI_MAX, mpi_proc_root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( wt_ct_l(n), wt_ct(n), 1, MPI_REAL, MPI_MAX, mpi_proc_root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( wt_cp_l(n), wt_cp(n), 1, MPI_REAL, MPI_MAX, mpi_proc_root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( wt_dir_l(n), wt_dir(n), 1, MPI_REAL, MPI_MAX, mpi_proc_root, MPI_COMM_WORLD, ierr )
        CALL MPI_REDUCE( fst_uinf_l(n), fst_uinf(n), 1, MPI_REAL, MPI_MAX, mpi_proc_root, MPI_COMM_WORLD, ierr )
      ENDIF
      IF( wtp_history_flag ) THEN
        CALL MPI_ALLREDUCE( fst_rho_ave_l(n), fst_rho_ave(n), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( fst_wdir_ave_l(n), fst_wdir_ave(n), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( wt_power_l(n), wt_power(n), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( wt_aif_l(n), wt_aif(n), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( wt_ct_l(n), wt_ct(n), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( wt_cp_l(n), wt_cp(n), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( wt_dir_l(n), wt_dir(n), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( fst_uinf_l(n), fst_uinf(n), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
      ENDIF
      IF( wtp_restart_flag ) THEN
        CALL MPI_ALLREDUCE( fst_i_l(n), fst_i(n), 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( fst_j_l(n), fst_j(n), 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( fst_rho_ave_l(n), fst_rho_ave(n), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( wt_dir_l(n), wt_dir(n), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( wt_yaw_error_l(n), wt_yaw_error(n), 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( fst_rho_ts_l(n, :), fst_rho_ts(n, :), n_ave_samples, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( fst_u_ts_l(n, :), fst_u_ts(n, :), n_ave_samples, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( fst_v_ts_l(n, :), fst_v_ts(n, :), n_ave_samples, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
        CALL MPI_ALLREDUCE( fst_wdir_ts_l(n, :), fst_wdir_ts(n, :), n_ave_samples, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
      ENDIF
#endif
      IF( opt_print_diag .EQ. 1 )THEN
        IF( wrf_dm_on_monitor() ) THEN
          IF( fst_uinf_l(n) .GE. wt_cutin( wt_type_id(n) ) .AND. fst_uinf_l(n) .LE. wt_cutout( wt_type_id(n) ) ) THEN
            WRITE(outstring, "(A, I4, A, A, F8.3, A, A, F7.3, A, A, F8.3, A, &
                               A, F7.3, A, A, F9.3, A, A, F5.3, A, A, F6.3,  &
                               A, A, F5.3, A, A, F10.3)")                    &
                             ' WT# ', n, ',',                                &
                             ' Wdirinf = ', fst_wdir_ave_l(n), ' [deg],',    &
                             ' Uinf = ', fst_uinf_l(n), ' [m/s],',           &
                             ' WTdir = ', wt_dir(n), ' [deg],',              &
                             ' Uhub = ', wt_vel_hub(n), ' [m/s],',           &
                             ' Power = ', wt_power(n)/1000.0, ' [kW],',      &
                             ' CP = ', wt_cp(n), ',',                        &
                             ' CT = ', wt_ct(n), ',',                        &
                             ' aif = ', wt_aif(n), ',',                      &
                             ' YawErr = ', wt_yaw_error(n)
            CALL wrf_message ( TRIM(outstring) )
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    rst_controller = .FALSE.
    read_from_rst_file = 1
  END SUBROUTINE calc_power_ct_sad
!=======================================================================
! This subroutine interpolates velocity components u, v, w to desired
! i, j, k levels.
!=======================================================================
  SUBROUTINE interp_scalar(            &
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
  END SUBROUTINE interp_scalar
!=======================================================================
! This subroutine interpolates velocity components u, v, w to desired
! i, j, k levels.
!=======================================================================
  SUBROUTINE interp_velocity(          &
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
  END SUBROUTINE interp_velocity
!=======================================================================
  FUNCTION vector_add(a, b)
! This function adds 2 vectors ( arrays real, dimension(3) )
!=======================================================================
  IMPLICIT NONE
    REAL, DIMENSION(3), INTENT(IN) :: a,b
    REAL, DIMENSION(3) :: vector_add
    vector_add(1) = a(1) + b(1)
    vector_add(2) = a(2) + b(2)
    vector_add(3) = a(3) + b(3)
    RETURN
  END FUNCTION vector_add
!=======================================================================
  FUNCTION vector_divide(a, b)
! This function divides one vector ( array real, dimension(3) by a
! number )
!=======================================================================
  IMPLICIT NONE
    REAL, DIMENSION(3), INTENT(IN) :: a
    REAL, INTENT(IN) :: b
    REAL, DIMENSION(3) :: vector_divide
    vector_divide(1) = a(1)/b
    vector_divide(2) = a(2)/b
    vector_divide(3) = a(3)/b
    RETURN
  END FUNCTION vector_divide
!=======================================================================
  FUNCTION vector_mag(a)
! This function calculates the magnitude of a vector
!=======================================================================
  IMPLICIT NONE
    REAL, DIMENSION(3), INTENT(IN) :: a
    REAL :: vector_mag
    vector_mag = ABS( SQRT( a(1)**2 + a(2)**2 + a(3)**2 ) )
    RETURN
  END FUNCTION vector_mag
!=======================================================================
  FUNCTION rotate_point(point_in, rotationPoint, axis, angle)
! This function performs rotation of a point with respect to an axis or
! rotation and a certain angle
!=======================================================================
  IMPLICIT NONE
    REAL, DIMENSION(3), INTENT(IN) :: point_in
    REAL, DIMENSION(3), INTENT(IN) :: rotationPoint
    REAL, DIMENSION(3), INTENT(IN) :: axis
    REAL, INTENT(IN) :: angle
    REAL, DIMENSION(3, 3) :: RM ! Rotation Matrix tensor
    REAL, DIMENSION(3) :: rotate_point, point
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
    point = vector_add(point, -rotationPoint)
    ! Perform rotation (multiplication matrix and vector)
    point = matrix_vector(RM, point)
    ! Return the rotated point to its new location relative to the rotation point
    rotate_point = point + rotationPoint
    RETURN
  END FUNCTION rotate_point
!=======================================================================
  FUNCTION matrix_vector(RM, point)
! This function multiplies a matrix and a vector
!=======================================================================
  IMPLICIT NONE
    REAL, DIMENSION(3, 3), INTENT(IN) :: RM ! Matrix
    REAL, DIMENSION(3), INTENT(IN) :: point ! vector point
    REAL, DIMENSION(3) :: matrix_vector
    ! Perform rotation
    matrix_vector(1) = RM(1, 1)*point(1) + RM(1, 2)*point(2) + RM(1, 3)*point(3)
    matrix_vector(2) = RM(2, 1)*point(1) + RM(2, 2)*point(2) + RM(2, 3)*point(3)
    matrix_vector(3) = RM(3, 1)*point(1) + RM(3, 2)*point(2) + RM(3, 3)*point(3)
    RETURN
  END FUNCTION matrix_vector
!=======================================================================
  SUBROUTINE compass_to_standard(dir)
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
  END SUBROUTINE compass_to_standard
!=======================================================================
! This subroutine is used to perform global to local rotation.
!=======================================================================
  SUBROUTINE rot_global_to_local(phi, global_variable_in, local_variable_out)
  IMPLICIT NONE
    REAL, INTENT(IN) :: phi
    REAL, DIMENSION(3), INTENT(IN) :: global_variable_in
    REAL, DIMENSION(3), INTENT(OUT) :: local_variable_out
! local
    REAL, DIMENSION(3, 3) :: transpose_yaw
    ! Negative rotation around z-axis
    transpose_yaw(1, 1:3) = (/   COS(phi*degrad), SIN(phi*degrad), 0.0 /)
    transpose_yaw(2, 1:3) = (/ - SIN(phi*degrad), COS(phi*degrad), 0.0 /)
    transpose_yaw(3, 1:3) = (/               0.0,             0.0, 1.0 /)
    local_variable_out = MATMUL(transpose_yaw, global_variable_in)
  END SUBROUTINE rot_global_to_local
!=======================================================================
! This subroutine calculates distance between vertical levels.
!=======================================================================
! wt_mod_hub_height, wt_rotor_radius are wind turbine hub height and
! radius  are height of full-level, which is assume to change very
! little, we can recalculate when yawing as well areas are between the
! vertical levels. Rs are the radius of the approximated segment at each
! half level; then normalized by dx.
  SUBROUTINE  calc_vert_area_fact(             &
        &  wt_mod_hub_height, wt_rotor_radius  &
        &, w_z, dx, dy, kts, ktf, k_bot, k_top &
        &, Fas, Rs                             &
        &)
  IMPLICIT NONE
    REAL, INTENT(IN) :: wt_mod_hub_height, wt_rotor_radius, dx, dy
    INTEGER, INTENT(IN) :: kts, ktf
    REAL, DIMENSION(kts:ktf),INTENT(IN) :: w_z
    INTEGER, INTENT(OUT) :: k_bot, k_top
    REAL, DIMENSION(kts:ktf-1), INTENT(OUT) :: Fas, Rs
! local 
    REAL :: wt_z_bot, wt_z_top, zr1, zr2, area1, area2
    REAL, DIMENSION(kts:ktf-1) :: z1s, z2s, areas
    INTEGER :: k     
    Fas = 0.0
    Rs = 0.0
    wt_z_bot = wt_mod_hub_height-wt_rotor_radius
    wt_z_top = wt_mod_hub_height+wt_rotor_radius
    k_bot = -1
    k_top = -1
    DO k = 1, ktf
      IF( (w_z(k) .GE. wt_z_bot) .AND. (k_bot .EQ. -1) ) k_bot = k
      IF( (w_z(k) .GE. wt_z_top) .AND. (k_top .EQ. -1) ) THEN
        k_top = k-1
        EXIT
      ENDIF
    ENDDO
    z1s = 0.0
    z2s = 0.0
     DO k = k_bot, k_top
        z1s(k) = w_z(k)
        z2s(k) = w_z(k+1)
        IF(k .EQ. k_bot) THEN
          z1s(k) = wt_z_bot
        ELSEIF (k .EQ. k_top) THEN
          z2s(k) = wt_z_top
        ENDIF
        area1 = circle_segment(wt_rotor_radius, z1s(k)-wt_mod_hub_height)
        area2 = circle_segment(wt_rotor_radius, z2s(k)-wt_mod_hub_height)
        areas(k) = area1-area2
        Fas(k) = areas(k)/( w_z(k+1)-w_z(k) )/(dx*dy)
        Rs(k) = areas(k)/( z2s(k)-z1s(k) )/2
     ENDDO 
     Rs = Rs/dx  !Normalized by dx
  END SUBROUTINE calc_vert_area_fact
!=======================================================================
! This subroutine updates wind turbine shaft direction with respect to
! incoming wind direction (i.e., wind turbine yaw motion).
!=======================================================================
  SUBROUTINE sad_yaw(                 &
        &  k_bot, k_top, wt_i, wt_j   &
        &, wt_Faz, wt_Rz, wt_dir      &
        &, sad_n_points, sad_i, sad_j &
        &, sad_k, sad_area_fact       &
        &)
  IMPLICIT NONE
    INTEGER, INTENT(IN) :: k_bot, k_top, wt_i, wt_j
    REAL, INTENT(IN) :: wt_dir
    REAL, DIMENSION(k_bot:k_top), INTENT(IN) :: wt_Faz, wt_Rz
    REAL, DIMENSION(1:sad_max_points), INTENT(OUT) :: sad_area_fact
    INTEGER, INTENT(OUT) :: sad_n_points
    INTEGER, DIMENSION(1:sad_max_points), INTENT(OUT) :: sad_i, sad_j, sad_k
! local
    INTEGER :: k, m, row_n_points
    INTEGER, DIMENSION(1:row_max_points) :: row_i, row_j
    REAL, DIMENSION(1:row_max_points) :: row_area_fact
    sad_n_points = 0
    sad_i = -999999
    sad_j = -999999
    sad_k = -999999
    sad_area_fact = -999999.999
    !Now calculate the 3D grid-based area factor
    DO k = k_bot, k_top
      CALL calc_horz_area_split(                    &
          &  wt_Faz(k)                              &
          &, REAL(wt_i)-wt_Rz(k)*SIN(wt_dir*degrad) &
          &, REAL(wt_i)+wt_Rz(k)*SIN(wt_dir*degrad) &
          &, REAL(wt_j)-wt_Rz(k)*COS(wt_dir*degrad) &
          &, REAL(wt_j)+wt_Rz(k)*COS(wt_dir*degrad) &
          &, row_max_points, row_n_points           &
          &, row_i, row_j, row_area_fact            &
          &)
      DO m = 1, row_n_points
        sad_n_points = sad_n_points+1
        IF(sad_n_points .GT. sad_max_points) THEN
          CALL wrf_error_fatal("***WTP-ERROR*** sad_n_points is too small. &
                                sad_n_points = %d", sad_max_points)
!          WRITE(0,*) "***WTP-ERROR*** sad_max_points is too small"
!          WRITE(0,*) "sad_max_points=", sad_max_points
!          STOP
        ENDIF
        sad_i(sad_n_points) = row_i(m)
        sad_j(sad_n_points) = row_j(m)
        sad_k(sad_n_points) = k
        sad_area_fact(sad_n_points) = row_area_fact(m)
      ENDDO 
    ENDDO
  END SUBROUTINE sad_yaw
!=======================================================================
! This subroutine discretizes vertical area into horizontal grid and
! returns the number of grid points (n_points), their associate i,j, and
! the splitted area.
!=======================================================================
  SUBROUTINE calc_horz_area_split(            &
        &  wt_Faz, wt_xs, wt_xe, wt_ys, wt_ye &
        &, row_max_points, row_n_points       &
        &, grid_i, grid_j, grid_area_fact     &
        &)
  IMPLICIT NONE
    REAL,INTENT(IN) :: wt_Faz, wt_xs, wt_xe, wt_ys, wt_ye
    INTEGER,INTENT(IN) :: row_max_points
    INTEGER,INTENT(OUT) :: row_n_points
    INTEGER,DIMENSION(1:row_max_points), INTENT(OUT) :: grid_i, grid_j
    REAL, DIMENSION(1:row_max_points), INTENT(OUT) :: grid_area_fact
! local
    INTEGER, PARAMETER :: n_seq = 50
    REAL, DIMENSION(1:n_seq) :: wt_i, wt_j
    INTEGER :: n, i, j, curr_i, curr_j, curr_n
    CALL spacing_linear(wt_xs, wt_xe, n_seq, wt_i)
    CALL spacing_linear(wt_ys, wt_ye, n_seq, wt_j)
    curr_i = -999999
    curr_j = -999999
    row_n_points = 0
    grid_area_fact = 0.0
    DO n = 1, n_seq
      i = INT( wt_i(n) )
      j = INT( wt_j(n) )
      IF(i .NE. curr_i .OR. j .NE. curr_j) THEN  ! a new cell
        row_n_points = row_n_points+1
        IF(row_n_points .GT. row_max_points) THEN
          CALL wrf_error_fatal("***WTP-ERROR*** row_max_points is too small. &
                                row_max_points = %d", row_max_points)
!          WRITE(0,*) "***WTP-ERROR*** row_max_points is too small"
!          WRITE(0,*) "row_max_points=", row_max_points
!          STOP
        ENDIF
        grid_i(row_n_points) = i
        grid_j(row_n_points) = j
        curr_i = i
        curr_j = j
      ENDIF
      grid_area_fact(row_n_points) = grid_area_fact(row_n_points) + (wt_Faz/n_seq)
    ENDDO
  END SUBROUTINE calc_horz_area_split
!=======================================================================
! This function calculates the circle segment area for z to the top of
! the cirle (z is relative to the circle origin).
!=======================================================================
  REAL FUNCTION circle_segment(r, z)
    IMPLICIT NONE
    REAL, INTENT(IN) :: r, z
! local
    REAL :: area
    REAL :: r1, angle
    IF(ABS(z) .GT. r) THEN
      area = 0.0
    ELSE
      r1 = SQRT(r**2-z**2)
      angle = ACOS(z/r)
      area = angle*r**2-z*r1
    END IF
    circle_segment = area
    RETURN 
  END FUNCTION circle_segment
!=======================================================================
! Simple linear spcaing.
!=======================================================================
 !n_point>=2
  SUBROUTINE spacing_linear(x_start, x_end, n_point, arr)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x_start, x_end
    INTEGER, INTENT(IN) :: n_point
    REAL, DIMENSION(1:n_point), INTENT(OUT) :: arr
! local
    INTEGER :: i
    REAL :: dx
    dx = (x_end-x_start)/(n_point-1)
    arr = (/(x_start+(i-1)*dx, i=1, n_point)/)
  END SUBROUTINE spacing_linear
!=======================================================================
! Simple linear interpolation.
!=======================================================================
  REAL FUNCTION interpolate_linear(y, x, n_x, x_0)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n_x
    REAL, DIMENSION(n_x), INTENT(IN) :: x, y
    REAL, INTENT(IN) :: x_0
    INTEGER :: i1, i2
    IF(x_0 .LE. x(1) ) THEN
      !interpolate_linear = y(1) 
      interpolate_linear = 0 ! cut-in speed
    ELSEIF(x_0 .GE. x(n_x) ) THEN
      !interpolate_linear = y(n_x)
      interpolate_linear = 0 ! cut out speed
    ELSE
      !Simple search, would be fast if using, for example binary search.
      DO i2 = 2, n_x
        IF(x(i2) .GE. x_0) THEN
          i1 = i2-1
          EXIT
        ENDIF
      ENDDO
      interpolate_linear = y(i1)+( x_0-x(i1) )*( y(i2)-y(i1) )/( x(i2)-x(i1) )
    ENDIF
    RETURN
  END FUNCTION interpolate_linear
!=======================================================================
! This function calculates wind turbine thrust coeffient from tabulated
! data.
!=======================================================================
  REAL FUNCTION calc_wt_ct(wt_type_id, wt_wsph, opt_aif)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: wt_type_id, opt_aif
    REAL, INTENT(IN) :: wt_wsph
! local
    REAL, DIMENSION(1:n_max_val) :: tab_wsph, tab_ct
    INTEGER :: n_lines
    tab_ct = wt_tab_ct(wt_type_id, :)
    IF(opt_aif .eq. 1) THEN
      tab_wsph = wt_tab_wsph(wt_type_id, :)  
    ELSEIF(opt_aif .eq. 2) THEN
      tab_wsph = wt_tab_wsph_1(wt_type_id, :)  
    ENDIF
    n_lines = wt_tab_n_lines(wt_type_id)
    ! We need to interpolate using these reference speed and coefficent:
    calc_wt_ct = interpolate_linear(tab_ct, tab_wsph, n_lines, wt_wsph)   
  END FUNCTION calc_wt_ct
!=======================================================================
! This function calculates wind turbine thrust coeffient from tabulated
! data (opt_aif=2 only).
!=======================================================================
  REAL FUNCTION calc_wt_aif(wt_type_id, wt_wsph)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: wt_type_id
    REAL, INTENT(IN) :: wt_wsph
! local
    REAL, DIMENSION(1:n_max_val) :: tab_wsph, tab_aif
    INTEGER :: n_lines
    tab_aif = wt_tab_aif(wt_type_id, :)
    ! Use hub height horizontal wind speed at rotor disk:
    tab_wsph = wt_tab_wsph_1(wt_type_id, :)
    n_lines = wt_tab_n_lines(wt_type_id)
    ! We need to interpolate using these reference speed and coefficent:   
    calc_wt_aif = interpolate_linear(tab_aif, tab_wsph, n_lines, wt_wsph)
  END FUNCTION calc_wt_aif
!=======================================================================
! This function calculates wind turbine power.
!=======================================================================
 REAL FUNCTION calc_wt_power(wt_type_id, wt_wsph)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: wt_type_id
    REAL, INTENT(IN) :: wt_wsph
! local
    REAL, DIMENSION(1:n_max_val) :: tab_wsph, tab_power
    INTEGER :: n_lines
    tab_power = wt_tab_power(wt_type_id, :)
    ! Use hub height horizontal wind speed at freestream:
    tab_wsph = wt_tab_wsph(wt_type_id, :)
    n_lines = wt_tab_n_lines(wt_type_id)
    ! We need to interpolate using these reference speed and coefficent:
    calc_wt_power = interpolate_linear(tab_power, tab_wsph, n_lines, wt_wsph)     
  END FUNCTION calc_wt_power
!=======================================================================
! This subroutine is used to initialize Simplified Actuator Disk module.
!=======================================================================
  SUBROUTINE init_module_wind_sad(        &
        &  grid, config_flags, restart    &
        &, xlong_u, xlat_v, ht, dx, dy    &
        &, ids, ide, jds, jde, kds, kde   &
        &, ims, ime, jms, jme, kms, kme   &
        &, its, ite, jts, jte, kts, kte   &
#ifdef DM_PARALLEL
        &, n_proc_x                       &                         
#endif
        &, dz8w, u_2, v_2, wt_dir         &
        &, plane_points_x, plane_points_y &
        &, plane_points_z &
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
    INTEGER, INTENT(IN) :: n_proc_x
#endif
!
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: xlong_u, xlat_v, ht
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: dz8w, u_2, v_2
!
    REAL, DIMENSION(1:config_flags%wind_wtp_pln_nTurbines), INTENT(IN) :: wt_dir
!
    REAL, DIMENSION(1:config_flags%wind_wtp_pln_nTurbines, &
                    1:config_flags%wind_wtp_pln_nz,        &
                    1:config_flags%wind_wtp_pln_nxy        ), INTENT(IN) :: plane_points_x
    REAL, DIMENSION(1:config_flags%wind_wtp_pln_nTurbines, &
                    1:config_flags%wind_wtp_pln_nz,        &
                    1:config_flags%wind_wtp_pln_nxy        ), INTENT(IN) :: plane_points_y
    REAL, DIMENSION(1:config_flags%wind_wtp_pln_nTurbines, &
                    1:config_flags%wind_wtp_pln_nz,        &
                    1:config_flags%wind_wtp_pln_nxy        ), INTENT(IN) :: plane_points_z
!
    LOGICAL, INTENT(IN) :: restart
!
    LOGICAL, EXTERNAL :: wrf_dm_on_monitor
! Local
    REAL :: ts_rx, ts_ry, known_lat, known_lon
    REAL :: dist_wt_to_plane
!
    INTEGER :: n_val, i, j, k, m, n, p, n_size
    INTEGER :: itf, jtf, ktf
!
    CHARACTER*256 :: input, wt_num
    CHARACTER*256 :: message_windturbine, outstring
!
    REAL, DIMENSION(:), ALLOCATABLE :: wt_dir_init, wt_yaw_angle
    REAL, DIMENSION(:), ALLOCATABLE :: wt_unit_vector_shaft_dir
    REAL, DIMENSION(:), ALLOCATABLE :: wt_min_x, wt_max_x, wt_min_y, wt_max_y
    REAL, DIMENSION(:), ALLOCATABLE :: wt_tab_wsph_tmp, wt_tab_power_tmp, wt_tab_ct_tmp
#ifdef DM_PARALLEL
    INTEGER, DIMENSION(:), ALLOCATABLE :: mpi_procs_comm_tmp, mpi_procs_tmp
!
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: mpi_procs, mpi_procs_comm
#endif
    INTEGER, DIMENSION(:), ALLOCATABLE :: arrayIN, arrayOUT
    INTEGER, DIMENSION(:), ALLOCATABLE :: fst_i_init, fst_j_init
!
    REAL, DIMENSION(its-1:ite+1, kts:kte, jts-1:jte+1) :: w_z ! z-coordinate of w-nodes (m) ! Grid point (node)
    REAL, DIMENSION(its-1:ite+1, kts:kte, jts-1:jte+1) :: dz  ! Thickness of each grid cell (m)
!=======================================================================
    IF( restart ) THEN
      rst_controller = .TRUE.
      rst_controller1 = .TRUE.
    ENDIF       
    ! namelist.input options for the WTP model:
    opt_yaw = config_flags%wind_wtp_yaw_opt
    opt_print_diag = config_flags%wind_wtp_printdiag_opt
    opt_tsout = config_flags%wind_wtp_sad_tsout_opt
    opt_aif = config_flags%wind_wtp_aif_opt
    n_wt = config_flags%wind_wtp_nTurbines
    n_ave_samples = config_flags%wind_wtp_nAvesamples
    spec_rho = config_flags%wind_wtp_spec_rho
    wt_ck_cst = config_flags%wind_wtp_tke_fact
!    max_radius = config_flags%wind_wtp_max_radius
!    min_dx = config_flags%wind_wtp_min_dx
!    min_dz = config_flags%wind_wtp_min_dz
!
    opt_tsout_plane = config_flags%wind_wtp_pln_tsout_opt
    n_pln_wt = config_flags%wind_wtp_pln_nTurbines
    n_horz_points = config_flags%wind_wtp_pln_nxy
    n_vert_points = config_flags%wind_wtp_pln_nz
    plane_width = config_flags%wind_wtp_pln_width
    plane_height = config_flags%wind_wtp_pln_height
    plane_distance = config_flags%wind_wtp_pln_dist
    IF( wrf_dm_on_monitor() ) THEN
! ... Open the file with the location of the wind turbines ...
      IF( config_flags%windfarm_ij .EQ. 1 ) THEN
        OPEN(70, FILE='windturbines-ij.txt', FORM='FORMATTED', STATUS='OLD')
      ELSE
        OPEN(70, FILE='windturbines.txt', FORM='FORMATTED', STATUS='OLD')
      ENDIF
 100  CONTINUE
      REWIND (70)
    ENDIF
! ... When doing spin-up runs in real configuration, real.exe sometimes initializes the finest domain
! ... with SAD twice. For that reason, first, deallocate all allocatable arrays to work around the
! ... allocated array issue.
    IF( ALLOCATED(spec_wt) ) DEALLOCATE(spec_wt) ! wind turbine IDs [-]
    IF( ALLOCATED(wt_x) ) DEALLOCATE(wt_x) ! wind turbine location in x-direction [m]
    IF( ALLOCATED(wt_y) ) DEALLOCATE(wt_y) ! wind turbine location in y-direction [m]
    IF( ALLOCATED(wt_z) ) DEALLOCATE(wt_z) ! wind turbine location in z-direction [m]
    IF( ALLOCATED(wt_lat) ) DEALLOCATE(wt_lat) ! wind turbine location in latitudes [deg]
    IF( ALLOCATED(wt_lon) ) DEALLOCATE(wt_lon) ! wind turbine location in longitudes [deg]
    IF( ALLOCATED(wt_type_id) ) DEALLOCATE(wt_type_id) ! wind turbine type id [-]
    IF( ALLOCATED(wt_i) ) DEALLOCATE(wt_i) ! wind turbine i-index [-]
    IF( ALLOCATED(wt_j) ) DEALLOCATE(wt_j) ! wind turbine j-index [-]
    IF( ALLOCATED(wt_tab_n_lines) ) DEALLOCATE(wt_tab_n_lines) ! tabulated data total number of lines [-]
    IF( ALLOCATED(wt_cutin) ) DEALLOCATE(wt_cutin) ! wind turbine cut-in wind speed [m/s]
    IF( ALLOCATED(wt_cutout) ) DEALLOCATE(wt_cutout) ! wind turbine cut-out wind speed [m/s]
    IF( ALLOCATED(wt_hub_height) ) DEALLOCATE(wt_hub_height) ! wind turbine hub height [m]
    IF( ALLOCATED(wt_rotor_diameter) ) DEALLOCATE(wt_rotor_diameter) ! wind turbine rotor diameter [m]
    IF( ALLOCATED(wt_std_ct) ) DEALLOCATE(wt_std_ct) ! wind turbine standing thrust coefficient [-]
    IF( ALLOCATED(wt_rated_power) ) DEALLOCATE(wt_rated_power) ! wind turbine rated power [MW]
    IF( ALLOCATED(wt_norm_radius) ) DEALLOCATE(wt_norm_radius) ! wind turbine normalized radius (R/dx) [-]
    IF( ALLOCATED(wt_area) ) DEALLOCATE(wt_area) ! wind turbine normalized area [-]
    IF( ALLOCATED(wt_tab_wsph) ) DEALLOCATE(wt_tab_wsph) ! tabulated freestream horizontal wind speed at hub height [m/s]
    IF( ALLOCATED(wt_tab_power) ) DEALLOCATE(wt_tab_power) ! tabulated power output [W]
    IF( ALLOCATED(wt_tab_ct) ) DEALLOCATE(wt_tab_ct) ! tabulated thrust coefficient [-]
    IF( ALLOCATED(wt_tab_wsph_1) ) DEALLOCATE(wt_tab_wsph_1) ! horizontal wind speed at the rotor disk calculated from tabulated data (for opt=1) [-]
    IF( ALLOCATED(wt_tab_aif) ) DEALLOCATE(wt_tab_aif) ! axial induction factor calculated from tabulated data [-]
    IF( ALLOCATED(wt_mod_hub_height) ) DEALLOCATE(wt_mod_hub_height) ! modified wind turbine hub height to account for wind turbine base height, wt_z(n) [m]
    IF( ALLOCATED(wt_min_x) ) DEALLOCATE(wt_min_x) ! minimum distance of the wind turbine circumference in x-direction [m]
    IF( ALLOCATED(wt_max_x) ) DEALLOCATE(wt_max_x) ! maximum distance of the wind turbine circumference in x-direction [m]
    IF( ALLOCATED(wt_min_y) ) DEALLOCATE(wt_min_y) ! minimum distance of the wind turbine circumference in y-direction [m]
    IF( ALLOCATED(wt_max_y) ) DEALLOCATE(wt_max_y) ! maximum distance of the wind turbine circumference in y-direction [m]
    IF( ALLOCATED(wt_yaw_angle) ) DEALLOCATE(wt_yaw_angle) ! wind turbine yaw angle [deg]
    IF( ALLOCATED(wt_delta_yaw_angle_l) ) DEALLOCATE(wt_delta_yaw_angle_l) ! [locally computed] wind turbine yaw angle increment at each timestep [deg]
    IF( ALLOCATED(wt_yaw_rate) ) DEALLOCATE(wt_yaw_rate) ! wind turbine yaw rate [deg/s]
    IF( ALLOCATED(wt_yaw_error_l) ) DEALLOCATE(wt_yaw_error_l) ! [locally computed] wind turbine accumulated yaw error [deg]
    IF( ALLOCATED(wt_yaw_error_threshold) ) DEALLOCATE(wt_yaw_error_threshold) ! yaw error threshold beyond which wind turbine starts to yaw [deg^2/s]
    IF( ALLOCATED(wt_ck) ) DEALLOCATE(wt_ck) ! turbulence kinetic energy coefficient, ck=ctke [-]
    IF( ALLOCATED(wt_do_yaw) ) DEALLOCATE(wt_do_yaw) ! switch to check whether wind turbine is yawing or not [-]
    IF( ALLOCATED(wt_do_yaw_l) ) DEALLOCATE(wt_do_yaw_l) ! [locally computed] switch to check whether wind turbine is yawing or not [-]
    IF( ALLOCATED(wt_dir_l) ) DEALLOCATE(wt_dir_l) ! [locally computed] wind turbine shaft direction at hub height [deg]
    IF( ALLOCATED(wt_unit_vector) ) DEALLOCATE(wt_unit_vector) ! wind turbine unit vector pointing from wind turbine base to the hub [-]
    IF( ALLOCATED(wt_unit_vector_shaft_dir) ) DEALLOCATE(wt_unit_vector_shaft_dir) ! wind turbine unit vector shaft direction at hub height [-]
    IF( ALLOCATED(wt_dir_init) ) DEALLOCATE(wt_dir_init) ! wind direction at hub height at wind turbine location [deg]
    IF( ALLOCATED(wt_k_bot) ) DEALLOCATE(wt_k_bot) ! blade tip k-index below hub height [-]
    IF( ALLOCATED(wt_k_hub) ) DEALLOCATE(wt_k_hub) ! blade tip k-index above hub height [-]
    IF( ALLOCATED(wt_k_top) ) DEALLOCATE(wt_k_top) ! k-index at hub height [-]
    IF( ALLOCATED(wt_vel_hub) ) DEALLOCATE(wt_vel_hub) ! hub height wind velocity at rotor disk [m/s]
    IF( ALLOCATED(wt_aif_l) ) DEALLOCATE(wt_aif_l) ! [locally computed] axial induction factor calculated from tabulated data [-]
    IF( ALLOCATED(wt_ct_l) ) DEALLOCATE(wt_ct_l) ! [locally computed] wind turbine thrust coefficient [-]
    IF( ALLOCATED(wt_cp_l) ) DEALLOCATE(wt_cp_l) ! [locally computed] wind turbine power coefficient [-]
    IF( ALLOCATED(wt_ck_l) ) DEALLOCATE(wt_ck_l) ! [locally computed] turbulence kinetic energy coefficient [-]
    IF( ALLOCATED(wt_power_l) ) DEALLOCATE(wt_power_l) ! [locally computed] wind turbine power [W]
    IF( ALLOCATED(wt_yaw_angle_l) ) DEALLOCATE(wt_yaw_angle_l) ! [locally computed] wind turbine yaw angle [deg]
!
    IF( ALLOCATED(dist_wt_base_to_hub) ) DEALLOCATE(dist_wt_base_to_hub) ! distance between wind turbine base and hub
!
    IF( ALLOCATED(dz_cst) ) DEALLOCATE(dz_cst) ! vertical grid distance at hub height [m]
!
    IF( ALLOCATED(sad_n_points) ) DEALLOCATE(sad_n_points) ! number of SAD points [-]
    IF( ALLOCATED(sad_n_points_l) ) DEALLOCATE(sad_n_points_l) ! [locally computed] number of SAD points [-]
    IF( ALLOCATED(sad_area_fact_l) ) DEALLOCATE(sad_area_fact_l) ! [locally computed] SAD grid area factor [-]
    IF( ALLOCATED(sad_i) ) DEALLOCATE(sad_i) ! SAD point i-index [-]
    IF( ALLOCATED(sad_j) ) DEALLOCATE(sad_j) ! SAD point j-index [-]
    IF( ALLOCATED(sad_k) ) DEALLOCATE(sad_k) ! SAD point k-index [-]
    IF( ALLOCATED(sad_i_l) ) DEALLOCATE(sad_i_l) ! [locally computed] SAD point i-index [-]
    IF( ALLOCATED(sad_j_l) ) DEALLOCATE(sad_j_l) ! [locally computed] SAD point j-index [-]
    IF( ALLOCATED(sad_k_l) ) DEALLOCATE(sad_k_l) ! [locally computed] SAD point k-index [-]
!
    IF( ALLOCATED(fst_i_l) ) DEALLOCATE(fst_i_l) ! [locally computed] i-index of freestream point [-]
    IF( ALLOCATED(fst_j_l) ) DEALLOCATE(fst_j_l) ! [locally computed] j-index of freestream point [-]
    IF( ALLOCATED(fst_i_init) ) DEALLOCATE(fst_i_init) ! i-index of freestream point at initialization step [-]
    IF( ALLOCATED(fst_j_init) ) DEALLOCATE(fst_j_init) ! j-index of freestream point at initialization step [-]
    IF( ALLOCATED(fst_dist_cst) ) DEALLOCATE(fst_dist_cst) ! a constant multiplier of rotor diameter to calculate distance between freestream (inflow) location
                                                           ! and wind turbine, e.g., fst_dist_cst*wt_rotor_diameter ==> fst_dist_cst = 1.0 [-]
    IF( ALLOCATED(fst_wdir_ave_l) ) DEALLOCATE(fst_wdir_ave_l) ! [locally computed] mean freestream wind direction at hub height [deg]
    IF( ALLOCATED(fst_rho_ave_l) ) DEALLOCATE(fst_rho_ave_l) ! [locally computed] mean freestream air density at hub height [kg/m3]
    IF( ALLOCATED(fst_uinf_l) ) DEALLOCATE(fst_uinf_l) ! [locally computed] mean freestream wind velocity [m/s]
    IF( ALLOCATED(fst_wdir_ts_l) ) DEALLOCATE(fst_wdir_ts_l) ! [locally computed] time series of mean freestream wind direction at hub height [deg]
    IF( ALLOCATED(fst_u_ts_l) ) DEALLOCATE(fst_u_ts_l) ! [locally computed] time series of mean freestream u-velocity component at hub height [m/s]
    IF( ALLOCATED(fst_v_ts_l) ) DEALLOCATE(fst_v_ts_l) ! [locally computed] time series of mean freestream v-velocity component at hub height [m/s]
    IF( ALLOCATED(fst_rho_ts_l) ) DEALLOCATE(fst_rho_ts_l) ! [locally computed] time series of mean freestream air density at hub height [m/s]
#ifdef DM_PARALLEL
    IF( ALLOCATED(mpi_sender) ) DEALLOCATE(mpi_sender) ! processor id for mpi_sender [-]
    IF( ALLOCATED(mpi_proc) ) DEALLOCATE(mpi_proc) ! processor for senders [-]
    IF( ALLOCATED(mpi_proc_id) ) DEALLOCATE(mpi_proc_id) ! processor ids for senders [-]
    IF( ALLOCATED(mpi_proc_fst_init) ) DEALLOCATE(mpi_proc_fst_init) ! initialized processor for senders at freestream location [-]
    IF( ALLOCATED(mpi_proc_id_fst_init) ) DEALLOCATE(mpi_proc_id_fst_init) ! initialized processor ids for senders at freestream location [-]
    IF( ALLOCATED(mpi_new_group) ) DEALLOCATE(mpi_new_group) ! group of processors [-]
    IF( ALLOCATED(mpi_new_comm) ) DEALLOCATE(mpi_new_comm) ! communicator for processors [-]
#endif
    IF( ALLOCATED(u_plane_l) ) DEALLOCATE(u_plane_l) ! [locally-computed] plane-element u-velocity component [m/s]
    IF( ALLOCATED(v_plane_l) ) DEALLOCATE(v_plane_l) ! [locally-computed] plane-element v-velocity component [m/s]
    IF( ALLOCATED(w_plane_l) ) DEALLOCATE(w_plane_l) ! [locally-computed] plane-element w-velocity component [m/s]
    IF( ALLOCATED(plane_points) ) DEALLOCATE(plane_points) ! reference plane-element locations [m]
! ... Allocate and initialize the configuration of the wind turbine(s) ...
    ALLOCATE( spec_wt(1:n_pln_wt) )
    ALLOCATE( wt_x(1:n_wt) )
    ALLOCATE( wt_y(1:n_wt) )
    ALLOCATE( wt_z(1:n_wt) )
    ALLOCATE( wt_type_id(1:n_wt) )
    ALLOCATE( wt_i(1:n_wt) )
    ALLOCATE( wt_j(1:n_wt) )
    ALLOCATE( wt_lat(1:n_wt) )
    ALLOCATE( wt_lon(1:n_wt) )
    spec_wt(:) = -999999
    wt_x(:) = -999999.999
    wt_y(:) = -999999.999
    wt_z(:) = -999999.999
    wt_type_id(:) = -999999
    wt_i(:) = -999999
    wt_j(:) = -999999
    wt_lat(:) = -999999.999
    wt_lon(:) = -999999.999
    IF( wrf_dm_on_monitor() ) THEN
      DO n = 1, n_wt
        IF( config_flags%windfarm_ij .EQ. 1 ) THEN
          READ(70, *) wt_x(n), wt_y(n), wt_z(n), wt_type_id(n)
          wt_i(n) = INT( (wt_x(n) + 0.5*dx)/dx )
          wt_j(n) = INT( (wt_y(n) + 0.5*dy)/dy )
        ELSE
          READ(70, *) wt_lat(n), wt_lon(n), wt_z(n), wt_type_id(n)
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
          CALL latlon_to_ij(ts_proj, wt_lat(n), wt_lon(n), ts_rx, ts_ry)
          wt_x(n) = ts_rx*dx
          wt_y(n) = ts_ry*dy
          wt_i(n) = INT( (wt_x(n) + 0.5*dx)/dx )
          wt_j(n) = INT( (wt_y(n) + 0.5*dy)/dy )
        ENDIF
      ENDDO
      CLOSE(70)
    ENDIF
    CALL wrf_dm_bcast_integer(wt_type_id, n_wt)
    n_wt_types = MAXVAL(wt_type_id)
    ALLOCATE( wt_tab_n_lines(1:n_wt_types) )
    ALLOCATE( wt_cutin(1:n_wt_types) )
    ALLOCATE( wt_cutout(1:n_wt_types) )
    ALLOCATE( wt_hub_height(1:n_wt_types) )
    ALLOCATE( wt_rotor_diameter(1:n_wt_types) )
    ALLOCATE( wt_std_ct(1:n_wt_types) )
    ALLOCATE( wt_rated_power(1:n_wt_types) )
    ALLOCATE( wt_norm_radius(1:n_wt_types) )
    ALLOCATE( wt_area(1:n_wt_types) )
    ALLOCATE( wt_tab_wsph(1:n_wt_types, 1:n_max_val) )
    ALLOCATE( wt_tab_power(1:n_wt_types, 1:n_max_val) )
    ALLOCATE( wt_tab_ct(1:n_wt_types, 1:n_max_val) )
    ALLOCATE( wt_tab_wsph_1(1:n_wt_types, 1:n_max_val) )
    ALLOCATE( wt_tab_aif(1:n_wt_types, 1:n_max_val) )
!=======================================================================
! Initialize allocated variables:
!=======================================================================
    wt_tab_n_lines(:) = -999999
    wt_cutin(:) = -999999.999
    wt_cutout(:) = -999999.999
    wt_hub_height(:) = -999999.999
    wt_rotor_diameter(:) = -999999.999
    wt_std_ct(:) = -999999.999
    wt_rated_power(:) = -999999.999
    wt_norm_radius(:) = -999999.999
    wt_area(:) = -999999.999
    wt_tab_wsph(:, :) = -999999.999
    wt_tab_power(:, :) = -999999.999
    wt_tab_ct(:, :) = -999999.999
    wt_tab_wsph_1(:, :) = -999999.999
    wt_tab_aif(:, :) = -999999.999
    IF( wrf_dm_on_monitor() ) THEN
      DO n = 1, n_wt_types
        WRITE(wt_num, *) n
        wt_num=TRIM( ADJUSTL(wt_num) )
        input = "wind-turbine-"//TRIM(wt_num)//".tbl"
        OPEN(FILE=TRIM(input), UNIT=15, FORM='FORMATTED', STATUS='OLD')
        READ(15, *, ERR=110) n_val
        READ(15, *, ERR=110) wt_hub_height(n), wt_rotor_diameter(n), &
                             wt_std_ct(n), wt_rated_power(n)
        wt_norm_radius(n) = wt_rotor_diameter(n)/dx
        wt_area(n) = 0.25*piconst*wt_rotor_diameter(n)**2
        ALLOCATE( wt_tab_wsph_tmp(1:n_val) )
        ALLOCATE( wt_tab_ct_tmp(1:n_val) )
        ALLOCATE( wt_tab_power_tmp(1:n_val) )
        wt_tab_wsph_tmp(:) = 0.0
        wt_tab_ct_tmp(:) = 0.0
        wt_tab_power_tmp(:) = 0.0
        DO m = 1, n_val
          READ(15, *, ERR=110) wt_tab_wsph_tmp(m), wt_tab_ct_tmp(m), wt_tab_power_tmp(m)
        ENDDO
        CLOSE(15)
        wt_tab_n_lines(n) = n_val
        wt_cutin(n) = wt_tab_wsph_tmp(1)
        wt_cutout(n) = wt_tab_wsph_tmp(n_val)
        wt_tab_wsph(n, 1:n_val) = wt_tab_wsph_tmp(1:n_val)
        wt_tab_ct(n, 1:n_val) = wt_tab_ct_tmp(1:n_val)
        wt_tab_power(n, 1:n_val) = wt_tab_power_tmp(1:n_val)
        ! Calculate wind speed at rotor disk from tabulated data:
        DO m = 1, n_val
          wt_tab_wsph_1(n, m) = 0.5*wt_tab_wsph(n, m)*(1 + SQRT( 1-wt_tab_ct(n, m) ) )
          wt_tab_aif(n, m) = 0.5*(1-SQRT( 1-wt_tab_ct(n, m) ) )
        ENDDO
        DEALLOCATE( wt_tab_wsph_tmp, wt_tab_ct_tmp, wt_tab_power_tmp)
 110    CONTINUE
      ENDDO
    ENDIF
!=======================================================================
! Broadcast data to all processors:  
#ifdef DM_PARALLEL
    CALL wrf_dm_bcast_integer(wt_i, n_wt)
    CALL wrf_dm_bcast_integer(wt_j, n_wt)
    CALL wrf_dm_bcast_integer(wt_tab_n_lines, n_wt_types)
!
    CALL wrf_dm_bcast_real(wt_x, n_wt)
    CALL wrf_dm_bcast_real(wt_y, n_wt)
    CALL wrf_dm_bcast_real(wt_z, n_wt)
    CALL wrf_dm_bcast_real(wt_lat, n_wt)
    CALL wrf_dm_bcast_real(wt_lon, n_wt)
    CALL wrf_dm_bcast_real(wt_hub_height, n_wt_types)
    CALL wrf_dm_bcast_real(wt_rotor_diameter, n_wt_types)
    CALL wrf_dm_bcast_real(wt_cutin, n_wt_types)
    CALL wrf_dm_bcast_real(wt_cutout, n_wt_types)
    CALL wrf_dm_bcast_real(wt_area, n_wt_types)
    CALL wrf_dm_bcast_real(wt_norm_radius, n_wt_types)
    CALL wrf_dm_bcast_real(wt_rated_power, n_wt_types)
    CALL wrf_dm_bcast_real(wt_std_ct, n_wt_types)
    CALL wrf_dm_bcast_real(wt_tab_wsph, n_wt_types*n_max_val)
    CALL wrf_dm_bcast_real(wt_tab_wsph_1, n_wt_types*n_max_val)
    CALL wrf_dm_bcast_real(wt_tab_aif, n_wt_types*n_max_val)
    CALL wrf_dm_bcast_real(wt_tab_ct, n_wt_types*n_max_val)
    CALL wrf_dm_bcast_real(wt_tab_power, n_wt_types*n_max_val)
#endif
! orig
!    sad_max_points = CEILING( (piconst*max_radius**2)/(min_dx*min_dz) )
!    row_max_points = CEILING(max_radius/min_dx)*2
! approximated
!    sad_max_points = CEILING( (piconst*(MAXVAL(wt_rotor_diameter)/2)**2)/(dx**2) )
!    row_max_points = CEILING((MAXVAL(wt_rotor_diameter)/2)/dx)*2
!=======================================================================
! Allocate arrays:
!=======================================================================
#ifdef DM_PARALLEL
    ALLOCATE( mpi_sender(1:n_wt) )
    ALLOCATE( mpi_proc(1:n_wt) )
    ALLOCATE( mpi_proc_id(1:n_wt) )
    ALLOCATE( mpi_proc_fst_init(1:n_wt) )
    ALLOCATE( mpi_proc_id_fst_init(1:n_wt) )
    ALLOCATE( mpi_new_group(1:n_wt) )
    ALLOCATE( mpi_new_comm(1:n_wt) )
#endif
!
!    ALLOCATE( sad_i(1:n_wt, 1:sad_max_points) )
!    ALLOCATE( sad_j(1:n_wt, 1:sad_max_points) )
!    ALLOCATE( sad_k(1:n_wt, 1:sad_max_points) )
    ALLOCATE( sad_n_points(1:n_wt) )
!
    ALLOCATE( wt_mod_hub_height(1:n_wt) )
    ALLOCATE( wt_min_x(1:n_wt) )
    ALLOCATE( wt_max_x(1:n_wt) )
    ALLOCATE( wt_min_y(1:n_wt) )
    ALLOCATE( wt_max_y(1:n_wt) )
    ALLOCATE( wt_yaw_rate(1:n_wt_types) )
    ALLOCATE( wt_yaw_error_threshold(1:n_wt_types) )
!
    ALLOCATE( wt_k_bot(1:n_wt) )
    ALLOCATE( wt_k_hub(1:n_wt) )
    ALLOCATE( wt_k_top(1:n_wt) )
    ALLOCATE( wt_unit_vector_shaft_dir(1:n_wt) )
    ALLOCATE( wt_dir_init(1:n_wt) )
    ALLOCATE( wt_do_yaw(1:n_wt) )
    ALLOCATE( wt_yaw_angle(1:n_wt) )
    ALLOCATE( wt_vel_hub(1:n_wt) )
    ALLOCATE( wt_ck(1:n_wt) )
!
    ALLOCATE( wt_unit_vector(1:n_wt, 1:3) )
!
    ALLOCATE( dist_wt_base_to_hub(1:n_wt, 1:3) )
!
    ALLOCATE( dz_cst(1:n_wt) )
!
    ALLOCATE( fst_i_init(1:n_wt) )
    ALLOCATE( fst_j_init(1:n_wt) )
    ALLOCATE( fst_dist_cst(1:n_wt) )   
!
    ALLOCATE( u_plane_l(1:n_pln_wt, 1:n_vert_points, 1:n_horz_points) )
    ALLOCATE( v_plane_l(1:n_pln_wt, 1:n_vert_points, 1:n_horz_points) )
    ALLOCATE( w_plane_l(1:n_pln_wt, 1:n_vert_points, 1:n_horz_points) )
!
    ALLOCATE( plane_points(1:n_pln_wt, 1:n_vert_points, 1:n_horz_points, 1:3) )
!=======================================================================
! Allocate [locally computed] arrays (reduced to all processors):
!=======================================================================
    ALLOCATE( wt_aif_l(1:n_wt) )
    ALLOCATE( wt_ct_l(1:n_wt) )
    ALLOCATE( wt_cp_l(1:n_wt) )
    ALLOCATE( wt_ck_l(1:n_wt) )
    ALLOCATE( wt_power_l(1:n_wt) )
    ALLOCATE( wt_yaw_angle_l(1:n_wt) )
    ALLOCATE( wt_yaw_error_l(1:n_wt) )
    ALLOCATE( wt_do_yaw_l(1:n_wt) )
    ALLOCATE( wt_delta_yaw_angle_l(1:n_wt) )
    ALLOCATE( wt_dir_l(1:n_wt) )
!
    ALLOCATE( sad_n_points_l(1:n_wt) )
!    ALLOCATE( sad_area_fact_l(1:n_wt, 1:sad_max_points) )
!    ALLOCATE( sad_i_l(1:n_wt, 1:sad_max_points) )
!    ALLOCATE( sad_j_l(1:n_wt, 1:sad_max_points) )
!    ALLOCATE( sad_k_l(1:n_wt, 1:sad_max_points) )
!
    ALLOCATE( fst_i_l(1:n_wt) )
    ALLOCATE( fst_j_l(1:n_wt) )
    ALLOCATE( fst_wdir_ave_l(1:n_wt) )
    ALLOCATE( fst_rho_ave_l(1:n_wt) )
    ALLOCATE( fst_uinf_l(1:n_wt) )
!
    ALLOCATE( fst_wdir_ts_l(1:n_wt, 1:n_ave_samples) )
    ALLOCATE( fst_u_ts_l(1:n_wt, 1:n_ave_samples) )
    ALLOCATE( fst_v_ts_l(1:n_wt, 1:n_ave_samples) )
    ALLOCATE( fst_rho_ts_l(1:n_wt, 1:n_ave_samples) )
!=======================================================================
! Initialize allocated arrays:
!=======================================================================
#ifdef DM_PARALLEL
    mpi_sender(:) = 0
    mpi_proc_id(:) = 0
    mpi_proc_id_fst_init(:) = 0
    mpi_proc(:) = 0
    mpi_proc_fst_init(:) = 0
    mpi_new_group(:) = 0
    mpi_new_comm(:) = 0
#endif
    sad_n_points(:) = -999999
!    sad_i(:, :) = -999999
!    sad_j(:, :) = -999999
!    sad_k(:, :) = -999999
!
    fst_i_init(:) = -999999
    fst_j_init(:) = -999999
    fst_dist_cst(:) = config_flags%wind_wtp_dist_cst
!
    wt_do_yaw(:) = .FALSE.
    wt_do_yaw_l(:) = .FALSE.
    wt_yaw_rate(:) = config_flags%wind_wtp_yaw_rate
    wt_yaw_error_threshold(:) = config_flags%wind_wtp_yaw_threshold
    wt_min_x(:) = -999999.999
    wt_max_x(:) = -999999.999
    wt_min_y(:) = -999999.999
    wt_max_y(:) = -999999.999
    wt_k_bot(:) = -1
    wt_k_hub(:) = -1
    wt_k_top(:) = -1
!
    dz_cst(:) = -999999.999
!
    wt_ck(:) = -999999.999
    wt_dir_init(:) = -999999.999
    wt_yaw_angle(:) = -999999.999
    wt_delta_yaw_angle_l(:) = -999999.999
    wt_vel_hub(:) = -999999.999
    wt_mod_hub_height(:) = -999999.999
    wt_unit_vector(:, :) = -999999.999
    wt_unit_vector_shaft_dir(:) = -999999.999
!
    dist_wt_base_to_hub(:, :) = -999999.999
!
    plane_points(:, :, :, :) = -999999.999
!=======================================================================
! Initialize [locally computed] arrays (reduced to all processors):
!=======================================================================
    wt_yaw_angle_l(:) = -999999.999
    wt_delta_yaw_angle_l(:) = -999999.999
    wt_yaw_error_l(:) = -999999.999
    wt_dir_l(:) = -999999.999
    wt_ct_l(:) = -999999.999
    wt_cp_l(:) = -999999.999
    wt_ck_l(:) = -999999.999
    wt_aif_l(:) = -999999.999
    wt_power_l(:) = -999999.999
!
    fst_i_l(:) = -999999
    fst_j_l(:) = -999999
    fst_rho_ave_l(:) = -999999.999
    fst_wdir_ave_l(:) = -999999.999
    fst_uinf_l(:) = -999999.999
    fst_wdir_ts_l(:, :) = 0.0
    fst_u_ts_l(:, :) = 0.0
    fst_v_ts_l(:, :) = 0.0
    fst_rho_ts_l(:, :) = 0.0
!
    sad_n_points_l(:) = -999999.999
!    sad_area_fact_l(:, :) = -999999.999
!    sad_i_l(:, :) = -999999
!    sad_j_l(:, :) = -999999
!    sad_k_l(:, :) = -999999
!
    u_plane_l(:, :, :) = -999999.999
    v_plane_l(:, :, :) = -999999.999
    w_plane_l(:, :, :) = -999999.999
!=======================================================================
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
        ENDDO
      ENDDO
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
#ifdef DM_PARALLEL
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size_comm, ierr)
#endif
    spec_wt = model_config_rec%wind_wtp_pln_turbine_ids(1:n_pln_wt)
    n_pln_wt_counter = 1
    DO n = 1, n_wt
      wt_min_x(n) = wt_x(n) - wt_rotor_diameter( wt_type_id(n) )/2
      wt_min_y(n) = wt_y(n) - wt_rotor_diameter( wt_type_id(n) )/2
      wt_max_x(n) = wt_x(n) + wt_rotor_diameter( wt_type_id(n) )/2
      wt_max_y(n) = wt_y(n) + wt_rotor_diameter( wt_type_id(n) )/2
      ! modified wind turbine hub height to account for wind turbine base height, wt_z(n)
      wt_mod_hub_height( wt_type_id(n) ) = wt_z(n) + wt_hub_height( wt_type_id(n) )
      IF(its .LE. wt_i(n) .AND. wt_i(n) .LE. itf .AND. jts .LE. wt_j(n) .AND. wt_j(n) .LE. jtf) THEN
        DO k = kts, ktf
          IF( ( w_z(wt_i(n), k, wt_j(n) ) .GE. &
              (wt_mod_hub_height( wt_type_id(n) )-wt_rotor_diameter( wt_type_id(n) )/2) ) .AND. &
              (wt_k_bot(n) .EQ. -1) ) wt_k_bot(n) = k
          IF( ( w_z(wt_i(n), k, wt_j(n) ) .GE. &
              (wt_mod_hub_height( wt_type_id(n) ) ) .AND. (wt_k_hub(n) .EQ. -1) ) ) wt_k_hub(n) = k
          IF( ( w_z(wt_i(n), k, wt_j(n) ) .GE. &
              (wt_mod_hub_height( wt_type_id(n) )+wt_rotor_diameter( wt_type_id(n) )/2) ) .AND. &
              (wt_k_top(n) .EQ. -1) ) THEN
              wt_k_top(n) = k-1
              EXIT
          ENDIF
        ENDDO
        dz_cst(n) = dz( wt_i(n), wt_k_hub(n), wt_j(n) )
        wt_dir_init(n) = ATAN2( u_2( wt_i(n), wt_k_hub(n), wt_j(n) ), v_2( wt_i(n), wt_k_hub(n), wt_j(n) ) )/degrad + 180.0
        IF( opt_error_check ) THEN
          ! Check division by zero for running averages:
          IF( n_ave_samples .LE. 0 ) THEN
            CALL wrf_error_fatal('WTP-ERROR ***** wind_wtp_sad_n_aave_samples &
                                  and wind_wtp_sad_n_series_tme &
                                  are not set in namelist.input properly. *****')
          ENDIF
          ! Check indices of the turbine dimensions:
          IF( ( wt_min_x(n) .LT. (ids-1)*dx .OR. wt_max_x(n) .GT. ide*dx ) .OR. &
              ( wt_min_y(n) .LT. (jds-1)*dy .OR. wt_max_x(N) .GT. jde*dy ) ) THEN
            CALL wrf_error_fatal('WTP-ERROR ***** Turbine dimensions are outside of the domain. &
                                  Relocate wind turbine(s) properly. *****')
          ENDIF
          ! Check non-positive indices in vertical direction:
          IF( wt_k_bot(n) .LE. 0 .OR. wt_k_hub(n) .LE. 0 .OR. wt_k_top(n) .LE. 0 ) THEN
            CALL wrf_error_fatal('WTP-ERROR ***** Non-positive "k" in vertical direction is found. &
                                  Check WRF-LES setup prior to the WTP-SAD simulation. *****')
          ENDIF
          DO k = wt_k_bot(n), wt_k_top(n)
            ! Check model grid spacing:
            IF( NINT(wt_rotor_diameter( wt_type_id(n) )/dx) .LT. 2 .OR. &
                NINT(wt_rotor_diameter( wt_type_id(n) )/dy) .LT. 2 .OR. &
                NINT(wt_rotor_diameter( wt_type_id(n) )/dz_cst( wt_type_id(n) ) ) .LT. 2 ) THEN
              WRITE(outstring, *)''
              CALL wrf_message ( TRIM(outstring) )
              WRITE(outstring, "(A, I3, A, I3, A, I3, A)") &
                               'WTP-WARNING ***** wt_k_bot = ', wt_k_bot( wt_type_id(n) ), &
                               ', wt_k_hub = ', wt_k_hub( wt_type_id(n) ), ', wt_k_top = ', &
                               wt_k_top( wt_type_id(n) ), ' *****'
              CALL wrf_message ( TRIM(outstring) )
              WRITE(outstring, "(A, F7.3, A, A, F7.3, A, A, F7.3, A, A, F7.3, A, I3, A)")                        &
                               'WTP-WARNING ***** rotor diameter = ', wt_rotor_diameter( wt_type_id(n) ), ' m,', &
                               ' dx = ', dx, ' m,', ' dy = ', dy, ' m,', ' dz = ', dz_cst( wt_type_id(n) ), &
                               ' m at k = ', k, ' *****'
              CALL wrf_message ( TRIM(outstring) )
              WRITE(outstring, "(A, I3, A, I3, A, I3, A)") &
                               'WTP-WARNING ***** (rotor diamater/dx) = ', NINT(wt_rotor_diameter( wt_type_id(n) )/dx), &
                               ' (rotor diameter/dy) = ', NINT(wt_rotor_diameter( wt_type_id(n) )/dy),                  &
                               ' and (rotor diameter/dz) = ', NINT(wt_rotor_diameter( wt_type_id(n) )/dz_cst( wt_type_id(n) ) ), &
                               ' must be equal to or greater than 2. *****'
              CALL wrf_message ( TRIM(outstring) )
              WRITE(outstring, *)''
              CALL wrf_message ( TRIM(outstring) )
              CALL wrf_error_fatal('WTP-ERROR ***** Increase spatial resolution for correct &
                                    modeling. dx, dy, dz values should be as close to each &
                                    other as possible. *****')
            ENDIF
          ENDDO
        ENDIF
#ifdef DM_PARALLEL
        CALL wrf_get_myproc( mpi_proc(n) )
#endif
      ENDIF
#ifdef DM_PARALLEL
      mpi_proc_id(n) = wrf_dm_max_int( mpi_proc(n) )
      CALL MPI_BCAST(dz_cst(n), 1, MPI_REAL, mpi_proc_id(n), MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(wt_dir_init(n), 1, MPI_REAL, mpi_proc_id(n), MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(wt_k_bot(n), 1, MPI_INT, mpi_proc_id(n), MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(wt_k_hub(n), 1, MPI_INT, mpi_proc_id(n), MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(wt_k_top(n), 1, MPI_INT, mpi_proc_id(n), MPI_COMM_WORLD, ierr)
#endif
      CALL compass_to_standard( wt_dir_init(n) )
      wt_yaw_angle(n) = wt_dir_init(n)
      ! Define vector aligned with the wind turbine base pointing from the ground to the hub
      dist_wt_base_to_hub(n, :) = (/wt_x(n), wt_y(n), wt_mod_hub_height( wt_type_id(n) )/)
      wt_unit_vector(n, :) = vector_add( dist_wt_base_to_hub(n, :), -(/wt_x(n), wt_y(n), wt_z(n)/) )
      wt_unit_vector(n, :) = vector_divide(wt_unit_vector(n, :), vector_mag( wt_unit_vector(n, :) ) )
      IF( opt_tsout_plane .AND. ANY(spec_wt .EQ. n) ) THEN
        IF( .NOT. restart ) THEN
          wt_unit_vector_shaft_dir(n) = -1.0
          dist_wt_to_plane = wt_unit_vector_shaft_dir(n)*plane_distance
          CALL init_wind_plane_tsout_locs(                        &
              &  n_pln_wt, n_pln_wt_counter, n_vert_points        &
              &, n_horz_points, dx, dy, plane_width, plane_height &
              &, dist_wt_to_plane, dist_wt_base_to_hub(n, :)      &
              &, wt_unit_vector(n, :), wt_yaw_angle(n), wt_x(n)   &
              &, wt_y(n), wt_z(n), plane_points                   &
              &)
          n_pln_wt_counter = n_pln_wt_counter + 1
        ELSE
          IF( (plane_points_x(n_pln_wt_counter, n_vert_points, n_horz_points) .EQ. 0.0) .AND. &
              (plane_points_y(n_pln_wt_counter, n_vert_points, n_horz_points) .EQ. 0.0) .AND. &
              (plane_points_z(n_pln_wt_counter, n_vert_points, n_horz_points) .EQ. 0.0) ) THEN
            wt_yaw_angle(n) = wt_dir(n)
            wt_unit_vector_shaft_dir(n) = -1.0
            dist_wt_to_plane = wt_unit_vector_shaft_dir(n)*plane_distance
            CALL init_wind_plane_tsout_locs(                        &
                &  n_pln_wt, n_pln_wt_counter, n_vert_points        &
                &, n_horz_points, dx, dy, plane_width, plane_height &
                &, dist_wt_to_plane, dist_wt_base_to_hub(n, :)      &
                &, wt_unit_vector(n, :), wt_yaw_angle(n), wt_x(n)   &
                &, wt_y(n), wt_z(n), plane_points                   &
                &)
            n_pln_wt_counter = n_pln_wt_counter + 1
          ELSE
            k = n_pln_wt_counter
            DO m = 1, n_horz_points
              DO p = 1, n_vert_points
                plane_points(k, p, m, 1) = plane_points_x(k, p, m)
                plane_points(k, p, m, 2) = plane_points_y(k, p, m)
                plane_points(k, p, m, 3) = plane_points_z(k, p, m)
              ENDDO
            ENDDO
            n_pln_wt_counter = n_pln_wt_counter + 1
          ENDIF
        ENDIF
      ENDIF
      fst_i_init(n) = INT( ( wt_x(n)-fst_dist_cst( wt_type_id(n) ) &
                    * wt_rotor_diameter( wt_type_id(n) )*COS(wt_yaw_angle(n)*degrad) + 0.5*dx )/dx )
      ! Set horizontal boundary index in x-direction to i_tile_start if fst_i_init(n)
      ! is less than its
      IF( fst_i_init(n) .LE. 0 ) fst_i_init(n) = 1
      ! Set horizontal boundary index in x-direction to i_domain_end-1 if fst_i_init(n)
      ! is greater than ide
      IF( fst_i_init(n) .GE. ide ) fst_i_init(n) = ide-1
      fst_j_init(n) = INT( ( wt_y(n)-fst_dist_cst( wt_type_id(n) ) &
                    * wt_rotor_diameter( wt_type_id(n) )*SIN(wt_yaw_angle(n)*degrad) + 0.5*dy )/dy )
      ! Set horizontal boundary index in y-direction to 1 if fst_j_init(n)
      ! is less than 0
      IF( fst_j_init(n) .LE. 0 ) fst_j_init(n) = 1
      ! Set horizontal boundary index in y-direction to j_domain_end-1 if fst_j_init(n)
      ! is greater than jde
      IF( fst_j_init(n) .GE. jde ) fst_j_init(n) = jde-1
      IF( ( its .LE. fst_i_init(n) .AND. fst_i_init(n) .LE. itf ) .AND. &
          ( jts .LE. fst_j_init(n) .AND. fst_j_init(n) .LE. jtf ) ) THEN
#ifdef DM_PARALLEL
        CALL wrf_get_myproc( mpi_proc_fst_init(n) )
#endif
      ENDIF
#ifdef DM_PARALLEL
      mpi_proc_id_fst_init(n) = wrf_dm_max_int( mpi_proc_fst_init(n) )
      IF( n_proc_x .GE. 3 .AND. mpi_size_comm .GT. 10 ) THEN
        procs_W = mpi_proc_id(n)-1
        procs_E = mpi_proc_id(n)+1
        procs_N = mpi_proc_id(n)+n_proc_x
        procs_S = mpi_proc_id(n)-n_proc_x
        procs_NW = procs_N-1
        procs_NE = procs_N+1
        procs_SW = procs_S-1
        procs_SE = procs_S+1
        IF( MOD(mpi_proc_id(n), n_proc_x ) .EQ. 0 ) THEN
          procs_W = -1
          procs_NW = -1
          procs_SW = -1
        ELSEIF( MOD(mpi_proc_id(n), n_proc_x ) .EQ. n_proc_x-1 ) THEN
          procs_E = -1
          procs_NE = -1
          procs_SE = -1
        ENDIF
        IF( mpi_proc_id(n) .EQ. mpi_proc_root .OR. mpi_proc_id(n) .EQ. mpi_proc_root+1 ) THEN
          ALLOCATE( mpi_procs(1:n_wt, 1:9) )
          mpi_procs(n, :) = (/procs_W, procs_E, procs_N, procs_S, procs_NW, procs_NE, procs_SW, procs_SE, mpi_proc_id(n)/)
        ELSE
          ALLOCATE( mpi_procs(1:n_wt, 1:10) )
          mpi_procs(n, :) = (/procs_W, procs_E, procs_N, procs_S, procs_NW, procs_NE, procs_SW, procs_SE, mpi_proc_root, mpi_proc_id(n)/)
        ENDIF
        IF( ANY( mpi_procs(n, :) .LT. mpi_proc_root ) .OR. ANY( mpi_procs(n, :) .GE. mpi_size_comm ) ) THEN
          mpi_procs_tmp = PACK( mpi_procs(n, :), MASK = ( mpi_procs(n, :) .GE. mpi_proc_root .AND. mpi_procs(n, :) .LT. mpi_size_comm ) )
          ALLOCATE( mpi_procs_comm_tmp( 1:SIZE(mpi_procs_tmp, DIM=1) ) )
          arrayIN = mpi_procs_tmp
          CALL removeDupsINT(arrayIN, arrayOUT, n_size)
          mpi_procs_comm_tmp = arrayOUT(1:n_size)
          ALLOCATE( mpi_procs_comm( 1:n_wt, 1:n_size ) )
          mpi_procs_comm(n, 1:n_size) = mpi_procs_comm_tmp(1:n_size)
        ELSE
          ALLOCATE( mpi_procs_comm_tmp( 1:SIZE(mpi_procs, DIM=2) ) )
          arrayIN = mpi_procs(n, :)
          CALL removeDupsINT(arrayIN, arrayOUT, n_size)
          mpi_procs_comm_tmp = arrayOUT(1:n_size)
          ALLOCATE( mpi_procs_comm( 1:n_wt, 1:n_size ) )
          mpi_procs_comm(n, 1:n_size) = mpi_procs_comm_tmp(1:n_size)
        ENDIF
        IF( opt_error_check ) THEN
          IF( FINDLOC(mpi_procs_comm(n, :), VALUE=mpi_proc_id_fst_init(n), DIM=1) .EQ. 0 ) THEN
            CALL wrf_error_fatal('WTP-SAD-ERROR ***** Inflow point cannot be placed outside &
                                  of processor boundaries within communicators. &
                                  Decrease inflow distance (wind_wtp_dist_cst) &
                                  or increase domain size. *****')
          ENDIF
        ENDIF
        ! Get the group or processes of the default communicator
        CALL MPI_COMM_GROUP(MPI_COMM_WORLD, mpi_world_group, ierr)
        ! Keep only the processes of "mpi_procs_comm" in the new group.
        CALL MPI_GROUP_INCL(mpi_world_group, SIZE( mpi_procs_comm(n, 1:n_size) ), mpi_procs_comm(n, 1:n_size), mpi_new_group(n), ierr)
        ! Create the new communicator from that group of processes.
        CALL MPI_COMM_CREATE(MPI_COMM_WORLD, mpi_new_group(n), mpi_new_comm(n), ierr)
        mpi_sender(n) = FINDLOC(mpi_procs_comm(n, 1:n_size), VALUE=mpi_proc_id(n), DIM=1) - 1
        DEALLOCATE(mpi_procs, mpi_procs_comm_tmp, mpi_procs_comm)
      ELSE
        mpi_sender(n) = mpi_proc_id(n)
        mpi_new_comm(n) = MPI_COMM_WORLD
      ENDIF
      IF ( wrf_dm_on_monitor() ) THEN
        IF ( config_flags%windfarm_ij .EQ. 1 ) THEN
          WRITE(message_windturbine, "(A, I5, A, F9.3, A, A, F9.3, A)")                  &
                                     'Wind Turbine #', n, ': LXrotor = ', wt_x(n), ' m', &
                                     ', LYrotor = ',  wt_y(n), ' m'
          CALL wrf_message(message_windturbine)
          WRITE(message_windturbine, "(A, I9, A, I9, A, I4)")                                              &
                                     '                     Irotor = ', wt_i(n), '  ,  Jrotor = ', wt_j(n), &
                                     '  , mpi_proc_id = ', mpi_proc_id(n)
          CALL wrf_message(message_windturbine)
        ELSE
          WRITE(message_windturbine, "(A, I5, A, F13.8, A, F13.8, A, I4, A, I4, A, I4)")                               &
                                     'Wind Turbine #', n ,': Lat = ', wt_lat(n), ', Lon = ', wt_lon(n), ', Irotor = ', &
                                      wt_i(n), ', Jrotor = ', wt_j(n), ', mpi_proc_id = ', mpi_proc_id(n)
          CALL wrf_message(message_windturbine)
        ENDIF
      ENDIF
#endif
    ENDDO
    sad_max_points = 2*CEILING( (piconst*(MAXVAL(wt_rotor_diameter)/2)**2)/( dx*MINVAL(dz_cst) ) )
    row_max_points = 2*CEILING( MAXVAL(wt_rotor_diameter)/dx )
!=======================================================================
! Allocate arrays:
!=======================================================================
    ALLOCATE( sad_i(1:n_wt, 1:sad_max_points) )
    ALLOCATE( sad_j(1:n_wt, 1:sad_max_points) )
    ALLOCATE( sad_k(1:n_wt, 1:sad_max_points) )
!=======================================================================
! Allocate [locally computed] arrays (reduced to all processors):
!=======================================================================
    ALLOCATE( sad_area_fact_l(1:n_wt, 1:sad_max_points) )
    ALLOCATE( sad_i_l(1:n_wt, 1:sad_max_points) )
    ALLOCATE( sad_j_l(1:n_wt, 1:sad_max_points) )
    ALLOCATE( sad_k_l(1:n_wt, 1:sad_max_points) )
!=======================================================================
! Initialize allocated arrays:
!=======================================================================
    sad_i(:, :) = -999999
    sad_j(:, :) = -999999
    sad_k(:, :) = -999999
!=======================================================================
! Initialize [locally computed] arrays (reduced to all processors):
!=======================================================================
    sad_area_fact_l(:, :) = -999999.999
    sad_i_l(:, :) = -999999
    sad_j_l(:, :) = -999999
    sad_k_l(:, :) = -999999
!=======================================================================
#ifdef DM_PARALLEL
    CALL wrf_dm_bcast_real(plane_points, 3*n_pln_wt*n_horz_points*n_vert_points)
#endif
  END SUBROUTINE init_module_wind_sad
!=======================================================================
  SUBROUTINE init_wind_plane_tsout_locs(               &
        &  n_pln_wt, n_pln_wt_counter, nz, nxy, dx, dy &
        &, plane_width, plane_height, dist_plane       &
        &, wt_base_to_hub, wt_unit_vector, yaw_angle   &
        &, Lx, Ly, Lz, plane_points                    &
        &)
!
  IMPLICIT NONE
!
    INTEGER, INTENT(IN) :: n_pln_wt, n_pln_wt_counter, nz, nxy
    REAL, INTENT(IN) :: plane_width, plane_height, dist_plane
    REAL, INTENT(IN) :: dx, dy, yaw_angle, Lx, Ly, Lz
    REAL, DIMENSION(1:3), INTENT(IN) :: wt_base_to_hub, wt_unit_vector
    REAL, DIMENSION(1:n_pln_wt, 1:nz, 1:nxy, 1:3), INTENT(INOUT) :: plane_points
! local
    INTEGER :: k, n, m
    REAL :: dist_z, dist_xy
    REAL :: dxy, dz
!
    dxy = plane_width/nxy
    dz = plane_height/nz
    dist_xy = 0.0
    k = n_pln_wt_counter
    DO m = 1, nxy
      dist_xy = dist_xy + 0.5*dxy
      dist_z = 0.0
      DO n = 1, nz
        dist_z = dist_z + 0.5*dz
        plane_points(k, n, m, 1) = Lx + dist_plane
        plane_points(k, n, m, 2) = Ly - (plane_width/2) + dist_xy
        plane_points(k, n, m, 3) = Lz + dist_z
        dist_z = dist_z + 0.5*dz
        plane_points(k, n, m, :) = rotate_point( plane_points(k, n, m, :), wt_base_to_hub, wt_unit_vector, yaw_angle )
      ENDDO
      dist_xy = dist_xy + 0.5*dxy
    ENDDO
    IF( opt_error_check ) THEN
    ! spec_bdy_width = 5
      IF( (Lx-dist_plane) .LE. (5*dx) .OR. ( Ly-(plane_width/2) ) .LE. (5*dy) ) THEN
        CALL wrf_error_fatal('WTP-ERROR ***** Plane must not be placed &
                              very close to domain boundaries. *****')
      ENDIF
    ENDIF
  END SUBROUTINE init_wind_plane_tsout_locs
!=======================================================================
  SUBROUTINE calc_wind_plane_tsout_vel(       &
        &  ims, ime, jms, jme, kms, kme       &
        &, its, itf, jts, jtf, kts, ktf       &
        &, n_pln_wt, n_pln_wt_counter, nz     &
        &, nxy, dx, dy, dz_cst, dz, wt_do_yaw &
        &, plane_points, dist_wt_base_to_hub  &
        &, wt_unit_vector, wt_delta_yaw       &
        &, c_z, w_z, u_2, v_2, w_2            &
        &, mpi_new_comm, ierr                 &
        &, u_plane_l, v_plane_l, w_plane_l    &
        &)
!
  IMPLICIT NONE
#if defined(DM_PARALLEL) && !defined(STUBMPI)
  INCLUDE "mpif.h"
#endif
    LOGICAL, INTENT(IN) :: wt_do_yaw
    INTEGER, INTENT(IN) :: ims, ime, jms, jme, kms, kme
    INTEGER, INTENT(IN) :: n_pln_wt, nz, nxy, its, itf, jts, jtf, kts, ktf
    INTEGER, INTENT(IN) :: n_pln_wt_counter, mpi_new_comm, ierr
    REAL, INTENT(IN) :: dx, dy, dz_cst, wt_delta_yaw
    REAL, DIMENSION(1:3), INTENT(IN) :: dist_wt_base_to_hub, wt_unit_vector
    REAL, DIMENSION(1:n_pln_wt, 1:nz, 1:nxy, 1:3), INTENT(INOUT) :: plane_points
    REAL, DIMENSION(its:itf, kts:ktf, jts:jtf), INTENT(IN) :: dz, c_z, w_z
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: u_2, v_2, w_2
    REAL, DIMENSION(1:n_pln_wt, 1:nz, 1:nxy), INTENT(OUT) :: u_plane_l, v_plane_l, w_plane_l
! local
    INTEGER :: k, m, n, ii, jj, kk
    REAL, DIMENSION(1:nz, 1:nxy) :: u_plane_les, v_plane_les, w_plane_les
    REAL, DIMENSION(1:nz, 1:nxy) :: u_plane_les_l, v_plane_les_l, w_plane_les_l
!
    k = n_pln_wt_counter
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
          plane_points(k, n, m, :) = rotate_point(plane_points(k, n, m, :), dist_wt_base_to_hub, wt_unit_vector, wt_delta_yaw)
        ENDIF
        ! Interpolation of the u-component:
        ii = INT( (plane_points(k, n, m, 1)/dx ) )
        jj = INT( (plane_points(k, n, m, 2) + 0.5*dy)/dy )
        IF( ( its .LE. ii .AND. ii .LE. itf ) .AND. &
            ( jts .LE. jj .AND. jj .LE. jtf ) ) THEN
          kk = MINLOC( ABS(c_z(ii, :, jj) + 0.5*dz_cst - plane_points(k, n, m, 3) ), 1)
          IF( ( c_z(ii, kk, jj) + 0.5*dz_cst ) .GT. plane_points(k, n, m, 3) ) kk = kk - 1
          IF( kk .EQ. 0 ) kk = kts
          CALL interp_velocity(                                     &
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
          CALL interp_velocity(                                     &
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
          CALL interp_velocity(                                     &
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
END MODULE module_wind_sad

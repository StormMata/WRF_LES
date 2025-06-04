!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine prints out the current value of variables at all specified
! time series locations that are within the current patch.
!
! Baris Kale, 27/04/2022 (reproduced from wrf_timeseries.F)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calc_ts_wtp_initialize( grid )
  USE module_domain, ONLY : domain
  USE module_configure, ONLY : model_config_rec, grid_config_rec_type, model_to_grid_config_rec
  USE module_state_description
  USE module_model_constants
  IMPLICIT NONE
    ! Arguments
    TYPE (domain), INTENT(INOUT) :: grid
    TYPE(grid_config_rec_type) :: config_flags
    ! Externals
    LOGICAL, EXTERNAL :: wrf_dm_on_monitor
    INTEGER, EXTERNAL :: get_unused_unit
    ! Local variables
    INTEGER :: ntsloc_extra_temp
    INTEGER :: i, j, k, iunit
    LOGICAL :: exists
    CHARACTER (LEN=132) :: message
    CHARACTER (LEN=35) :: ts_filename_wtp
    CHARACTER (LEN=19) simulation_start_date
    INTEGER :: simulation_start_year   , &
               simulation_start_month  , &
               simulation_start_day    , &
               simulation_start_hour   , &
               simulation_start_minute , &
               simulation_start_second
    CALL model_to_grid_config_rec(grid%id, model_config_rec, config_flags)
    IF(config_flags%wind_wtp_opt .GT. 0 .AND. config_flags%wind_wtp_opt .LT. 4) THEN
#if ((EM_CORE == 1) && (DA_CORE != 1))
      IF( grid%dfi_stage == DFI_FST ) THEN
#endif
        grid%next_ts_time_wtp = 1
        ! Determine simulation start time
        CALL nl_get_simulation_start_year   ( 1 , simulation_start_year   )
        CALL nl_get_simulation_start_month  ( 1 , simulation_start_month  )
        CALL nl_get_simulation_start_day    ( 1 , simulation_start_day    )
        CALL nl_get_simulation_start_hour   ( 1 , simulation_start_hour   )
        CALL nl_get_simulation_start_minute ( 1 , simulation_start_minute )
        CALL nl_get_simulation_start_second ( 1 , simulation_start_second )
        WRITE(simulation_start_date , FMT = '(I4.4,"-",I2.2,"-",I2.2,"_",I2.2,":",I2.2,":",I2.2)') &
              simulation_start_year,simulation_start_month,simulation_start_day,simulation_start_hour, &
              simulation_start_minute,simulation_start_second
        WRITE(message, '(A55,I3,A15,A19)') 'Computing time series for the GAD/GAL model for domain ', grid%id, &
        ' starting from ',simulation_start_date
        CALL wrf_message(message)
        CALL wind_wtp_tsout_dir(exists)
        IF(exists) THEN
          IF( wrf_dm_on_monitor() ) THEN
            DO j = 1, config_flags%wind_wtp_nTurbines
            WRITE(ts_filename_wtp,'(A)') 'turbine.0000'
              ! Create the output files for the vertical profiles, one file for each variable
              iunit = get_unused_unit()
              IF( iunit .LE. 0 ) THEN
                CALL wrf_error_fatal('Error in calc_ts_wtp_initialize: could not find a free Fortran unit.')
              ENDIF
              i = LEN_TRIM(ts_filename_wtp)
              IF(j .LT. 10  ) WRITE(ts_filename_wtp(i  :i),'(I1)') j
              IF(j .GE. 10  ) WRITE(ts_filename_wtp(i-1:i),'(I2)') j
              IF(j .GE. 100 ) WRITE(ts_filename_wtp(i-2:i),'(I3)') j
              IF(j .GE. 1000) WRITE(ts_filename_wtp(i-3:i),'(I4)') j
              WRITE(ts_filename_wtp,'(A)') 'tsout/'//TRIM(ts_filename_wtp)
              OPEN(UNIT=iunit, FILE=TRIM(ts_filename_wtp), FORM='FORMATTED', STATUS='REPLACE')
#if (EM_CORE == 1)
                WRITE(UNIT=iunit, FMT='(A11,A16,A12,A11,A9,A12,A11,A16,A17, &
                A13,A12,A14,A15,A15,A15,A15,A15,A15,A12,A12,A12,A16,A16,A16)') &
                'timestep[-]', 'time[s]', 'dt[s]', 'rho[kg/m3]',            &
                'TSR[-]', 'rotSpd[rpm]', 'pitch[deg]', 'fstWindSpd[m/s]',   &
                'fstWindDir[deg]', 'turbDir[deg]', 'yawAng[deg]',           &
                'wsRotor[m/s]', 'rotThr[N]', 'rotTrq[Nm]', 'genTrq[Nm]',    &
                'rotPwr[W]', 'mechPwr[W]', 'genPwr[W]', 'thrCoef[-]',       &
                'trqCoef[-]', 'pwrCoef[-]', 'rotApexPos_x[m]',              &
                'rotApexPos_y[m]', 'rotApexPos_z[m]'
#endif
              CLOSE(UNIT=iunit)
              IF(config_flags%wind_wtp_opt .EQ. 3) THEN
                WRITE(ts_filename_wtp,'(A)') 'turbine.0000.disk'
                ! Create the output files for the vertical profiles, one file for each variable
                iunit = get_unused_unit()
                IF( iunit .LE. 0 ) THEN
                  CALL wrf_error_fatal('Error in calc_ts_wtp_initialize: could not find a free Fortran unit.')
                ENDIF
                i = LEN_TRIM(ts_filename_wtp)
                IF(j .LT. 10  ) WRITE(ts_filename_wtp(i-5:i-5),'(I1)') j
                IF(j .GE. 10  ) WRITE(ts_filename_wtp(i-6:i-5),'(I2)') j
                IF(j .GE. 100 ) WRITE(ts_filename_wtp(i-7:i-6),'(I3)') j
                IF(j .GE. 1000) WRITE(ts_filename_wtp(i-8:i-7),'(I4)') j
                WRITE(ts_filename_wtp,'(A)') 'tsout/'//TRIM(ts_filename_wtp)
                OPEN(UNIT=iunit, FILE=TRIM(ts_filename_wtp), FORM='FORMATTED', STATUS='REPLACE')
#if (EM_CORE == 1)
                WRITE(UNIT=iunit, FMT='(A11,A16,A12,A11,A13,A17,A13,A13, &
                A13,A13,A13,A13,A14,A14,A12,A12,A12,A15,A14,A15,A15,A11  &
                )')                                                      &
                'timestep[-]', 'time[s]', 'dt[s]', 'blEl[-]',            &
                'V1[m/s]', 'Vtheta[m/s]', 'Vt[m/s]', 'Vr[m/s]',          &
                'Vz[m/s]', 'u[m/s]', 'v[m/s]', 'w[m/s]', 'Phi[deg]',     &
                'AoA[deg]', 'Cl[-]', 'Cd[-]', 'Cm[-]', 'L[N/m]',         &
                'D[N/m]', 'Ft[N/m]', 'Fn[N/m]', 'F[-]'
#endif
                CLOSE(UNIT=iunit)
              ELSE
                WRITE(ts_filename_wtp,'(A)') 'turbine.0000.blade.000'
                ! Create the output files for the vertical profiles, one file for each variable
                iunit = get_unused_unit()
                IF( iunit .LE. 0 ) THEN
                  CALL wrf_error_fatal('Error in calc_ts_wtp_initialize: could not find a free Fortran unit.')
                ENDIF
                i = LEN_TRIM(ts_filename_wtp)
                IF(j .LT. 10  ) WRITE(ts_filename_wtp(i-10:i-10),'(I1)') j
                IF(j .GE. 10  ) WRITE(ts_filename_wtp(i-11:i-10),'(I2)') j
                IF(j .GE. 100 ) WRITE(ts_filename_wtp(i-12:i-11),'(I3)') j
                IF(j .GE. 1000) WRITE(ts_filename_wtp(i-13:i-12),'(I4)') j
                DO k = 1, config_flags%wind_wtp_nSections
                  WRITE(ts_filename_wtp,'(A)') 'turbine.0000.blade.000'
                  i = LEN_TRIM(ts_filename_wtp)
                  IF(j .LT. 10  ) WRITE(ts_filename_wtp(i-10:i-10),'(I1)') j
                  IF(j .GE. 10  ) WRITE(ts_filename_wtp(i-11:i-10),'(I2)') j
                  IF(j .GE. 100 ) WRITE(ts_filename_wtp(i-12:i-11),'(I3)') j
                  IF(j .GE. 1000) WRITE(ts_filename_wtp(i-13:i-12),'(I4)') j
                  IF(k .LT. 10  ) WRITE(ts_filename_wtp(i  :i),'(I1)') k
                  IF(k .GE. 10  ) WRITE(ts_filename_wtp(i-1:i),'(I2)') k
                  IF(k .GE. 100 ) WRITE(ts_filename_wtp(i-2:i),'(I3)') k
                  WRITE(ts_filename_wtp,'(A)') 'tsout/'//TRIM(ts_filename_wtp)
                  OPEN(UNIT=iunit, FILE=TRIM(ts_filename_wtp), FORM='FORMATTED', STATUS='REPLACE')
#if (EM_CORE == 1)
                  WRITE(UNIT=iunit, FMT='(A11,A16,A12,A11,A13,A17,A13,A13, &
                  A13,A13,A13,A13,A14,A14,A12,A12,A12,A15,A14,A15,A15,A11, &
                  A15,A15,A15)')                                           &
                  'timestep[-]', 'time[s]', 'dt[s]', 'blEl[-]',            &
                  'V1[m/s]', 'Vtheta[m/s]', 'Vt[m/s]', 'Vr[m/s]',          &
                  'Vz[m/s]', 'u[m/s]', 'v[m/s]', 'w[m/s]', 'Phi[deg]',     &
                  'AoA[deg]', 'Cl[-]', 'Cd[-]', 'Cm[-]', 'L[N/m]',         &
                  'D[N/m]', 'Ft[N/m]', 'Fn[N/m]', 'F[-]',                  &
                  'blPts_x[m]', 'blPts_y[m]', 'blPts_z[m]' 
#endif
                  CLOSE(UNIT=iunit)
                ENDDO
              ENDIF
            ENDDO
          ENDIF
        ELSE
          CALL wrf_error_fatal('Error in calc_ts_wtp_initialize: No tsout directory is found. &
                                Create "./tsout/" directory to save time series.')
        ENDIF
      ENDIF
    ENDIF
END SUBROUTINE calc_ts_wtp_initialize
!
SUBROUTINE calc_ts_wtp_advance( grid )
  USE module_domain
  USE module_configure, ONLY : model_config_rec, grid_config_rec_type, model_to_grid_config_rec
  USE module_model_constants
  IMPLICIT NONE
    ! Arguments
    TYPE(domain), INTENT(INOUT) :: grid
    LOGICAL, EXTERNAL :: wrf_dm_on_monitor
    ! Local variables
    INTEGER :: i, k, l, m, n
    REAL :: xtime_minutes
    TYPE(grid_config_rec_type) :: config_flags
      REAL, ALLOCATABLE, DIMENSION(:, :) :: v1_ave_temp, vtheta_ave_temp, vt_ave_temp
      REAL, ALLOCATABLE, DIMENSION(:, :) :: vrel_ave_temp, vz_ave_temp
      REAL, ALLOCATABLE, DIMENSION(:, :) :: u_ave_temp, v_ave_temp, w_ave_temp
      REAL, ALLOCATABLE, DIMENSION(:, :) :: phi_ave_temp, alpha_ave_temp, cl_ave_temp
      REAL, ALLOCATABLE, DIMENSION(:, :) :: cd_ave_temp, cm_ave_temp, L_ave_temp
      REAL, ALLOCATABLE, DIMENSION(:, :) :: D_ave_temp, fn_ave_temp, ft_ave_temp
      REAL, ALLOCATABLE, DIMENSION(:, :) :: F_ave_temp, bladePoints_x_ave_temp
      REAL, ALLOCATABLE, DIMENSION(:, :, :) :: v1_temp, vtheta_temp, vt_temp
      REAL, ALLOCATABLE, DIMENSION(:, :, :) :: vrel_temp, vz_temp
      REAL, ALLOCATABLE, DIMENSION(:, :, :) :: u_temp, v_temp, w_temp
      REAL, ALLOCATABLE, DIMENSION(:, :, :) :: phi_temp, alpha_temp, cl_temp
      REAL, ALLOCATABLE, DIMENSION(:, :, :) :: cd_temp, cm_temp, L_temp
      REAL, ALLOCATABLE, DIMENSION(:, :, :) :: D_temp, fn_temp, ft_temp
      REAL, ALLOCATABLE, DIMENSION(:, :, :) :: F_temp, bladePoints_x_temp
      REAL, ALLOCATABLE, DIMENSION(:, :, :) :: bladePoints_y_temp, bladePoints_z_temp
    CALL model_to_grid_config_rec(grid%id, model_config_rec, config_flags)
#if ((EM_CORE == 1) && (DA_CORE != 1))
    IF( grid%dfi_opt /= DFI_NODFI .AND. grid%dfi_stage /= DFI_FST ) RETURN
#endif
    n = grid%next_ts_time_wtp
    IF(config_flags%wind_wtp_opt .GT. 0 .AND. config_flags%wind_wtp_opt .LT. 4) THEN
#if ( EM_CORE == 1 )
      IF(config_flags%wind_wtp_opt .EQ. 3) THEN
        ALLOCATE( v1_ave_temp(config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements) )
        ALLOCATE( vtheta_ave_temp(config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements) )
        ALLOCATE( vt_ave_temp(config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements) )
        ALLOCATE( vrel_ave_temp(config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements) )
        ALLOCATE( vz_ave_temp(config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements) )
        ALLOCATE( u_ave_temp(config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements) )
        ALLOCATE( v_ave_temp(config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements) )
        ALLOCATE( w_ave_temp(config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements) )
        ALLOCATE( phi_ave_temp(config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements) )
        ALLOCATE( alpha_ave_temp(config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements) )
        ALLOCATE( cl_ave_temp(config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements) )
        ALLOCATE( cd_ave_temp(config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements) )
        ALLOCATE( cm_ave_temp(config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements) )
        ALLOCATE( L_ave_temp(config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements) )
        ALLOCATE( D_ave_temp(config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements) )
        ALLOCATE( fn_ave_temp(config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements) )
        ALLOCATE( ft_ave_temp(config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements) )
        ALLOCATE( F_ave_temp(config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements) )
                                     
        v1_ave_temp = RESHAPE( grid%wtp_v1_ave, (/ config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements /) )
        vtheta_ave_temp = RESHAPE( grid%wtp_vtheta_ave, (/ config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements /) )
        vt_ave_temp = RESHAPE( grid%wtp_vt_ave, (/ config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements /) )
        vrel_ave_temp = RESHAPE( grid%wtp_vrel_ave, (/ config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements /) )
        vz_ave_temp = RESHAPE( grid%wtp_vz_ave, (/ config_flags%wind_wtp_nTurbines,config_flags%wind_wtp_nElements /) )
        u_ave_temp = RESHAPE( grid%wtp_u_ave, (/ config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements /) )
        v_ave_temp = RESHAPE( grid%wtp_v_ave, (/ config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements /) )
        w_ave_temp = RESHAPE( grid%wtp_w_ave, (/ config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements /) )
        phi_ave_temp = RESHAPE( grid%wtp_phi_ave, (/ config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements /) )
        alpha_ave_temp = RESHAPE( grid%wtp_alpha_ave, (/ config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements /) )
        cl_ave_temp = RESHAPE( grid%wtp_cl_ave, (/ config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements /) )
        cd_ave_temp = RESHAPE( grid%wtp_cd_ave, (/ config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements /) )
        cm_ave_temp = RESHAPE( grid%wtp_cm_ave, (/ config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements /) )
        L_ave_temp = RESHAPE( grid%wtp_L_ave, (/ config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements /) )
        D_ave_temp = RESHAPE( grid%wtp_D_ave, (/ config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements /) )
        fn_ave_temp = RESHAPE( grid%wtp_fn_ave, (/ config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements /) )
        ft_ave_temp = RESHAPE( grid%wtp_ft_ave, (/ config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements /) )
        F_ave_temp = RESHAPE( grid%wtp_F_ave, (/ config_flags%wind_wtp_nTurbines, config_flags%wind_wtp_nElements /) )
      ELSE
        ALLOCATE( v1_temp(config_flags%wind_wtp_nTurbines, &
                          config_flags%wind_wtp_nSections, &
                          config_flags%wind_wtp_nElements) )
        ALLOCATE( vtheta_temp(config_flags%wind_wtp_nTurbines, &
                              config_flags%wind_wtp_nSections, &
                              config_flags%wind_wtp_nElements) )
        ALLOCATE( vt_temp(config_flags%wind_wtp_nTurbines, &
                          config_flags%wind_wtp_nSections, &
                          config_flags%wind_wtp_nElements) )
        ALLOCATE( vrel_temp(config_flags%wind_wtp_nTurbines, &
                            config_flags%wind_wtp_nSections, &
                            config_flags%wind_wtp_nElements) )
        ALLOCATE( vz_temp(config_flags%wind_wtp_nTurbines, &
                          config_flags%wind_wtp_nSections, &
                          config_flags%wind_wtp_nElements) )
        ALLOCATE( u_temp(config_flags%wind_wtp_nTurbines, &
                         config_flags%wind_wtp_nSections, &
                         config_flags%wind_wtp_nElements) )
        ALLOCATE( v_temp(config_flags%wind_wtp_nTurbines, &
                         config_flags%wind_wtp_nSections, &
                         config_flags%wind_wtp_nElements) )
        ALLOCATE( w_temp(config_flags%wind_wtp_nTurbines, &
                         config_flags%wind_wtp_nSections, &
                         config_flags%wind_wtp_nElements) )
        ALLOCATE( phi_temp(config_flags%wind_wtp_nTurbines, &
                           config_flags%wind_wtp_nSections, &
                           config_flags%wind_wtp_nElements) )
        ALLOCATE( alpha_temp(config_flags%wind_wtp_nTurbines, &
                             config_flags%wind_wtp_nSections, &
                             config_flags%wind_wtp_nElements) )
        ALLOCATE( cl_temp(config_flags%wind_wtp_nTurbines, &
                          config_flags%wind_wtp_nSections, &
                          config_flags%wind_wtp_nElements) )
        ALLOCATE( cd_temp(config_flags%wind_wtp_nTurbines, &
                          config_flags%wind_wtp_nSections, &
                          config_flags%wind_wtp_nElements) )
        ALLOCATE( cm_temp(config_flags%wind_wtp_nTurbines, &
                          config_flags%wind_wtp_nSections, &
                          config_flags%wind_wtp_nElements) )
        ALLOCATE( L_temp(config_flags%wind_wtp_nTurbines, &
                         config_flags%wind_wtp_nSections, &
                         config_flags%wind_wtp_nElements) )
        ALLOCATE( D_temp(config_flags%wind_wtp_nTurbines, &
                         config_flags%wind_wtp_nSections, &
                         config_flags%wind_wtp_nElements) )
        ALLOCATE( fn_temp(config_flags%wind_wtp_nTurbines, &
                          config_flags%wind_wtp_nSections, &
                          config_flags%wind_wtp_nElements) )
        ALLOCATE( ft_temp(config_flags%wind_wtp_nTurbines, &
                          config_flags%wind_wtp_nSections, &
                          config_flags%wind_wtp_nElements) )
        ALLOCATE( F_temp(config_flags%wind_wtp_nTurbines, &
                         config_flags%wind_wtp_nSections, &
                         config_flags%wind_wtp_nElements) )
        ALLOCATE( bladePoints_x_temp(config_flags%wind_wtp_nTurbines, &
                                     config_flags%wind_wtp_nSections, &
                                     config_flags%wind_wtp_nElements) )
        ALLOCATE( bladePoints_y_temp(config_flags%wind_wtp_nTurbines, &
                                     config_flags%wind_wtp_nSections, &
                                     config_flags%wind_wtp_nElements) )
        ALLOCATE( bladePoints_z_temp(config_flags%wind_wtp_nTurbines, &
                                     config_flags%wind_wtp_nSections, &
                                     config_flags%wind_wtp_nElements) )
                                     
        v1_temp = RESHAPE( grid%wtp_v1, (/ config_flags%wind_wtp_nTurbines, &
                                           config_flags%wind_wtp_nSections, &
                                           config_flags%wind_wtp_nElements /) )
        vtheta_temp = RESHAPE( grid%wtp_vtheta, (/ config_flags%wind_wtp_nTurbines, &
                                                   config_flags%wind_wtp_nSections, &
                                                   config_flags%wind_wtp_nElements /) )
        vt_temp = RESHAPE( grid%wtp_vt, (/ config_flags%wind_wtp_nTurbines, &
                                           config_flags%wind_wtp_nSections, &
                                           config_flags%wind_wtp_nElements /) )
        vrel_temp = RESHAPE( grid%wtp_vrel, (/ config_flags%wind_wtp_nTurbines, &
                                               config_flags%wind_wtp_nSections, &
                                               config_flags%wind_wtp_nElements /) )
        vz_temp = RESHAPE( grid%wtp_vz, (/ config_flags%wind_wtp_nTurbines, &
                                           config_flags%wind_wtp_nSections, &
                                           config_flags%wind_wtp_nElements /) )
        u_temp = RESHAPE( grid%wtp_u, (/ config_flags%wind_wtp_nTurbines, &
                                         config_flags%wind_wtp_nSections, &
                                         config_flags%wind_wtp_nElements /) )
        v_temp = RESHAPE( grid%wtp_v, (/ config_flags%wind_wtp_nTurbines, &
                                         config_flags%wind_wtp_nSections, &
                                         config_flags%wind_wtp_nElements /) )
        w_temp = RESHAPE( grid%wtp_w, (/ config_flags%wind_wtp_nTurbines, &
                                         config_flags%wind_wtp_nSections, &
                                         config_flags%wind_wtp_nElements /) )
        phi_temp = RESHAPE( grid%wtp_phi, (/ config_flags%wind_wtp_nTurbines, &
                                             config_flags%wind_wtp_nSections, &
                                             config_flags%wind_wtp_nElements /) )
        alpha_temp = RESHAPE( grid%wtp_alpha, (/ config_flags%wind_wtp_nTurbines, &
                                                 config_flags%wind_wtp_nSections, &
                                                 config_flags%wind_wtp_nElements /) )
        cl_temp = RESHAPE( grid%wtp_cl, (/ config_flags%wind_wtp_nTurbines, &
                                           config_flags%wind_wtp_nSections, &
                                           config_flags%wind_wtp_nElements /) )
        cd_temp = RESHAPE( grid%wtp_cd, (/ config_flags%wind_wtp_nTurbines, &
                                           config_flags%wind_wtp_nSections, &
                                           config_flags%wind_wtp_nElements /) )
        cm_temp = RESHAPE( grid%wtp_cm, (/ config_flags%wind_wtp_nTurbines, &
                                           config_flags%wind_wtp_nSections, &
                                           config_flags%wind_wtp_nElements /) )
        L_temp = RESHAPE( grid%wtp_L, (/ config_flags%wind_wtp_nTurbines, &
                                         config_flags%wind_wtp_nSections, &
                                         config_flags%wind_wtp_nElements /) )
        D_temp = RESHAPE( grid%wtp_D, (/ config_flags%wind_wtp_nTurbines, &
                                         config_flags%wind_wtp_nSections, &
                                         config_flags%wind_wtp_nElements /) )
        fn_temp = RESHAPE( grid%wtp_fn, (/ config_flags%wind_wtp_nTurbines, &
                                           config_flags%wind_wtp_nSections, &
                                           config_flags%wind_wtp_nElements /) )
        ft_temp = RESHAPE( grid%wtp_ft, (/ config_flags%wind_wtp_nTurbines, &
                                           config_flags%wind_wtp_nSections, &
                                           config_flags%wind_wtp_nElements /) )
        F_temp = RESHAPE( grid%wtp_F, (/ config_flags%wind_wtp_nTurbines, &
                                         config_flags%wind_wtp_nSections, &
                                         config_flags%wind_wtp_nElements /) )
        bladePoints_x_temp = RESHAPE( grid%wtp_bladePoints_x, (/ config_flags%wind_wtp_nTurbines, &
                                                                 config_flags%wind_wtp_nSections, &
                                                                 config_flags%wind_wtp_nElements /) )
        bladePoints_y_temp = RESHAPE( grid%wtp_bladePoints_y, (/ config_flags%wind_wtp_nTurbines, &
                                                                 config_flags%wind_wtp_nSections, &
                                                                 config_flags%wind_wtp_nElements /) )
        bladePoints_z_temp = RESHAPE( grid%wtp_bladePoints_z, (/ config_flags%wind_wtp_nTurbines, &
                                                                 config_flags%wind_wtp_nSections, &
                                                                 config_flags%wind_wtp_nElements /) )
      ENDIF
#endif
      CALL domain_clock_get( grid, minutesSinceSimulationStart=xtime_minutes )
      grid%wtp_ts_sec(n) = xtime_minutes*60.
      grid%wtp_ts_dt(n) = grid%dt
      grid%wtp_ts_timestep(n) = grid%itimestep
      DO k = 1, config_flags%wind_wtp_nTurbines
        grid%wtp_ts_rho_fst_ave(n,k) = grid%wtp_rho_fst_ave(k)
        grid%wtp_ts_tsr(n,k) = grid%wtp_tsr(k)
        grid%wtp_ts_omega(n,k) = grid%wtp_omega(k)*(30.0/piconst)
        grid%wtp_ts_pitch(n,k) = grid%wtp_pitch(k)
        grid%wtp_ts_v0_fst_ave(n,k) = grid%wtp_v0_fst_ave(k)
        grid%wtp_ts_wind_dir_fst_ave(n,k) = grid%wtp_wind_dir_fst_ave(k)
        grid%wtp_ts_turbine_dir(n,k) = grid%wtp_turbine_dir(k)
        grid%wtp_ts_capital_phi(n,k) = grid%wtp_capital_phi(k)
        grid%wtp_ts_wsRotor(n,k) = grid%wtp_wsRotor(k)
        grid%wtp_ts_thrust(n,k) = grid%wtp_thrust(k)
        grid%wtp_ts_torque(n,k) = grid%wtp_torque(k)
        grid%wtp_ts_power(n,k) = grid%wtp_power(k)
        grid%wtp_ts_trqGenCommandedOld(n,k) = grid%wtp_trqGenCommandedOld(k)
        grid%wtp_ts_power_mech(n,k) = grid%wtp_power_mech(k)
        grid%wtp_ts_power_gen(n,k) = grid%wtp_power_gen(k)
        grid%wtp_ts_thrust_coeff(n,k) = grid%wtp_thrust_coeff(k)
        grid%wtp_ts_torque_coeff(n,k) = grid%wtp_torque_coeff(k)
        grid%wtp_ts_power_coeff(n,k) = grid%wtp_power_coeff(k)
        grid%wtp_ts_rotorApex_x(n,k) = grid%wtp_rotorApex_x(k)
        grid%wtp_ts_rotorApex_y(n,k) = grid%wtp_rotorApex_y(k)
        grid%wtp_ts_rotorApex_z(n,k) = grid%wtp_rotorApex_z(k)
        IF(config_flags%wind_wtp_opt .EQ. 3) THEN
          DO m = 1, config_flags%wind_wtp_nElements
#if (EM_CORE == 1)
            grid%wtp_ts_v1_ave(n,k,m) = v1_ave_temp(k,m)
            grid%wtp_ts_vtheta_ave(n,k,m) = vtheta_ave_temp(k,m)
            grid%wtp_ts_vt_ave(n,k,m) = vt_ave_temp(k,m)
            grid%wtp_ts_vrel_ave(n,k,m) = vrel_ave_temp(k,m)
            grid%wtp_ts_vz_ave(n,k,m) = vz_ave_temp(k,m)
            grid%wtp_ts_u_ave(n,k,m) = u_ave_temp(k,m)
            grid%wtp_ts_v_ave(n,k,m) = v_ave_temp(k,m)
            grid%wtp_ts_w_ave(n,k,m) = w_ave_temp(k,m)
            grid%wtp_ts_phi_ave(n,k,m) = phi_ave_temp(k,m)
            grid%wtp_ts_alpha_ave(n,k,m) = alpha_ave_temp(k,m)
            grid%wtp_ts_cl_ave(n,k,m) = cl_ave_temp(k,m)
            grid%wtp_ts_cd_ave(n,k,m) = cd_ave_temp(k,m)
            grid%wtp_ts_cm_ave(n,k,m) = cm_ave_temp(k,m)
            grid%wtp_ts_L_ave(n,k,m) = L_ave_temp(k,m)
            grid%wtp_ts_D_ave(n,k,m) = D_ave_temp(k,m)
            grid%wtp_ts_fn_ave(n,k,m) = fn_ave_temp(k,m)
            grid%wtp_ts_ft_ave(n,k,m) = ft_ave_temp(k,m)
            grid%wtp_ts_F_ave(n,k,m) = F_ave_temp(k,m)
#endif
          ENDDO
        ELSE
          DO l = 1, config_flags%wind_wtp_nSections
            DO m = 1, config_flags%wind_wtp_nElements
#if (EM_CORE == 1)
              grid%wtp_ts_v1(n,k,l,m) = v1_temp(k,l,m)
              grid%wtp_ts_vtheta(n,k,l,m) = vtheta_temp(k,l,m)
              grid%wtp_ts_vt(n,k,l,m) = vt_temp(k,l,m)
              grid%wtp_ts_vrel(n,k,l,m) = vrel_temp(k,l,m)
              grid%wtp_ts_vz(n,k,l,m) = vz_temp(k,l,m)
              grid%wtp_ts_u(n,k,l,m) = u_temp(k,l,m)
              grid%wtp_ts_v(n,k,l,m) = v_temp(k,l,m)
              grid%wtp_ts_w(n,k,l,m) = w_temp(k,l,m)
              grid%wtp_ts_phi(n,k,l,m) = phi_temp(k,l,m)
              grid%wtp_ts_alpha(n,k,l,m) = alpha_temp(k,l,m)
              grid%wtp_ts_cl(n,k,l,m) = cl_temp(k,l,m)
              grid%wtp_ts_cd(n,k,l,m) = cd_temp(k,l,m)
              grid%wtp_ts_cm(n,k,l,m) = cm_temp(k,l,m)
              grid%wtp_ts_L(n,k,l,m) = L_temp(k,l,m)
              grid%wtp_ts_D(n,k,l,m) = D_temp(k,l,m)
              grid%wtp_ts_fn(n,k,l,m) = fn_temp(k,l,m)
              grid%wtp_ts_ft(n,k,l,m) = ft_temp(k,l,m)
              grid%wtp_ts_F(n,k,l,m) = F_temp(k,l,m)
              grid%wtp_ts_bladePoints_x(n,k,l,m) = bladePoints_x_temp(k,l,m)
              grid%wtp_ts_bladePoints_y(n,k,l,m) = bladePoints_y_temp(k,l,m)
              grid%wtp_ts_bladePoints_z(n,k,l,m) = bladePoints_z_temp(k,l,m)
#endif
            ENDDO
          ENDDO
        ENDIF
      ENDDO
#if (EM_CORE == 1)
      IF(config_flags%wind_wtp_opt .EQ. 3) THEN
        DEALLOCATE(v1_ave_temp, vtheta_ave_temp, vt_ave_temp, vrel_ave_temp, vz_ave_temp)
        DEALLOCATE(phi_ave_temp, u_ave_temp, v_ave_temp, w_ave_temp)
        DEALLOCATE(alpha_ave_temp, cl_ave_temp, cd_ave_temp, cm_ave_temp)
        DEALLOCATE(L_ave_temp, D_ave_temp, fn_ave_temp, ft_ave_temp, F_ave_temp)
      ELSE
        DEALLOCATE(v1_temp, vtheta_temp, vt_temp, vrel_temp, vz_temp)
        DEALLOCATE(phi_temp, u_temp, v_temp, w_temp)
        DEALLOCATE(alpha_temp, cl_temp, cd_temp, cm_temp)
        DEALLOCATE(L_temp, D_temp, fn_temp, ft_temp, F_temp)
        DEALLOCATE(bladePoints_x_temp)
        DEALLOCATE(bladePoints_y_temp)
        DEALLOCATE(bladePoints_z_temp)
      ENDIF
#endif
      grid%next_ts_time_wtp = grid%next_ts_time_wtp + 1
      IF(grid%next_ts_time_wtp .GT. config_flags%wind_wtp_ts_buf_size) CALL write_ts_wtp(grid)
    ENDIF
END SUBROUTINE calc_ts_wtp_advance
!!
SUBROUTINE write_ts_wtp( grid )
  USE module_domain, ONLY : domain
  USE module_dm, ONLY : wrf_dm_min_reals
  USE module_state_description
  USE module_configure, ONLY : model_config_rec, grid_config_rec_type, model_to_grid_config_rec
  IMPLICIT NONE
    ! Arguments
    TYPE (domain), INTENT(INOUT) :: grid
    TYPE(grid_config_rec_type) :: config_flags
    LOGICAL, EXTERNAL :: wrf_dm_on_monitor
    INTEGER, EXTERNAL :: get_unused_unit
    ! Local variables
    INTEGER :: i, j, k, l, m, n, iunit
    CHARACTER (LEN=35) :: ts_filename_wtp
    CALL model_to_grid_config_rec(grid%id, model_config_rec, config_flags)
    IF(config_flags%wind_wtp_opt .GT. 0 .AND. config_flags%wind_wtp_opt .LT. 4) THEN
#if ((EM_CORE == 1) && (DA_CORE != 1))
      IF( grid%dfi_opt /= DFI_NODFI .AND. grid%dfi_stage /= DFI_FST ) RETURN
#endif
      IF( wrf_dm_on_monitor() ) THEN
        DO k = 1, config_flags%wind_wtp_nTurbines
          WRITE(ts_filename_wtp,'(A)') 'turbine.0000'
          ! Create the output files for the vertical profiles, one file for each variable
          iunit = get_unused_unit()
          IF( iunit .LE. 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts_wtp: could not find a free Fortran unit.')
          ENDIF
          i = LEN_TRIM(ts_filename_wtp)
          IF(k .LT. 10  ) WRITE(ts_filename_wtp(i  :i),'(I1)') k
          IF(k .GE. 10  ) WRITE(ts_filename_wtp(i-1:i),'(I2)') k
          IF(k .GE. 100 ) WRITE(ts_filename_wtp(i-2:i),'(I3)') k
          IF(k .GE. 1000) WRITE(ts_filename_wtp(i-3:i),'(I4)') k
          WRITE(ts_filename_wtp,'(A)') 'tsout/'//TRIM(ts_filename_wtp)
          OPEN(UNIT=iunit, FILE=ts_filename_wtp, STATUS='unknown', POSITION='append', FORM='formatted')
          DO n = 1, grid%next_ts_time_wtp - 1
#if (EM_CORE == 1)
            WRITE(UNIT=iunit,FMT='(I11,3x,f13.6,3x,f9.6,3x,f8.3,3x,                     &
            f6.3,3x,f9.3,3x,f8.3,3x,f13.3,3x,f14.3,3x,f10.3,3x,f9.3,3x,f11.3,3x,        &
            f12.3,3x,f12.3,3x,f12.3,3x,f12.3,3x,f12.3,3x,f12.3,3x,f9.3,3x,              &
            f9.3,3x,f9.3,3x,f13.3,3x,f13.3,3x,f13.3)')                                  &
            grid%wtp_ts_timestep(n), grid%wtp_ts_sec(n), grid%wtp_ts_dt(n),             &
            grid%wtp_ts_rho_fst_ave(n,k), grid%wtp_ts_tsr(n,k),                         &
            grid%wtp_ts_omega(n,k), grid%wtp_ts_pitch(n,k),                             &
            grid%wtp_ts_v0_fst_ave(n,k), grid%wtp_ts_wind_dir_fst_ave(n,k),             &
            grid%wtp_ts_turbine_dir(n,k), grid%wtp_ts_capital_phi(n,k),                 &
            grid%wtp_ts_wsRotor(n,k), grid%wtp_ts_thrust(n,k), grid%wtp_ts_torque(n,k), &
            grid%wtp_ts_trqGenCommandedOld(n,k), grid%wtp_ts_power(n,k),                &
            grid%wtp_ts_power_mech(n,k), grid%wtp_ts_power_gen(n,k),                    &
            grid%wtp_ts_thrust_coeff(n,k), grid%wtp_ts_torque_coeff(n,k),               &
            grid%wtp_ts_power_coeff(n,k), grid%wtp_ts_rotorApex_x(n,k),                 &
            grid%wtp_ts_rotorApex_y(n,k), grid%wtp_ts_rotorApex_z(n,k)
#endif
          ENDDO
          CLOSE(UNIT=iunit)
          IF(config_flags%wind_wtp_opt .EQ. 3) THEN
            WRITE(ts_filename_wtp,'(A)') 'turbine.0000.disk'
            ! Create the output files for the vertical profiles, one file for each variable
            iunit = get_unused_unit()
            IF( iunit .LE. 0 ) THEN
              CALL wrf_error_fatal('Error in write_ts_wtp: could not find a free Fortran unit.')
            ENDIF
            i = LEN_TRIM(ts_filename_wtp)
            IF(k .LT. 10  ) WRITE(ts_filename_wtp(i-5:i-5),'(I1)') k
            IF(k .GE. 10  ) WRITE(ts_filename_wtp(i-6:i-5),'(I2)') k
            IF(k .GE. 100 ) WRITE(ts_filename_wtp(i-7:i-6),'(I3)') k
            IF(k .GE. 1000) WRITE(ts_filename_wtp(i-8:i-7),'(I4)') k
            WRITE(ts_filename_wtp,'(A)') 'tsout/'//TRIM(ts_filename_wtp)
            OPEN(UNIT=iunit, FILE=ts_filename_wtp, STATUS='unknown', POSITION='append', FORM='formatted')
            DO n = 1, grid%next_ts_time_wtp - 1
#if (EM_CORE == 1)
              WRITE(UNIT=iunit,FMT='(I11,3x,f13.6,3x,f9.6,3x)')  &
                         grid%wtp_ts_timestep(n), grid%wtp_ts_sec(n), grid%wtp_ts_dt(n)
#endif
#if (EM_CORE == 1)
              DO m = 1, config_flags%wind_wtp_nElements
                WRITE(UNIT=iunit,FMT='(I50,3x,f10.6,3x,f14.6,3x,f10.6,3x,f10.6,3x,     &
                f10.6,3x,f10.6,3x,f10.6,3x,f10.6,3x,f11.6,3x,f11.6,3x,f9.6,3x,f9.6,3x, &
                f9.6,3x,f12.6,3x,f11.6,3x,f12.6,3x,f12.6,3x,f8.6,3x)')                 &
                m, grid%wtp_ts_v1_ave(n,k,m), grid%wtp_ts_vtheta_ave(n,k,m),           &
                grid%wtp_ts_vt_ave(n,k,m), grid%wtp_ts_vrel_ave(n,k,m),                &
                grid%wtp_ts_vz_ave(n,k,m), grid%wtp_ts_u_ave(n,k,m),                   &
                grid%wtp_ts_v_ave(n,k,m), grid%wtp_ts_w_ave(n,k,m),                    &
                grid%wtp_ts_phi_ave(n,k,m), grid%wtp_ts_alpha_ave(n,k,m),              &
                grid%wtp_ts_cl_ave(n,k,m), grid%wtp_ts_cd_ave(n,k,m),                  &
                grid%wtp_ts_cm_ave(n,k,m), grid%wtp_ts_L_ave(n,k,m),                   &
                grid%wtp_ts_D_ave(n,k,m), grid%wtp_ts_ft_ave(n,k,m),                   &
                grid%wtp_ts_fn_ave(n,k,m), grid%wtp_ts_F_ave(n,k,m)
#endif
              ENDDO
            ENDDO
            CLOSE(UNIT=iunit)
          ELSE
            WRITE(ts_filename_wtp,'(A)') 'turbine.0000.blade.000'
            ! Create the output files for the vertical profiles, one file for each variable
            iunit = get_unused_unit()
            IF( iunit .LE. 0 ) THEN
              CALL wrf_error_fatal('Error in write_ts_wtp: could not find a free Fortran unit.')
            ENDIF
            i = LEN_TRIM(ts_filename_wtp)
            IF(k .LT. 10  ) WRITE(ts_filename_wtp(i-10:i-10),'(I1)') k
            IF(k .GE. 10  ) WRITE(ts_filename_wtp(i-11:i-10),'(I2)') k
            IF(k .GE. 100 ) WRITE(ts_filename_wtp(i-12:i-11),'(I3)') k
            IF(k .GE. 1000) WRITE(ts_filename_wtp(i-13:i-12),'(I4)') k
            DO l = 1, config_flags%wind_wtp_nSections
              WRITE(ts_filename_wtp,'(A)') 'turbine.0000.blade.000'
              i = LEN_TRIM(ts_filename_wtp)
              IF(k .LT. 10  ) WRITE(ts_filename_wtp(i-10:i-10),'(I1)') k
              IF(k .GE. 10  ) WRITE(ts_filename_wtp(i-11:i-10),'(I2)') k
              IF(k .GE. 100 ) WRITE(ts_filename_wtp(i-12:i-11),'(I3)') k
              IF(k .GE. 1000) WRITE(ts_filename_wtp(i-13:i-12),'(I4)') k
              IF(l .LT. 10  ) WRITE(ts_filename_wtp(i  :i),'(I1)') l
              IF(l .GE. 10  ) WRITE(ts_filename_wtp(i-1:i),'(I2)') l
              IF(l .GE. 100 ) WRITE(ts_filename_wtp(i-2:i),'(I3)') l
              WRITE(ts_filename_wtp,'(A)') 'tsout/'//TRIM(ts_filename_wtp)
              OPEN(UNIT=iunit, FILE=ts_filename_wtp, STATUS='unknown', POSITION='append', FORM='formatted')
              DO n = 1, grid%next_ts_time_wtp - 1
#if (EM_CORE == 1)
                WRITE(UNIT=iunit,FMT='(I11,3x,f13.6,3x,f9.6,3x)')  &
                                  grid%wtp_ts_timestep(n), grid%wtp_ts_sec(n), grid%wtp_ts_dt(n)
#endif
#if (EM_CORE == 1)
                DO m = 1, config_flags%wind_wtp_nElements
                  WRITE(UNIT=iunit,FMT='(I50,3x,f10.6,3x,f14.6,3x,f10.6,3x,f10.6,3x,     &
                  f10.6,3x,f10.6,3x,f10.6,3x,f10.6,3x,f11.6,3x,f11.6,3x,f9.6,3x,f9.6,3x, &
                  f9.6,3x,f12.6,3x,f11.6,3x,f12.6,3x,f12.6,3x,f8.6,3x,f12.6,3x,f12.6,3x, &
                  f12.6,3x)')                                                            &
                  m, grid%wtp_ts_v1(n,k,l,m), grid%wtp_ts_vtheta(n,k,l,m),               &
                  grid%wtp_ts_vt(n,k,l,m), grid%wtp_ts_vrel(n,k,l,m),                    &
                  grid%wtp_ts_vz(n,k,l,m), grid%wtp_ts_u(n,k,l,m),                       &
                  grid%wtp_ts_v(n,k,l,m), grid%wtp_ts_w(n,k,l,m),                        &
                  grid%wtp_ts_phi(n,k,l,m), grid%wtp_ts_alpha(n,k,l,m),                  &
                  grid%wtp_ts_cl(n,k,l,m), grid%wtp_ts_cd(n,k,l,m),                      &
                  grid%wtp_ts_cm(n,k,l,m), grid%wtp_ts_L(n,k,l,m),                       &
                  grid%wtp_ts_D(n,k,l,m), grid%wtp_ts_ft(n,k,l,m),                       &
                  grid%wtp_ts_fn(n,k,l,m), grid%wtp_ts_F(n,k,l,m),                       &
                  grid%wtp_ts_bladePoints_x(n,k,l,m),                                    &
                  grid%wtp_ts_bladePoints_y(n,k,l,m),                                    &
                  grid%wtp_ts_bladePoints_z(n,k,l,m)
#endif
                ENDDO
              ENDDO
              CLOSE(UNIT=iunit)
            ENDDO
          ENDIF
        ENDDO
      ENDIF
    grid%next_ts_time_wtp = 1
    ENDIF
END SUBROUTINE write_ts_wtp
!!
SUBROUTINE wind_wtp_tsout_dir(exists)
   IMPLICIT NONE
   LOGICAL , INTENT(OUT) :: exists
   INQUIRE(FILE='./tsout/', EXIST=exists)
END SUBROUTINE wind_wtp_tsout_dir

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine prints out the current value of variables at all specified
! time series locations that are within the current patch.
!
! Baris Kale, 01/12/2023 (reproduced from wrf_timeseries.F)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calc_ts_wtp_sad_initialize( grid )
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
    CHARACTER (LEN=35) :: ts_filename_wtp_sad
    CHARACTER (LEN=19) simulation_start_date
    INTEGER :: simulation_start_year   , &
               simulation_start_month  , &
               simulation_start_day    , &
               simulation_start_hour   , &
               simulation_start_minute , &
               simulation_start_second
    CALL model_to_grid_config_rec(grid%id, model_config_rec, config_flags)
    IF(config_flags%wind_wtp_opt .EQ. 4) THEN
#if ((EM_CORE == 1) && (DA_CORE != 1))
      IF( grid%dfi_stage == DFI_FST ) THEN
#endif
        grid%next_ts_time_wtp_sad = 1
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
        WRITE(message, '(A51,I3,A15,A19)') 'Computing time series for the SAD model for domain ', grid%id, &
        ' starting from ',simulation_start_date
        CALL wrf_message(message)
        CALL wind_wtp_sad_tsout_dir(exists)
        IF(exists) THEN
          IF( wrf_dm_on_monitor() ) THEN
            DO j = 1, config_flags%wind_wtp_nTurbines
            WRITE(ts_filename_wtp_sad,'(A)') 'turbine.0000'
              ! Create the output files for the vertical profiles, one file for each variable
              iunit = get_unused_unit()
              IF( iunit .LE. 0 ) THEN
                CALL wrf_error_fatal('Error in calc_ts_wtp_sad_initialize: could not find a free Fortran unit.')
              ENDIF
              i = LEN_TRIM(ts_filename_wtp_sad)
              IF(j .LT. 10  ) WRITE(ts_filename_wtp_sad(i  :i),'(I1)') j
              IF(j .GE. 10  ) WRITE(ts_filename_wtp_sad(i-1:i),'(I2)') j
              IF(j .GE. 100 ) WRITE(ts_filename_wtp_sad(i-2:i),'(I3)') j
              IF(j .GE. 1000) WRITE(ts_filename_wtp_sad(i-3:i),'(I4)') j
              WRITE(ts_filename_wtp_sad,'(A)') 'tsout/'//TRIM(ts_filename_wtp_sad)
              OPEN(UNIT=iunit, FILE=TRIM(ts_filename_wtp_sad), FORM='FORMATTED', STATUS='REPLACE')
#if (EM_CORE == 1)
                WRITE(UNIT=iunit, FMT='(A11,A16,A12,A11,A16,A17,A13,A15,A12,A12,A16)')      &
                'timestep[-]', 'time[s]', 'dt[s]', 'rho[kg/m3]', 'fstWindSpd[m/s]',         &
                'fstWindDir[deg]', 'turbDir[deg]', 'rotPwr[W]', 'thrCoef[-]', 'pwrCoef[-]', &
                'axialIndFact[-]'
#endif
              CLOSE(UNIT=iunit)
            ENDDO
          ENDIF
        ELSE
          CALL wrf_error_fatal('Error in calc_ts_wtp_sad_initialize: No tsout directory is found. &
                                Create "./tsout/" directory to save time series.')
        ENDIF
      ENDIF
    ENDIF
END SUBROUTINE calc_ts_wtp_sad_initialize
!
SUBROUTINE calc_ts_wtp_sad_advance( grid )
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
    CALL model_to_grid_config_rec(grid%id, model_config_rec, config_flags)
#if ((EM_CORE == 1) && (DA_CORE != 1))
    IF( grid%dfi_opt /= DFI_NODFI .AND. grid%dfi_stage /= DFI_FST ) RETURN
#endif
    n = grid%next_ts_time_wtp_sad
    IF(config_flags%wind_wtp_opt .EQ. 4) THEN
      CALL domain_clock_get( grid, minutesSinceSimulationStart=xtime_minutes )
      grid%wtp_ts_sec(n) = xtime_minutes*60.
      grid%wtp_ts_dt(n) = grid%dt
      grid%wtp_ts_timestep(n) = grid%itimestep
      DO k = 1, config_flags%wind_wtp_nTurbines
        grid%wtp_ts_rho_fst_ave(n,k) = grid%wtp_rho_fst_ave(k)
        grid%wtp_ts_v0_fst_ave(n,k) = grid%wtp_v0_fst_ave(k)
        grid%wtp_ts_wind_dir_fst_ave(n,k) = grid%wtp_wind_dir_fst_ave(k)
        grid%wtp_ts_turbine_dir(n,k) = grid%wtp_turbine_dir(k)
        grid%wtp_ts_power(n,k) = grid%wtp_power(k)
        grid%wtp_ts_thrust_coeff(n,k) = grid%wtp_thrust_coeff(k)
        grid%wtp_ts_power_coeff(n,k) = grid%wtp_power_coeff(k)
        grid%wtp_ts_aif(n,k) = grid%wtp_aif(k)
      ENDDO
      grid%next_ts_time_wtp_sad = grid%next_ts_time_wtp_sad + 1
      IF( grid%next_ts_time_wtp_sad .GT. config_flags%wind_wtp_ts_buf_size ) CALL write_ts_wtp_sad(grid)
    ENDIF
END SUBROUTINE calc_ts_wtp_sad_advance
!!
SUBROUTINE write_ts_wtp_sad( grid )
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
    CHARACTER (LEN=35) :: ts_filename_wtp_sad
    CALL model_to_grid_config_rec(grid%id, model_config_rec, config_flags)
    IF(config_flags%wind_wtp_opt .EQ. 4) THEN
#if ((EM_CORE == 1) && (DA_CORE != 1))
      IF( grid%dfi_opt /= DFI_NODFI .AND. grid%dfi_stage /= DFI_FST ) RETURN
#endif
      IF( wrf_dm_on_monitor() ) THEN
        DO k = 1, config_flags%wind_wtp_nTurbines
          WRITE(ts_filename_wtp_sad,'(A)') 'turbine.0000'
          ! Create the output files for the vertical profiles, one file for each variable
          iunit = get_unused_unit()
          IF( iunit .LE. 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts_wtp_sad: could not find a free Fortran unit.')
          ENDIF
          i = LEN_TRIM(ts_filename_wtp_sad)
          IF(k .LT. 10  ) WRITE(ts_filename_wtp_sad(i  :i),'(I1)') k
          IF(k .GE. 10  ) WRITE(ts_filename_wtp_sad(i-1:i),'(I2)') k
          IF(k .GE. 100 ) WRITE(ts_filename_wtp_sad(i-2:i),'(I3)') k
          IF(k .GE. 1000) WRITE(ts_filename_wtp_sad(i-3:i),'(I4)') k
          WRITE(ts_filename_wtp_sad,'(A)') 'tsout/'//TRIM(ts_filename_wtp_sad)
          OPEN(UNIT=iunit, FILE=ts_filename_wtp_sad, STATUS='unknown', POSITION='append', FORM='formatted')
          DO n = 1, grid%next_ts_time_wtp_sad - 1
#if (EM_CORE == 1)
            WRITE(UNIT=iunit,FMT='(I11,3x,f13.6,3x,f9.6,3x,f8.3,3x,          &
            f13.3,3x,f14.3,3x,f10.3,3x,f12.3,3x,f9.3,3x,f9.3,3x,f13.3,3x)')  &
            grid%wtp_ts_timestep(n), grid%wtp_ts_sec(n), grid%wtp_ts_dt(n),  &
            grid%wtp_ts_rho_fst_ave(n,k), grid%wtp_ts_v0_fst_ave(n,k),       &
            grid%wtp_ts_wind_dir_fst_ave(n,k), grid%wtp_ts_turbine_dir(n,k), &
            grid%wtp_ts_power(n,k), grid%wtp_ts_thrust_coeff(n,k),           &
            grid%wtp_ts_power_coeff(n,k), grid%wtp_ts_aif(n,k)
#endif
          ENDDO
          CLOSE(UNIT=iunit)
        ENDDO
      ENDIF
    grid%next_ts_time_wtp_sad = 1
    ENDIF
END SUBROUTINE write_ts_wtp_sad
!!
SUBROUTINE wind_wtp_sad_tsout_dir(exists)
   IMPLICIT NONE
   LOGICAL , INTENT(OUT) :: exists
   INQUIRE(FILE='./tsout/', EXIST=exists)
END SUBROUTINE wind_wtp_sad_tsout_dir

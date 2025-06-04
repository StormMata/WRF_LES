!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine prints out the current value of variables at all specified
! time series locations that are within the current patch.
!
! Baris Kale, 23/11/2023 (reproduced from wrf_timeseries_wtp.F)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calc_ts_wtp_pln_initialize( grid )
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
    INTEGER :: i, j, k, ij, iunit
    LOGICAL :: exists
    CHARACTER (LEN=132) :: message
    CHARACTER (LEN=35) :: ts_filename_wtp_pln
    CHARACTER (LEN=19) simulation_start_date
    INTEGER :: simulation_start_year   , &
               simulation_start_month  , &
               simulation_start_day    , &
               simulation_start_hour   , &
               simulation_start_minute , &
               simulation_start_second
    CALL model_to_grid_config_rec(grid%id, model_config_rec, config_flags)
    IF( config_flags%wind_wtp_opt .GT. 0 ) THEN
      IF( config_flags%wind_wtp_pln_tsout_opt ) THEN
#if ((EM_CORE == 1) && (DA_CORE != 1))
        IF( grid%dfi_stage == DFI_FST ) THEN
#endif
          grid%next_ts_time_wtp_pln = 1
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
          WRITE(message, '(A55,I3,A15,A19)') 'Computing time series of velocity components in &
                          a vertical plane for the GAD model for domain ', grid%id, &
          ' starting from ',simulation_start_date
          CALL wrf_message(message)
          CALL wind_wtp_pln_tsout_dir(exists)
          IF(exists) THEN
            IF( wrf_dm_on_monitor() ) THEN
              DO k = 1, config_flags%wind_wtp_pln_nTurbines
                j = model_config_rec%wind_wtp_pln_turbine_ids(k)
                WRITE(ts_filename_wtp_pln,'(A)') 'turbine.0000.plane.0000'
                ! Create the output files for the vertical profiles, one file for each variable
                iunit = get_unused_unit()
                IF( iunit .LE. 0 ) THEN
                  CALL wrf_error_fatal('Error in calc_ts_wtp_pln_initialize: could not find a free Fortran unit.')
                ENDIF
                i = LEN_TRIM(ts_filename_wtp_pln)
                IF(j .LT. 10  ) WRITE(ts_filename_wtp_pln(i-11:i-11),'(I1)') j
                IF(j .GE. 10  ) WRITE(ts_filename_wtp_pln(i-12:i-11),'(I2)') j
                IF(j .GE. 100 ) WRITE(ts_filename_wtp_pln(i-13:i-12),'(I3)') j
                IF(j .GE. 1000) WRITE(ts_filename_wtp_pln(i-14:i-13),'(I4)') j
                DO ij = 1, config_flags%wind_wtp_pln_nxy
                  WRITE(ts_filename_wtp_pln,'(A)') 'turbine.0000.plane.0000'
                  i = LEN_TRIM(ts_filename_wtp_pln)
                  IF(j .LT. 10  ) WRITE(ts_filename_wtp_pln(i-11:i-11),'(I1)') j
                  IF(j .GE. 10  ) WRITE(ts_filename_wtp_pln(i-12:i-11),'(I2)') j
                  IF(j .GE. 100 ) WRITE(ts_filename_wtp_pln(i-13:i-12),'(I3)') j
                  IF(j .GE. 1000) WRITE(ts_filename_wtp_pln(i-14:i-13),'(I4)') j
                  IF(ij .LT. 10  ) WRITE(ts_filename_wtp_pln(i  :i),'(I1)') ij
                  IF(ij .GE. 10  ) WRITE(ts_filename_wtp_pln(i-1:i),'(I2)') ij
                  IF(ij .GE. 100 ) WRITE(ts_filename_wtp_pln(i-2:i),'(I3)') ij
                  IF(ij .GE. 1000) WRITE(ts_filename_wtp_pln(i-3:i),'(I4)') ij
                  WRITE(ts_filename_wtp_pln,'(A)') 'tsout/'//TRIM(ts_filename_wtp_pln)
                  OPEN(UNIT=iunit, FILE=TRIM(ts_filename_wtp_pln), FORM='FORMATTED', STATUS='REPLACE')
#if (EM_CORE == 1)
                  WRITE(UNIT=iunit, FMT='(A11,A16,A12,A11,A13,A13,A13,A15,A15,A15)') &
                  'timestep[-]', 'time[s]', 'dt[s]', 'plEl[-]', 'u[m/s]', 'v[m/s]',  &
                  'w[m/s]', 'plPts_x[m]', 'plPts_y[m]', 'plPts_z[m]'
#endif
                  CLOSE(UNIT=iunit)
                ENDDO
              ENDDO
            ENDIF
          ELSE
            CALL wrf_error_fatal('Error in calc_ts_wtp_pln_initialize: No tsout directory is found. &
                                  Create "./tsout/" directory to save time series.')
          ENDIF
        ENDIF
      ENDIF
    ENDIF
END SUBROUTINE calc_ts_wtp_pln_initialize
!
SUBROUTINE calc_ts_wtp_pln_advance( grid )
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
    REAL, ALLOCATABLE, DIMENSION(:, :, :) :: u_temp, v_temp, w_temp
    REAL, ALLOCATABLE, DIMENSION(:, :, :) :: planePoints_x_temp
    REAL, ALLOCATABLE, DIMENSION(:, :, :) :: planePoints_y_temp
    REAL, ALLOCATABLE, DIMENSION(:, :, :) :: planePoints_z_temp
    TYPE(grid_config_rec_type) :: config_flags
    CALL model_to_grid_config_rec(grid%id, model_config_rec, config_flags)
#if ((EM_CORE == 1) && (DA_CORE != 1))
    IF( grid%dfi_opt /= DFI_NODFI .AND. grid%dfi_stage /= DFI_FST ) RETURN
#endif
    n = grid%next_ts_time_wtp_pln
    IF( config_flags%wind_wtp_opt .GT. 0 ) THEN
      IF( config_flags%wind_wtp_pln_tsout_opt ) THEN
#if ( EM_CORE == 1 )
        ALLOCATE( u_temp(config_flags%wind_wtp_pln_nTurbines, &
                         config_flags%wind_wtp_pln_nz, &
                         config_flags%wind_wtp_pln_nxy) )
        ALLOCATE( v_temp(config_flags%wind_wtp_pln_nTurbines, &
                         config_flags%wind_wtp_pln_nz, &
                         config_flags%wind_wtp_pln_nxy) )
        ALLOCATE( w_temp(config_flags%wind_wtp_pln_nTurbines, &
                         config_flags%wind_wtp_pln_nz, &
                         config_flags%wind_wtp_pln_nxy) )
        ALLOCATE( planePoints_x_temp(config_flags%wind_wtp_pln_nTurbines, &
                         config_flags%wind_wtp_pln_nz, &
                                     config_flags%wind_wtp_pln_nxy) )
        ALLOCATE( planePoints_y_temp(config_flags%wind_wtp_pln_nTurbines, &
                         config_flags%wind_wtp_pln_nz, &
                                     config_flags%wind_wtp_pln_nxy) )
        ALLOCATE( planePoints_z_temp(config_flags%wind_wtp_pln_nTurbines, &
                                     config_flags%wind_wtp_pln_nz, &
                                     config_flags%wind_wtp_pln_nxy) )
        u_temp = RESHAPE( grid%wtp_u_plane, (/ config_flags%wind_wtp_pln_nTurbines, &
                                               config_flags%wind_wtp_pln_nz, &
                                               config_flags%wind_wtp_pln_nxy /) )
        v_temp = RESHAPE( grid%wtp_v_plane, (/ config_flags%wind_wtp_pln_nTurbines, &
                                               config_flags%wind_wtp_pln_nz, &
                                               config_flags%wind_wtp_pln_nxy /) )
        w_temp = RESHAPE( grid%wtp_w_plane, (/ config_flags%wind_wtp_pln_nTurbines, &
                                               config_flags%wind_wtp_pln_nz, &
                                               config_flags%wind_wtp_pln_nxy /) )
        planePoints_x_temp = RESHAPE( grid%wtp_planePoints_x, (/ config_flags%wind_wtp_pln_nTurbines, &
                                                                 config_flags%wind_wtp_pln_nz, &
                                                                 config_flags%wind_wtp_pln_nxy /) )
        planePoints_y_temp = RESHAPE( grid%wtp_planePoints_y, (/ config_flags%wind_wtp_pln_nTurbines, &
                                                                 config_flags%wind_wtp_pln_nz, &
                                                                 config_flags%wind_wtp_pln_nxy /) )
        planePoints_z_temp = RESHAPE( grid%wtp_planePoints_z, (/ config_flags%wind_wtp_pln_nTurbines, &
                                                                 config_flags%wind_wtp_pln_nz, &
                                                                 config_flags%wind_wtp_pln_nxy /) )
#endif
        CALL domain_clock_get( grid, minutesSinceSimulationStart=xtime_minutes )
        grid%wtp_ts_sec(n) = xtime_minutes*60.
        grid%wtp_ts_dt(n) = grid%dt
        grid%wtp_ts_timestep(n) = grid%itimestep
        DO k = 1, config_flags%wind_wtp_pln_nTurbines
          DO l = 1, config_flags%wind_wtp_pln_nxy
            DO m = 1, config_flags%wind_wtp_pln_nz
#if (EM_CORE == 1)
            grid%wtp_ts_u_plane(n,k,m,l) = u_temp(k,m,l)
            grid%wtp_ts_v_plane(n,k,m,l) = v_temp(k,m,l)
            grid%wtp_ts_w_plane(n,k,m,l) = w_temp(k,m,l)
            grid%wtp_ts_planePoints_x(n,k,m,l) = planePoints_x_temp(k,m,l)
            grid%wtp_ts_planePoints_y(n,k,m,l) = planePoints_y_temp(k,m,l)
            grid%wtp_ts_planePoints_z(n,k,m,l) = planePoints_z_temp(k,m,l)
#endif
            ENDDO
          ENDDO
        ENDDO
#if (EM_CORE == 1) 
        DEALLOCATE(u_temp, v_temp, w_temp)
        DEALLOCATE(planePoints_x_temp)
        DEALLOCATE(planePoints_y_temp)
        DEALLOCATE(planePoints_z_temp)
#endif
        grid%next_ts_time_wtp_pln = grid%next_ts_time_wtp_pln + 1
        IF( grid%next_ts_time_wtp_pln .GT. config_flags%wind_wtp_pln_ts_buf_size ) CALL write_ts_wtp_pln(grid)
      ENDIF
    ENDIF
END SUBROUTINE calc_ts_wtp_pln_advance
!!
SUBROUTINE write_ts_wtp_pln( grid )
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
    CHARACTER (LEN=35) :: ts_filename_wtp_pln
    CALL model_to_grid_config_rec(grid%id, model_config_rec, config_flags)
    IF( config_flags%wind_wtp_opt .GT. 0 ) THEN
      IF( config_flags%wind_wtp_pln_tsout_opt ) THEN
#if ((EM_CORE == 1) && (DA_CORE != 1))
        IF( grid%dfi_opt /= DFI_NODFI .AND. grid%dfi_stage /= DFI_FST ) RETURN
#endif
        IF( wrf_dm_on_monitor() ) THEN
          DO k = 1, config_flags%wind_wtp_pln_nTurbines
            j = model_config_rec%wind_wtp_pln_turbine_ids(k)
            WRITE(ts_filename_wtp_pln,'(A)') 'turbine.0000.plane.0000'
            ! Create the output files for the vertical profiles, one file for each variable
            iunit = get_unused_unit()
            IF( iunit .LE. 0 ) THEN
              CALL wrf_error_fatal('Error in write_ts_wtp_pln: could not find a free Fortran unit.')
            ENDIF
            i = LEN_TRIM(ts_filename_wtp_pln)
            IF(j .LT. 10  ) WRITE(ts_filename_wtp_pln(i-11:i-11),'(I1)') j
            IF(j .GE. 10  ) WRITE(ts_filename_wtp_pln(i-12:i-11),'(I2)') j
            IF(j .GE. 100 ) WRITE(ts_filename_wtp_pln(i-13:i-12),'(I3)') j
            IF(j .GE. 1000) WRITE(ts_filename_wtp_pln(i-14:i-13),'(I4)') j
            DO l = 1, config_flags%wind_wtp_pln_nxy
              WRITE(ts_filename_wtp_pln,'(A)') 'turbine.0000.plane.0000'
              i = LEN_TRIM(ts_filename_wtp_pln)
              IF(j .LT. 10  ) WRITE(ts_filename_wtp_pln(i-11:i-11),'(I1)') j
              IF(j .GE. 10  ) WRITE(ts_filename_wtp_pln(i-12:i-11),'(I2)') j
              IF(j .GE. 100 ) WRITE(ts_filename_wtp_pln(i-13:i-12),'(I3)') j
              IF(j .GE. 1000) WRITE(ts_filename_wtp_pln(i-14:i-13),'(I4)') j
              IF(l .LT. 10  ) WRITE(ts_filename_wtp_pln(i  :i),'(I1)') l
              IF(l .GE. 10  ) WRITE(ts_filename_wtp_pln(i-1:i),'(I2)') l
              IF(l .GE. 100 ) WRITE(ts_filename_wtp_pln(i-2:i),'(I3)') l
              IF(l .GE. 1000) WRITE(ts_filename_wtp_pln(i-3:i),'(I3)') l
              WRITE(ts_filename_wtp_pln,'(A)') 'tsout/'//TRIM(ts_filename_wtp_pln)
              OPEN(UNIT=iunit, FILE=ts_filename_wtp_pln, STATUS='unknown', POSITION='append', FORM='formatted')
              DO n = 1, grid%next_ts_time_wtp_pln - 1
#if (EM_CORE == 1)
                WRITE(UNIT=iunit,FMT='(I11,3x,f13.6,3x,f9.6,3x)')  &
                                  grid%wtp_ts_timestep(n), grid%wtp_ts_sec(n), grid%wtp_ts_dt(n)
#endif
#if (EM_CORE == 1)
                DO m = 1, config_flags%wind_wtp_pln_nz
                  WRITE(UNIT=iunit,FMT='(I50,3x,f10.6,3x,f10.6,3x,f10.6,3x,f12.6,3x,f12.6,3x,f12.6,3x)') &
                  m, grid%wtp_ts_u_plane(n,k,m,l), grid%wtp_ts_v_plane(n,k,m,l), grid%wtp_ts_w_plane(n,k,m,l), &
                  grid%wtp_ts_planePoints_x(n,k,m,l), grid%wtp_ts_planePoints_y(n,k,m,l), grid%wtp_ts_planePoints_z(n,k,m,l)
#endif
                ENDDO
              ENDDO
              CLOSE(UNIT=iunit)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
    grid%next_ts_time_wtp_pln = 1
    ENDIF
END SUBROUTINE write_ts_wtp_pln
!!
SUBROUTINE wind_wtp_pln_tsout_dir(exists)
   IMPLICIT NONE
   LOGICAL , INTENT(OUT) :: exists
   INQUIRE(FILE='./tsout/', EXIST=exists)
END SUBROUTINE wind_wtp_pln_tsout_dir


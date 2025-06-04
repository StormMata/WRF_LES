!WRF:MEDIATION_LAYER:PHYSICS
!
MODULE module_wind_driver
  CONTAINS
  SUBROUTINE wind_driver(                                         &
        &  grid, wind_wtp_opt, nTrb, nSct, nElm, nVertPts         &    
        &, nHorizPts, nNacPts, nTowPts, nAveSmpls, ru_tendf       &
        &, rv_tendf, rw_tendf, tke_tend, muu, muv, mut, ht        &
        &, u_2, v_2, w_2, dt, dx, dy, rdzw                        &
        &, ids,ide, jds,jde, kds,kde                              &
        &, ims,ime, jms,jme, kms,kme                              &
        &, its,ite, jts,jte, kts,kte                              &
        &, itimestep, restart, iFst, jFst, windDirFstAve          &
        &, rhoFstTS, uFstTS, vFstTS, wFstTS, windDirFstTS         &
        &, wsRotor, turbDir, capPhi, eYaw, rho, rhoFstAve         &
        &, tsr, thr, trq, pow, mechPow, genPow, thrCoeff          &
        &, trqCoeff, powCoeff, fn, ft, F, phi, alpha, cl, cd, cm  &
        &, L, D, v1, vTheta, vT, vRel, vZ, v0Fst, uLES, vLES      &
        &, wLES, aif, rotSpeed, psi, azimuth, deltaAzimuth        &
        &, trqGenCommandedOld, genSpeedOld, genSpeedFOld          &
        &, pitch, speedError, intSpeedError, rotorApexX           &
        &, rotorApexY, rotorApexZ, bladePointsX, bladePointsY     &
        &, bladePointsZ, towerPointsX, towerPointsY, towerPointsZ &
        &, nacellePointsX, nacellePointsY, nacellePointsZ         &
        &, readFromRSTFile, nPlnTrb, uPlane, vPlane, wPlane       &
        &, planePointsX, planePointsY, planePointsZ, fnAve, ftAve &
        &, FAve, phiAve, alphaAve, clAve, cdAve, cmAve,LAve, DAve &
        &, v1Ave, vThetaAve, vTAve, vRelAve, vZAve, uLESAve       &
        &, vLESAve, wLESAve, shapiroM                             )
! Driver layer modules
  USE module_domain
  USE module_configure
  USE module_tiles
  USE module_dm
  USE module_machine
  USE module_state_description
! Model layer modules
  USE module_bc
! *** add new modules of schemes here
#if (EM_CORE==1)
  USE module_wind_gal
  USE module_wind_gad
  USE module_wind_gad_rs
  USE module_wind_sad
#endif
!=======================================================================
!  This driver calls GAL/GAD/GADrs/SAD model subroutines.
!  GAL scheme   : Generalized Actuator Line
!  GAD scheme   : Generalized Actuator Disk
!  GADrs scheme : Generalized Actuator Disk ring-segment
!  SAD scheme   : Simplified Actuator Disk
!=======================================================================
  IMPLICIT NONE
    TYPE(domain), TARGET         ::                              grid
    INTEGER,    INTENT(IN)       ::                              nTrb, &
                                                                 nSct, &
                                                                 nElm, &
                                                             nVertPts, &
                                                            nHorizPts, &
                                                              nNacPts, &
                                                              nTowPts, &
                                                            nAveSmpls, &
                                                              nPlnTrb
    INTEGER,    INTENT(IN)    ::              wind_wtp_opt, itimestep
    INTEGER,    INTENT(IN)    ::            ids,ide, jds,jde, kds,kde, &
                                            ims,ime, jms,jme, kms,kme, &
                                            its,ite, jts,jte, kts,kte 
    REAL,       INTENT(IN)    ::                           DX, DY, DT
    REAL,       DIMENSION( ims:ime, kms:kme, jms:jme ),                &
                INTENT(IN)              ::                       rdzw, &
                                                   rho, u_2, v_2, w_2
    REAL,       DIMENSION( ims:ime, jms:jme ),                         &
                INTENT(IN)    ::                                  muu, &
                                                                  muv, &
                                                                  mut, &
                                                                  ht
    REAL,       DIMENSION( ims:ime, kms:kme, jms:jme ),                &
                INTENT(INOUT)    ::                          ru_tendf, &
                                                             rv_tendf, &
                                                             rw_tendf, &
                                                            tke_tend
    LOGICAL,    INTENT(IN)    ::                              restart
    INTEGER,    INTENT(INOUT)    ::                   readFromRSTFile
    INTEGER,    DIMENSION(1:nTrb),                                     &
                INTENT(INOUT) ::                           iFst, jFst
    REAL,       DIMENSION(1:nTrb),                                     &
                INTENT(INOUT) ::                              turbDir, &
                                                              wsRotor, &
                                                                 eYaw, &
                                                                  tsr, &
                                                         deltaAzimuth, &
                                                             rotSpeed, &
                                                   trqGenCommandedOld, &
                                                          genSpeedOld, &
                                                         genSpeedFOld, &
                                                                pitch, &
                                                           speedError, &
                                                        intSpeedError, &
                                                                  psi, &
                                                              azimuth, &
                                                           rotorApexX, &
                                                           rotorApexY, &
                                                           rotorApexZ
    REAL,       DIMENSION(1:nTrb),                                     &
                INTENT(OUT) ::                          windDirFstAve, &
                                                               capPhi, &
                                                            rhoFstAve, &
                                                                v0Fst, &
                                                                  aif, &
                                                                  thr, &
                                                                  trq, &
                                                                  pow, &
                                                              mechPow, &
                                                             thrCoeff, &
                                                             trqCoeff, &
                                                             powCoeff, &
                                                             shapiroM, &
                                                               genPow
    REAL,       DIMENSION(1:nTrb, 1:nElm),                             &
                INTENT(OUT) ::                                  fnAve, &
                                                                ftAve, &
                                                                 FAve, &
                                                               phiAve, &
                                                             alphaAve, &
                                                                clAve, &
                                                                cdAve, &
                                                                cmAve, &
                                                                 LAve, &
                                                                 DAve, &
                                                                v1Ave, &
                                                            vThetaAve, &
                                                                vTAve, &
                                                              vRelAve, &
                                                                vZAve, &
                                                              uLESAve, &
                                                              vLESAve, &
                                                              wLESAve
                                                         
    REAL,       DIMENSION(1:nTrb, 1:nTowPts),                          &
                INTENT(OUT) ::                           towerPointsX, &
                                                         towerPointsY, &
                                                         towerPointsZ
    REAL,       DIMENSION(1:nTrb, 1:nNacPts),                          &
                INTENT(OUT) ::                         nacellePointsX, &
                                                       nacellePointsY, &
                                                       nacellePointsZ
    REAL,       DIMENSION(1:nTrb, 1:nAveSmpls),                        &
                INTENT(INOUT) ::                             rhoFstTS, &
                                                               uFstTS, &
                                                               vFstTS, &
                                                               wFstTS, &
                                                         windDirFstTS
    REAL,       DIMENSION(1:nTrb, 1:nSct, 1:nElm),                     &
                INTENT(OUT) ::                                     fn, &
                                                                   ft, &
                                                                    F, &
                                                                  phi, &
                                                                alpha, &
                                                                   cl, &
                                                                   cd, &
                                                                   cm, &
                                                                    L, &
                                                                    D, &
                                                                   v1, &
                                                               vTheta, &
                                                                   vT, &
                                                                 vRel, &
                                                                   vZ, &
                                                                 uLES, &
                                                                 vLES, &
                                                                 wLES, &
                                                         bladePointsX, &
                                                         bladePointsY, &
                                                         bladePointsZ
    REAL,       DIMENSION(1:nPlnTrb, 1:nVertPts, 1:nHorizPts),         &
                INTENT(OUT) ::                                 uPlane, &
                                                               vPlane, &
                                                               wPlane
    REAL,       DIMENSION(1: nPlnTrb, 1:nVertPts, 1:nHorizPts),        &
                INTENT(INOUT) ::                         planePointsX, &
                                                         planePointsY, &
                                                         planePointsZ
! local
    LOGICAL       :: wtp_history_flag, wtp_restart_flag
!=======================================================================
    IF( wind_wtp_opt .GT. 0 ) THEN
      wtp_history_flag = .FALSE.
      IF( Is_alarm_tstep( grid%domain_clock, grid%alarms(HISTORY_ALARM) ) ) THEN
        wtp_history_flag = .TRUE.
      ENDIF
      wtp_restart_flag = .FALSE.
      IF( Is_alarm_tstep( grid%domain_clock, grid%alarms(RESTART_ALARM) ) ) THEN
        wtp_restart_flag = .TRUE.
      ENDIF
#if (EM_CORE==1)
      wind_wtpar_select: SELECT CASE(wind_wtp_opt)
! ... BK: Generalized Actuator Line parameterization ...
      CASE(galscheme)
        CALL wrf_debug(100,'in phys/module_wind_gal.F')
        CALL calc_forces_gal(                                &
            & dx, dy, rdzw, dt, muu, muv, mut                &
            &,ht, ru_tendf, rv_tendf, rw_tendf               &
            &,u_2, v_2, w_2                                  &
            &,ids,ide, jds,jde, kds,kde                      &
            &,ims,ime, jms,jme, kms,kme                      &
            &,its,ite, jts,jte, kts,kte                      &
            &,itimestep, restart, iFst, jFst                 &
            &,windDirFstAve, rhoFstAve, v0Fst                &
            &,rhoFstTS, uFstTS, vFstTS, wFstTS, windDirFstTS &
            &,wsRotor, turbDir, deltaAzimuth, capPhi         &
            &,eYaw, rho, tsr, thr, trq, pow, mechPow, genPow &
            &,thrCoeff, trqCoeff, powCoeff, fn, ft, F, phi   &
            &,alpha, cl, cd, cm, L, D, v1, vTheta            &
            &,vT, vRel, vZ, uLES, vLES, wLES                 &
            &,rotSpeed, psi, azimuth, trqGenCommandedOld     &
            &,genSpeedOld, genSpeedFOld                      &
            &,pitch, speedError, intSpeedError, rotorApexX   &
            &,rotorApexY, rotorApexZ, bladePointsX           &
            &,bladePointsY, bladePointsZ, towerPointsX       &
            &,towerPointsY, towerPointsZ, nacellePointsX     &
            &,nacellePointsY, nacellePointsZ                 &
            &,wtp_history_flag, wtp_restart_flag             &
            &,readFromRSTFile, uPlane, vPlane, wPlane        &
            &,planePointsX, planePointsY, planePointsZ       )
! ... BK: Generalized Actuator Disk parameterization ...
      CASE(gadscheme)
        CALL wrf_debug(100,'in phys/module_wind_gad.F')       
        CALL calc_forces_gad(                                &
            & dx, dy, rdzw, dt, muu, muv, mut                &
            &,ht, ru_tendf, rv_tendf, rw_tendf               &
            &,u_2, v_2, w_2                                  &
            &,ids,ide, jds,jde, kds,kde                      &
            &,ims,ime, jms,jme, kms,kme                      &
            &,its,ite, jts,jte, kts,kte                      &
            &,itimestep, restart, iFst, jFst                 &
            &,windDirFstAve, rhoFstAve, v0Fst                &
            &,rhoFstTS, uFstTS, vFstTS, wFstTS               &
            &,windDirFstTS, wsRotor, turbDir, capPhi         &
            &,eYaw, rho, tsr, thr, trq, pow, mechPow, genPow &
            &,thrCoeff, trqCoeff, powCoeff, fn, ft, F, phi   &
            &,alpha, cl, cd, cm, L, D, v1, vTheta            &
            &,vT, vRel, vZ, uLES, vLES, wLES                 &
            &,rotSpeed, trqGenCommandedOld                   &
            &,genSpeedOld, genSpeedFOld                      &
            &,pitch, speedError, intSpeedError, rotorApexX   &
            &,rotorApexY, rotorApexZ, bladePointsX           &
            &,bladePointsY, bladePointsZ, towerPointsX       &
            &,towerPointsY, towerPointsZ, nacellePointsX     &
            &,nacellePointsY, nacellePointsZ                 &
            &,wtp_history_flag, wtp_restart_flag             &
            &,readFromRSTFile, uPlane, vPlane, wPlane        &
            &,planePointsX, planePointsY, planePointsZ       &
            &,shapiroM                                       )
! ... BK: Generalized Actuator Disk ring-segment parameterization ...
      CASE(gadrsscheme)
        CALL wrf_debug(100,'in phys/module_wind_gad_rs.F')       
        CALL calc_forces_gad_rs(                             &
            & dx, dy, rdzw, dt, muu, muv, mut                &
            &,ht, ru_tendf, rv_tendf, rw_tendf               &
            &,u_2, v_2, w_2                                  &
            &,ids,ide, jds,jde, kds,kde                      &
            &,ims,ime, jms,jme, kms,kme                      &
            &,its,ite, jts,jte, kts,kte                      &
            &,itimestep, restart, iFst, jFst, windDirFstAve  &
            &,rhoFstAve, v0Fst, rhoFstTS, uFstTS, vFstTS     &
            &,wFstTS, windDirFstTS, wsRotor, turbDir, capPhi &
            &,eYaw, rho, tsr, thr, trq, pow, mechPow, genPow &
            &,thrCoeff, trqCoeff, powCoeff, fnAve, ftAve     &
            &,FAve, phiAve, alphaAve, clAve, cdAve, cmAve    &
            &,LAve, DAve, v1Ave, vThetaAve, vTAve, vRelAve   &
            &,vZAve, uLESAve, vLESAve, wLESAve, rotSpeed     &
            &,trqGenCommandedOld, genSpeedOld, genSpeedFOld  &
            &,pitch, speedError, intSpeedError, rotorApexX   &
            &,rotorApexY, rotorApexZ, bladePointsX           &
            &,bladePointsY, bladePointsZ, towerPointsX       &
            &,towerPointsY, towerPointsZ, nacellePointsX     &
            &,nacellePointsY, nacellePointsZ                 &
            &,wtp_history_flag, wtp_restart_flag             &
            &,readFromRSTFile, uPlane, vPlane, wPlane        &
            &,planePointsX, planePointsY, planePointsZ       )
! ... BK: Simplified Actuator Disk parameterization ...
      CASE(sadscheme)
        CALL wrf_debug(100,'in phys/module_wind_sad.F')       
        CALL calc_power_ct_sad(                              &
            & dx, dy, rdzw, dt, muu, muv, mut                &
            &,ht, ru_tendf, rv_tendf, tke_tend               &
            &,u_2, v_2, w_2                                  &
            &,ids,ide, jds,jde, kds,kde                      &
            &,ims,ime, jms,jme, kms,kme                      &
            &,its,ite, jts,jte, kts,kte                      &
            &,iFst, jFst, windDirFstAve, rhoFstAve, v0Fst    &
            &,rhoFstTS, uFstTS, vFstTS, windDirFstTS         &
            &,turbDir, eYaw, rho, pow, powCoeff, thrCoeff    &
            &,aif, wtp_history_flag, wtp_restart_flag        &
            &,readFromRSTFile, uPlane, vPlane, wPlane        &
            &,planePointsX, planePointsY, planePointsZ       )
      END SELECT wind_wtpar_select
#endif
    ENDIF
  END SUBROUTINE wind_driver
END MODULE module_wind_driver

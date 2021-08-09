












! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_config
!
!> \brief MPAS ocean specific config
!> \details
!>  This module contains config specific to the ocean model.
!
!-----------------------------------------------------------------------

module ocn_config

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_kind_types

   implicit none
   public
   save

      character (len=StrKIND), pointer :: config_ocean_run_mode

      logical, pointer :: config_do_restart
      character (len=StrKIND), pointer :: config_restart_timestamp_name
      character (len=StrKIND), pointer :: config_start_time
      character (len=StrKIND), pointer :: config_stop_time
      character (len=StrKIND), pointer :: config_run_duration
      character (len=StrKIND), pointer :: config_calendar_type

      logical, pointer :: config_write_output_on_startup
      integer, pointer :: config_pio_num_iotasks
      integer, pointer :: config_pio_stride

      integer, pointer :: config_num_halos
      character (len=StrKIND), pointer :: config_block_decomp_file_prefix
      integer, pointer :: config_number_of_blocks
      logical, pointer :: config_explicit_proc_decomp
      character (len=StrKIND), pointer :: config_proc_decomp_file_prefix

      character (len=StrKIND), pointer :: config_dt
      character (len=StrKIND), pointer :: config_time_integrator

      logical, pointer :: config_hmix_scaleWithMesh
      real (kind=RKIND), pointer :: config_maxMeshDensity
      logical, pointer :: config_hmix_use_ref_cell_width
      real (kind=RKIND), pointer :: config_hmix_ref_cell_width
      real (kind=RKIND), pointer :: config_apvm_scale_factor

      logical, pointer :: config_use_mom_del2
      real (kind=RKIND), pointer :: config_mom_del2
      logical, pointer :: config_use_tracer_del2
      real (kind=RKIND), pointer :: config_tracer_del2

      logical, pointer :: config_use_mom_del4
      real (kind=RKIND), pointer :: config_mom_del4
      real (kind=RKIND), pointer :: config_mom_del4_div_factor
      logical, pointer :: config_use_tracer_del4
      real (kind=RKIND), pointer :: config_tracer_del4

      logical, pointer :: config_use_Leith_del2
      real (kind=RKIND), pointer :: config_Leith_parameter
      real (kind=RKIND), pointer :: config_Leith_dx
      real (kind=RKIND), pointer :: config_Leith_visc2_max

      character (len=StrKIND), pointer :: config_eddying_resolution_taper
      real (kind=RKIND), pointer :: config_eddying_resolution_ramp_min
      real (kind=RKIND), pointer :: config_eddying_resolution_ramp_max

      logical, pointer :: config_use_Redi
      character (len=StrKIND), pointer :: config_Redi_closure
      real (kind=RKIND), pointer :: config_Redi_constant_kappa
      real (kind=RKIND), pointer :: config_Redi_maximum_slope
      logical, pointer :: config_Redi_use_slope_taper
      logical, pointer :: config_Redi_use_surface_taper
      logical, pointer :: config_Redi_N2_based_taper_enable
      real (kind=RKIND), pointer :: config_Redi_N2_based_taper_min
      logical, pointer :: config_Redi_N2_based_taper_limit_term1

      logical, pointer :: config_use_GM
      character (len=StrKIND), pointer :: config_GM_closure
      real (kind=RKIND), pointer :: config_GM_constant_kappa
      real (kind=RKIND), pointer :: config_GM_constant_gravWaveSpeed
      real (kind=RKIND), pointer :: config_GM_spatially_variable_min_kappa
      real (kind=RKIND), pointer :: config_GM_spatially_variable_max_kappa
      real (kind=RKIND), pointer :: config_GM_spatially_variable_baroclinic_mode
      real (kind=RKIND), pointer :: config_GM_Visbeck_alpha
      real (kind=RKIND), pointer :: config_GM_Visbeck_max_depth
      real (kind=RKIND), pointer :: config_GM_EG_riMin
      real (kind=RKIND), pointer :: config_GM_EG_kappa_factor
      real (kind=RKIND), pointer :: config_GM_EG_Rossby_factor
      real (kind=RKIND), pointer :: config_GM_EG_Rhines_factor

      logical, pointer :: config_Rayleigh_friction
      real (kind=RKIND), pointer :: config_Rayleigh_damping_coeff
      logical, pointer :: config_Rayleigh_damping_depth_variable
      logical, pointer :: config_Rayleigh_bottom_friction
      real (kind=RKIND), pointer :: config_Rayleigh_bottom_damping_coeff

      logical, pointer :: config_use_cvmix
      real (kind=RKIND), pointer :: config_cvmix_prandtl_number
      character (len=StrKIND), pointer :: config_cvmix_background_scheme
      real (kind=RKIND), pointer :: config_cvmix_background_diffusion
      real (kind=RKIND), pointer :: config_cvmix_background_viscosity
      real (kind=RKIND), pointer :: config_cvmix_BryanLewis_bl1
      real (kind=RKIND), pointer :: config_cvmix_BryanLewis_bl2
      real (kind=RKIND), pointer :: config_cvmix_BryanLewis_transitionDepth
      real (kind=RKIND), pointer :: config_cvmix_BryanLewis_transitionWidth
      logical, pointer :: config_use_cvmix_convection
      real (kind=RKIND), pointer :: config_cvmix_convective_diffusion
      real (kind=RKIND), pointer :: config_cvmix_convective_viscosity
      logical, pointer :: config_cvmix_convective_basedOnBVF
      real (kind=RKIND), pointer :: config_cvmix_convective_triggerBVF
      logical, pointer :: config_use_cvmix_shear
      integer, pointer :: config_cvmix_num_ri_smooth_loops
      logical, pointer :: config_cvmix_use_BLD_smoothing
      character (len=StrKIND), pointer :: config_cvmix_shear_mixing_scheme
      real (kind=RKIND), pointer :: config_cvmix_shear_PP_nu_zero
      real (kind=RKIND), pointer :: config_cvmix_shear_PP_alpha
      real (kind=RKIND), pointer :: config_cvmix_shear_PP_exp
      real (kind=RKIND), pointer :: config_cvmix_shear_KPP_nu_zero
      real (kind=RKIND), pointer :: config_cvmix_shear_KPP_Ri_zero
      real (kind=RKIND), pointer :: config_cvmix_shear_KPP_exp
      logical, pointer :: config_use_cvmix_tidal_mixing
      logical, pointer :: config_use_cvmix_double_diffusion
      logical, pointer :: config_use_cvmix_kpp
      logical, pointer :: config_use_cvmix_fixed_boundary_layer
      real (kind=RKIND), pointer :: config_cvmix_kpp_boundary_layer_depth
      real (kind=RKIND), pointer :: config_cvmix_kpp_criticalBulkRichardsonNumber
      character (len=StrKIND), pointer :: config_cvmix_kpp_matching
      logical, pointer :: config_cvmix_kpp_EkmanOBL
      logical, pointer :: config_cvmix_kpp_MonObOBL
      character (len=StrKIND), pointer :: config_cvmix_kpp_interpolationOMLType
      real (kind=RKIND), pointer :: config_cvmix_kpp_surface_layer_extent
      real (kind=RKIND), pointer :: config_cvmix_kpp_surface_layer_averaging
      real (kind=RKIND), pointer :: configure_cvmix_kpp_minimum_OBL_under_sea_ice
      real (kind=RKIND), pointer :: config_cvmix_kpp_stop_OBL_search
      logical, pointer :: config_cvmix_kpp_use_enhanced_diff
      logical, pointer :: config_cvmix_kpp_nonlocal_with_implicit_mix
      logical, pointer :: config_cvmix_kpp_use_theory_wave
      character (len=StrKIND), pointer :: config_cvmix_kpp_langmuir_mixing_opt
      character (len=StrKIND), pointer :: config_cvmix_kpp_langmuir_entrainment_opt

      logical, pointer :: config_use_gotm
      character (len=StrKIND), pointer :: config_gotm_namelist_file
      real (kind=RKIND), pointer :: config_gotm_constant_surface_roughness_length
      real (kind=RKIND), pointer :: config_gotm_constant_bottom_roughness_length
      real (kind=RKIND), pointer :: config_gotm_constant_bottom_drag_coeff

      logical, pointer :: config_use_variable_drag
      logical, pointer :: config_use_bulk_wind_stress
      logical, pointer :: config_use_bulk_thickness_flux
      real (kind=RKIND), pointer :: config_flux_attenuation_coefficient
      real (kind=RKIND), pointer :: config_flux_attenuation_coefficient_runoff

      logical, pointer :: config_use_time_varying_atmospheric_forcing
      character (len=StrKIND), pointer :: config_time_varying_atmospheric_forcing_type
      character (len=StrKIND), pointer :: config_time_varying_atmospheric_forcing_start_time
      character (len=StrKIND), pointer :: config_time_varying_atmospheric_forcing_reference_time
      character (len=StrKIND), pointer :: config_time_varying_atmospheric_forcing_cycle_start
      character (len=StrKIND), pointer :: config_time_varying_atmospheric_forcing_cycle_duration
      character (len=StrKIND), pointer :: config_time_varying_atmospheric_forcing_interval
      real (kind=RKIND), pointer :: config_time_varying_atmospheric_forcing_ramp
      real (kind=RKIND), pointer :: config_time_varying_atmospheric_forcing_ramp_delay
      logical, pointer :: config_use_time_varying_land_ice_forcing
      character (len=StrKIND), pointer :: config_time_varying_land_ice_forcing_start_time
      character (len=StrKIND), pointer :: config_time_varying_land_ice_forcing_reference_time
      character (len=StrKIND), pointer :: config_time_varying_land_ice_forcing_cycle_start
      character (len=StrKIND), pointer :: config_time_varying_land_ice_forcing_cycle_duration
      character (len=StrKIND), pointer :: config_time_varying_land_ice_forcing_interval

      real (kind=RKIND), pointer :: config_ssh_grad_relax_timescale
      logical, pointer :: config_remove_AIS_coupler_runoff

      character (len=StrKIND), pointer :: config_sw_absorption_type
      integer, pointer :: config_jerlov_water_type
      real (kind=RKIND), pointer :: config_surface_buoyancy_depth
      logical, pointer :: config_enable_shortwave_energy_fixer

      logical, pointer :: config_use_tidal_forcing
      real (kind=RKIND), pointer :: config_use_tidal_forcing_tau
      character (len=StrKIND), pointer :: config_tidal_forcing_type
      character (len=StrKIND), pointer :: config_tidal_forcing_model
      real (kind=RKIND), pointer :: config_tidal_forcing_monochromatic_amp
      real (kind=RKIND), pointer :: config_tidal_forcing_monochromatic_period
      real (kind=RKIND), pointer :: config_tidal_forcing_monochromatic_phaseLag
      real (kind=RKIND), pointer :: config_tidal_forcing_monochromatic_baseline

      logical, pointer :: config_use_tidal_potential_forcing
      character (len=StrKIND), pointer :: config_tidal_potential_reference_time
      logical, pointer :: config_use_tidal_potential_forcing_M2
      logical, pointer :: config_use_tidal_potential_forcing_S2
      logical, pointer :: config_use_tidal_potential_forcing_N2
      logical, pointer :: config_use_tidal_potential_forcing_K2
      logical, pointer :: config_use_tidal_potential_forcing_K1
      logical, pointer :: config_use_tidal_potential_forcing_O1
      logical, pointer :: config_use_tidal_potential_forcing_Q1
      logical, pointer :: config_use_tidal_potential_forcing_P1
      real (kind=RKIND), pointer :: config_tidal_potential_ramp
      real (kind=RKIND), pointer :: config_self_attraction_and_loading_beta

      logical, pointer :: config_use_vegetation_drag
      logical, pointer :: config_use_vegetation_manning_equation
      real (kind=RKIND), pointer :: config_vegetation_drag_coefficient

      logical, pointer :: config_use_frazil_ice_formation
      logical, pointer :: config_frazil_in_open_ocean
      logical, pointer :: config_frazil_under_land_ice
      real (kind=RKIND), pointer :: config_frazil_heat_of_fusion
      real (kind=RKIND), pointer :: config_frazil_ice_density
      real (kind=RKIND), pointer :: config_frazil_fractional_thickness_limit
      real (kind=RKIND), pointer :: config_specific_heat_sea_water
      real (kind=RKIND), pointer :: config_frazil_maximum_depth
      real (kind=RKIND), pointer :: config_frazil_sea_ice_reference_salinity
      real (kind=RKIND), pointer :: config_frazil_land_ice_reference_salinity
      real (kind=RKIND), pointer :: config_frazil_maximum_freezing_temperature
      logical, pointer :: config_frazil_use_surface_pressure

      character (len=StrKIND), pointer :: config_land_ice_flux_mode
      character (len=StrKIND), pointer :: config_land_ice_flux_formulation
      logical, pointer :: config_land_ice_flux_useHollandJenkinsAdvDiff
      real (kind=RKIND), pointer :: config_land_ice_flux_attenuation_coefficient
      real (kind=RKIND), pointer :: config_land_ice_flux_boundaryLayerThickness
      real (kind=RKIND), pointer :: config_land_ice_flux_boundaryLayerNeighborWeight
      real (kind=RKIND), pointer :: config_land_ice_flux_cp_ice
      real (kind=RKIND), pointer :: config_land_ice_flux_rho_ice
      real (kind=RKIND), pointer :: config_land_ice_flux_topDragCoeff
      real (kind=RKIND), pointer :: config_land_ice_flux_ISOMIP_gammaT
      real (kind=RKIND), pointer :: config_land_ice_flux_rms_tidal_velocity
      real (kind=RKIND), pointer :: config_land_ice_flux_jenkins_heat_transfer_coefficient
      real (kind=RKIND), pointer :: config_land_ice_flux_jenkins_salt_transfer_coefficient

      character (len=StrKIND), pointer :: config_vert_tracer_adv
      integer, pointer :: config_vert_tracer_adv_order
      integer, pointer :: config_horiz_tracer_adv_order
      real (kind=RKIND), pointer :: config_coef_3rd_order
      logical, pointer :: config_monotonic

      logical, pointer :: config_use_implicit_bottom_drag
      real (kind=RKIND), pointer :: config_implicit_bottom_drag_coeff
      logical, pointer :: config_use_implicit_bottom_roughness
      logical, pointer :: config_use_implicit_bottom_drag_variable
      logical, pointer :: config_use_implicit_bottom_drag_variable_mannings
      logical, pointer :: config_use_explicit_bottom_drag
      real (kind=RKIND), pointer :: config_explicit_bottom_drag_coeff
      logical, pointer :: config_use_topographic_wave_drag
      real (kind=RKIND), pointer :: config_topographic_wave_drag_coeff

      logical, pointer :: config_use_wetting_drying
      logical, pointer :: config_prevent_drying
      real (kind=RKIND), pointer :: config_drying_min_cell_height
      logical, pointer :: config_zero_drying_velocity
      logical, pointer :: config_verify_not_dry
      character (len=StrKIND), pointer :: config_thickness_flux_type
      real (kind=RKIND), pointer :: config_drying_safety_height

      real (kind=RKIND), pointer :: config_density0

      character (len=StrKIND), pointer :: config_pressure_gradient_type
      real (kind=RKIND), pointer :: config_common_level_weight
      real (kind=RKIND), pointer :: config_zonal_ssh_grad
      real (kind=RKIND), pointer :: config_meridional_ssh_grad

      character (len=StrKIND), pointer :: config_eos_type
      real (kind=RKIND), pointer :: config_open_ocean_freezing_temperature_coeff_0
      real (kind=RKIND), pointer :: config_open_ocean_freezing_temperature_coeff_S
      real (kind=RKIND), pointer :: config_open_ocean_freezing_temperature_coeff_p
      real (kind=RKIND), pointer :: config_open_ocean_freezing_temperature_coeff_pS
      real (kind=RKIND), pointer :: config_open_ocean_freezing_temperature_coeff_mushy_az1_liq
      real (kind=RKIND), pointer :: config_land_ice_cavity_freezing_temperature_coeff_0
      real (kind=RKIND), pointer :: config_land_ice_cavity_freezing_temperature_coeff_S
      real (kind=RKIND), pointer :: config_land_ice_cavity_freezing_temperature_coeff_p
      real (kind=RKIND), pointer :: config_land_ice_cavity_freezing_temperature_coeff_pS

      real (kind=RKIND), pointer :: config_eos_linear_alpha
      real (kind=RKIND), pointer :: config_eos_linear_beta
      real (kind=RKIND), pointer :: config_eos_linear_Tref
      real (kind=RKIND), pointer :: config_eos_linear_Sref
      real (kind=RKIND), pointer :: config_eos_linear_densityref

      real (kind=RKIND), pointer :: config_eos_wright_ref_pressure

      integer, pointer :: config_n_ts_iter
      integer, pointer :: config_n_bcl_iter_beg
      integer, pointer :: config_n_bcl_iter_mid
      integer, pointer :: config_n_bcl_iter_end

      character (len=StrKIND), pointer :: config_btr_dt
      integer, pointer :: config_n_btr_cor_iter
      logical, pointer :: config_vel_correction
      integer, pointer :: config_btr_subcycle_loop_factor
      real (kind=RKIND), pointer :: config_btr_gam1_velWt1
      real (kind=RKIND), pointer :: config_btr_gam2_SSHWt1
      real (kind=RKIND), pointer :: config_btr_gam3_velWt2
      logical, pointer :: config_btr_solve_SSH2

      character (len=StrKIND), pointer :: config_btr_si_preconditioner
      real (kind=RKIND), pointer :: config_btr_si_tolerance
      integer, pointer :: config_n_btr_si_outer_iter
      logical, pointer :: config_btr_si_partition_match_mode

      character (len=StrKIND), pointer :: config_vert_coord_movement
      character (len=StrKIND), pointer :: config_ALE_thickness_proportionality
      real (kind=RKIND), pointer :: config_vert_taper_weight_depth_1
      real (kind=RKIND), pointer :: config_vert_taper_weight_depth_2
      logical, pointer :: config_use_min_max_thickness
      real (kind=RKIND), pointer :: config_min_thickness
      real (kind=RKIND), pointer :: config_max_thickness_factor
      logical, pointer :: config_dzdk_positive

      logical, pointer :: config_use_freq_filtered_thickness
      real (kind=RKIND), pointer :: config_thickness_filter_timescale
      logical, pointer :: config_use_highFreqThick_restore
      real (kind=RKIND), pointer :: config_highFreqThick_restore_time
      logical, pointer :: config_use_highFreqThick_del2
      real (kind=RKIND), pointer :: config_highFreqThick_del2

      logical, pointer :: config_check_zlevel_consistency
      logical, pointer :: config_check_ssh_consistency
      logical, pointer :: config_filter_btr_mode
      logical, pointer :: config_prescribe_velocity
      logical, pointer :: config_prescribe_thickness
      logical, pointer :: config_include_KE_vertex
      logical, pointer :: config_check_tracer_monotonicity
      logical, pointer :: config_compute_active_tracer_budgets
      logical, pointer :: config_disable_thick_all_tend
      logical, pointer :: config_disable_thick_hadv
      logical, pointer :: config_disable_thick_vadv
      logical, pointer :: config_disable_thick_sflux
      logical, pointer :: config_disable_vel_all_tend
      logical, pointer :: config_disable_vel_coriolis
      logical, pointer :: config_disable_vel_pgrad
      logical, pointer :: config_disable_vel_hmix
      logical, pointer :: config_disable_vel_surface_stress
      logical, pointer :: config_disable_vel_topographic_wave_drag
      logical, pointer :: config_disable_vel_explicit_bottom_drag
      logical, pointer :: config_disable_vel_vmix
      logical, pointer :: config_disable_vel_vadv
      logical, pointer :: config_disable_tr_all_tend
      logical, pointer :: config_disable_tr_adv
      logical, pointer :: config_disable_tr_hmix
      logical, pointer :: config_disable_tr_vmix
      logical, pointer :: config_disable_tr_sflux
      logical, pointer :: config_disable_tr_nonlocalflux
      logical, pointer :: config_disable_redi_k33
      logical, pointer :: config_read_nearest_restart

      logical, pointer :: config_conduct_tests
      logical, pointer :: config_test_tensors
      character (len=StrKIND), pointer :: config_tensor_test_function

      integer, pointer :: config_vert_levels

      logical, pointer :: config_use_activeTracers
      logical, pointer :: config_use_activeTracers_surface_bulk_forcing
      logical, pointer :: config_use_activeTracers_surface_restoring
      logical, pointer :: config_use_activeTracers_interior_restoring
      logical, pointer :: config_use_activeTracers_exponential_decay
      logical, pointer :: config_use_activeTracers_idealAge_forcing
      logical, pointer :: config_use_activeTracers_ttd_forcing
      logical, pointer :: config_use_surface_salinity_monthly_restoring
      character (len=StrKIND), pointer :: config_surface_salinity_monthly_restoring_compute_interval
      real (kind=RKIND), pointer :: config_salinity_restoring_constant_piston_velocity
      real (kind=RKIND), pointer :: config_salinity_restoring_max_difference
      logical, pointer :: config_salinity_restoring_under_sea_ice

      logical, pointer :: config_use_debugTracers
      logical, pointer :: config_reset_debugTracers_near_surface
      integer, pointer :: config_reset_debugTracers_top_nLayers
      logical, pointer :: config_use_debugTracers_surface_bulk_forcing
      logical, pointer :: config_use_debugTracers_surface_restoring
      logical, pointer :: config_use_debugTracers_interior_restoring
      logical, pointer :: config_use_debugTracers_exponential_decay
      logical, pointer :: config_use_debugTracers_idealAge_forcing
      logical, pointer :: config_use_debugTracers_ttd_forcing

      logical, pointer :: config_use_ecosysTracers
      character (len=StrKIND), pointer :: config_ecosys_atm_co2_option
      character (len=StrKIND), pointer :: config_ecosys_atm_alt_co2_option
      logical, pointer :: config_ecosys_atm_alt_co2_use_eco
      real (kind=RKIND), pointer :: config_ecosys_atm_co2_constant_value
      logical, pointer :: config_use_ecosysTracers_surface_bulk_forcing
      logical, pointer :: config_use_ecosysTracers_surface_restoring
      logical, pointer :: config_use_ecosysTracers_interior_restoring
      logical, pointer :: config_use_ecosysTracers_exponential_decay
      logical, pointer :: config_use_ecosysTracers_idealAge_forcing
      logical, pointer :: config_use_ecosysTracers_ttd_forcing
      logical, pointer :: config_use_ecosysTracers_surface_value
      logical, pointer :: config_use_ecosysTracers_sea_ice_coupling
      logical, pointer :: config_ecosysTracers_diagnostic_fields_level1
      logical, pointer :: config_ecosysTracers_diagnostic_fields_level2
      logical, pointer :: config_ecosysTracers_diagnostic_fields_level3
      logical, pointer :: config_ecosysTracers_diagnostic_fields_level4
      logical, pointer :: config_ecosysTracers_diagnostic_fields_level5

      logical, pointer :: config_use_DMSTracers
      logical, pointer :: config_use_DMSTracers_surface_bulk_forcing
      logical, pointer :: config_use_DMSTracers_surface_restoring
      logical, pointer :: config_use_DMSTracers_interior_restoring
      logical, pointer :: config_use_DMSTracers_exponential_decay
      logical, pointer :: config_use_DMSTracers_idealAge_forcing
      logical, pointer :: config_use_DMSTracers_ttd_forcing
      logical, pointer :: config_use_DMSTracers_surface_value
      logical, pointer :: config_use_DMSTracers_sea_ice_coupling

      logical, pointer :: config_use_MacroMoleculesTracers
      logical, pointer :: config_use_MacroMoleculesTracers_surface_bulk_forcing
      logical, pointer :: config_use_MacroMoleculesTracers_surface_restoring
      logical, pointer :: config_use_MacroMoleculesTracers_interior_restoring
      logical, pointer :: config_use_MacroMoleculesTracers_exponential_decay
      logical, pointer :: config_use_MacroMoleculesTracers_idealAge_forcing
      logical, pointer :: config_use_MacroMoleculesTracers_ttd_forcing
      logical, pointer :: config_use_MacroMoleculesTracers_surface_value
      logical, pointer :: config_use_MacroMoleculesTracers_sea_ice_coupling

      logical, pointer :: config_AM_globalStats_enable
      character (len=StrKIND), pointer :: config_AM_globalStats_compute_interval
      logical, pointer :: config_AM_globalStats_compute_on_startup
      logical, pointer :: config_AM_globalStats_write_on_startup
      logical, pointer :: config_AM_globalStats_text_file
      character (len=StrKIND), pointer :: config_AM_globalStats_directory
      character (len=StrKIND), pointer :: config_AM_globalStats_output_stream

      logical, pointer :: config_AM_surfaceAreaWeightedAverages_enable
      logical, pointer :: config_AM_surfaceAreaWeightedAverages_compute_on_startup
      logical, pointer :: config_AM_surfaceAreaWeightedAverages_write_on_startup
      character (len=StrKIND), pointer :: config_AM_surfaceAreaWeightedAverages_compute_interval
      character (len=StrKIND), pointer :: config_AM_surfaceAreaWeightedAverages_output_stream

      logical, pointer :: config_AM_waterMassCensus_enable
      character (len=StrKIND), pointer :: config_AM_waterMassCensus_compute_interval
      character (len=StrKIND), pointer :: config_AM_waterMassCensus_output_stream
      logical, pointer :: config_AM_waterMassCensus_compute_on_startup
      logical, pointer :: config_AM_waterMassCensus_write_on_startup
      real (kind=RKIND), pointer :: config_AM_waterMassCensus_minTemperature
      real (kind=RKIND), pointer :: config_AM_waterMassCensus_maxTemperature
      real (kind=RKIND), pointer :: config_AM_waterMassCensus_minSalinity
      real (kind=RKIND), pointer :: config_AM_waterMassCensus_maxSalinity
      logical, pointer :: config_AM_waterMassCensus_compute_predefined_regions
      character (len=StrKIND), pointer :: config_AM_waterMassCensus_region_group

      logical, pointer :: config_AM_layerVolumeWeightedAverage_enable
      character (len=StrKIND), pointer :: config_AM_layerVolumeWeightedAverage_compute_interval
      logical, pointer :: config_AM_layerVolumeWeightedAverage_compute_on_startup
      logical, pointer :: config_AM_layerVolumeWeightedAverage_write_on_startup
      character (len=StrKIND), pointer :: config_AM_layerVolumeWeightedAverage_output_stream

      logical, pointer :: config_AM_zonalMean_enable
      logical, pointer :: config_AM_zonalMean_compute_on_startup
      logical, pointer :: config_AM_zonalMean_write_on_startup
      character (len=StrKIND), pointer :: config_AM_zonalMean_compute_interval
      character (len=StrKIND), pointer :: config_AM_zonalMean_output_stream
      integer, pointer :: config_AM_zonalMean_num_bins
      real (kind=RKIND), pointer :: config_AM_zonalMean_min_bin
      real (kind=RKIND), pointer :: config_AM_zonalMean_max_bin

      logical, pointer :: config_AM_okuboWeiss_enable
      logical, pointer :: config_AM_okuboWeiss_compute_on_startup
      logical, pointer :: config_AM_okuboWeiss_write_on_startup
      character (len=StrKIND), pointer :: config_AM_okuboWeiss_compute_interval
      character (len=StrKIND), pointer :: config_AM_okuboWeiss_output_stream
      character (len=StrKIND), pointer :: config_AM_okuboWeiss_directory
      real (kind=RKIND), pointer :: config_AM_okuboWeiss_threshold_value
      real (kind=RKIND), pointer :: config_AM_okuboWeiss_normalization
      real (kind=RKIND), pointer :: config_AM_okuboWeiss_lambda2_normalization
      logical, pointer :: config_AM_okuboWeiss_use_lat_lon_coords
      logical, pointer :: config_AM_okuboWeiss_compute_eddy_census
      integer, pointer :: config_AM_okuboWeiss_eddy_min_cells

      logical, pointer :: config_AM_meridionalHeatTransport_enable
      character (len=StrKIND), pointer :: config_AM_meridionalHeatTransport_compute_interval
      logical, pointer :: config_AM_meridionalHeatTransport_compute_on_startup
      logical, pointer :: config_AM_meridionalHeatTransport_write_on_startup
      character (len=StrKIND), pointer :: config_AM_meridionalHeatTransport_output_stream
      integer, pointer :: config_AM_meridionalHeatTransport_num_bins
      real (kind=RKIND), pointer :: config_AM_meridionalHeatTransport_min_bin
      real (kind=RKIND), pointer :: config_AM_meridionalHeatTransport_max_bin
      character (len=StrKIND), pointer :: config_AM_meridionalHeatTransport_region_group

      logical, pointer :: config_AM_testComputeInterval_enable
      character (len=StrKIND), pointer :: config_AM_testComputeInterval_compute_interval
      logical, pointer :: config_AM_testComputeInterval_compute_on_startup
      logical, pointer :: config_AM_testComputeInterval_write_on_startup
      character (len=StrKIND), pointer :: config_AM_testComputeInterval_output_stream

      logical, pointer :: config_AM_highFrequencyOutput_enable
      character (len=StrKIND), pointer :: config_AM_highFrequencyOutput_compute_interval
      character (len=StrKIND), pointer :: config_AM_highFrequencyOutput_output_stream
      logical, pointer :: config_AM_highFrequencyOutput_compute_on_startup
      logical, pointer :: config_AM_highFrequencyOutput_write_on_startup

      logical, pointer :: config_AM_timeFilters_enable
      character (len=StrKIND), pointer :: config_AM_timeFilters_compute_interval
      character (len=StrKIND), pointer :: config_AM_timeFilters_output_stream
      character (len=StrKIND), pointer :: config_AM_timeFilters_restart_stream
      logical, pointer :: config_AM_timeFilters_compute_on_startup
      logical, pointer :: config_AM_timeFilters_write_on_startup
      logical, pointer :: config_AM_timeFilters_initialize_filters
      character (len=StrKIND), pointer :: config_AM_timeFilters_tau
      logical, pointer :: config_AM_timeFilters_compute_cell_centered_values

      logical, pointer :: config_AM_lagrPartTrack_enable
      character (len=StrKIND), pointer :: config_AM_lagrPartTrack_compute_interval
      logical, pointer :: config_AM_lagrPartTrack_compute_on_startup
      character (len=StrKIND), pointer :: config_AM_lagrPartTrack_output_stream
      character (len=StrKIND), pointer :: config_AM_lagrPartTrack_restart_stream
      character (len=StrKIND), pointer :: config_AM_lagrPartTrack_input_stream
      logical, pointer :: config_AM_lagrPartTrack_write_on_startup
      integer, pointer :: config_AM_lagrPartTrack_filter_number
      integer, pointer :: config_AM_lagrPartTrack_timeIntegration
      character (len=StrKIND), pointer :: config_AM_lagrPartTrack_reset_criteria
      character (len=StrKIND), pointer :: config_AM_lagrPartTrack_reset_global_timestamp
      character (len=StrKIND), pointer :: config_AM_lagrPartTrack_region_stream
      logical, pointer :: config_AM_lagrPartTrack_reset_if_outside_region
      logical, pointer :: config_AM_lagrPartTrack_reset_if_inside_region
      logical, pointer :: config_AM_lagrPartTrack_sample_horizontal_interp
      logical, pointer :: config_AM_lagrPartTrack_sample_temperature
      logical, pointer :: config_AM_lagrPartTrack_sample_salinity
      logical, pointer :: config_AM_lagrPartTrack_sample_DIC
      logical, pointer :: config_AM_lagrPartTrack_sample_ALK
      logical, pointer :: config_AM_lagrPartTrack_sample_PO4
      logical, pointer :: config_AM_lagrPartTrack_sample_NO3
      logical, pointer :: config_AM_lagrPartTrack_sample_SiO3
      logical, pointer :: config_AM_lagrPartTrack_sample_NH4
      logical, pointer :: config_AM_lagrPartTrack_sample_Fe
      logical, pointer :: config_AM_lagrPartTrack_sample_O2

      logical, pointer :: config_AM_eliassenPalm_enable
      character (len=StrKIND), pointer :: config_AM_eliassenPalm_compute_interval
      character (len=StrKIND), pointer :: config_AM_eliassenPalm_output_stream
      character (len=StrKIND), pointer :: config_AM_eliassenPalm_restart_stream
      logical, pointer :: config_AM_eliassenPalm_compute_on_startup
      logical, pointer :: config_AM_eliassenPalm_write_on_startup
      logical, pointer :: config_AM_eliassenPalm_debug
      integer, pointer :: config_AM_eliassenPalm_nBuoyancyLayers
      real (kind=RKIND), pointer :: config_AM_eliassenPalm_rhomin_buoycoor
      real (kind=RKIND), pointer :: config_AM_eliassenPalm_rhomax_buoycoor

      logical, pointer :: config_AM_mixedLayerDepths_enable
      character (len=StrKIND), pointer :: config_AM_mixedLayerDepths_compute_interval
      character (len=StrKIND), pointer :: config_AM_mixedLayerDepths_output_stream
      logical, pointer :: config_AM_mixedLayerDepths_write_on_startup
      logical, pointer :: config_AM_mixedLayerDepths_compute_on_startup
      logical, pointer :: config_AM_mixedLayerDepths_Tthreshold
      logical, pointer :: config_AM_mixedLayerDepths_Dthreshold
      real (kind=RKIND), pointer :: config_AM_mixedLayerDepths_crit_temp_threshold
      real (kind=RKIND), pointer :: config_AM_mixedLayerDepths_crit_dens_threshold
      real (kind=RKIND), pointer :: config_AM_mixedLayerDepths_reference_pressure
      logical, pointer :: config_AM_mixedLayerDepths_Tgradient
      logical, pointer :: config_AM_mixedLayerDepths_Dgradient
      real (kind=RKIND), pointer :: config_AM_mixedLayerDepths_temp_gradient_threshold
      real (kind=RKIND), pointer :: config_AM_mixedLayerDepths_den_gradient_threshold
      integer, pointer :: config_AM_mixedLayerDepths_interp_method

      logical, pointer :: config_AM_regionalStatsDaily_enable
      logical, pointer :: config_AM_regionalStatsDaily_compute_on_startup
      logical, pointer :: config_AM_regionalStatsDaily_write_on_startup
      character (len=StrKIND), pointer :: config_AM_regionalStatsDaily_compute_interval
      character (len=StrKIND), pointer :: config_AM_regionalStatsDaily_output_stream
      character (len=StrKIND), pointer :: config_AM_regionalStatsDaily_restart_stream
      character (len=StrKIND), pointer :: config_AM_regionalStatsDaily_input_stream
      character (len=StrKIND), pointer :: config_AM_regionalStatsDaily_operation
      character (len=StrKIND), pointer :: config_AM_regionalStatsDaily_region_type
      character (len=StrKIND), pointer :: config_AM_regionalStatsDaily_region_group
      character (len=StrKIND), pointer :: config_AM_regionalStatsDaily_1d_weighting_function
      character (len=StrKIND), pointer :: config_AM_regionalStatsDaily_2d_weighting_function
      character (len=StrKIND), pointer :: config_AM_regionalStatsDaily_1d_weighting_field
      character (len=StrKIND), pointer :: config_AM_regionalStatsDaily_2d_weighting_field
      character (len=StrKIND), pointer :: config_AM_regionalStatsDaily_vertical_mask
      character (len=StrKIND), pointer :: config_AM_regionalStatsDaily_vertical_dimension

      logical, pointer :: config_AM_regionalStatsWeekly_enable
      logical, pointer :: config_AM_regionalStatsWeekly_compute_on_startup
      logical, pointer :: config_AM_regionalStatsWeekly_write_on_startup
      character (len=StrKIND), pointer :: config_AM_regionalStatsWeekly_compute_interval
      character (len=StrKIND), pointer :: config_AM_regionalStatsWeekly_output_stream
      character (len=StrKIND), pointer :: config_AM_regionalStatsWeekly_restart_stream
      character (len=StrKIND), pointer :: config_AM_regionalStatsWeekly_input_stream
      character (len=StrKIND), pointer :: config_AM_regionalStatsWeekly_operation
      character (len=StrKIND), pointer :: config_AM_regionalStatsWeekly_region_type
      character (len=StrKIND), pointer :: config_AM_regionalStatsWeekly_region_group
      character (len=StrKIND), pointer :: config_AM_regionalStatsWeekly_1d_weighting_function
      character (len=StrKIND), pointer :: config_AM_regionalStatsWeekly_2d_weighting_function
      character (len=StrKIND), pointer :: config_AM_regionalStatsWeekly_1d_weighting_field
      character (len=StrKIND), pointer :: config_AM_regionalStatsWeekly_2d_weighting_field
      character (len=StrKIND), pointer :: config_AM_regionalStatsWeekly_vertical_mask
      character (len=StrKIND), pointer :: config_AM_regionalStatsWeekly_vertical_dimension

      logical, pointer :: config_AM_regionalStatsMonthly_enable
      logical, pointer :: config_AM_regionalStatsMonthly_compute_on_startup
      logical, pointer :: config_AM_regionalStatsMonthly_write_on_startup
      character (len=StrKIND), pointer :: config_AM_regionalStatsMonthly_compute_interval
      character (len=StrKIND), pointer :: config_AM_regionalStatsMonthly_output_stream
      character (len=StrKIND), pointer :: config_AM_regionalStatsMonthly_restart_stream
      character (len=StrKIND), pointer :: config_AM_regionalStatsMonthly_input_stream
      character (len=StrKIND), pointer :: config_AM_regionalStatsMonthly_operation
      character (len=StrKIND), pointer :: config_AM_regionalStatsMonthly_region_type
      character (len=StrKIND), pointer :: config_AM_regionalStatsMonthly_region_group
      character (len=StrKIND), pointer :: config_AM_regionalStatsMonthly_1d_weighting_function
      character (len=StrKIND), pointer :: config_AM_regionalStatsMonthly_2d_weighting_function
      character (len=StrKIND), pointer :: config_AM_regionalStatsMonthly_1d_weighting_field
      character (len=StrKIND), pointer :: config_AM_regionalStatsMonthly_2d_weighting_field
      character (len=StrKIND), pointer :: config_AM_regionalStatsMonthly_vertical_mask
      character (len=StrKIND), pointer :: config_AM_regionalStatsMonthly_vertical_dimension

      logical, pointer :: config_AM_regionalStatsCustom_enable
      logical, pointer :: config_AM_regionalStatsCustom_compute_on_startup
      logical, pointer :: config_AM_regionalStatsCustom_write_on_startup
      character (len=StrKIND), pointer :: config_AM_regionalStatsCustom_compute_interval
      character (len=StrKIND), pointer :: config_AM_regionalStatsCustom_output_stream
      character (len=StrKIND), pointer :: config_AM_regionalStatsCustom_restart_stream
      character (len=StrKIND), pointer :: config_AM_regionalStatsCustom_input_stream
      character (len=StrKIND), pointer :: config_AM_regionalStatsCustom_operation
      character (len=StrKIND), pointer :: config_AM_regionalStatsCustom_region_type
      character (len=StrKIND), pointer :: config_AM_regionalStatsCustom_region_group
      character (len=StrKIND), pointer :: config_AM_regionalStatsCustom_1d_weighting_function
      character (len=StrKIND), pointer :: config_AM_regionalStatsCustom_2d_weighting_function
      character (len=StrKIND), pointer :: config_AM_regionalStatsCustom_1d_weighting_field
      character (len=StrKIND), pointer :: config_AM_regionalStatsCustom_2d_weighting_field
      character (len=StrKIND), pointer :: config_AM_regionalStatsCustom_vertical_mask
      character (len=StrKIND), pointer :: config_AM_regionalStatsCustom_vertical_dimension

      logical, pointer :: config_AM_timeSeriesStatsDaily_enable
      logical, pointer :: config_AM_timeSeriesStatsDaily_compute_on_startup
      logical, pointer :: config_AM_timeSeriesStatsDaily_write_on_startup
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsDaily_compute_interval
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsDaily_output_stream
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsDaily_restart_stream
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsDaily_operation
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsDaily_reference_times
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsDaily_duration_intervals
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsDaily_repeat_intervals
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsDaily_reset_intervals
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsDaily_backward_output_offset

      logical, pointer :: config_AM_timeSeriesStatsMonthly_enable
      logical, pointer :: config_AM_timeSeriesStatsMonthly_compute_on_startup
      logical, pointer :: config_AM_timeSeriesStatsMonthly_write_on_startup
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthly_compute_interval
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthly_output_stream
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthly_restart_stream
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthly_operation
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthly_reference_times
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthly_duration_intervals
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthly_repeat_intervals
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthly_reset_intervals
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthly_backward_output_offset

      logical, pointer :: config_AM_timeSeriesStatsClimatology_enable
      logical, pointer :: config_AM_timeSeriesStatsClimatology_compute_on_startup
      logical, pointer :: config_AM_timeSeriesStatsClimatology_write_on_startup
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsClimatology_compute_interval
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsClimatology_output_stream
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsClimatology_restart_stream
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsClimatology_operation
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsClimatology_reference_times
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsClimatology_duration_intervals
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsClimatology_repeat_intervals
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsClimatology_reset_intervals
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsClimatology_backward_output_offset

      logical, pointer :: config_AM_timeSeriesStatsMonthlyMax_enable
      logical, pointer :: config_AM_timeSeriesStatsMonthlyMax_compute_on_startup
      logical, pointer :: config_AM_timeSeriesStatsMonthlyMax_write_on_startup
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthlyMax_compute_interval
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthlyMax_output_stream
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthlyMax_restart_stream
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthlyMax_operation
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthlyMax_reference_times
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthlyMax_duration_intervals
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthlyMax_repeat_intervals
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthlyMax_reset_intervals
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthlyMax_backward_output_offset

      logical, pointer :: config_AM_timeSeriesStatsMonthlyMin_enable
      logical, pointer :: config_AM_timeSeriesStatsMonthlyMin_compute_on_startup
      logical, pointer :: config_AM_timeSeriesStatsMonthlyMin_write_on_startup
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthlyMin_compute_interval
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthlyMin_output_stream
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthlyMin_restart_stream
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthlyMin_operation
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthlyMin_reference_times
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthlyMin_duration_intervals
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthlyMin_repeat_intervals
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthlyMin_reset_intervals
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsMonthlyMin_backward_output_offset

      logical, pointer :: config_AM_timeSeriesStatsCustom_enable
      logical, pointer :: config_AM_timeSeriesStatsCustom_compute_on_startup
      logical, pointer :: config_AM_timeSeriesStatsCustom_write_on_startup
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsCustom_compute_interval
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsCustom_output_stream
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsCustom_restart_stream
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsCustom_operation
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsCustom_reference_times
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsCustom_duration_intervals
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsCustom_repeat_intervals
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsCustom_reset_intervals
      character (len=StrKIND), pointer :: config_AM_timeSeriesStatsCustom_backward_output_offset

      logical, pointer :: config_AM_pointwiseStats_enable
      character (len=StrKIND), pointer :: config_AM_pointwiseStats_compute_interval
      character (len=StrKIND), pointer :: config_AM_pointwiseStats_output_stream
      logical, pointer :: config_AM_pointwiseStats_compute_on_startup
      logical, pointer :: config_AM_pointwiseStats_write_on_startup

      logical, pointer :: config_AM_debugDiagnostics_enable
      character (len=StrKIND), pointer :: config_AM_debugDiagnostics_compute_interval
      character (len=StrKIND), pointer :: config_AM_debugDiagnostics_output_stream
      logical, pointer :: config_AM_debugDiagnostics_compute_on_startup
      logical, pointer :: config_AM_debugDiagnostics_write_on_startup
      logical, pointer :: config_AM_debugDiagnostics_check_state

      logical, pointer :: config_AM_rpnCalculator_enable
      logical, pointer :: config_AM_rpnCalculator_compute_on_startup
      logical, pointer :: config_AM_rpnCalculator_write_on_startup
      character (len=StrKIND), pointer :: config_AM_rpnCalculator_compute_interval
      character (len=StrKIND), pointer :: config_AM_rpnCalculator_output_stream
      character (len=StrKIND), pointer :: config_AM_rpnCalculator_variable_a
      character (len=StrKIND), pointer :: config_AM_rpnCalculator_variable_b
      character (len=StrKIND), pointer :: config_AM_rpnCalculator_variable_c
      character (len=StrKIND), pointer :: config_AM_rpnCalculator_variable_d
      character (len=StrKIND), pointer :: config_AM_rpnCalculator_variable_e
      character (len=StrKIND), pointer :: config_AM_rpnCalculator_variable_f
      character (len=StrKIND), pointer :: config_AM_rpnCalculator_variable_g
      character (len=StrKIND), pointer :: config_AM_rpnCalculator_variable_h
      character (len=StrKIND), pointer :: config_AM_rpnCalculator_expression_1
      character (len=StrKIND), pointer :: config_AM_rpnCalculator_expression_2
      character (len=StrKIND), pointer :: config_AM_rpnCalculator_expression_3
      character (len=StrKIND), pointer :: config_AM_rpnCalculator_expression_4
      character (len=StrKIND), pointer :: config_AM_rpnCalculator_output_name_1
      character (len=StrKIND), pointer :: config_AM_rpnCalculator_output_name_2
      character (len=StrKIND), pointer :: config_AM_rpnCalculator_output_name_3
      character (len=StrKIND), pointer :: config_AM_rpnCalculator_output_name_4

      logical, pointer :: config_AM_transectTransport_enable
      character (len=StrKIND), pointer :: config_AM_transectTransport_compute_interval
      character (len=StrKIND), pointer :: config_AM_transectTransport_output_stream
      logical, pointer :: config_AM_transectTransport_compute_on_startup
      logical, pointer :: config_AM_transectTransport_write_on_startup
      character (len=StrKIND), pointer :: config_AM_transectTransport_transect_group

      logical, pointer :: config_AM_eddyProductVariables_enable
      character (len=StrKIND), pointer :: config_AM_eddyProductVariables_compute_interval
      character (len=StrKIND), pointer :: config_AM_eddyProductVariables_output_stream
      logical, pointer :: config_AM_eddyProductVariables_compute_on_startup
      logical, pointer :: config_AM_eddyProductVariables_write_on_startup

      logical, pointer :: config_AM_mocStreamfunction_enable
      character (len=StrKIND), pointer :: config_AM_mocStreamfunction_compute_interval
      character (len=StrKIND), pointer :: config_AM_mocStreamfunction_output_stream
      logical, pointer :: config_AM_mocStreamfunction_compute_on_startup
      logical, pointer :: config_AM_mocStreamfunction_write_on_startup
      real (kind=RKIND), pointer :: config_AM_mocStreamfunction_min_bin
      real (kind=RKIND), pointer :: config_AM_mocStreamfunction_max_bin
      integer, pointer :: config_AM_mocStreamfunction_num_bins
      character (len=StrKIND), pointer :: config_AM_mocStreamfunction_region_group
      character (len=StrKIND), pointer :: config_AM_mocStreamfunction_transect_group

      logical, pointer :: config_AM_oceanHeatContent_enable
      character (len=StrKIND), pointer :: config_AM_oceanHeatContent_compute_interval
      character (len=StrKIND), pointer :: config_AM_oceanHeatContent_output_stream
      logical, pointer :: config_AM_oceanHeatContent_compute_on_startup
      logical, pointer :: config_AM_oceanHeatContent_write_on_startup

      logical, pointer :: config_AM_mixedLayerHeatBudget_enable
      character (len=StrKIND), pointer :: config_AM_mixedLayerHeatBudget_compute_interval
      character (len=StrKIND), pointer :: config_AM_mixedLayerHeatBudget_output_stream
      logical, pointer :: config_AM_mixedLayerHeatBudget_compute_on_startup
      logical, pointer :: config_AM_mixedLayerHeatBudget_write_on_startup

      logical, pointer :: config_AM_sedimentFluxIndex_enable
      logical, pointer :: config_AM_sedimentFluxIndex_compute_on_startup
      logical, pointer :: config_AM_sedimentFluxIndex_write_on_startup
      character (len=StrKIND), pointer :: config_AM_sedimentFluxIndex_compute_interval
      character (len=StrKIND), pointer :: config_AM_sedimentFluxIndex_output_stream
      character (len=StrKIND), pointer :: config_AM_sedimentFluxIndex_directory
      logical, pointer :: config_AM_sedimentFluxIndex_use_lat_lon_coords

      logical, pointer :: config_AM_sedimentTransport_enable
      logical, pointer :: config_AM_sedimentTransport_compute_on_startup
      logical, pointer :: config_AM_sedimentTransport_write_on_startup
      character (len=StrKIND), pointer :: config_AM_sedimentTransport_compute_interval
      character (len=StrKIND), pointer :: config_AM_sedimentTransport_output_stream
      character (len=StrKIND), pointer :: config_AM_sedimentTransport_directory
      real (kind=RKIND), pointer :: config_AM_sedimentTransport_grain_size
      character (len=StrKIND), pointer :: config_AM_sedimentTransport_ws_formula
      character (len=StrKIND), pointer :: config_AM_sedimentTransport_bedld_formula
      character (len=StrKIND), pointer :: config_AM_sedimentTransport_SSC_ref_formula
      real (kind=RKIND), pointer :: config_AM_sedimentTransport_drag_coefficient
      real (kind=RKIND), pointer :: config_AM_sedimentTransport_erate
      real (kind=RKIND), pointer :: config_AM_sedimentTransport_tau_ce
      real (kind=RKIND), pointer :: config_AM_sedimentTransport_tau_cd
      real (kind=RKIND), pointer :: config_AM_sedimentTransport_Manning_coef
      real (kind=RKIND), pointer :: config_AM_sedimentTransport_grain_porosity
      real (kind=RKIND), pointer :: config_AM_sedimentTransport_water_density
      real (kind=RKIND), pointer :: config_AM_sedimentTransport_grain_density
      real (kind=RKIND), pointer :: config_AM_sedimentTransport_alpha
      real (kind=RKIND), pointer :: config_AM_sedimentTransport_kinematic_viscosity
      real (kind=RKIND), pointer :: config_AM_sedimentTransport_vertical_diffusion_coefficient
      logical, pointer :: config_AM_sedimentTransport_bedload
      logical, pointer :: config_AM_sedimentTransport_suspended
      logical, pointer :: config_AM_sedimentTransport_use_lat_lon_coords

      logical, pointer :: config_AM_harmonicAnalysis_enable
      character (len=StrKIND), pointer :: config_AM_harmonicAnalysis_compute_interval
      character (len=StrKIND), pointer :: config_AM_harmonicAnalysis_start
      character (len=StrKIND), pointer :: config_AM_harmonicAnalysis_end
      character (len=StrKIND), pointer :: config_AM_harmonicAnalysis_output_stream
      character (len=StrKIND), pointer :: config_AM_harmonicAnalysis_restart_stream
      logical, pointer :: config_AM_harmonicAnalysis_compute_on_startup
      logical, pointer :: config_AM_harmonicAnalysis_write_on_startup
      logical, pointer :: config_AM_harmonicAnalysis_use_M2
      logical, pointer :: config_AM_harmonicAnalysis_use_S2
      logical, pointer :: config_AM_harmonicAnalysis_use_N2
      logical, pointer :: config_AM_harmonicAnalysis_use_K2
      logical, pointer :: config_AM_harmonicAnalysis_use_K1
      logical, pointer :: config_AM_harmonicAnalysis_use_O1
      logical, pointer :: config_AM_harmonicAnalysis_use_Q1
      logical, pointer :: config_AM_harmonicAnalysis_use_P1

      integer, pointer :: config_baroclinic_channel_vert_levels
      logical, pointer :: config_baroclinic_channel_use_distances
      real (kind=RKIND), pointer :: config_baroclinic_channel_surface_temperature
      real (kind=RKIND), pointer :: config_baroclinic_channel_bottom_temperature
      real (kind=RKIND), pointer :: config_baroclinic_channel_temperature_difference
      real (kind=RKIND), pointer :: config_baroclinic_channel_gradient_width_frac
      real (kind=RKIND), pointer :: config_baroclinic_channel_gradient_width_dist
      real (kind=RKIND), pointer :: config_baroclinic_channel_bottom_depth
      real (kind=RKIND), pointer :: config_baroclinic_channel_salinity
      real (kind=RKIND), pointer :: config_baroclinic_channel_coriolis_parameter

      integer, pointer :: config_lock_exchange_vert_levels
      real (kind=RKIND), pointer :: config_lock_exchange_bottom_depth
      real (kind=RKIND), pointer :: config_lock_exchange_cold_temperature
      real (kind=RKIND), pointer :: config_lock_exchange_warm_temperature
      character (len=StrKIND), pointer :: config_lock_exchange_direction
      real (kind=RKIND), pointer :: config_lock_exchange_salinity
      character (len=StrKIND), pointer :: config_lock_exchange_layer_type
      real (kind=RKIND), pointer :: config_lock_exchange_isopycnal_min_thickness

      integer, pointer :: config_internal_waves_vert_levels
      logical, pointer :: config_internal_waves_use_distances
      real (kind=RKIND), pointer :: config_internal_waves_surface_temperature
      real (kind=RKIND), pointer :: config_internal_waves_bottom_temperature
      real (kind=RKIND), pointer :: config_internal_waves_temperature_difference
      real (kind=RKIND), pointer :: config_internal_waves_amplitude_width_frac
      real (kind=RKIND), pointer :: config_internal_waves_amplitude_width_dist
      real (kind=RKIND), pointer :: config_internal_waves_bottom_depth
      real (kind=RKIND), pointer :: config_internal_waves_salinity
      character (len=StrKIND), pointer :: config_internal_waves_layer_type
      real (kind=RKIND), pointer :: config_internal_waves_isopycnal_displacement

      integer, pointer :: config_overflow_vert_levels
      logical, pointer :: config_overflow_use_distances
      real (kind=RKIND), pointer :: config_overflow_bottom_depth
      real (kind=RKIND), pointer :: config_overflow_ridge_depth
      real (kind=RKIND), pointer :: config_overflow_plug_temperature
      real (kind=RKIND), pointer :: config_overflow_domain_temperature
      real (kind=RKIND), pointer :: config_overflow_salinity
      real (kind=RKIND), pointer :: config_overflow_plug_width_frac
      real (kind=RKIND), pointer :: config_overflow_slope_center_frac
      real (kind=RKIND), pointer :: config_overflow_slope_width_frac
      real (kind=RKIND), pointer :: config_overflow_plug_width_dist
      real (kind=RKIND), pointer :: config_overflow_slope_center_dist
      real (kind=RKIND), pointer :: config_overflow_slope_width_dist
      character (len=StrKIND), pointer :: config_overflow_layer_type
      real (kind=RKIND), pointer :: config_overflow_isopycnal_min_thickness

      integer, pointer :: config_dam_break_vert_levels
      real (kind=RKIND), pointer :: config_dam_break_eta0
      real (kind=RKIND), pointer :: config_dam_break_dc
      real (kind=RKIND), pointer :: config_dam_break_R0
      real (kind=RKIND), pointer :: config_dam_break_Xl
      real (kind=RKIND), pointer :: config_dam_break_Yl
      real (kind=RKIND), pointer :: config_dam_break_Inlet

      real (kind=RKIND), pointer :: config_global_ocean_minimum_depth
      character (len=StrKIND), pointer :: config_global_ocean_depth_file
      character (len=StrKIND), pointer :: config_global_ocean_depth_dimname
      character (len=StrKIND), pointer :: config_global_ocean_depth_varname
      real (kind=RKIND), pointer :: config_global_ocean_depth_conversion_factor
      character (len=StrKIND), pointer :: config_global_ocean_temperature_file
      character (len=StrKIND), pointer :: config_global_ocean_salinity_file
      character (len=StrKIND), pointer :: config_global_ocean_tracer_nlat_dimname
      character (len=StrKIND), pointer :: config_global_ocean_tracer_nlon_dimname
      character (len=StrKIND), pointer :: config_global_ocean_tracer_ndepth_dimname
      real (kind=RKIND), pointer :: config_global_ocean_tracer_depth_conversion_factor
      integer, pointer :: config_global_ocean_tracer_vert_levels
      character (len=StrKIND), pointer :: config_global_ocean_temperature_varname
      character (len=StrKIND), pointer :: config_global_ocean_salinity_varname
      logical, pointer :: config_global_ocean_tracer_latlon_degrees
      character (len=StrKIND), pointer :: config_global_ocean_tracer_lat_varname
      character (len=StrKIND), pointer :: config_global_ocean_tracer_lon_varname
      character (len=StrKIND), pointer :: config_global_ocean_tracer_depth_varname
      character (len=StrKIND), pointer :: config_global_ocean_tracer_method
      integer, pointer :: config_global_ocean_smooth_TS_iterations
      character (len=StrKIND), pointer :: config_global_ocean_swData_file
      character (len=StrKIND), pointer :: config_global_ocean_swData_nlat_dimname
      character (len=StrKIND), pointer :: config_global_ocean_swData_nlon_dimname
      character (len=StrKIND), pointer :: config_global_ocean_swData_lat_varname
      character (len=StrKIND), pointer :: config_global_ocean_swData_lon_varname
      logical, pointer :: config_global_ocean_swData_latlon_degrees
      character (len=StrKIND), pointer :: config_global_ocean_swData_method
      character (len=StrKIND), pointer :: config_global_ocean_chlorophyll_varname
      character (len=StrKIND), pointer :: config_global_ocean_zenithAngle_varname
      character (len=StrKIND), pointer :: config_global_ocean_clearSky_varname
      real (kind=RKIND), pointer :: config_global_ocean_piston_velocity
      real (kind=RKIND), pointer :: config_global_ocean_interior_restore_rate
      character (len=StrKIND), pointer :: config_global_ocean_topography_file
      character (len=StrKIND), pointer :: config_global_ocean_topography_nlat_dimname
      character (len=StrKIND), pointer :: config_global_ocean_topography_nlon_dimname
      logical, pointer :: config_global_ocean_topography_latlon_degrees
      character (len=StrKIND), pointer :: config_global_ocean_topography_lat_varname
      character (len=StrKIND), pointer :: config_global_ocean_topography_lon_varname
      character (len=StrKIND), pointer :: config_global_ocean_topography_varname
      logical, pointer :: config_global_ocean_topography_has_ocean_frac
      character (len=StrKIND), pointer :: config_global_ocean_topography_ocean_frac_varname
      character (len=StrKIND), pointer :: config_global_ocean_topography_method
      logical, pointer :: config_global_ocean_fill_bathymetry_holes
      integer, pointer :: config_global_ocean_topography_smooth_iterations
      real (kind=RKIND), pointer :: config_global_ocean_topography_smooth_weight
      logical, pointer :: config_global_ocean_deepen_critical_passages
      logical, pointer :: config_global_ocean_depress_by_land_ice
      character (len=StrKIND), pointer :: config_global_ocean_land_ice_topo_file
      character (len=StrKIND), pointer :: config_global_ocean_land_ice_topo_nlat_dimname
      character (len=StrKIND), pointer :: config_global_ocean_land_ice_topo_nlon_dimname
      logical, pointer :: config_global_ocean_land_ice_topo_latlon_degrees
      character (len=StrKIND), pointer :: config_global_ocean_land_ice_topo_lat_varname
      character (len=StrKIND), pointer :: config_global_ocean_land_ice_topo_lon_varname
      character (len=StrKIND), pointer :: config_global_ocean_land_ice_topo_thickness_varname
      character (len=StrKIND), pointer :: config_global_ocean_land_ice_topo_draft_varname
      character (len=StrKIND), pointer :: config_global_ocean_land_ice_topo_ice_frac_varname
      character (len=StrKIND), pointer :: config_global_ocean_land_ice_topo_grounded_frac_varname
      logical, pointer :: config_global_ocean_use_constant_land_ice_cavity_temperature
      real (kind=RKIND), pointer :: config_global_ocean_constant_land_ice_cavity_temperature
      logical, pointer :: config_global_ocean_cull_inland_seas
      character (len=StrKIND), pointer :: config_global_ocean_windstress_file
      character (len=StrKIND), pointer :: config_global_ocean_windstress_nlat_dimname
      character (len=StrKIND), pointer :: config_global_ocean_windstress_nlon_dimname
      logical, pointer :: config_global_ocean_windstress_latlon_degrees
      character (len=StrKIND), pointer :: config_global_ocean_windstress_lat_varname
      character (len=StrKIND), pointer :: config_global_ocean_windstress_lon_varname
      character (len=StrKIND), pointer :: config_global_ocean_windstress_zonal_varname
      character (len=StrKIND), pointer :: config_global_ocean_windstress_meridional_varname
      character (len=StrKIND), pointer :: config_global_ocean_windstress_method
      real (kind=RKIND), pointer :: config_global_ocean_windstress_conversion_factor
      character (len=StrKIND), pointer :: config_global_ocean_ecosys_file
      character (len=StrKIND), pointer :: config_global_ocean_ecosys_forcing_file
      character (len=StrKIND), pointer :: config_global_ocean_ecosys_nlat_dimname
      character (len=StrKIND), pointer :: config_global_ocean_ecosys_nlon_dimname
      character (len=StrKIND), pointer :: config_global_ocean_ecosys_ndepth_dimname
      real (kind=RKIND), pointer :: config_global_ocean_ecosys_depth_conversion_factor
      integer, pointer :: config_global_ocean_ecosys_vert_levels
      character (len=StrKIND), pointer :: config_global_ocean_ecosys_lat_varname
      character (len=StrKIND), pointer :: config_global_ocean_ecosys_lon_varname
      character (len=StrKIND), pointer :: config_global_ocean_ecosys_depth_varname
      logical, pointer :: config_global_ocean_ecosys_latlon_degrees
      character (len=StrKIND), pointer :: config_global_ocean_ecosys_method
      character (len=StrKIND), pointer :: config_global_ocean_ecosys_forcing_time_dimname
      integer, pointer :: config_global_ocean_smooth_ecosys_iterations

      integer, pointer :: config_cvmix_WSwSBF_vert_levels
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_surface_temperature
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_surface_salinity
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_surface_restoring_temperature
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_surface_restoring_salinity
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_temperature_piston_velocity
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_salinity_piston_velocity
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_sensible_heat_flux
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_latent_heat_flux
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_shortwave_heat_flux
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_rain_flux
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_evaporation_flux
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_interior_temperature_restoring_rate
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_interior_salinity_restoring_rate
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_temperature_gradient
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_salinity_gradient
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_temperature_gradient_mixed_layer
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_salinity_gradient_mixed_layer
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_mixed_layer_depth_temperature
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_mixed_layer_depth_salinity
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_mixed_layer_temperature_change
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_mixed_layer_salinity_change
      character (len=StrKIND), pointer :: config_cvmix_WSwSBF_vertical_grid
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_bottom_depth
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_max_windstress
      real (kind=RKIND), pointer :: config_cvmix_WSwSBF_coriolis_parameter

      integer, pointer :: config_iso_vert_levels
      real (kind=RKIND), pointer :: config_iso_main_channel_depth
      real (kind=RKIND), pointer :: config_iso_north_wall_lat
      real (kind=RKIND), pointer :: config_iso_south_wall_lat
      logical, pointer :: config_iso_ridge_flag
      real (kind=RKIND), pointer :: config_iso_ridge_center_lon
      real (kind=RKIND), pointer :: config_iso_ridge_height
      real (kind=RKIND), pointer :: config_iso_ridge_width
      logical, pointer :: config_iso_plateau_flag
      real (kind=RKIND), pointer :: config_iso_plateau_center_lon
      real (kind=RKIND), pointer :: config_iso_plateau_center_lat
      real (kind=RKIND), pointer :: config_iso_plateau_height
      real (kind=RKIND), pointer :: config_iso_plateau_radius
      real (kind=RKIND), pointer :: config_iso_plateau_slope_width
      logical, pointer :: config_iso_shelf_flag
      real (kind=RKIND), pointer :: config_iso_shelf_depth
      real (kind=RKIND), pointer :: config_iso_shelf_width
      logical, pointer :: config_iso_cont_slope_flag
      real (kind=RKIND), pointer :: config_iso_max_cont_slope
      logical, pointer :: config_iso_embayment_flag
      real (kind=RKIND), pointer :: config_iso_embayment_center_lon
      real (kind=RKIND), pointer :: config_iso_embayment_center_lat
      real (kind=RKIND), pointer :: config_iso_embayment_radius
      real (kind=RKIND), pointer :: config_iso_embayment_depth
      logical, pointer :: config_iso_depression_flag
      real (kind=RKIND), pointer :: config_iso_depression_center_lon
      real (kind=RKIND), pointer :: config_iso_depression_south_lat
      real (kind=RKIND), pointer :: config_iso_depression_north_lat
      real (kind=RKIND), pointer :: config_iso_depression_width
      real (kind=RKIND), pointer :: config_iso_depression_depth
      real (kind=RKIND), pointer :: config_iso_salinity
      real (kind=RKIND), pointer :: config_iso_wind_stress_max
      real (kind=RKIND), pointer :: config_iso_acc_wind
      real (kind=RKIND), pointer :: config_iso_asf_wind
      real (kind=RKIND), pointer :: config_iso_wind_trans
      real (kind=RKIND), pointer :: config_iso_heat_flux_south
      real (kind=RKIND), pointer :: config_iso_heat_flux_middle
      real (kind=RKIND), pointer :: config_iso_heat_flux_north
      real (kind=RKIND), pointer :: config_iso_heat_flux_lat_ss
      real (kind=RKIND), pointer :: config_iso_heat_flux_lat_sm
      real (kind=RKIND), pointer :: config_iso_heat_flux_lat_mn
      real (kind=RKIND), pointer :: config_iso_region1_center_lon
      real (kind=RKIND), pointer :: config_iso_region1_center_lat
      real (kind=RKIND), pointer :: config_iso_region2_center_lon
      real (kind=RKIND), pointer :: config_iso_region2_center_lat
      real (kind=RKIND), pointer :: config_iso_region3_center_lon
      real (kind=RKIND), pointer :: config_iso_region3_center_lat
      real (kind=RKIND), pointer :: config_iso_region4_center_lon
      real (kind=RKIND), pointer :: config_iso_region4_center_lat
      logical, pointer :: config_iso_heat_flux_region1_flag
      real (kind=RKIND), pointer :: config_iso_heat_flux_region1
      real (kind=RKIND), pointer :: config_iso_heat_flux_region1_radius
      logical, pointer :: config_iso_heat_flux_region2_flag
      real (kind=RKIND), pointer :: config_iso_heat_flux_region2
      real (kind=RKIND), pointer :: config_iso_heat_flux_region2_radius
      real (kind=RKIND), pointer :: config_iso_surface_temperature_piston_velocity
      real (kind=RKIND), pointer :: config_iso_initial_temp_t1
      real (kind=RKIND), pointer :: config_iso_initial_temp_t2
      real (kind=RKIND), pointer :: config_iso_initial_temp_h0
      real (kind=RKIND), pointer :: config_iso_initial_temp_h1
      real (kind=RKIND), pointer :: config_iso_initial_temp_mt
      real (kind=RKIND), pointer :: config_iso_initial_temp_latS
      real (kind=RKIND), pointer :: config_iso_initial_temp_latN
      real (kind=RKIND), pointer :: config_iso_temperature_sponge_t1
      real (kind=RKIND), pointer :: config_iso_temperature_sponge_h1
      real (kind=RKIND), pointer :: config_iso_temperature_sponge_l1
      real (kind=RKIND), pointer :: config_iso_temperature_sponge_tau1
      logical, pointer :: config_iso_temperature_restore_region1_flag
      real (kind=RKIND), pointer :: config_iso_temperature_restore_t1
      real (kind=RKIND), pointer :: config_iso_temperature_restore_lcx1
      real (kind=RKIND), pointer :: config_iso_temperature_restore_lcy1
      logical, pointer :: config_iso_temperature_restore_region2_flag
      real (kind=RKIND), pointer :: config_iso_temperature_restore_t2
      real (kind=RKIND), pointer :: config_iso_temperature_restore_lcx2
      real (kind=RKIND), pointer :: config_iso_temperature_restore_lcy2
      logical, pointer :: config_iso_temperature_restore_region3_flag
      real (kind=RKIND), pointer :: config_iso_temperature_restore_t3
      real (kind=RKIND), pointer :: config_iso_temperature_restore_lcx3
      real (kind=RKIND), pointer :: config_iso_temperature_restore_lcy3
      logical, pointer :: config_iso_temperature_restore_region4_flag
      real (kind=RKIND), pointer :: config_iso_temperature_restore_t4
      real (kind=RKIND), pointer :: config_iso_temperature_restore_lcx4
      real (kind=RKIND), pointer :: config_iso_temperature_restore_lcy4

      integer, pointer :: config_soma_vert_levels
      real (kind=RKIND), pointer :: config_soma_domain_width
      real (kind=RKIND), pointer :: config_soma_center_latitude
      real (kind=RKIND), pointer :: config_soma_center_longitude
      real (kind=RKIND), pointer :: config_soma_phi
      real (kind=RKIND), pointer :: config_soma_bottom_depth
      real (kind=RKIND), pointer :: config_soma_shelf_width
      real (kind=RKIND), pointer :: config_soma_shelf_depth
      real (kind=RKIND), pointer :: config_soma_ref_density
      real (kind=RKIND), pointer :: config_soma_density_difference
      real (kind=RKIND), pointer :: config_soma_thermocline_depth
      real (kind=RKIND), pointer :: config_soma_density_difference_linear
      real (kind=RKIND), pointer :: config_soma_surface_temperature
      real (kind=RKIND), pointer :: config_soma_surface_salinity
      logical, pointer :: config_soma_use_surface_temp_restoring
      real (kind=RKIND), pointer :: config_soma_surface_temp_restoring_at_center_latitude
      real (kind=RKIND), pointer :: config_soma_surface_temp_restoring_latitude_gradient
      real (kind=RKIND), pointer :: config_soma_restoring_temp_piston_vel

      integer, pointer :: config_ziso_vert_levels
      logical, pointer :: config_ziso_add_easterly_wind_stress_ASF
      real (kind=RKIND), pointer :: config_ziso_wind_transition_position
      real (kind=RKIND), pointer :: config_ziso_antarctic_shelf_front_width
      real (kind=RKIND), pointer :: config_ziso_wind_stress_shelf_front_max
      logical, pointer :: config_ziso_use_slopping_bathymetry
      real (kind=RKIND), pointer :: config_ziso_meridional_extent
      real (kind=RKIND), pointer :: config_ziso_zonal_extent
      real (kind=RKIND), pointer :: config_ziso_bottom_depth
      real (kind=RKIND), pointer :: config_ziso_shelf_depth
      real (kind=RKIND), pointer :: config_ziso_slope_half_width
      real (kind=RKIND), pointer :: config_ziso_slope_center_position
      real (kind=RKIND), pointer :: config_ziso_reference_coriolis
      real (kind=RKIND), pointer :: config_ziso_coriolis_gradient
      real (kind=RKIND), pointer :: config_ziso_wind_stress_max
      real (kind=RKIND), pointer :: config_ziso_mean_restoring_temp
      real (kind=RKIND), pointer :: config_ziso_restoring_temp_dev_ta
      real (kind=RKIND), pointer :: config_ziso_restoring_temp_dev_tb
      real (kind=RKIND), pointer :: config_ziso_restoring_temp_tau
      real (kind=RKIND), pointer :: config_ziso_restoring_temp_piston_vel
      real (kind=RKIND), pointer :: config_ziso_restoring_temp_ze
      real (kind=RKIND), pointer :: config_ziso_restoring_sponge_l
      real (kind=RKIND), pointer :: config_ziso_initial_temp_t1
      real (kind=RKIND), pointer :: config_ziso_initial_temp_t2
      real (kind=RKIND), pointer :: config_ziso_initial_temp_h1
      real (kind=RKIND), pointer :: config_ziso_initial_temp_mt
      logical, pointer :: config_ziso_frazil_enable
      real (kind=RKIND), pointer :: config_ziso_frazil_temperature_anomaly

      integer, pointer :: config_sub_ice_shelf_2D_vert_levels
      real (kind=RKIND), pointer :: config_sub_ice_shelf_2D_bottom_depth
      real (kind=RKIND), pointer :: config_sub_ice_shelf_2D_cavity_thickness
      real (kind=RKIND), pointer :: config_sub_ice_shelf_2D_slope_height
      real (kind=RKIND), pointer :: config_sub_ice_shelf_2D_edge_width
      real (kind=RKIND), pointer :: config_sub_ice_shelf_2D_y1
      real (kind=RKIND), pointer :: config_sub_ice_shelf_2D_y2
      real (kind=RKIND), pointer :: config_sub_ice_shelf_2D_temperature
      real (kind=RKIND), pointer :: config_sub_ice_shelf_2D_surface_salinity
      real (kind=RKIND), pointer :: config_sub_ice_shelf_2D_bottom_salinity

      integer, pointer :: config_periodic_planar_vert_levels
      real (kind=RKIND), pointer :: config_periodic_planar_bottom_depth
      real (kind=RKIND), pointer :: config_periodic_planar_velocity_strength

      integer, pointer :: config_ecosys_column_vert_levels
      character (len=StrKIND), pointer :: config_ecosys_column_vertical_grid
      character (len=StrKIND), pointer :: config_ecosys_column_TS_filename
      character (len=StrKIND), pointer :: config_ecosys_column_ecosys_filename
      real (kind=RKIND), pointer :: config_ecosys_column_bottom_depth

      integer, pointer :: config_sea_mount_vert_levels
      character (len=StrKIND), pointer :: config_sea_mount_layer_type
      character (len=StrKIND), pointer :: config_sea_mount_stratification_type
      real (kind=RKIND), pointer :: config_sea_mount_density_coef_linear
      real (kind=RKIND), pointer :: config_sea_mount_density_coef_exp
      real (kind=RKIND), pointer :: config_sea_mount_density_gradient_linear
      real (kind=RKIND), pointer :: config_sea_mount_density_gradient_exp
      real (kind=RKIND), pointer :: config_sea_mount_density_depth_linear
      real (kind=RKIND), pointer :: config_sea_mount_density_depth_exp
      real (kind=RKIND), pointer :: config_sea_mount_density_ref
      real (kind=RKIND), pointer :: config_sea_mount_density_Tref
      real (kind=RKIND), pointer :: config_sea_mount_density_alpha
      real (kind=RKIND), pointer :: config_sea_mount_bottom_depth
      real (kind=RKIND), pointer :: config_sea_mount_height
      real (kind=RKIND), pointer :: config_sea_mount_radius
      real (kind=RKIND), pointer :: config_sea_mount_width
      real (kind=RKIND), pointer :: config_sea_mount_salinity
      real (kind=RKIND), pointer :: config_sea_mount_coriolis_parameter

      integer, pointer :: config_isomip_vert_levels
      character (len=StrKIND), pointer :: config_isomip_vertical_level_distribution
      real (kind=RKIND), pointer :: config_isomip_bottom_depth
      real (kind=RKIND), pointer :: config_isomip_temperature
      real (kind=RKIND), pointer :: config_isomip_salinity
      real (kind=RKIND), pointer :: config_isomip_restoring_temperature
      real (kind=RKIND), pointer :: config_isomip_temperature_piston_velocity
      real (kind=RKIND), pointer :: config_isomip_restoring_salinity
      real (kind=RKIND), pointer :: config_isomip_salinity_piston_velocity
      real (kind=RKIND), pointer :: config_isomip_coriolis_parameter
      real (kind=RKIND), pointer :: config_isomip_southern_boundary
      real (kind=RKIND), pointer :: config_isomip_northern_boundary
      real (kind=RKIND), pointer :: config_isomip_western_boundary
      real (kind=RKIND), pointer :: config_isomip_eastern_boundary
      real (kind=RKIND), pointer :: config_isomip_y1
      real (kind=RKIND), pointer :: config_isomip_z1
      real (kind=RKIND), pointer :: config_isomip_ice_fraction1
      real (kind=RKIND), pointer :: config_isomip_y2
      real (kind=RKIND), pointer :: config_isomip_z2
      real (kind=RKIND), pointer :: config_isomip_ice_fraction2
      real (kind=RKIND), pointer :: config_isomip_y3
      real (kind=RKIND), pointer :: config_isomip_z3
      real (kind=RKIND), pointer :: config_isomip_ice_fraction3

      integer, pointer :: config_isomip_plus_vert_levels
      character (len=StrKIND), pointer :: config_isomip_plus_vertical_level_distribution
      real (kind=RKIND), pointer :: config_isomip_plus_max_bottom_depth
      integer, pointer :: config_isomip_plus_minimum_levels
      real (kind=RKIND), pointer :: config_isomip_plus_min_column_thickness
      real (kind=RKIND), pointer :: config_isomip_plus_min_ocean_fraction
      character (len=StrKIND), pointer :: config_isomip_plus_topography_file
      real (kind=RKIND), pointer :: config_isomip_plus_init_top_temp
      real (kind=RKIND), pointer :: config_isomip_plus_init_bot_temp
      real (kind=RKIND), pointer :: config_isomip_plus_init_top_sal
      real (kind=RKIND), pointer :: config_isomip_plus_init_bot_sal
      real (kind=RKIND), pointer :: config_isomip_plus_restore_top_temp
      real (kind=RKIND), pointer :: config_isomip_plus_restore_bot_temp
      real (kind=RKIND), pointer :: config_isomip_plus_restore_top_sal
      real (kind=RKIND), pointer :: config_isomip_plus_restore_bot_sal
      real (kind=RKIND), pointer :: config_isomip_plus_restore_rate
      real (kind=RKIND), pointer :: config_isomip_plus_restore_evap_rate
      real (kind=RKIND), pointer :: config_isomip_plus_restore_xMin
      real (kind=RKIND), pointer :: config_isomip_plus_restore_xMax
      real (kind=RKIND), pointer :: config_isomip_plus_coriolis_parameter
      real (kind=RKIND), pointer :: config_isomip_plus_effective_density

      integer, pointer :: config_hurricane_vert_levels
      real (kind=RKIND), pointer :: config_hurricane_min_depth
      real (kind=RKIND), pointer :: config_hurricane_max_depth
      real (kind=RKIND), pointer :: config_hurricane_gaussian_hump_amplitude
      logical, pointer :: config_hurricane_use_gaussian_hump
      real (kind=RKIND), pointer :: config_hurricane_gaussian_lon_center
      real (kind=RKIND), pointer :: config_hurricane_gaussian_lat_center
      real (kind=RKIND), pointer :: config_hurricane_gaussian_width
      real (kind=RKIND), pointer :: config_hurricane_gaussian_amplitude
      real (kind=RKIND), pointer :: config_hurricane_gaussian_slr_amp
      real (kind=RKIND), pointer :: config_hurricane_land_z_limit
      real (kind=RKIND), pointer :: config_hurricane_marsh_z_limit
      real (kind=RKIND), pointer :: config_hurricane_land_drag
      real (kind=RKIND), pointer :: config_hurricane_marsh_drag
      real (kind=RKIND), pointer :: config_hurricane_channel_drag
      real (kind=RKIND), pointer :: config_hurricane_sea_level_rise_adjustment

      integer, pointer :: config_tidal_boundary_vert_levels
      integer, pointer :: config_tidal_boundary_min_vert_levels
      character (len=StrKIND), pointer :: config_tidal_boundary_layer_type
      real (kind=RKIND), pointer :: config_tidal_boundary_right_bottom_depth
      logical, pointer :: config_tidal_start_dry
      logical, pointer :: config_tidal_boundary_use_distances
      real (kind=RKIND), pointer :: config_tidal_boundary_left_value
      real (kind=RKIND), pointer :: config_tidal_boundary_right_value
      real (kind=RKIND), pointer :: config_tidal_boundary_left_bottom_depth
      real (kind=RKIND), pointer :: config_tidal_boundary_salinity
      real (kind=RKIND), pointer :: config_tidal_boundary_domain_temperature
      real (kind=RKIND), pointer :: config_tidal_boundary_plug_temperature
      real (kind=RKIND), pointer :: config_tidal_boundary_plug_width_frac
      real (kind=RKIND), pointer :: config_tidal_forcing_left_Cd_or_n
      real (kind=RKIND), pointer :: config_tidal_forcing_right_Cd_or_n
      logical, pointer :: config_use_idealized_transect
      real (kind=RKIND), pointer :: config_idealized_transect_Lshore
      real (kind=RKIND), pointer :: config_idealized_transect_Sshore
      real (kind=RKIND), pointer :: config_idealized_transect_Lcoast
      real (kind=RKIND), pointer :: config_idealized_transect_Scoast
      real (kind=RKIND), pointer :: config_idealized_transect_Lmarsh
      real (kind=RKIND), pointer :: config_idealized_transect_Smarsh
      real (kind=RKIND), pointer :: config_idealized_transect_roughness
      real (kind=RKIND), pointer :: config_idealized_transect_roughness_marsh
      real (kind=RKIND), pointer :: config_idealized_vegetation_diameter
      real (kind=RKIND), pointer :: config_idealized_vegetation_height
      real (kind=RKIND), pointer :: config_idealized_vegetation_density

      real (kind=RKIND), pointer :: config_cosine_bell_temperature
      real (kind=RKIND), pointer :: config_cosine_bell_salinity
      real (kind=RKIND), pointer :: config_cosine_bell_lat_center
      real (kind=RKIND), pointer :: config_cosine_bell_lon_center
      real (kind=RKIND), pointer :: config_cosine_bell_psi0
      real (kind=RKIND), pointer :: config_cosine_bell_radius
      real (kind=RKIND), pointer :: config_cosine_bell_vel_pd

      integer, pointer :: config_mixed_layer_eddy_vert_levels
      real (kind=RKIND), pointer :: config_mixed_layer_eddy_bottom_depth
      real (kind=RKIND), pointer :: config_mixed_layer_eddy_mixed_layer_depth
      real (kind=RKIND), pointer :: config_mixed_layer_eddy_base_temperature
      real (kind=RKIND), pointer :: config_mixed_layer_eddy_temperature_stratification_mixed_layer
      real (kind=RKIND), pointer :: config_mixed_layer_eddy_temperature_stratification_interior
      real (kind=RKIND), pointer :: config_mixed_layer_eddy_temperature_horizontal_gradient
      real (kind=RKIND), pointer :: config_mixed_layer_eddy_temperature_front_width
      real (kind=RKIND), pointer :: config_mixed_layer_eddy_temperature_perturbation_magnitude
      real (kind=RKIND), pointer :: config_mixed_layer_eddy_salinity
      logical, pointer :: config_mixed_layer_eddy_two_fronts
      real (kind=RKIND), pointer :: config_mixed_layer_eddy_restoring_width
      real (kind=RKIND), pointer :: config_mixed_layer_eddy_restoring_tau
      real (kind=RKIND), pointer :: config_mixed_layer_eddy_heat_flux
      real (kind=RKIND), pointer :: config_mixed_layer_eddy_evaporation_flux
      real (kind=RKIND), pointer :: config_mixed_layer_eddy_wind_stress_zonal
      real (kind=RKIND), pointer :: config_mixed_layer_eddy_wind_stress_meridional
      real (kind=RKIND), pointer :: config_mixed_layer_eddy_coriolis_parameter

      logical, pointer :: config_alter_ICs_for_pcs
      character (len=StrKIND), pointer :: config_pc_alteration_type
      real (kind=RKIND), pointer :: config_min_pc_fraction

      character (len=StrKIND), pointer :: config_init_configuration
      logical, pointer :: config_expand_sphere
      logical, pointer :: config_realistic_coriolis_parameter
      logical, pointer :: config_write_cull_cell_mask
      character (len=StrKIND), pointer :: config_vertical_grid

      real (kind=RKIND), pointer :: config_1dCVTgenerator_stretch1
      real (kind=RKIND), pointer :: config_1dCVTgenerator_stretch2
      real (kind=RKIND), pointer :: config_1dCVTgenerator_dzSeed

      character (len=StrKIND), pointer :: config_init_vertical_grid_type

      integer, pointer :: config_rx1_outer_iter_count
      integer, pointer :: config_rx1_inner_iter_count
      real (kind=RKIND), pointer :: config_rx1_init_inner_weight
      real (kind=RKIND), pointer :: config_rx1_max
      real (kind=RKIND), pointer :: config_rx1_horiz_smooth_weight
      real (kind=RKIND), pointer :: config_rx1_vert_smooth_weight
      real (kind=RKIND), pointer :: config_rx1_slope_weight
      real (kind=RKIND), pointer :: config_rx1_zstar_weight
      integer, pointer :: config_rx1_horiz_smooth_open_ocean_cells
      integer, pointer :: config_rx1_min_levels
      real (kind=RKIND), pointer :: config_rx1_min_layer_thickness


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_config_init
!
!> \brief   Initializes the ocean config
!> \details
!>  This routine sets up config for use in the ocean model.
!
!-----------------------------------------------------------------------
   subroutine ocn_config_init(configPool)!{{{
       type (mpas_pool_type), pointer :: configPool

      call mpas_pool_get_config(configPool, 'config_ocean_run_mode', config_ocean_run_mode)

      call mpas_pool_get_config(configPool, 'config_do_restart', config_do_restart)
      call mpas_pool_get_config(configPool, 'config_restart_timestamp_name', config_restart_timestamp_name)
      call mpas_pool_get_config(configPool, 'config_start_time', config_start_time)
      call mpas_pool_get_config(configPool, 'config_stop_time', config_stop_time)
      call mpas_pool_get_config(configPool, 'config_run_duration', config_run_duration)
      call mpas_pool_get_config(configPool, 'config_calendar_type', config_calendar_type)

      call mpas_pool_get_config(configPool, 'config_write_output_on_startup', config_write_output_on_startup)
      call mpas_pool_get_config(configPool, 'config_pio_num_iotasks', config_pio_num_iotasks)
      call mpas_pool_get_config(configPool, 'config_pio_stride', config_pio_stride)

      call mpas_pool_get_config(configPool, 'config_num_halos', config_num_halos)
      call mpas_pool_get_config(configPool, 'config_block_decomp_file_prefix', config_block_decomp_file_prefix)
      call mpas_pool_get_config(configPool, 'config_number_of_blocks', config_number_of_blocks)
      call mpas_pool_get_config(configPool, 'config_explicit_proc_decomp', config_explicit_proc_decomp)
      call mpas_pool_get_config(configPool, 'config_proc_decomp_file_prefix', config_proc_decomp_file_prefix)

      call mpas_pool_get_config(configPool, 'config_dt', config_dt)
      call mpas_pool_get_config(configPool, 'config_time_integrator', config_time_integrator)

      call mpas_pool_get_config(configPool, 'config_hmix_scaleWithMesh', config_hmix_scaleWithMesh)
      call mpas_pool_get_config(configPool, 'config_maxMeshDensity', config_maxMeshDensity)
      call mpas_pool_get_config(configPool, 'config_hmix_use_ref_cell_width', config_hmix_use_ref_cell_width)
      call mpas_pool_get_config(configPool, 'config_hmix_ref_cell_width', config_hmix_ref_cell_width)
      call mpas_pool_get_config(configPool, 'config_apvm_scale_factor', config_apvm_scale_factor)

      call mpas_pool_get_config(configPool, 'config_use_mom_del2', config_use_mom_del2)
      call mpas_pool_get_config(configPool, 'config_mom_del2', config_mom_del2)
      call mpas_pool_get_config(configPool, 'config_use_tracer_del2', config_use_tracer_del2)
      call mpas_pool_get_config(configPool, 'config_tracer_del2', config_tracer_del2)

      call mpas_pool_get_config(configPool, 'config_use_mom_del4', config_use_mom_del4)
      call mpas_pool_get_config(configPool, 'config_mom_del4', config_mom_del4)
      call mpas_pool_get_config(configPool, 'config_mom_del4_div_factor', config_mom_del4_div_factor)
      call mpas_pool_get_config(configPool, 'config_use_tracer_del4', config_use_tracer_del4)
      call mpas_pool_get_config(configPool, 'config_tracer_del4', config_tracer_del4)

      call mpas_pool_get_config(configPool, 'config_use_Leith_del2', config_use_Leith_del2)
      call mpas_pool_get_config(configPool, 'config_Leith_parameter', config_Leith_parameter)
      call mpas_pool_get_config(configPool, 'config_Leith_dx', config_Leith_dx)
      call mpas_pool_get_config(configPool, 'config_Leith_visc2_max', config_Leith_visc2_max)

      call mpas_pool_get_config(configPool, 'config_eddying_resolution_taper', config_eddying_resolution_taper)
      call mpas_pool_get_config(configPool, 'config_eddying_resolution_ramp_min', config_eddying_resolution_ramp_min)
      call mpas_pool_get_config(configPool, 'config_eddying_resolution_ramp_max', config_eddying_resolution_ramp_max)

      call mpas_pool_get_config(configPool, 'config_use_Redi', config_use_Redi)
      call mpas_pool_get_config(configPool, 'config_Redi_closure', config_Redi_closure)
      call mpas_pool_get_config(configPool, 'config_Redi_constant_kappa', config_Redi_constant_kappa)
      call mpas_pool_get_config(configPool, 'config_Redi_maximum_slope', config_Redi_maximum_slope)
      call mpas_pool_get_config(configPool, 'config_Redi_use_slope_taper', config_Redi_use_slope_taper)
      call mpas_pool_get_config(configPool, 'config_Redi_use_surface_taper', config_Redi_use_surface_taper)
      call mpas_pool_get_config(configPool, 'config_Redi_N2_based_taper_enable', config_Redi_N2_based_taper_enable)
      call mpas_pool_get_config(configPool, 'config_Redi_N2_based_taper_min', config_Redi_N2_based_taper_min)
      call mpas_pool_get_config(configPool, 'config_Redi_N2_based_taper_limit_term1', config_Redi_N2_based_taper_limit_term1)

      call mpas_pool_get_config(configPool, 'config_use_GM', config_use_GM)
      call mpas_pool_get_config(configPool, 'config_GM_closure', config_GM_closure)
      call mpas_pool_get_config(configPool, 'config_GM_constant_kappa', config_GM_constant_kappa)
      call mpas_pool_get_config(configPool, 'config_GM_constant_gravWaveSpeed', config_GM_constant_gravWaveSpeed)
      call mpas_pool_get_config(configPool, 'config_GM_spatially_variable_min_kappa', config_GM_spatially_variable_min_kappa)
      call mpas_pool_get_config(configPool, 'config_GM_spatially_variable_max_kappa', config_GM_spatially_variable_max_kappa)
      call mpas_pool_get_config(configPool, 'config_GM_spatially_variable_baroclinic_mode', &
config_GM_spatially_variable_baroclinic_mode)
      call mpas_pool_get_config(configPool, 'config_GM_Visbeck_alpha', config_GM_Visbeck_alpha)
      call mpas_pool_get_config(configPool, 'config_GM_Visbeck_max_depth', config_GM_Visbeck_max_depth)
      call mpas_pool_get_config(configPool, 'config_GM_EG_riMin', config_GM_EG_riMin)
      call mpas_pool_get_config(configPool, 'config_GM_EG_kappa_factor', config_GM_EG_kappa_factor)
      call mpas_pool_get_config(configPool, 'config_GM_EG_Rossby_factor', config_GM_EG_Rossby_factor)
      call mpas_pool_get_config(configPool, 'config_GM_EG_Rhines_factor', config_GM_EG_Rhines_factor)

      call mpas_pool_get_config(configPool, 'config_Rayleigh_friction', config_Rayleigh_friction)
      call mpas_pool_get_config(configPool, 'config_Rayleigh_damping_coeff', config_Rayleigh_damping_coeff)
      call mpas_pool_get_config(configPool, 'config_Rayleigh_damping_depth_variable', config_Rayleigh_damping_depth_variable)
      call mpas_pool_get_config(configPool, 'config_Rayleigh_bottom_friction', config_Rayleigh_bottom_friction)
      call mpas_pool_get_config(configPool, 'config_Rayleigh_bottom_damping_coeff', config_Rayleigh_bottom_damping_coeff)

      call mpas_pool_get_config(configPool, 'config_use_cvmix', config_use_cvmix)
      call mpas_pool_get_config(configPool, 'config_cvmix_prandtl_number', config_cvmix_prandtl_number)
      call mpas_pool_get_config(configPool, 'config_cvmix_background_scheme', config_cvmix_background_scheme)
      call mpas_pool_get_config(configPool, 'config_cvmix_background_diffusion', config_cvmix_background_diffusion)
      call mpas_pool_get_config(configPool, 'config_cvmix_background_viscosity', config_cvmix_background_viscosity)
      call mpas_pool_get_config(configPool, 'config_cvmix_BryanLewis_bl1', config_cvmix_BryanLewis_bl1)
      call mpas_pool_get_config(configPool, 'config_cvmix_BryanLewis_bl2', config_cvmix_BryanLewis_bl2)
      call mpas_pool_get_config(configPool, 'config_cvmix_BryanLewis_transitionDepth', config_cvmix_BryanLewis_transitionDepth)
      call mpas_pool_get_config(configPool, 'config_cvmix_BryanLewis_transitionWidth', config_cvmix_BryanLewis_transitionWidth)
      call mpas_pool_get_config(configPool, 'config_use_cvmix_convection', config_use_cvmix_convection)
      call mpas_pool_get_config(configPool, 'config_cvmix_convective_diffusion', config_cvmix_convective_diffusion)
      call mpas_pool_get_config(configPool, 'config_cvmix_convective_viscosity', config_cvmix_convective_viscosity)
      call mpas_pool_get_config(configPool, 'config_cvmix_convective_basedOnBVF', config_cvmix_convective_basedOnBVF)
      call mpas_pool_get_config(configPool, 'config_cvmix_convective_triggerBVF', config_cvmix_convective_triggerBVF)
      call mpas_pool_get_config(configPool, 'config_use_cvmix_shear', config_use_cvmix_shear)
      call mpas_pool_get_config(configPool, 'config_cvmix_num_ri_smooth_loops', config_cvmix_num_ri_smooth_loops)
      call mpas_pool_get_config(configPool, 'config_cvmix_use_BLD_smoothing', config_cvmix_use_BLD_smoothing)
      call mpas_pool_get_config(configPool, 'config_cvmix_shear_mixing_scheme', config_cvmix_shear_mixing_scheme)
      call mpas_pool_get_config(configPool, 'config_cvmix_shear_PP_nu_zero', config_cvmix_shear_PP_nu_zero)
      call mpas_pool_get_config(configPool, 'config_cvmix_shear_PP_alpha', config_cvmix_shear_PP_alpha)
      call mpas_pool_get_config(configPool, 'config_cvmix_shear_PP_exp', config_cvmix_shear_PP_exp)
      call mpas_pool_get_config(configPool, 'config_cvmix_shear_KPP_nu_zero', config_cvmix_shear_KPP_nu_zero)
      call mpas_pool_get_config(configPool, 'config_cvmix_shear_KPP_Ri_zero', config_cvmix_shear_KPP_Ri_zero)
      call mpas_pool_get_config(configPool, 'config_cvmix_shear_KPP_exp', config_cvmix_shear_KPP_exp)
      call mpas_pool_get_config(configPool, 'config_use_cvmix_tidal_mixing', config_use_cvmix_tidal_mixing)
      call mpas_pool_get_config(configPool, 'config_use_cvmix_double_diffusion', config_use_cvmix_double_diffusion)
      call mpas_pool_get_config(configPool, 'config_use_cvmix_kpp', config_use_cvmix_kpp)
      call mpas_pool_get_config(configPool, 'config_use_cvmix_fixed_boundary_layer', config_use_cvmix_fixed_boundary_layer)
      call mpas_pool_get_config(configPool, 'config_cvmix_kpp_boundary_layer_depth', config_cvmix_kpp_boundary_layer_depth)
      call mpas_pool_get_config(configPool, 'config_cvmix_kpp_criticalBulkRichardsonNumber', &
config_cvmix_kpp_criticalBulkRichardsonNumber)
      call mpas_pool_get_config(configPool, 'config_cvmix_kpp_matching', config_cvmix_kpp_matching)
      call mpas_pool_get_config(configPool, 'config_cvmix_kpp_EkmanOBL', config_cvmix_kpp_EkmanOBL)
      call mpas_pool_get_config(configPool, 'config_cvmix_kpp_MonObOBL', config_cvmix_kpp_MonObOBL)
      call mpas_pool_get_config(configPool, 'config_cvmix_kpp_interpolationOMLType', config_cvmix_kpp_interpolationOMLType)
      call mpas_pool_get_config(configPool, 'config_cvmix_kpp_surface_layer_extent', config_cvmix_kpp_surface_layer_extent)
      call mpas_pool_get_config(configPool, 'config_cvmix_kpp_surface_layer_averaging', config_cvmix_kpp_surface_layer_averaging)
      call mpas_pool_get_config(configPool, 'configure_cvmix_kpp_minimum_OBL_under_sea_ice', &
configure_cvmix_kpp_minimum_OBL_under_sea_ice)
      call mpas_pool_get_config(configPool, 'config_cvmix_kpp_stop_OBL_search', config_cvmix_kpp_stop_OBL_search)
      call mpas_pool_get_config(configPool, 'config_cvmix_kpp_use_enhanced_diff', config_cvmix_kpp_use_enhanced_diff)
      call mpas_pool_get_config(configPool, 'config_cvmix_kpp_nonlocal_with_implicit_mix', &
config_cvmix_kpp_nonlocal_with_implicit_mix)
      call mpas_pool_get_config(configPool, 'config_cvmix_kpp_use_theory_wave', config_cvmix_kpp_use_theory_wave)
      call mpas_pool_get_config(configPool, 'config_cvmix_kpp_langmuir_mixing_opt', config_cvmix_kpp_langmuir_mixing_opt)
      call mpas_pool_get_config(configPool, 'config_cvmix_kpp_langmuir_entrainment_opt', config_cvmix_kpp_langmuir_entrainment_opt)

      call mpas_pool_get_config(configPool, 'config_use_gotm', config_use_gotm)
      call mpas_pool_get_config(configPool, 'config_gotm_namelist_file', config_gotm_namelist_file)
      call mpas_pool_get_config(configPool, 'config_gotm_constant_surface_roughness_length', &
config_gotm_constant_surface_roughness_length)
      call mpas_pool_get_config(configPool, 'config_gotm_constant_bottom_roughness_length', &
config_gotm_constant_bottom_roughness_length)
      call mpas_pool_get_config(configPool, 'config_gotm_constant_bottom_drag_coeff', config_gotm_constant_bottom_drag_coeff)

      call mpas_pool_get_config(configPool, 'config_use_variable_drag', config_use_variable_drag)
      call mpas_pool_get_config(configPool, 'config_use_bulk_wind_stress', config_use_bulk_wind_stress)
      call mpas_pool_get_config(configPool, 'config_use_bulk_thickness_flux', config_use_bulk_thickness_flux)
      call mpas_pool_get_config(configPool, 'config_flux_attenuation_coefficient', config_flux_attenuation_coefficient)
      call mpas_pool_get_config(configPool, 'config_flux_attenuation_coefficient_runoff', &
config_flux_attenuation_coefficient_runoff)

      call mpas_pool_get_config(configPool, 'config_use_time_varying_atmospheric_forcing', &
config_use_time_varying_atmospheric_forcing)
      call mpas_pool_get_config(configPool, 'config_time_varying_atmospheric_forcing_type', &
config_time_varying_atmospheric_forcing_type)
      call mpas_pool_get_config(configPool, 'config_time_varying_atmospheric_forcing_start_time', &
config_time_varying_atmospheric_forcing_start_time)
      call mpas_pool_get_config(configPool, 'config_time_varying_atmospheric_forcing_reference_time', &
config_time_varying_atmospheric_forcing_reference_time)
      call mpas_pool_get_config(configPool, 'config_time_varying_atmospheric_forcing_cycle_start', &
config_time_varying_atmospheric_forcing_cycle_start)
      call mpas_pool_get_config(configPool, 'config_time_varying_atmospheric_forcing_cycle_duration', &
config_time_varying_atmospheric_forcing_cycle_duration)
      call mpas_pool_get_config(configPool, 'config_time_varying_atmospheric_forcing_interval', &
config_time_varying_atmospheric_forcing_interval)
      call mpas_pool_get_config(configPool, 'config_time_varying_atmospheric_forcing_ramp', &
config_time_varying_atmospheric_forcing_ramp)
      call mpas_pool_get_config(configPool, 'config_time_varying_atmospheric_forcing_ramp_delay', &
config_time_varying_atmospheric_forcing_ramp_delay)
      call mpas_pool_get_config(configPool, 'config_use_time_varying_land_ice_forcing', config_use_time_varying_land_ice_forcing)
      call mpas_pool_get_config(configPool, 'config_time_varying_land_ice_forcing_start_time', &
config_time_varying_land_ice_forcing_start_time)
      call mpas_pool_get_config(configPool, 'config_time_varying_land_ice_forcing_reference_time', &
config_time_varying_land_ice_forcing_reference_time)
      call mpas_pool_get_config(configPool, 'config_time_varying_land_ice_forcing_cycle_start', &
config_time_varying_land_ice_forcing_cycle_start)
      call mpas_pool_get_config(configPool, 'config_time_varying_land_ice_forcing_cycle_duration', &
config_time_varying_land_ice_forcing_cycle_duration)
      call mpas_pool_get_config(configPool, 'config_time_varying_land_ice_forcing_interval', &
config_time_varying_land_ice_forcing_interval)

      call mpas_pool_get_config(configPool, 'config_ssh_grad_relax_timescale', config_ssh_grad_relax_timescale)
      call mpas_pool_get_config(configPool, 'config_remove_AIS_coupler_runoff', config_remove_AIS_coupler_runoff)

      call mpas_pool_get_config(configPool, 'config_sw_absorption_type', config_sw_absorption_type)
      call mpas_pool_get_config(configPool, 'config_jerlov_water_type', config_jerlov_water_type)
      call mpas_pool_get_config(configPool, 'config_surface_buoyancy_depth', config_surface_buoyancy_depth)
      call mpas_pool_get_config(configPool, 'config_enable_shortwave_energy_fixer', config_enable_shortwave_energy_fixer)

      call mpas_pool_get_config(configPool, 'config_use_tidal_forcing', config_use_tidal_forcing)
      call mpas_pool_get_config(configPool, 'config_use_tidal_forcing_tau', config_use_tidal_forcing_tau)
      call mpas_pool_get_config(configPool, 'config_tidal_forcing_type', config_tidal_forcing_type)
      call mpas_pool_get_config(configPool, 'config_tidal_forcing_model', config_tidal_forcing_model)
      call mpas_pool_get_config(configPool, 'config_tidal_forcing_monochromatic_amp', config_tidal_forcing_monochromatic_amp)
      call mpas_pool_get_config(configPool, 'config_tidal_forcing_monochromatic_period', config_tidal_forcing_monochromatic_period)
      call mpas_pool_get_config(configPool, 'config_tidal_forcing_monochromatic_phaseLag', &
config_tidal_forcing_monochromatic_phaseLag)
      call mpas_pool_get_config(configPool, 'config_tidal_forcing_monochromatic_baseline', &
config_tidal_forcing_monochromatic_baseline)

      call mpas_pool_get_config(configPool, 'config_use_tidal_potential_forcing', config_use_tidal_potential_forcing)
      call mpas_pool_get_config(configPool, 'config_tidal_potential_reference_time', config_tidal_potential_reference_time)
      call mpas_pool_get_config(configPool, 'config_use_tidal_potential_forcing_M2', config_use_tidal_potential_forcing_M2)
      call mpas_pool_get_config(configPool, 'config_use_tidal_potential_forcing_S2', config_use_tidal_potential_forcing_S2)
      call mpas_pool_get_config(configPool, 'config_use_tidal_potential_forcing_N2', config_use_tidal_potential_forcing_N2)
      call mpas_pool_get_config(configPool, 'config_use_tidal_potential_forcing_K2', config_use_tidal_potential_forcing_K2)
      call mpas_pool_get_config(configPool, 'config_use_tidal_potential_forcing_K1', config_use_tidal_potential_forcing_K1)
      call mpas_pool_get_config(configPool, 'config_use_tidal_potential_forcing_O1', config_use_tidal_potential_forcing_O1)
      call mpas_pool_get_config(configPool, 'config_use_tidal_potential_forcing_Q1', config_use_tidal_potential_forcing_Q1)
      call mpas_pool_get_config(configPool, 'config_use_tidal_potential_forcing_P1', config_use_tidal_potential_forcing_P1)
      call mpas_pool_get_config(configPool, 'config_tidal_potential_ramp', config_tidal_potential_ramp)
      call mpas_pool_get_config(configPool, 'config_self_attraction_and_loading_beta', config_self_attraction_and_loading_beta)

      call mpas_pool_get_config(configPool, 'config_use_vegetation_drag', config_use_vegetation_drag)
      call mpas_pool_get_config(configPool, 'config_use_vegetation_manning_equation', config_use_vegetation_manning_equation)
      call mpas_pool_get_config(configPool, 'config_vegetation_drag_coefficient', config_vegetation_drag_coefficient)

      call mpas_pool_get_config(configPool, 'config_use_frazil_ice_formation', config_use_frazil_ice_formation)
      call mpas_pool_get_config(configPool, 'config_frazil_in_open_ocean', config_frazil_in_open_ocean)
      call mpas_pool_get_config(configPool, 'config_frazil_under_land_ice', config_frazil_under_land_ice)
      call mpas_pool_get_config(configPool, 'config_frazil_heat_of_fusion', config_frazil_heat_of_fusion)
      call mpas_pool_get_config(configPool, 'config_frazil_ice_density', config_frazil_ice_density)
      call mpas_pool_get_config(configPool, 'config_frazil_fractional_thickness_limit', config_frazil_fractional_thickness_limit)
      call mpas_pool_get_config(configPool, 'config_specific_heat_sea_water', config_specific_heat_sea_water)
      call mpas_pool_get_config(configPool, 'config_frazil_maximum_depth', config_frazil_maximum_depth)
      call mpas_pool_get_config(configPool, 'config_frazil_sea_ice_reference_salinity', config_frazil_sea_ice_reference_salinity)
      call mpas_pool_get_config(configPool, 'config_frazil_land_ice_reference_salinity', config_frazil_land_ice_reference_salinity)
      call mpas_pool_get_config(configPool, 'config_frazil_maximum_freezing_temperature', &
config_frazil_maximum_freezing_temperature)
      call mpas_pool_get_config(configPool, 'config_frazil_use_surface_pressure', config_frazil_use_surface_pressure)

      call mpas_pool_get_config(configPool, 'config_land_ice_flux_mode', config_land_ice_flux_mode)
      call mpas_pool_get_config(configPool, 'config_land_ice_flux_formulation', config_land_ice_flux_formulation)
      call mpas_pool_get_config(configPool, 'config_land_ice_flux_useHollandJenkinsAdvDiff', &
config_land_ice_flux_useHollandJenkinsAdvDiff)
      call mpas_pool_get_config(configPool, 'config_land_ice_flux_attenuation_coefficient', &
config_land_ice_flux_attenuation_coefficient)
      call mpas_pool_get_config(configPool, 'config_land_ice_flux_boundaryLayerThickness', &
config_land_ice_flux_boundaryLayerThickness)
      call mpas_pool_get_config(configPool, 'config_land_ice_flux_boundaryLayerNeighborWeight', &
config_land_ice_flux_boundaryLayerNeighborWeight)
      call mpas_pool_get_config(configPool, 'config_land_ice_flux_cp_ice', config_land_ice_flux_cp_ice)
      call mpas_pool_get_config(configPool, 'config_land_ice_flux_rho_ice', config_land_ice_flux_rho_ice)
      call mpas_pool_get_config(configPool, 'config_land_ice_flux_topDragCoeff', config_land_ice_flux_topDragCoeff)
      call mpas_pool_get_config(configPool, 'config_land_ice_flux_ISOMIP_gammaT', config_land_ice_flux_ISOMIP_gammaT)
      call mpas_pool_get_config(configPool, 'config_land_ice_flux_rms_tidal_velocity', config_land_ice_flux_rms_tidal_velocity)
      call mpas_pool_get_config(configPool, 'config_land_ice_flux_jenkins_heat_transfer_coefficient', &
config_land_ice_flux_jenkins_heat_transfer_coefficient)
      call mpas_pool_get_config(configPool, 'config_land_ice_flux_jenkins_salt_transfer_coefficient', &
config_land_ice_flux_jenkins_salt_transfer_coefficient)

      call mpas_pool_get_config(configPool, 'config_vert_tracer_adv', config_vert_tracer_adv)
      call mpas_pool_get_config(configPool, 'config_vert_tracer_adv_order', config_vert_tracer_adv_order)
      call mpas_pool_get_config(configPool, 'config_horiz_tracer_adv_order', config_horiz_tracer_adv_order)
      call mpas_pool_get_config(configPool, 'config_coef_3rd_order', config_coef_3rd_order)
      call mpas_pool_get_config(configPool, 'config_monotonic', config_monotonic)

      call mpas_pool_get_config(configPool, 'config_use_implicit_bottom_drag', config_use_implicit_bottom_drag)
      call mpas_pool_get_config(configPool, 'config_implicit_bottom_drag_coeff', config_implicit_bottom_drag_coeff)
      call mpas_pool_get_config(configPool, 'config_use_implicit_bottom_roughness', config_use_implicit_bottom_roughness)
      call mpas_pool_get_config(configPool, 'config_use_implicit_bottom_drag_variable', config_use_implicit_bottom_drag_variable)
      call mpas_pool_get_config(configPool, 'config_use_implicit_bottom_drag_variable_mannings', &
config_use_implicit_bottom_drag_variable_mannings)
      call mpas_pool_get_config(configPool, 'config_use_explicit_bottom_drag', config_use_explicit_bottom_drag)
      call mpas_pool_get_config(configPool, 'config_explicit_bottom_drag_coeff', config_explicit_bottom_drag_coeff)
      call mpas_pool_get_config(configPool, 'config_use_topographic_wave_drag', config_use_topographic_wave_drag)
      call mpas_pool_get_config(configPool, 'config_topographic_wave_drag_coeff', config_topographic_wave_drag_coeff)

      call mpas_pool_get_config(configPool, 'config_use_wetting_drying', config_use_wetting_drying)
      call mpas_pool_get_config(configPool, 'config_prevent_drying', config_prevent_drying)
      call mpas_pool_get_config(configPool, 'config_drying_min_cell_height', config_drying_min_cell_height)
      call mpas_pool_get_config(configPool, 'config_zero_drying_velocity', config_zero_drying_velocity)
      call mpas_pool_get_config(configPool, 'config_verify_not_dry', config_verify_not_dry)
      call mpas_pool_get_config(configPool, 'config_thickness_flux_type', config_thickness_flux_type)
      call mpas_pool_get_config(configPool, 'config_drying_safety_height', config_drying_safety_height)

      call mpas_pool_get_config(configPool, 'config_density0', config_density0)

      call mpas_pool_get_config(configPool, 'config_pressure_gradient_type', config_pressure_gradient_type)
      call mpas_pool_get_config(configPool, 'config_common_level_weight', config_common_level_weight)
      call mpas_pool_get_config(configPool, 'config_zonal_ssh_grad', config_zonal_ssh_grad)
      call mpas_pool_get_config(configPool, 'config_meridional_ssh_grad', config_meridional_ssh_grad)

      call mpas_pool_get_config(configPool, 'config_eos_type', config_eos_type)
      call mpas_pool_get_config(configPool, 'config_open_ocean_freezing_temperature_coeff_0', &
config_open_ocean_freezing_temperature_coeff_0)
      call mpas_pool_get_config(configPool, 'config_open_ocean_freezing_temperature_coeff_S', &
config_open_ocean_freezing_temperature_coeff_S)
      call mpas_pool_get_config(configPool, 'config_open_ocean_freezing_temperature_coeff_p', &
config_open_ocean_freezing_temperature_coeff_p)
      call mpas_pool_get_config(configPool, 'config_open_ocean_freezing_temperature_coeff_pS', &
config_open_ocean_freezing_temperature_coeff_pS)
      call mpas_pool_get_config(configPool, 'config_open_ocean_freezing_temperature_coeff_mushy_az1_liq', &
config_open_ocean_freezing_temperature_coeff_mushy_az1_liq)
      call mpas_pool_get_config(configPool, 'config_land_ice_cavity_freezing_temperature_coeff_0', &
config_land_ice_cavity_freezing_temperature_coeff_0)
      call mpas_pool_get_config(configPool, 'config_land_ice_cavity_freezing_temperature_coeff_S', &
config_land_ice_cavity_freezing_temperature_coeff_S)
      call mpas_pool_get_config(configPool, 'config_land_ice_cavity_freezing_temperature_coeff_p', &
config_land_ice_cavity_freezing_temperature_coeff_p)
      call mpas_pool_get_config(configPool, 'config_land_ice_cavity_freezing_temperature_coeff_pS', &
config_land_ice_cavity_freezing_temperature_coeff_pS)

      call mpas_pool_get_config(configPool, 'config_eos_linear_alpha', config_eos_linear_alpha)
      call mpas_pool_get_config(configPool, 'config_eos_linear_beta', config_eos_linear_beta)
      call mpas_pool_get_config(configPool, 'config_eos_linear_Tref', config_eos_linear_Tref)
      call mpas_pool_get_config(configPool, 'config_eos_linear_Sref', config_eos_linear_Sref)
      call mpas_pool_get_config(configPool, 'config_eos_linear_densityref', config_eos_linear_densityref)

      call mpas_pool_get_config(configPool, 'config_eos_wright_ref_pressure', config_eos_wright_ref_pressure)

      call mpas_pool_get_config(configPool, 'config_n_ts_iter', config_n_ts_iter)
      call mpas_pool_get_config(configPool, 'config_n_bcl_iter_beg', config_n_bcl_iter_beg)
      call mpas_pool_get_config(configPool, 'config_n_bcl_iter_mid', config_n_bcl_iter_mid)
      call mpas_pool_get_config(configPool, 'config_n_bcl_iter_end', config_n_bcl_iter_end)

      call mpas_pool_get_config(configPool, 'config_btr_dt', config_btr_dt)
      call mpas_pool_get_config(configPool, 'config_n_btr_cor_iter', config_n_btr_cor_iter)
      call mpas_pool_get_config(configPool, 'config_vel_correction', config_vel_correction)
      call mpas_pool_get_config(configPool, 'config_btr_subcycle_loop_factor', config_btr_subcycle_loop_factor)
      call mpas_pool_get_config(configPool, 'config_btr_gam1_velWt1', config_btr_gam1_velWt1)
      call mpas_pool_get_config(configPool, 'config_btr_gam2_SSHWt1', config_btr_gam2_SSHWt1)
      call mpas_pool_get_config(configPool, 'config_btr_gam3_velWt2', config_btr_gam3_velWt2)
      call mpas_pool_get_config(configPool, 'config_btr_solve_SSH2', config_btr_solve_SSH2)

      call mpas_pool_get_config(configPool, 'config_btr_si_preconditioner', config_btr_si_preconditioner)
      call mpas_pool_get_config(configPool, 'config_btr_si_tolerance', config_btr_si_tolerance)
      call mpas_pool_get_config(configPool, 'config_n_btr_si_outer_iter', config_n_btr_si_outer_iter)
      call mpas_pool_get_config(configPool, 'config_btr_si_partition_match_mode', config_btr_si_partition_match_mode)

      call mpas_pool_get_config(configPool, 'config_vert_coord_movement', config_vert_coord_movement)
      call mpas_pool_get_config(configPool, 'config_ALE_thickness_proportionality', config_ALE_thickness_proportionality)
      call mpas_pool_get_config(configPool, 'config_vert_taper_weight_depth_1', config_vert_taper_weight_depth_1)
      call mpas_pool_get_config(configPool, 'config_vert_taper_weight_depth_2', config_vert_taper_weight_depth_2)
      call mpas_pool_get_config(configPool, 'config_use_min_max_thickness', config_use_min_max_thickness)
      call mpas_pool_get_config(configPool, 'config_min_thickness', config_min_thickness)
      call mpas_pool_get_config(configPool, 'config_max_thickness_factor', config_max_thickness_factor)
      call mpas_pool_get_config(configPool, 'config_dzdk_positive', config_dzdk_positive)

      call mpas_pool_get_config(configPool, 'config_use_freq_filtered_thickness', config_use_freq_filtered_thickness)
      call mpas_pool_get_config(configPool, 'config_thickness_filter_timescale', config_thickness_filter_timescale)
      call mpas_pool_get_config(configPool, 'config_use_highFreqThick_restore', config_use_highFreqThick_restore)
      call mpas_pool_get_config(configPool, 'config_highFreqThick_restore_time', config_highFreqThick_restore_time)
      call mpas_pool_get_config(configPool, 'config_use_highFreqThick_del2', config_use_highFreqThick_del2)
      call mpas_pool_get_config(configPool, 'config_highFreqThick_del2', config_highFreqThick_del2)

      call mpas_pool_get_config(configPool, 'config_check_zlevel_consistency', config_check_zlevel_consistency)
      call mpas_pool_get_config(configPool, 'config_check_ssh_consistency', config_check_ssh_consistency)
      call mpas_pool_get_config(configPool, 'config_filter_btr_mode', config_filter_btr_mode)
      call mpas_pool_get_config(configPool, 'config_prescribe_velocity', config_prescribe_velocity)
      call mpas_pool_get_config(configPool, 'config_prescribe_thickness', config_prescribe_thickness)
      call mpas_pool_get_config(configPool, 'config_include_KE_vertex', config_include_KE_vertex)
      call mpas_pool_get_config(configPool, 'config_check_tracer_monotonicity', config_check_tracer_monotonicity)
      call mpas_pool_get_config(configPool, 'config_compute_active_tracer_budgets', config_compute_active_tracer_budgets)
      call mpas_pool_get_config(configPool, 'config_disable_thick_all_tend', config_disable_thick_all_tend)
      call mpas_pool_get_config(configPool, 'config_disable_thick_hadv', config_disable_thick_hadv)
      call mpas_pool_get_config(configPool, 'config_disable_thick_vadv', config_disable_thick_vadv)
      call mpas_pool_get_config(configPool, 'config_disable_thick_sflux', config_disable_thick_sflux)
      call mpas_pool_get_config(configPool, 'config_disable_vel_all_tend', config_disable_vel_all_tend)
      call mpas_pool_get_config(configPool, 'config_disable_vel_coriolis', config_disable_vel_coriolis)
      call mpas_pool_get_config(configPool, 'config_disable_vel_pgrad', config_disable_vel_pgrad)
      call mpas_pool_get_config(configPool, 'config_disable_vel_hmix', config_disable_vel_hmix)
      call mpas_pool_get_config(configPool, 'config_disable_vel_surface_stress', config_disable_vel_surface_stress)
      call mpas_pool_get_config(configPool, 'config_disable_vel_topographic_wave_drag', config_disable_vel_topographic_wave_drag)
      call mpas_pool_get_config(configPool, 'config_disable_vel_explicit_bottom_drag', config_disable_vel_explicit_bottom_drag)
      call mpas_pool_get_config(configPool, 'config_disable_vel_vmix', config_disable_vel_vmix)
      call mpas_pool_get_config(configPool, 'config_disable_vel_vadv', config_disable_vel_vadv)
      call mpas_pool_get_config(configPool, 'config_disable_tr_all_tend', config_disable_tr_all_tend)
      call mpas_pool_get_config(configPool, 'config_disable_tr_adv', config_disable_tr_adv)
      call mpas_pool_get_config(configPool, 'config_disable_tr_hmix', config_disable_tr_hmix)
      call mpas_pool_get_config(configPool, 'config_disable_tr_vmix', config_disable_tr_vmix)
      call mpas_pool_get_config(configPool, 'config_disable_tr_sflux', config_disable_tr_sflux)
      call mpas_pool_get_config(configPool, 'config_disable_tr_nonlocalflux', config_disable_tr_nonlocalflux)
      call mpas_pool_get_config(configPool, 'config_disable_redi_k33', config_disable_redi_k33)
      call mpas_pool_get_config(configPool, 'config_read_nearest_restart', config_read_nearest_restart)

      call mpas_pool_get_config(configPool, 'config_conduct_tests', config_conduct_tests)
      call mpas_pool_get_config(configPool, 'config_test_tensors', config_test_tensors)
      call mpas_pool_get_config(configPool, 'config_tensor_test_function', config_tensor_test_function)

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)

      call mpas_pool_get_config(configPool, 'config_use_activeTracers', config_use_activeTracers)
      call mpas_pool_get_config(configPool, 'config_use_activeTracers_surface_bulk_forcing', &
config_use_activeTracers_surface_bulk_forcing)
      call mpas_pool_get_config(configPool, 'config_use_activeTracers_surface_restoring', &
config_use_activeTracers_surface_restoring)
      call mpas_pool_get_config(configPool, 'config_use_activeTracers_interior_restoring', &
config_use_activeTracers_interior_restoring)
      call mpas_pool_get_config(configPool, 'config_use_activeTracers_exponential_decay', &
config_use_activeTracers_exponential_decay)
      call mpas_pool_get_config(configPool, 'config_use_activeTracers_idealAge_forcing', config_use_activeTracers_idealAge_forcing)
      call mpas_pool_get_config(configPool, 'config_use_activeTracers_ttd_forcing', config_use_activeTracers_ttd_forcing)
      call mpas_pool_get_config(configPool, 'config_use_surface_salinity_monthly_restoring', &
config_use_surface_salinity_monthly_restoring)
      call mpas_pool_get_config(configPool, 'config_surface_salinity_monthly_restoring_compute_interval', &
config_surface_salinity_monthly_restoring_compute_interval)
      call mpas_pool_get_config(configPool, 'config_salinity_restoring_constant_piston_velocity', &
config_salinity_restoring_constant_piston_velocity)
      call mpas_pool_get_config(configPool, 'config_salinity_restoring_max_difference', config_salinity_restoring_max_difference)
      call mpas_pool_get_config(configPool, 'config_salinity_restoring_under_sea_ice', config_salinity_restoring_under_sea_ice)

      call mpas_pool_get_config(configPool, 'config_use_debugTracers', config_use_debugTracers)
      call mpas_pool_get_config(configPool, 'config_reset_debugTracers_near_surface', config_reset_debugTracers_near_surface)
      call mpas_pool_get_config(configPool, 'config_reset_debugTracers_top_nLayers', config_reset_debugTracers_top_nLayers)
      call mpas_pool_get_config(configPool, 'config_use_debugTracers_surface_bulk_forcing', &
config_use_debugTracers_surface_bulk_forcing)
      call mpas_pool_get_config(configPool, 'config_use_debugTracers_surface_restoring', config_use_debugTracers_surface_restoring)
      call mpas_pool_get_config(configPool, 'config_use_debugTracers_interior_restoring', &
config_use_debugTracers_interior_restoring)
      call mpas_pool_get_config(configPool, 'config_use_debugTracers_exponential_decay', config_use_debugTracers_exponential_decay)
      call mpas_pool_get_config(configPool, 'config_use_debugTracers_idealAge_forcing', config_use_debugTracers_idealAge_forcing)
      call mpas_pool_get_config(configPool, 'config_use_debugTracers_ttd_forcing', config_use_debugTracers_ttd_forcing)

      call mpas_pool_get_config(configPool, 'config_use_ecosysTracers', config_use_ecosysTracers)
      call mpas_pool_get_config(configPool, 'config_ecosys_atm_co2_option', config_ecosys_atm_co2_option)
      call mpas_pool_get_config(configPool, 'config_ecosys_atm_alt_co2_option', config_ecosys_atm_alt_co2_option)
      call mpas_pool_get_config(configPool, 'config_ecosys_atm_alt_co2_use_eco', config_ecosys_atm_alt_co2_use_eco)
      call mpas_pool_get_config(configPool, 'config_ecosys_atm_co2_constant_value', config_ecosys_atm_co2_constant_value)
      call mpas_pool_get_config(configPool, 'config_use_ecosysTracers_surface_bulk_forcing', &
config_use_ecosysTracers_surface_bulk_forcing)
      call mpas_pool_get_config(configPool, 'config_use_ecosysTracers_surface_restoring', &
config_use_ecosysTracers_surface_restoring)
      call mpas_pool_get_config(configPool, 'config_use_ecosysTracers_interior_restoring', &
config_use_ecosysTracers_interior_restoring)
      call mpas_pool_get_config(configPool, 'config_use_ecosysTracers_exponential_decay', &
config_use_ecosysTracers_exponential_decay)
      call mpas_pool_get_config(configPool, 'config_use_ecosysTracers_idealAge_forcing', config_use_ecosysTracers_idealAge_forcing)
      call mpas_pool_get_config(configPool, 'config_use_ecosysTracers_ttd_forcing', config_use_ecosysTracers_ttd_forcing)
      call mpas_pool_get_config(configPool, 'config_use_ecosysTracers_surface_value', config_use_ecosysTracers_surface_value)
      call mpas_pool_get_config(configPool, 'config_use_ecosysTracers_sea_ice_coupling', config_use_ecosysTracers_sea_ice_coupling)
      call mpas_pool_get_config(configPool, 'config_ecosysTracers_diagnostic_fields_level1', &
config_ecosysTracers_diagnostic_fields_level1)
      call mpas_pool_get_config(configPool, 'config_ecosysTracers_diagnostic_fields_level2', &
config_ecosysTracers_diagnostic_fields_level2)
      call mpas_pool_get_config(configPool, 'config_ecosysTracers_diagnostic_fields_level3', &
config_ecosysTracers_diagnostic_fields_level3)
      call mpas_pool_get_config(configPool, 'config_ecosysTracers_diagnostic_fields_level4', &
config_ecosysTracers_diagnostic_fields_level4)
      call mpas_pool_get_config(configPool, 'config_ecosysTracers_diagnostic_fields_level5', &
config_ecosysTracers_diagnostic_fields_level5)

      call mpas_pool_get_config(configPool, 'config_use_DMSTracers', config_use_DMSTracers)
      call mpas_pool_get_config(configPool, 'config_use_DMSTracers_surface_bulk_forcing', &
config_use_DMSTracers_surface_bulk_forcing)
      call mpas_pool_get_config(configPool, 'config_use_DMSTracers_surface_restoring', config_use_DMSTracers_surface_restoring)
      call mpas_pool_get_config(configPool, 'config_use_DMSTracers_interior_restoring', config_use_DMSTracers_interior_restoring)
      call mpas_pool_get_config(configPool, 'config_use_DMSTracers_exponential_decay', config_use_DMSTracers_exponential_decay)
      call mpas_pool_get_config(configPool, 'config_use_DMSTracers_idealAge_forcing', config_use_DMSTracers_idealAge_forcing)
      call mpas_pool_get_config(configPool, 'config_use_DMSTracers_ttd_forcing', config_use_DMSTracers_ttd_forcing)
      call mpas_pool_get_config(configPool, 'config_use_DMSTracers_surface_value', config_use_DMSTracers_surface_value)
      call mpas_pool_get_config(configPool, 'config_use_DMSTracers_sea_ice_coupling', config_use_DMSTracers_sea_ice_coupling)

      call mpas_pool_get_config(configPool, 'config_use_MacroMoleculesTracers', config_use_MacroMoleculesTracers)
      call mpas_pool_get_config(configPool, 'config_use_MacroMoleculesTracers_surface_bulk_forcing', &
config_use_MacroMoleculesTracers_surface_bulk_forcing)
      call mpas_pool_get_config(configPool, 'config_use_MacroMoleculesTracers_surface_restoring', &
config_use_MacroMoleculesTracers_surface_restoring)
      call mpas_pool_get_config(configPool, 'config_use_MacroMoleculesTracers_interior_restoring', &
config_use_MacroMoleculesTracers_interior_restoring)
      call mpas_pool_get_config(configPool, 'config_use_MacroMoleculesTracers_exponential_decay', &
config_use_MacroMoleculesTracers_exponential_decay)
      call mpas_pool_get_config(configPool, 'config_use_MacroMoleculesTracers_idealAge_forcing', &
config_use_MacroMoleculesTracers_idealAge_forcing)
      call mpas_pool_get_config(configPool, 'config_use_MacroMoleculesTracers_ttd_forcing', &
config_use_MacroMoleculesTracers_ttd_forcing)
      call mpas_pool_get_config(configPool, 'config_use_MacroMoleculesTracers_surface_value', &
config_use_MacroMoleculesTracers_surface_value)
      call mpas_pool_get_config(configPool, 'config_use_MacroMoleculesTracers_sea_ice_coupling', &
config_use_MacroMoleculesTracers_sea_ice_coupling)

      call mpas_pool_get_config(configPool, 'config_AM_globalStats_enable', config_AM_globalStats_enable)
      call mpas_pool_get_config(configPool, 'config_AM_globalStats_compute_interval', config_AM_globalStats_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_globalStats_compute_on_startup', config_AM_globalStats_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_globalStats_write_on_startup', config_AM_globalStats_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_globalStats_text_file', config_AM_globalStats_text_file)
      call mpas_pool_get_config(configPool, 'config_AM_globalStats_directory', config_AM_globalStats_directory)
      call mpas_pool_get_config(configPool, 'config_AM_globalStats_output_stream', config_AM_globalStats_output_stream)

      call mpas_pool_get_config(configPool, 'config_AM_surfaceAreaWeightedAverages_enable', &
config_AM_surfaceAreaWeightedAverages_enable)
      call mpas_pool_get_config(configPool, 'config_AM_surfaceAreaWeightedAverages_compute_on_startup', &
config_AM_surfaceAreaWeightedAverages_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_surfaceAreaWeightedAverages_write_on_startup', &
config_AM_surfaceAreaWeightedAverages_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_surfaceAreaWeightedAverages_compute_interval', &
config_AM_surfaceAreaWeightedAverages_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_surfaceAreaWeightedAverages_output_stream', &
config_AM_surfaceAreaWeightedAverages_output_stream)

      call mpas_pool_get_config(configPool, 'config_AM_waterMassCensus_enable', config_AM_waterMassCensus_enable)
      call mpas_pool_get_config(configPool, 'config_AM_waterMassCensus_compute_interval', &
config_AM_waterMassCensus_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_waterMassCensus_output_stream', config_AM_waterMassCensus_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_waterMassCensus_compute_on_startup', &
config_AM_waterMassCensus_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_waterMassCensus_write_on_startup', &
config_AM_waterMassCensus_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_waterMassCensus_minTemperature', config_AM_waterMassCensus_minTemperature)
      call mpas_pool_get_config(configPool, 'config_AM_waterMassCensus_maxTemperature', config_AM_waterMassCensus_maxTemperature)
      call mpas_pool_get_config(configPool, 'config_AM_waterMassCensus_minSalinity', config_AM_waterMassCensus_minSalinity)
      call mpas_pool_get_config(configPool, 'config_AM_waterMassCensus_maxSalinity', config_AM_waterMassCensus_maxSalinity)
      call mpas_pool_get_config(configPool, 'config_AM_waterMassCensus_compute_predefined_regions', &
config_AM_waterMassCensus_compute_predefined_regions)
      call mpas_pool_get_config(configPool, 'config_AM_waterMassCensus_region_group', config_AM_waterMassCensus_region_group)

      call mpas_pool_get_config(configPool, 'config_AM_layerVolumeWeightedAverage_enable', &
config_AM_layerVolumeWeightedAverage_enable)
      call mpas_pool_get_config(configPool, 'config_AM_layerVolumeWeightedAverage_compute_interval', &
config_AM_layerVolumeWeightedAverage_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_layerVolumeWeightedAverage_compute_on_startup', &
config_AM_layerVolumeWeightedAverage_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_layerVolumeWeightedAverage_write_on_startup', &
config_AM_layerVolumeWeightedAverage_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_layerVolumeWeightedAverage_output_stream', &
config_AM_layerVolumeWeightedAverage_output_stream)

      call mpas_pool_get_config(configPool, 'config_AM_zonalMean_enable', config_AM_zonalMean_enable)
      call mpas_pool_get_config(configPool, 'config_AM_zonalMean_compute_on_startup', config_AM_zonalMean_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_zonalMean_write_on_startup', config_AM_zonalMean_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_zonalMean_compute_interval', config_AM_zonalMean_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_zonalMean_output_stream', config_AM_zonalMean_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_zonalMean_num_bins', config_AM_zonalMean_num_bins)
      call mpas_pool_get_config(configPool, 'config_AM_zonalMean_min_bin', config_AM_zonalMean_min_bin)
      call mpas_pool_get_config(configPool, 'config_AM_zonalMean_max_bin', config_AM_zonalMean_max_bin)

      call mpas_pool_get_config(configPool, 'config_AM_okuboWeiss_enable', config_AM_okuboWeiss_enable)
      call mpas_pool_get_config(configPool, 'config_AM_okuboWeiss_compute_on_startup', config_AM_okuboWeiss_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_okuboWeiss_write_on_startup', config_AM_okuboWeiss_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_okuboWeiss_compute_interval', config_AM_okuboWeiss_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_okuboWeiss_output_stream', config_AM_okuboWeiss_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_okuboWeiss_directory', config_AM_okuboWeiss_directory)
      call mpas_pool_get_config(configPool, 'config_AM_okuboWeiss_threshold_value', config_AM_okuboWeiss_threshold_value)
      call mpas_pool_get_config(configPool, 'config_AM_okuboWeiss_normalization', config_AM_okuboWeiss_normalization)
      call mpas_pool_get_config(configPool, 'config_AM_okuboWeiss_lambda2_normalization', &
config_AM_okuboWeiss_lambda2_normalization)
      call mpas_pool_get_config(configPool, 'config_AM_okuboWeiss_use_lat_lon_coords', config_AM_okuboWeiss_use_lat_lon_coords)
      call mpas_pool_get_config(configPool, 'config_AM_okuboWeiss_compute_eddy_census', config_AM_okuboWeiss_compute_eddy_census)
      call mpas_pool_get_config(configPool, 'config_AM_okuboWeiss_eddy_min_cells', config_AM_okuboWeiss_eddy_min_cells)

      call mpas_pool_get_config(configPool, 'config_AM_meridionalHeatTransport_enable', config_AM_meridionalHeatTransport_enable)
      call mpas_pool_get_config(configPool, 'config_AM_meridionalHeatTransport_compute_interval', &
config_AM_meridionalHeatTransport_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_meridionalHeatTransport_compute_on_startup', &
config_AM_meridionalHeatTransport_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_meridionalHeatTransport_write_on_startup', &
config_AM_meridionalHeatTransport_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_meridionalHeatTransport_output_stream', &
config_AM_meridionalHeatTransport_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_meridionalHeatTransport_num_bins', &
config_AM_meridionalHeatTransport_num_bins)
      call mpas_pool_get_config(configPool, 'config_AM_meridionalHeatTransport_min_bin', config_AM_meridionalHeatTransport_min_bin)
      call mpas_pool_get_config(configPool, 'config_AM_meridionalHeatTransport_max_bin', config_AM_meridionalHeatTransport_max_bin)
      call mpas_pool_get_config(configPool, 'config_AM_meridionalHeatTransport_region_group', &
config_AM_meridionalHeatTransport_region_group)

      call mpas_pool_get_config(configPool, 'config_AM_testComputeInterval_enable', config_AM_testComputeInterval_enable)
      call mpas_pool_get_config(configPool, 'config_AM_testComputeInterval_compute_interval', &
config_AM_testComputeInterval_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_testComputeInterval_compute_on_startup', &
config_AM_testComputeInterval_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_testComputeInterval_write_on_startup', &
config_AM_testComputeInterval_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_testComputeInterval_output_stream', &
config_AM_testComputeInterval_output_stream)

      call mpas_pool_get_config(configPool, 'config_AM_highFrequencyOutput_enable', config_AM_highFrequencyOutput_enable)
      call mpas_pool_get_config(configPool, 'config_AM_highFrequencyOutput_compute_interval', &
config_AM_highFrequencyOutput_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_highFrequencyOutput_output_stream', &
config_AM_highFrequencyOutput_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_highFrequencyOutput_compute_on_startup', &
config_AM_highFrequencyOutput_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_highFrequencyOutput_write_on_startup', &
config_AM_highFrequencyOutput_write_on_startup)

      call mpas_pool_get_config(configPool, 'config_AM_timeFilters_enable', config_AM_timeFilters_enable)
      call mpas_pool_get_config(configPool, 'config_AM_timeFilters_compute_interval', config_AM_timeFilters_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_timeFilters_output_stream', config_AM_timeFilters_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_timeFilters_restart_stream', config_AM_timeFilters_restart_stream)
      call mpas_pool_get_config(configPool, 'config_AM_timeFilters_compute_on_startup', config_AM_timeFilters_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_timeFilters_write_on_startup', config_AM_timeFilters_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_timeFilters_initialize_filters', config_AM_timeFilters_initialize_filters)
      call mpas_pool_get_config(configPool, 'config_AM_timeFilters_tau', config_AM_timeFilters_tau)
      call mpas_pool_get_config(configPool, 'config_AM_timeFilters_compute_cell_centered_values', &
config_AM_timeFilters_compute_cell_centered_values)

      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_enable', config_AM_lagrPartTrack_enable)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_compute_interval', config_AM_lagrPartTrack_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_compute_on_startup', &
config_AM_lagrPartTrack_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_output_stream', config_AM_lagrPartTrack_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_restart_stream', config_AM_lagrPartTrack_restart_stream)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_input_stream', config_AM_lagrPartTrack_input_stream)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_write_on_startup', config_AM_lagrPartTrack_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_filter_number', config_AM_lagrPartTrack_filter_number)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_timeIntegration', config_AM_lagrPartTrack_timeIntegration)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_reset_criteria', config_AM_lagrPartTrack_reset_criteria)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_reset_global_timestamp', &
config_AM_lagrPartTrack_reset_global_timestamp)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_region_stream', config_AM_lagrPartTrack_region_stream)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_reset_if_outside_region', &
config_AM_lagrPartTrack_reset_if_outside_region)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_reset_if_inside_region', &
config_AM_lagrPartTrack_reset_if_inside_region)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_sample_horizontal_interp', &
config_AM_lagrPartTrack_sample_horizontal_interp)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_sample_temperature', &
config_AM_lagrPartTrack_sample_temperature)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_sample_salinity', config_AM_lagrPartTrack_sample_salinity)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_sample_DIC', config_AM_lagrPartTrack_sample_DIC)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_sample_ALK', config_AM_lagrPartTrack_sample_ALK)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_sample_PO4', config_AM_lagrPartTrack_sample_PO4)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_sample_NO3', config_AM_lagrPartTrack_sample_NO3)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_sample_SiO3', config_AM_lagrPartTrack_sample_SiO3)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_sample_NH4', config_AM_lagrPartTrack_sample_NH4)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_sample_Fe', config_AM_lagrPartTrack_sample_Fe)
      call mpas_pool_get_config(configPool, 'config_AM_lagrPartTrack_sample_O2', config_AM_lagrPartTrack_sample_O2)

      call mpas_pool_get_config(configPool, 'config_AM_eliassenPalm_enable', config_AM_eliassenPalm_enable)
      call mpas_pool_get_config(configPool, 'config_AM_eliassenPalm_compute_interval', config_AM_eliassenPalm_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_eliassenPalm_output_stream', config_AM_eliassenPalm_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_eliassenPalm_restart_stream', config_AM_eliassenPalm_restart_stream)
      call mpas_pool_get_config(configPool, 'config_AM_eliassenPalm_compute_on_startup', config_AM_eliassenPalm_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_eliassenPalm_write_on_startup', config_AM_eliassenPalm_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_eliassenPalm_debug', config_AM_eliassenPalm_debug)
      call mpas_pool_get_config(configPool, 'config_AM_eliassenPalm_nBuoyancyLayers', config_AM_eliassenPalm_nBuoyancyLayers)
      call mpas_pool_get_config(configPool, 'config_AM_eliassenPalm_rhomin_buoycoor', config_AM_eliassenPalm_rhomin_buoycoor)
      call mpas_pool_get_config(configPool, 'config_AM_eliassenPalm_rhomax_buoycoor', config_AM_eliassenPalm_rhomax_buoycoor)

      call mpas_pool_get_config(configPool, 'config_AM_mixedLayerDepths_enable', config_AM_mixedLayerDepths_enable)
      call mpas_pool_get_config(configPool, 'config_AM_mixedLayerDepths_compute_interval', &
config_AM_mixedLayerDepths_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_mixedLayerDepths_output_stream', config_AM_mixedLayerDepths_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_mixedLayerDepths_write_on_startup', &
config_AM_mixedLayerDepths_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_mixedLayerDepths_compute_on_startup', &
config_AM_mixedLayerDepths_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_mixedLayerDepths_Tthreshold', config_AM_mixedLayerDepths_Tthreshold)
      call mpas_pool_get_config(configPool, 'config_AM_mixedLayerDepths_Dthreshold', config_AM_mixedLayerDepths_Dthreshold)
      call mpas_pool_get_config(configPool, 'config_AM_mixedLayerDepths_crit_temp_threshold', &
config_AM_mixedLayerDepths_crit_temp_threshold)
      call mpas_pool_get_config(configPool, 'config_AM_mixedLayerDepths_crit_dens_threshold', &
config_AM_mixedLayerDepths_crit_dens_threshold)
      call mpas_pool_get_config(configPool, 'config_AM_mixedLayerDepths_reference_pressure', &
config_AM_mixedLayerDepths_reference_pressure)
      call mpas_pool_get_config(configPool, 'config_AM_mixedLayerDepths_Tgradient', config_AM_mixedLayerDepths_Tgradient)
      call mpas_pool_get_config(configPool, 'config_AM_mixedLayerDepths_Dgradient', config_AM_mixedLayerDepths_Dgradient)
      call mpas_pool_get_config(configPool, 'config_AM_mixedLayerDepths_temp_gradient_threshold', &
config_AM_mixedLayerDepths_temp_gradient_threshold)
      call mpas_pool_get_config(configPool, 'config_AM_mixedLayerDepths_den_gradient_threshold', &
config_AM_mixedLayerDepths_den_gradient_threshold)
      call mpas_pool_get_config(configPool, 'config_AM_mixedLayerDepths_interp_method', config_AM_mixedLayerDepths_interp_method)

      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsDaily_enable', config_AM_regionalStatsDaily_enable)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsDaily_compute_on_startup', &
config_AM_regionalStatsDaily_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsDaily_write_on_startup', &
config_AM_regionalStatsDaily_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsDaily_compute_interval', &
config_AM_regionalStatsDaily_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsDaily_output_stream', &
config_AM_regionalStatsDaily_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsDaily_restart_stream', &
config_AM_regionalStatsDaily_restart_stream)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsDaily_input_stream', config_AM_regionalStatsDaily_input_stream)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsDaily_operation', config_AM_regionalStatsDaily_operation)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsDaily_region_type', config_AM_regionalStatsDaily_region_type)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsDaily_region_group', config_AM_regionalStatsDaily_region_group)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsDaily_1d_weighting_function', &
config_AM_regionalStatsDaily_1d_weighting_function)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsDaily_2d_weighting_function', &
config_AM_regionalStatsDaily_2d_weighting_function)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsDaily_1d_weighting_field', &
config_AM_regionalStatsDaily_1d_weighting_field)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsDaily_2d_weighting_field', &
config_AM_regionalStatsDaily_2d_weighting_field)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsDaily_vertical_mask', &
config_AM_regionalStatsDaily_vertical_mask)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsDaily_vertical_dimension', &
config_AM_regionalStatsDaily_vertical_dimension)

      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsWeekly_enable', config_AM_regionalStatsWeekly_enable)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsWeekly_compute_on_startup', &
config_AM_regionalStatsWeekly_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsWeekly_write_on_startup', &
config_AM_regionalStatsWeekly_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsWeekly_compute_interval', &
config_AM_regionalStatsWeekly_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsWeekly_output_stream', &
config_AM_regionalStatsWeekly_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsWeekly_restart_stream', &
config_AM_regionalStatsWeekly_restart_stream)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsWeekly_input_stream', &
config_AM_regionalStatsWeekly_input_stream)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsWeekly_operation', config_AM_regionalStatsWeekly_operation)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsWeekly_region_type', config_AM_regionalStatsWeekly_region_type)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsWeekly_region_group', &
config_AM_regionalStatsWeekly_region_group)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsWeekly_1d_weighting_function', &
config_AM_regionalStatsWeekly_1d_weighting_function)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsWeekly_2d_weighting_function', &
config_AM_regionalStatsWeekly_2d_weighting_function)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsWeekly_1d_weighting_field', &
config_AM_regionalStatsWeekly_1d_weighting_field)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsWeekly_2d_weighting_field', &
config_AM_regionalStatsWeekly_2d_weighting_field)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsWeekly_vertical_mask', &
config_AM_regionalStatsWeekly_vertical_mask)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsWeekly_vertical_dimension', &
config_AM_regionalStatsWeekly_vertical_dimension)

      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsMonthly_enable', config_AM_regionalStatsMonthly_enable)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsMonthly_compute_on_startup', &
config_AM_regionalStatsMonthly_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsMonthly_write_on_startup', &
config_AM_regionalStatsMonthly_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsMonthly_compute_interval', &
config_AM_regionalStatsMonthly_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsMonthly_output_stream', &
config_AM_regionalStatsMonthly_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsMonthly_restart_stream', &
config_AM_regionalStatsMonthly_restart_stream)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsMonthly_input_stream', &
config_AM_regionalStatsMonthly_input_stream)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsMonthly_operation', config_AM_regionalStatsMonthly_operation)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsMonthly_region_type', &
config_AM_regionalStatsMonthly_region_type)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsMonthly_region_group', &
config_AM_regionalStatsMonthly_region_group)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsMonthly_1d_weighting_function', &
config_AM_regionalStatsMonthly_1d_weighting_function)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsMonthly_2d_weighting_function', &
config_AM_regionalStatsMonthly_2d_weighting_function)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsMonthly_1d_weighting_field', &
config_AM_regionalStatsMonthly_1d_weighting_field)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsMonthly_2d_weighting_field', &
config_AM_regionalStatsMonthly_2d_weighting_field)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsMonthly_vertical_mask', &
config_AM_regionalStatsMonthly_vertical_mask)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsMonthly_vertical_dimension', &
config_AM_regionalStatsMonthly_vertical_dimension)

      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsCustom_enable', config_AM_regionalStatsCustom_enable)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsCustom_compute_on_startup', &
config_AM_regionalStatsCustom_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsCustom_write_on_startup', &
config_AM_regionalStatsCustom_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsCustom_compute_interval', &
config_AM_regionalStatsCustom_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsCustom_output_stream', &
config_AM_regionalStatsCustom_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsCustom_restart_stream', &
config_AM_regionalStatsCustom_restart_stream)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsCustom_input_stream', &
config_AM_regionalStatsCustom_input_stream)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsCustom_operation', config_AM_regionalStatsCustom_operation)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsCustom_region_type', config_AM_regionalStatsCustom_region_type)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsCustom_region_group', &
config_AM_regionalStatsCustom_region_group)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsCustom_1d_weighting_function', &
config_AM_regionalStatsCustom_1d_weighting_function)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsCustom_2d_weighting_function', &
config_AM_regionalStatsCustom_2d_weighting_function)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsCustom_1d_weighting_field', &
config_AM_regionalStatsCustom_1d_weighting_field)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsCustom_2d_weighting_field', &
config_AM_regionalStatsCustom_2d_weighting_field)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsCustom_vertical_mask', &
config_AM_regionalStatsCustom_vertical_mask)
      call mpas_pool_get_config(configPool, 'config_AM_regionalStatsCustom_vertical_dimension', &
config_AM_regionalStatsCustom_vertical_dimension)

      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsDaily_enable', config_AM_timeSeriesStatsDaily_enable)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsDaily_compute_on_startup', &
config_AM_timeSeriesStatsDaily_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsDaily_write_on_startup', &
config_AM_timeSeriesStatsDaily_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsDaily_compute_interval', &
config_AM_timeSeriesStatsDaily_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsDaily_output_stream', &
config_AM_timeSeriesStatsDaily_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsDaily_restart_stream', &
config_AM_timeSeriesStatsDaily_restart_stream)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsDaily_operation', config_AM_timeSeriesStatsDaily_operation)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsDaily_reference_times', &
config_AM_timeSeriesStatsDaily_reference_times)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsDaily_duration_intervals', &
config_AM_timeSeriesStatsDaily_duration_intervals)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsDaily_repeat_intervals', &
config_AM_timeSeriesStatsDaily_repeat_intervals)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsDaily_reset_intervals', &
config_AM_timeSeriesStatsDaily_reset_intervals)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsDaily_backward_output_offset', &
config_AM_timeSeriesStatsDaily_backward_output_offset)

      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthly_enable', config_AM_timeSeriesStatsMonthly_enable)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthly_compute_on_startup', &
config_AM_timeSeriesStatsMonthly_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthly_write_on_startup', &
config_AM_timeSeriesStatsMonthly_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthly_compute_interval', &
config_AM_timeSeriesStatsMonthly_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthly_output_stream', &
config_AM_timeSeriesStatsMonthly_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthly_restart_stream', &
config_AM_timeSeriesStatsMonthly_restart_stream)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthly_operation', &
config_AM_timeSeriesStatsMonthly_operation)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthly_reference_times', &
config_AM_timeSeriesStatsMonthly_reference_times)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthly_duration_intervals', &
config_AM_timeSeriesStatsMonthly_duration_intervals)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthly_repeat_intervals', &
config_AM_timeSeriesStatsMonthly_repeat_intervals)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthly_reset_intervals', &
config_AM_timeSeriesStatsMonthly_reset_intervals)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthly_backward_output_offset', &
config_AM_timeSeriesStatsMonthly_backward_output_offset)

      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsClimatology_enable', &
config_AM_timeSeriesStatsClimatology_enable)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsClimatology_compute_on_startup', &
config_AM_timeSeriesStatsClimatology_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsClimatology_write_on_startup', &
config_AM_timeSeriesStatsClimatology_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsClimatology_compute_interval', &
config_AM_timeSeriesStatsClimatology_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsClimatology_output_stream', &
config_AM_timeSeriesStatsClimatology_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsClimatology_restart_stream', &
config_AM_timeSeriesStatsClimatology_restart_stream)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsClimatology_operation', &
config_AM_timeSeriesStatsClimatology_operation)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsClimatology_reference_times', &
config_AM_timeSeriesStatsClimatology_reference_times)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsClimatology_duration_intervals', &
config_AM_timeSeriesStatsClimatology_duration_intervals)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsClimatology_repeat_intervals', &
config_AM_timeSeriesStatsClimatology_repeat_intervals)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsClimatology_reset_intervals', &
config_AM_timeSeriesStatsClimatology_reset_intervals)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsClimatology_backward_output_offset', &
config_AM_timeSeriesStatsClimatology_backward_output_offset)

      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_enable', &
config_AM_timeSeriesStatsMonthlyMax_enable)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_compute_on_startup', &
config_AM_timeSeriesStatsMonthlyMax_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_write_on_startup', &
config_AM_timeSeriesStatsMonthlyMax_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_compute_interval', &
config_AM_timeSeriesStatsMonthlyMax_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_output_stream', &
config_AM_timeSeriesStatsMonthlyMax_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_restart_stream', &
config_AM_timeSeriesStatsMonthlyMax_restart_stream)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_operation', &
config_AM_timeSeriesStatsMonthlyMax_operation)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_reference_times', &
config_AM_timeSeriesStatsMonthlyMax_reference_times)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_duration_intervals', &
config_AM_timeSeriesStatsMonthlyMax_duration_intervals)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_repeat_intervals', &
config_AM_timeSeriesStatsMonthlyMax_repeat_intervals)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_reset_intervals', &
config_AM_timeSeriesStatsMonthlyMax_reset_intervals)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMax_backward_output_offset', &
config_AM_timeSeriesStatsMonthlyMax_backward_output_offset)

      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_enable', &
config_AM_timeSeriesStatsMonthlyMin_enable)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_compute_on_startup', &
config_AM_timeSeriesStatsMonthlyMin_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_write_on_startup', &
config_AM_timeSeriesStatsMonthlyMin_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_compute_interval', &
config_AM_timeSeriesStatsMonthlyMin_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_output_stream', &
config_AM_timeSeriesStatsMonthlyMin_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_restart_stream', &
config_AM_timeSeriesStatsMonthlyMin_restart_stream)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_operation', &
config_AM_timeSeriesStatsMonthlyMin_operation)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_reference_times', &
config_AM_timeSeriesStatsMonthlyMin_reference_times)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_duration_intervals', &
config_AM_timeSeriesStatsMonthlyMin_duration_intervals)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_repeat_intervals', &
config_AM_timeSeriesStatsMonthlyMin_repeat_intervals)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_reset_intervals', &
config_AM_timeSeriesStatsMonthlyMin_reset_intervals)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsMonthlyMin_backward_output_offset', &
config_AM_timeSeriesStatsMonthlyMin_backward_output_offset)

      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsCustom_enable', config_AM_timeSeriesStatsCustom_enable)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsCustom_compute_on_startup', &
config_AM_timeSeriesStatsCustom_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsCustom_write_on_startup', &
config_AM_timeSeriesStatsCustom_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsCustom_compute_interval', &
config_AM_timeSeriesStatsCustom_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsCustom_output_stream', &
config_AM_timeSeriesStatsCustom_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsCustom_restart_stream', &
config_AM_timeSeriesStatsCustom_restart_stream)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsCustom_operation', config_AM_timeSeriesStatsCustom_operation)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsCustom_reference_times', &
config_AM_timeSeriesStatsCustom_reference_times)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsCustom_duration_intervals', &
config_AM_timeSeriesStatsCustom_duration_intervals)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsCustom_repeat_intervals', &
config_AM_timeSeriesStatsCustom_repeat_intervals)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsCustom_reset_intervals', &
config_AM_timeSeriesStatsCustom_reset_intervals)
      call mpas_pool_get_config(configPool, 'config_AM_timeSeriesStatsCustom_backward_output_offset', &
config_AM_timeSeriesStatsCustom_backward_output_offset)

      call mpas_pool_get_config(configPool, 'config_AM_pointwiseStats_enable', config_AM_pointwiseStats_enable)
      call mpas_pool_get_config(configPool, 'config_AM_pointwiseStats_compute_interval', config_AM_pointwiseStats_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_pointwiseStats_output_stream', config_AM_pointwiseStats_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_pointwiseStats_compute_on_startup', &
config_AM_pointwiseStats_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_pointwiseStats_write_on_startup', config_AM_pointwiseStats_write_on_startup)

      call mpas_pool_get_config(configPool, 'config_AM_debugDiagnostics_enable', config_AM_debugDiagnostics_enable)
      call mpas_pool_get_config(configPool, 'config_AM_debugDiagnostics_compute_interval', &
config_AM_debugDiagnostics_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_debugDiagnostics_output_stream', config_AM_debugDiagnostics_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_debugDiagnostics_compute_on_startup', &
config_AM_debugDiagnostics_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_debugDiagnostics_write_on_startup', &
config_AM_debugDiagnostics_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_debugDiagnostics_check_state', config_AM_debugDiagnostics_check_state)

      call mpas_pool_get_config(configPool, 'config_AM_rpnCalculator_enable', config_AM_rpnCalculator_enable)
      call mpas_pool_get_config(configPool, 'config_AM_rpnCalculator_compute_on_startup', &
config_AM_rpnCalculator_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_rpnCalculator_write_on_startup', config_AM_rpnCalculator_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_rpnCalculator_compute_interval', config_AM_rpnCalculator_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_rpnCalculator_output_stream', config_AM_rpnCalculator_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_rpnCalculator_variable_a', config_AM_rpnCalculator_variable_a)
      call mpas_pool_get_config(configPool, 'config_AM_rpnCalculator_variable_b', config_AM_rpnCalculator_variable_b)
      call mpas_pool_get_config(configPool, 'config_AM_rpnCalculator_variable_c', config_AM_rpnCalculator_variable_c)
      call mpas_pool_get_config(configPool, 'config_AM_rpnCalculator_variable_d', config_AM_rpnCalculator_variable_d)
      call mpas_pool_get_config(configPool, 'config_AM_rpnCalculator_variable_e', config_AM_rpnCalculator_variable_e)
      call mpas_pool_get_config(configPool, 'config_AM_rpnCalculator_variable_f', config_AM_rpnCalculator_variable_f)
      call mpas_pool_get_config(configPool, 'config_AM_rpnCalculator_variable_g', config_AM_rpnCalculator_variable_g)
      call mpas_pool_get_config(configPool, 'config_AM_rpnCalculator_variable_h', config_AM_rpnCalculator_variable_h)
      call mpas_pool_get_config(configPool, 'config_AM_rpnCalculator_expression_1', config_AM_rpnCalculator_expression_1)
      call mpas_pool_get_config(configPool, 'config_AM_rpnCalculator_expression_2', config_AM_rpnCalculator_expression_2)
      call mpas_pool_get_config(configPool, 'config_AM_rpnCalculator_expression_3', config_AM_rpnCalculator_expression_3)
      call mpas_pool_get_config(configPool, 'config_AM_rpnCalculator_expression_4', config_AM_rpnCalculator_expression_4)
      call mpas_pool_get_config(configPool, 'config_AM_rpnCalculator_output_name_1', config_AM_rpnCalculator_output_name_1)
      call mpas_pool_get_config(configPool, 'config_AM_rpnCalculator_output_name_2', config_AM_rpnCalculator_output_name_2)
      call mpas_pool_get_config(configPool, 'config_AM_rpnCalculator_output_name_3', config_AM_rpnCalculator_output_name_3)
      call mpas_pool_get_config(configPool, 'config_AM_rpnCalculator_output_name_4', config_AM_rpnCalculator_output_name_4)

      call mpas_pool_get_config(configPool, 'config_AM_transectTransport_enable', config_AM_transectTransport_enable)
      call mpas_pool_get_config(configPool, 'config_AM_transectTransport_compute_interval', &
config_AM_transectTransport_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_transectTransport_output_stream', config_AM_transectTransport_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_transectTransport_compute_on_startup', &
config_AM_transectTransport_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_transectTransport_write_on_startup', &
config_AM_transectTransport_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_transectTransport_transect_group', &
config_AM_transectTransport_transect_group)

      call mpas_pool_get_config(configPool, 'config_AM_eddyProductVariables_enable', config_AM_eddyProductVariables_enable)
      call mpas_pool_get_config(configPool, 'config_AM_eddyProductVariables_compute_interval', &
config_AM_eddyProductVariables_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_eddyProductVariables_output_stream', &
config_AM_eddyProductVariables_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_eddyProductVariables_compute_on_startup', &
config_AM_eddyProductVariables_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_eddyProductVariables_write_on_startup', &
config_AM_eddyProductVariables_write_on_startup)

      call mpas_pool_get_config(configPool, 'config_AM_mocStreamfunction_enable', config_AM_mocStreamfunction_enable)
      call mpas_pool_get_config(configPool, 'config_AM_mocStreamfunction_compute_interval', &
config_AM_mocStreamfunction_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_mocStreamfunction_output_stream', config_AM_mocStreamfunction_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_mocStreamfunction_compute_on_startup', &
config_AM_mocStreamfunction_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_mocStreamfunction_write_on_startup', &
config_AM_mocStreamfunction_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_mocStreamfunction_min_bin', config_AM_mocStreamfunction_min_bin)
      call mpas_pool_get_config(configPool, 'config_AM_mocStreamfunction_max_bin', config_AM_mocStreamfunction_max_bin)
      call mpas_pool_get_config(configPool, 'config_AM_mocStreamfunction_num_bins', config_AM_mocStreamfunction_num_bins)
      call mpas_pool_get_config(configPool, 'config_AM_mocStreamfunction_region_group', config_AM_mocStreamfunction_region_group)
      call mpas_pool_get_config(configPool, 'config_AM_mocStreamfunction_transect_group', &
config_AM_mocStreamfunction_transect_group)

      call mpas_pool_get_config(configPool, 'config_AM_oceanHeatContent_enable', config_AM_oceanHeatContent_enable)
      call mpas_pool_get_config(configPool, 'config_AM_oceanHeatContent_compute_interval', &
config_AM_oceanHeatContent_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_oceanHeatContent_output_stream', config_AM_oceanHeatContent_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_oceanHeatContent_compute_on_startup', &
config_AM_oceanHeatContent_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_oceanHeatContent_write_on_startup', &
config_AM_oceanHeatContent_write_on_startup)

      call mpas_pool_get_config(configPool, 'config_AM_mixedLayerHeatBudget_enable', config_AM_mixedLayerHeatBudget_enable)
      call mpas_pool_get_config(configPool, 'config_AM_mixedLayerHeatBudget_compute_interval', &
config_AM_mixedLayerHeatBudget_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_mixedLayerHeatBudget_output_stream', &
config_AM_mixedLayerHeatBudget_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_mixedLayerHeatBudget_compute_on_startup', &
config_AM_mixedLayerHeatBudget_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_mixedLayerHeatBudget_write_on_startup', &
config_AM_mixedLayerHeatBudget_write_on_startup)

      call mpas_pool_get_config(configPool, 'config_AM_sedimentFluxIndex_enable', config_AM_sedimentFluxIndex_enable)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentFluxIndex_compute_on_startup', &
config_AM_sedimentFluxIndex_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentFluxIndex_write_on_startup', &
config_AM_sedimentFluxIndex_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentFluxIndex_compute_interval', &
config_AM_sedimentFluxIndex_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentFluxIndex_output_stream', config_AM_sedimentFluxIndex_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentFluxIndex_directory', config_AM_sedimentFluxIndex_directory)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentFluxIndex_use_lat_lon_coords', &
config_AM_sedimentFluxIndex_use_lat_lon_coords)

      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_enable', config_AM_sedimentTransport_enable)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_compute_on_startup', &
config_AM_sedimentTransport_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_write_on_startup', &
config_AM_sedimentTransport_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_compute_interval', &
config_AM_sedimentTransport_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_output_stream', config_AM_sedimentTransport_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_directory', config_AM_sedimentTransport_directory)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_grain_size', config_AM_sedimentTransport_grain_size)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_ws_formula', config_AM_sedimentTransport_ws_formula)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_bedld_formula', config_AM_sedimentTransport_bedld_formula)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_SSC_ref_formula', &
config_AM_sedimentTransport_SSC_ref_formula)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_drag_coefficient', &
config_AM_sedimentTransport_drag_coefficient)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_erate', config_AM_sedimentTransport_erate)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_tau_ce', config_AM_sedimentTransport_tau_ce)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_tau_cd', config_AM_sedimentTransport_tau_cd)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_Manning_coef', config_AM_sedimentTransport_Manning_coef)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_grain_porosity', &
config_AM_sedimentTransport_grain_porosity)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_water_density', config_AM_sedimentTransport_water_density)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_grain_density', config_AM_sedimentTransport_grain_density)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_alpha', config_AM_sedimentTransport_alpha)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_kinematic_viscosity', &
config_AM_sedimentTransport_kinematic_viscosity)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_vertical_diffusion_coefficient', &
config_AM_sedimentTransport_vertical_diffusion_coefficient)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_bedload', config_AM_sedimentTransport_bedload)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_suspended', config_AM_sedimentTransport_suspended)
      call mpas_pool_get_config(configPool, 'config_AM_sedimentTransport_use_lat_lon_coords', &
config_AM_sedimentTransport_use_lat_lon_coords)

      call mpas_pool_get_config(configPool, 'config_AM_harmonicAnalysis_enable', config_AM_harmonicAnalysis_enable)
      call mpas_pool_get_config(configPool, 'config_AM_harmonicAnalysis_compute_interval', &
config_AM_harmonicAnalysis_compute_interval)
      call mpas_pool_get_config(configPool, 'config_AM_harmonicAnalysis_start', config_AM_harmonicAnalysis_start)
      call mpas_pool_get_config(configPool, 'config_AM_harmonicAnalysis_end', config_AM_harmonicAnalysis_end)
      call mpas_pool_get_config(configPool, 'config_AM_harmonicAnalysis_output_stream', config_AM_harmonicAnalysis_output_stream)
      call mpas_pool_get_config(configPool, 'config_AM_harmonicAnalysis_restart_stream', config_AM_harmonicAnalysis_restart_stream)
      call mpas_pool_get_config(configPool, 'config_AM_harmonicAnalysis_compute_on_startup', &
config_AM_harmonicAnalysis_compute_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_harmonicAnalysis_write_on_startup', &
config_AM_harmonicAnalysis_write_on_startup)
      call mpas_pool_get_config(configPool, 'config_AM_harmonicAnalysis_use_M2', config_AM_harmonicAnalysis_use_M2)
      call mpas_pool_get_config(configPool, 'config_AM_harmonicAnalysis_use_S2', config_AM_harmonicAnalysis_use_S2)
      call mpas_pool_get_config(configPool, 'config_AM_harmonicAnalysis_use_N2', config_AM_harmonicAnalysis_use_N2)
      call mpas_pool_get_config(configPool, 'config_AM_harmonicAnalysis_use_K2', config_AM_harmonicAnalysis_use_K2)
      call mpas_pool_get_config(configPool, 'config_AM_harmonicAnalysis_use_K1', config_AM_harmonicAnalysis_use_K1)
      call mpas_pool_get_config(configPool, 'config_AM_harmonicAnalysis_use_O1', config_AM_harmonicAnalysis_use_O1)
      call mpas_pool_get_config(configPool, 'config_AM_harmonicAnalysis_use_Q1', config_AM_harmonicAnalysis_use_Q1)
      call mpas_pool_get_config(configPool, 'config_AM_harmonicAnalysis_use_P1', config_AM_harmonicAnalysis_use_P1)

      call mpas_pool_get_config(configPool, 'config_baroclinic_channel_vert_levels', config_baroclinic_channel_vert_levels)
      call mpas_pool_get_config(configPool, 'config_baroclinic_channel_use_distances', config_baroclinic_channel_use_distances)
      call mpas_pool_get_config(configPool, 'config_baroclinic_channel_surface_temperature', &
config_baroclinic_channel_surface_temperature)
      call mpas_pool_get_config(configPool, 'config_baroclinic_channel_bottom_temperature', &
config_baroclinic_channel_bottom_temperature)
      call mpas_pool_get_config(configPool, 'config_baroclinic_channel_temperature_difference', &
config_baroclinic_channel_temperature_difference)
      call mpas_pool_get_config(configPool, 'config_baroclinic_channel_gradient_width_frac', &
config_baroclinic_channel_gradient_width_frac)
      call mpas_pool_get_config(configPool, 'config_baroclinic_channel_gradient_width_dist', &
config_baroclinic_channel_gradient_width_dist)
      call mpas_pool_get_config(configPool, 'config_baroclinic_channel_bottom_depth', config_baroclinic_channel_bottom_depth)
      call mpas_pool_get_config(configPool, 'config_baroclinic_channel_salinity', config_baroclinic_channel_salinity)
      call mpas_pool_get_config(configPool, 'config_baroclinic_channel_coriolis_parameter', &
config_baroclinic_channel_coriolis_parameter)

      call mpas_pool_get_config(configPool, 'config_lock_exchange_vert_levels', config_lock_exchange_vert_levels)
      call mpas_pool_get_config(configPool, 'config_lock_exchange_bottom_depth', config_lock_exchange_bottom_depth)
      call mpas_pool_get_config(configPool, 'config_lock_exchange_cold_temperature', config_lock_exchange_cold_temperature)
      call mpas_pool_get_config(configPool, 'config_lock_exchange_warm_temperature', config_lock_exchange_warm_temperature)
      call mpas_pool_get_config(configPool, 'config_lock_exchange_direction', config_lock_exchange_direction)
      call mpas_pool_get_config(configPool, 'config_lock_exchange_salinity', config_lock_exchange_salinity)
      call mpas_pool_get_config(configPool, 'config_lock_exchange_layer_type', config_lock_exchange_layer_type)
      call mpas_pool_get_config(configPool, 'config_lock_exchange_isopycnal_min_thickness', &
config_lock_exchange_isopycnal_min_thickness)

      call mpas_pool_get_config(configPool, 'config_internal_waves_vert_levels', config_internal_waves_vert_levels)
      call mpas_pool_get_config(configPool, 'config_internal_waves_use_distances', config_internal_waves_use_distances)
      call mpas_pool_get_config(configPool, 'config_internal_waves_surface_temperature', config_internal_waves_surface_temperature)
      call mpas_pool_get_config(configPool, 'config_internal_waves_bottom_temperature', config_internal_waves_bottom_temperature)
      call mpas_pool_get_config(configPool, 'config_internal_waves_temperature_difference', &
config_internal_waves_temperature_difference)
      call mpas_pool_get_config(configPool, 'config_internal_waves_amplitude_width_frac', &
config_internal_waves_amplitude_width_frac)
      call mpas_pool_get_config(configPool, 'config_internal_waves_amplitude_width_dist', &
config_internal_waves_amplitude_width_dist)
      call mpas_pool_get_config(configPool, 'config_internal_waves_bottom_depth', config_internal_waves_bottom_depth)
      call mpas_pool_get_config(configPool, 'config_internal_waves_salinity', config_internal_waves_salinity)
      call mpas_pool_get_config(configPool, 'config_internal_waves_layer_type', config_internal_waves_layer_type)
      call mpas_pool_get_config(configPool, 'config_internal_waves_isopycnal_displacement', &
config_internal_waves_isopycnal_displacement)

      call mpas_pool_get_config(configPool, 'config_overflow_vert_levels', config_overflow_vert_levels)
      call mpas_pool_get_config(configPool, 'config_overflow_use_distances', config_overflow_use_distances)
      call mpas_pool_get_config(configPool, 'config_overflow_bottom_depth', config_overflow_bottom_depth)
      call mpas_pool_get_config(configPool, 'config_overflow_ridge_depth', config_overflow_ridge_depth)
      call mpas_pool_get_config(configPool, 'config_overflow_plug_temperature', config_overflow_plug_temperature)
      call mpas_pool_get_config(configPool, 'config_overflow_domain_temperature', config_overflow_domain_temperature)
      call mpas_pool_get_config(configPool, 'config_overflow_salinity', config_overflow_salinity)
      call mpas_pool_get_config(configPool, 'config_overflow_plug_width_frac', config_overflow_plug_width_frac)
      call mpas_pool_get_config(configPool, 'config_overflow_slope_center_frac', config_overflow_slope_center_frac)
      call mpas_pool_get_config(configPool, 'config_overflow_slope_width_frac', config_overflow_slope_width_frac)
      call mpas_pool_get_config(configPool, 'config_overflow_plug_width_dist', config_overflow_plug_width_dist)
      call mpas_pool_get_config(configPool, 'config_overflow_slope_center_dist', config_overflow_slope_center_dist)
      call mpas_pool_get_config(configPool, 'config_overflow_slope_width_dist', config_overflow_slope_width_dist)
      call mpas_pool_get_config(configPool, 'config_overflow_layer_type', config_overflow_layer_type)
      call mpas_pool_get_config(configPool, 'config_overflow_isopycnal_min_thickness', config_overflow_isopycnal_min_thickness)

      call mpas_pool_get_config(configPool, 'config_dam_break_vert_levels', config_dam_break_vert_levels)
      call mpas_pool_get_config(configPool, 'config_dam_break_eta0', config_dam_break_eta0)
      call mpas_pool_get_config(configPool, 'config_dam_break_dc', config_dam_break_dc)
      call mpas_pool_get_config(configPool, 'config_dam_break_R0', config_dam_break_R0)
      call mpas_pool_get_config(configPool, 'config_dam_break_Xl', config_dam_break_Xl)
      call mpas_pool_get_config(configPool, 'config_dam_break_Yl', config_dam_break_Yl)
      call mpas_pool_get_config(configPool, 'config_dam_break_Inlet', config_dam_break_Inlet)

      call mpas_pool_get_config(configPool, 'config_global_ocean_minimum_depth', config_global_ocean_minimum_depth)
      call mpas_pool_get_config(configPool, 'config_global_ocean_depth_file', config_global_ocean_depth_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_depth_dimname', config_global_ocean_depth_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_depth_varname', config_global_ocean_depth_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_depth_conversion_factor', &
config_global_ocean_depth_conversion_factor)
      call mpas_pool_get_config(configPool, 'config_global_ocean_temperature_file', config_global_ocean_temperature_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_salinity_file', config_global_ocean_salinity_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_tracer_nlat_dimname', config_global_ocean_tracer_nlat_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_tracer_nlon_dimname', config_global_ocean_tracer_nlon_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_tracer_ndepth_dimname', config_global_ocean_tracer_ndepth_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_tracer_depth_conversion_factor', &
config_global_ocean_tracer_depth_conversion_factor)
      call mpas_pool_get_config(configPool, 'config_global_ocean_tracer_vert_levels', config_global_ocean_tracer_vert_levels)
      call mpas_pool_get_config(configPool, 'config_global_ocean_temperature_varname', config_global_ocean_temperature_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_salinity_varname', config_global_ocean_salinity_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_tracer_latlon_degrees', config_global_ocean_tracer_latlon_degrees)
      call mpas_pool_get_config(configPool, 'config_global_ocean_tracer_lat_varname', config_global_ocean_tracer_lat_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_tracer_lon_varname', config_global_ocean_tracer_lon_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_tracer_depth_varname', config_global_ocean_tracer_depth_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_tracer_method', config_global_ocean_tracer_method)
      call mpas_pool_get_config(configPool, 'config_global_ocean_smooth_TS_iterations', config_global_ocean_smooth_TS_iterations)
      call mpas_pool_get_config(configPool, 'config_global_ocean_swData_file', config_global_ocean_swData_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_swData_nlat_dimname', config_global_ocean_swData_nlat_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_swData_nlon_dimname', config_global_ocean_swData_nlon_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_swData_lat_varname', config_global_ocean_swData_lat_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_swData_lon_varname', config_global_ocean_swData_lon_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_swData_latlon_degrees', config_global_ocean_swData_latlon_degrees)
      call mpas_pool_get_config(configPool, 'config_global_ocean_swData_method', config_global_ocean_swData_method)
      call mpas_pool_get_config(configPool, 'config_global_ocean_chlorophyll_varname', config_global_ocean_chlorophyll_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_zenithAngle_varname', config_global_ocean_zenithAngle_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_clearSky_varname', config_global_ocean_clearSky_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_piston_velocity', config_global_ocean_piston_velocity)
      call mpas_pool_get_config(configPool, 'config_global_ocean_interior_restore_rate', config_global_ocean_interior_restore_rate)
      call mpas_pool_get_config(configPool, 'config_global_ocean_topography_file', config_global_ocean_topography_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_topography_nlat_dimname', &
config_global_ocean_topography_nlat_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_topography_nlon_dimname', &
config_global_ocean_topography_nlon_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_topography_latlon_degrees', &
config_global_ocean_topography_latlon_degrees)
      call mpas_pool_get_config(configPool, 'config_global_ocean_topography_lat_varname', &
config_global_ocean_topography_lat_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_topography_lon_varname', &
config_global_ocean_topography_lon_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_topography_varname', config_global_ocean_topography_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_topography_has_ocean_frac', &
config_global_ocean_topography_has_ocean_frac)
      call mpas_pool_get_config(configPool, 'config_global_ocean_topography_ocean_frac_varname', &
config_global_ocean_topography_ocean_frac_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_topography_method', config_global_ocean_topography_method)
      call mpas_pool_get_config(configPool, 'config_global_ocean_fill_bathymetry_holes', config_global_ocean_fill_bathymetry_holes)
      call mpas_pool_get_config(configPool, 'config_global_ocean_topography_smooth_iterations', &
config_global_ocean_topography_smooth_iterations)
      call mpas_pool_get_config(configPool, 'config_global_ocean_topography_smooth_weight', &
config_global_ocean_topography_smooth_weight)
      call mpas_pool_get_config(configPool, 'config_global_ocean_deepen_critical_passages', &
config_global_ocean_deepen_critical_passages)
      call mpas_pool_get_config(configPool, 'config_global_ocean_depress_by_land_ice', config_global_ocean_depress_by_land_ice)
      call mpas_pool_get_config(configPool, 'config_global_ocean_land_ice_topo_file', config_global_ocean_land_ice_topo_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_land_ice_topo_nlat_dimname', &
config_global_ocean_land_ice_topo_nlat_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_land_ice_topo_nlon_dimname', &
config_global_ocean_land_ice_topo_nlon_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_land_ice_topo_latlon_degrees', &
config_global_ocean_land_ice_topo_latlon_degrees)
      call mpas_pool_get_config(configPool, 'config_global_ocean_land_ice_topo_lat_varname', &
config_global_ocean_land_ice_topo_lat_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_land_ice_topo_lon_varname', &
config_global_ocean_land_ice_topo_lon_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_land_ice_topo_thickness_varname', &
config_global_ocean_land_ice_topo_thickness_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_land_ice_topo_draft_varname', &
config_global_ocean_land_ice_topo_draft_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_land_ice_topo_ice_frac_varname', &
config_global_ocean_land_ice_topo_ice_frac_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_land_ice_topo_grounded_frac_varname', &
config_global_ocean_land_ice_topo_grounded_frac_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_use_constant_land_ice_cavity_temperature', &
config_global_ocean_use_constant_land_ice_cavity_temperature)
      call mpas_pool_get_config(configPool, 'config_global_ocean_constant_land_ice_cavity_temperature', &
config_global_ocean_constant_land_ice_cavity_temperature)
      call mpas_pool_get_config(configPool, 'config_global_ocean_cull_inland_seas', config_global_ocean_cull_inland_seas)
      call mpas_pool_get_config(configPool, 'config_global_ocean_windstress_file', config_global_ocean_windstress_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_windstress_nlat_dimname', &
config_global_ocean_windstress_nlat_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_windstress_nlon_dimname', &
config_global_ocean_windstress_nlon_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_windstress_latlon_degrees', &
config_global_ocean_windstress_latlon_degrees)
      call mpas_pool_get_config(configPool, 'config_global_ocean_windstress_lat_varname', &
config_global_ocean_windstress_lat_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_windstress_lon_varname', &
config_global_ocean_windstress_lon_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_windstress_zonal_varname', &
config_global_ocean_windstress_zonal_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_windstress_meridional_varname', &
config_global_ocean_windstress_meridional_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_windstress_method', config_global_ocean_windstress_method)
      call mpas_pool_get_config(configPool, 'config_global_ocean_windstress_conversion_factor', &
config_global_ocean_windstress_conversion_factor)
      call mpas_pool_get_config(configPool, 'config_global_ocean_ecosys_file', config_global_ocean_ecosys_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_ecosys_forcing_file', config_global_ocean_ecosys_forcing_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_ecosys_nlat_dimname', config_global_ocean_ecosys_nlat_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_ecosys_nlon_dimname', config_global_ocean_ecosys_nlon_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_ecosys_ndepth_dimname', config_global_ocean_ecosys_ndepth_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_ecosys_depth_conversion_factor', &
config_global_ocean_ecosys_depth_conversion_factor)
      call mpas_pool_get_config(configPool, 'config_global_ocean_ecosys_vert_levels', config_global_ocean_ecosys_vert_levels)
      call mpas_pool_get_config(configPool, 'config_global_ocean_ecosys_lat_varname', config_global_ocean_ecosys_lat_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_ecosys_lon_varname', config_global_ocean_ecosys_lon_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_ecosys_depth_varname', config_global_ocean_ecosys_depth_varname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_ecosys_latlon_degrees', config_global_ocean_ecosys_latlon_degrees)
      call mpas_pool_get_config(configPool, 'config_global_ocean_ecosys_method', config_global_ocean_ecosys_method)
      call mpas_pool_get_config(configPool, 'config_global_ocean_ecosys_forcing_time_dimname', &
config_global_ocean_ecosys_forcing_time_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_smooth_ecosys_iterations', &
config_global_ocean_smooth_ecosys_iterations)

      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_vert_levels', config_cvmix_WSwSBF_vert_levels)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_surface_temperature', config_cvmix_WSwSBF_surface_temperature)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_surface_salinity', config_cvmix_WSwSBF_surface_salinity)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_surface_restoring_temperature', &
config_cvmix_WSwSBF_surface_restoring_temperature)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_surface_restoring_salinity', &
config_cvmix_WSwSBF_surface_restoring_salinity)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_temperature_piston_velocity', &
config_cvmix_WSwSBF_temperature_piston_velocity)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_salinity_piston_velocity', &
config_cvmix_WSwSBF_salinity_piston_velocity)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_sensible_heat_flux', config_cvmix_WSwSBF_sensible_heat_flux)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_latent_heat_flux', config_cvmix_WSwSBF_latent_heat_flux)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_shortwave_heat_flux', config_cvmix_WSwSBF_shortwave_heat_flux)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_rain_flux', config_cvmix_WSwSBF_rain_flux)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_evaporation_flux', config_cvmix_WSwSBF_evaporation_flux)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_interior_temperature_restoring_rate', &
config_cvmix_WSwSBF_interior_temperature_restoring_rate)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_interior_salinity_restoring_rate', &
config_cvmix_WSwSBF_interior_salinity_restoring_rate)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_temperature_gradient', config_cvmix_WSwSBF_temperature_gradient)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_salinity_gradient', config_cvmix_WSwSBF_salinity_gradient)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_temperature_gradient_mixed_layer', &
config_cvmix_WSwSBF_temperature_gradient_mixed_layer)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_salinity_gradient_mixed_layer', &
config_cvmix_WSwSBF_salinity_gradient_mixed_layer)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_mixed_layer_depth_temperature', &
config_cvmix_WSwSBF_mixed_layer_depth_temperature)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_mixed_layer_depth_salinity', &
config_cvmix_WSwSBF_mixed_layer_depth_salinity)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_mixed_layer_temperature_change', &
config_cvmix_WSwSBF_mixed_layer_temperature_change)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_mixed_layer_salinity_change', &
config_cvmix_WSwSBF_mixed_layer_salinity_change)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_vertical_grid', config_cvmix_WSwSBF_vertical_grid)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_bottom_depth', config_cvmix_WSwSBF_bottom_depth)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_max_windstress', config_cvmix_WSwSBF_max_windstress)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_coriolis_parameter', config_cvmix_WSwSBF_coriolis_parameter)

      call mpas_pool_get_config(configPool, 'config_iso_vert_levels', config_iso_vert_levels)
      call mpas_pool_get_config(configPool, 'config_iso_main_channel_depth', config_iso_main_channel_depth)
      call mpas_pool_get_config(configPool, 'config_iso_north_wall_lat', config_iso_north_wall_lat)
      call mpas_pool_get_config(configPool, 'config_iso_south_wall_lat', config_iso_south_wall_lat)
      call mpas_pool_get_config(configPool, 'config_iso_ridge_flag', config_iso_ridge_flag)
      call mpas_pool_get_config(configPool, 'config_iso_ridge_center_lon', config_iso_ridge_center_lon)
      call mpas_pool_get_config(configPool, 'config_iso_ridge_height', config_iso_ridge_height)
      call mpas_pool_get_config(configPool, 'config_iso_ridge_width', config_iso_ridge_width)
      call mpas_pool_get_config(configPool, 'config_iso_plateau_flag', config_iso_plateau_flag)
      call mpas_pool_get_config(configPool, 'config_iso_plateau_center_lon', config_iso_plateau_center_lon)
      call mpas_pool_get_config(configPool, 'config_iso_plateau_center_lat', config_iso_plateau_center_lat)
      call mpas_pool_get_config(configPool, 'config_iso_plateau_height', config_iso_plateau_height)
      call mpas_pool_get_config(configPool, 'config_iso_plateau_radius', config_iso_plateau_radius)
      call mpas_pool_get_config(configPool, 'config_iso_plateau_slope_width', config_iso_plateau_slope_width)
      call mpas_pool_get_config(configPool, 'config_iso_shelf_flag', config_iso_shelf_flag)
      call mpas_pool_get_config(configPool, 'config_iso_shelf_depth', config_iso_shelf_depth)
      call mpas_pool_get_config(configPool, 'config_iso_shelf_width', config_iso_shelf_width)
      call mpas_pool_get_config(configPool, 'config_iso_cont_slope_flag', config_iso_cont_slope_flag)
      call mpas_pool_get_config(configPool, 'config_iso_max_cont_slope', config_iso_max_cont_slope)
      call mpas_pool_get_config(configPool, 'config_iso_embayment_flag', config_iso_embayment_flag)
      call mpas_pool_get_config(configPool, 'config_iso_embayment_center_lon', config_iso_embayment_center_lon)
      call mpas_pool_get_config(configPool, 'config_iso_embayment_center_lat', config_iso_embayment_center_lat)
      call mpas_pool_get_config(configPool, 'config_iso_embayment_radius', config_iso_embayment_radius)
      call mpas_pool_get_config(configPool, 'config_iso_embayment_depth', config_iso_embayment_depth)
      call mpas_pool_get_config(configPool, 'config_iso_depression_flag', config_iso_depression_flag)
      call mpas_pool_get_config(configPool, 'config_iso_depression_center_lon', config_iso_depression_center_lon)
      call mpas_pool_get_config(configPool, 'config_iso_depression_south_lat', config_iso_depression_south_lat)
      call mpas_pool_get_config(configPool, 'config_iso_depression_north_lat', config_iso_depression_north_lat)
      call mpas_pool_get_config(configPool, 'config_iso_depression_width', config_iso_depression_width)
      call mpas_pool_get_config(configPool, 'config_iso_depression_depth', config_iso_depression_depth)
      call mpas_pool_get_config(configPool, 'config_iso_salinity', config_iso_salinity)
      call mpas_pool_get_config(configPool, 'config_iso_wind_stress_max', config_iso_wind_stress_max)
      call mpas_pool_get_config(configPool, 'config_iso_acc_wind', config_iso_acc_wind)
      call mpas_pool_get_config(configPool, 'config_iso_asf_wind', config_iso_asf_wind)
      call mpas_pool_get_config(configPool, 'config_iso_wind_trans', config_iso_wind_trans)
      call mpas_pool_get_config(configPool, 'config_iso_heat_flux_south', config_iso_heat_flux_south)
      call mpas_pool_get_config(configPool, 'config_iso_heat_flux_middle', config_iso_heat_flux_middle)
      call mpas_pool_get_config(configPool, 'config_iso_heat_flux_north', config_iso_heat_flux_north)
      call mpas_pool_get_config(configPool, 'config_iso_heat_flux_lat_ss', config_iso_heat_flux_lat_ss)
      call mpas_pool_get_config(configPool, 'config_iso_heat_flux_lat_sm', config_iso_heat_flux_lat_sm)
      call mpas_pool_get_config(configPool, 'config_iso_heat_flux_lat_mn', config_iso_heat_flux_lat_mn)
      call mpas_pool_get_config(configPool, 'config_iso_region1_center_lon', config_iso_region1_center_lon)
      call mpas_pool_get_config(configPool, 'config_iso_region1_center_lat', config_iso_region1_center_lat)
      call mpas_pool_get_config(configPool, 'config_iso_region2_center_lon', config_iso_region2_center_lon)
      call mpas_pool_get_config(configPool, 'config_iso_region2_center_lat', config_iso_region2_center_lat)
      call mpas_pool_get_config(configPool, 'config_iso_region3_center_lon', config_iso_region3_center_lon)
      call mpas_pool_get_config(configPool, 'config_iso_region3_center_lat', config_iso_region3_center_lat)
      call mpas_pool_get_config(configPool, 'config_iso_region4_center_lon', config_iso_region4_center_lon)
      call mpas_pool_get_config(configPool, 'config_iso_region4_center_lat', config_iso_region4_center_lat)
      call mpas_pool_get_config(configPool, 'config_iso_heat_flux_region1_flag', config_iso_heat_flux_region1_flag)
      call mpas_pool_get_config(configPool, 'config_iso_heat_flux_region1', config_iso_heat_flux_region1)
      call mpas_pool_get_config(configPool, 'config_iso_heat_flux_region1_radius', config_iso_heat_flux_region1_radius)
      call mpas_pool_get_config(configPool, 'config_iso_heat_flux_region2_flag', config_iso_heat_flux_region2_flag)
      call mpas_pool_get_config(configPool, 'config_iso_heat_flux_region2', config_iso_heat_flux_region2)
      call mpas_pool_get_config(configPool, 'config_iso_heat_flux_region2_radius', config_iso_heat_flux_region2_radius)
      call mpas_pool_get_config(configPool, 'config_iso_surface_temperature_piston_velocity', &
config_iso_surface_temperature_piston_velocity)
      call mpas_pool_get_config(configPool, 'config_iso_initial_temp_t1', config_iso_initial_temp_t1)
      call mpas_pool_get_config(configPool, 'config_iso_initial_temp_t2', config_iso_initial_temp_t2)
      call mpas_pool_get_config(configPool, 'config_iso_initial_temp_h0', config_iso_initial_temp_h0)
      call mpas_pool_get_config(configPool, 'config_iso_initial_temp_h1', config_iso_initial_temp_h1)
      call mpas_pool_get_config(configPool, 'config_iso_initial_temp_mt', config_iso_initial_temp_mt)
      call mpas_pool_get_config(configPool, 'config_iso_initial_temp_latS', config_iso_initial_temp_latS)
      call mpas_pool_get_config(configPool, 'config_iso_initial_temp_latN', config_iso_initial_temp_latN)
      call mpas_pool_get_config(configPool, 'config_iso_temperature_sponge_t1', config_iso_temperature_sponge_t1)
      call mpas_pool_get_config(configPool, 'config_iso_temperature_sponge_h1', config_iso_temperature_sponge_h1)
      call mpas_pool_get_config(configPool, 'config_iso_temperature_sponge_l1', config_iso_temperature_sponge_l1)
      call mpas_pool_get_config(configPool, 'config_iso_temperature_sponge_tau1', config_iso_temperature_sponge_tau1)
      call mpas_pool_get_config(configPool, 'config_iso_temperature_restore_region1_flag', &
config_iso_temperature_restore_region1_flag)
      call mpas_pool_get_config(configPool, 'config_iso_temperature_restore_t1', config_iso_temperature_restore_t1)
      call mpas_pool_get_config(configPool, 'config_iso_temperature_restore_lcx1', config_iso_temperature_restore_lcx1)
      call mpas_pool_get_config(configPool, 'config_iso_temperature_restore_lcy1', config_iso_temperature_restore_lcy1)
      call mpas_pool_get_config(configPool, 'config_iso_temperature_restore_region2_flag', &
config_iso_temperature_restore_region2_flag)
      call mpas_pool_get_config(configPool, 'config_iso_temperature_restore_t2', config_iso_temperature_restore_t2)
      call mpas_pool_get_config(configPool, 'config_iso_temperature_restore_lcx2', config_iso_temperature_restore_lcx2)
      call mpas_pool_get_config(configPool, 'config_iso_temperature_restore_lcy2', config_iso_temperature_restore_lcy2)
      call mpas_pool_get_config(configPool, 'config_iso_temperature_restore_region3_flag', &
config_iso_temperature_restore_region3_flag)
      call mpas_pool_get_config(configPool, 'config_iso_temperature_restore_t3', config_iso_temperature_restore_t3)
      call mpas_pool_get_config(configPool, 'config_iso_temperature_restore_lcx3', config_iso_temperature_restore_lcx3)
      call mpas_pool_get_config(configPool, 'config_iso_temperature_restore_lcy3', config_iso_temperature_restore_lcy3)
      call mpas_pool_get_config(configPool, 'config_iso_temperature_restore_region4_flag', &
config_iso_temperature_restore_region4_flag)
      call mpas_pool_get_config(configPool, 'config_iso_temperature_restore_t4', config_iso_temperature_restore_t4)
      call mpas_pool_get_config(configPool, 'config_iso_temperature_restore_lcx4', config_iso_temperature_restore_lcx4)
      call mpas_pool_get_config(configPool, 'config_iso_temperature_restore_lcy4', config_iso_temperature_restore_lcy4)

      call mpas_pool_get_config(configPool, 'config_soma_vert_levels', config_soma_vert_levels)
      call mpas_pool_get_config(configPool, 'config_soma_domain_width', config_soma_domain_width)
      call mpas_pool_get_config(configPool, 'config_soma_center_latitude', config_soma_center_latitude)
      call mpas_pool_get_config(configPool, 'config_soma_center_longitude', config_soma_center_longitude)
      call mpas_pool_get_config(configPool, 'config_soma_phi', config_soma_phi)
      call mpas_pool_get_config(configPool, 'config_soma_bottom_depth', config_soma_bottom_depth)
      call mpas_pool_get_config(configPool, 'config_soma_shelf_width', config_soma_shelf_width)
      call mpas_pool_get_config(configPool, 'config_soma_shelf_depth', config_soma_shelf_depth)
      call mpas_pool_get_config(configPool, 'config_soma_ref_density', config_soma_ref_density)
      call mpas_pool_get_config(configPool, 'config_soma_density_difference', config_soma_density_difference)
      call mpas_pool_get_config(configPool, 'config_soma_thermocline_depth', config_soma_thermocline_depth)
      call mpas_pool_get_config(configPool, 'config_soma_density_difference_linear', config_soma_density_difference_linear)
      call mpas_pool_get_config(configPool, 'config_soma_surface_temperature', config_soma_surface_temperature)
      call mpas_pool_get_config(configPool, 'config_soma_surface_salinity', config_soma_surface_salinity)
      call mpas_pool_get_config(configPool, 'config_soma_use_surface_temp_restoring', config_soma_use_surface_temp_restoring)
      call mpas_pool_get_config(configPool, 'config_soma_surface_temp_restoring_at_center_latitude', &
config_soma_surface_temp_restoring_at_center_latitude)
      call mpas_pool_get_config(configPool, 'config_soma_surface_temp_restoring_latitude_gradient', &
config_soma_surface_temp_restoring_latitude_gradient)
      call mpas_pool_get_config(configPool, 'config_soma_restoring_temp_piston_vel', config_soma_restoring_temp_piston_vel)

      call mpas_pool_get_config(configPool, 'config_ziso_vert_levels', config_ziso_vert_levels)
      call mpas_pool_get_config(configPool, 'config_ziso_add_easterly_wind_stress_ASF', config_ziso_add_easterly_wind_stress_ASF)
      call mpas_pool_get_config(configPool, 'config_ziso_wind_transition_position', config_ziso_wind_transition_position)
      call mpas_pool_get_config(configPool, 'config_ziso_antarctic_shelf_front_width', config_ziso_antarctic_shelf_front_width)
      call mpas_pool_get_config(configPool, 'config_ziso_wind_stress_shelf_front_max', config_ziso_wind_stress_shelf_front_max)
      call mpas_pool_get_config(configPool, 'config_ziso_use_slopping_bathymetry', config_ziso_use_slopping_bathymetry)
      call mpas_pool_get_config(configPool, 'config_ziso_meridional_extent', config_ziso_meridional_extent)
      call mpas_pool_get_config(configPool, 'config_ziso_zonal_extent', config_ziso_zonal_extent)
      call mpas_pool_get_config(configPool, 'config_ziso_bottom_depth', config_ziso_bottom_depth)
      call mpas_pool_get_config(configPool, 'config_ziso_shelf_depth', config_ziso_shelf_depth)
      call mpas_pool_get_config(configPool, 'config_ziso_slope_half_width', config_ziso_slope_half_width)
      call mpas_pool_get_config(configPool, 'config_ziso_slope_center_position', config_ziso_slope_center_position)
      call mpas_pool_get_config(configPool, 'config_ziso_reference_coriolis', config_ziso_reference_coriolis)
      call mpas_pool_get_config(configPool, 'config_ziso_coriolis_gradient', config_ziso_coriolis_gradient)
      call mpas_pool_get_config(configPool, 'config_ziso_wind_stress_max', config_ziso_wind_stress_max)
      call mpas_pool_get_config(configPool, 'config_ziso_mean_restoring_temp', config_ziso_mean_restoring_temp)
      call mpas_pool_get_config(configPool, 'config_ziso_restoring_temp_dev_ta', config_ziso_restoring_temp_dev_ta)
      call mpas_pool_get_config(configPool, 'config_ziso_restoring_temp_dev_tb', config_ziso_restoring_temp_dev_tb)
      call mpas_pool_get_config(configPool, 'config_ziso_restoring_temp_tau', config_ziso_restoring_temp_tau)
      call mpas_pool_get_config(configPool, 'config_ziso_restoring_temp_piston_vel', config_ziso_restoring_temp_piston_vel)
      call mpas_pool_get_config(configPool, 'config_ziso_restoring_temp_ze', config_ziso_restoring_temp_ze)
      call mpas_pool_get_config(configPool, 'config_ziso_restoring_sponge_l', config_ziso_restoring_sponge_l)
      call mpas_pool_get_config(configPool, 'config_ziso_initial_temp_t1', config_ziso_initial_temp_t1)
      call mpas_pool_get_config(configPool, 'config_ziso_initial_temp_t2', config_ziso_initial_temp_t2)
      call mpas_pool_get_config(configPool, 'config_ziso_initial_temp_h1', config_ziso_initial_temp_h1)
      call mpas_pool_get_config(configPool, 'config_ziso_initial_temp_mt', config_ziso_initial_temp_mt)
      call mpas_pool_get_config(configPool, 'config_ziso_frazil_enable', config_ziso_frazil_enable)
      call mpas_pool_get_config(configPool, 'config_ziso_frazil_temperature_anomaly', config_ziso_frazil_temperature_anomaly)

      call mpas_pool_get_config(configPool, 'config_sub_ice_shelf_2D_vert_levels', config_sub_ice_shelf_2D_vert_levels)
      call mpas_pool_get_config(configPool, 'config_sub_ice_shelf_2D_bottom_depth', config_sub_ice_shelf_2D_bottom_depth)
      call mpas_pool_get_config(configPool, 'config_sub_ice_shelf_2D_cavity_thickness', config_sub_ice_shelf_2D_cavity_thickness)
      call mpas_pool_get_config(configPool, 'config_sub_ice_shelf_2D_slope_height', config_sub_ice_shelf_2D_slope_height)
      call mpas_pool_get_config(configPool, 'config_sub_ice_shelf_2D_edge_width', config_sub_ice_shelf_2D_edge_width)
      call mpas_pool_get_config(configPool, 'config_sub_ice_shelf_2D_y1', config_sub_ice_shelf_2D_y1)
      call mpas_pool_get_config(configPool, 'config_sub_ice_shelf_2D_y2', config_sub_ice_shelf_2D_y2)
      call mpas_pool_get_config(configPool, 'config_sub_ice_shelf_2D_temperature', config_sub_ice_shelf_2D_temperature)
      call mpas_pool_get_config(configPool, 'config_sub_ice_shelf_2D_surface_salinity', config_sub_ice_shelf_2D_surface_salinity)
      call mpas_pool_get_config(configPool, 'config_sub_ice_shelf_2D_bottom_salinity', config_sub_ice_shelf_2D_bottom_salinity)

      call mpas_pool_get_config(configPool, 'config_periodic_planar_vert_levels', config_periodic_planar_vert_levels)
      call mpas_pool_get_config(configPool, 'config_periodic_planar_bottom_depth', config_periodic_planar_bottom_depth)
      call mpas_pool_get_config(configPool, 'config_periodic_planar_velocity_strength', config_periodic_planar_velocity_strength)

      call mpas_pool_get_config(configPool, 'config_ecosys_column_vert_levels', config_ecosys_column_vert_levels)
      call mpas_pool_get_config(configPool, 'config_ecosys_column_vertical_grid', config_ecosys_column_vertical_grid)
      call mpas_pool_get_config(configPool, 'config_ecosys_column_TS_filename', config_ecosys_column_TS_filename)
      call mpas_pool_get_config(configPool, 'config_ecosys_column_ecosys_filename', config_ecosys_column_ecosys_filename)
      call mpas_pool_get_config(configPool, 'config_ecosys_column_bottom_depth', config_ecosys_column_bottom_depth)

      call mpas_pool_get_config(configPool, 'config_sea_mount_vert_levels', config_sea_mount_vert_levels)
      call mpas_pool_get_config(configPool, 'config_sea_mount_layer_type', config_sea_mount_layer_type)
      call mpas_pool_get_config(configPool, 'config_sea_mount_stratification_type', config_sea_mount_stratification_type)
      call mpas_pool_get_config(configPool, 'config_sea_mount_density_coef_linear', config_sea_mount_density_coef_linear)
      call mpas_pool_get_config(configPool, 'config_sea_mount_density_coef_exp', config_sea_mount_density_coef_exp)
      call mpas_pool_get_config(configPool, 'config_sea_mount_density_gradient_linear', config_sea_mount_density_gradient_linear)
      call mpas_pool_get_config(configPool, 'config_sea_mount_density_gradient_exp', config_sea_mount_density_gradient_exp)
      call mpas_pool_get_config(configPool, 'config_sea_mount_density_depth_linear', config_sea_mount_density_depth_linear)
      call mpas_pool_get_config(configPool, 'config_sea_mount_density_depth_exp', config_sea_mount_density_depth_exp)
      call mpas_pool_get_config(configPool, 'config_sea_mount_density_ref', config_sea_mount_density_ref)
      call mpas_pool_get_config(configPool, 'config_sea_mount_density_Tref', config_sea_mount_density_Tref)
      call mpas_pool_get_config(configPool, 'config_sea_mount_density_alpha', config_sea_mount_density_alpha)
      call mpas_pool_get_config(configPool, 'config_sea_mount_bottom_depth', config_sea_mount_bottom_depth)
      call mpas_pool_get_config(configPool, 'config_sea_mount_height', config_sea_mount_height)
      call mpas_pool_get_config(configPool, 'config_sea_mount_radius', config_sea_mount_radius)
      call mpas_pool_get_config(configPool, 'config_sea_mount_width', config_sea_mount_width)
      call mpas_pool_get_config(configPool, 'config_sea_mount_salinity', config_sea_mount_salinity)
      call mpas_pool_get_config(configPool, 'config_sea_mount_coriolis_parameter', config_sea_mount_coriolis_parameter)

      call mpas_pool_get_config(configPool, 'config_isomip_vert_levels', config_isomip_vert_levels)
      call mpas_pool_get_config(configPool, 'config_isomip_vertical_level_distribution', config_isomip_vertical_level_distribution)
      call mpas_pool_get_config(configPool, 'config_isomip_bottom_depth', config_isomip_bottom_depth)
      call mpas_pool_get_config(configPool, 'config_isomip_temperature', config_isomip_temperature)
      call mpas_pool_get_config(configPool, 'config_isomip_salinity', config_isomip_salinity)
      call mpas_pool_get_config(configPool, 'config_isomip_restoring_temperature', config_isomip_restoring_temperature)
      call mpas_pool_get_config(configPool, 'config_isomip_temperature_piston_velocity', config_isomip_temperature_piston_velocity)
      call mpas_pool_get_config(configPool, 'config_isomip_restoring_salinity', config_isomip_restoring_salinity)
      call mpas_pool_get_config(configPool, 'config_isomip_salinity_piston_velocity', config_isomip_salinity_piston_velocity)
      call mpas_pool_get_config(configPool, 'config_isomip_coriolis_parameter', config_isomip_coriolis_parameter)
      call mpas_pool_get_config(configPool, 'config_isomip_southern_boundary', config_isomip_southern_boundary)
      call mpas_pool_get_config(configPool, 'config_isomip_northern_boundary', config_isomip_northern_boundary)
      call mpas_pool_get_config(configPool, 'config_isomip_western_boundary', config_isomip_western_boundary)
      call mpas_pool_get_config(configPool, 'config_isomip_eastern_boundary', config_isomip_eastern_boundary)
      call mpas_pool_get_config(configPool, 'config_isomip_y1', config_isomip_y1)
      call mpas_pool_get_config(configPool, 'config_isomip_z1', config_isomip_z1)
      call mpas_pool_get_config(configPool, 'config_isomip_ice_fraction1', config_isomip_ice_fraction1)
      call mpas_pool_get_config(configPool, 'config_isomip_y2', config_isomip_y2)
      call mpas_pool_get_config(configPool, 'config_isomip_z2', config_isomip_z2)
      call mpas_pool_get_config(configPool, 'config_isomip_ice_fraction2', config_isomip_ice_fraction2)
      call mpas_pool_get_config(configPool, 'config_isomip_y3', config_isomip_y3)
      call mpas_pool_get_config(configPool, 'config_isomip_z3', config_isomip_z3)
      call mpas_pool_get_config(configPool, 'config_isomip_ice_fraction3', config_isomip_ice_fraction3)

      call mpas_pool_get_config(configPool, 'config_isomip_plus_vert_levels', config_isomip_plus_vert_levels)
      call mpas_pool_get_config(configPool, 'config_isomip_plus_vertical_level_distribution', &
config_isomip_plus_vertical_level_distribution)
      call mpas_pool_get_config(configPool, 'config_isomip_plus_max_bottom_depth', config_isomip_plus_max_bottom_depth)
      call mpas_pool_get_config(configPool, 'config_isomip_plus_minimum_levels', config_isomip_plus_minimum_levels)
      call mpas_pool_get_config(configPool, 'config_isomip_plus_min_column_thickness', config_isomip_plus_min_column_thickness)
      call mpas_pool_get_config(configPool, 'config_isomip_plus_min_ocean_fraction', config_isomip_plus_min_ocean_fraction)
      call mpas_pool_get_config(configPool, 'config_isomip_plus_topography_file', config_isomip_plus_topography_file)
      call mpas_pool_get_config(configPool, 'config_isomip_plus_init_top_temp', config_isomip_plus_init_top_temp)
      call mpas_pool_get_config(configPool, 'config_isomip_plus_init_bot_temp', config_isomip_plus_init_bot_temp)
      call mpas_pool_get_config(configPool, 'config_isomip_plus_init_top_sal', config_isomip_plus_init_top_sal)
      call mpas_pool_get_config(configPool, 'config_isomip_plus_init_bot_sal', config_isomip_plus_init_bot_sal)
      call mpas_pool_get_config(configPool, 'config_isomip_plus_restore_top_temp', config_isomip_plus_restore_top_temp)
      call mpas_pool_get_config(configPool, 'config_isomip_plus_restore_bot_temp', config_isomip_plus_restore_bot_temp)
      call mpas_pool_get_config(configPool, 'config_isomip_plus_restore_top_sal', config_isomip_plus_restore_top_sal)
      call mpas_pool_get_config(configPool, 'config_isomip_plus_restore_bot_sal', config_isomip_plus_restore_bot_sal)
      call mpas_pool_get_config(configPool, 'config_isomip_plus_restore_rate', config_isomip_plus_restore_rate)
      call mpas_pool_get_config(configPool, 'config_isomip_plus_restore_evap_rate', config_isomip_plus_restore_evap_rate)
      call mpas_pool_get_config(configPool, 'config_isomip_plus_restore_xMin', config_isomip_plus_restore_xMin)
      call mpas_pool_get_config(configPool, 'config_isomip_plus_restore_xMax', config_isomip_plus_restore_xMax)
      call mpas_pool_get_config(configPool, 'config_isomip_plus_coriolis_parameter', config_isomip_plus_coriolis_parameter)
      call mpas_pool_get_config(configPool, 'config_isomip_plus_effective_density', config_isomip_plus_effective_density)

      call mpas_pool_get_config(configPool, 'config_hurricane_vert_levels', config_hurricane_vert_levels)
      call mpas_pool_get_config(configPool, 'config_hurricane_min_depth', config_hurricane_min_depth)
      call mpas_pool_get_config(configPool, 'config_hurricane_max_depth', config_hurricane_max_depth)
      call mpas_pool_get_config(configPool, 'config_hurricane_gaussian_hump_amplitude', config_hurricane_gaussian_hump_amplitude)
      call mpas_pool_get_config(configPool, 'config_hurricane_use_gaussian_hump', config_hurricane_use_gaussian_hump)
      call mpas_pool_get_config(configPool, 'config_hurricane_gaussian_lon_center', config_hurricane_gaussian_lon_center)
      call mpas_pool_get_config(configPool, 'config_hurricane_gaussian_lat_center', config_hurricane_gaussian_lat_center)
      call mpas_pool_get_config(configPool, 'config_hurricane_gaussian_width', config_hurricane_gaussian_width)
      call mpas_pool_get_config(configPool, 'config_hurricane_gaussian_amplitude', config_hurricane_gaussian_amplitude)
      call mpas_pool_get_config(configPool, 'config_hurricane_gaussian_slr_amp', config_hurricane_gaussian_slr_amp)
      call mpas_pool_get_config(configPool, 'config_hurricane_land_z_limit', config_hurricane_land_z_limit)
      call mpas_pool_get_config(configPool, 'config_hurricane_marsh_z_limit', config_hurricane_marsh_z_limit)
      call mpas_pool_get_config(configPool, 'config_hurricane_land_drag', config_hurricane_land_drag)
      call mpas_pool_get_config(configPool, 'config_hurricane_marsh_drag', config_hurricane_marsh_drag)
      call mpas_pool_get_config(configPool, 'config_hurricane_channel_drag', config_hurricane_channel_drag)
      call mpas_pool_get_config(configPool, 'config_hurricane_sea_level_rise_adjustment', &
config_hurricane_sea_level_rise_adjustment)

      call mpas_pool_get_config(configPool, 'config_tidal_boundary_vert_levels', config_tidal_boundary_vert_levels)
      call mpas_pool_get_config(configPool, 'config_tidal_boundary_min_vert_levels', config_tidal_boundary_min_vert_levels)
      call mpas_pool_get_config(configPool, 'config_tidal_boundary_layer_type', config_tidal_boundary_layer_type)
      call mpas_pool_get_config(configPool, 'config_tidal_boundary_right_bottom_depth', config_tidal_boundary_right_bottom_depth)
      call mpas_pool_get_config(configPool, 'config_tidal_start_dry', config_tidal_start_dry)
      call mpas_pool_get_config(configPool, 'config_tidal_boundary_use_distances', config_tidal_boundary_use_distances)
      call mpas_pool_get_config(configPool, 'config_tidal_boundary_left_value', config_tidal_boundary_left_value)
      call mpas_pool_get_config(configPool, 'config_tidal_boundary_right_value', config_tidal_boundary_right_value)
      call mpas_pool_get_config(configPool, 'config_tidal_boundary_left_bottom_depth', config_tidal_boundary_left_bottom_depth)
      call mpas_pool_get_config(configPool, 'config_tidal_boundary_salinity', config_tidal_boundary_salinity)
      call mpas_pool_get_config(configPool, 'config_tidal_boundary_domain_temperature', config_tidal_boundary_domain_temperature)
      call mpas_pool_get_config(configPool, 'config_tidal_boundary_plug_temperature', config_tidal_boundary_plug_temperature)
      call mpas_pool_get_config(configPool, 'config_tidal_boundary_plug_width_frac', config_tidal_boundary_plug_width_frac)
      call mpas_pool_get_config(configPool, 'config_tidal_forcing_left_Cd_or_n', config_tidal_forcing_left_Cd_or_n)
      call mpas_pool_get_config(configPool, 'config_tidal_forcing_right_Cd_or_n', config_tidal_forcing_right_Cd_or_n)
      call mpas_pool_get_config(configPool, 'config_use_idealized_transect', config_use_idealized_transect)
      call mpas_pool_get_config(configPool, 'config_idealized_transect_Lshore', config_idealized_transect_Lshore)
      call mpas_pool_get_config(configPool, 'config_idealized_transect_Sshore', config_idealized_transect_Sshore)
      call mpas_pool_get_config(configPool, 'config_idealized_transect_Lcoast', config_idealized_transect_Lcoast)
      call mpas_pool_get_config(configPool, 'config_idealized_transect_Scoast', config_idealized_transect_Scoast)
      call mpas_pool_get_config(configPool, 'config_idealized_transect_Lmarsh', config_idealized_transect_Lmarsh)
      call mpas_pool_get_config(configPool, 'config_idealized_transect_Smarsh', config_idealized_transect_Smarsh)
      call mpas_pool_get_config(configPool, 'config_idealized_transect_roughness', config_idealized_transect_roughness)
      call mpas_pool_get_config(configPool, 'config_idealized_transect_roughness_marsh', config_idealized_transect_roughness_marsh)
      call mpas_pool_get_config(configPool, 'config_idealized_vegetation_diameter', config_idealized_vegetation_diameter)
      call mpas_pool_get_config(configPool, 'config_idealized_vegetation_height', config_idealized_vegetation_height)
      call mpas_pool_get_config(configPool, 'config_idealized_vegetation_density', config_idealized_vegetation_density)

      call mpas_pool_get_config(configPool, 'config_cosine_bell_temperature', config_cosine_bell_temperature)
      call mpas_pool_get_config(configPool, 'config_cosine_bell_salinity', config_cosine_bell_salinity)
      call mpas_pool_get_config(configPool, 'config_cosine_bell_lat_center', config_cosine_bell_lat_center)
      call mpas_pool_get_config(configPool, 'config_cosine_bell_lon_center', config_cosine_bell_lon_center)
      call mpas_pool_get_config(configPool, 'config_cosine_bell_psi0', config_cosine_bell_psi0)
      call mpas_pool_get_config(configPool, 'config_cosine_bell_radius', config_cosine_bell_radius)
      call mpas_pool_get_config(configPool, 'config_cosine_bell_vel_pd', config_cosine_bell_vel_pd)

      call mpas_pool_get_config(configPool, 'config_mixed_layer_eddy_vert_levels', config_mixed_layer_eddy_vert_levels)
      call mpas_pool_get_config(configPool, 'config_mixed_layer_eddy_bottom_depth', config_mixed_layer_eddy_bottom_depth)
      call mpas_pool_get_config(configPool, 'config_mixed_layer_eddy_mixed_layer_depth', config_mixed_layer_eddy_mixed_layer_depth)
      call mpas_pool_get_config(configPool, 'config_mixed_layer_eddy_base_temperature', config_mixed_layer_eddy_base_temperature)
      call mpas_pool_get_config(configPool, 'config_mixed_layer_eddy_temperature_stratification_mixed_layer', &
config_mixed_layer_eddy_temperature_stratification_mixed_layer)
      call mpas_pool_get_config(configPool, 'config_mixed_layer_eddy_temperature_stratification_interior', &
config_mixed_layer_eddy_temperature_stratification_interior)
      call mpas_pool_get_config(configPool, 'config_mixed_layer_eddy_temperature_horizontal_gradient', &
config_mixed_layer_eddy_temperature_horizontal_gradient)
      call mpas_pool_get_config(configPool, 'config_mixed_layer_eddy_temperature_front_width', &
config_mixed_layer_eddy_temperature_front_width)
      call mpas_pool_get_config(configPool, 'config_mixed_layer_eddy_temperature_perturbation_magnitude', &
config_mixed_layer_eddy_temperature_perturbation_magnitude)
      call mpas_pool_get_config(configPool, 'config_mixed_layer_eddy_salinity', config_mixed_layer_eddy_salinity)
      call mpas_pool_get_config(configPool, 'config_mixed_layer_eddy_two_fronts', config_mixed_layer_eddy_two_fronts)
      call mpas_pool_get_config(configPool, 'config_mixed_layer_eddy_restoring_width', config_mixed_layer_eddy_restoring_width)
      call mpas_pool_get_config(configPool, 'config_mixed_layer_eddy_restoring_tau', config_mixed_layer_eddy_restoring_tau)
      call mpas_pool_get_config(configPool, 'config_mixed_layer_eddy_heat_flux', config_mixed_layer_eddy_heat_flux)
      call mpas_pool_get_config(configPool, 'config_mixed_layer_eddy_evaporation_flux', config_mixed_layer_eddy_evaporation_flux)
      call mpas_pool_get_config(configPool, 'config_mixed_layer_eddy_wind_stress_zonal', config_mixed_layer_eddy_wind_stress_zonal)
      call mpas_pool_get_config(configPool, 'config_mixed_layer_eddy_wind_stress_meridional', &
config_mixed_layer_eddy_wind_stress_meridional)
      call mpas_pool_get_config(configPool, 'config_mixed_layer_eddy_coriolis_parameter', &
config_mixed_layer_eddy_coriolis_parameter)

      call mpas_pool_get_config(configPool, 'config_alter_ICs_for_pcs', config_alter_ICs_for_pcs)
      call mpas_pool_get_config(configPool, 'config_pc_alteration_type', config_pc_alteration_type)
      call mpas_pool_get_config(configPool, 'config_min_pc_fraction', config_min_pc_fraction)

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)
      call mpas_pool_get_config(configPool, 'config_expand_sphere', config_expand_sphere)
      call mpas_pool_get_config(configPool, 'config_realistic_coriolis_parameter', config_realistic_coriolis_parameter)
      call mpas_pool_get_config(configPool, 'config_write_cull_cell_mask', config_write_cull_cell_mask)
      call mpas_pool_get_config(configPool, 'config_vertical_grid', config_vertical_grid)

      call mpas_pool_get_config(configPool, 'config_1dCVTgenerator_stretch1', config_1dCVTgenerator_stretch1)
      call mpas_pool_get_config(configPool, 'config_1dCVTgenerator_stretch2', config_1dCVTgenerator_stretch2)
      call mpas_pool_get_config(configPool, 'config_1dCVTgenerator_dzSeed', config_1dCVTgenerator_dzSeed)

      call mpas_pool_get_config(configPool, 'config_init_vertical_grid_type', config_init_vertical_grid_type)

      call mpas_pool_get_config(configPool, 'config_rx1_outer_iter_count', config_rx1_outer_iter_count)
      call mpas_pool_get_config(configPool, 'config_rx1_inner_iter_count', config_rx1_inner_iter_count)
      call mpas_pool_get_config(configPool, 'config_rx1_init_inner_weight', config_rx1_init_inner_weight)
      call mpas_pool_get_config(configPool, 'config_rx1_max', config_rx1_max)
      call mpas_pool_get_config(configPool, 'config_rx1_horiz_smooth_weight', config_rx1_horiz_smooth_weight)
      call mpas_pool_get_config(configPool, 'config_rx1_vert_smooth_weight', config_rx1_vert_smooth_weight)
      call mpas_pool_get_config(configPool, 'config_rx1_slope_weight', config_rx1_slope_weight)
      call mpas_pool_get_config(configPool, 'config_rx1_zstar_weight', config_rx1_zstar_weight)
      call mpas_pool_get_config(configPool, 'config_rx1_horiz_smooth_open_ocean_cells', config_rx1_horiz_smooth_open_ocean_cells)
      call mpas_pool_get_config(configPool, 'config_rx1_min_levels', config_rx1_min_levels)
      call mpas_pool_get_config(configPool, 'config_rx1_min_layer_thickness', config_rx1_min_layer_thickness)


   end subroutine ocn_config_init!}}}

!***********************************************************************

end module ocn_config

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker

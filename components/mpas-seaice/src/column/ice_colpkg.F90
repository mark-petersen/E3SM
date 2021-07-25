!  SVN:$Id: ice_colpkg.F90 1175 2017-03-02 19:53:26Z akt $
!=========================================================================
!
! flags and interface routines for the column package
!
! authors: Elizabeth C. Hunke, LANL

      module ice_colpkg

      use ice_kinds_mod
      use ice_colpkg_shared ! namelist and other parameters
      use ice_warnings, only: add_warning

      implicit none

      private

      ! initialization
      public :: &
           colpkg_init_itd, &
           colpkg_init_itd_hist, &
           colpkg_init_thermo, &
           colpkg_init_orbit, &
           colpkg_init_trcr, &
           colpkg_init_bgc, &
           colpkg_init_zbgc, &
           colpkg_init_hbrine, &
           colpkg_init_zsalinity, &
           colpkg_init_ocean_conc, &
           colpkg_init_OceanConcArray, &
           colpkg_init_bgc_trcr, &
           colpkg_init_parameters, &
           colpkg_init_tracer_flags, &
           colpkg_init_tracer_indices, &
           colpkg_init_tracer_numbers

      ! time stepping
      public :: &
           colpkg_step_snow, &
           colpkg_step_therm1, &
           colpkg_biogeochemistry, &
           colpkg_step_therm2, &
           colpkg_prep_radiation, &
           colpkg_step_radiation, &
           colpkg_step_ridge

      ! other column routines
      public :: &
           colpkg_aggregate, &
           colpkg_ice_strength, &
           colpkg_atm_boundary, &
           colpkg_ocn_mixed_layer

      ! temperature inquiry functions
      public :: &
           colpkg_ice_temperature, &
           colpkg_snow_temperature, &
           colpkg_liquidus_temperature, &
           colpkg_sea_freezing_temperature, &
           colpkg_enthalpy_ice, &
           colpkg_enthalpy_snow, &
           colpkg_salinity_profile

      ! warning messages
      public :: &
           colpkg_clear_warnings, &
           colpkg_get_warnings, &
           colpkg_print_warnings


!=======================================================================

      contains

!=======================================================================
!     Initialization routines
!=======================================================================

! Initialize area fraction and thickness boundaries for the itd model
!
! authors: William H. Lipscomb and Elizabeth C. Hunke, LANL
!          C. M. Bitz, UW

      subroutine colpkg_init_itd(ncat, hin_max_R4, l_stop, stop_label)

      use ice_colpkg_shared, only: kcatbound, kitd
      use ice_therm_shared, only: hi_min
      use ice_constants_colpkg, only: p01, p1, c0, c1, c2, c3, c15, c25, c100

      integer (kind=int_kind), intent(in) :: &
           ncat ! number of thickness categories

      real (kind=real_kind), intent(out) :: &
           hin_max_R4(0:ncat)  ! category limits (m)

      real (kind=dbl_kind) :: &
           hin_max(0:ncat)  ! category limits (m)
      logical (kind=log_kind), intent(inout) :: &
         l_stop          ! if true, print diagnostics and abort model

      character (len=*), intent(out) :: &
         stop_label   ! abort error message

      ! local variables

      integer (kind=int_kind) :: &
           n    ! thickness category index

      real (kind=dbl_kind) :: &
           cc1, cc2, cc3, & ! parameters for kcatbound = 0
           x1           , &
           rn           , & ! real(n)
           rncat        , & ! real(ncat)
           d1           , & ! parameters for kcatbound = 1 (m)
           d2           , & !
           b1           , & ! parameters for kcatbound = 3
           b2           , & !
           b3

      real (kind=dbl_kind), dimension(5) :: wmo5 ! data for wmo itd
      real (kind=dbl_kind), dimension(6) :: wmo6 ! data for wmo itd
      real (kind=dbl_kind), dimension(7) :: wmo7 ! data for wmo itd

      hin_max = real(hin_max_R4, kind=dbl_kind)

      l_stop = .false.

      rncat = real(ncat, kind=dbl_kind)
      d1 = 3.0_dbl_kind / rncat
      d2 = 0.5_dbl_kind / rncat
      b1 = p1         ! asymptotic category width (m)
      b2 = c3         ! thickness for which participation function is small (m)
      b3 = max(rncat*(rncat-1), c2*b2/b1)

      hi_min = p01    ! minimum ice thickness allowed (m) for thermo
                      ! note hi_min is reset to 0.1 for kitd=0, below

      !-----------------------------------------------------------------
      ! Choose category boundaries based on one of four options.
      !
      ! The first formula (kcatbound = 0) was used in Lipscomb (2001) 
      !  and in CICE versions 3.0 and 3.1.
      !
      ! The second formula is more user-friendly in the sense that it
      !  is easy to obtain round numbers for category boundaries:
      !
      !    H(n) = n * [d1 + d2*(n-1)] 
      ! 
      ! Default values are d1 = 300/ncat, d2 = 50/ncat.
      ! For ncat = 5, boundaries in cm are 60, 140, 240, 360, which are 
      !  close to the standard values given by the first formula.
      ! For ncat = 10, boundaries in cm are 30, 70, 120, 180, 250, 330,
      !  420, 520, 630.    
      !
      ! The third option provides support for World Meteorological
      !  Organization classification based on thickness.  The full
      !  WMO thickness distribution is used if ncat = 7;  if ncat=5 
      !  or ncat = 6, some of the thinner categories are combined.
      ! For ncat = 5,  boundaries are         30, 70, 120, 200, >200 cm.
      ! For ncat = 6,  boundaries are     15, 30, 70, 120, 200, >200 cm.
      ! For ncat = 7,  boundaries are 10, 15, 30, 70, 120, 200, >200 cm.
      !
      ! The fourth formula asymptotes to a particular category width as
      ! the number of categories increases, given by the parameter b1.
      ! The parameter b3 is computed so that the category boundaries
      ! are even numbers.
      !
      !    H(n) = b1 * [n + b3*n*(n+1)/(2*N*(N-1))] for N=ncat
      !
      ! kcatbound=-1 is available only for 1-category runs, with
      ! boundaries 0 and 100 m.
      !-----------------------------------------------------------------

      if (kcatbound == -1) then ! single category
         hin_max(0) = c0
         hin_max(1) = c100

      elseif (kcatbound == 0) then   ! original scheme

         if (kitd == 1) then
            ! linear remapping itd category limits
            cc1 = c3/rncat
            cc2 = c15*cc1
            cc3 = c3

            hin_max(0) = c0     ! minimum ice thickness, m
         else
            ! delta function itd category limits
#ifndef CCSMCOUPLED
            hi_min = p1    ! minimum ice thickness allowed (m) for thermo
#endif
            cc1 = max(1.1_dbl_kind/rncat,hi_min)
            cc2 = c25*cc1
            cc3 = 2.25_dbl_kind

            ! hin_max(0) should not be zero
            ! use some caution in making it less than 0.10
            hin_max(0) = hi_min ! minimum ice thickness, m
         endif                  ! kitd

         do n = 1, ncat
            x1 = real(n-1,kind=dbl_kind) / rncat
            hin_max(n) = hin_max(n-1) &
                       + cc1 + cc2*(c1 + tanh(cc3*(x1-c1)))
         enddo

      elseif (kcatbound == 1) then  ! new scheme

         hin_max(0) = c0
         do n = 1, ncat
            rn = real(n, kind=dbl_kind)
            hin_max(n) = rn * (d1 + (rn-c1)*d2)
         enddo

      elseif (kcatbound == 2) then  ! WMO standard

        if (ncat == 5) then
         ! thinnest 3 categories combined
         data wmo5 / 0.30_dbl_kind, 0.70_dbl_kind, &
                    1.20_dbl_kind, 2.00_dbl_kind,  &
                    999._dbl_kind  /
         hin_max(0) = c0
         do n = 1, ncat
            hin_max(n) = wmo5(n)
         enddo
       elseif (ncat == 6) then
         ! thinnest 2 categories combined
         data wmo6 / 0.15_dbl_kind, &
                    0.30_dbl_kind, 0.70_dbl_kind,  &
                    1.20_dbl_kind, 2.00_dbl_kind,  &
                    999._dbl_kind /
!echmod wmo6a
!         data wmo6 /0.30_dbl_kind, 0.70_dbl_kind,  &
!                    1.20_dbl_kind, 2.00_dbl_kind,  &
!                    4.56729_dbl_kind, &
!                    999._dbl_kind /

         hin_max(0) = c0
         do n = 1, ncat
            hin_max(n) = wmo6(n)
         enddo
       elseif (ncat == 7) then
         ! all thickness categories 
         data wmo7 / 0.10_dbl_kind, 0.15_dbl_kind, &
                    0.30_dbl_kind, 0.70_dbl_kind,  &
                    1.20_dbl_kind, 2.00_dbl_kind,  &
                    999._dbl_kind  /
         hin_max(0) = c0
         do n = 1, ncat
            hin_max(n) = wmo7(n)
         enddo
       else
         stop_label = 'kcatbound=2 (WMO) must have ncat=5, 6 or 7'
         l_stop = .true. 
         return
       endif

      elseif (kcatbound == 3) then  ! asymptotic scheme

         hin_max(0) = c0
         do n = 1, ncat
            rn = real(n, kind=dbl_kind)
            hin_max(n) = b1 * (rn + b3*rn*(rn+c1)/(c2*rncat*(rncat-c1)))
         enddo

      endif ! kcatbound

      if (kitd == 1) then
         hin_max(ncat) = 999.9_dbl_kind ! arbitrary big number
      endif
      hin_max_R4 = real(hin_max, kind=real_kind)

      end subroutine colpkg_init_itd

!=======================================================================

! Initialize area fraction and thickness boundaries for the itd model
!
! authors: William H. Lipscomb and Elizabeth C. Hunke, LANL
!          C. M. Bitz, UW

      subroutine colpkg_init_itd_hist (ncat, hin_max, c_hi_range)

      use ice_colpkg_shared, only: kcatbound, kitd
      use ice_constants_colpkg, only: p01, p1, c2, c3, c15, c25, c100

      integer (kind=int_kind), intent(in) :: &
           ncat ! number of thickness categories

      real (kind=dbl_kind), intent(in) :: &
           hin_max(0:ncat)  ! category limits (m)

      character (len=35), intent(out) :: &
           c_hi_range(ncat) ! string for history output

      ! local variables

      integer (kind=int_kind) :: &
           n    ! thickness category index

      character(len=8) :: c_hinmax1,c_hinmax2
      character(len=2) :: c_nc

      character(len=char_len_long) :: &
           warning ! warning message

         write(warning,*) ' '
         call add_warning(warning)
         write(warning,*) 'hin_max(n-1) < Cat n < hin_max(n)'
         call add_warning(warning)
         do n = 1, ncat
            write(warning,*) hin_max(n-1),' < Cat ',n, ' < ',hin_max(n)
            call add_warning(warning)
            ! Write integer n to character string
            write (c_nc, '(i2)') n    

            ! Write hin_max to character string
            write (c_hinmax1, '(f6.3)') hin_max(n-1)
            write (c_hinmax2, '(f6.3)') hin_max(n)

            ! Save character string to write to history file
            c_hi_range(n)=c_hinmax1//'m < hi Cat '//c_nc//' < '//c_hinmax2//'m'
         enddo

         write(warning,*) ' '
         call add_warning(warning)

      end subroutine colpkg_init_itd_hist

!=======================================================================
!
! Initialize the vertical profile of ice salinity and melting temperature.
!
! authors: C. M. Bitz, UW
!          William H. Lipscomb, LANL

      subroutine colpkg_init_thermo(nilyr, sprofile)

      use ice_colpkg_shared, only: saltmax, ktherm, heat_capacity, &
          min_salin
      use ice_constants_colpkg, only: p5, c0, c1, c2, pi
      use ice_therm_shared, only: l_brine

      integer (kind=int_kind), intent(in) :: &
         nilyr                            ! number of ice layers

      real (kind=dbl_kind), dimension(:), intent(out) :: &
         sprofile                         ! vertical salinity profile

      real (kind=dbl_kind), parameter :: &
         nsal    = 0.407_dbl_kind, &
         msal    = 0.573_dbl_kind

      integer (kind=int_kind) :: k        ! ice layer index
      real (kind=dbl_kind)    :: zn       ! normalized ice thickness

      !-----------------------------------------------------------------
      ! Determine l_brine based on saltmax.
      ! Thermodynamic solver will not converge if l_brine is true and
      !  saltmax is close to zero.
      ! Set l_brine to false for zero layer thermodynamics
      !-----------------------------------------------------------------

      heat_capacity = .true.      
      if (ktherm == 0) heat_capacity = .false. ! 0-layer thermodynamics

      l_brine = .false.
      if (saltmax > min_salin .and. heat_capacity) l_brine = .true.

      !-----------------------------------------------------------------
      ! Prescibe vertical profile of salinity and melting temperature.
      ! Note this profile is only used for BL99 thermodynamics.
      !-----------------------------------------------------------------

      if (l_brine) then
         do k = 1, nilyr
            zn = (real(k,kind=dbl_kind)-p5) /  &
                  real(nilyr,kind=dbl_kind)
            sprofile(k)=(saltmax/c2)*(c1-cos(pi*zn**(nsal/(msal+zn))))
            sprofile(k) = max(sprofile(k), min_salin)
         enddo ! k
         sprofile(nilyr+1) = saltmax

      else ! .not. l_brine
         do k = 1, nilyr+1
            sprofile(k) = c0
         enddo
      endif ! l_brine

      end subroutine colpkg_init_thermo

!=======================================================================
! Initial salinity profile
!
! authors: C. M. Bitz, UW
!          William H. Lipscomb, LANL

      function colpkg_salinity_profile(zn) result(salinity)

        use ice_colpkg_shared, only: saltmax
        use ice_constants_colpkg, only: c1, c2, pi

        real(kind=dbl_kind), intent(in) :: &
             zn ! depth

        real(kind=dbl_kind) :: &
             salinity ! initial salinity profile

        real (kind=dbl_kind), parameter :: &
             nsal    = 0.407_dbl_kind, &
             msal    = 0.573_dbl_kind

        salinity = (saltmax/c2)*(c1-cos(pi*zn**(nsal/(msal+zn))))

      end function colpkg_salinity_profile

!=======================================================================
! Compute orbital parameters for the specified date.
!
! author:  Bruce P. Briegleb, NCAR 

      subroutine colpkg_init_orbit(l_stop, stop_label)

      use ice_constants_colpkg, only: iyear_AD, eccen, obliqr, lambm0, &
         mvelpp, obliq, mvelp, decln, eccf, log_print

#ifdef CCSMCOUPLED
      use shr_orb_mod, only: shr_orb_params
#else
      use ice_orbital, only: shr_orb_params
#endif
      !  real (kind=8)    :: obliq_R8
      !  real (kind=8)    :: mvelp_R8
      !  real (kind=8)    :: eccen_R8
      !  real (kind=8)    :: mvelpp_R8
      !  real (kind=8)    :: lambm0_R8
      !  real (kind=8)    :: obliqr_R8
      !  real (kind=8)    :: decln_R8
      !  real (kind=8)    :: eccf_R8


      logical (kind=log_kind), intent(out) :: &
         l_stop          ! if true, abort the model

      character (len=*), intent(out) :: stop_label

      ! eccen_R8 = real(eccen, kind = 8)
      ! mvelpp_R8 = real(mvelpp, kind = 8)
      ! mvelp_R8 = real(mvelp, kind = 8)
      ! obliq_R8  = real(obliq, kind =8)
      ! lambm0_R8  = real(lambm0, kind =8)
      ! obliqr_R8  = real(obliqr, kind =8)
      !decln_R8  = real(decln, kind = R8KIND)
      !eccf_R8  = real(eccf, kind = R8KIND)
      l_stop = .false.      ! initialized for CCSMCOUPLED
      stop_label = ''       ! initialized for CCSMCOUPLED
      iyear_AD  = 1950
      log_print = .false.   ! if true, write out orbital parameters

#ifdef CCSMCOUPLED
      call shr_orb_params( iyear_AD, eccen , obliq , mvelp , &
                           obliqr  , lambm0, mvelpp, log_print)
      ! eccen = real(eccen_R8, kind=dbl_kind) 
      ! obliq = real(obliq_R8, kind=dbl_kind)
      ! obliqr = real(obliqr_R8, kind=dbl_kind)
      ! mvelp = real(mvelp_R8, kind=dbl_kind)
      ! mvelpp = real(mvelpp_R8, kind=dbl_kind)
      ! lambm0 = real(lambm0_R8, kind=dbl_kind)
#else
      call shr_orb_params( iyear_AD, eccen , obliq , mvelp    , &
                           obliqr  , lambm0, mvelpp, log_print, &
                           l_stop, stop_label)
#endif

      end subroutine colpkg_init_orbit
 
!=======================================================================

      subroutine colpkg_init_trcr(Tair_R4,     Tf_R4,       &
                                  Sprofile_R4, Tprofile_R4, &
                                  Tsfc_R4,               &
                                  nilyr,    nslyr,    &
                                  qin_R4,      qsn_R4)

      use ice_colpkg_shared, only: calc_Tsfc
      use ice_constants_colpkg, only: Tsmelt, Tffresh, p5, cp_ice, cp_ocn, &
          Lfresh, rhoi, rhos, c0, c1
      use ice_mushy_physics, only: enthalpy_mush

      integer (kind=int_kind), intent(in) :: &
         nilyr, &    ! number of ice layers
         nslyr       ! number of snow layers

      real (kind=real_kind), intent(in) :: &
         Tair_R4, &     ! air temperature (C)
         Tf_R4          ! freezing temperature (C)

      real (kind=real_kind), dimension(:), intent(in) :: &
         Sprofile_R4, & ! vertical salinity profile (ppt)
         Tprofile_R4    ! vertical temperature profile (C)

      real (kind=real_kind), intent(out) :: &
         Tsfc_R4        ! surface temperature (C)

      real (kind=real_kind), dimension(:), intent(out) :: &
         qin_R4, &      ! ice enthalpy profile (J/m3)
         qsn_R4         ! snow enthalpy profile (J/m3)

      ! local variables

      integer (kind=int_kind) :: k

      real (kind=dbl_kind) :: & 
         Tair, &     ! air temperature (C)
         Tf          ! freezing temperature (C)
      real (kind=dbl_kind), dimension(:), allocatable :: &
         Sprofile, & ! vertical salinity profile (ppt)
         Tprofile    ! vertical temperature profile (C)
      real (kind=dbl_kind) :: &
         Tsfc        ! surface temperature (C)
      real (kind=dbl_kind), dimension(:), allocatable :: &
         qin, &      ! ice enthalpy profile (J/m3)
         qsn         ! snow enthalpy profile (J/m3)

      real (kind=dbl_kind) :: &
         slope, Ti

      Tair = real(Tair_R4, kind = dbl_kind)
      Tf = real(Tf_R4, kind = dbl_kind)          
      Tsfc = real(Tsfc_R4, kind = dbl_kind)        

     Sprofile = real(Sprofile_R4, kind = dbl_kind)
     Tprofile = real(Tprofile_R4, kind = dbl_kind)    
     qin = real(qin_R4, kind = dbl_kind)
     qsn = real(qsn_R4, kind = dbl_kind)         
      
            ! surface temperature
            Tsfc = Tf ! default
            if (calc_Tsfc) Tsfc = min(Tsmelt, Tair - Tffresh) ! deg C

            if (heat_capacity) then

               ! ice enthalpy
               do k = 1, nilyr
                  ! assume linear temp profile and compute enthalpy
                  slope = Tf - Tsfc
                  Ti = Tsfc + slope*(real(k,kind=dbl_kind)-p5) &
                                    /real(nilyr,kind=dbl_kind)
                  if (ktherm == 2) then
                     qin(k) = enthalpy_mush(Ti, Sprofile(k))
                  else
                     qin(k) = -(rhoi * (cp_ice*(Tprofile(k)-Ti) &
                         + Lfresh*(c1-Tprofile(k)/Ti) - cp_ocn*Tprofile(k)))
                  endif
               enddo               ! nilyr

               ! snow enthalpy
               do k = 1, nslyr
                  Ti = min(c0, Tsfc)
                  qsn(k) = -rhos*(Lfresh - cp_ice*Ti)
               enddo               ! nslyr

            else  ! one layer with zero heat capacity

               ! ice energy
               qin(1) = -rhoi * Lfresh 

               ! snow energy
               qsn(1) = -rhos * Lfresh 

            endif               ! heat_capacity


      Tsfc_R4 = real(Tsfc, kind = real_kind)        
      qin_R4 = real(qin, kind = real_kind)
      qsn_R4 = real(qsn, kind = real_kind) 

      deallocate(Sprofile)
      deallocate(Tprofile)
      deallocate(qin)
      deallocate(qsn)
      end subroutine colpkg_init_trcr

!=======================================================================

      subroutine colpkg_init_bgc(dt_R4, ncat, nblyr, nilyr, ntrcr_o, cgrid_R4, igrid_R4, &
         restart_bgc, ntrcr, nbtrcr, sicen_R4, trcrn_R4, &
         sss_R4, nit_R4, amm_R4, sil_R4, dmsp_R4, dms_R4, algalN_R4, &
         doc_R4, don_R4, dic_R4, fed_R4, fep_R4, zaeros_R4, hum_R4,  &
         ocean_bio_all_R4, &
         max_algae, max_doc, max_dic, max_don,  max_fe, max_nbtrcr, max_aero, &
         l_stop, stop_label)

      use ice_constants_colpkg, only: c0, c1, c2, p1, p15, p5
      use ice_zbgc_shared, only: R_S2N, zbgc_frac_init, zbgc_init_frac, remap_zbgc

      ! column package includes
      use ice_colpkg_tracers, only: nt_fbri, nt_bgc_S, nt_sice, nt_zbgc_frac, &
         bio_index_o,  bio_index  
      use ice_colpkg_shared, only: solve_zsal, ktherm, hs_ssl,  &
         skl_bgc, scale_bgc, grid_o_t,  fe_data_type, &
         R_C2N, R_chl2N

      real (kind=real_kind), intent(in) :: &
         dt_R4        ! time step

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         nilyr , & ! number of ice layers
         nblyr , & ! number of bio layers
         ntrcr_o, & ! number of tracers not including bgc
         ntrcr , & ! number of tracers in use
         nbtrcr, & ! number of bio tracers in use
         max_algae, &
         max_doc, &
         max_dic, &
         max_don, &
         max_fe, &
         max_nbtrcr, &
         max_aero
 
      logical (kind=log_kind), intent(in) :: & 
         restart_bgc ! if .true., read bgc restart file



      real (kind=real_kind), dimension (nblyr+1), intent(inout) :: &
         igrid_R4     ! biology vertical interface points
 
      real (kind=real_kind), dimension (nilyr+1), intent(inout) :: &
         cgrid_R4     ! CICE vertical coordinate   

      real (kind=real_kind), dimension(nilyr, ncat), intent(in) :: &
         sicen_R4     ! salinity on the cice grid

      real (kind=real_kind), dimension (:,:), intent(inout) :: &
         trcrn_R4     ! subset of tracer array (only bgc) 

      real (kind=real_kind), intent(in) :: &
         sss_R4       ! sea surface salinity (ppt)

      real (kind=real_kind), intent(inout) :: &
         nit_R4   , & ! ocean nitrate (mmol/m^3)          
         amm_R4   , & ! ammonia/um (mmol/m^3)
         sil_R4   , & ! silicate (mmol/m^3)
         dmsp_R4  , & ! dmsp (mmol/m^3)
         dms_R4   , & ! dms (mmol/m^3)
         hum_R4       ! hum (mmol/m^3)

      real (kind=real_kind), dimension (max_algae), intent(inout) :: &
         algalN_R4    ! ocean algal nitrogen (mmol/m^3) (diatoms, pico, phaeocystis)

      real (kind=real_kind), dimension (max_doc), intent(inout) :: &
         doc_R4       ! ocean doc (mmol/m^3)  (proteins, EPS, lipid)

      real (kind=real_kind), dimension (max_don), intent(inout) :: &
         don_R4       ! ocean don (mmol/m^3) 

      real (kind=real_kind), dimension (max_dic), intent(inout) :: &
         dic_R4       ! ocean dic (mmol/m^3) 

      real (kind=real_kind), dimension (max_fe), intent(inout) :: &
         fed_R4, fep_R4  ! ocean disolved and particulate fe (nM) 

      real (kind=real_kind), dimension (max_aero), intent(inout) :: &
         zaeros_R4    ! ocean aerosols (mmol/m^3) 

      real (kind=real_kind), dimension (:), intent(inout) :: &
         ocean_bio_all_R4   ! fixed order, all values even for tracers false
      logical (kind=log_kind), intent(inout) :: &
         l_stop            ! if true, print diagnostics and abort on return

      character (len=*), intent(inout) :: stop_label

      ! local variables

      integer (kind=int_kind) :: &
         k     , & ! vertical index 
         n     , & ! category index 
         mm    , & ! bio tracer index
         ki    , & ! loop index
         ks    , & ! 
         ntrcr_bgc

      real (kind=dbl_kind), dimension (ntrcr+2) :: & 
         trtmp     ! temporary, remapped tracers   
      
      real (kind=dbl_kind), dimension (nblyr+1) :: &
         zspace    ! vertical grid spacing

      real (kind=dbl_kind) :: & 
         dvssl , & ! volume of snow surface layer (m)
         dvint , & ! volume of snow interior      (m)
         nit_dum, &!
         sil_dum

      real (kind=dbl_kind) :: &
         dt        ! time step
      real (kind=dbl_kind), dimension (nblyr+1)  :: &
         igrid     ! biology vertical interface points
      real (kind=dbl_kind), dimension (nilyr+1)  :: &
         cgrid     ! CICE vertical coordinate   
      real (kind=dbl_kind), dimension(nilyr, ncat) :: &
         sicen     ! salinity on the cice grid
      real (kind=dbl_kind), dimension (:,:), allocatable  :: &
         trcrn     ! subset of tracer array (only bgc) 
      real (kind=dbl_kind) :: &
         sss       ! sea surface salinity (ppt)
      real (kind=dbl_kind) :: &
         nit,  &
         amm,  &
         sil,  &
         dmsp,  &
         dms,  &
         hum       ! hum (mmol/m^3)
      real (kind=dbl_kind), dimension (max_algae) :: &
         algalN    ! ocean algal nitrogen (mmol/m^3) (diatoms :: &
      real (kind=dbl_kind), dimension (max_doc) :: &
         doc       ! ocean doc (mmol/m^3)  (proteins :: &
      real (kind=dbl_kind), dimension (max_don) :: &
         don       ! ocean don (mmol/m^3) 
      real (kind=dbl_kind), dimension (max_dic) :: &
         dic       ! ocean dic (mmol/m^3) 
      real (kind=dbl_kind), dimension (max_fe) :: &
         fed, fep
      real (kind=dbl_kind), dimension (max_aero) :: &
         zaeros    ! ocean aerosols (mmol/m^3) 
      real (kind=dbl_kind), dimension (:), allocatable :: &
         ocean_bio_all   ! fixed order :: &

     ocean_bio_all =real(ocean_bio_all_R4, kind = dbl_kind)
     trcrn =real(trcrn_R4, kind = dbl_kind)

      dt = real(dt_R4, kind = dbl_kind)         
      igrid = real(igrid_R4, kind = dbl_kind)      
      cgrid = real(cgrid_R4, kind = dbl_kind)      
      sicen = real(sicen_R4, kind = dbl_kind)      
      sss = real(sss_R4, kind = dbl_kind)        
      nit = real(nit_R4, kind = dbl_kind)     
      amm = real(amm_R4, kind = dbl_kind)     
      sil = real(sil_R4, kind = dbl_kind)     
      dmsp = real(dmsp_R4, kind = dbl_kind)    
      dms = real(dms_R4, kind = dbl_kind)     
      hum = real(hum_R4, kind = dbl_kind)        
      algalN = real(algalN_R4, kind = dbl_kind)     
      doc = real(doc_R4, kind = dbl_kind)        
      don = real(don_R4, kind = dbl_kind)        
      dic = real(dic_R4, kind = dbl_kind)        
      fed = real(fed_R4, kind = dbl_kind)  
      fep = real(fep_R4, kind = dbl_kind) 
      zaeros = real(zaeros_R4, kind = dbl_kind)     
      ! ocean_bio_all = real(ocean_bio_all_R4, kind = dbl_kind)   


      zspace(:)       = c1/real(nblyr,kind=dbl_kind)
      zspace(1)       = p5*zspace(1)
      zspace(nblyr+1) = p5*zspace(nblyr+1)
      ntrcr_bgc       = ntrcr-ntrcr_o

      call colpkg_init_OceanConcArray_double(max_nbtrcr,                &
                                 max_algae, max_don,  max_doc,   &
                                 max_dic,   max_aero, max_fe,    &
                                 nit,       amm,      sil,       &
                                 dmsp,      dms,      algalN,    &
                                 doc,       don,      dic,       &  
                                 fed,       fep,      zaeros,    &
                                 ocean_bio_all,       hum)

      if (.not. restart_bgc) then  ! not restarting

      !-----------------------------------------------------------------------------   
      !     Skeletal Layer Model
      !  All bgc tracers are Bulk quantities in units of mmol or mg per m^3
      !  The skeletal layer model assumes a constant 
      !  layer depth (sk_l) and porosity (phi_sk)
      !-----------------------------------------------------------------------------   
         if (skl_bgc) then
       
            do  n = 1,ncat
            do mm = 1,nbtrcr
               ! bulk concentration (mmol or mg per m^3, or 10^-3 mmol/m^3)
               trcrn(bio_index(mm)-ntrcr_o, n) = ocean_bio_all(bio_index_o(mm))
            enddo       ! nbtrcr
            enddo       ! n 

      !-----------------------------------------------------------------------------   
      !    zbgc Model
      !  All bgc tracers are Bulk quantities in units of mmol or mg per m^3
      !  The vertical layer model uses prognosed porosity and layer depth
      !-----------------------------------------------------------------------------   

         else   ! not skl_bgc

            if (scale_bgc .and. solve_zsal) then ! bulk concentration (mmol or mg per m^3)
               do n = 1,ncat
               do mm = 1,nbtrcr
                  do k = 2, nblyr
                     trcrn(bio_index(mm)+k-1-ntrcr_o,n) = &
                          (p5*(trcrn(nt_bgc_S+k-1-ntrcr_o,n)+ trcrn(nt_bgc_S+k-2-ntrcr_o,n)) &
                         / sss*ocean_bio_all(bio_index_o(mm))) 
                  enddo  !k
                  trcrn(nt_zbgc_frac-1+mm-ntrcr_o,n) = zbgc_frac_init(mm)
                  trcrn(bio_index(mm)-ntrcr_o,n) = (trcrn(nt_bgc_S-ntrcr_o,n) &
                                         / sss*ocean_bio_all(bio_index_o(mm))) 
                  trcrn(bio_index(mm)+nblyr-ntrcr_o,n) = (trcrn(nt_bgc_S+nblyr-1-ntrcr_o,n) &
                                               / sss*ocean_bio_all(bio_index_o(mm)))
                  trcrn(bio_index(mm)+nblyr+1-ntrcr_o:bio_index(mm)+nblyr+2-ntrcr_o,n) = c0 ! snow
               enddo ! mm
               enddo ! n 
    
            elseif (scale_bgc .and. ktherm == 2) then
               trtmp(:) = c0
               do n = 1,ncat     
                  call remap_zbgc(nilyr,            nilyr,    &
                                  1,                          &
                                  sicen(:,n),       trtmp,    &
                                  0,                nblyr+1,  &
                                  c1,               c1,       &
                                  cgrid(2:nilyr+1),           &
                                  igrid(1:nblyr+1),           &
                                  sicen(1,n),                 &
                                  l_stop,           stop_label)
                  if (l_stop) return

                  do mm = 1,nbtrcr
                  do k = 1, nblyr + 1            
                     trcrn(bio_index(mm)+k-1-ntrcr_o,n) =   &
                          (trtmp(k)/sss*ocean_bio_all(bio_index_o(mm)))
                     trcrn(bio_index(mm)+nblyr+1-ntrcr_o:bio_index(mm)+nblyr+2-ntrcr_o,n) = c0 ! snow
                  enddo  ! k
                  enddo  ! mm
               enddo     ! n 

            elseif (nbtrcr > 0 .and. nt_fbri > 0) then ! not scale_bgc         
     
               do n = 1,ncat
               do mm = 1,nbtrcr
               do k = 1, nblyr+1
                  trcrn(bio_index(mm)+k-1-ntrcr_o,n) = ocean_bio_all(bio_index_o(mm)) &
                                             * zbgc_init_frac(mm) 
                  trcrn(bio_index(mm)+nblyr+1-ntrcr_o:bio_index(mm)+nblyr+2-ntrcr_o,n) = c0 ! snow
               enddo    ! k
               trcrn(nt_zbgc_frac-1+mm-ntrcr_o,n) = zbgc_frac_init(mm)
               enddo    ! mm
               enddo    ! n 
              
            endif  ! scale_bgc
         endif     ! skl_bgc
      endif        ! restart
    
      igrid_R4 = real(igrid, kind = real_kind)      
      cgrid_R4 = real(cgrid, kind = real_kind)      
      trcrn_R4 = real(trcrn, kind = real_kind)      
      nit_R4 = real(nit, kind = real_kind)     
      amm_R4 = real(amm, kind = real_kind)     
      sil_R4 = real(sil, kind = real_kind)     
      dmsp_R4 = real(dmsp, kind = real_kind)    
      dms_R4 = real(dms, kind = real_kind)     
      hum_R4 = real(hum, kind = real_kind)        
      algalN_R4 = real(algalN, kind = real_kind)     
      doc_R4 = real(doc, kind = real_kind)        
      don_R4 = real(don, kind = real_kind)        
      dic_R4 = real(dic, kind = real_kind)        
      fed_R4 = real(fed, kind = real_kind)  
      fep_R4 = real(fep, kind = real_kind) 
      zaeros_R4 = real(zaeros, kind = real_kind)     
      ocean_bio_all_R4 = real(ocean_bio_all, kind = real_kind)   

      deallocate(ocean_bio_all)
      deallocate(trcrn)
      end subroutine colpkg_init_bgc

   
      subroutine colpkg_init_bgc_double(dt, ncat, nblyr, nilyr, ntrcr_o, cgrid, igrid, &
         restart_bgc, ntrcr, nbtrcr, sicen, trcrn, &
         sss, nit, amm, sil, dmsp, dms, algalN, &
         doc, don, dic, fed, fep, zaeros, hum,  &
         ocean_bio_all, &
         max_algae, max_doc, max_dic, max_don,  max_fe, max_nbtrcr, max_aero, &
         l_stop, stop_label)

      use ice_constants_colpkg, only: c0, c1, c2, p1, p15, p5
      use ice_zbgc_shared, only: R_S2N, zbgc_frac_init, zbgc_init_frac, remap_zbgc

      ! column package includes
      use ice_colpkg_tracers, only: nt_fbri, nt_bgc_S, nt_sice, nt_zbgc_frac, &
         bio_index_o,  bio_index  
      use ice_colpkg_shared, only: solve_zsal, ktherm, hs_ssl,  &
         skl_bgc, scale_bgc, grid_o_t,  fe_data_type, &
         R_C2N, R_chl2N

      real (kind=dbl_kind), intent(in) :: &
         dt        ! time step

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         nilyr , & ! number of ice layers
         nblyr , & ! number of bio layers
         ntrcr_o, & ! number of tracers not including bgc
         ntrcr , & ! number of tracers in use
         nbtrcr, & ! number of bio tracers in use
         max_algae, &
         max_doc, &
         max_dic, &
         max_don, &
         max_fe, &
         max_nbtrcr, &
         max_aero
 
      real (kind=dbl_kind), dimension (nblyr+1), intent(inout) :: &
         igrid     ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), intent(inout) :: &
         cgrid     ! CICE vertical coordinate   

      logical (kind=log_kind), intent(in) :: & 
         restart_bgc ! if .true., read bgc restart file

      real (kind=dbl_kind), dimension(nilyr, ncat), intent(in) :: &
         sicen     ! salinity on the cice grid

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn     ! subset of tracer array (only bgc) 

      real (kind=dbl_kind), intent(in) :: &
         sss       ! sea surface salinity (ppt)

      real (kind=dbl_kind), intent(inout) :: &
         nit   , & ! ocean nitrate (mmol/m^3)          
         amm   , & ! ammonia/um (mmol/m^3)
         sil   , & ! silicate (mmol/m^3)
         dmsp  , & ! dmsp (mmol/m^3)
         dms   , & ! dms (mmol/m^3)
         hum       ! hum (mmol/m^3)

      real (kind=dbl_kind), dimension (max_algae), intent(inout) :: &
         algalN    ! ocean algal nitrogen (mmol/m^3) (diatoms, pico, phaeocystis)

      real (kind=dbl_kind), dimension (max_doc), intent(inout) :: &
         doc       ! ocean doc (mmol/m^3)  (proteins, EPS, lipid)

      real (kind=dbl_kind), dimension (max_don), intent(inout) :: &
         don       ! ocean don (mmol/m^3) 

      real (kind=dbl_kind), dimension (max_dic), intent(inout) :: &
         dic       ! ocean dic (mmol/m^3) 

      real (kind=dbl_kind), dimension (max_fe), intent(inout) :: &
         fed, fep  ! ocean disolved and particulate fe (nM) 

      real (kind=dbl_kind), dimension (max_aero), intent(inout) :: &
         zaeros    ! ocean aerosols (mmol/m^3) 

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         ocean_bio_all   ! fixed order, all values even for tracers false

      logical (kind=log_kind), intent(inout) :: &
         l_stop            ! if true, print diagnostics and abort on return

      character (len=*), intent(inout) :: stop_label

      ! local variables

      integer (kind=int_kind) :: &
         k     , & ! vertical index 
         n     , & ! category index 
         mm    , & ! bio tracer index
         ki    , & ! loop index
         ks    , & ! 
         ntrcr_bgc

      real (kind=dbl_kind), dimension (ntrcr+2) :: & 
         trtmp     ! temporary, remapped tracers   
      
      real (kind=dbl_kind), dimension (nblyr+1) :: &
         zspace    ! vertical grid spacing

      real (kind=dbl_kind) :: & 
         dvssl , & ! volume of snow surface layer (m)
         dvint , & ! volume of snow interior      (m)
         nit_dum, &!
         sil_dum

      zspace(:)       = c1/real(nblyr,kind=dbl_kind)
      zspace(1)       = p5*zspace(1)
      zspace(nblyr+1) = p5*zspace(nblyr+1)
      ntrcr_bgc       = ntrcr-ntrcr_o

      call colpkg_init_OceanConcArray_double(max_nbtrcr,                &
                                 max_algae, max_don,  max_doc,   &
                                 max_dic,   max_aero, max_fe,    &
                                 nit,       amm,      sil,       &
                                 dmsp,      dms,      algalN,    &
                                 doc,       don,      dic,       &  
                                 fed,       fep,      zaeros,    &
                                 ocean_bio_all,       hum)

      if (.not. restart_bgc) then  ! not restarting

      !-----------------------------------------------------------------------------   
      !     Skeletal Layer Model
      !  All bgc tracers are Bulk quantities in units of mmol or mg per m^3
      !  The skeletal layer model assumes a constant 
      !  layer depth (sk_l) and porosity (phi_sk)
      !-----------------------------------------------------------------------------   
         if (skl_bgc) then
       
            do  n = 1,ncat
            do mm = 1,nbtrcr
               ! bulk concentration (mmol or mg per m^3, or 10^-3 mmol/m^3)
               trcrn(bio_index(mm)-ntrcr_o, n) = ocean_bio_all(bio_index_o(mm))
            enddo       ! nbtrcr
            enddo       ! n 

      !-----------------------------------------------------------------------------   
      !    zbgc Model
      !  All bgc tracers are Bulk quantities in units of mmol or mg per m^3
      !  The vertical layer model uses prognosed porosity and layer depth
      !-----------------------------------------------------------------------------   

         else   ! not skl_bgc

            if (scale_bgc .and. solve_zsal) then ! bulk concentration (mmol or mg per m^3)
               do n = 1,ncat
               do mm = 1,nbtrcr
                  do k = 2, nblyr
                     trcrn(bio_index(mm)+k-1-ntrcr_o,n) = &
                          (p5*(trcrn(nt_bgc_S+k-1-ntrcr_o,n)+ trcrn(nt_bgc_S+k-2-ntrcr_o,n)) &
                         / sss*ocean_bio_all(bio_index_o(mm))) 
                  enddo  !k
                  trcrn(nt_zbgc_frac-1+mm-ntrcr_o,n) = zbgc_frac_init(mm)
                  trcrn(bio_index(mm)-ntrcr_o,n) = (trcrn(nt_bgc_S-ntrcr_o,n) &
                                         / sss*ocean_bio_all(bio_index_o(mm))) 
                  trcrn(bio_index(mm)+nblyr-ntrcr_o,n) = (trcrn(nt_bgc_S+nblyr-1-ntrcr_o,n) &
                                               / sss*ocean_bio_all(bio_index_o(mm)))
                  trcrn(bio_index(mm)+nblyr+1-ntrcr_o:bio_index(mm)+nblyr+2-ntrcr_o,n) = c0 ! snow
               enddo ! mm
               enddo ! n 
    
            elseif (scale_bgc .and. ktherm == 2) then
               trtmp(:) = c0
               do n = 1,ncat     
                  call remap_zbgc(nilyr,            nilyr,    &
                                  1,                          &
                                  sicen(:,n),       trtmp,    &
                                  0,                nblyr+1,  &
                                  c1,               c1,       &
                                  cgrid(2:nilyr+1),           &
                                  igrid(1:nblyr+1),           &
                                  sicen(1,n),                 &
                                  l_stop,           stop_label)
                  if (l_stop) return

                  do mm = 1,nbtrcr
                  do k = 1, nblyr + 1            
                     trcrn(bio_index(mm)+k-1-ntrcr_o,n) =   &
                          (trtmp(k)/sss*ocean_bio_all(bio_index_o(mm)))
                     trcrn(bio_index(mm)+nblyr+1-ntrcr_o:bio_index(mm)+nblyr+2-ntrcr_o,n) = c0 ! snow
                  enddo  ! k
                  enddo  ! mm
               enddo     ! n 

            elseif (nbtrcr > 0 .and. nt_fbri > 0) then ! not scale_bgc         
     
               do n = 1,ncat
               do mm = 1,nbtrcr
               do k = 1, nblyr+1
                  trcrn(bio_index(mm)+k-1-ntrcr_o,n) = ocean_bio_all(bio_index_o(mm)) &
                                             * zbgc_init_frac(mm) 
                  trcrn(bio_index(mm)+nblyr+1-ntrcr_o:bio_index(mm)+nblyr+2-ntrcr_o,n) = c0 ! snow
               enddo    ! k
               trcrn(nt_zbgc_frac-1+mm-ntrcr_o,n) = zbgc_frac_init(mm)
               enddo    ! mm
               enddo    ! n 
              
            endif  ! scale_bgc
         endif     ! skl_bgc
      endif        ! restart

      end subroutine colpkg_init_bgc_double

!=======================================================================
!=======================================================================

      subroutine colpkg_init_zbgc (nblyr, nilyr, nslyr, &
                 n_algae, n_zaero, n_doc, n_dic, n_don, n_fed, n_fep, &
                 trcr_base_R4, trcr_depend, n_trcr_strata, nt_strata, nbtrcr_sw, &
                 tr_brine, nt_fbri, ntrcr, nbtrcr, nt_bgc_Nit, nt_bgc_Am, &
                 nt_bgc_Sil, nt_bgc_DMS, nt_bgc_PON, nt_bgc_S, nt_bgc_N, &
                 nt_bgc_C, nt_bgc_chl, nt_bgc_DOC, nt_bgc_DON, nt_bgc_DIC, & 
                 nt_zaero, nt_bgc_DMSPp, nt_bgc_DMSPd, nt_bgc_Fed, nt_bgc_Fep, &
                 nt_zbgc_frac, tr_bgc_Nit, tr_bgc_Am, tr_bgc_Sil, tr_bgc_DMS, &
                 tr_bgc_PON, tr_bgc_S, tr_bgc_N, tr_bgc_C, tr_bgc_chl, &
                 tr_bgc_DON, tr_bgc_Fe, tr_zaero, nlt_zaero_sw, nlt_chl_sw, &
                 nlt_bgc_N, nlt_bgc_Nit, nlt_bgc_Am, nlt_bgc_Sil, &
                 nlt_bgc_DMS, nlt_bgc_DMSPp, nlt_bgc_DMSPd, &
                 nlt_bgc_C, nlt_bgc_chl, nlt_bgc_DIC, nlt_bgc_DOC, &
                 nlt_bgc_PON, nlt_bgc_DON, nlt_bgc_Fed, nlt_bgc_Fep, &
                 nlt_zaero, &
                 nt_bgc_hum, nlt_bgc_hum, tr_bgc_hum, solve_zsal, &
                 skl_bgc, z_tracers, dEdd_algae, solve_zbgc, &
                 frazil_scav, initbio_frac, bio_index_o, bio_index, ntrcr_o, &
                 max_algae, max_doc, max_dic, max_don, max_fe, &
                 ratio_Si2N_diatoms, ratio_Si2N_sp, ratio_Si2N_phaeo, &
                 ratio_S2N_diatoms, ratio_S2N_sp, ratio_S2N_phaeo, &
                 ratio_Fe2C_diatoms, ratio_Fe2C_sp, ratio_Fe2C_phaeo, &
                 ratio_Fe2N_diatoms, ratio_Fe2N_sp, ratio_Fe2N_phaeo, &
                 ratio_Fe2DON, ratio_Fe2DOC_s,  ratio_Fe2DOC_l, & 
                 chlabs_diatoms, chlabs_sp, chlabs_phaeo, &    
                 alpha2max_low_diatoms, alpha2max_low_sp, alpha2max_low_phaeo, &  
                 beta2max_diatoms, beta2max_sp, beta2max_phaeo, &    
                 mu_max_diatoms, mu_max_sp, mu_max_phaeo, &      
                 grow_Tdep_diatoms, grow_Tdep_sp, grow_Tdep_phaeo, &      
                 fr_graze_diatoms, fr_graze_sp, fr_graze_phaeo, &    
                 mort_pre_diatoms, mort_pre_sp, mort_pre_phaeo, &        
                 mort_Tdep_diatoms, mort_Tdep_sp, mort_Tdep_phaeo, &
                 k_exude_diatoms, k_exude_sp, k_exude_phaeo, &   
                 K_Nit_diatoms, K_Nit_sp, K_Nit_phaeo, &     
                 K_Am_diatoms, K_Am_sp, K_Am_phaeo, &     
                 K_Sil_diatoms, K_Sil_sp, K_Sil_phaeo, &     
                 K_Fe_diatoms, K_Fe_sp, K_Fe_phaeo, &  
                 f_don_protein, kn_bac_protein, & 
                 f_don_Am_protein ,f_doc_s, f_doc_l, &
                 f_exude_s, f_exude_l, k_bac_s,  k_bac_l, &
                 algaltype_diatoms, algaltype_sp, algaltype_phaeo, &
                 doctype_s, doctype_l, dictype_1, dontype_protein, &
                 fedtype_1, feptype_1, zaerotype_bc1, zaerotype_bc2, &
                 zaerotype_dust1, zaerotype_dust2, zaerotype_dust3, &
                 zaerotype_dust4, &
                 ratio_C2N_diatoms, ratio_C2N_sp, ratio_C2N_phaeo, &
                 ratio_chl2N_diatoms, ratio_chl2N_sp, ratio_chl2N_phaeo, &
                 F_abs_chl_diatoms, F_abs_chl_sp, F_abs_chl_phaeo, &
                 ratio_C2N_proteins, &
                 nitratetype, ammoniumtype, dmspptype, dmspdtype, &
                 silicatetype, humtype, tau_min, tau_max)
                    
      use ice_constants_colpkg, only: c1, p5, c0, c2

      use ice_colpkg_shared, only: &
         algaltype, doctype, dictype, dontype, fedtype, feptype, zaerotype, &
         R_C2N, R_chl2N, F_abs_chl, R_C2N_DON

      use ice_zbgc_shared, only: zbgc_init_frac, &
         bgc_tracer_type, zbgc_frac_init, &
         tau_ret, tau_rel, R_Si2N, R_S2N, R_Fe2C, &
         R_Fe2N, R_Fe2DON, R_Fe2DOC, &
         chlabs, alpha2max_low, beta2max, &
         mu_max, grow_Tdep, fr_graze, &
         mort_pre, mort_Tdep, k_exude, &
         K_Nit, K_Am, K_Sil, K_Fe, &
         f_don, kn_bac, f_don_Am, &
         f_doc, f_exude, k_bac
         

      integer (kind=int_kind), intent(in) :: &
         nblyr     , & ! number of bio/brine layers per category 
         nilyr     , & ! number of ice layers per category
         nslyr     , & ! number of snow layers per category
         n_zaero   , & ! number of z aerosols in use 
         n_algae   , & ! number of algae in use 
         n_doc     , & ! number of DOC pools in use
         n_dic     , & ! number of DIC pools in use
         n_don     , & ! number of DON pools in use
         n_fed     , & ! number of Fe  pools in use dissolved Fe
         n_fep     , & ! number of Fe  pools in use particulate Fe
         max_algae , &
         max_doc   , &
         max_dic   , &
         max_don   , &
         max_fe

      integer (kind=int_kind), intent(inout) :: &
         ntrcr_o,     & ! number of non-bio tracers in use
         ntrcr,       & ! number of tracers
         nbtrcr,      & ! number of bgc tracers in use
         nbtrcr_sw      ! size of shorwave tracer vector

      integer (kind=int_kind), dimension (:), intent(inout) :: &
         trcr_depend   ! = 0 for ice area tracers
                       ! = 1 for ice volume tracers
                       ! = 2 for snow volume tracers

      integer (kind=int_kind), dimension (:), intent(inout) :: &
         n_trcr_strata ! number of underlying tracer layers

      integer (kind=int_kind), dimension (:,:), intent(inout) :: &
         nt_strata     ! indices of underlying tracer layers

      real (kind=real_kind), dimension (:,:), intent(inout) :: &
         trcr_base_R4     ! = 0 or 1 depending on tracer dependency
                       ! argument 2:  (1) aice, (2) vice, (3) vsno

      logical (kind=log_kind), intent(in) :: &
         tr_brine,       & ! if .true., brine height differs from ice thickness
         tr_bgc_S,       & ! if .true., use zsalinity
         tr_zaero,       & ! if .true., black carbon is tracers  (n_zaero)
         tr_bgc_Nit,     & ! if .true. Nitrate tracer in ice 
         tr_bgc_N,       & ! if .true., algal nitrogen tracers  (n_algae)
         tr_bgc_DON,     & ! if .true., DON pools are tracers  (n_don)
         tr_bgc_C,       & ! if .true., algal carbon tracers + DOC and DIC 
         tr_bgc_chl,     & ! if .true., algal chlorophyll tracers 
         tr_bgc_Am,      & ! if .true., ammonia/um as nutrient tracer 
         tr_bgc_Sil,     & ! if .true., silicon as nutrient tracer 
         tr_bgc_DMS,     & ! if .true., DMS as  tracer 
         tr_bgc_Fe,      & ! if .true., Fe as  tracer 
         tr_bgc_PON,     & ! if .true., PON as tracer 
         tr_bgc_hum,     & ! if .true., humic material as tracer
         solve_zsal,     & ! if true, update salinity profile from solve_S_dt
         z_tracers,      & ! if .true., bgc or aerosol tracers are vertically resolved
         solve_zbgc,     & ! if .true., solve vertical biochemistry portion of code
         dEdd_algae,     & ! if .true., algal absorption of Shortwave is computed in the 
         skl_bgc           ! if true, solve skeletal biochemistry

       integer (kind=int_kind), intent(out) :: &
         nt_fbri,      & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
         nt_bgc_Nit,   & ! nutrients  
         nt_bgc_Am,    & ! 
         nt_bgc_Sil,   & !
         nt_bgc_DMSPp, & ! trace gases (skeletal layer)
         nt_bgc_DMSPd, & ! 
         nt_bgc_DMS,   & ! 
         nt_bgc_PON,   & ! zooplankton and detritus 
         nt_bgc_hum,   & ! humic material 
                         ! bio layer indicess
         nlt_bgc_Nit,  & ! nutrients  
         nlt_bgc_Am,   & ! 
         nlt_bgc_Sil,  & !
         nlt_bgc_DMSPp,& ! trace gases (skeletal layer)
         nlt_bgc_DMSPd,& ! 
         nlt_bgc_DMS,  & ! 
         nlt_bgc_PON,  & ! zooplankton and detritus 
         nlt_bgc_hum,  & ! humic material 
         nlt_chl_sw,   & ! points to total chla in trcrn_sw
         nt_zbgc_frac, & ! fraction of tracer in the mobile phase
         nt_bgc_S        ! Bulk salinity in fraction ice with dynamic salinity (Bio grid)

      integer (kind=int_kind), dimension(:), intent(out) :: &  
         nt_bgc_N , & ! diatoms, phaeocystis, pico/small   
         nt_bgc_C , & ! diatoms, phaeocystis, pico/small   
         nt_bgc_chl,& ! diatoms, phaeocystis, pico/small  
         nlt_bgc_N ,& ! diatoms, phaeocystis, pico/small   
         nlt_bgc_C ,& ! diatoms, phaeocystis, pico/small   
         nlt_bgc_chl   ! diatoms, phaeocystis, pico/small 

      integer (kind=int_kind), dimension(:), intent(out) :: &  
         nt_bgc_DOC,   & !  dissolved organic carbon  
         nlt_bgc_DOC     !  dissolved organic carbon

      integer (kind=int_kind), dimension(:), intent(out) :: & 
         nt_bgc_DON,   & !  dissolved organic nitrogen
         nlt_bgc_DON     !  dissolved organic nitrogen

      integer (kind=int_kind), dimension(:), intent(out) :: &  
         nt_bgc_DIC,    & !  dissolved inorganic carbon
         nlt_bgc_DIC      !  dissolved inorganic carbon

      integer (kind=int_kind), dimension(:), intent(out) :: & 
         nt_bgc_Fed,     & !  dissolved iron
         nt_bgc_Fep,     & !  particulate iron
         nlt_bgc_Fed,    & !  dissolved iron
         nlt_bgc_Fep       !  particulate iron

      integer (kind=int_kind), dimension(:), intent(out) :: &  
         nt_zaero,    & !  black carbon and other aerosols 
         nlt_zaero,   & !  black carbon and other aerosols
         nlt_zaero_sw   
    
      integer (kind=int_kind), dimension(:), intent(out) :: &   
         bio_index_o , & ! nlt  to appropriate value in ocean data array
         bio_index       ! nlt to nt

      real (kind=real_kind), intent(in) :: &
         initbio_frac, & ! fraction of ocean tracer concentration used to initialize tracer 
         frazil_scav     ! multiple of ocean tracer concentration due to frazil scavenging

      real (kind=real_kind), intent(in) :: &
        ratio_Si2N_diatoms, &   ! algal Si to N (mol/mol)
        ratio_Si2N_sp     , &
        ratio_Si2N_phaeo  , &
        ratio_S2N_diatoms , &   ! algal S  to N (mol/mol)
        ratio_S2N_sp      , &
        ratio_S2N_phaeo   , &
        ratio_Fe2C_diatoms, &   ! algal Fe to C  (umol/mol)
        ratio_Fe2C_sp     , &
        ratio_Fe2C_phaeo  , &
        ratio_Fe2N_diatoms, &   ! algal Fe to N  (umol/mol)
        ratio_Fe2N_sp     , &
        ratio_Fe2N_phaeo  , &
        ratio_Fe2DON      , &   ! Fe to N of DON (nmol/umol)
        ratio_Fe2DOC_s    , &   ! Fe to C of DOC (nmol/umol) saccharids
        ratio_Fe2DOC_l    , &   ! Fe to C of DOC (nmol/umol) lipids 
        tau_min           , &   ! rapid mobile to stationary exchanges (s) = 1.5 hours
        tau_max           , &   ! long time mobile to stationary exchanges (s) = 2 days
        chlabs_diatoms   , & ! chl absorption (1/m/(mg/m^3))
        chlabs_sp        , & !
        chlabs_phaeo     , & !
        alpha2max_low_diatoms , & ! light limitation (1/(W/m^2))  
        alpha2max_low_sp      , & 
        alpha2max_low_phaeo   , & 
        beta2max_diatoms , & ! light inhibition (1/(W/m^2))  
        beta2max_sp      , & 
        beta2max_phaeo   , & 
        mu_max_diatoms   , & ! maximum growth rate (1/day)       
        mu_max_sp        , & 
        mu_max_phaeo     , & 
        grow_Tdep_diatoms, & ! Temperature dependence of growth (1/C)
        grow_Tdep_sp     , & 
        grow_Tdep_phaeo  , & 
        fr_graze_diatoms , & ! Fraction grazed
        fr_graze_sp      , & 
        fr_graze_phaeo   , & 
        mort_pre_diatoms , & ! Mortality (1/day)
        mort_pre_sp      , & 
        mort_pre_phaeo   , & 
        mort_Tdep_diatoms, & ! T dependence of mortality (1/C)
        mort_Tdep_sp     , &  
        mort_Tdep_phaeo  , &  
        k_exude_diatoms  , & ! algal exudation (1/d)
        k_exude_sp       , &  
        k_exude_phaeo    , &  
        K_Nit_diatoms    , & ! nitrate half saturation (mmol/m^3)
        K_Nit_sp        , &  
        K_Nit_phaeo      , &  
        K_Am_diatoms     , & ! ammonium half saturation (mmol/m^3)
        K_Am_sp         , &   
        K_Am_phaeo       , &   
        K_Sil_diatoms    , & ! silicate half saturation (mmol/m^3)
        K_Sil_sp        , &   
        K_Sil_phaeo      , &   
        K_Fe_diatoms     , & ! iron half saturation (nM)
        K_Fe_sp         , &   
        K_Fe_phaeo       , &    
        f_don_protein    , & ! fraction of spilled grazing to proteins          
        kn_bac_protein   , & ! Bacterial degredation of DON (1/d)               
        f_don_Am_protein , & ! fraction of remineralized DON to ammonium        
        f_doc_s         , & ! fraction of mortality to DOC 
        f_doc_l         , &   
        f_exude_s        , & ! fraction of exudation to DOC
        f_exude_l        , & 
        k_bac_s         , & ! Bacterial degredation of DOC (1/d)
        k_bac_l         , & 
        algaltype_diatoms  , & ! mobility type
        algaltype_sp       , & !
        algaltype_phaeo    , & !
        nitratetype        , & !
        ammoniumtype       , & !
        silicatetype       , & !
        dmspptype         , & !
        dmspdtype         , & !
        humtype           , & !
        doctype_s         , & !
        doctype_l         , & !
        dictype_1         , & !
        dontype_protein    , & !
        fedtype_1         , & !
        feptype_1         , & !
        zaerotype_bc1      , & !
        zaerotype_bc2      , & !
        zaerotype_dust1    , & !
        zaerotype_dust2    , & !
        zaerotype_dust3    , & !
        zaerotype_dust4    , & !
        ratio_C2N_diatoms  , & ! algal C to N ratio (mol/mol)
        ratio_C2N_sp       , & !
        ratio_C2N_phaeo    , & !
        ratio_chl2N_diatoms, & ! algal chlorophyll to N ratio (mg/mmol)
        ratio_chl2N_sp     , & !
        ratio_chl2N_phaeo  , & !
        F_abs_chl_diatoms  , & ! scales absorbed radiation for dEdd
        F_abs_chl_sp       , & !
        F_abs_chl_phaeo    , & !
        ratio_C2N_proteins     ! ratio of C to N in proteins (mol/mol)   

      ! local variables
      real (kind=dbl_kind) :: &
         nitratetype_R8, &
         ammoniumtype_R8, &
         silicatetype_R8, &
         dmspptype_R8, &
         dmspdtype_R8, &
         humtype_R8
      real (kind=dbl_kind), dimension (:,:), allocatable :: &
         trcr_base     ! = 0 or 1 depending on tracer dependency
                       ! argument 2:  (1) aice, (2) vice, (3) vsno
      integer (kind=int_kind) :: &
        k, mm    , & ! loop index  
        ntd      , & ! for tracer dependency calculation
        nk       , & !
        nt_depend

     trcr_base =real(trcr_base_R4, kind=dbl_kind)
      nitratetype_R8 = real(nitratetype, kind = dbl_kind)
      ammoniumtype_R8 = real(ammoniumtype, kind = dbl_kind)
      silicatetype_R8 = real(silicatetype, kind = dbl_kind)
      dmspptype_R8 = real(dmspptype, kind = dbl_kind)
      dmspdtype_R8 = real(dmspdtype, kind = dbl_kind)
      humtype_R8 = real(humtype, kind = dbl_kind)
      ntrcr_o = ntrcr
      nt_fbri = 0
      if (tr_brine) then
          nt_fbri = ntrcr + 1   ! ice volume fraction with salt
          ntrcr = ntrcr + 1
          trcr_depend(nt_fbri)   = 1   ! volume-weighted
          trcr_base  (nt_fbri,1) = c0  ! volume-weighted
          trcr_base  (nt_fbri,2) = c1  ! volume-weighted
          trcr_base  (nt_fbri,3) = c0  ! volume-weighted
          n_trcr_strata(nt_fbri) = 0
          nt_strata  (nt_fbri,1) = 0
          nt_strata  (nt_fbri,2) = 0
      endif

      ntd = 0                    ! if nt_fbri /= 0 then use fbri dependency
      if (nt_fbri == 0) ntd = -1 ! otherwise make tracers depend on ice volume

      if (solve_zsal) then       ! .true. only if tr_brine = .true.
          nt_bgc_S = ntrcr + 1
          ntrcr = ntrcr + nblyr
          do k = 1,nblyr
             trcr_depend(nt_bgc_S + k - 1) = 2 + nt_fbri + ntd
             trcr_base  (nt_bgc_S,1) = c0  ! default: ice area
             trcr_base  (nt_bgc_S,2) = c1 
             trcr_base  (nt_bgc_S,3) = c0  
             n_trcr_strata(nt_bgc_S) = 1
             nt_strata(nt_bgc_S,1) = nt_fbri
             nt_strata(nt_bgc_S,2) = 0
          enddo
      endif 

      !-----------------------------------------------------------------
      ! biogeochemistry
      !-----------------------------------------------------------------

      nbtrcr = 0
      nbtrcr_sw = 0

      ! vectors of size max_algae
      nlt_bgc_N(:) = 0
      nlt_bgc_C(:) = 0
      nlt_bgc_chl(:) = 0
      nt_bgc_N(:) = 0
      nt_bgc_C(:) = 0
      nt_bgc_chl(:) = 0

      ! vectors of size max_dic
      nlt_bgc_DIC(:) = 0
      nt_bgc_DIC(:) = 0

      ! vectors of size max_doc
      nlt_bgc_DOC(:) = 0
      nt_bgc_DOC(:) = 0

      ! vectors of size max_don
      nlt_bgc_DON(:) = 0
      nt_bgc_DON(:) = 0

      ! vectors of size max_fe 
      nlt_bgc_Fed(:) = 0
      nlt_bgc_Fep(:) = 0
      nt_bgc_Fed(:) = 0
      nt_bgc_Fep(:) = 0

      ! vectors of size max_aero
      nlt_zaero(:) = 0
      nlt_zaero_sw(:) = 0
      nt_zaero(:) = 0

      nlt_bgc_Nit    = 0
      nlt_bgc_Am     = 0
      nlt_bgc_Sil    = 0
      nlt_bgc_DMSPp  = 0
      nlt_bgc_DMSPd  = 0
      nlt_bgc_DMS    = 0
      nlt_bgc_PON    = 0
      nlt_bgc_hum    = 0
      nlt_chl_sw     = 0
      bio_index(:)   = 0
      bio_index_o(:) = 0

      nt_bgc_Nit    = 0
      nt_bgc_Am     = 0
      nt_bgc_Sil    = 0
      nt_bgc_DMSPp  = 0
      nt_bgc_DMSPd  = 0
      nt_bgc_DMS    = 0
      nt_bgc_PON    = 0
      nt_bgc_hum    = 0

      !-----------------------------------------------------------------
      ! Define array parameters
      !-----------------------------------------------------------------
      R_Si2N(1) = real(ratio_Si2N_diatoms, kind = dbl_kind)
      R_Si2N(2) = real(ratio_Si2N_sp, kind = dbl_kind)
      R_Si2N(3) = real(ratio_Si2N_phaeo, kind = dbl_kind)

      R_S2N(1) = real(ratio_S2N_diatoms, kind = dbl_kind)
      R_S2N(2) = real(ratio_S2N_sp, kind = dbl_kind)
      R_S2N(3) = real(ratio_S2N_phaeo, kind = dbl_kind)

      R_Fe2C(1) = real(ratio_Fe2C_diatoms, kind = dbl_kind)
      R_Fe2C(2) = real(ratio_Fe2C_sp, kind = dbl_kind)
      R_Fe2C(3) = real(ratio_Fe2C_phaeo, kind = dbl_kind)

      R_Fe2N(1) = real(ratio_Fe2N_diatoms, kind = dbl_kind)
      R_Fe2N(2) = real(ratio_Fe2N_sp, kind = dbl_kind)
      R_Fe2N(3) = real(ratio_Fe2N_phaeo, kind = dbl_kind)

      R_C2N(1) = real(ratio_C2N_diatoms, kind = dbl_kind)
      R_C2N(2) = real(ratio_C2N_sp, kind = dbl_kind)
      R_C2N(3) = real(ratio_C2N_phaeo, kind = dbl_kind)

      R_chl2N(1) = real(ratio_chl2N_diatoms, kind = dbl_kind)
      R_chl2N(2) = real(ratio_chl2N_sp, kind = dbl_kind)
      R_chl2N(3) = real(ratio_chl2N_phaeo, kind = dbl_kind)

      F_abs_chl(1) = real(F_abs_chl_diatoms, kind = dbl_kind)
      F_abs_chl(2) = real(F_abs_chl_sp, kind = dbl_kind)
      F_abs_chl(3) = real(F_abs_chl_phaeo, kind = dbl_kind)

      R_Fe2DON(1) = real(ratio_Fe2DON, kind = dbl_kind)
      R_C2N_DON(1) = real(ratio_C2N_proteins, kind = dbl_kind)
     
      R_Fe2DOC(1) = real(ratio_Fe2DOC_s, kind = dbl_kind)
      R_Fe2DOC(2) = real(ratio_Fe2DOC_l, kind = dbl_kind)
      R_Fe2DOC(3) = real(c0, kind = dbl_kind)

      chlabs(1) = real(chlabs_diatoms, kind = dbl_kind)
      chlabs(2) = real(chlabs_sp, kind = dbl_kind)
      chlabs(3) = real(chlabs_phaeo, kind = dbl_kind)

      alpha2max_low(1) = real(alpha2max_low_diatoms, kind = dbl_kind)
      alpha2max_low(2) = real(alpha2max_low_sp, kind = dbl_kind)
      alpha2max_low(3) = real(alpha2max_low_phaeo, kind = dbl_kind)

      beta2max(1) = real(beta2max_diatoms, kind = dbl_kind)
      beta2max(2) = real(beta2max_sp, kind = dbl_kind)
      beta2max(3) = real(beta2max_phaeo, kind = dbl_kind)

      mu_max(1) = real(mu_max_diatoms, kind = dbl_kind)
      mu_max(2) = real(mu_max_sp, kind = dbl_kind)
      mu_max(3) = real(mu_max_phaeo, kind = dbl_kind)

      grow_Tdep(1) = real(grow_Tdep_diatoms, kind = dbl_kind)
      grow_Tdep(2) = real(grow_Tdep_sp, kind = dbl_kind)
      grow_Tdep(3) = real(grow_Tdep_phaeo, kind = dbl_kind)

      fr_graze(1) = real(fr_graze_diatoms, kind = dbl_kind)
      fr_graze(2) = real(fr_graze_sp, kind = dbl_kind)
      fr_graze(3) = real(fr_graze_phaeo, kind = dbl_kind)

      mort_pre(1) = real(mort_pre_diatoms, kind = dbl_kind)
      mort_pre(2) = real(mort_pre_sp, kind = dbl_kind)
      mort_pre(3) = real(mort_pre_phaeo, kind = dbl_kind)

      mort_Tdep(1) = real(mort_Tdep_diatoms, kind = dbl_kind)
      mort_Tdep(2) = real(mort_Tdep_sp, kind = dbl_kind)
      mort_Tdep(3) = real(mort_Tdep_phaeo, kind = dbl_kind)

      k_exude(1) = real(k_exude_diatoms, kind = dbl_kind)
      k_exude(2) = real(k_exude_sp, kind = dbl_kind)
      k_exude(3) = real(k_exude_phaeo, kind = dbl_kind)

      K_Nit(1) = real(K_Nit_diatoms, kind = dbl_kind)
      K_Nit(2) = real(K_Nit_sp, kind = dbl_kind)
      K_Nit(3) = real(K_Nit_phaeo, kind = dbl_kind)

      K_Am(1) = real(K_Am_diatoms, kind = dbl_kind)
      K_Am(2) = real(K_Am_sp, kind = dbl_kind)
      K_Am(3) = real(K_Am_phaeo, kind = dbl_kind)

      K_Sil(1) = real(K_Sil_diatoms, kind = dbl_kind)
      K_Sil(2) = real(K_Sil_sp, kind = dbl_kind)
      K_Sil(3) = real(K_Sil_phaeo, kind = dbl_kind)

      K_Fe(1) = real(K_Fe_diatoms, kind = dbl_kind)
      K_Fe(2) = real(K_Fe_sp, kind = dbl_kind)
      K_Fe(3) = real(K_Fe_phaeo, kind = dbl_kind)

      f_doc(1) = real(f_doc_s, kind = dbl_kind)
      f_doc(2) = real(f_doc_l, kind = dbl_kind)

      f_don(1) = real(f_don_protein, kind = dbl_kind)
      kn_bac(1) = real(kn_bac_protein, kind = dbl_kind)
      f_don_Am(1) = real(f_don_Am_protein, kind = dbl_kind)

      f_exude(1) = real(f_exude_s, kind = dbl_kind)
      f_exude(2) = real(f_exude_l, kind = dbl_kind)
      k_bac(1) = real(k_bac_s, kind = dbl_kind)
      k_bac(2) = real(k_bac_l, kind = dbl_kind)

      algaltype(1) = real(algaltype_diatoms, kind = dbl_kind)
      algaltype(2) = real(algaltype_sp, kind = dbl_kind)
      algaltype(3) = real(algaltype_phaeo, kind = dbl_kind)

      doctype(1) = real(doctype_s, kind = dbl_kind)
      doctype(2) = real(doctype_l, kind = dbl_kind)
 
      dictype(1) = real(dictype_1, kind = dbl_kind)

      dontype(1) = real(dontype_protein, kind = dbl_kind)

      fedtype(1) = real(fedtype_1, kind = dbl_kind)
      feptype(1) = real(feptype_1, kind = dbl_kind)

      zaerotype(1) = real(zaerotype_bc1, kind = dbl_kind)
      zaerotype(2) = real(zaerotype_bc2, kind = dbl_kind)
      zaerotype(3) = real(zaerotype_dust1, kind = dbl_kind)
      zaerotype(4) = real(zaerotype_dust2, kind = dbl_kind)
      zaerotype(5) = real(zaerotype_dust3, kind = dbl_kind)
      zaerotype(6) = real(zaerotype_dust4, kind = dbl_kind)     

      if (skl_bgc) then

         nk = 1
         nt_depend = 0

         if (dEdd_algae) then 
           nlt_chl_sw = 1
           nbtrcr_sw = nilyr+nslyr+2  ! only the bottom layer 
                                                 ! will be nonzero    
         endif  
         
      elseif (z_tracers) then ! defined on nblyr+1 in ice 
                              ! and 2 snow layers (snow surface + interior)

         nk = nblyr + 1
         nt_depend = 2 + nt_fbri + ntd 

         if (tr_bgc_N) then
            if (dEdd_algae) then
               nlt_chl_sw = 1
               nbtrcr_sw =  nilyr+nslyr+2
            endif
         endif ! tr_bgc_N

      endif ! skl_bgc or z_tracers

      if (skl_bgc .or. z_tracers) then

      !-----------------------------------------------------------------
      ! assign tracer indices and dependencies
      ! bgc_tracer_type: < 0  purely mobile , >= 0 stationary 
      !------------------------------------------------------------------

      if (tr_bgc_N) then
         do mm = 1, n_algae      
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_N(mm),    nlt_bgc_N(mm), &
                                      algaltype(mm),   nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_N(mm)) = mm
         enddo   ! mm
      endif ! tr_bgc_N

      if (tr_bgc_Nit) then
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_Nit,      nlt_bgc_Nit,   &
                                      nitratetype_R8,     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Nit) = max_algae + 1
      endif ! tr_bgc_Nit
         
      if (tr_bgc_C) then
       !
       ! Algal C is not yet distinct from algal N
       ! * Reqires exudation and/or changing C:N ratios
       ! for implementation
       !
       !  do mm = 1,n_algae      
       !     call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
       !                               nt_bgc_C(mm),    nlt_bgc_C(mm), &
       !                               algaltype(mm),   nt_depend,     &
       !                               ntrcr,           nbtrcr,        &
       !                               bgc_tracer_type, trcr_depend,   &
       !                               trcr_base,       n_trcr_strata, &
       !                               nt_strata,       bio_index)
       !     bio_index_o(nlt_bgc_C(mm)) = max_algae + 1 + mm
       !  enddo   ! mm

         do mm = 1, n_doc
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_DOC(mm),  nlt_bgc_DOC(mm), &
                                      doctype(mm),     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DOC(mm)) = max_algae + 1 + mm
         enddo   ! mm
         do mm = 1, n_dic
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_DIC(mm),  nlt_bgc_DIC(mm), &
                                      dictype(mm),     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DIC(mm)) = max_algae + max_doc + 1 + mm
         enddo   ! mm
      endif      ! tr_bgc_C

      if (tr_bgc_chl) then
         do mm = 1, n_algae
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_chl(mm),  nlt_bgc_chl(mm), &
                                      algaltype(mm),   nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_chl(mm)) = max_algae + 1 + max_doc + max_dic + mm
         enddo   ! mm
      endif      ! tr_bgc_chl

      if (tr_bgc_Am) then
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_Am,       nlt_bgc_Am,    &
                                      ammoniumtype_R8,    nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Am) = 2*max_algae + max_doc + max_dic + 2
      endif    
      if (tr_bgc_Sil) then
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_Sil,      nlt_bgc_Sil,   &
                                      silicatetype_R8,    nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Sil) = 2*max_algae + max_doc + max_dic + 3
      endif    
      if (tr_bgc_DMS) then   ! all together
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_DMSPp,    nlt_bgc_DMSPp, &
                                      dmspptype_R8,       nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DMSPp) = 2*max_algae + max_doc + max_dic + 4

            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_DMSPd,    nlt_bgc_DMSPd, &
                                      dmspdtype_R8,       nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DMSPd) = 2*max_algae + max_doc + max_dic + 5

            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_DMS,      nlt_bgc_DMS,   &
                                      dmspdtype_R8,       nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DMS) = 2*max_algae + max_doc + max_dic + 6
      endif    
      if (tr_bgc_PON) then
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_PON,      nlt_bgc_PON, &
                                      nitratetype_R8,     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_PON) =  2*max_algae + max_doc + max_dic + 7
      endif
      if (tr_bgc_DON) then
         do mm = 1, n_don
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_DON(mm),  nlt_bgc_DON(mm), &
                                      dontype(mm),     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DON(mm)) = 2*max_algae + max_doc + max_dic + 7 + mm
         enddo   ! mm
      endif      ! tr_bgc_DON
      if (tr_bgc_Fe) then
         do mm = 1, n_fed
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_Fed(mm),  nlt_bgc_Fed(mm), &
                                      fedtype(mm),     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Fed(mm)) = 2*max_algae + max_doc + max_dic &
                                         + max_don + 7 + mm
         enddo   ! mm
         do mm = 1, n_fep
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_Fep(mm),  nlt_bgc_Fep(mm), &
                                      feptype(mm),     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Fep(mm)) = 2*max_algae + max_doc + max_dic &
                                         + max_don + max_fe + 7 + mm
         enddo   ! mm
      endif      ! tr_bgc_Fe 
  
      if (tr_bgc_hum) then
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_hum,      nlt_bgc_hum,   &
                                      humtype_R8,         nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_hum) =   2*max_algae + max_doc + 8 + max_dic &
                                         + max_don + 2*max_fe + max_aero 
      endif
      endif  ! skl_bgc or z_tracers

      if (z_tracers) then ! defined on nblyr+1 in ice 
                              ! and 2 snow layers (snow surface + interior)

         nk = nblyr + 1
         nt_depend = 2 + nt_fbri + ntd 

         ! z layer aerosols
         if (tr_zaero) then
            do mm = 1, n_zaero
               if (dEdd_algae) then
                  nlt_zaero_sw(mm) = nbtrcr_sw + 1
                  nbtrcr_sw = nbtrcr_sw + nilyr + nslyr+2
               endif
               call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                         nt_zaero(mm),    nlt_zaero(mm), &
                                         zaerotype(mm),   nt_depend,     &
                                         ntrcr,           nbtrcr,        &
                                         bgc_tracer_type, trcr_depend,   &
                                         trcr_base,       n_trcr_strata, &
                                         nt_strata,       bio_index)
               bio_index_o(nlt_zaero(mm)) = 2*max_algae + max_doc + max_dic &
                                          + max_don + 2*max_fe + 7 + mm
            enddo   ! mm
         endif      ! tr_zaero

         nt_zbgc_frac = 0
         if (nbtrcr > 0) then
            nt_zbgc_frac = ntrcr + 1
            ntrcr = ntrcr + nbtrcr
            do k = 1,nbtrcr 
               zbgc_frac_init(k) = c1   
               trcr_depend(nt_zbgc_frac+k-1) =  2+nt_fbri 
               trcr_base(nt_zbgc_frac+ k - 1,1)  = c0
               trcr_base(nt_zbgc_frac+ k - 1,2)  = c1
               trcr_base(nt_zbgc_frac+ k - 1,3)  = c0
               n_trcr_strata(nt_zbgc_frac+ k - 1)= 1  
               nt_strata(nt_zbgc_frac+ k - 1,1)  = nt_fbri
               nt_strata(nt_zbgc_frac+ k - 1,2)  = 0 
               tau_ret(k) = c1
               tau_rel(k) = c1
               if (bgc_tracer_type(k) >=  c0 .and. bgc_tracer_type(k) < p5) then
                  tau_ret(k) = tau_min
                  tau_rel(k) = tau_max
                  zbgc_frac_init(k) = c1
               elseif (bgc_tracer_type(k) >= p5 .and. bgc_tracer_type(k) < c1) then
                  tau_ret(k) = tau_min
                  tau_rel(k) = tau_min
                  zbgc_frac_init(k) = c1
               elseif (bgc_tracer_type(k) >= c1 .and. bgc_tracer_type(k) < c2) then
                  tau_ret(k) = tau_max
                  tau_rel(k) = tau_min
                  zbgc_frac_init(k) = c1
               elseif (bgc_tracer_type(k) >= c2 ) then
                  tau_ret(k) = tau_max
                  tau_rel(k) = tau_max
                  zbgc_frac_init(k) = c1
               endif
            enddo
         endif

      endif ! z_tracers

      do k = 1, nbtrcr
         zbgc_init_frac(k) = real(frazil_scav, kind = dbl_kind)
         if (bgc_tracer_type(k) < c0)  zbgc_init_frac(k) = real(initbio_frac, kind=dbl_kind)
      enddo  
      trcr_base_R4 = real(trcr_base, kind=real_kind)
      deallocate(trcr_base)
      end subroutine colpkg_init_zbgc

!=======================================================================

      subroutine colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc,          nlt_bgc,       &
                                      bgctype,         nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)

      use ice_constants_colpkg, only: c0, c1

      integer (kind=int_kind), intent(in) :: &
         nk           , & ! counter
         nt_depend    , & ! tracer dependency index
         nt_fbri

      integer (kind=int_kind), intent(inout) :: &
         ntrcr        , & ! number of tracers
         nbtrcr       , & ! number of bio tracers
         nt_bgc       , & ! tracer index
         nlt_bgc          ! bio tracer index

      integer (kind=int_kind), dimension(:), intent(inout) :: &
         trcr_depend  , & ! tracer dependencies
         n_trcr_strata, & ! number of underlying tracer layers
         bio_index        !

      integer (kind=int_kind), dimension(:,:), intent(inout) :: &
         nt_strata        ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         trcr_base        ! = 0 or 1 depending on tracer dependency
                          ! argument 2:  (1) aice, (2) vice, (3) vsno

      real (kind=dbl_kind) :: &
         bgctype          ! bio tracer transport type (mobile vs stationary)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         bgc_tracer_type  ! bio tracer transport type array

      ! local variables

      integer (kind=int_kind) :: &
         k         , & ! loop index
         n_strata  , & ! temporary values
         nt_strata1, & ! 
         nt_strata2

      real (kind=dbl_kind) :: &
         trcr_base1, & ! temporary values
         trcr_base2, &
         trcr_base3
         nt_bgc = ntrcr + 1 
         nbtrcr = nbtrcr + 1
         nlt_bgc = nbtrcr
         bgc_tracer_type(nbtrcr) = bgctype
   
         if (nk > 1) then 
            ! include vertical bgc in snow
            do k = nk, nk+1
               ntrcr = ntrcr + 1
               trcr_depend  (nt_bgc + k  ) = 2 ! snow volume
               trcr_base    (nt_bgc + k,1) = c0
               trcr_base    (nt_bgc + k,2) = c0
               trcr_base    (nt_bgc + k,3) = c1
               n_trcr_strata(nt_bgc + k  ) = 0
               nt_strata    (nt_bgc + k,1) = 0
               nt_strata    (nt_bgc + k,2) = 0
            enddo

            trcr_base1 = c0      
            trcr_base2 = c1     
            trcr_base3 = c0
            n_strata = 1    
            nt_strata1 = nt_fbri
            nt_strata2 = 0
         else  ! nk = 1
            trcr_base1 = c1
            trcr_base2 = c0
            trcr_base3 = c0
            n_strata = 0
            nt_strata1 = 0
            nt_strata2 = 0
         endif ! nk

         do k = 1, nk     !in ice
            ntrcr = ntrcr + 1
            trcr_depend  (nt_bgc + k - 1  ) = nt_depend
            trcr_base    (nt_bgc + k - 1,1) = trcr_base1
            trcr_base    (nt_bgc + k - 1,2) = trcr_base2
            trcr_base    (nt_bgc + k - 1,3) = trcr_base3
            n_trcr_strata(nt_bgc + k - 1  ) = n_strata
            nt_strata    (nt_bgc + k - 1,1) = nt_strata1
            nt_strata    (nt_bgc + k - 1,2) = nt_strata2
         enddo

         bio_index (nlt_bgc) = nt_bgc

      end subroutine colpkg_init_bgc_trcr

!=======================================================================
!     Temperature functions
!=======================================================================

      function colpkg_liquidus_temperature(Sin) result(Tmlt)

        use ice_colpkg_shared, only: ktherm
        use ice_constants_colpkg, only: depressT
        use ice_mushy_physics, only: liquidus_temperature_mush

        real(dbl_kind), intent(in) :: Sin
        real(dbl_kind) :: Tmlt

        if (ktherm == 2) then

           Tmlt = liquidus_temperature_mush(Sin)

        else

           Tmlt = -depressT * Sin

        endif

      end function colpkg_liquidus_temperature

!=======================================================================

      function colpkg_sea_freezing_temperature(sss) result(Tf)

        use ice_colpkg_shared, only: tfrz_option
        use ice_constants_colpkg, only: depressT, Tocnfrz

        real(dbl_kind), intent(in) :: sss
        real(dbl_kind) :: Tf

        if (trim(tfrz_option) == 'mushy') then

           Tf = colpkg_liquidus_temperature(sss) ! deg C
           
        elseif (trim(tfrz_option) == 'linear_salt') then

           Tf = -depressT * sss ! deg C

        else

           Tf = Tocnfrz 

        endif

      end function colpkg_sea_freezing_temperature

!=======================================================================

      function colpkg_ice_temperature(qin, Sin) result(Tin)

        use ice_colpkg_shared, only: ktherm
        use ice_constants_colpkg, only: depressT
        use ice_mushy_physics, only: temperature_mush
        use ice_therm_shared, only: calculate_Tin_from_qin

        real(kind=dbl_kind), intent(in) :: qin, Sin
        real(kind=dbl_kind) :: Tin

        real(kind=dbl_kind) :: Tmlts

        if (ktherm == 2) then

           Tin = temperature_mush(qin, Sin)

        else

           Tmlts = -depressT * Sin
           Tin = calculate_Tin_from_qin(qin,Tmlts)

        endif

      end function colpkg_ice_temperature

   

!=======================================================================

      function colpkg_snow_temperature(qin) result(Tsn)

        use ice_colpkg_shared, only: ktherm
        use ice_mushy_physics, only: temperature_snow
        use ice_constants_colpkg, only: Lfresh, rhos, cp_ice

        real(kind=dbl_kind), intent(in) :: qin
        real(kind=dbl_kind) :: Tsn

        if (ktherm == 2) then

           Tsn = temperature_snow(qin)

        else

           Tsn = (Lfresh + qin/rhos)/cp_ice

        endif

      end function colpkg_snow_temperature

!=======================================================================

      function colpkg_enthalpy_ice(zTin, zSin) result(qin)

        use ice_colpkg_shared, only: ktherm
        use ice_mushy_physics, only: enthalpy_mush
        use ice_constants_colpkg, only: depressT, rhoi, cp_ice, Lfresh, cp_ocn, c1

        real(kind=dbl_kind), intent(in) :: zTin
        real(kind=dbl_kind), intent(in) :: zSin
        real(kind=dbl_kind) :: qin

        real(kind=dbl_kind) :: Tmlt

        if (ktherm == 2) then

           qin = enthalpy_mush(zTin, zSin)

        else

           Tmlt = -zSin*depressT
           qin = -(rhoi * (cp_ice*(Tmlt-zTin) + Lfresh*(c1-Tmlt/zTin) - cp_ocn*Tmlt))

        endif

      end function colpkg_enthalpy_ice

!=======================================================================

      function colpkg_enthalpy_snow(zTsn) result(qsn)

        use ice_mushy_physics, only: enthalpy_snow

        real(kind=dbl_kind), intent(in) :: zTsn
        real(kind=dbl_kind) :: qsn

        qsn = enthalpy_snow(zTsn)

      end function colpkg_enthalpy_snow

!=======================================================================
!     Time-stepping routines
!=======================================================================

! Driver for thermodynamic changes not needed for coupling:
! transport in thickness space, lateral growth and melting.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine colpkg_step_therm1(dt_R4, ncat, nilyr, nslyr, n_aero, &
                                    aice0_R4       ,               &
                                    aicen_init_R4  ,               &
                                    vicen_init_R4  , vsnon_init_R4  , &
                                    aice_R4        , aicen_R4       , &
                                    vice_R4        , vicen_R4       , &
                                    vsno_R4        , vsnon_R4       , &
                                    uvel_R4        , vvel_R4        , &
                                    Tsfc_R4        , zqsn_R4        , &
                                    zqin_R4        , zSin_R4        , &
                                    smice_R4       , smliq_R4       , &
                                    alvl_R4        , vlvl_R4        , &
                                    apnd_R4        , hpnd_R4        , &
                                    ipnd_R4        ,                  &
                                    iage_R4        , FY_R4          , &
                                    rsnw_R4        , use_smliq_pnd  , &
                                    aerosno_R4     , aeroice_R4     , &
                                    uatm_R4        , vatm_R4        , &
                                    wind_R4        , zlvl_R4        , &
                                    Qa_R4          , rhoa_R4        , &
                                    Tair_R4        , Tref_R4        , &
                                    Qref_R4        , Uref_R4        , &
                                    Cdn_atm_ratio_R4,                 &
                                    Cdn_ocn_R4     , Cdn_ocn_skin_R4, &
                                    Cdn_ocn_floe_R4, Cdn_ocn_keel_R4, &
                                    Cdn_atm_R4     , Cdn_atm_skin_R4, &
                                    Cdn_atm_floe_R4, Cdn_atm_pond_R4, &
                                    Cdn_atm_rdg_R4 , hfreebd_R4     , &
                                    hdraft_R4      , hridge_R4      , &
                                    distrdg_R4     , hkeel_R4       , &
                                    dkeel_R4       , lfloe_R4       , &
                                    dfloe_R4       ,                  &
                                    strax_R4       , stray_R4       , &
                                    strairxT_R4    , strairyT_R4    , &
                                    potT_R4        , sst_R4         , &
                                    sss_R4         , Tf_R4          , &
                                    strocnxT_R4    , strocnyT_R4    , &
                                    fbot_R4        ,               &
                                    frzmlt_R4      , rside_R4       , &
                                    fsnow_R4       , frain_R4       , &
                                    fpond_R4       , fsloss_R4      , &
                                    fsurf_R4       , fsurfn_R4      , &
                                    fcondtop_R4    , fcondtopn_R4   , &
                                    fswsfcn_R4     , fswintn_R4     , &
                                    fswthrun_R4    , fswabs_R4      , &
                                    flwout_R4      ,               &
                                    Sswabsn_R4     , Iswabsn_R4     , &
                                    flw_R4         , coszen_R4      , & 
                                    fsens_R4       , fsensn_R4      , &
                                    flat_R4        , flatn_R4       , &
                                    evap_R4        ,               &
                                    fresh_R4       , fsalt_R4       , &
                                    fhocn_R4       , fswthru_R4     , &
                                    flatn_f_R4     , fsensn_f_R4    , &
                                    fsurfn_f_R4    , fcondtopn_f_R4 , &
                                    faero_atm_R4   , faero_ocn_R4   , &
                                    dhsn_R4        , ffracn_R4      , &
                                    meltt_R4       , melttn_R4      , &
                                    meltb_R4       , meltbn_R4      , &
                                    meltl_R4       ,               &
                                    melts_R4       , meltsn_R4      , &
                                    meltsliq_R4   , meltsliqn_R4   , &
                                    congel_R4      , congeln_R4     , &
                                    snoice_R4      , snoicen_R4     , &
                                    dsnown_R4      , frazil_R4      , &
                                    lmask_n     , lmask_s     , &
                                    mlt_onset_R4   , frz_onset_R4   , &
                                    yday_R4        , l_stop      , &
                                    stop_label  , prescribed_ice)

      use ice_aerosol, only: update_aerosol
      use ice_atmo, only: neutral_drag_coeffs
      use ice_age, only: increment_age
      use ice_constants_colpkg, only: rhofresh, rhoi, rhos, c0, c1, puny, &
          snwlvlfac
      use ice_firstyear, only: update_FYarea
      use ice_flux_colpkg, only: set_sfcflux, merge_fluxes
      use ice_meltpond_cesm, only: compute_ponds_cesm
      use ice_meltpond_lvl, only: compute_ponds_lvl
      use ice_meltpond_topo, only: compute_ponds_topo
      use ice_snow, only: drain_snow
      use ice_therm_shared, only: hi_min
      use ice_therm_vertical, only: frzmlt_bottom_lateral, thermo_vertical
      use ice_colpkg_tracers, only: tr_iage, tr_FY, tr_aero, tr_pond, &
          tr_pond_cesm, tr_pond_lvl, tr_pond_topo, tr_snow, tr_rsnw

      integer (kind=int_kind), intent(in) :: &
         ncat    , & ! number of thickness categories
         nilyr   , & ! number of ice layers
         nslyr   , & ! number of snow layers
         n_aero      ! number of aerosol tracers in use

      real (kind=real_kind), intent(in) :: &
         dt_R4          , & ! time step
         uvel_R4        , & ! x-component of velocity (m/s)
         vvel_R4        , & ! y-component of velocity (m/s)
         strax_R4       , & ! wind stress components (N/m^2)
         stray_R4       , & ! 
         yday_R4            ! day of year
      
      logical (kind=log_kind), intent(in) :: &
         lmask_n     , & ! northern hemisphere mask
         lmask_s     , & ! southern hemisphere mask
         use_smliq_pnd   ! if true, use snow liquid tracer for ponds

      logical (kind=log_kind), intent(in), optional :: &
         prescribed_ice  ! if .true., use prescribed ice instead of computed
      
      real (kind=real_kind), intent(inout) :: &
         aice0_R4       , & ! open water fraction
         aice_R4       , & ! sea ice concentration
         vice_R4        , & ! volume per unit area of ice          (m)
         vsno_R4        , & ! volume per unit area of snow         (m)
         zlvl_R4        , & ! atm level height (m)
         uatm_R4        , & ! wind velocity components (m/s)
         vatm_R4        , &
         wind_R4        , & ! wind speed (m/s)
         potT_R4        , & ! air potential temperature  (K)
         Tair_R4        , & ! air temperature  (K)
         Qa_R4          , & ! specific humidity (kg/kg)
         rhoa_R4        , & ! air density (kg/m^3)
         frain_R4       , & ! rainfall rate (kg/m^2 s)
         fsnow_R4       , & ! snowfall rate (kg/m^2 s)
         fsloss_R4      , & ! blowing snow loss to leads (kg/m^2/s)
         fpond_R4       , & ! fresh water flux to ponds (kg/m^2/s)
         fresh_R4       , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt_R4       , & ! salt flux to ocean (kg/m^2/s)
         fhocn_R4      , & ! net heat flux to ocean (W/m^2)
         fswthru_R4     , & ! shortwave penetrating to ocean (W/m^2)
         fsurf_R4       , & ! net surface heat flux (excluding fcondtop)(W/m^2)
         fcondtop_R4    , & ! top surface conductive flux        (W/m^2)
         fsens_R4       , & ! sensible heat flux (W/m^2)
         flat_R4        , & ! latent heat flux   (W/m^2)
         fswabs_R4      , & ! shortwave flux absorbed in ice and ocean (W/m^2)
         coszen_R4      , & ! cosine solar zenith angle, < 0 for sun below horizon 
         flw_R4         , & ! incoming longwave radiation (W/m^2)
         flwout_R4      , & ! outgoing longwave radiation (W/m^2)
         evap_R4        , & ! evaporative water flux (kg/m^2/s)
         congel_R4      , & ! basal ice growth         (m/step-->cm/day)
         frazil_R4      , & ! frazil ice growth        (m/step-->cm/day)
         snoice_R4      , & ! snow-ice formation       (m/step-->cm/day)
         Tref_R4        , & ! 2m atm reference temperature (K)
         Qref_R4        , & ! 2m atm reference spec humidity (kg/kg)
         Uref_R4        , & ! 10m atm reference wind speed (m/s)
         Cdn_atm_R4     , & ! atm drag coefficient
         Cdn_ocn_R4     , & ! ocn drag coefficient
         hfreebd_R4     , & ! freeboard (m)
         hdraft_R4      , & ! draft of ice + snow column (Stoessel1993)
         hridge_R4      , & ! ridge height
         distrdg_R4     , & ! distance between ridges
         hkeel_R4       , & ! keel depth
         dkeel_R4       , & ! distance between keels
         lfloe_R4       , & ! floe length
         dfloe_R4       , & ! distance between floes
         Cdn_atm_skin_R4, & ! neutral skin drag coefficient
         Cdn_atm_floe_R4, & ! neutral floe edge drag coefficient
         Cdn_atm_pond_R4, & ! neutral pond edge drag coefficient
         Cdn_atm_rdg_R4 , & ! neutral ridge drag coefficient
         Cdn_ocn_skin_R4, & ! skin drag coefficient
         Cdn_ocn_floe_R4, & ! floe edge drag coefficient
         Cdn_ocn_keel_R4, & ! keel drag coefficient
         Cdn_atm_ratio_R4,& ! ratio drag atm / neutral drag atm
         strairxT_R4    , & ! stress on ice by air, x-direction
         strairyT_R4    , & ! stress on ice by air, y-direction
         strocnxT_R4    , & ! ice-ocean stress, x-direction
         strocnyT_R4    , & ! ice-ocean stress, y-direction
         fbot_R4        , & ! ice-ocean heat flux at bottom surface (W/m^2)
         frzmlt_R4      , & ! freezing/melting potential (W/m^2)
         rside_R4       , & ! fraction of ice that melts laterally
         sst_R4         , & ! sea surface temperature (C)
         Tf_R4          , & ! freezing temperature (C)
         sss_R4         , & ! sea surface salinity (ppt)
         meltt_R4       , & ! top ice melt             (m/step-->cm/day)
         melts_R4       , & ! snow melt                (m/step-->cm/day)
         meltsliq_R4    , & ! snow melt mass           (kg/m^2/step-->kg/m^2/day)
         meltb_R4       , & ! basal ice melt           (m/step-->cm/day)
         meltl_R4       , & ! lateral ice melt         (m/step-->cm/day)
         mlt_onset_R4   , & ! day of year that sfc melting begins
         frz_onset_R4       ! day of year that freezing begins (congel or frazil)

      real (kind=real_kind), dimension(:), intent(inout) :: &
         aicen_init_R4  , & ! fractional area of ice
         vicen_init_R4  , & ! volume per unit area of ice (m)
         vsnon_init_R4  , & ! volume per unit area of snow (m)
         aicen_R4       , & ! concentration of ice
         vicen_R4       , & ! volume per unit area of ice          (m)
         vsnon_R4       , & ! volume per unit area of snow         (m)
         Tsfc_R4        , & ! ice/snow surface temperature, Tsfcn
         alvl_R4        , & ! level ice area fraction
         vlvl_R4        , & ! level ice volume fraction
         apnd_R4        , & ! melt pond area fraction
         hpnd_R4        , & ! melt pond depth (m)
         ipnd_R4        , & ! melt pond refrozen lid thickness (m)
         iage_R4        , & ! volume-weighted ice age
         FY_R4          , & ! area-weighted first-year ice area
         fsurfn_R4      , & ! net flux to top surface, excluding fcondtop
         fcondtopn_R4   , & ! downward cond flux at top surface (W m-2)
         flatn_R4       , & ! latent heat flux (W m-2)
         fsensn_R4      , & ! sensible heat flux (W m-2)
         fsurfn_f_R4    , & ! net flux to top surface, excluding fcondtop
         fcondtopn_f_R4 , & ! downward cond flux at top surface (W m-2)
         flatn_f_R4     , & ! latent heat flux (W m-2)
         fsensn_f_R4    , & ! sensible heat flux (W m-2)
         fswsfcn_R4     , & ! SW absorbed at ice/snow surface (W m-2)
         fswthrun_R4    , & ! SW through ice to ocean            (W/m^2)
         fswintn_R4     , & ! SW absorbed in ice interior, below surface (W m-2)
         faero_atm_R4   , & ! aerosol deposition rate (kg/m^2 s)
         faero_ocn_R4   , & ! aerosol flux to ocean  (kg/m^2/s)
         dhsn_R4        , & ! depth difference for snow on sea ice and pond ice
         ffracn_R4      , & ! fraction of fsurfn used to melt ipond
         meltsn_R4      , & ! snow melt                       (m)
         meltsliqn_R4   , & ! snow melt mass                  (kg/m^2)
         melttn_R4      , & ! top ice melt                    (m)
         meltbn_R4      , & ! bottom ice melt                 (m)
         congeln_R4     , & ! congelation ice growth          (m)
         snoicen_R4     , & ! snow-ice growth                 (m)
         dsnown_R4          ! change in snow thickness (m/step-->cm/day)



      real (kind=real_kind), dimension(:,:), intent(inout) :: &
         zqsn_R4        , & ! snow layer enthalpy (J m-3)
         zqin_R4        , & ! ice layer enthalpy (J m-3)
         zSin_R4        , & ! internal ice layer salinities
         smice_R4       , & ! ice mass tracer in snow (kg/m^3)
         smliq_R4       , & ! liquid water mass tracer in snow (kg/m^3)
         Sswabsn_R4     , & ! SW radiation absorbed in snow layers (W m-2)
         Iswabsn_R4     , & ! SW radiation absorbed in ice layers (W m-2)
         rsnw_R4           ! snow grain radius (10^-6 m) in snow layers 

      real (kind=real_kind), dimension(:,:,:), intent(inout) :: &
         aerosno_R4    , &  ! snow aerosol tracer (kg/m^2)
         aeroice_R4         ! ice aerosol tracer (kg/m^2)



      logical (kind=log_kind), intent(out) :: &
         l_stop          ! if true, abort model

      character (len=*), intent(out) :: &
         stop_label      ! abort error message

      ! local variables

      integer (kind=int_kind) :: &
         n               ! category index

      real (kind=dbl_kind) :: &
         worka, workb    ! temporary variables

      ! 2D coupler variables (computed for each category, then aggregated)
      real (kind=dbl_kind) :: &
         fswabsn     , & ! shortwave absorbed by ice          (W/m^2)
         flwoutn     , & ! upward LW at surface               (W/m^2)
         evapn       , & ! flux of vapor, atmos to ice   (kg m-2 s-1)
         freshn      , & ! flux of water, ice to ocean     (kg/m^2/s)
         fsaltn      , & ! flux of salt, ice to ocean      (kg/m^2/s)
         fhocnn      , & ! fbot corrected for leftover energy (W/m^2)
         strairxn    , & ! air/ice zonal  stress,             (N/m^2)
         strairyn    , & ! air/ice meridional stress,         (N/m^2)
         Cdn_atm_ratio_n, & ! drag coefficient ratio
         Trefn       , & ! air tmp reference level                (K)
         Urefn       , & ! air speed reference level            (m/s)
         Qrefn       , & ! air sp hum reference level         (kg/kg)
         Tbot        , & ! ice bottom surface temperature (deg C)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef      , & ! transfer coefficient for latent heat
         rfrac           ! water fraction retained for melt ponds

      real (kind=dbl_kind) :: &
         raice       , & ! 1/aice
         pond            ! water retained in ponds (m)
      !---------------------------------------------------------------
      ! Double precision temp vars
      !---------------------------------------------------------------
      ! in variables
      real (kind=dbl_kind) :: &
         dt          , & ! time step
         uvel        , & ! x-component of velocity (m/s)
         vvel        , & ! y-component of velocity (m/s)
         strax       , & ! wind stress components (N/m^2)
         stray       , & ! 
         yday            ! day of year
      !out and in/out variables must be copied at end of subroutine
      real (kind=dbl_kind) :: &
         aice0       , & ! open water fraction
         aice        , & ! sea ice concentration
         vice        , & ! volume per unit area of ice          (m)
         vsno        , & ! volume per unit area of snow         (m)
         zlvl        , & ! atm level height (m)
         uatm        , & ! wind velocity components (m/s)
         vatm        , &
         wind        , & ! wind speed (m/s)
         potT        , & ! air potential temperature  (K)
         Tair        , & ! air temperature  (K)
         Qa          , & ! specific humidity (kg/kg)
         rhoa        , & ! air density (kg/m^3)
         frain       , & ! rainfall rate (kg/m^2 s)
         fsnow       , & ! snowfall rate (kg/m^2 s)
         fsloss      , & ! blowing snow loss to leads (kg/m^2/s)
         fpond       , & ! fresh water flux to ponds (kg/m^2/s)
         fresh       , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt       , & ! salt flux to ocean (kg/m^2/s)
         fhocn       , & ! net heat flux to ocean (W/m^2)
         fswthru     , & ! shortwave penetrating to ocean (W/m^2)
         fsurf       , & ! net surface heat flux (excluding fcondtop)(W/m^2)
         fcondtop    , & ! top surface conductive flux        (W/m^2)
         fsens       , & ! sensible heat flux (W/m^2)
         flat        , & ! latent heat flux   (W/m^2)
         fswabs      , & ! shortwave flux absorbed in ice and ocean (W/m^2)
         coszen      , & ! cosine solar zenith angle, < 0 for sun below horizon 
         flw         , & ! incoming longwave radiation (W/m^2)
         flwout      , & ! outgoing longwave radiation (W/m^2)
         evap        , & ! evaporative water flux (kg/m^2/s)
         congel      , & ! basal ice growth         (m/step-->cm/day)
         frazil      , & ! frazil ice growth        (m/step-->cm/day)
         snoice      , & ! snow-ice formation       (m/step-->cm/day)
         Tref        , & ! 2m atm reference temperature (K)
         Qref        , & ! 2m atm reference spec humidity (kg/kg)
         Uref        , & ! 10m atm reference wind speed (m/s)
         Cdn_atm     , & ! atm drag coefficient
         Cdn_ocn     , & ! ocn drag coefficient
         hfreebd     , & ! freeboard (m)
         hdraft      , & ! draft of ice + snow column (Stoessel1993)
         hridge      , & ! ridge height
         distrdg     , & ! distance between ridges
         hkeel       , & ! keel depth
         dkeel       , & ! distance between keels
         lfloe       , & ! floe length
         dfloe       , & ! distance between floes
         Cdn_atm_skin, & ! neutral skin drag coefficient
         Cdn_atm_floe, & ! neutral floe edge drag coefficient
         Cdn_atm_pond, & ! neutral pond edge drag coefficient
         Cdn_atm_rdg , & ! neutral ridge drag coefficient
         Cdn_ocn_skin, & ! skin drag coefficient
         Cdn_ocn_floe, & ! floe edge drag coefficient
         Cdn_ocn_keel, & ! keel drag coefficient
         Cdn_atm_ratio,& ! ratio drag atm / neutral drag atm
         strairxT    , & ! stress on ice by air, x-direction
         strairyT    , & ! stress on ice by air, y-direction
         strocnxT    , & ! ice-ocean stress, x-direction
         strocnyT    , & ! ice-ocean stress, y-direction
         fbot        , & ! ice-ocean heat flux at bottom surface (W/m^2)
         frzmlt      , & ! freezing/melting potential (W/m^2)
         rside       , & ! fraction of ice that melts laterally
         sst         , & ! sea surface temperature (C)
         Tf          , & ! freezing temperature (C)
         sss         , & ! sea surface salinity (ppt)
         meltt       , & ! top ice melt             (m/step-->cm/day)
         melts       , & ! snow melt                (m/step-->cm/day)
         meltsliq    , & ! snow melt mass           (kg/m^2/step-->kg/m^2/day)
         meltb       , & ! basal ice melt           (m/step-->cm/day)
         meltl       , & ! lateral ice melt         (m/step-->cm/day)
         mlt_onset   , & ! day of year that sfc melting begins
         frz_onset       ! day of year that freezing begins (congel or frazil)

      real (kind=dbl_kind), dimension(:), allocatable :: &
         aicen_init  , & ! fractional area of ice
         vicen_init  , & ! volume per unit area of ice (m)
         vsnon_init  , & ! volume per unit area of snow (m)
         aicen       , & ! concentration of ice
         vicen       , & ! volume per unit area of ice          (m)
         vsnon       , & ! volume per unit area of snow         (m)
         Tsfc        , & ! ice/snow surface temperature, Tsfcn
         alvl        , & ! level ice area fraction
         vlvl        , & ! level ice volume fraction
         apnd        , & ! melt pond area fraction
         hpnd        , & ! melt pond depth (m)
         ipnd        , & ! melt pond refrozen lid thickness (m)
         iage        , & ! volume-weighted ice age
         FY          , & ! area-weighted first-year ice area
         fsurfn      , & ! net flux to top surface, excluding fcondtop
         fcondtopn   , & ! downward cond flux at top surface (W m-2)
         flatn       , & ! latent heat flux (W m-2)
         fsensn      , & ! sensible heat flux (W m-2)
         fsurfn_f    , & ! net flux to top surface, excluding fcondtop
         fcondtopn_f , & ! downward cond flux at top surface (W m-2)
         flatn_f     , & ! latent heat flux (W m-2)
         fsensn_f    , & ! sensible heat flux (W m-2)
         fswsfcn     , & ! SW absorbed at ice/snow surface (W m-2)
         fswthrun    , & ! SW through ice to ocean            (W/m^2)
         fswintn     , & ! SW absorbed in ice interior, below surface (W m-2)
         faero_atm   , & ! aerosol deposition rate (kg/m^2 s)
         faero_ocn   , & ! aerosol flux to ocean  (kg/m^2/s)
         dhsn        , & ! depth difference for snow on sea ice and pond ice
         ffracn      , & ! fraction of fsurfn used to melt ipond
         meltsn      , & ! snow melt                       (m)
         meltsliqn   , & ! snow melt mass                  (kg/m^2)
         melttn      , & ! top ice melt                    (m)
         meltbn      , & ! bottom ice melt                 (m)
         congeln     , & ! congelation ice growth          (m)
         snoicen     , & ! snow-ice growth                 (m)
         dsnown          ! change in snow thickness (m/step-->cm/day)

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         zqsn        , & ! snow layer enthalpy (J m-3)
         zqin        , & ! ice layer enthalpy (J m-3)
         zSin        , & ! internal ice layer salinities
         smice       , & ! ice mass tracer in snow (kg/m^3)
         smliq       , & ! liquid water mass tracer in snow (kg/m^3)
         Sswabsn     , & ! SW radiation absorbed in snow layers (W m-2)
         Iswabsn     , & ! SW radiation absorbed in ice layers (W m-2)
         rsnw            ! snow grain radius (10^-6 m) in snow layers 
      real (kind=dbl_kind), dimension(:,:,:), allocatable :: &
         aerosno    , &  ! snow aerosol tracer (kg/m^2)
         aeroice         ! ice aerosol tracer (kg/m^2)

      !---------------------------------------------------------------
      ! copy values to doubles
      !---------------------------------------------------------------

         dt            = real(dt_R4           , kind= dbl_kind ) 
         uvel          = real(uvel_R4         , kind= dbl_kind ) 
         vvel          = real(vvel_R4         , kind= dbl_kind ) 
         strax         = real(strax_R4        , kind= dbl_kind ) 
         stray         = real(stray_R4        , kind= dbl_kind ) 
         yday          = real(yday_R4         , kind= dbl_kind )   
        
         aice0         = real(aice0_R4        , kind= dbl_kind ) 
         aice          = real(aice_R4         , kind= dbl_kind ) 
         vice          = real(vice_R4         , kind= dbl_kind ) 
         vsno          = real(vsno_R4         , kind= dbl_kind ) 
         zlvl          = real(zlvl_R4         , kind= dbl_kind ) 
         uatm          = real(uatm_R4         , kind= dbl_kind ) 
         vatm          = real(vatm_R4         , kind= dbl_kind )         
         wind          = real(wind_R4         , kind= dbl_kind ) 
         potT          = real(potT_R4         , kind= dbl_kind ) 
         Tair          = real(Tair_R4         , kind= dbl_kind ) 
         Qa            = real(Qa_R4           , kind= dbl_kind ) 
         rhoa          = real(rhoa_R4         , kind= dbl_kind ) 
         frain         = real(frain_R4        , kind= dbl_kind ) 
         fsnow         = real(fsnow_R4        , kind= dbl_kind ) 
         fsloss        = real(fsloss_R4       , kind= dbl_kind ) 
         fpond         = real(fpond_R4        , kind= dbl_kind ) 
         fresh         = real(fresh_R4        , kind= dbl_kind ) 
         fsalt         = real(fsalt_R4        , kind= dbl_kind ) 
         fhocn         = real(fhocn_R4        , kind= dbl_kind ) 
         fswthru       = real(fswthru_R4      , kind= dbl_kind ) 
         fsurf         = real(fsurf_R4        , kind= dbl_kind ) 
         fcondtop      = real(fcondtop_R4     , kind= dbl_kind ) 
         fsens         = real(fsens_R4        , kind= dbl_kind ) 
         flat          = real(flat_R4         , kind= dbl_kind ) 
         fswabs        = real(fswabs_R4       , kind= dbl_kind ) 
         coszen        = real(coszen_R4       , kind= dbl_kind ) 
         flw           = real(flw_R4          , kind= dbl_kind ) 
         flwout        = real(flwout_R4       , kind= dbl_kind ) 
         evap          = real(evap_R4         , kind= dbl_kind ) 
         congel        = real(congel_R4       , kind= dbl_kind ) 
         frazil        = real(frazil_R4       , kind= dbl_kind ) 
         snoice        = real(snoice_R4       , kind= dbl_kind ) 
         Tref          = real(Tref_R4         , kind= dbl_kind ) 
         Qref          = real(Qref_R4         , kind= dbl_kind ) 
         Uref          = real(Uref_R4         , kind= dbl_kind ) 
         Cdn_atm       = real(Cdn_atm_R4      , kind= dbl_kind ) 
         Cdn_ocn       = real(Cdn_ocn_R4      , kind= dbl_kind ) 
         hfreebd       = real(hfreebd_R4      , kind= dbl_kind ) 
         hdraft        = real(hdraft_R4       , kind= dbl_kind ) 
         hridge        = real(hridge_R4       , kind= dbl_kind ) 
         distrdg       = real(distrdg_R4      , kind= dbl_kind ) 
         hkeel         = real(hkeel_R4        , kind= dbl_kind ) 
         dkeel         = real(dkeel_R4        , kind= dbl_kind ) 
         lfloe         = real(lfloe_R4        , kind= dbl_kind ) 
         dfloe         = real(dfloe_R4        , kind= dbl_kind ) 
         Cdn_atm_skin  = real(Cdn_atm_skin_R4 , kind= dbl_kind ) 
         Cdn_atm_floe  = real(Cdn_atm_floe_R4 , kind= dbl_kind ) 
         Cdn_atm_pond  = real(Cdn_atm_pond_R4 , kind= dbl_kind ) 
         Cdn_atm_rdg   = real(Cdn_atm_rdg_R4  , kind= dbl_kind ) 
         Cdn_ocn_skin  = real(Cdn_ocn_skin_R4 , kind= dbl_kind ) 
         Cdn_ocn_floe  = real(Cdn_ocn_floe_R4 , kind= dbl_kind ) 
         Cdn_ocn_keel  = real(Cdn_ocn_keel_R4 , kind= dbl_kind ) 
         Cdn_atm_ratio = real(Cdn_atm_ratio_R4, kind= dbl_kind )  
         strairxT      = real(strairxT_R4     , kind= dbl_kind ) 
         strairyT      = real(strairyT_R4     , kind= dbl_kind ) 
         strocnxT      = real(strocnxT_R4     , kind= dbl_kind ) 
         strocnyT      = real(strocnyT_R4     , kind= dbl_kind ) 
         fbot          = real(fbot_R4         , kind= dbl_kind ) 
         frzmlt        = real(frzmlt_R4       , kind= dbl_kind ) 
         rside         = real(rside_R4        , kind= dbl_kind ) 
         sst           = real(sst_R4          , kind= dbl_kind ) 
         Tf            = real(Tf_R4           , kind= dbl_kind ) 
         sss           = real(sss_R4          , kind= dbl_kind ) 
         meltt         = real(meltt_R4        , kind= dbl_kind ) 
         melts         = real(melts_R4        , kind= dbl_kind ) 
         meltsliq      = real(meltsliq_R4     , kind= dbl_kind ) 
         meltb         = real(meltb_R4        , kind= dbl_kind ) 
         meltl         = real(meltl_R4        , kind= dbl_kind ) 
         mlt_onset     = real(mlt_onset_R4    , kind= dbl_kind ) 
         frz_onset     = real(frz_onset_R4    , kind= dbl_kind )   

        
        aicen_init   = real(aicen_init_R4   , kind= dbl_kind ) 
        vicen_init   = real(vicen_init_R4   , kind= dbl_kind ) 
        vsnon_init   = real(vsnon_init_R4   , kind= dbl_kind ) 
        aicen        = real(aicen_R4        , kind= dbl_kind ) 
        vicen        = real(vicen_R4        , kind= dbl_kind ) 
        vsnon        = real(vsnon_R4        , kind= dbl_kind ) 
        Tsfc         = real(Tsfc_R4         , kind= dbl_kind ) 
        alvl         = real(alvl_R4         , kind= dbl_kind ) 
        vlvl         = real(vlvl_R4         , kind= dbl_kind ) 
        apnd         = real(apnd_R4         , kind= dbl_kind ) 
        hpnd         = real(hpnd_R4         , kind= dbl_kind ) 
        ipnd         = real(ipnd_R4         , kind= dbl_kind ) 
        iage         = real(iage_R4         , kind= dbl_kind ) 
        FY           = real(FY_R4           , kind= dbl_kind ) 
        fsurfn       = real(fsurfn_R4       , kind= dbl_kind ) 
        fcondtopn    = real(fcondtopn_R4    , kind= dbl_kind ) 
        flatn        = real(flatn_R4        , kind= dbl_kind ) 
        fsensn       = real(fsensn_R4       , kind= dbl_kind ) 
        fsurfn_f     = real(fsurfn_f_R4     , kind= dbl_kind ) 
        fcondtopn_f  = real(fcondtopn_f_R4  , kind= dbl_kind ) 
        flatn_f      = real(flatn_f_R4      , kind= dbl_kind ) 
        fsensn_f     = real(fsensn_f_R4     , kind= dbl_kind ) 
        fswsfcn      = real(fswsfcn_R4      , kind= dbl_kind ) 
        fswthrun     = real(fswthrun_R4     , kind= dbl_kind ) 
        fswintn      = real(fswintn_R4      , kind= dbl_kind ) 
        faero_atm    = real(faero_atm_R4    , kind= dbl_kind ) 
        faero_ocn    = real(faero_ocn_R4    , kind= dbl_kind ) 
        dhsn         = real(dhsn_R4         , kind= dbl_kind ) 
        ffracn       = real(ffracn_R4       , kind= dbl_kind ) 
        meltsn       = real(meltsn_R4       , kind= dbl_kind ) 
        meltsliqn    = real(meltsliqn_R4    , kind= dbl_kind ) 
        melttn       = real(melttn_R4       , kind= dbl_kind ) 
        meltbn       = real(meltbn_R4       , kind= dbl_kind ) 
        congeln      = real(congeln_R4      , kind= dbl_kind ) 
        snoicen      = real(snoicen_R4      , kind= dbl_kind ) 
        dsnown       = real(dsnown_R4       , kind= dbl_kind )   
        zqsn         = real(zqsn_R4         , kind= dbl_kind ) 
        zqin         = real(zqin_R4         , kind= dbl_kind ) 
        zSin         = real(zSin_R4         , kind= dbl_kind ) 
        smice        = real(smice_R4        , kind= dbl_kind ) 
        smliq        = real(smliq_R4        , kind= dbl_kind ) 
        Sswabsn      = real(Sswabsn_R4      , kind= dbl_kind ) 
        Iswabsn      = real(Iswabsn_R4      , kind= dbl_kind ) 
        rsnw         = real(rsnw_R4         , kind= dbl_kind )   
        aerosno      = real(aerosno_R4      , kind= dbl_kind )
        aeroice      = real(aeroice_R4      , kind= dbl_kind )  
      !---------------------------------------------------------------
      ! Initialize rate of snow loss to leads
      !---------------------------------------------------------------

      fsloss = fsnow*aice0

      !---------------------------------------------------------------
      ! 30% rule for snow redistribution: precip factor
      !---------------------------------------------------------------

      if (trim(snwredist) == '30percent') then
         worka = c0      
         do n = 1, ncat
            worka = worka + alvl(n)
         enddo
         worka  = worka * snwlvlfac/(c1+snwlvlfac)
         fsloss = fsloss + fsnow*(c1-worka)
         fsnow  =          fsnow*    worka
      endif ! snwredist

      !-----------------------------------------------------------------
      ! Adjust frzmlt to account for ice-ocean heat fluxes since last
      !  call to coupler.
      ! Compute lateral and bottom heat fluxes.
      !-----------------------------------------------------------------

      call frzmlt_bottom_lateral (dt,        ncat,      &
                                  nilyr,     nslyr,     &
                                  aice,      frzmlt,    &
                                  vicen,     vsnon,     &
                                  zqin,      zqsn,      &
                                  sst,       Tf,        &
                                  ustar_min,            &
                                  fbot_xfer_type,       &
                                  strocnxT,  strocnyT,  &
                                  Tbot,      fbot,      &
                                  rside,     Cdn_ocn)
      
      !-----------------------------------------------------------------
      ! Update the neutral drag coefficients to account for form drag
      ! Oceanic and atmospheric drag coefficients
      !-----------------------------------------------------------------

      if (formdrag) then
         call neutral_drag_coeffs (apnd         , &
                                   hpnd        , ipnd         , &
                                   alvl        , vlvl         , &
                                   aice        , vice,          &
                                   vsno        , aicen        , &
                                   vicen       , vsnon        , &
                                   Cdn_ocn     , Cdn_ocn_skin, &
                                   Cdn_ocn_floe, Cdn_ocn_keel, &
                                   Cdn_atm     , Cdn_atm_skin, &
                                   Cdn_atm_floe, Cdn_atm_pond, &
                                   Cdn_atm_rdg , hfreebd     , &
                                   hdraft      , hridge      , &
                                   distrdg     , hkeel       , &
                                   dkeel       , lfloe       , &
                                   dfloe       , ncat)
      endif

      do n = 1, ncat

         meltsn (n) = c0
         meltsliqn(n) = c0
         melttn (n) = c0
         meltbn (n) = c0
         congeln(n) = c0
         snoicen(n) = c0
         dsnown (n) = c0

         Trefn  = c0
         Qrefn  = c0
         Urefn  = c0
         lhcoef = c0
         shcoef = c0
         worka  = c0
         workb  = c0

         if (aicen_init(n) > puny) then

            if (calc_Tsfc .or. calc_strair) then 

      !-----------------------------------------------------------------
      ! Atmosphere boundary layer calculation; compute coefficients
      ! for sensible and latent heat fluxes.
      !
      ! NOTE: The wind stress is computed here for later use if 
      !       calc_strair = .true.   Otherwise, the wind stress
      !       components are set to the data values.
      !-----------------------------------------------------------------

               call colpkg_atm_boundary_double( 'ice',                  &
                                        Tsfc(n),  potT,          &
                                        uatm,     vatm,          &
                                        wind,     zlvl,          &
                                        Qa,       rhoa,          &
                                        strairxn, strairyn,      &
                                        Trefn,    Qrefn,         &
                                        worka,    workb,         &
                                        lhcoef,   shcoef,        &
                                        Cdn_atm,                 &
                                        Cdn_atm_ratio_n,         &
                                        uvel,     vvel,          &
                                        Uref=Urefn)

            endif   ! calc_Tsfc or calc_strair

            if (.not.(calc_strair)) then
#ifndef CICE_IN_NEMO
               ! Set to data values (on T points)
               strairxn = strax
               strairyn = stray
#else
               ! NEMO wind stress is supplied on u grid, multipied 
               ! by ice concentration and set directly in evp, so
               ! strairxT/yT = 0. Zero u-components here for safety.
               strairxn = c0
               strairyn = c0
#endif
            endif

      !-----------------------------------------------------------------
      ! Update ice age
      ! This is further adjusted for freezing in the thermodynamics.
      ! Melting does not alter the ice age.
      !-----------------------------------------------------------------

            if (tr_iage) call increment_age (dt, iage(n))
            if (tr_FY)   call update_FYarea (dt,               &
                                             lmask_n, lmask_s, &
                                             yday,    FY(n))

      !-----------------------------------------------------------------
      ! Vertical thermodynamics: Heat conduction, growth and melting.
      !----------------------------------------------------------------- 

            if (.not.(calc_Tsfc)) then

               ! If not calculating surface temperature and fluxes, set 
               ! surface fluxes (flatn, fsurfn, and fcondtopn) to be used 
               ! in thickness_changes
 
               ! hadgem routine sets fluxes to default values in ice-only mode
               call set_sfcflux(aicen      (n),                 &
                                flatn_f    (n), fsensn_f   (n), &
                                fcondtopn_f(n),                 &
                                fsurfn_f   (n),                 &
                                flatn      (n), fsensn     (n), &
                                fsurfn     (n),                 &
                                fcondtopn  (n))
            endif

            call thermo_vertical(nilyr,        nslyr,        &
                                 dt,           aicen    (n), &
                                 vicen    (n), vsnon    (n), &
                                 Tsfc     (n), zSin   (:,n), &
                                 zqin   (:,n), zqsn   (:,n), &
                                 smice  (:,n), smliq  (:,n), &
                                 tr_snow,      apnd     (n), &
                                 hpnd     (n), iage     (n), &
                                 tr_pond_topo,               &
                                 flw,          potT,         &
                                 Qa,           rhoa,         &
                                 fsnow,        fpond,        &
                                 fbot,         Tbot,         &
                                 sss,          rsnw   (:,n), &
                                 lhcoef,       shcoef,       &
                                 fswsfcn  (n), fswintn  (n), &
                                 Sswabsn(:,n), Iswabsn(:,n), &
                                 fsurfn   (n), fcondtopn(n), &
                                 fsensn   (n), flatn    (n), &
                                 flwoutn,      evapn,        &
                                 freshn,       fsaltn,       &
                                 fhocnn,       frain,        &
                                 melttn   (n), meltsn   (n), &
                                 meltbn   (n), meltsliqn(n), &
                                 congeln  (n), snoicen  (n), &
                                 mlt_onset,    frz_onset,    &
                                 yday,         dsnown   (n), &
                                 tr_rsnw,                    &
                                 l_stop,       stop_label,   &
                                 prescribed_ice)
               
            if (l_stop) then
               stop_label = 'ice: Vertical thermo error: '//trim(stop_label)
               return
            endif
               
      !-----------------------------------------------------------------
      ! Total absorbed shortwave radiation
      !-----------------------------------------------------------------

            fswabsn = fswsfcn(n) + fswintn(n) + fswthrun(n)

      !-----------------------------------------------------------------
      ! Aerosol update
      !-----------------------------------------------------------------

            if (tr_aero) then
               call update_aerosol (dt,                             &
                                    nilyr, nslyr, n_aero,           &
                                    melttn     (n), meltsn     (n), &
                                    meltbn     (n), congeln    (n), &
                                    snoicen    (n), fsnow,          &
                                    aerosno(:,:,n), aeroice(:,:,n), &
                                    aicen_init (n), vicen_init (n), &
                                    vsnon_init (n),                 &
                                    vicen      (n), vsnon      (n), &
                                    aicen      (n),                 &
                                    faero_atm     ,  faero_ocn)
            endif

         endif   ! aicen_init

      !-----------------------------------------------------------------
      ! Transport liquid water in snow between layers and 
      ! compute the meltpond contribution
      !-----------------------------------------------------------------

       if (use_smliq_pnd) then
         call drain_snow (dt,            nslyr,        &
                          vsnon   (n) ,  aicen    (n), &
                          smice  (:,n),  smliq  (:,n), &
                          meltsliqn(n))
       endif

      !-----------------------------------------------------------------
      ! Melt ponds
      ! If using tr_pond_cesm, the full calculation is performed here.
      ! If using tr_pond_topo, the rest of the calculation is done after
      ! the surface fluxes are merged, below.
      !-----------------------------------------------------------------

         !call ice_timer_start(timer_ponds)
         if (tr_pond) then
            if (tr_pond_cesm) then
               rfrac = rfracmin + (rfracmax-rfracmin) * aicen(n) 
               call compute_ponds_cesm(dt,        hi_min,    &
                                       pndaspect, rfrac,     &
                                       melttn(n), meltsn(n), &
                                       frain,                &
                                       aicen (n), vicen (n), &
                                       vsnon (n), Tsfc  (n), &
                                       apnd  (n), hpnd  (n), &
                                       meltsliqn(n),         &
                                       use_smliq_pnd)
                  
            elseif (tr_pond_lvl) then
               rfrac = rfracmin + (rfracmax-rfracmin) * aicen(n)
               call compute_ponds_lvl(dt,        nilyr,     &
                                      ktherm,               &
                                      hi_min,               &
                                      dpscale,   frzpnd,    &
                                      pndaspect, rfrac,     &
                                      melttn(n), meltsn(n), &
                                      frain,     Tair,      &
                                      fsurfn(n),            &
                                      dhsn  (n), ffracn(n), &
                                      aicen (n), vicen (n), &
                                      vsnon (n),            &
                                      zqin(:,n), zSin(:,n), &
                                      Tsfc  (n), alvl  (n), &
                                      apnd  (n), hpnd  (n), &
                                      ipnd  (n),            &
                                      meltsliqn(n),         &
                                      use_smliq_pnd)
                  
            elseif (tr_pond_topo) then
               if (aicen_init(n) > puny) then
                     
                  ! collect liquid water in ponds
                  ! assume salt still runs off
                  rfrac = rfracmin + (rfracmax-rfracmin) * aicen(n)
                  if (use_smliq_pnd) then
                    pond = rfrac/rhofresh * (melttn(n)*rhoi &
                         +                   meltsliqn(n))
                  else
                    pond = rfrac/rhofresh * (melttn(n)*rhoi &
                         +                   meltsn(n)*rhos &
                         +                   frain *dt)
                  endif
                  ! if pond does not exist, create new pond over full ice area
                  ! otherwise increase pond depth without changing pond area
                  if (apnd(n) < puny) then
                     hpnd(n) = c0
                     apnd(n) = c1
                  endif
                  hpnd(n) = (pond + hpnd(n)*apnd(n)) / apnd(n)
                  fpond = fpond + pond * aicen(n) ! m
               endif ! aicen_init
            endif

         endif ! tr_pond
         !call ice_timer_stop(timer_ponds)

      !-----------------------------------------------------------------
      ! Increment area-weighted fluxes.
      !-----------------------------------------------------------------

         if (aicen_init(n) > puny) &
            call merge_fluxes (aicen_init(n),            &
                               flw,        coszen,       & 
                               strairxn,   strairyn,     &
                               Cdn_atm_ratio_n,          &
                               fsurfn(n),  fcondtopn(n), &
                               fsensn(n),  flatn(n),     &
                               fswabsn,    flwoutn,      &
                               evapn,                    &
                               Trefn,      Qrefn,        &
                               freshn,     fsaltn,       &
                               fhocnn,     fswthrun(n),  &
                               strairxT,   strairyT,     &
                               Cdn_atm_ratio,            &
                               fsurf,      fcondtop,     &
                               fsens,      flat,         &
                               fswabs,     flwout,       &
                               evap,                     &
                               Tref,       Qref,         &
                               fresh,      fsalt,        &
                               fhocn,      fswthru,      &
                               melttn (n), meltsn(n),    &
                               meltbn (n), congeln(n),   &
                               snoicen(n), meltsliqn(n), &
                               meltt,      melts,        &
                               meltb,      congel,       &
                               snoice,     meltsliq,     &
                               Uref,       Urefn)

      enddo                  ! ncat

      !-----------------------------------------------------------------
      ! Calculate ponds from the topographic scheme
      !-----------------------------------------------------------------
      !call ice_timer_start(timer_ponds)
      if (tr_pond_topo) then
         call compute_ponds_topo(dt,       ncat,      nilyr,     &
                                 ktherm,   heat_capacity,        &
                                 aice,     aicen,                &
                                 vice,     vicen,                &
                                 vsno,     vsnon,                &
                                 potT,     meltt,                &
                                 fsurf,    fpond,                &
                                 Tsfc,     Tf,                   &
                                 zqin,     zSin,                 &
                                 apnd,     hpnd,      ipnd,      &
                                 l_stop,   stop_label)
      endif
      !call ice_timer_stop(timer_ponds)

      !-----------------------------------------------------------------
      ! Copy values back to og variables
      !-----------------------------------------------------------------
        
        
      aice0_R4         = real(aice0        , kind= real_kind ) 
      aice_R4          = real(aice         , kind= real_kind ) 
      vice_R4          = real(vice         , kind= real_kind ) 
      vsno_R4          = real(vsno         , kind= real_kind ) 
      zlvl_R4           = real(zlvl         , kind= real_kind ) 
      uatm_R4          = real(uatm         , kind= real_kind ) 
      vatm_R4          = real(vatm         , kind= real_kind )         
      wind_R4          = real(wind         , kind= real_kind ) 
      potT_R4          = real(potT         , kind= real_kind ) 
      Tair_R4          = real(Tair         , kind= real_kind ) 
      Qa_R4            = real(Qa           , kind= real_kind ) 
      rhoa_R4          = real(rhoa         , kind= real_kind ) 
      frain_R4         = real(frain        , kind= real_kind ) 
      fsnow_R4         = real(fsnow        , kind= real_kind ) 
      fsloss_R4        = real(fsloss       , kind= real_kind ) 
      fpond_R4         = real(fpond        , kind= real_kind ) 
      fresh_R4         = real(fresh        , kind= real_kind ) 
      fsalt_R4          = real(fsalt        , kind= real_kind ) 
      fhocn_R4         = real(fhocn        , kind= real_kind ) 
      fswthru_R4       = real(fswthru      , kind= real_kind ) 
      fsurf_R4         = real(fsurf        , kind= real_kind ) 
      fcondtop_R4      = real(fcondtop     , kind= real_kind ) 
      fsens_R4         = real(fsens        , kind= real_kind ) 
      flat_R4           = real(flat         , kind= real_kind ) 
      fswabs_R4        = real(fswabs       , kind= real_kind ) 
      coszen_R4        = real(coszen       , kind= real_kind ) 
      flw_R4           = real(flw          , kind= real_kind ) 
      flwout_R4        = real(flwout       , kind= real_kind ) 
      evap_R4          = real(evap         , kind= real_kind ) 
      congel_R4        = real(congel       , kind= real_kind ) 
      frazil_R4        = real(frazil       , kind= real_kind ) 
      snoice_R4        = real(snoice       , kind= real_kind ) 
      Tref_R4          = real(Tref         , kind= real_kind ) 
      Qref_R4          = real(Qref         , kind= real_kind ) 
      Uref_R4          = real(Uref         , kind= real_kind ) 
      Cdn_atm_R4       = real(Cdn_atm      , kind= real_kind ) 
      Cdn_ocn_R4       = real(Cdn_ocn      , kind= real_kind ) 
      hfreebd_R4       = real(hfreebd      , kind= real_kind ) 
      hdraft_R4        = real(hdraft       , kind= real_kind ) 
      hridge_R4        = real(hridge       , kind= real_kind ) 
      distrdg_R4       = real(distrdg      , kind= real_kind ) 
      hkeel_R4         = real(hkeel        , kind= real_kind ) 
      dkeel_R4         = real(dkeel        , kind= real_kind ) 
      lfloe_R4         = real(lfloe        , kind= real_kind ) 
      dfloe_R4         = real(dfloe        , kind= real_kind ) 
      Cdn_atm_skin_R4  = real(Cdn_atm_skin , kind= real_kind ) 
      Cdn_atm_floe_R4  = real(Cdn_atm_floe , kind= real_kind ) 
      Cdn_atm_pond_R4  = real(Cdn_atm_pond , kind= real_kind ) 
      Cdn_atm_rdg_R4   = real(Cdn_atm_rdg  , kind= real_kind ) 
      Cdn_ocn_skin_R4  = real(Cdn_ocn_skin , kind= real_kind ) 
      Cdn_ocn_floe_R4  = real(Cdn_ocn_floe , kind= real_kind ) 
      Cdn_ocn_keel_R4  = real(Cdn_ocn_keel , kind= real_kind ) 
      Cdn_atm_ratio_R4 = real(Cdn_atm_ratio, kind= real_kind )  
      strairxT_R4      = real(strairxT     , kind= real_kind ) 
      strairyT_R4      = real(strairyT     , kind= real_kind ) 
      strocnxT_R4      = real(strocnxT     , kind= real_kind ) 
      strocnyT_R4      = real(strocnyT     , kind= real_kind ) 
      fbot_R4           = real(fbot         , kind= real_kind ) 
      frzmlt_R4        = real(frzmlt       , kind= real_kind ) 
      rside_R4         = real(rside        , kind= real_kind ) 
      sst_R4           = real(sst          , kind= real_kind ) 
      Tf_R4            = real(Tf           , kind= real_kind ) 
      sss_R4           = real(sss          , kind= real_kind ) 
      meltt_R4          = real(meltt        , kind= real_kind ) 
      melts_R4         = real(melts        , kind= real_kind ) 
      meltsliq_R4      = real(meltsliq     , kind= real_kind ) 
      meltb_R4         = real(meltb        , kind= real_kind ) 
      meltl_R4         = real(meltl        , kind= real_kind ) 
      mlt_onset_R4     = real(mlt_onset    , kind= real_kind ) 
      frz_onset_R4     = real(frz_onset    , kind= real_kind )   

      
      aicen_init_R4    = real(aicen_init   , kind= real_kind ) 
      vicen_init_R4    = real(vicen_init   , kind= real_kind ) 
      vsnon_init_R4    = real(vsnon_init   , kind= real_kind ) 
      aicen_R4         = real(aicen        , kind= real_kind ) 
      vicen_R4         = real(vicen        , kind= real_kind ) 
      vsnon_R4         = real(vsnon        , kind= real_kind ) 
      Tsfc_R4          = real(Tsfc         , kind= real_kind ) 
      alvl_R4           = real(alvl         , kind= real_kind ) 
      vlvl_R4           = real(vlvl         , kind= real_kind ) 
      apnd_R4          = real(apnd         , kind= real_kind ) 
      hpnd_R4          = real(hpnd         , kind= real_kind ) 
      ipnd_R4          = real(ipnd         , kind= real_kind ) 
      iage_R4          = real(iage         , kind= real_kind ) 
      FY_R4            = real(FY           , kind= real_kind ) 
      fsurfn_R4        = real(fsurfn       , kind= real_kind ) 
      fcondtopn_R4     = real(fcondtopn    , kind= real_kind ) 
      flatn_R4         = real(flatn        , kind= real_kind ) 
      fsensn_R4        = real(fsensn       , kind= real_kind ) 
      fsurfn_f_R4      = real(fsurfn_f     , kind= real_kind ) 
      fcondtopn_f_R4   = real(fcondtopn_f  , kind= real_kind ) 
      flatn_f_R4       = real(flatn_f      , kind= real_kind ) 
      fsensn_f_R4      = real(fsensn_f     , kind= real_kind ) 
      fswsfcn_R4       = real(fswsfcn      , kind= real_kind ) 
      fswthrun_R4      = real(fswthrun     , kind= real_kind ) 
      fswintn_R4       = real(fswintn      , kind= real_kind ) 
      faero_atm_R4     = real(faero_atm    , kind= real_kind ) 
      faero_ocn_R4     = real(faero_ocn    , kind= real_kind ) 
      dhsn_R4          = real(dhsn         , kind= real_kind ) 
      ffracn_R4        = real(ffracn       , kind= real_kind ) 
      meltsn_R4        = real(meltsn       , kind= real_kind ) 
      meltsliqn_R4     = real(meltsliqn    , kind= real_kind ) 
      melttn_R4        = real(melttn       , kind= real_kind ) 
      meltbn_R4        = real(meltbn       , kind= real_kind ) 
      congeln_R4       = real(congeln      , kind= real_kind ) 
      snoicen_R4       = real(snoicen      , kind= real_kind ) 
      dsnown_R4        = real(dsnown       , kind= real_kind )   

      
      zqsn_R4          = real(zqsn         , kind= real_kind ) 
      zqin_R4          = real(zqin         , kind= real_kind ) 
      zSin_R4          = real(zSin         , kind= real_kind ) 
      smice_R4         = real(smice        , kind= real_kind ) 
      smliq_R4         = real(smliq        , kind= real_kind ) 
      Sswabsn_R4       = real(Sswabsn      , kind= real_kind ) 
      Iswabsn_R4       = real(Iswabsn      , kind= real_kind ) 
      rsnw_R4          = real(rsnw         , kind= real_kind )   
      
      aerosno_R4       = real(aerosno      , kind= real_kind )
      aeroice_R4       = real(aeroice      , kind= real_kind )     
      
      deallocate(aicen_init)
      deallocate(vicen_init)
      deallocate(vsnon_init)
      deallocate(aicen)
      deallocate(vicen)
      deallocate(vsnon)
      deallocate(Tsfc)
      deallocate(alvl)
      deallocate(vlvl)
      deallocate(apnd)
      deallocate(hpnd)
      deallocate(ipnd)
      deallocate(iage)
      deallocate(FY)
      deallocate(fsurfn)
      deallocate(fcondtopn)
      deallocate(flatn)
      deallocate(fsensn)
      deallocate(fsurfn_f)
      deallocate(fcondtopn_f)
      deallocate(flatn_f)
      deallocate(fsensn_f)
      deallocate(fswsfcn)
      deallocate(fswthrun)
      deallocate(fswintn)
      deallocate(faero_atm)
      deallocate(faero_ocn)
      deallocate(dhsn)
      deallocate(ffracn)
      deallocate(meltsn)
      deallocate(meltsliqn)
      deallocate(melttn)
      deallocate(meltbn)
      deallocate(congeln)
      deallocate(snoicen)
      deallocate(dsnown)
      deallocate(zqsn)
      deallocate(zqin)
      deallocate(zSin)
      deallocate(smice)
      deallocate(smliq)
      deallocate(Sswabsn)
      deallocate(Iswabsn)
      deallocate(rsnw)
      deallocate(aerosno)
      deallocate(aeroice)
      end subroutine colpkg_step_therm1

!=======================================================================
! Driver for thermodynamic changes not needed for coupling:
! transport in thickness space, lateral growth and melting.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine colpkg_step_therm2 (dt_R4, ncat, n_aero, nbtrcr,    &
                                     nilyr,           nslyr,         &
                                     hin_max_R4,      nblyr,         &
                                     aicen_R4,                       &
                                     vicen_R4,        vsnon_R4,      &
                                     aicen_init_R4,   vicen_init_R4, &
                                     trcrn_R4,                       &
                                     aice0_R4,        aice_R4,       &
                                     trcr_depend,                    &
                                     trcr_base_R4,    n_trcr_strata, &
                                     nt_strata,                      &
                                     Tf_R4,           sss_R4,        &
                                     salinz_R4,                      &
                                     rside_R4,        meltl_R4,      &
                                     frzmlt_R4,       frazil_R4,     &
                                     frain_R4,        fpond_R4,      &
                                     fresh_R4,        fsalt_R4,      &
                                     fhocn_R4,        update_ocn_f,  &
                                     bgrid_R4,        cgrid_R4,      &
                                     igrid_R4,        faero_ocn_R4,     &
                                     first_ice,       fzsal_R4,      &
                                     flux_bio_R4,     ocean_bio_R4,  &
                                     l_stop,          stop_label,    &
                                     frazil_diag_R4,                 &
                                     frz_onset_R4,    yday_R4)


      use ice_constants_colpkg, only: puny, c0
      use ice_itd, only: aggregate_area, reduce_area, cleanup_itd
      use ice_therm_itd, only: linear_itd, add_new_ice, lateral_melt
      use ice_colpkg_tracers, only: ntrcr, tr_aero, tr_pond_topo, tr_brine, nt_fbri, bio_index

      integer (kind=int_kind), intent(in) :: &
         ncat     , & ! number of thickness categories
         nbtrcr   , & ! number of zbgc tracers
         nblyr    , & ! number of bio layers
         nilyr    , & ! number of ice layers
         nslyr    , & ! number of snow layers
         n_aero       ! number of aerosol tracers

      logical (kind=log_kind), intent(in) :: &
         update_ocn_f     ! if true, update fresh water and salt fluxes

      real (kind=real_kind), dimension(0:ncat), intent(in) :: &
         hin_max_R4      ! category boundaries (m)

      real (kind=real_kind), intent(in) :: &
         dt_R4       , & ! time step
         Tf_R4       , & ! freezing temperature (C)
         sss_R4      , & ! sea surface salinity (ppt)
         rside_R4    , & ! fraction of ice that melts laterally
         frzmlt_R4       ! freezing/melting potential (W/m^2)

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=real_kind), dimension (:,:), intent(in) :: &
         trcr_base_R4      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      real (kind=real_kind), dimension (nblyr+2), intent(in) :: &
         bgrid_R4              ! biology nondimensional vertical grid points

      real (kind=real_kind), dimension (nblyr+1), intent(in) :: &
         igrid_R4              ! biology vertical interface points
 
      real (kind=real_kind), dimension (nilyr+1), intent(in) :: &
         cgrid_R4              ! CICE vertical coordinate   

      real (kind=real_kind), dimension(:), intent(in) :: &
         salinz_R4   , & ! initial salinity profile
         ocean_bio_R4    ! ocean concentration of biological tracer

      real (kind=real_kind), intent(inout) :: &
         aice_R4     , & ! sea ice concentration
         aice0_R4    , & ! concentration of open water
         frain_R4    , & ! rainfall rate (kg/m^2 s)
         fpond_R4    , & ! fresh water flux to ponds (kg/m^2/s)
         fresh_R4    , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt_R4    , & ! salt flux to ocean (kg/m^2/s)
         fhocn_R4    , & ! net heat flux to ocean (W/m^2)
         fzsal_R4    , & ! salt flux to ocean from zsalinity (kg/m^2/s)
         meltl_R4    , & ! lateral ice melt         (m/step-->cm/day)
         frazil_R4   , & ! frazil ice growth        (m/step-->cm/day)
         frazil_diag_R4  ! frazil ice growth diagnostic (m/step-->cm/day)

      real (kind=real_kind), dimension(:), intent(inout) :: &
         aicen_init_R4,& ! initial concentration of ice
         vicen_init_R4,& ! initial volume per unit area of ice          (m)
         aicen_R4    , & ! concentration of ice
         vicen_R4    , & ! volume per unit area of ice          (m)
         vsnon_R4    , & ! volume per unit area of snow         (m)
         faero_ocn_R4, & ! aerosol flux to ocean  (kg/m^2/s)
         flux_bio_R4     ! all bio fluxes to ocean

      real (kind=real_kind), dimension(:,:), intent(inout) :: &
         trcrn_R4        ! tracers
 
      logical (kind=log_kind), dimension(:), intent(inout) :: &
         first_ice      ! true until ice forms

      logical (kind=log_kind), intent(out) :: &
         l_stop         ! if true, abort model

      character (len=*), intent(out) :: stop_label

      real (kind=real_kind), intent(inout), optional :: &
         frz_onset_R4    ! day of year that freezing begins (congel or frazil)

      real (kind=real_kind), intent(in), optional :: &
         yday_R4         ! day of year


      !-----------------------------------------------------------------
      ! Local double variable copies 
      !-----------------------------------------------------------------


      !intent in variables 
      real (kind=dbl_kind), dimension(0:ncat) :: &
         hin_max      ! category boundaries (m)

      real (kind=dbl_kind) :: &
         dt       , & ! time step
         Tf       , & ! freezing temperature (C)
         sss      , & ! sea surface salinity (ppt)
         rside    , & ! fraction of ice that melts laterally
         frzmlt       ! freezing/melting potential (W/m^2)

      real (kind=dbl_kind), dimension (:,:), allocatable :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      real (kind=dbl_kind), dimension (nblyr+2) :: &
         bgrid              ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1) :: &
         igrid              ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1) :: &
         cgrid              ! CICE vertical coordinate   

      real (kind=dbl_kind), dimension(:), allocatable :: &
         salinz   , & ! initial salinity profile
         ocean_bio    ! ocean concentration of biological tracer
      real (kind=dbl_kind) :: yday         ! day of year
         !nothing
      !intent in/out variables
      real (kind=dbl_kind) :: &
         aice     , & ! sea ice concentration
         aice0    , & ! concentration of open water
         frain    , & ! rainfall rate (kg/m^2 s)
         fpond    , & ! fresh water flux to ponds (kg/m^2/s)
         fresh    , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt    , & ! salt flux to ocean (kg/m^2/s)
         fhocn    , & ! net heat flux to ocean (W/m^2)
         fzsal    , & ! salt flux to ocean from zsalinity (kg/m^2/s)
         meltl    , & ! lateral ice melt         (m/step-->cm/day)
         frazil   , & ! frazil ice growth        (m/step-->cm/day)
         frazil_diag  ! frazil ice growth diagnostic (m/step-->cm/day)

      real (kind=dbl_kind), dimension(:), allocatable :: &
         aicen_init,& ! initial concentration of ice
         vicen_init,& ! initial volume per unit area of ice          (m)
         aicen    , & ! concentration of ice
         vicen    , & ! volume per unit area of ice          (m)
         vsnon    , & ! volume per unit area of snow         (m)
         faero_ocn, & ! aerosol flux to ocean  (kg/m^2/s)
         flux_bio     ! all bio fluxes to ocean

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         trcrn        ! tracers
      real (kind=dbl_kind) :: &
         frz_onset    ! day of year that freezing begins (congel or frazil)
      hin_max  = real(      hin_max_R4, kind = dbl_kind)
      dt       = real(      dt_R4     , kind = dbl_kind)
      Tf       = real(      Tf_R4     , kind = dbl_kind)
      sss      = real(      sss_R4    , kind = dbl_kind)
      rside    = real(      rside_R4  , kind = dbl_kind)
      frzmlt       = real(  frzmlt_R4     , kind = dbl_kind)
     trcr_base     = real(  trcr_base_R4  , kind = dbl_kind)
      bgrid        = real(  bgrid_R4      , kind = dbl_kind)      
      igrid        = real(  igrid_R4      , kind = dbl_kind)      
      cgrid        = real(  cgrid_R4      , kind = dbl_kind)      

      if(present(yday_R4)) then
         yday         = real(  yday_R4       , kind = dbl_kind)
      endif

      aice     = real(      aice_R4   , kind = dbl_kind)
      aice0    = real(      aice0_R4  , kind = dbl_kind)
      frain    = real(      frain_R4  , kind = dbl_kind)
      fpond    = real(      fpond_R4  , kind = dbl_kind)
      fresh    = real(      fresh_R4  , kind = dbl_kind)
      fsalt    = real(      fsalt_R4  , kind = dbl_kind)
      fhocn    = real(      fhocn_R4  , kind = dbl_kind)
      fzsal    = real(      fzsal_R4  , kind = dbl_kind)
      meltl    = real(      meltl_R4  , kind = dbl_kind)
      frazil   = real(      frazil_R4 , kind = dbl_kind)
      frazil_diag  = real(  frazil_diag_R4, kind = dbl_kind)

     salinz    = real(      salinz_R4 , kind = dbl_kind)
     ocean_bio     = real(  ocean_bio_R4  , kind = dbl_kind)

     aicen_init = real(      aicen_init_R4, kind = dbl_kind)
     vicen_init  = real(    vicen_init_R4, kind = dbl_kind)
     aicen     = real(      aicen_R4  , kind = dbl_kind)
     vicen     = real(      vicen_R4  , kind = dbl_kind)
     vsnon     = real(      vsnon_R4  , kind = dbl_kind)
     faero_ocn = real(      faero_ocn_R4, kind = dbl_kind)
     flux_bio      = real(  flux_bio_R4   , kind = dbl_kind)
     trcrn         = real(  trcrn_R4      , kind = dbl_kind)
      if (present(frz_onset_R4)) then
         frz_onset    = real(  frz_onset_R4  , kind = dbl_kind)
      endif




      l_stop = .false.

      !-----------------------------------------------------------------
      ! Let rain drain through to the ocean.
      !-----------------------------------------------------------------

      fresh  = fresh + frain * aice

      !-----------------------------------------------------------------
      ! Given thermodynamic growth rates, transport ice between
      ! thickness categories.
      !-----------------------------------------------------------------

!      call ice_timer_start(timer_catconv)    ! category conversions

      !-----------------------------------------------------------------
      ! Compute fractional ice area in each grid cell.
      !-----------------------------------------------------------------

      flux_bio(:) = c0

      call aggregate_area (ncat, aicen, aice, aice0)

      if (kitd == 1) then

      !-----------------------------------------------------------------
      ! Identify grid cells with ice.
      !-----------------------------------------------------------------

         if (aice > puny) then

            call linear_itd (ncat,     hin_max,     &
                             nilyr,    nslyr,       &
                             ntrcr,    trcr_depend, &
                             trcr_base,             & 
                             n_trcr_strata,         &
                             nt_strata,   Tf,       &
                             aicen_init,            &
                             vicen_init,            &
                             aicen,                 &
                             trcrn,                 & 
                             vicen,                 &
                             vsnon,                 &
                             aice,                  &
                             aice0,                 &
                             fpond,       l_stop,   &
                             stop_label)

            if (l_stop) return

         endif ! aice > puny

      endif  ! kitd = 1

!      call ice_timer_stop(timer_catconv)    ! category conversions

      !-----------------------------------------------------------------
      ! Add frazil ice growing in leads.
      !-----------------------------------------------------------------

      ! identify ice-ocean cells

         call add_new_ice (ncat,          nilyr,        &
                           nblyr,                       &
                           n_aero,        dt,           &
                           ntrcr,         nbtrcr,       &
                           hin_max,       ktherm,       &
                           aicen,         trcrn,        &
                           vicen,         vsnon(1),     &
                           aice0,         aice,         &
                           frzmlt,        frazil,       &
                           frz_onset,     yday,         &
                           update_ocn_f,                &
                           fresh,         fsalt,        &
                           Tf,            sss,          &
                           salinz,        phi_init,     &
                           dSin0_frazil,  bgrid,        &
                           cgrid,         igrid,        &
                           flux_bio,                    &
                           ocean_bio,     fzsal,        &
                           frazil_diag,                 &
                           l_stop,        stop_label)

         if (l_stop) return

      !-----------------------------------------------------------------
      ! Melt ice laterally.
      !-----------------------------------------------------------------

      call lateral_melt (dt,        ncat,          &
                         nilyr,     nslyr,         &
                         n_aero,    fpond,         &
                         fresh,     fsalt,         &
                         fhocn,     faero_ocn,     &
                         rside,     meltl,         &
                         aicen,     vicen,         &
                         vsnon,     trcrn,         &
                         fzsal,     flux_bio,      &
                         nbtrcr,    nblyr)

      !-----------------------------------------------------------------
      ! For the special case of a single category, adjust the area and
      ! volume (assuming that half the volume change decreases the
      ! thickness, and the other half decreases the area).  
      !-----------------------------------------------------------------

!echmod: test this
      if (ncat==1) &
          call reduce_area (hin_max   (0),                &
                            aicen     (1), vicen     (1), &
                            aicen_init(1), vicen_init(1))

      !-----------------------------------------------------------------
      ! ITD cleanup: Rebin thickness categories if necessary, and remove
      !  categories with very small areas.
      !-----------------------------------------------------------------

      call cleanup_itd (dt,                   Tf,               &
                        ntrcr,                                  &
                        nilyr,                nslyr,            &
                        ncat,                 hin_max,          &
                        aicen,                trcrn(1:ntrcr,:), &
                        vicen,                vsnon,            &
                        aice0,                aice,             & 
                        n_aero,                                 &
                        nbtrcr,               nblyr,            &
                        l_stop,               stop_label,       &
                        tr_aero,                                &
                        tr_pond_topo,         heat_capacity,    &
                        first_ice,                              &
                        trcr_depend,          trcr_base,        &
                        n_trcr_strata,        nt_strata,        &
                        fpond,                fresh,            &
                        fsalt,                fhocn,            &
                        faero_ocn,            fzsal,            &
                        flux_bio)   
      
      aice_R4     = real(      aice   , kind = real_kind)
      aice0_R4    = real(      aice0  , kind = real_kind)
      frain_R4    = real(      frain  , kind = real_kind)
      fpond_R4    = real(      fpond  , kind = real_kind)
      fresh_R4    = real(      fresh  , kind = real_kind)
      fsalt_R4    = real(      fsalt  , kind = real_kind)
      fhocn_R4    = real(      fhocn  , kind = real_kind)
      fzsal_R4    = real(      fzsal  , kind = real_kind)
      meltl_R4    = real(      meltl  , kind = real_kind)
      frazil_R4   = real(      frazil , kind = real_kind)
      frazil_diag_R4  = real(  frazil_diag, kind = real_kind)
      aicen_init_R4= real(      aicen_init, kind = real_kind)
      vicen_init_R4 = real(    vicen_init, kind = real_kind)
      aicen_R4    = real(      aicen  , kind = real_kind)
      vicen_R4    = real(      vicen  , kind = real_kind)
      vsnon_R4    = real(      vsnon  , kind = real_kind)
      faero_ocn_R4= real(      faero_ocn, kind = real_kind)
      flux_bio_R4     = real(  flux_bio   , kind = real_kind)
      trcrn_R4        = real(  trcrn      , kind = real_kind)
      if (present(frz_onset_R4)) then
         frz_onset_R4    = real(  frz_onset  , kind = real_kind)
      endif      

      deallocate(trcr_base)
      deallocate(salinz)
      deallocate(ocean_bio)
      deallocate(aicen_init)
      deallocate(vicen_init)
      deallocate(aicen)
      deallocate(vicen)
      deallocate(vsnon)
      deallocate(faero_ocn)
      deallocate(flux_bio)
      deallocate(trcrn)
      end subroutine colpkg_step_therm2

!=======================================================================
!
! Scales radiation fields computed on the previous time step.
!
! authors: Elizabeth Hunke, LANL

      subroutine colpkg_prep_radiation (ncat, nilyr, nslyr,          &
                                        aice_R4,        aicen_R4,    &
                                        swvdr_R4,       swvdf_R4,    &
                                        swidr_R4,       swidf_R4,    &
                                        alvdr_ai_R4,    alvdf_ai_R4, &
                                        alidr_ai_R4,    alidf_ai_R4, &
                                        scale_factor_R4,             &
                                        fswsfcn_R4,     fswintn_R4,  &
                                        fswthrun_R4,    fswpenln_R4, &
                                        Sswabsn_R4,     Iswabsn_R4)

      use ice_constants_colpkg, only: c0, c1, puny

      integer (kind=int_kind), intent(in) :: &
         ncat    , & ! number of ice thickness categories
         nilyr   , & ! number of ice layers
         nslyr       ! number of snow layers

      real (kind=real_kind), intent(in) :: &
         aice_R4        , & ! ice area fraction
         swvdr_R4       , & ! sw down, visible, direct  (W/m^2)
         swvdf_R4       , & ! sw down, visible, diffuse (W/m^2)
         swidr_R4       , & ! sw down, near IR, direct  (W/m^2)
         swidf_R4       , & ! sw down, near IR, diffuse (W/m^2)
         ! grid-box-mean albedos aggregated over categories (if calc_Tsfc)
         alvdr_ai_R4    , & ! visible, direct   (fraction)
         alidr_ai_R4    , & ! near-ir, direct   (fraction)
         alvdf_ai_R4    , & ! visible, diffuse  (fraction)
         alidf_ai_R4        ! near-ir, diffuse  (fraction)

      real (kind=real_kind), dimension(:), intent(in) :: &
         aicen_R4           ! ice area fraction in each category

      real (kind=real_kind), intent(inout) :: &
         scale_factor_R4    ! shortwave scaling factor, ratio new:old

      real (kind=real_kind), dimension(:), intent(inout) :: &
         fswsfcn_R4     , & ! SW absorbed at ice/snow surface (W m-2)
         fswintn_R4     , & ! SW absorbed in ice interior, below surface (W m-2)
         fswthrun_R4        ! SW through ice to ocean (W/m^2)

      real (kind=real_kind), dimension(:,:), intent(inout) :: &
         fswpenln_R4    , & ! visible SW entering ice layers (W m-2)
         Iswabsn_R4     , & ! SW radiation absorbed in ice layers (W m-2)
         Sswabsn_R4         ! SW radiation absorbed in snow layers (W m-2)

      ! local variables
      integer (kind=int_kind) :: &
         k , & ! vertical index       
         n               ! thickness category index
      real (kind=dbl_kind) :: netsw 
      
      !in variables 
      real (kind=dbl_kind) :: & 
         aice        , & ! ice area fraction
         swvdr       , & ! sw down, visible, direct  (W/m^2)
         swvdf       , & ! sw down, visible, diffuse (W/m^2)
         swidr       , & ! sw down, near IR, direct  (W/m^2)
         swidf       , & ! sw down, near IR, diffuse (W/m^2)
         ! grid-box-mean albedos aggregated over categories (if calc_Tsfc)
         alvdr_ai    , & ! visible, direct   (fraction)
         alidr_ai    , & ! near-ir, direct   (fraction)
         alvdf_ai    , & ! visible, diffuse  (fraction)
         alidf_ai        ! near-ir, diffuse  (fraction)

      real (kind=dbl_kind), dimension(:), allocatable :: & 
         aicen           ! ice area fraction in each category
      !out variables 
      real (kind=dbl_kind) :: & 
         scale_factor    ! shortwave scaling factor, ratio new:old

      real (kind=dbl_kind), dimension(:), allocatable :: & 
         fswsfcn     , & ! SW absorbed at ice/snow surface (W m-2)
         fswintn     , & ! SW absorbed in ice interior, below surface (W m-2)
         fswthrun        ! SW through ice to ocean (W/m^2)

      real (kind=dbl_kind), dimension(:,:), allocatable :: & 
         fswpenln    , & ! visible SW entering ice layers (W m-2)
         Iswabsn     , & ! SW radiation absorbed in ice layers (W m-2)
         Sswabsn         ! SW radiation absorbed in snow layers (W m-2)

      !copying R4 to R8
      aice         = real(aice_R4        , kind=dbl_kind)
      swvdr        = real(swvdr_R4       , kind=dbl_kind)
      swvdf        = real(swvdf_R4       , kind=dbl_kind)
      swidr        = real(swidr_R4       , kind=dbl_kind)
      swidf        = real(swidf_R4       , kind=dbl_kind)
      
      alvdr_ai     = real(alvdr_ai_R4    , kind=dbl_kind)
      alidr_ai     = real(alidr_ai_R4    , kind=dbl_kind)
      alvdf_ai     = real(alvdf_ai_R4    , kind=dbl_kind)
      alidf_ai     = real(alidf_ai_R4    , kind=dbl_kind)    

     aicen     = real(aicen_R4       , kind=dbl_kind)    
     fswsfcn   = real(fswsfcn_R4     , kind=dbl_kind)
     fswintn   = real(fswintn_R4     , kind=dbl_kind)
     fswthrun  = real(fswthrun_R4    , kind=dbl_kind)    
     fswpenln  = real(fswpenln_R4    , kind=dbl_kind)
     Iswabsn   = real(Iswabsn_R4     , kind=dbl_kind)
     Sswabsn   = real(Sswabsn_R4     , kind=dbl_kind)    

      scale_factor = real(scale_factor_R4, kind=dbl_kind)    

      !-----------------------------------------------------------------
      ! Compute netsw scaling factor (new netsw / old netsw)
      !-----------------------------------------------------------------

         if (aice > c0 .and. scale_factor > puny) then
            netsw = swvdr*(c1 - alvdr_ai) &
                  + swvdf*(c1 - alvdf_ai) &
                  + swidr*(c1 - alidr_ai) &
                  + swidf*(c1 - alidf_ai)
            scale_factor = netsw / scale_factor
         else
            scale_factor = c1
         endif

         do n = 1, ncat

            if (aicen(n) > puny) then

      !-----------------------------------------------------------------
      ! Scale absorbed solar radiation for change in net shortwave
      !-----------------------------------------------------------------

               fswsfcn(n)  = scale_factor*fswsfcn (n)
               fswintn(n)  = scale_factor*fswintn (n)
               fswthrun(n) = scale_factor*fswthrun(n)
               do k = 1,nilyr+1
                  fswpenln(k,n) = scale_factor*fswpenln(k,n)
               enddo       !k
               do k=1,nslyr
                  Sswabsn(k,n) = scale_factor*Sswabsn(k,n)
               enddo
               do k=1,nilyr
                  Iswabsn(k,n) = scale_factor*Iswabsn(k,n)
               enddo

            endif
         enddo                  ! ncat

         scale_factor_R4 = real(scale_factor, kind=real_kind)    
         fswsfcn_R4      = real(fswsfcn     , kind=real_kind)
         fswintn_R4      = real(fswintn     , kind=real_kind)
         fswthrun_R4     = real(fswthrun    , kind=real_kind)    
         fswpenln_R4     = real(fswpenln    , kind=real_kind)
         Iswabsn_R4      = real(Iswabsn     , kind=real_kind)
         Sswabsn_R4      = real(Sswabsn     , kind=real_kind)  

         deallocate(aicen)
         deallocate(fswsfcn)
         deallocate(fswintn)
         deallocate(fswthrun)
         deallocate(fswpenln)
         deallocate(Iswabsn)
         deallocate(Sswabsn)

      end subroutine colpkg_prep_radiation

!=======================================================================
!
! Computes radiation fields
!
! authors: William H. Lipscomb, LANL
!          David Bailey, NCAR
!          Elizabeth C. Hunke, LANL

      subroutine colpkg_step_radiation (dt_R4,       ncat,      & 
                                        n_algae,             &
                                        nblyr,    ntrcr,     &
                                        nbtrcr,   nbtrcr_sw, &
                                        nilyr,    nslyr,     &
                                        n_aero,   n_zaero,   &
                                        dEdd_algae,          &
                                        nlt_chl_sw,          &
                                        nlt_zaero_sw,        &
                                        swgrid_R4,   igrid_R4,     &
                                        fbri_R4,                &
                                        aicen_R4,    vicen_R4,     &
                                        vsnon_R4,    Tsfcn_R4,     &
                                        alvln_R4,    apndn_R4,     &
                                        hpndn_R4,    ipndn_R4,     &
                                        snwredist,           &
                                        rsnow_R4,               &
                                        aeron_R4,               &
                                        zbion_R4,               &
                                        trcrn_R4,               &
                                        TLAT_R4,     TLON_R4,      &
                                        calendar_type,       &
                                        days_per_year,       &
                                        nextsw_cday_R4,         &
                                        yday_R4,     sec,       &
                                        kaer_tab_R4, waer_tab_R4,  &
                                        gaer_tab_R4,            &
                                        kaer_bc_tab_R4,         &
                                        waer_bc_tab_R4,         &
                                        gaer_bc_tab_R4,         &
                                        bcenh_R4,               &
                                        modal_aero,          &
                                        swvdr_R4,    swvdf_R4,     &
                                        swidr_R4,    swidf_R4,     &
                                        coszen_R4,   fsnow_R4,     &
                                        alvdrn_R4,   alvdfn_R4,    &
                                        alidrn_R4,   alidfn_R4,    &
                                        fswsfcn_R4,  fswintn_R4,   &
                                        fswthrun_R4, fswpenln_R4,  &
                                        Sswabsn_R4,  Iswabsn_R4,   &
                                        albicen_R4,  albsnon_R4,   &
                                        albpndn_R4,  apeffn_R4,    &
                                        snowfracn_R4,           &
                                        dhsn_R4,     ffracn_R4,    &
                                        l_print_point,       &
                                        initonly,           &
                                        asm_prm_ice_drc_R4,     &
                                        asm_prm_ice_dfs_R4,     &
                                        ss_alb_ice_drc_R4,      &
                                        ss_alb_ice_dfs_R4,      &
                                        ext_cff_mss_ice_drc_R4, &
                                        ext_cff_mss_ice_dfs_R4, &
                                        kaer_tab_5bd_R4,          &
                                        waer_tab_5bd_R4,          &
                                        gaer_tab_5bd_R4,          &
                                        kaer_bc_tab_5bd_R4,       &
                                        waer_bc_tab_5bd_R4,       &
                                        gaer_bc_tab_5bd_R4,       &
                                        bcenh_5bd_R4,             &
                                        rsnw_dEddn_R4)

      use ice_constants_colpkg, only: c0, puny
      use ice_shortwave, only: run_dEdd, shortwave_ccsm3, compute_shortwave_trcr
      use ice_colpkg_tracers, only: tr_pond_cesm, tr_pond_lvl, tr_pond_topo, &
                                    tr_bgc_N, tr_aero, tr_rsnw, tr_zaero

      use ice_colpkg_shared, only:  z_tracers, skl_bgc

      integer (kind=int_kind), intent(in) :: &
         ncat      , & ! number of ice thickness categories
         nilyr     , & ! number of ice layers
         nslyr     , & ! number of snow layers
         n_aero    , & ! number of aerosols
         n_zaero   , & ! number of zaerosols 
         nlt_chl_sw, & ! index for chla
         nblyr     , &
         ntrcr     , &
         nbtrcr    , &
         nbtrcr_sw , &
         n_algae

      integer (kind=int_kind), dimension(:), intent(in) :: &
        nlt_zaero_sw   ! index for zaerosols

      real (kind=real_kind), intent(in) :: &
         dt_R4        , & ! time step (s)
         swvdr_R4     , & ! sw down, visible, direct  (W/m^2)
         swvdf_R4     , & ! sw down, visible, diffuse (W/m^2)
         swidr_R4     , & ! sw down, near IR, direct  (W/m^2)
         swidf_R4     , & ! sw down, near IR, diffuse (W/m^2)
         fsnow_R4     , & ! snowfall rate (kg/m^2 s)
         TLAT_R4, TLON_R4    ! latitude and longitude (radian)

      character (len=char_len), intent(in) :: &
         calendar_type       ! differentiates Gregorian from other calendars

      integer (kind=int_kind), intent(in) :: &
         days_per_year, &    ! number of days in one year
         sec                 ! elapsed seconds into date

      real (kind=real_kind), intent(in) :: &
         nextsw_cday_R4     , & ! julian day of next shortwave calculation
         yday_R4                ! day of the year

      real (kind=real_kind), intent(inout) :: &
         coszen_R4        ! cosine solar zenith angle, < 0 for sun below horizon_R4 

      real (kind=real_kind), dimension (:), intent(in) :: &
         igrid_R4              ! biology vertical interface points
 
      real (kind=real_kind), dimension (:), intent(in) :: &
         swgrid_R4                ! grid for ice tracers used in dEdd scheme
        
      real (kind=real_kind), dimension(:,:), intent(in) :: & 
         kaer_tab_R4, & ! aerosol mass extinction cross section (m2/kg)
         waer_tab_R4, & ! aerosol single scatter albedo (fraction)
         gaer_tab_R4, & ! aerosol asymmetry parameter (cos(theta))
         rsnow_R4       ! snow grain radius tracer (10^-6 m)

      real (kind=real_kind), dimension(:,:), intent(in) :: & 
         kaer_bc_tab_R4, & ! aerosol mass extinction cross section (m2/kg)
         waer_bc_tab_R4, & ! aerosol single scatter albedo (fraction)
         gaer_bc_tab_R4    ! aerosol asymmetry parameter (cos(theta))

      real (kind=real_kind), dimension(:,:,:), intent(in) :: & 
         bcenh_R4 

      real (kind=real_kind), dimension(:), intent(in) :: &
         aicen_R4     , & ! ice area fraction in each category
         vicen_R4     , & ! ice volume in each category (m)
         vsnon_R4     , & ! snow volume in each category (m)
         Tsfcn_R4     , & ! surface temperature (deg C)
         alvln_R4     , & ! level-ice area fraction
         apndn_R4     , & ! pond area fraction
         hpndn_R4     , & ! pond depth (m)
         ipndn_R4     , & ! pond refrozen lid thickness (m)
         fbri_R4           ! brine fraction 

      character(len=char_len), intent(in) :: & 
         snwredist                ! type of snow redistribution

      real(kind=real_kind), dimension(:,:), intent(in) :: &
         aeron_R4     , & ! aerosols (kg/m^3)
         trcrn_R4         ! tracers

      real(kind=real_kind), dimension(:,:), intent(inout) :: &
         zbion_R4        ! zaerosols (kg/m^3) and chla (mg/m^3)

      real (kind=real_kind), dimension(:), intent(inout) :: &
         alvdrn_R4    , & ! visible, direct  albedo (fraction)
         alidrn_R4    , & ! near-ir, direct   (fraction)
         alvdfn_R4    , & ! visible, diffuse  (fraction)
         alidfn_R4    , & ! near-ir, diffuse  (fraction)
         fswsfcn_R4   , & ! SW absorbed at ice/snow surface (W m-2)
         fswintn_R4   , & ! SW absorbed in ice interior, below surface (W m-2)
         fswthrun_R4  , & ! SW through ice to ocean (W/m^2)
         snowfracn_R4, & ! snow fraction on each category
         dhsn_R4      , & ! depth difference for snow on sea ice and pond ice_R4
         ffracn_R4    , & ! fraction of fsurfn used to melt ipond
                       ! albedo components for history
         albicen_R4   , & ! bare ice 
         albsnon_R4   , & ! snow 
         albpndn_R4   , & ! pond 
         rsnw_dEddn_R4, & ! snow grain radius (um)
         apeffn_R4        ! effective pond area used for radiation calculation_R4

      real (kind=real_kind), dimension(:,:), intent(inout) :: &
         fswpenln_R4  , & ! visible SW entering ice layers (W m-2)
         Iswabsn_R4   , & ! SW radiation absorbed in ice layers (W m-2)
         Sswabsn_R4       ! SW radiation absorbed in snow layers (W m-2)

      logical (kind=log_kind), intent(in) :: &
         l_print_point, & ! flag for printing diagnostics
         dEdd_algae   , & ! .true. use prognostic chla in dEdd
         modal_aero       ! .true. use modal aerosol optical treatment

      logical (kind=log_kind), optional :: &
         initonly         ! flag to indicate init only, default is false


      ! snow grain single-scattering properties for
      ! direct (drc) and diffuse (dfs) shortwave incidents
      real (kind=real_kind), dimension(:,:), intent(in) :: & ! Model SNICAR snow SSP
        asm_prm_ice_drc_R4,     & ! snow asymmetry factor (cos(theta))
        asm_prm_ice_dfs_R4,     & ! snow asymmetry factor (cos(theta))
        ss_alb_ice_drc_R4,      & ! snow single scatter albedo (fraction)
        ss_alb_ice_dfs_R4,      & ! snow single scatter albedo (fraction)
        ext_cff_mss_ice_drc_R4, & ! snow mass extinction cross section (m2/kg)
        ext_cff_mss_ice_dfs_R4    ! snow mass extinction cross section (m2/kg)

      real (kind=real_kind), dimension(:,:), intent(in) :: &
        kaer_tab_5bd_R4, & ! aerosol mass extinction cross section (m2/kg)
        waer_tab_5bd_R4, & ! aerosol single scatter albedo (fraction)
        gaer_tab_5bd_R4    ! aerosol asymmetry parameter (cos(theta))

      real (kind=real_kind), dimension(:,:), intent(in) :: & ! Modal aerosol treatment
        kaer_bc_tab_5bd_R4, & ! aerosol mass extinction cross section (m2/kg)
        waer_bc_tab_5bd_R4, & ! aerosol single scatter albedo (fraction)
        gaer_bc_tab_5bd_R4    ! aerosol asymmetry parameter (cos(theta))

      real (kind=real_kind), dimension(:,:,:), intent(in) :: & ! Modal aerosol treatment
        bcenh_5bd_R4         ! BC absorption enhancement factor

      ! local variables

      integer (kind=int_kind) :: &
         n                  ! thickness category index

      logical (kind=log_kind) :: &
         l_stop      ,&  ! if true, abort the model
         linitonly       ! local flag for initonly

      character (char_len) :: stop_label

      real(kind=dbl_kind) :: &
        hin,         & ! Ice thickness (m)
        hbri           ! brine thickness (m)


      real (kind=dbl_kind) :: &
         dt        , & ! time step (s)
         swvdr     , & ! sw down, visible, direct  (W/m^2)
         swvdf     , & ! sw down, visible, diffuse (W/m^2)
         swidr     , & ! sw down, near IR, direct  (W/m^2)
         swidf     , & ! sw down, near IR, diffuse (W/m^2)
         fsnow     , & ! snowfall rate (kg/m^2 s)
         TLAT, TLON    ! latitude and longitude (radian)

      !!!!!!!!
      ! local variables copies 
      !!!!!!!!

      real (kind=dbl_kind) :: &
         nextsw_cday     , & ! julian day of next shortwave calculation
         yday                ! day of the year

      real (kind=dbl_kind):: &
         coszen        ! cosine solar zenith angle, < 0 for sun below horizon 

      real (kind=dbl_kind), dimension (:), allocatable :: &
         igrid              ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (:), allocatable :: &
         swgrid                ! grid for ice tracers used in dEdd scheme
        
      real (kind=dbl_kind), dimension(:,:), allocatable :: & 
         kaer_tab, & ! aerosol mass extinction cross section (m2/kg)
         waer_tab, & ! aerosol single scatter albedo (fraction)
         gaer_tab, & ! aerosol asymmetry parameter (cos(theta))
         rsnow       ! snow grain radius tracer (10^-6 m)

      real (kind=dbl_kind), dimension(:,:), allocatable :: & 
         kaer_bc_tab, & ! aerosol mass extinction cross section (m2/kg)
         waer_bc_tab, & ! aerosol single scatter albedo (fraction)
         gaer_bc_tab    ! aerosol asymmetry parameter (cos(theta))

      real (kind=dbl_kind), dimension(:,:,:), allocatable :: & 
         bcenh 

      real (kind=dbl_kind), dimension(:), allocatable :: &
         aicen     , & ! ice area fraction in each category
         vicen     , & ! ice volume in each category (m)
         vsnon     , & ! snow volume in each category (m)
         Tsfcn     , & ! surface temperature (deg C)
         alvln     , & ! level-ice area fraction
         apndn     , & ! pond area fraction
         hpndn     , & ! pond depth (m)
         ipndn     , & ! pond refrozen lid thickness (m)
         fbri           ! brine fraction 

      real(kind=dbl_kind), dimension(:,:), allocatable :: &
         aeron     , & ! aerosols (kg/m^3)
         trcrn         ! tracers

      real(kind=dbl_kind), dimension(:,:), allocatable :: &
         zbion        ! zaerosols (kg/m^3) and chla (mg/m^3)

      real (kind=dbl_kind), dimension(:), allocatable:: &
         alvdrn    , & ! visible, direct  albedo (fraction)
         alidrn    , & ! near-ir, direct   (fraction)
         alvdfn    , & ! visible, diffuse  (fraction)
         alidfn    , & ! near-ir, diffuse  (fraction)
         fswsfcn   , & ! SW absorbed at ice/snow surface (W m-2)
         fswintn   , & ! SW absorbed in ice interior, below surface (W m-2)
         fswthrun  , & ! SW through ice to ocean (W/m^2)
         snowfracn, & ! snow fraction on each category
         dhsn      , & ! depth difference for snow on sea ice and pond ice
         ffracn    , & ! fraction of fsurfn used to melt ipond
                       ! albedo components for history
         albicen   , & ! bare ice 
         albsnon   , & ! snow 
         albpndn   , & ! pond 
         rsnw_dEddn, & ! snow grain radius (um)
         apeffn        ! effective pond area used for radiation calculation

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         fswpenln  , & ! visible SW entering ice layers (W m-2)
         Iswabsn   , & ! SW radiation absorbed in ice layers (W m-2)
         Sswabsn       ! SW radiation absorbed in snow layers (W m-2)

      ! snow grain single-scattering properties for
      ! direct (drc) and diffuse (dfs) shortwave incidents
      real (kind=dbl_kind), dimension(:,:), allocatable :: & ! Model SNICAR snow SSP
        asm_prm_ice_drc,     & ! snow asymmetry factor (cos(theta))
        asm_prm_ice_dfs,     & ! snow asymmetry factor (cos(theta))
        ss_alb_ice_drc,      & ! snow single scatter albedo (fraction)
        ss_alb_ice_dfs,      & ! snow single scatter albedo (fraction)
        ext_cff_mss_ice_drc, & ! snow mass extinction cross section (m2/kg)
        ext_cff_mss_ice_dfs    ! snow mass extinction cross section (m2/kg)

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
        kaer_tab_5bd, & ! aerosol mass extinction cross section (m2/kg)
        waer_tab_5bd, & ! aerosol single scatter albedo (fraction)
        gaer_tab_5bd    ! aerosol asymmetry parameter (cos(theta))

      real (kind=dbl_kind), dimension(:,:), allocatable :: & ! Modal aerosol treatment
        kaer_bc_tab_5bd, & ! aerosol mass extinction cross section (m2/kg)
        waer_bc_tab_5bd, & ! aerosol single scatter albedo (fraction)
        gaer_bc_tab_5bd    ! aerosol asymmetry parameter (cos(theta))

      real (kind=dbl_kind), dimension(:,:,:), allocatable :: & ! Modal aerosol treatment
        bcenh_5bd         ! BC absorption enhancement factor

      swgrid =real(swgrid_R4, kind=dbl_kind)
      igrid =real(igrid_R4  , kind=dbl_kind)

      dt    = real(dt_R4, kind = dbl_kind)     
      swvdr = real(swvdr_R4, kind = dbl_kind)     
      swvdf = real(swvdf_R4, kind = dbl_kind)     
      swidr = real(swidr_R4, kind = dbl_kind)     
      swidf = real(swidf_R4, kind = dbl_kind)     
      fsnow = real(fsnow_R4, kind = dbl_kind)     
      TLAT  = real(TLAT_R4, kind = dbl_kind)
      TLON  = real(TLON_R4, kind = dbl_kind)
      nextsw_cday = real(nextsw_cday_R4, kind=dbl_kind)     
      yday        = real(yday_R4       , kind=dbl_kind)         
   
      coszen      = real(coszen_R4     , kind=dbl_kind)   
   
     kaer_tab     = real(kaer_tab_R4   , kind=dbl_kind)
     waer_tab     = real(waer_tab_R4   , kind=dbl_kind)
     gaer_tab     = real(gaer_tab_R4   , kind=dbl_kind) 
     rsnow        = real(rsnow_R4      , kind=dbl_kind) 

     kaer_bc_tab  = real(kaer_bc_tab_R4, kind = dbl_kind)
     waer_bc_tab  = real(waer_bc_tab_R4, kind = dbl_kind)
     gaer_bc_tab  = real(gaer_bc_tab_R4, kind = dbl_kind)    
   
     bcenh  = real(bcenh_R4, kind = dbl_kind) 
   
     aicen  = real(aicen_R4, kind = dbl_kind)     
     vicen  = real(vicen_R4, kind = dbl_kind)     
     vsnon  = real(vsnon_R4, kind = dbl_kind)     
     Tsfcn  = real(Tsfcn_R4, kind = dbl_kind)     
     alvln  = real(alvln_R4, kind = dbl_kind)     
     apndn  = real(apndn_R4, kind = dbl_kind)     
     hpndn  = real(hpndn_R4, kind = dbl_kind)     
     ipndn  = real(ipndn_R4, kind = dbl_kind)     
     fbri  = real(fbri_R4, kind = dbl_kind)          

     aeron     = real(aeron_R4    , kind = dbl_kind)       
     trcrn     = real(trcrn_R4    , kind = dbl_kind)         

     zbion     = real(zbion_R4    , kind = dbl_kind)       

     alvdrn    = real(alvdrn_R4   , kind = dbl_kind)      
     alidrn    = real(alidrn_R4   , kind = dbl_kind)       
     alvdfn    = real(alvdfn_R4   , kind = dbl_kind)    
     alidfn    = real(alidfn_R4   , kind = dbl_kind)      
     fswsfcn   = real(fswsfcn_R4  , kind = dbl_kind)  
     fswintn   = real(fswintn_R4  , kind = dbl_kind)  
     fswthrun  = real(fswthrun_R4 , kind = dbl_kind)  
     snowfracn  = real(snowfracn_R4, kind = dbl_kind) 
     dhsn  = real(dhsn_R4     , kind = dbl_kind)    
     ffracn    = real(ffracn_R4   , kind = dbl_kind)  
     albicen     = real(albicen_R4   , kind = dbl_kind)
     albsnon     = real(albsnon_R4   , kind = dbl_kind)
     albpndn     = real(albpndn_R4   , kind = dbl_kind)
     rsnw_dEddn  = real(rsnw_dEddn_R4, kind = dbl_kind)
     apeffn      = real(apeffn_R4    , kind = dbl_kind)    
   
     fswpenln    = real(fswpenln_R4, kind = dbl_kind)  
     Iswabsn     = real(Iswabsn_R4 , kind = dbl_kind)  
     Sswabsn     = real(Sswabsn_R4 , kind = dbl_kind)      

     asm_prm_ice_drc  = real(asm_prm_ice_drc_R4, kind = dbl_kind)
     asm_prm_ice_dfs  = real(asm_prm_ice_dfs_R4, kind = dbl_kind)
     ss_alb_ice_drc  = real(ss_alb_ice_drc_R4, kind = dbl_kind)
     ss_alb_ice_dfs  = real(ss_alb_ice_dfs_R4, kind = dbl_kind)
     ext_cff_mss_ice_drc  = real(ext_cff_mss_ice_drc_R4, kind = dbl_kind)
     ext_cff_mss_ice_dfs  = real(ext_cff_mss_ice_dfs_R4, kind = dbl_kind)    

     kaer_tab_5bd  = real(kaer_tab_5bd_R4, kind = dbl_kind)
     waer_tab_5bd  = real(waer_tab_5bd_R4, kind = dbl_kind)
     gaer_tab_5bd  = real(gaer_tab_5bd_R4, kind = dbl_kind)    

     kaer_bc_tab_5bd  = real(kaer_bc_tab_5bd_R4, kind = dbl_kind)
     waer_bc_tab_5bd  = real(waer_bc_tab_5bd_R4, kind = dbl_kind)
     gaer_bc_tab_5bd  = real(gaer_bc_tab_5bd_R4, kind = dbl_kind)    

     bcenh_5bd  = real(bcenh_5bd_R4, kind = dbl_kind)         


      !!!!!!!!
      ! local variables copies end
      !!!!!!!!
        hin = c0
        hbri = c0
        linitonly = .false.
        if (present(initonly)) then
           linitonly = initonly
        endif

         ! Initialize
         do n = 1, ncat
            alvdrn  (n) = c0
            alidrn  (n) = c0
            alvdfn  (n) = c0
            alidfn  (n) = c0
            fswsfcn (n) = c0
            fswintn (n) = c0
            fswthrun(n) = c0
         enddo   ! ncat
         fswpenln (:,:) = c0
         Iswabsn  (:,:) = c0
         Sswabsn  (:,:) = c0
         zbion(:,:) = c0

         ! Interpolate z-shortwave tracers to shortwave grid
         if (dEdd_algae) then
         do n = 1, ncat
              if (aicen(n) .gt. puny) then
                 hin = vicen(n)/aicen(n)
                 hbri= fbri(n)*hin
                 call compute_shortwave_trcr(n_algae, nslyr,  &
                                     trcrn(1:ntrcr,n),        &
                                     zbion(1:nbtrcr_sw,n),    &
                                     swgrid,       hin,       &
                                     hbri,         ntrcr,     &
                                     nilyr,        nblyr,     &
                                     igrid,                   &
                                     nbtrcr_sw,    n_zaero,   &
                                     skl_bgc,      z_tracers, &
                                     l_stop,       stop_label)
              endif
         enddo
         endif

         if (calc_Tsfc) then
         if (trim(shortwave) == 'dEdd') then ! delta Eddington
            
            call run_dEdd(dt,           tr_aero,        &
                          tr_pond_cesm,                 &
                          tr_pond_lvl,                  &
                          tr_pond_topo,                 &
                          ncat,         n_aero,         &
                          n_zaero,      dEdd_algae,     &
                          nlt_chl_sw,   nlt_zaero_sw,   &
                          tr_bgc_N,     tr_zaero,       &
                          nilyr,        nslyr,          &
                          aicen,        vicen,          &
                          vsnon,        Tsfcn,          &
                          alvln,        apndn,          &
                          hpndn,        ipndn,          &
                          snwredist,                    &
                          rsnow,        tr_rsnw,        &
                          aeron,        kalg,           &
                          zbion,                        &
                          heat_capacity,                &
                          TLAT,         TLON,           &
                          calendar_type,days_per_year,  &
                          nextsw_cday,  yday,           &
                          sec,          R_ice,          &
                          R_pnd,        R_snw,          &
                          dT_mlt,       rsnw_mlt,       &
                          hs0,          hs1,            &
                          hp1,          pndaspect,      &
                          kaer_tab,     waer_tab,       &
                          gaer_tab,                     &
                          kaer_bc_tab,                  &
                          waer_bc_tab,                  &
                          gaer_bc_tab,                  &
                          bcenh,                        &
                          modal_aero,                   &
                          swvdr,        swvdf,          &
                          swidr,        swidf,          &
                          coszen,       fsnow,          &
                          alvdrn,       alvdfn,         &
                          alidrn,       alidfn,         &
                          fswsfcn,      fswintn,        &
                          fswthrun,     fswpenln,       &
                          Sswabsn,      Iswabsn,        &
                          albicen,      albsnon,        &
                          albpndn,      apeffn,         &
                          snowfracn,                    &
                          dhsn,         ffracn,         &
                          rsnw_dEddn,                   &
                          l_print_point,                &
                          linitonly,                    &
                          use_snicar,                   &
                          asm_prm_ice_drc,              &
                          asm_prm_ice_dfs,              &
                          ss_alb_ice_drc,               &
                          ss_alb_ice_dfs,               &
                          ext_cff_mss_ice_drc,          &
                          ext_cff_mss_ice_dfs,          &
                          kaer_tab_5bd,                 &
                          waer_tab_5bd,                 &
                          gaer_tab_5bd,                 &
                          kaer_bc_tab_5bd,              &
                          waer_bc_tab_5bd,              &
                          gaer_bc_tab_5bd,              &
                          bcenh_5bd)
 
         else  ! .not. dEdd

            call shortwave_ccsm3(aicen,      vicen,      &
                                 vsnon,                  &
                                 Tsfcn,                  &
                                 swvdr,      swvdf,      &
                                 swidr,      swidf,      &
                                 heat_capacity,          &
                                 albedo_type,            &
                                 albicev,    albicei,    &
                                 albsnowv,   albsnowi,   &
                                 ahmax,                  &
                                 alvdrn,     alidrn,     &
                                 alvdfn,     alidfn,     &
                                 fswsfcn,    fswintn,    &
                                 fswthrun,               &
                                 fswpenln,               &
                                 Iswabsn,                &
                                 Sswabsn,                &
                                 albicen,    albsnon,    &
                                 coszen,     ncat)

         endif   ! shortwave

      else    ! .not. calc_Tsfc

      ! Calculate effective pond area for HadGEM

         if (tr_pond_topo) then
            do n = 1, ncat
               apeffn(n) = c0 
               if (aicen(n) > puny) then
               ! Lid effective if thicker than hp1
                 if (apndn(n)*aicen(n) > puny .and. ipndn(n) < hp1) then
                    apeffn(n) = apndn(n)
                 else
                    apeffn(n) = c0
                 endif
                 if (apndn(n) < puny) apeffn(n) = c0
               endif
            enddo  ! ncat
 
         endif ! tr_pond_topo

         ! Initialize for safety
         do n = 1, ncat
            alvdrn(n) = c0
            alidrn(n) = c0
            alvdfn(n) = c0
            alidfn(n) = c0
            fswsfcn(n) = c0
            fswintn(n) = c0
            fswthrun(n) = c0
         enddo   ! ncat
         Iswabsn(:,:) = c0
         Sswabsn(:,:) = c0

      endif    ! calc_Tsfc

      zbion_R4    = real(zbion    , kind = real_kind)       

      alvdrn_R4   = real(alvdrn   , kind = real_kind)      
      alidrn_R4   = real(alidrn   , kind = real_kind)       
      alvdfn_R4   = real(alvdfn   , kind = real_kind)    
      alidfn_R4   = real(alidfn   , kind = real_kind)      
      fswsfcn_R4  = real(fswsfcn  , kind = real_kind)  
      fswintn_R4  = real(fswintn  , kind = real_kind)  
      fswthrun_R4 = real(fswthrun , kind = real_kind)  
      snowfracn_R4 = real(snowfracn, kind = real_kind) 
      dhsn_R4     = real(dhsn     , kind = real_kind)    
      ffracn_R4   = real(ffracn   , kind = real_kind)  

      albicen_R4    = real(albicen   , kind = real_kind)
      albsnon_R4    = real(albsnon   , kind = real_kind)
      albpndn_R4    = real(albpndn   , kind = real_kind)
      rsnw_dEddn_R4 = real(rsnw_dEddn, kind = real_kind)
      apeffn_R4     = real(apeffn    , kind = real_kind)    
   
      fswpenln_R4   = real(fswpenln, kind = real_kind)  
      Iswabsn_R4    = real(Iswabsn , kind = real_kind)  
      Sswabsn_R4    = real(Sswabsn , kind = real_kind)      



      ! asm_prm_ice_drc_R4 = real(asm_prm_ice_drc, kind = real_kind)
      ! asm_prm_ice_dfs_R4 = real(asm_prm_ice_dfs, kind = real_kind)
      ! ss_alb_ice_drc_R4 = real(ss_alb_ice_drc, kind = real_kind)
      ! ss_alb_ice_dfs_R4 = real(ss_alb_ice_dfs, kind = real_kind)
      ! ext_cff_mss_ice_drc_R4 = real(ext_cff_mss_ice_drc, kind = real_kind)
      ! ext_cff_mss_ice_dfs_R4 = real(ext_cff_mss_ice_dfs, kind = real_kind)    

   
      ! kaer_tab_5bd_R4 = real(kaer_tab_5bd, kind = real_kind)
      ! waer_tab_5bd_R4 = real(waer_tab_5bd, kind = real_kind)
      ! gaer_tab_5bd_R4 = real(gaer_tab_5bd, kind = real_kind)    

   
      ! kaer_bc_tab_5bd_R4 = real(kaer_bc_tab_5bd, kind = real_kind)
      ! waer_bc_tab_5bd_R4 = real(waer_bc_tab_5bd, kind = real_kind)
      ! gaer_bc_tab_5bd_R4 = real(gaer_bc_tab_5bd, kind = real_kind)    

   
      ! bcenh_5bd_R4 = real(bcenh_5bd, kind = real_kind)       

      deallocate (swgrid)
      deallocate (igrid)
      deallocate(kaer_tab)
      deallocate(waer_tab)
      deallocate(gaer_tab)
      deallocate(rsnow)
      deallocate(kaer_bc_tab)
      deallocate(waer_bc_tab)
      deallocate(gaer_bc_tab)
      deallocate(bcenh)
      deallocate(aicen)
      deallocate(vicen)
      deallocate(vsnon)
      deallocate(Tsfcn)
      deallocate(alvln)
      deallocate(apndn)
      deallocate(hpndn)
      deallocate(ipndn)
      deallocate(fbri)
      deallocate(aeron)
      deallocate(trcrn)
      deallocate(zbion)
      deallocate(alvdrn)
      deallocate(alidrn)
      deallocate(alvdfn)
      deallocate(alidfn)
      deallocate(fswsfcn)
      deallocate(fswintn)
      deallocate(fswthrun)
      deallocate(snowfracn)
      deallocate(dhsn)
      deallocate(ffracn)
      deallocate(albicen)
      deallocate(albsnon)
      deallocate(albpndn)
      deallocate(rsnw_dEddn)
      deallocate(apeffn)
      deallocate(fswpenln)
      deallocate(Iswabsn)
      deallocate(Sswabsn)
      deallocate(asm_prm_ice_drc)
      deallocate(asm_prm_ice_dfs)
      deallocate(ss_alb_ice_drc)
      deallocate(ss_alb_ice_dfs)
      deallocate(ext_cff_mss_ice_drc)
      deallocate(ext_cff_mss_ice_dfs)
      deallocate(kaer_tab_5bd)
      deallocate(waer_tab_5bd)
      deallocate(gaer_tab_5bd)
      deallocate(kaer_bc_tab_5bd)
      deallocate(waer_bc_tab_5bd)
      deallocate(gaer_bc_tab_5bd)
      deallocate(bcenh_5bd)
      end subroutine colpkg_step_radiation

!=======================================================================
!
! Computes sea ice mechanical deformation
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine colpkg_step_ridge (dt_R4,           ndtd,          &
                                    nilyr,        nslyr,         &
                                    nblyr,                       &
                                    ncat,         hin_max_R4,       &
                                    rdg_conv_R4,     rdg_shear_R4,     &
                                    Tf_R4,                          &
                                    aicen_R4,                       &
                                    trcrn_R4,                       &
                                    vicen_R4,        vsnon_R4,         &
                                    aice0_R4,        trcr_depend,   &
                                    trcr_base_R4,    n_trcr_strata, &
                                    nt_strata,                   &
                                    dardg1dt_R4,     dardg2dt_R4,      &
                                    dvirdgdt_R4,     opening_R4,       &
                                    fpond_R4,                       &
                                    fresh_R4,        fhocn_R4,         &
                                    n_aero,                      &
                                    faero_ocn_R4,                   &
                                    aparticn_R4,     krdgn_R4,         &
                                    aredistn_R4,     vredistn_R4,      &
                                    dardg1ndt_R4,    dardg2ndt_R4,     &
                                    dvirdgndt_R4,                   &
                                    araftn_R4,       vraftn_R4,        &
                                    aice_R4,         fsalt_R4,         &
                                    first_ice,    fzsal_R4,         &
                                    flux_bio_R4,                    &
                                    l_stop,       stop_label)

      use ice_mechred, only: ridge_ice
      use ice_itd, only: cleanup_itd
      use ice_colpkg_tracers, only: tr_pond_topo, tr_aero, tr_brine, ntrcr, nbtrcr

      real (kind=real_kind), intent(in) :: &
         dt_R4    , & ! time step
         Tf_R4        ! ocean freezing temperature

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         ndtd  , & ! number of dynamics supercycles
         nblyr , & ! number of bio layers
         nilyr , & ! number of ice layers
         nslyr , & ! number of snow layers
         n_aero    ! number of aerosol tracers

      real (kind=real_kind), dimension(0:ncat), intent(inout) :: &
         hin_max_R4   ! category limits (m)

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=real_kind), dimension (:,:), intent(in) :: &
         trcr_base_R4      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      real (kind=real_kind), intent(inout) :: &
         aice_R4     , & ! sea ice concentration
         aice0_R4    , & ! concentration of open water
         rdg_conv_R4 , & ! convergence term for ridging (1/s)
         rdg_shear_R4, & ! shear term for ridging (1/s)
         dardg1dt_R4 , & ! rate of area loss by ridging ice (1/s)
         dardg2dt_R4 , & ! rate of area gain by new ridges (1/s)
         dvirdgdt_R4 , & ! rate of ice volume ridged (m/s)
         opening_R4  , & ! rate of opening due to divergence/shear (1/s)
         fpond_R4    , & ! fresh water flux to ponds (kg/m^2/s)
         fresh_R4    , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt_R4    , & ! salt flux to ocean (kg/m^2/s)
         fhocn_R4    , & ! net heat flux to ocean (W/m^2)
         fzsal_R4        ! zsalinity flux to ocean(kg/m^2/s)

      real (kind=real_kind), dimension(:), intent(inout) :: &
         aicen_R4    , & ! concentration of ice
         vicen_R4    , & ! volume per unit area of ice          (m)
         vsnon_R4    , & ! volume per unit area of snow         (m)
         dardg1ndt_R4, & ! rate of area loss by ridging ice (1/s)
         dardg2ndt_R4, & ! rate of area gain by new ridges (1/s)
         dvirdgndt_R4, & ! rate of ice volume ridged (m/s)
         aparticn_R4 , & ! participation function
         krdgn_R4    , & ! mean ridge thickness/thickness of ridging ice
         araftn_R4   , & ! rafting ice area
         vraftn_R4   , & ! rafting ice volume 
         aredistn_R4 , & ! redistribution function: fraction of new ridge area_R4
         vredistn_R4 , & ! redistribution function: fraction of new ridge volume_R4
         faero_ocn_R4, & ! aerosol flux to ocean  (kg/m^2/s)
         flux_bio_R4     ! all bio fluxes to ocean

      real (kind=real_kind), dimension(:,:), intent(inout) :: &
         trcrn_R4        ! tracers

      !logical (kind=log_kind), intent(in) :: &
         !tr_pond_topo,& ! if .true., use explicit topography-based ponds
         !tr_aero     ,& ! if .true., use aerosol tracers
         !tr_brine    !,& ! if .true., brine height differs from ice thickness
         !heat_capacity  ! if true, ice has nonzero heat capacity

      logical (kind=log_kind), dimension(:), intent(inout) :: &
         first_ice    ! true until ice forms

      logical (kind=log_kind), intent(out) :: &
         l_stop       ! if true, abort the model

      character (len=*), intent(out) :: &
         stop_label   ! diagnostic information for abort

      ! local variables

      real (kind=dbl_kind) :: &
         dtt      , & ! thermo time step
         atmp     , & ! temporary ice area
         atmp0        ! temporary open water area

      real (kind=dbl_kind) :: &
         dt    , & ! time step
         Tf        ! ocean freezing temperature

      real (kind=dbl_kind), dimension(0:ncat) :: &
         hin_max   ! category limits (m)

      real (kind=dbl_kind), dimension (:,:), allocatable :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      real (kind=dbl_kind) :: &
         aice     , & ! sea ice concentration
         aice0    , & ! concentration of open water
         rdg_conv , & ! convergence term for ridging (1/s)
         rdg_shear, & ! shear term for ridging (1/s)
         dardg1dt , & ! rate of area loss by ridging ice (1/s)
         dardg2dt , & ! rate of area gain by new ridges (1/s)
         dvirdgdt , & ! rate of ice volume ridged (m/s)
         opening  , & ! rate of opening due to divergence/shear (1/s)
         fpond    , & ! fresh water flux to ponds (kg/m^2/s)
         fresh    , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt    , & ! salt flux to ocean (kg/m^2/s)
         fhocn    , & ! net heat flux to ocean (W/m^2)
         fzsal        ! zsalinity flux to ocean(kg/m^2/s)

      real (kind=dbl_kind), dimension(:), allocatable:: &
         aicen    , & ! concentration of ice
         vicen    , & ! volume per unit area of ice          (m)
         vsnon    , & ! volume per unit area of snow         (m)
         dardg1ndt, & ! rate of area loss by ridging ice (1/s)
         dardg2ndt, & ! rate of area gain by new ridges (1/s)
         dvirdgndt, & ! rate of ice volume ridged (m/s)
         aparticn , & ! participation function
         krdgn    , & ! mean ridge thickness/thickness of ridging ice
         araftn   , & ! rafting ice area
         vraftn   , & ! rafting ice volume 
         aredistn , & ! redistribution function: fraction of new ridge area
         vredistn , & ! redistribution function: fraction of new ridge volume
         faero_ocn, & ! aerosol flux to ocean  (kg/m^2/s)
         flux_bio     ! all bio fluxes to ocean

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         trcrn        ! tracers

      dt = real(dt_R4, kind = dbl_kind)    
      Tf = real(Tf_R4, kind = dbl_kind)        

      hin_max = real(hin_max_R4, kind = dbl_kind)   

     trcr_base  = real(trcr_base_R4, kind = dbl_kind)     

      aice = real(aice_R4, kind = dbl_kind)     
      aice0 = real(aice0_R4, kind = dbl_kind)    
      rdg_conv = real(rdg_conv_R4, kind = dbl_kind) 
      rdg_shear = real(rdg_shear_R4, kind = dbl_kind)
      dardg1dt = real(dardg1dt_R4, kind = dbl_kind) 
      dardg2dt = real(dardg2dt_R4, kind = dbl_kind) 
      dvirdgdt = real(dvirdgdt_R4, kind = dbl_kind) 
      opening = real(opening_R4, kind = dbl_kind)  
      fpond = real(fpond_R4, kind = dbl_kind)    
      fresh = real(fresh_R4, kind = dbl_kind)    
      fsalt = real(fsalt_R4, kind = dbl_kind)    
      fhocn = real(fhocn_R4, kind = dbl_kind)    
      fzsal = real(fzsal_R4, kind = dbl_kind)        

     aicen  = real(aicen_R4, kind = dbl_kind)    
     vicen  = real(vicen_R4, kind = dbl_kind)    
     vsnon  = real(vsnon_R4, kind = dbl_kind)    
     dardg1ndt  = real(dardg1ndt_R4, kind = dbl_kind)
     dardg2ndt  = real(dardg2ndt_R4, kind = dbl_kind)
     dvirdgndt  = real(dvirdgndt_R4, kind = dbl_kind)
     aparticn  = real(aparticn_R4, kind = dbl_kind) 
     krdgn  = real(krdgn_R4, kind = dbl_kind)    
     araftn  = real(araftn_R4, kind = dbl_kind)   
     vraftn  = real(vraftn_R4, kind = dbl_kind)   
     aredistn  = real(aredistn_R4, kind = dbl_kind) 
     vredistn  = real(vredistn_R4, kind = dbl_kind) 
     faero_ocn  = real(faero_ocn_R4, kind = dbl_kind)
     flux_bio  = real(flux_bio_R4, kind = dbl_kind)     

     trcrn  = real(trcrn_R4, kind = dbl_kind)        


      l_stop = .false.

      !-----------------------------------------------------------------
      ! Identify ice-ocean cells.
      ! Note:  We can not limit the loop here using aice>puny because
      !        aice has not yet been updated since the transport (and
      !        it may be out of whack, which the ridging helps fix).-ECH
      !-----------------------------------------------------------------
           
         call ridge_ice (dt,           ndtd,           &
                         ncat,         n_aero,         &
                         nilyr,        nslyr,          &
                         ntrcr,        hin_max,        &
                         rdg_conv,     rdg_shear,      &
                         aicen,                        &
                         trcrn,                        &
                         vicen,        vsnon,          &
                         aice0,                        &
                         trcr_depend,                  &
                         trcr_base,                    &
                         n_trcr_strata,                &
                         nt_strata,                    &
                         l_stop,                       &
                         stop_label,                   &
                         krdg_partic, krdg_redist, &
                         mu_rdg,                   &
                         dardg1dt,     dardg2dt,       &
                         dvirdgdt,     opening,        &
                         fpond,                        &
                         fresh,        fhocn,          &
                         tr_brine,     faero_ocn,      &
                         aparticn,     krdgn,          &
                         aredistn,     vredistn,       &
                         dardg1ndt,    dardg2ndt,      &
                         dvirdgndt,                    &
                         araftn,       vraftn,         &
                         Tf)        

         if (l_stop) return

      !-----------------------------------------------------------------
      ! ITD cleanup: Rebin thickness categories if necessary, and remove
      !  categories with very small areas.
      !-----------------------------------------------------------------

      dtt = dt * ndtd  ! for proper averaging over thermo timestep
      call cleanup_itd (dtt,                  Tf,               &
                        ntrcr,                                  &
                        nilyr,                nslyr,            &
                        ncat,                 hin_max,          &
                        aicen,                trcrn,            &
                        vicen,                vsnon,            &
                        aice0,                aice,             &          
                        n_aero,                                 &
                        nbtrcr,               nblyr,            &
                        l_stop,               stop_label,       &
                        tr_aero,                                &
                        tr_pond_topo,         heat_capacity,    &  
                        first_ice,                              &                
                        trcr_depend,          trcr_base,        &
                        n_trcr_strata,        nt_strata,        &
                        fpond,                fresh,            &
                        fsalt,                fhocn,            &
                        faero_ocn,            fzsal,            &
                        flux_bio)

      if (l_stop) then
         stop_label = 'ice: ITD cleanup error in colpkg_step_ridge'
      endif

      ! dt_R4 = real(dt, kind = real_kind)    
      ! Tf_R4 = real(Tf, kind = real_kind)        

      hin_max_R4 = real(hin_max, kind = real_kind)   

      ! trcr_base_R4 = real(trcr_base, kind = real_kind)      

      aice_R4 = real(aice, kind = real_kind)     
      aice0_R4 = real(aice0, kind = real_kind)    
      rdg_conv_R4 = real(rdg_conv, kind = real_kind) 
      rdg_shear_R4 = real(rdg_shear, kind = real_kind)
      dardg1dt_R4 = real(dardg1dt, kind = real_kind) 
      dardg2dt_R4 = real(dardg2dt, kind = real_kind) 
      dvirdgdt_R4 = real(dvirdgdt, kind = real_kind) 
      opening_R4 = real(opening, kind = real_kind)  
      fpond_R4 = real(fpond, kind = real_kind)    
      fresh_R4 = real(fresh, kind = real_kind)    
      fsalt_R4 = real(fsalt, kind = real_kind)    
      fhocn_R4 = real(fhocn, kind = real_kind)    
      fzsal_R4 = real(fzsal, kind = real_kind)        

      aicen_R4 = real(aicen, kind = real_kind)    
      vicen_R4 = real(vicen, kind = real_kind)    
      vsnon_R4 = real(vsnon, kind = real_kind)    
      dardg1ndt_R4 = real(dardg1ndt, kind = real_kind)
      dardg2ndt_R4 = real(dardg2ndt, kind = real_kind)
      dvirdgndt_R4 = real(dvirdgndt, kind = real_kind)
      aparticn_R4 = real(aparticn, kind = real_kind) 
      krdgn_R4 = real(krdgn, kind = real_kind)    
      araftn_R4 = real(araftn, kind = real_kind)   
      vraftn_R4 = real(vraftn, kind = real_kind)   
      aredistn_R4 = real(aredistn, kind = real_kind) 
      vredistn_R4 = real(vredistn, kind = real_kind) 
      faero_ocn_R4 = real(faero_ocn, kind = real_kind)
      flux_bio_R4 = real(flux_bio, kind = real_kind)     

      trcrn_R4 = real(trcrn, kind = real_kind)        

      deallocate(trcr_base)
      deallocate(aicen)
      deallocate(vicen)
      deallocate(vsnon)
      deallocate(dardg1ndt)
      deallocate(dardg2ndt)
      deallocate(dvirdgndt)
      deallocate(aparticn)
      deallocate(krdgn)
      deallocate(araftn)
      deallocate(vraftn)
      deallocate(aredistn)
      deallocate(vredistn)
      deallocate(faero_ocn)
      deallocate(flux_bio)
      deallocate(trcrn)
      end subroutine colpkg_step_ridge

!=======================================================================

! Aggregate ice state variables over thickness categories.
!
! authors: C. M. Bitz, UW
!          W. H. Lipscomb, LANL

      subroutine colpkg_aggregate (ncat,     Tf_R4,       &
                                   aicen_R4,    trcrn_R4,    &
                                   vicen_R4,    vsnon_R4,    &
                                   aice_R4,     trcr_R4,     &
                                   vice_R4,     vsno_R4,     &
                                   aice0_R4,              &
                                   ntrcr,              &
                                   trcr_depend,        &
                                   trcr_base_R4,          & 
                                   n_trcr_strata,      &
                                   nt_strata)

      use ice_constants_colpkg, only: c0, c1
      use ice_colpkg_tracers, only: colpkg_compute_tracers

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         ntrcr     ! number of tracers in use

      real (kind=real_kind), intent(in) :: &
         Tf_R4        ! ocean freezing temperature           (Celsius)

      real (kind=real_kind), dimension (:), intent(in) :: &
         aicen_R4 , & ! concentration of ice
         vicen_R4 , & ! volume per unit area of ice          (m)
         vsnon_R4     ! volume per unit area of snow         (m)

      real (kind=real_kind), dimension (:,:), &
         intent(inout) :: &
         trcrn_R4     ! ice tracers

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=real_kind), dimension (:,:), intent(in) :: &
         trcr_base_R4      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      real (kind=real_kind), intent(out) :: &
         aice_R4  , & ! concentration of ice
         vice_R4  , & ! volume per unit area of ice          (m)
         vsno_R4  , & ! volume per unit area of snow         (m)
         aice0_R4     ! concentration of open water

      real (kind=real_kind), dimension (:),  &
         intent(out) :: &
         trcr_R4      ! ice tracers

      ! local variables

      integer (kind=int_kind) :: &
         n, it, itl, & ! loop indices
         ntr           ! tracer index

      real (kind=dbl_kind), dimension (:), allocatable :: &
         atrcr     ! sum of aicen*trcrn or vicen*trcrn or vsnon*trcrn

      real (kind=dbl_kind) :: &
         atrcrn    ! category value
      !-----------------------------------------------------------------
      ! local double copies
      !-----------------------------------------------------------------

      real (kind=dbl_kind) :: &
         Tf        ! ocean freezing temperature           (Celsius)

      real (kind=dbl_kind), dimension (:), allocatable :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:), allocatable :: &
         trcrn     ! ice tracers

      real (kind=dbl_kind), dimension (:,:), allocatable :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno
      !intent out vars 
      real (kind=dbl_kind) :: &
         aice  , & ! concentration of ice
         vice  , & ! volume per unit area of ice          (m)
         vsno  , & ! volume per unit area of snow         (m)
         aice0     ! concentration of open water

      real (kind=dbl_kind), dimension (:), allocatable :: &
         trcr      ! ice tracers

      Tf = real(Tf_R4, kind = dbl_kind)            
   
      trcrn  = real(trcrn_R4, kind = dbl_kind)      
      trcr_base  = real(trcr_base_R4, kind = dbl_kind)       
      vicen =real(vicen_R4, kind = dbl_kind)
      vsnon =real(vsnon_R4, kind = dbl_kind) 
      trcr =real(trcr_R4, kind = dbl_kind) 

      aice = real(aice_R4, kind = dbl_kind)   
      vice = real(vice_R4, kind = dbl_kind)   
      vsno = real(vsno_R4, kind = dbl_kind)   
      aice0 = real(aice0_R4, kind = dbl_kind)      
   
      trcr = real(trcr_R4, kind = dbl_kind)       
      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      aice0 = c1
      aice  = c0
      vice  = c0
      vsno  = c0

      allocate (atrcr(ntrcr))

      aicen =real(aicen_R4, kind = dbl_kind)


      ! aicen = real(aicen_R4, kind = dbl_kind)  
      ! vicen = real(vicen_R4, kind = dbl_kind)  
      ! vsnon = real(vsnon_R4, kind = dbl_kind)   
      !-----------------------------------------------------------------
      ! Aggregate
      !-----------------------------------------------------------------

      atrcr(:) = c0

      do n = 1, ncat

            aice = aice + aicen(n)
            vice = vice + vicen(n)
            vsno = vsno + vsnon(n)

         do it = 1, ntrcr
            atrcrn = trcrn(it,n)*(trcr_base(it,1) * aicen(n) &
                                + trcr_base(it,2) * vicen(n) &
                                + trcr_base(it,3) * vsnon(n))
            if (n_trcr_strata(it) > 0) then  ! additional tracer layers
               do itl = 1, n_trcr_strata(it)
                  ntr = nt_strata(it,itl)
                  atrcrn = atrcrn * trcrn(ntr,n)
               enddo
            endif
            atrcr(it) = atrcr(it) + atrcrn
         enddo                  ! ntrcr
      enddo                     ! ncat

      ! Open water fraction
      aice0 = max (c1 - aice, c0)

      ! Tracers
      call colpkg_compute_tracers (ntrcr,     trcr_depend,   &
                                   atrcr,     aice,          &
                                   vice ,     vsno,          &
                                   trcr_base, n_trcr_strata, &
                                   nt_strata, trcr,          &
                                   Tf)   

      deallocate (atrcr)

   

      trcrn_R4 = real(trcrn, kind = real_kind)      
      aice_R4 = real(aice, kind = real_kind)   
      vice_R4 = real(vice, kind = real_kind)   
      vsno_R4 = real(vsno, kind = real_kind)   
      aice0_R4 = real(aice0, kind = real_kind)      
      trcr_R4 = real(trcr, kind = real_kind)     

      deallocate (aicen) 
      deallocate(trcrn)
      deallocate(trcr_base)
      deallocate (vicen)
      deallocate (vsnon)  
      deallocate (trcr)  

      end subroutine colpkg_aggregate

!=======================================================================

! Compute the strength of the ice pack, defined as the energy (J m-2)
! dissipated per unit area removed from the ice pack under compression,
! and assumed proportional to the change in potential energy caused
! by ridging.
!
! See Rothrock (1975) and Hibler (1980).
!
! For simpler strength parameterization, see this reference:
! Hibler, W. D. III, 1979: A dynamic-thermodynamic sea ice model,
!  J. Phys. Oceanog., 9, 817-846.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine colpkg_ice_strength (ncat,               &
                                      aice_R4,     vice_R4,     &
                                      aice0_R4,    aicen_R4,    &
                                      vicen_R4,    &
                                      strength_R4)

      use ice_constants_colpkg, only: p333, c0, c1, c2, Cf, Cp, Pstar, Cstar, &
          rhoi, puny
      use ice_mechred, only: asum_ridging, ridge_itd

      integer (kind=int_kind), intent(in) :: & 
         ncat       ! number of thickness categories

      real (kind=real_kind), intent(in) :: &
         aice_R4   , & ! concentration of ice
         vice_R4  , & ! volume per unit area of ice  (m)
         aice0_R4      ! concentration of open water

      real (kind=real_kind), dimension(:), intent(in) :: &
         aicen_R4  , & ! concentration of ice
         vicen_R4      ! volume per unit area of ice  (m)

      real (kind=real_kind), intent(inout) :: &
         strength_R4   ! ice strength (N/m)

      ! local variables

      real (kind=dbl_kind) :: &
         aice   , & ! concentration of ice
         vice   , & ! volume per unit area of ice  (m)
         aice0      ! concentration of open water

      real (kind=dbl_kind), dimension(:), allocatable :: &
         aicen  , & ! concentration of ice
         vicen      ! volume per unit area of ice  (m)

      real (kind=dbl_kind)  :: &
         strength   ! ice strength (N/m)

      real (kind=dbl_kind) :: &
         asum   , & ! sum of ice and open water area
         aksum      ! ratio of area removed to area ridged

      real (kind=dbl_kind), dimension (0:ncat) :: &
         apartic    ! participation function; fraction of ridging
                    ! and closing associated w/ category n

      real (kind=dbl_kind), dimension (ncat) :: &
         hrmin  , & ! minimum ridge thickness
         hrmax  , & ! maximum ridge thickness (krdg_redist = 0)
         hrexp  , & ! ridge e-folding thickness (krdg_redist = 1) 
         krdg       ! mean ridge thickness/thickness of ridging ice

      integer (kind=int_kind) :: &
         n          ! thickness category index

      real (kind=dbl_kind) :: &
         hi     , & ! ice thickness (m)
         h2rdg  , & ! mean value of h^2 for new ridge
         dh2rdg     ! change in mean value of h^2 per unit area
                    ! consumed by ridging 
      aice = real(aice_R4, kind = dbl_kind)
      vice = real(vice_R4, kind = dbl_kind)
      aice0 = real(aice0_R4, kind = dbl_kind)
      aicen  = real(aicen_R4, kind = dbl_kind)
      vicen  = real(vicen_R4, kind = dbl_kind)
      strength = real(strength_R4, kind = dbl_kind)

      if (kstrength == 1) then  ! Rothrock '75 formulation

      !-----------------------------------------------------------------
      ! Compute thickness distribution of ridging and ridged ice.
      !-----------------------------------------------------------------

         call asum_ridging (ncat, aicen, aice0, asum)

         call ridge_itd (ncat,     aice0,      &
                         aicen,    vicen,      &
                         krdg_partic, krdg_redist, &
                         mu_rdg,                   &
                         aksum,    apartic,    &
                         hrmin,    hrmax,      &
                         hrexp,    krdg)   

      !-----------------------------------------------------------------
      ! Compute ice strength based on change in potential energy,
      ! as in Rothrock (1975)
      !-----------------------------------------------------------------

         if (krdg_redist==0) then ! Hibler 1980 formulation

            do n = 1, ncat
               if (aicen(n) > puny .and. apartic(n) > c0)then
                  hi = vicen(n) / aicen(n)
                  h2rdg = p333 * (hrmax(n)**3 - hrmin(n)**3)  &
                               / (hrmax(n) - hrmin(n)) 
                  dh2rdg = -hi*hi + h2rdg/krdg(n)
                  strength = strength + apartic(n) * dh2rdg
               endif         ! aicen > puny
            enddo               ! n

         elseif (krdg_redist==1) then ! exponential formulation

            do n = 1, ncat
               if (aicen(n) > puny .and. apartic(n) > c0) then
                  hi = vicen(n) / aicen(n)
                  h2rdg =    hrmin(n)*hrmin(n) &
                        + c2*hrmin(n)*hrexp(n) &
                        + c2*hrexp(n)*hrexp(n)
                  dh2rdg = -hi*hi + h2rdg/krdg(n)
                  strength = strength + apartic(n) * dh2rdg
               endif
            enddo               ! n

         endif                  ! krdg_redist

         strength = Cf * Cp * strength / aksum
                       ! Cp = (g/2)*(rhow-rhoi)*(rhoi/rhow)
                       ! Cf accounts for frictional dissipation

      else                      ! kstrength /= 1:  Hibler (1979) form

      !-----------------------------------------------------------------
      ! Compute ice strength as in Hibler (1979)
      !-----------------------------------------------------------------

         strength = Pstar*vice*exp(-Cstar*(c1-aice))

      endif                     ! kstrength
      deallocate(aicen)
      deallocate(vicen)
      strength_R4 = real(strength, kind = real_kind)

      end subroutine colpkg_ice_strength

!=======================================================================

      subroutine colpkg_atm_boundary(sfctype,                    &
                                     Tsf_R4,         potT_R4,          &
                                     uatm_R4,        vatm_R4,          &
                                     wind_R4,        zlvl_R4,          &
                                     Qa_R4,          rhoa_R4,          &
                                     strx_R4,        stry_R4,          &
                                     Tref_R4,        Qref_R4,          &
                                     delt_R4,        delq_R4,          &
                                     lhcoef_R4,      shcoef_R4,        &
                                     Cdn_atm_R4,                    &
                                     Cdn_atm_ratio_n_R4,            &
                                     uvel_R4,        vvel_R4,          &
                                     Uref_R4)

      use ice_atmo, only: atmo_boundary_const, atmo_boundary_layer
      use ice_constants_colpkg, only: c0

      character (len=3), intent(in) :: &
         sfctype      ! ice or ocean
      real (kind=real_kind), intent(in) :: &
         Tsf_R4      , & ! surface temperature of ice or ocean
         potT_R4     , & ! air potential temperature  (K)
         uatm_R4     , & ! x-direction wind speed (m/s)
         vatm_R4     , & ! y-direction wind speed (m/s)
         wind_R4     , & ! wind speed (m/s)
         zlvl_R4     , & ! atm level height (m)
         Qa_R4       , & ! specific humidity (kg/kg)
         rhoa_R4         ! air density (kg/m^3)

      real (kind=real_kind), intent(inout) :: &
         Cdn_atm_R4  , &    ! neutral drag coefficient
         Cdn_atm_ratio_n_R4 ! ratio drag coeff / neutral drag coeff

      real (kind=real_kind), &
         intent(inout) :: &
         strx_R4     , & ! x surface stress (N)
         stry_R4        ! y surface stress (N)

      real (kind=real_kind), intent(inout) :: &
         Tref_R4     , & ! reference height temperature  (K)
         Qref_R4     , & ! reference height specific humidity (kg/kg)
         delt_R4     , & ! potential T difference   (K)
         delq_R4     , & ! humidity difference      (kg/kg)
         shcoef_R4   , & ! transfer coefficient for sensible heat
         lhcoef_R4       ! transfer coefficient for latent heat

      real (kind=real_kind), optional, intent(in) :: &
         uvel_R4     , & ! x-direction ice speed (m/s)
         vvel_R4         ! y-direction ice speed (m/s)

      real (kind=real_kind), optional, intent(out) :: &
         Uref_R4         ! reference height wind speed (m/s)

      real (kind=dbl_kind) :: &
         worku, workv, workr

      real (kind=dbl_kind) :: &
         Tsf      , & ! surface temperature of ice or ocean
         potT     , & ! air potential temperature  (K)
         uatm     , & ! x-direction wind speed (m/s)
         vatm     , & ! y-direction wind speed (m/s)
         wind     , & ! wind speed (m/s)
         zlvl     , & ! atm level height (m)
         Qa       , & ! specific humidity (kg/kg)
         rhoa         ! air density (kg/m^3)

      real (kind=dbl_kind) :: &
         Cdn_atm  , &    ! neutral drag coefficient
         Cdn_atm_ratio_n ! ratio drag coeff / neutral drag coeff

      real (kind=dbl_kind) :: &
         strx     , & ! x surface stress (N)
         stry        ! y surface stress (N)

      real (kind=dbl_kind) :: &
         Tref     , & ! reference height temperature  (K)
         Qref     , & ! reference height specific humidity (kg/kg)
         delt     , & ! potential T difference   (K)
         delq     , & ! humidity difference      (kg/kg)
         shcoef   , & ! transfer coefficient for sensible heat
         lhcoef       ! transfer coefficient for latent heat

      real (kind=dbl_kind) :: &
         uvel     , & ! x-direction ice speed (m/s)
         vvel         ! y-direction ice speed (m/s)

      real (kind=dbl_kind) :: &
         Uref         ! reference height wind speed (m/s)     

      
      Tsf = real(Tsf_R4, kind = dbl_kind)      
      potT = real(potT_R4, kind = dbl_kind)     
      uatm = real(uatm_R4, kind = dbl_kind)     
      vatm = real(vatm_R4, kind = dbl_kind)     
      wind = real(wind_R4, kind = dbl_kind)     
      zlvl = real(zlvl_R4, kind = dbl_kind)     
      Qa = real(Qa_R4, kind = dbl_kind)       
      rhoa = real(rhoa_R4, kind = dbl_kind)         

   
      Cdn_atm = real(Cdn_atm_R4, kind = dbl_kind)  
      Cdn_atm_ratio_n = real(Cdn_atm_ratio_n_R4, kind = dbl_kind) 

   
      strx = real(strx_R4, kind = dbl_kind)     
      stry = real(stry_R4, kind = dbl_kind)        

   
      Tref = real(Tref_R4, kind = dbl_kind)     
      Qref = real(Qref_R4, kind = dbl_kind)     
      delt = real(delt_R4, kind = dbl_kind)     
      delq = real(delq_R4, kind = dbl_kind)     
      shcoef = real(shcoef_R4, kind = dbl_kind)   
      lhcoef = real(lhcoef_R4, kind = dbl_kind)       

      
      uvel = real(uvel_R4, kind = dbl_kind)     
      vvel = real(vvel_R4, kind = dbl_kind)         

   
      Uref = real(Uref_R4, kind = dbl_kind)         

      worku = c0
      workv = c0
      workr = c0
      if (present(uvel_R4)) then
         uvel = real(uvel_R4, kind = dbl_kind)     
      endif
      ! should this be for vvel,workv?
      if (present(vvel_R4)) then
         vvel = real(vvel_R4, kind = dbl_kind)         
      endif
      if (present(Uref_R4)) then
         Uref = real(Uref_R4, kind = dbl_kind)         
      endif

      if (present(uvel_R4)) then
         worku = uvel_R4
      endif
      ! should this be for vvel,workv?
      if (present(uvel_R4)) then
         worku = uvel_R4
      endif

               if (trim(atmbndy) == 'constant') then
                  call atmo_boundary_const (sfctype,  calc_strair, &
                                            uatm,     vatm,     &
                                            wind,     rhoa,     &
                                            strx,     stry,     &
                                            Tsf,      potT,     &
                                            Qa,                 &
                                            delt,     delq,     &
                                            lhcoef,   shcoef,   &
                                            Cdn_atm)
               else ! default
                  call atmo_boundary_layer (sfctype,                 &
                                            calc_strair, formdrag,   &
                                            highfreq, natmiter,      &
                                            Tsf,      potT,          &
                                            uatm,     vatm,          &
                                            wind,     zlvl,          &
                                            Qa,       rhoa,          &
                                            strx,     stry,          &
                                            Tref,     Qref,          &
                                            delt,     delq,          &
                                            lhcoef,   shcoef,        &
                                            Cdn_atm,                 &
                                            Cdn_atm_ratio_n,         &
                                            worku,    workv,         &
                                            workr)
               endif ! atmbndy

      if (present(Uref_R4)) then
         Uref_R4 = real(workr, kind=real_kind)
      endif

      Cdn_atm_R4 = real(Cdn_atm, kind = real_kind)  
      Cdn_atm_ratio_n_R4 = real(Cdn_atm_ratio_n, kind = real_kind) 

      strx_R4 = real(strx, kind = real_kind)     
      stry_R4 = real(stry, kind = real_kind)        

      Tref_R4 = real(Tref, kind = real_kind)     
      Qref_R4 = real(Qref, kind = real_kind)     
      delt_R4 = real(delt, kind = real_kind)     
      delq_R4 = real(delq, kind = real_kind)     
      shcoef_R4 = real(shcoef, kind = real_kind)   
      lhcoef_R4 = real(lhcoef, kind = real_kind)      
   
      if (present(Uref_R4)) then
         Uref_R4 = real(Uref, kind = real_kind)         
      endif
      end subroutine colpkg_atm_boundary
!=======================================================================

      subroutine colpkg_atm_boundary_double(sfctype,                    &
         Tsf,         potT,          &
         uatm,        vatm,          &
         wind,        zlvl,          &
         Qa,          rhoa,          &
         strx,        stry,          &
         Tref,        Qref,          &
         delt,        delq,          &
         lhcoef,      shcoef,        &
         Cdn_atm,                    &
         Cdn_atm_ratio_n,            &
         uvel,        vvel,          &
         Uref)

      use ice_atmo, only: atmo_boundary_const, atmo_boundary_layer
      use ice_constants_colpkg, only: c0

      character (len=3), intent(in) :: &
      sfctype      ! ice or ocean

      real (kind=dbl_kind), intent(in) :: &
      Tsf      , & ! surface temperature of ice or ocean
      potT     , & ! air potential temperature  (K)
      uatm     , & ! x-direction wind speed (m/s)
      vatm     , & ! y-direction wind speed (m/s)
      wind     , & ! wind speed (m/s)
      zlvl     , & ! atm level height (m)
      Qa       , & ! specific humidity (kg/kg)
      rhoa         ! air density (kg/m^3)

      real (kind=dbl_kind), intent(inout) :: &
      Cdn_atm  , &    ! neutral drag coefficient
      Cdn_atm_ratio_n ! ratio drag coeff / neutral drag coeff

      real (kind=dbl_kind), &
      intent(inout) :: &
      strx     , & ! x surface stress (N)
      stry         ! y surface stress (N)

      real (kind=dbl_kind), intent(inout) :: &
      Tref     , & ! reference height temperature  (K)
      Qref     , & ! reference height specific humidity (kg/kg)
      delt     , & ! potential T difference   (K)
      delq     , & ! humidity difference      (kg/kg)
      shcoef   , & ! transfer coefficient for sensible heat
      lhcoef       ! transfer coefficient for latent heat

      real (kind=dbl_kind), optional, intent(in) :: &
      uvel     , & ! x-direction ice speed (m/s)
      vvel         ! y-direction ice speed (m/s)

      real (kind=dbl_kind), optional, intent(out) :: &
      Uref         ! reference height wind speed (m/s)

      real (kind=dbl_kind) :: &
      worku, workv, workr

      worku = c0
      workv = c0
      workr = c0
      if (present(uvel)) then
      worku = uvel
      endif
      ! should this be for vvel,workv?
      if (present(uvel)) then
      worku = uvel
      endif

      if (trim(atmbndy) == 'constant') then
      call atmo_boundary_const (sfctype,  calc_strair, &
                     uatm,     vatm,     &
                     wind,     rhoa,     &
                     strx,     stry,     &
                     Tsf,      potT,     &
                     Qa,                 &
                     delt,     delq,     &
                     lhcoef,   shcoef,   &
                     Cdn_atm)
      else ! default
      call atmo_boundary_layer (sfctype,                 &
                     calc_strair, formdrag,   &
                     highfreq, natmiter,      &
                     Tsf,      potT,          &
                     uatm,     vatm,          &
                     wind,     zlvl,          &
                     Qa,       rhoa,          &
                     strx,     stry,          &
                     Tref,     Qref,          &
                     delt,     delq,          &
                     lhcoef,   shcoef,        &
                     Cdn_atm,                 &
                     Cdn_atm_ratio_n,         &
                     worku,    workv,         &
                     workr)
      endif ! atmbndy

      if (present(Uref)) then
      Uref = workr
      endif

      end subroutine colpkg_atm_boundary_double

!=======================================================================

! Compute the energy available to freeze or melt ice.
! NOTE: SST changes due to fluxes through the ice are computed in
!       ice_therm_vertical.

      subroutine colpkg_ocn_mixed_layer (alvdr_ocn_R4, swvdr_R4,      &
                                         alidr_ocn_R4, swidr_R4,      &
                                         alvdf_ocn_R4, swvdf_R4,      &
                                         alidf_ocn_R4, swidf_R4,      &
                                         sst_R4,       flwout_ocn_R4, &
                                         fsens_ocn_R4, shcoef_R4,     &
                                         flat_ocn_R4,  lhcoef_R4,     &
                                         evap_ocn_R4,  flw_R4,        &
                                         delt_R4,      delq_R4,       &
                                         aice_R4,      fhocn_R4,      &
                                         fswthru_R4,   hmix_R4,       &
                                         Tf_R4,        qdp_R4,        &
                                         frzmlt_R4,    dt_R4)

      use ice_constants_colpkg, only: c0, c1, c1000, &
          cp_ocn, Tffresh, stefan_boltzmann, Lvap, cprho

      real (kind=real_kind), intent(in) :: &
         alvdr_ocn_R4 , & ! visible, direct   (fraction)
         alidr_ocn_R4 , & ! near-ir, direct   (fraction)
         alvdf_ocn_R4 , & ! visible, diffuse  (fraction)
         alidf_ocn_R4 , & ! near-ir, diffuse  (fraction)
         swvdr_R4     , & ! sw down, visible, direct  (W/m^2)
         swvdf_R4     , & ! sw down, visible, diffuse (W/m^2)
         swidr_R4     , & ! sw down, near IR, direct  (W/m^2)
         swidf_R4     , & ! sw down, near IR, diffuse (W/m^2)
         flw_R4       , & ! incoming longwave radiation (W/m^2)
         Tf_R4        , & ! freezing temperature (C)
         hmix_R4      , & ! mixed layer depth (m)
         delt_R4      , & ! potential temperature difference   (K)
         delq_R4      , & ! specific humidity difference   (kg/kg)
         shcoef_R4    , & ! transfer coefficient for sensible heat
         lhcoef_R4    , & ! transfer coefficient for latent heat
         fhocn_R4     , & ! net heat flux to ocean (W/m^2)
         fswthru_R4   , & ! shortwave penetrating to ocean (W/m^2)
         aice_R4      , & ! ice area fraction
         dt_R4            ! time step (s)

      real (kind=real_kind), intent(inout) :: &
         flwout_ocn_R4, & ! outgoing longwave radiation (W/m^2)
         fsens_ocn_R4 , & ! sensible heat flux (W/m^2)
         flat_ocn_R4  , & ! latent heat flux   (W/m^2)
         evap_ocn_R4  , & ! evaporative water flux (kg/m^2/s)
         qdp_R4       , & ! deep ocean heat flux (W/m^2), negative upward
         sst_R4       , & ! sea surface temperature (C)
         frzmlt_R4        ! freezing/melting potential (W/m^2)

      ! local variables

      real (kind=dbl_kind), parameter :: &
         frzmlt_max = c1000   ! max magnitude of frzmlt (W/m^2)

      real (kind=dbl_kind) :: &
         TsfK , & ! surface temperature (K)
         swabs    ! surface absorbed shortwave heat flux (W/m^2)

      real (kind=dbl_kind) :: &
         alvdr_ocn , & ! visible, direct   (fraction)
         alidr_ocn , & ! near-ir, direct   (fraction)
         alvdf_ocn , & ! visible, diffuse  (fraction)
         alidf_ocn , & ! near-ir, diffuse  (fraction)
         swvdr     , & ! sw down, visible, direct  (W/m^2)
         swvdf     , & ! sw down, visible, diffuse (W/m^2)
         swidr     , & ! sw down, near IR, direct  (W/m^2)
         swidf     , & ! sw down, near IR, diffuse (W/m^2)
         flw       , & ! incoming longwave radiation (W/m^2)
         Tf        , & ! freezing temperature (C)
         hmix      , & ! mixed layer depth (m)
         delt      , & ! potential temperature difference   (K)
         delq      , & ! specific humidity difference   (kg/kg)
         shcoef    , & ! transfer coefficient for sensible heat
         lhcoef    , & ! transfer coefficient for latent heat
         fhocn     , & ! net heat flux to ocean (W/m^2)
         fswthru   , & ! shortwave penetrating to ocean (W/m^2)
         aice      , & ! ice area fraction
         dt            ! time step (s)

      real (kind=dbl_kind) :: &
         flwout_ocn, & ! outgoing longwave radiation (W/m^2)
         fsens_ocn , & ! sensible heat flux (W/m^2)
         flat_ocn  , & ! latent heat flux   (W/m^2)
         evap_ocn  , & ! evaporative water flux (kg/m^2/s)
         qdp       , & ! deep ocean heat flux (W/m^2), negative upward
         sst       , & ! sea surface temperature (C)
         frzmlt        ! freezing/melting potential (W/m^2)


      alvdr_ocn = real(alvdr_ocn_R4, kind = dbl_kind) 
      alidr_ocn = real(alidr_ocn_R4, kind = dbl_kind) 
      alvdf_ocn = real(alvdf_ocn_R4, kind = dbl_kind) 
      alidf_ocn = real(alidf_ocn_R4, kind = dbl_kind) 
      swvdr = real(swvdr_R4, kind = dbl_kind)     
      swvdf = real(swvdf_R4, kind = dbl_kind)     
      swidr = real(swidr_R4, kind = dbl_kind)     
      swidf = real(swidf_R4, kind = dbl_kind)     
      flw = real(flw_R4, kind = dbl_kind)       
      Tf = real(Tf_R4, kind = dbl_kind)        
      hmix = real(hmix_R4, kind = dbl_kind)      
      delt = real(delt_R4, kind = dbl_kind)      
      delq = real(delq_R4, kind = dbl_kind)      
      shcoef = real(shcoef_R4, kind = dbl_kind)    
      lhcoef = real(lhcoef_R4, kind = dbl_kind)    
      fhocn = real(fhocn_R4, kind = dbl_kind)     
      fswthru = real(fswthru_R4, kind = dbl_kind)   
      aice = real(aice_R4, kind = dbl_kind)      
      dt = real(dt_R4, kind = dbl_kind)            

      flwout_ocn = real(flwout_ocn_R4, kind = dbl_kind)
      fsens_ocn = real(fsens_ocn_R4, kind = dbl_kind) 
      flat_ocn = real(flat_ocn_R4, kind = dbl_kind)  
      evap_ocn = real(evap_ocn_R4, kind = dbl_kind)  
      qdp = real(qdp_R4, kind = dbl_kind)       
      sst = real(sst_R4, kind = dbl_kind)       
      frzmlt = real(frzmlt_R4, kind = dbl_kind)        

      ! shortwave radiative flux
      swabs = (c1-alvdr_ocn) * swvdr + (c1-alidr_ocn) * swidr &
            + (c1-alvdf_ocn) * swvdf + (c1-alidf_ocn) * swidf 

      ! ocean surface temperature in Kelvin
      TsfK = sst + Tffresh

      ! longwave radiative flux
      flwout_ocn = -stefan_boltzmann * TsfK**4

      ! downward latent and sensible heat fluxes
      fsens_ocn =  shcoef * delt
      flat_ocn  =  lhcoef * delq
      evap_ocn  = -flat_ocn / Lvap

      ! Compute sst change due to exchange with atm/ice above
      sst = sst + dt * ( &
            (fsens_ocn + flat_ocn + flwout_ocn + flw + swabs) * (c1-aice) &
          + fhocn + fswthru)         &  ! these are *aice
          / (cprho*hmix)

      ! adjust qdp if cooling of mixed layer would occur when sst <= Tf
      if (sst <= Tf .and. qdp > c0) qdp = c0

      ! computed T change due to exchange with deep layers:
      sst = sst - qdp*dt/(cprho*hmix)

      ! compute potential to freeze or melt ice
      frzmlt = (Tf-sst)*cprho*hmix/dt
      frzmlt = min(max(frzmlt,-frzmlt_max),frzmlt_max)

      ! if sst is below freezing, reset sst to Tf
      if (sst <= Tf) sst = Tf


      flwout_ocn_R4 = real(flwout_ocn, kind = real_kind)
      fsens_ocn_R4 = real(fsens_ocn, kind = real_kind) 
      flat_ocn_R4 = real(flat_ocn, kind = real_kind)  
      evap_ocn_R4 = real(evap_ocn, kind = real_kind)  
      qdp_R4 = real(qdp, kind = real_kind)       
      sst_R4 = real(sst, kind = real_kind)       
      frzmlt_R4 = real(frzmlt, kind = real_kind)        

      end subroutine colpkg_ocn_mixed_layer

!=======================================================================
!
! Updates snow tracers
!
! authors: Elizabeth C. Hunke, LANL
!          Nicole Jeffery, LANL

      subroutine colpkg_step_snow (dt_R4,       wind_R4,   &
                                   nilyr,              &
                                   nslyr,    ncat,    &
                                   aice_R4,     aicen_R4,   &
                                   vicen_R4,    vsnon_R4,   &
                                   alvl_R4,     vlvl_R4,    &
                                   smice_R4,    smliq_R4,   &
                                   rhos_effn_R4,rhos_eff_R4,&
                                   rhos_cmpn_R4,rhos_cmp_R4,&
                                   rsnw_R4,     zqin1_R4,   &
                                   zSin1_R4,    Tsfc_R4,    &
                                   zqsn_R4,               &
                                   fresh_R4,    fhocn_R4,   &
                                   fsloss_R4,   fsnow_R4,   &
                                   rhosnew_R4,  rhosmax_R4, &
                                   windmin_R4,  drhosdwind_R4,&
                                   snowage_tau_R4,&
                                   snowage_kappa_R4,&
                                   snowage_drdt0_R4,&
                                   idx_T_max,&
                                   idx_Tgrd_max,&
                                   idx_rhos_max,&
                                   l_stop,&
                                   stop_label)

      use ice_colpkg_tracers, only: tr_snow, tr_rsnw
      use ice_constants_colpkg, only: c0, puny, rhos
      use ice_snow, only: snow_effective_density, update_snow_radius, &
                          snow_redist

      integer (kind=int_kind), intent(in) :: & 
         nslyr, & ! number of snow layers
         nilyr, & ! number of ice  layers
         ncat, &  ! number of thickness categories
         idx_T_max, & ! dimensions of snow parameter matrix
         idx_Tgrd_max, &
         idx_rhos_max

      real (kind=real_kind), intent(in) :: &
         dt_R4     , & ! time step
         wind_R4   , & ! wind speed (m/s)
         fsnow_R4  , & ! snowfall rate (kg m-2 s-1)
         aice_R4   , & ! ice area fraction
         rhosnew_R4, & ! new snow density (kg/m^3)
         rhosmax_R4, & ! maximum snow density (kg/m^3)
         windmin_R4, & ! minimum wind speed to compact snow (m/s)
         drhosdwind_R4 ! wind compaction factor (kg s/m^4)

      real (kind=real_kind), dimension(:), intent(in) :: &
         aicen_R4, & ! ice area fraction
         vicen_R4, & ! ice volume (m)
         Tsfc_R4 , & ! surface temperature (C)
         zqin1_R4, & ! ice upper layer enthalpy
         zSin1_R4, & ! ice upper layer salinity
         alvl_R4,  & ! level ice area tracer
         vlvl_R4     ! level ice volume tracer

      real (kind=real_kind), intent(inout) :: &
         fresh_R4    , & ! fresh water flux to ocean (kg/m^2/s)
         fhocn_R4    , & ! net heat flux to ocean (W/m^2)
         fsloss_R4       ! snow loss to leads (kg/m^2/s)

      real (kind=real_kind), dimension(:), intent(inout) :: &
         vsnon_R4    ! snow volume (m)

      real (kind=real_kind), dimension(:,:), intent(inout) :: &
         zqsn_R4     , & ! snow enthalpy (J/m^3)
         smice_R4    , & ! mass of ice in snow (kg/m^3)
         smliq_R4    , & ! mass of liquid in snow (kg/m^3)
         rsnw_R4     , & ! snow grain radius (10^-6 m)
         rhos_effn_R4, & ! effective snow density: content (kg/m^3)
         rhos_cmpn_R4    ! effective snow density: compaction (kg/m^3)

      real (kind=real_kind), intent(inout) :: &
         rhos_eff_R4 , & ! mean effective snow density: content (kg/m^3)
         rhos_cmp_R4     ! mean effective snow density: compaction (kg/m^3)

      ! dry snow aging parameters
      real (kind=real_kind), dimension(idx_rhos_max,idx_Tgrd_max,idx_T_max), intent(in) :: &  
         snowage_tau_R4,   & ! (10^-6 m)
         snowage_kappa_R4, & ! 
         snowage_drdt0_R4    ! (10^-6 m/hr)

      logical (kind=log_kind), intent(inout) :: &
         l_stop          ! if true, print diagnostics and abort model

      character (len=*), intent(out) :: &
         stop_label   ! abort error message

      ! local temporary variables

      integer (kind=int_kind) :: n

      real (kind=dbl_kind), dimension(ncat) :: &
         zTin,  & ! ice upper layer temperature (oC)
         hsn ,  & ! snow thickness (m)
         hin      ! ice thickness

      real (kind=dbl_kind) :: &
         vsno,  & ! snow volume (m)
         tmp1, tmp2

      character(len=char_len_long) :: &
           warning ! warning message

      real (kind=dbl_kind) :: &
         dt     , & ! time step
         wind   , & ! wind speed (m/s)
         fsnow  , & ! snowfall rate (kg m-2 s-1)
         aice   , & ! ice area fraction
         rhosnew, & ! new snow density (kg/m^3)
         rhosmax, & ! maximum snow density (kg/m^3)
         windmin, & ! minimum wind speed to compact snow (m/s)
         drhosdwind ! wind compaction factor (kg s/m^4)

      real (kind=dbl_kind), dimension(:), allocatable :: &
         aicen, & ! ice area fraction
         vicen, & ! ice volume (m)
         Tsfc , & ! surface temperature (C)
         zqin1, & ! ice upper layer enthalpy
         zSin1, & ! ice upper layer salinity
         alvl,  & ! level ice area tracer
         vlvl     ! level ice volume tracer

      real (kind=dbl_kind) :: &
         fresh    , & ! fresh water flux to ocean (kg/m^2/s)
         fhocn    , & ! net heat flux to ocean (W/m^2)
         fsloss       ! snow loss to leads (kg/m^2/s)

      real (kind=dbl_kind), dimension(:), allocatable :: &
         vsnon    ! snow volume (m)

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         zqsn     , & ! snow enthalpy (J/m^3)
         smice    , & ! mass of ice in snow (kg/m^3)
         smliq    , & ! mass of liquid in snow (kg/m^3)
         rsnw     , & ! snow grain radius (10^-6 m)
         rhos_effn, & ! effective snow density: content (kg/m^3)
         rhos_cmpn    ! effective snow density: compaction (kg/m^3)

      real (kind=dbl_kind) :: &
         rhos_eff , & ! mean effective snow density: content (kg/m^3)
         rhos_cmp     ! mean effective snow density: compaction (kg/m^3)

      ! dry snow aging parameters
      real (kind=dbl_kind), dimension(idx_rhos_max,idx_Tgrd_max,idx_T_max) :: &
         snowage_tau,   & ! (10^-6 m)
         snowage_kappa, & ! 
         snowage_drdt0    ! (10^-6 m/hr)


      dt      = real(dt_R4     , kind = dbl_kind)
      wind    = real(wind_R4   , kind = dbl_kind)
      fsnow   = real(fsnow_R4  , kind = dbl_kind)
      aice    = real(aice_R4   , kind = dbl_kind)
      rhosnew = real(rhosnew_R4, kind = dbl_kind)
      rhosmax = real(rhosmax_R4, kind = dbl_kind)
      windmin = real(windmin_R4, kind = dbl_kind)
      drhosdwind  = real(drhosdwind_R4 , kind = dbl_kind)


      fresh     = real(fresh_R4    , kind = dbl_kind)
      fhocn     = real(fhocn_R4    , kind = dbl_kind)
      fsloss        = real(fsloss_R4       , kind = dbl_kind)

     aicen  = real(aicen_R4, kind = dbl_kind)
     vicen  = real(vicen_R4, kind = dbl_kind)
     Tsfc   = real(Tsfc_R4 , kind = dbl_kind)
     zqin1  = real(zqin1_R4, kind = dbl_kind)
     zSin1  = real(zSin1_R4, kind = dbl_kind)
     alvl  = real(alvl_R4, kind = dbl_kind)
     vlvl       = real(vlvl_R4     , kind = dbl_kind)

     vsnon      = real(vsnon_R4    , kind = dbl_kind)

     zqsn       = real(zqsn_R4     , kind = dbl_kind)
     smice      = real(smice_R4    , kind = dbl_kind)
     smliq      = real(smliq_R4    , kind = dbl_kind)
     rsnw       = real(rsnw_R4     , kind = dbl_kind)
     rhos_effn  = real(rhos_effn_R4, kind = dbl_kind)
     rhos_cmpn      = real(rhos_cmpn_R4    , kind = dbl_kind)

      rhos_eff  = real(rhos_eff_R4 , kind = dbl_kind)
      rhos_cmp      = real(rhos_cmp_R4     , kind = dbl_kind)


      snowage_tau = real(snowage_tau_R4, kind = dbl_kind)
      snowage_kappa = real(snowage_kappa_R4, kind = dbl_kind)
      snowage_drdt0 = real(snowage_drdt0_R4, kind = dbl_kind)    

      l_stop = .false.
      stop_label = ''

      if (tr_snow) then

      !-----------------------------------------------------------------
      ! Compute effective density of snow
      !-----------------------------------------------------------------

      vsno = c0
      do n = 1, ncat
         vsno = vsno + vsnon(n)
      enddo

      call snow_effective_density(nslyr,     ncat,     &
                                  vsnon,     vsno,     &
                                  smice,     smliq,    &
                                  rhosnew,             &
                                  rhos_effn, rhos_eff, &
                                  rhos_cmpn, rhos_cmp)

      !-----------------------------------------------------------------
      ! Redistribute snow based on wind
      !-----------------------------------------------------------------

      tmp1 = rhos*vsno + fresh*dt

      if (snwredist(1:3) == 'ITD' .and. aice > puny) then
         call snow_redist(dt,                  &
                          nslyr,    ncat,      &
                          wind,     aicen(:),  &
                          vicen(:), vsnon(:),  &
                          zqsn(:,:),snwredist, &
                          alvl(:),  vlvl(:),   &
                          fresh,    fhocn,     &
                          fsloss,   rhos_cmpn, &
                          fsnow,    rhosmax,   &
                          windmin,  drhosdwind, &
                          l_stop,   stop_label)
      endif

      vsno = c0
      do n = 1, ncat
         vsno = vsno + vsnon(n)
      enddo
      tmp2 = rhos*vsno + fresh*dt
      if (abs(tmp1-tmp2)>puny) then
        write(warning,*) ' '
        call add_warning(warning)
        write(warning,*)'tmp1 ne tmp2',tmp1, tmp2
        call add_warning(warning)
        stop_label ='snow redistribution error' 
        l_stop = .true.
      endif

      endif ! tr_snow

      !-----------------------------------------------------------------
      ! Adjust snow grain radius
      !-----------------------------------------------------------------

      if (tr_rsnw) then
       do n = 1, ncat
          zTin(n)= c0
          hsn(n) = c0
          hin(n) = c0
          if (aicen(n) > puny) then
              zTin(n)  = colpkg_ice_temperature(zqin1(n),zSin1(n))
              hsn(n)   = vsnon(n)/aicen(n)
              hin(n)   = vicen(n)/aicen(n)
          endif
       enddo

       call update_snow_radius       (dt,         ncat,  &
                                      nslyr,      nilyr, &
                                      rsnw,       hin,   &
                                      Tsfc,       zTin,  &
                                      hsn,        zqsn,  &
                                      smice,      smliq, &
                                      rsnw_fall,  rsnw_tmax, &
                                      snowage_tau, &
                                      snowage_kappa, &
                                      snowage_drdt0, &
                                      idx_T_max, &
                                      idx_Tgrd_max, &
                                      idx_rhos_max)
      endif

      ! dt_R4      = real(dt     , kind = real_kind)
      ! wind_R4    = real(wind   , kind = real_kind)
      ! fsnow_R4   = real(fsnow  , kind = real_kind)
      ! aice_R4    = real(aice   , kind = real_kind)
      ! rhosnew_R4 = real(rhosnew, kind = real_kind)
      ! rhosmax_R4 = real(rhosmax, kind = real_kind)
      ! windmin_R4 = real(windmin, kind = real_kind)
      ! drhosdwind_R4  = real(drhosdwind , kind = real_kind)

      ! aicen_R4 = real(aicen, kind = real_kind)
      ! vicen_R4 = real(vicen, kind = real_kind)
      ! Tsfc_R4  = real(Tsfc , kind = real_kind)
      ! zqin1_R4 = real(zqin1, kind = real_kind)
      ! zSin1_R4 = real(zSin1, kind = real_kind)
      ! alvl_R4 = real(alvl, kind = real_kind)
      ! vlvl_R4      = real(vlvl     , kind = real_kind)

      fresh_R4     = real(fresh    , kind = real_kind)
      fhocn_R4     = real(fhocn    , kind = real_kind)
      fsloss_R4        = real(fsloss       , kind = real_kind)

      vsnon_R4     = real(vsnon    , kind = real_kind)

      zqsn_R4      = real(zqsn     , kind = real_kind)
      smice_R4     = real(smice    , kind = real_kind)
      smliq_R4     = real(smliq    , kind = real_kind)
      rsnw_R4      = real(rsnw     , kind = real_kind)
      rhos_effn_R4 = real(rhos_effn, kind = real_kind)
      rhos_cmpn_R4     = real(rhos_cmpn    , kind = real_kind)

      rhos_eff_R4  = real(rhos_eff , kind = real_kind)
      rhos_cmp_R4      = real(rhos_cmp     , kind = real_kind)


      ! snowage_tau_R4 = real(snowage_tau, kind = real_kind)
      ! snowage_kappa_R4 = real(snowage_kappa, kind = real_kind)
      ! snowage_drdt0_R4 = real(snowage_drdt0, kind = real_kind)    

      deallocate(aicen)
      deallocate(vicen)
      deallocate(Tsfc)
      deallocate(zqin1)
      deallocate(zSin1)
      deallocate(alvl)
      deallocate(vlvl)
      deallocate(vsnon)
      deallocate(zqsn)
      deallocate(smice)
      deallocate(smliq)
      deallocate(rsnw)
      deallocate(rhos_effn)
      deallocate(rhos_cmpn)

      end subroutine colpkg_step_snow

!=======================================================================
! subroutine to set the column package internal parameters

      subroutine colpkg_init_parameters(&
           ktherm_in, &
           conduct_in, &
           fbot_xfer_type_in, &
           calc_Tsfc_in, &
           ustar_min_in, &
           a_rapid_mode_in, &
           Rac_rapid_mode_in, &
           aspect_rapid_mode_in, &
           dSdt_slow_mode_in, &
           phi_c_slow_mode_in, &
           phi_i_mushy_in, &
           shortwave_in, &
           use_snicar_in, &
           albedo_type_in, &
           albicev_in, &
           albicei_in, &
           albsnowv_in, &
           albsnowi_in, &
           ahmax_in, &
           R_ice_in, &
           R_pnd_in, &
           R_snw_in, &
           dT_mlt_in, &
           rsnw_mlt_in, &
           kalg_in, &
           kstrength_in, &
           krdg_partic_in, &
           krdg_redist_in, &
           mu_rdg_in, &
           Cf_in, &
           atmbndy_in, &
           calc_strair_in, &
           formdrag_in, &
           highfreq_in, &
           natmiter_in, &
           oceanmixed_ice_in, &
           tfrz_option_in, &
           kitd_in, &
           kcatbound_in, &
           hs0_in, &
           frzpnd_in, &
           dpscale_in, &
           rfracmin_in, &
           rfracmax_in, &
           pndaspect_in, &
           hs1_in, &
           hp1_in, &
         ! bgc_data_dir_in, &
         ! sil_data_type_in, &
         ! nit_data_type_in, &
         ! fe_data_type_in, &
           bgc_flux_type_in, &
           z_tracers_in, &
           scale_bgc_in, &
           solve_zbgc_in, &
           dEdd_algae_in, &
           modal_aero_in, &
           skl_bgc_in, &
           solve_zsal_in, &
           grid_o_in, &
           l_sk_in, &
           grid_o_t_in, &
           initbio_frac_in, &
           frazil_scav_in, &
           grid_oS_in, &
           l_skS_in, &
           phi_snow_in, &
           ratio_Si2N_diatoms_in, &
           ratio_Si2N_sp_in, &
           ratio_Si2N_phaeo_in, &
           ratio_S2N_diatoms_in, &
           ratio_S2N_sp_in, &      
           ratio_S2N_phaeo_in, &   
           ratio_Fe2C_diatoms_in, & 
           ratio_Fe2C_sp_in, &     
           ratio_Fe2C_phaeo_in, &  
           ratio_Fe2N_diatoms_in, & 
           ratio_Fe2N_sp_in, &     
           ratio_Fe2N_phaeo_in, &  
           ratio_Fe2DON_in, &       
           ratio_Fe2DOC_s_in, &     
           ratio_Fe2DOC_l_in, &     
           fr_resp_in, &            
           tau_min_in, &            
           tau_max_in, &            
           algal_vel_in, &          
           R_dFe2dust_in, &         
           dustFe_sol_in, &         
           chlabs_diatoms_in, &    
           chlabs_sp_in, &         
           chlabs_phaeo_in, &      
           alpha2max_low_diatoms_in, &  
           alpha2max_low_sp_in, &       
           alpha2max_low_phaeo_in, &    
           beta2max_diatoms_in, & 
           beta2max_sp_in, &       
           beta2max_phaeo_in, &    
           mu_max_diatoms_in, &   
           mu_max_sp_in, &         
           mu_max_phaeo_in, &      
           grow_Tdep_diatoms_in, &
           grow_Tdep_sp_in, &      
           grow_Tdep_phaeo_in, &   
           fr_graze_diatoms_in, & 
           fr_graze_sp_in, &       
           fr_graze_phaeo_in, &    
           mort_pre_diatoms_in, & 
           mort_pre_sp_in, &       
           mort_pre_phaeo_in, &    
           mort_Tdep_diatoms_in, &
           mort_Tdep_sp_in, &       
           mort_Tdep_phaeo_in, &    
           k_exude_diatoms_in, &  
           k_exude_sp_in, &         
           k_exude_phaeo_in, &      
           K_Nit_diatoms_in, &    
           K_Nit_sp_in, &           
           K_Nit_phaeo_in, &        
           K_Am_diatoms_in, &     
           K_Am_sp_in, &             
           K_Am_phaeo_in, &          
           K_Sil_diatoms_in, &    
           K_Sil_sp_in, &            
           K_Sil_phaeo_in, &         
           K_Fe_diatoms_in, &     
           K_Fe_sp_in, &             
           K_Fe_phaeo_in, &           
           f_don_protein_in, &    
           kn_bac_protein_in, &   
           f_don_Am_protein_in, & 
           f_doc_s_in, &            
           f_doc_l_in, &               
           f_exude_s_in, &          
           f_exude_l_in, &           
           k_bac_s_in, &            
           k_bac_l_in, &             
           T_max_in, &              
           fsal_in, &               
           op_dep_min_in, &         
           fr_graze_s_in, &         
           fr_graze_e_in, &         
           fr_mort2min_in, &        
           fr_dFe_in, &             
           k_nitrif_in, &           
           t_iron_conv_in, &        
           max_loss_in, &           
           max_dfe_doc1_in, &       
           fr_resp_s_in, &          
           y_sk_DMS_in, &           
           t_sk_conv_in, &          
           t_sk_ox_in, &             
           algaltype_diatoms_in, &   
           algaltype_sp_in, &       
           algaltype_phaeo_in, &    
           nitratetype_in, &        
           ammoniumtype_in, &       
           silicatetype_in, &       
           dmspptype_in, &          
           dmspdtype_in, &          
           humtype_in, &            
           doctype_s_in, &          
           doctype_l_in, &
           dictype_1_in, &
           dontype_protein_in, &     
           fedtype_1_in, &           
           feptype_1_in, &           
           zaerotype_bc1_in, &       
           zaerotype_bc2_in, &       
           zaerotype_dust1_in, &     
           zaerotype_dust2_in, &     
           zaerotype_dust3_in, &     
           zaerotype_dust4_in, &     
           ratio_C2N_diatoms_in, &   
           ratio_C2N_sp_in, &        
           ratio_C2N_phaeo_in, &     
           ratio_chl2N_diatoms_in, & 
           ratio_chl2N_sp_in, &      
           ratio_chl2N_phaeo_in, &   
           F_abs_chl_diatoms_in, &   
           F_abs_chl_sp_in, &        
           F_abs_chl_phaeo_in, &
           ratio_C2N_proteins_in, &
           snwredist_in, &
           use_smliq_pnd_in, &
           rsnw_fall_in, &
           rsnw_tmax_in, &
           rhosnew_in, &
           rhosmax_in, &
           windmin_in, &
           drhosdwind_in)
           !restore_bgc_in)

        use ice_colpkg_shared, only: &
             ktherm, &
             conduct, &
             fbot_xfer_type, &
             calc_Tsfc, &
             ustar_min, &
             a_rapid_mode, &
             Rac_rapid_mode, &
             aspect_rapid_mode, &
             dSdt_slow_mode, &
             phi_c_slow_mode, &
             phi_i_mushy, &
             shortwave, &
             use_snicar, &
             albedo_type, &
             albicev, &
             albicei, &
             albsnowv, &
             albsnowi, &
             ahmax, &
             R_ice, &
             R_pnd, &
             R_snw, &
             dT_mlt, &
             rsnw_mlt, &
             kalg, &
             kstrength, &
             krdg_partic, &
             krdg_redist, &
             mu_rdg, &
             Cf, &
             atmbndy, &
             calc_strair, &
             formdrag, &
             highfreq, &
             natmiter, &
             oceanmixed_ice, &
             tfrz_option, &
             kitd, &
             kcatbound, &
             hs0, &
             frzpnd, &
             dpscale, &
             rfracmin, &
             rfracmax, &
             pndaspect, &
             hs1, &
             hp1, &
           ! bgc_data_dir, &
           ! sil_data_type, &
           ! nit_data_type, &
           ! fe_data_type, &
             bgc_flux_type, &
             z_tracers, &
             scale_bgc, &
             solve_zbgc, &
             dEdd_algae, &
             modal_aero, &
             skl_bgc, &
             solve_zsal, &
             grid_o, &
             l_sk, &
             grid_o_t, &
             initbio_frac, &
             frazil_scav, &
             grid_oS, &
             l_skS, &
             phi_snow, &
             ratio_Si2N_diatoms, & 
             ratio_Si2N_sp     , &
             ratio_Si2N_phaeo  , &
             ratio_S2N_diatoms , & 
             ratio_S2N_sp      , &
             ratio_S2N_phaeo   , &
             ratio_Fe2C_diatoms, & 
             ratio_Fe2C_sp     , &
             ratio_Fe2C_phaeo  , &
             ratio_Fe2N_diatoms, & 
             ratio_Fe2N_sp     , &
             ratio_Fe2N_phaeo  , &
             ratio_Fe2DON      , & 
             ratio_Fe2DOC_s    , & 
             ratio_Fe2DOC_l    , & 
             fr_resp           , & 
             tau_min           , & 
             tau_max           , & 
             algal_vel         , & 
             R_dFe2dust        , & 
             dustFe_sol        , & 
             chlabs_diatoms    , &
             chlabs_sp         , &
             chlabs_phaeo      , &
             alpha2max_low_diatoms , & 
             alpha2max_low_sp      , & 
             alpha2max_low_phaeo   , & 
             beta2max_diatoms , &
             beta2max_sp      , & 
             beta2max_phaeo   , & 
             mu_max_diatoms   , &
             mu_max_sp        , & 
             mu_max_phaeo     , & 
             grow_Tdep_diatoms, &
             grow_Tdep_sp     , & 
             grow_Tdep_phaeo  , & 
             fr_graze_diatoms , &
             fr_graze_sp      , & 
             fr_graze_phaeo   , & 
             mort_pre_diatoms , &
             mort_pre_sp      , & 
             mort_pre_phaeo   , & 
             mort_Tdep_diatoms, &
             mort_Tdep_sp     , &  
             mort_Tdep_phaeo  , &  
             k_exude_diatoms  , &
             k_exude_sp       , &  
             k_exude_phaeo    , &  
             K_Nit_diatoms    , &
             K_Nit_sp         , &  
             K_Nit_phaeo      , &  
             K_Am_diatoms     , &
             K_Am_sp          , &   
             K_Am_phaeo       , &   
             K_Sil_diatoms    , &
             K_Sil_sp         , &   
             K_Sil_phaeo      , &   
             K_Fe_diatoms     , &
             K_Fe_sp          , &   
             K_Fe_phaeo       , &    
             f_don_protein    , &
             kn_bac_protein   , &
             f_don_Am_protein , &
             f_doc_s            , &
             f_doc_l            , &   
             f_exude_s          , &
             f_exude_l          , & 
             k_bac_s            , &
             k_bac_l            , & 
             T_max              , &
             fsal               , &
             op_dep_min         , &
             fr_graze_s         , &
             fr_graze_e         , &
             fr_mort2min        , &
             fr_dFe             , &
             k_nitrif           , &
             t_iron_conv        , &
             max_loss           , &
             max_dfe_doc1       , &
             fr_resp_s          , &
             y_sk_DMS           , &
             t_sk_conv          , &
             t_sk_ox            , & 
             algaltype_diatoms  , & 
             algaltype_sp       , &
             algaltype_phaeo    , &
             nitratetype        , &
             ammoniumtype       , &
             silicatetype       , &
             dmspptype          , &
             dmspdtype          , &
             humtype            , &
             doctype_s          , &
             doctype_l          , &
             dictype_1          , &
             dontype_protein    , & 
             fedtype_1          , & 
             feptype_1          , & 
             zaerotype_bc1      , & 
             zaerotype_bc2      , & 
             zaerotype_dust1    , & 
             zaerotype_dust2    , & 
             zaerotype_dust3    , & 
             zaerotype_dust4    , & 
             ratio_C2N_diatoms  , & 
             ratio_C2N_sp       , & 
             ratio_C2N_phaeo    , & 
             ratio_chl2N_diatoms, & 
             ratio_chl2N_sp     , & 
             ratio_chl2N_phaeo  , & 
             F_abs_chl_diatoms  , & 
             F_abs_chl_sp       , & 
             F_abs_chl_phaeo    , & 
             ratio_C2N_proteins , &
             snwredist, &
             use_smliq_pnd, &
             rsnw_fall, &
             rsnw_tmax, &
             rhosnew, &
             rhosmax, &
             windmin, &
             drhosdwind
            !restore_bgc

!-----------------------------------------------------------------------
! Parameters for thermodynamics
!-----------------------------------------------------------------------

        integer (kind=int_kind), intent(in) :: &
             ktherm_in          ! type of thermodynamics
                                ! 0 = 0-layer approximation
                                ! 1 = Bitz and Lipscomb 1999
                                ! 2 = mushy layer theory

        character (char_len), intent(in) :: &
             conduct_in, &      ! 'MU71' or 'bubbly'
             fbot_xfer_type_in  ! transfer coefficient type for ice-ocean heat flux
        
        logical (kind=log_kind), intent(in) :: &
             calc_Tsfc_in       ! if true, calculate surface temperature
                                ! if false, Tsfc is computed elsewhere and
                                ! atmos-ice fluxes are provided to CICE

        real (kind=real_kind), intent(in) :: &
             ustar_min_in       ! minimum friction velocity for ice-ocean heat flux
 
        ! mushy thermo
        real(kind=real_kind), intent(in) :: &
             a_rapid_mode_in      , & ! channel radius for rapid drainage mode (m)
             Rac_rapid_mode_in    , & ! critical Rayleigh number for rapid drainage mode
             aspect_rapid_mode_in , & ! aspect ratio for rapid drainage mode (larger=wider)
             dSdt_slow_mode_in    , & ! slow mode drainage strength (m s-1 K-1)
             phi_c_slow_mode_in   , & ! liquid fraction porosity cutoff for slow mode
             phi_i_mushy_in           ! liquid fraction of congelation ice
        
!-----------------------------------------------------------------------
! Parameters for radiation
!-----------------------------------------------------------------------

        character (len=char_len), intent(in) :: &
             shortwave_in, & ! shortwave method, 'default' ('ccsm3') or 'dEdd'
             albedo_type_in  ! albedo parameterization, 'default' ('ccsm3') or 'constant'
                             ! shortwave='dEdd' overrides this parameter

        ! baseline albedos for ccsm3 shortwave, set in namelist
        real (kind=real_kind), intent(in) :: &
             albicev_in  , & ! visible ice albedo for h > ahmax
             albicei_in  , & ! near-ir ice albedo for h > ahmax
             albsnowv_in , & ! cold snow albedo, visible
             albsnowi_in , & ! cold snow albedo, near IR
             ahmax_in        ! thickness above which ice albedo is constant (m)
        
        ! dEdd tuning parameters, set in namelist
        real (kind=real_kind), intent(in) :: &
             R_ice_in    , & ! sea ice tuning parameter; +1 > 1sig increase in albedo
             R_pnd_in    , & ! ponded ice tuning parameter; +1 > 1sig increase in albedo
             R_snw_in    , & ! snow tuning parameter; +1 > ~.01 change in broadband albedo
             dT_mlt_in   , & ! change in temp for non-melt to melt snow grain 
                             ! radius change (C)
             rsnw_mlt_in , & ! maximum melting snow grain radius (10^-6 m)
             kalg_in         ! algae absorption coefficient for 0.5 m thick layer

        ! snicar 5 band system, set in namelist
        logical (kind=log_kind), intent(in) :: &
             use_snicar_in ! if true, use 5-band snicar IOPs for
                              ! shortwave radiative calculation of
                              ! snow-coverd sea ice

!-----------------------------------------------------------------------
! Parameters for ridging and strength
!-----------------------------------------------------------------------

        integer (kind=int_kind), intent(in) :: & ! defined in namelist 
             kstrength_in  , & ! 0 for simple Hibler (1979) formulation 
                               ! 1 for Rothrock (1975) pressure formulation 
             krdg_partic_in, & ! 0 for Thorndike et al. (1975) formulation 
                               ! 1 for exponential participation function 
             krdg_redist_in    ! 0 for Hibler (1980) formulation 
                               ! 1 for exponential redistribution function 
 
        real (kind=real_kind), intent(in) :: &  
             mu_rdg_in, &      ! gives e-folding scale of ridged ice (m^.5) 
                               ! (krdg_redist = 1) 
             Cf_in             ! ratio of ridging work to PE change in ridging (kstrength = 1)

!-----------------------------------------------------------------------
! Parameters for atmosphere
!-----------------------------------------------------------------------

        character (len=char_len), intent(in) :: &
             atmbndy_in ! atmo boundary method, 'default' ('ccsm3') or 'constant'
        
        logical (kind=log_kind), intent(in) :: &
             calc_strair_in, &  ! if true, calculate wind stress components
             formdrag_in,    &  ! if true, calculate form drag
             highfreq_in        ! if true, use high frequency coupling
        
        integer (kind=int_kind), intent(in) :: &
             natmiter_in        ! number of iterations for boundary layer calculations
        
!-----------------------------------------------------------------------
! Parameters for ocean
!-----------------------------------------------------------------------

        logical (kind=log_kind), intent(in) :: &
             oceanmixed_ice_in           ! if true, use ocean mixed layer
        
        character(len=char_len), intent(in) :: &
             tfrz_option_in              ! form of ocean freezing temperature
                                         ! 'minus1p8' = -1.8 C
                                         ! 'linear_salt' = -depressT * sss
                                         ! 'mushy' conforms with ktherm=2

!-----------------------------------------------------------------------
! Parameters for the ice thickness distribution
!-----------------------------------------------------------------------

        integer (kind=int_kind), intent(in) :: &
             kitd_in        , & ! type of itd conversions
                                !   0 = delta function
                                !   1 = linear remap
             kcatbound_in       !   0 = old category boundary formula
                                !   1 = new formula giving round numbers
                                !   2 = WMO standard
                                !   3 = asymptotic formula

!-----------------------------------------------------------------------
! Parameters for biogeochemistry
!-----------------------------------------------------------------------

   !  character(char_len_long), intent(in) :: & 
   !     bgc_data_dir_in   ! directory for biogeochemistry data

     character(char_len), intent(in) :: &     
        bgc_flux_type_in    ! type of ocean-ice piston velocity 
                            ! 'constant', 'Jin2006'      
    !    sil_data_type_in  , & ! 'default', 'clim'
    !    nit_data_type_in  , & ! 'default', 'clim'   
    !    fe_data_type_in   , & ! 'default', 'clim'      

      logical (kind=log_kind), intent(in) :: &
         z_tracers_in,      & ! if .true., bgc or aerosol tracers are vertically resolved
         scale_bgc_in,      & ! if .true., initialize bgc tracers proportionally with salinity
         solve_zbgc_in,     & ! if .true., solve vertical biochemistry portion of code
         dEdd_algae_in,     & ! if .true., algal absorptionof Shortwave is computed in the
         modal_aero_in        ! if .true., use modal aerosol formulation in shortwave
        
      logical (kind=log_kind), intent(in) :: & 
         skl_bgc_in,        &   ! if true, solve skeletal biochemistry
         solve_zsal_in          ! if true, update salinity profile from solve_S_dt

      real (kind=real_kind), intent(in) :: & 
         grid_o_in      , & ! for bottom flux        
         l_sk_in        , & ! characteristic diffusive scale (zsalinity) (m)
         grid_o_t_in    , & ! top grid point length scale 
         initbio_frac_in, & ! fraction of ocean tracer concentration used to initialize tracer 
         frazil_scav_in , & ! multiple of ocean tracer concentration due to frazil scavenging
         phi_snow_in        ! snow porosity at the ice/snow interface 

      real (kind=real_kind), intent(in) :: & 
         grid_oS_in     , & ! for bottom flux (zsalinity)
         l_skS_in           ! 0.02 characteristic skeletal layer thickness (m) (zsalinity)
      real (kind=real_kind), intent(in) :: &
         ratio_Si2N_diatoms_in, &   ! algal Si to N (mol/mol)
         ratio_Si2N_sp_in     , &
         ratio_Si2N_phaeo_in  , &
         ratio_S2N_diatoms_in , &   ! algal S  to N (mol/mol)
         ratio_S2N_sp_in      , &
         ratio_S2N_phaeo_in   , &
         ratio_Fe2C_diatoms_in, &   ! algal Fe to C  (umol/mol)
         ratio_Fe2C_sp_in     , &
         ratio_Fe2C_phaeo_in  , &
         ratio_Fe2N_diatoms_in, &   ! algal Fe to N  (umol/mol)
         ratio_Fe2N_sp_in     , &
         ratio_Fe2N_phaeo_in  , &
         ratio_Fe2DON_in      , &   ! Fe to N of DON (nmol/umol)
         ratio_Fe2DOC_s_in    , &   ! Fe to C of DOC (nmol/umol) saccharids
         ratio_Fe2DOC_l_in    , &   ! Fe to C of DOC (nmol/umol) lipids
         fr_resp_in           , &   ! fraction of algal growth lost due to respiration
         tau_min_in           , &   ! rapid mobile to stationary exchanges (s) = 1.5 hours
         tau_max_in           , &   ! long time mobile to stationary exchanges (s) = 2 days
         algal_vel_in         , &   ! 0.5 cm/d(m/s) Lavoie 2005  1.5 cm/day
         R_dFe2dust_in        , &   !  g/g (3.5% content) Tagliabue 2009
         dustFe_sol_in        , &   ! solubility fraction
         chlabs_diatoms_in   , & ! chl absorption (1/m/(mg/m^3))
         chlabs_sp_in        , & !
         chlabs_phaeo_in     , & !
         alpha2max_low_diatoms_in , & ! light limitation (1/(W/m^2))  
         alpha2max_low_sp_in      , & 
         alpha2max_low_phaeo_in   , & 
         beta2max_diatoms_in , & ! light inhibition (1/(W/m^2))  
         beta2max_sp_in      , & 
         beta2max_phaeo_in   , & 
         mu_max_diatoms_in   , & ! maximum growth rate (1/day)       
         mu_max_sp_in        , & 
         mu_max_phaeo_in     , & 
         grow_Tdep_diatoms_in, & ! Temperature dependence of growth (1/C)
         grow_Tdep_sp_in     , & 
         grow_Tdep_phaeo_in  , & 
         fr_graze_diatoms_in , & ! Fraction grazed
         fr_graze_sp_in      , & 
         fr_graze_phaeo_in   , & 
         mort_pre_diatoms_in , & ! Mortality (1/day)
         mort_pre_sp_in      , & 
         mort_pre_phaeo_in   , & 
         mort_Tdep_diatoms_in, & ! T dependence of mortality (1/C) 
         mort_Tdep_sp_in     , &  
         mort_Tdep_phaeo_in  , &  
         k_exude_diatoms_in  , & ! algal exudation (1/d)
         k_exude_sp_in       , &  
         k_exude_phaeo_in    , &  
         K_Nit_diatoms_in    , & ! nitrate half saturation (mmol/m^3)
         K_Nit_sp_in         , &  
         K_Nit_phaeo_in      , &  
         K_Am_diatoms_in     , & ! ammonium half saturation (mmol/m^3)
         K_Am_sp_in          , &   
         K_Am_phaeo_in       , &   
         K_Sil_diatoms_in    , & ! silicate half saturation (mmol/m^3)
         K_Sil_sp_in         , &   
         K_Sil_phaeo_in      , &   
         K_Fe_diatoms_in     , & ! iron half saturation (nM)
         K_Fe_sp_in          , &   
         K_Fe_phaeo_in       , &    
         f_don_protein_in    , & ! fraction of spilled grazing to proteins            
         kn_bac_protein_in   , & ! Bacterial degredation of DON (1/d)                  
         f_don_Am_protein_in , & ! fraction of remineralized DON to ammonium          
         f_doc_s_in          , & ! fraction of mortality to DOC 
         f_doc_l_in          , &   
         f_exude_s_in        , & ! fraction of exudation to DOC
         f_exude_l_in        , & 
         k_bac_s_in          , & ! Bacterial degredation of DOC (1/d)
         k_bac_l_in          , & 
         T_max_in            , & ! maximum temperature (C)
         fsal_in             , & ! Salinity limitation (ppt)
         op_dep_min_in       , & ! Light attenuates for optical depths exceeding min
         fr_graze_s_in       , & ! fraction of grazing spilled or slopped
         fr_graze_e_in       , & ! fraction of assimilation excreted 
         fr_mort2min_in      , & ! fractionation of mortality to Am
         fr_dFe_in           , & ! fraction of remineralized nitrogen 
                                    ! (in units of algal iron)
         k_nitrif_in         , & ! nitrification rate (1/day)            
         t_iron_conv_in      , & ! desorption loss pFe to dFe (day)
         max_loss_in         , & ! restrict uptake to % of remaining value 
         max_dfe_doc1_in     , & ! max ratio of dFe to saccharides in the ice 
                                    ! (nM Fe/muM C)    
         fr_resp_s_in        , & ! DMSPd fraction of respiration loss as DMSPd
         y_sk_DMS_in         , & ! fraction conversion given high yield
         t_sk_conv_in        , & ! Stefels conversion time (d)
         t_sk_ox_in          , &   ! DMS oxidation time (d)
         algaltype_diatoms_in  , & ! mobility type
         algaltype_sp_in       , & !
         algaltype_phaeo_in    , & !
         nitratetype_in        , & !
         ammoniumtype_in       , & !
         silicatetype_in       , & !
         dmspptype_in          , & !
         dmspdtype_in          , & !
         humtype_in            , & !
         doctype_s_in          , & !
         doctype_l_in          , & !
         dictype_1_in          , & !
         dontype_protein_in    , & !
         fedtype_1_in          , & !
         feptype_1_in          , & !
         zaerotype_bc1_in      , & !
         zaerotype_bc2_in      , & !
         zaerotype_dust1_in    , & !
         zaerotype_dust2_in    , & !
         zaerotype_dust3_in    , & !
         zaerotype_dust4_in    , & !
         ratio_C2N_diatoms_in  , & ! algal C to N ratio (mol/mol)
         ratio_C2N_sp_in       , & !
         ratio_C2N_phaeo_in    , & !
         ratio_chl2N_diatoms_in, & ! algal chlorophyll to N ratio (mg/mmol)
         ratio_chl2N_sp_in     , & !
         ratio_chl2N_phaeo_in  , & !
         F_abs_chl_diatoms_in  , & ! scales absorbed radiation for dEdd
         F_abs_chl_sp_in       , & !
         F_abs_chl_phaeo_in    , & !
         ratio_C2N_proteins_in     ! ratio of C to N in proteins (mol/mol)       

     !logical (kind=log_kind), intent(in) :: & 
     !   restore_bgc_in      ! if true, restore nitrate

!-----------------------------------------------------------------------
! Parameters for melt ponds
!-----------------------------------------------------------------------

        real (kind=real_kind), intent(in) :: &
             hs0_in             ! snow depth for transition to bare sea ice (m)
        
        ! level-ice ponds
        character (len=char_len), intent(in) :: &
             frzpnd_in          ! pond refreezing parameterization
        
        real (kind=real_kind), intent(in) :: &
             dpscale_in, &      ! alter e-folding time scale for flushing 
             rfracmin_in, &     ! minimum retained fraction of meltwater
             rfracmax_in, &     ! maximum retained fraction of meltwater
             pndaspect_in, &    ! ratio of pond depth to pond fraction
             hs1_in             ! tapering parameter for snow on pond ice
        
        ! topo ponds
        real (kind=real_kind), intent(in) :: &
             hp1_in             ! critical parameter for pond ice thickness
        
!-----------------------------------------------------------------------
! Parameters for snow
!-----------------------------------------------------------------------

      ! snow metamorphism parameters, set in namelist
      real (kind=real_kind), intent(in) :: &
         rsnw_fall_in , & ! fallen snow grain radius (10^-6 m))  54.5 um CLM **
                       ! 30 um is minimum for defined mie properties 
         rsnw_tmax_in , & ! maximum dry metamorphism snow grain radius (10^-6 m)
                       ! 1500 um is maximum for defined mie properties
         rhosnew_in   , & ! new snow density (kg/m^3)
         rhosmax_in   , & ! maximum snow density (kg/m^3)
         windmin_in   , & ! minimum wind speed to compact snow (m/s)
         drhosdwind_in    ! wind compaction factor (kg s/m^4)

      character(len=char_len), intent(in) :: & 
         snwredist_in     ! type of snow redistribution
                       ! '30percent' = 30% rule, precip only
                       ! '30percentsw' = 30% rule with shortwave
                       ! 'ITDsd' = Lecomte PhD, 2014
                       ! 'ITDrdg' = like ITDsd but use level/ridged ice
                       ! 'default' or 'none' = none

      logical (kind=log_kind), intent(in) :: &
         use_smliq_pnd_in ! if true, use snow liquid tracer for ponds

        ktherm = ktherm_in
        conduct = conduct_in
        fbot_xfer_type = fbot_xfer_type_in
        calc_Tsfc = calc_Tsfc_in
        ustar_min = ustar_min_in
        a_rapid_mode = a_rapid_mode_in
        Rac_rapid_mode = Rac_rapid_mode_in
        aspect_rapid_mode = aspect_rapid_mode_in
        dSdt_slow_mode = dSdt_slow_mode_in
        phi_c_slow_mode = phi_c_slow_mode_in
        phi_i_mushy = phi_i_mushy_in
        shortwave = shortwave_in
        use_snicar = use_snicar_in
        albedo_type = albedo_type_in
        albicev = albicev_in
        albicei = albicei_in
        albsnowv = albsnowv_in
        albsnowi = albsnowi_in
        ahmax = ahmax_in
        R_ice = R_ice_in
        R_pnd = R_pnd_in
        R_snw = R_snw_in
        dT_mlt = dT_mlt_in
        rsnw_mlt = rsnw_mlt_in
        kalg = kalg_in
        kstrength = kstrength_in
        krdg_partic = krdg_partic_in
        krdg_redist = krdg_redist_in
        mu_rdg = mu_rdg_in
        Cf = Cf_in
        atmbndy = atmbndy_in
        calc_strair = calc_strair_in
        formdrag = formdrag_in
        highfreq = highfreq_in
        natmiter = natmiter_in
        oceanmixed_ice = oceanmixed_ice_in
        tfrz_option = tfrz_option_in
        kitd = kitd_in
        kcatbound = kcatbound_in
        hs0 = hs0_in
        frzpnd = frzpnd_in
        dpscale = dpscale_in
        rfracmin = rfracmin_in
        rfracmax = rfracmax_in
        pndaspect = pndaspect_in
        hs1 = hs1_in
        hp1 = hp1_in
     !  bgc_data_dir = bgc_data_dir_in
     !  sil_data_type= sil_data_type_in
     !  nit_data_type = nit_data_type_in
     !  fe_data_type = fe_data_type_in
        bgc_flux_type = bgc_flux_type_in
        z_tracers = z_tracers_in
        scale_bgc = scale_bgc_in
        solve_zbgc = solve_zbgc_in
        dEdd_algae = dEdd_algae_in
        skl_bgc = skl_bgc_in
        grid_o = grid_o_in
        l_sk = l_sk_in
        grid_o_t = grid_o_t_in
        initbio_frac = initbio_frac_in
        frazil_scav = frazil_scav_in
        grid_oS = grid_oS_in
        l_skS = l_skS_in
        phi_snow = phi_snow_in
     !  restore_bgc = restore_bgc_in
        ratio_Si2N_diatoms= ratio_Si2N_diatoms_in 
        ratio_Si2N_sp     = ratio_Si2N_sp_in
        ratio_Si2N_phaeo  = ratio_Si2N_phaeo_in
        ratio_S2N_diatoms = ratio_S2N_diatoms_in
        ratio_S2N_sp      = ratio_S2N_sp_in
        ratio_S2N_phaeo   = ratio_S2N_phaeo_in
        ratio_Fe2C_diatoms= ratio_Fe2C_diatoms_in 
        ratio_Fe2C_sp     = ratio_Fe2C_sp_in
        ratio_Fe2C_phaeo  = ratio_Fe2C_phaeo_in
        ratio_Fe2N_diatoms= ratio_Fe2N_diatoms_in 
        ratio_Fe2N_sp     = ratio_Fe2N_sp_in
        ratio_Fe2N_phaeo  = ratio_Fe2N_phaeo_in
        ratio_Fe2DON      = ratio_Fe2DON_in
        ratio_Fe2DOC_s    = ratio_Fe2DOC_s_in
        ratio_Fe2DOC_l    = ratio_Fe2DOC_l_in
        fr_resp           = fr_resp_in
        tau_min           = tau_min_in
        tau_max           = tau_max_in
        algal_vel         = algal_vel_in
        R_dFe2dust        = R_dFe2dust_in
        dustFe_sol        = dustFe_sol_in
        chlabs_diatoms    = chlabs_diatoms_in
        chlabs_sp         = chlabs_sp_in
        chlabs_phaeo      = chlabs_phaeo_in
        alpha2max_low_diatoms = alpha2max_low_diatoms_in
        alpha2max_low_sp      = alpha2max_low_sp_in
        alpha2max_low_phaeo   = alpha2max_low_phaeo_in
        beta2max_diatoms = beta2max_diatoms_in
        beta2max_sp      = beta2max_sp_in
        beta2max_phaeo   = beta2max_phaeo_in
        mu_max_diatoms   = mu_max_diatoms_in
        mu_max_sp        = mu_max_sp_in
        mu_max_phaeo     = mu_max_phaeo_in
        grow_Tdep_diatoms= grow_Tdep_diatoms_in
        grow_Tdep_sp     = grow_Tdep_sp_in
        grow_Tdep_phaeo  = grow_Tdep_phaeo_in
        fr_graze_diatoms = fr_graze_diatoms_in
        fr_graze_sp      = fr_graze_sp_in
        fr_graze_phaeo   = fr_graze_phaeo_in
        mort_pre_diatoms = mort_pre_diatoms_in
        mort_pre_sp      = mort_pre_sp_in
        mort_pre_phaeo   = mort_pre_phaeo_in
        mort_Tdep_diatoms= mort_Tdep_diatoms_in
        mort_Tdep_sp     = mort_Tdep_sp_in
        mort_Tdep_phaeo  = mort_Tdep_phaeo_in
        k_exude_diatoms  = k_exude_diatoms_in
        k_exude_sp       = k_exude_sp_in
        k_exude_phaeo    = k_exude_phaeo_in
        K_Nit_diatoms    = K_Nit_diatoms_in
        K_Nit_sp         = K_Nit_sp_in
        K_Nit_phaeo      = K_Nit_phaeo_in
        K_Am_diatoms     = K_Am_diatoms_in
        K_Am_sp          = K_Am_sp_in
        K_Am_phaeo       = K_Am_phaeo_in
        K_Sil_diatoms    = K_Sil_diatoms_in
        K_Sil_sp         = K_Sil_sp_in
        K_Sil_phaeo      = K_Sil_phaeo_in
        K_Fe_diatoms     = K_Fe_diatoms_in
        K_Fe_sp          = K_Fe_sp_in
        K_Fe_phaeo       = K_Fe_phaeo_in
        f_don_protein    = f_don_protein_in
        kn_bac_protein   = kn_bac_protein_in
        f_don_Am_protein = f_don_Am_protein_in
        f_doc_s          = f_doc_s_in
        f_doc_l          = f_doc_l_in
        f_exude_s        = f_exude_s_in
        f_exude_l        = f_exude_l_in
        k_bac_s          = k_bac_s_in
        k_bac_l          = k_bac_l_in
        T_max            = T_max_in
        fsal             = fsal_in
        op_dep_min       = op_dep_min_in
        fr_graze_s       = fr_graze_s_in
        fr_graze_e       = fr_graze_e_in
        fr_mort2min      = fr_mort2min_in
        fr_dFe           = fr_dFe_in
        k_nitrif         = k_nitrif_in
        t_iron_conv      = t_iron_conv_in
        max_loss         = max_loss_in
        max_dfe_doc1     = max_dfe_doc1_in
        fr_resp_s        = fr_resp_s_in
        y_sk_DMS         = y_sk_DMS_in
        t_sk_conv        = t_sk_conv_in
        t_sk_ox          = t_sk_ox_in
        algaltype_diatoms  = algaltype_diatoms_in
        algaltype_sp       = algaltype_sp_in
        algaltype_phaeo    = algaltype_phaeo_in
        nitratetype        = nitratetype_in
        ammoniumtype       = ammoniumtype_in
        silicatetype       = silicatetype_in
        dmspptype          = dmspptype_in
        dmspdtype          = dmspdtype_in
        humtype            = humtype_in
        doctype_s          = doctype_s_in
        doctype_l          = doctype_l_in
        dictype_1          = dictype_1_in
        dontype_protein    = dontype_protein_in
        fedtype_1          = fedtype_1_in
        feptype_1          = feptype_1_in
        zaerotype_bc1      = zaerotype_bc1_in
        zaerotype_bc2      = zaerotype_bc2_in
        zaerotype_dust1    = zaerotype_dust1_in
        zaerotype_dust2    = zaerotype_dust2_in
        zaerotype_dust3    = zaerotype_dust3_in
        zaerotype_dust4    = zaerotype_dust4_in
        ratio_C2N_diatoms  = ratio_C2N_diatoms_in
        ratio_C2N_sp       = ratio_C2N_sp_in
        ratio_C2N_phaeo    = ratio_C2N_phaeo_in
        ratio_chl2N_diatoms= ratio_chl2N_diatoms_in
        ratio_chl2N_sp     = ratio_chl2N_sp_in
        ratio_chl2N_phaeo  = ratio_chl2N_phaeo_in
        F_abs_chl_diatoms  = F_abs_chl_diatoms_in
        F_abs_chl_sp       = F_abs_chl_sp_in
        F_abs_chl_phaeo    = F_abs_chl_phaeo_in
        ratio_C2N_proteins = ratio_C2N_proteins_in
        snwredist = snwredist_in
        use_smliq_pnd = use_smliq_pnd_in
        rsnw_fall = rsnw_fall_in
        rsnw_tmax = rsnw_tmax_in
        rhosnew = rhosnew_in
        rhosmax = rhosmax_in
        windmin = windmin_in
        drhosdwind = drhosdwind_in

      end subroutine colpkg_init_parameters

!=======================================================================
! set tracer active flags

      subroutine colpkg_init_tracer_flags(&
           tr_iage_in      , & ! if .true., use age tracer
           tr_FY_in        , & ! if .true., use first-year area tracer
           tr_lvl_in       , & ! if .true., use level ice tracer
           tr_pond_in      , & ! if .true., use melt pond tracer
           tr_pond_cesm_in , & ! if .true., use cesm pond tracer
           tr_pond_lvl_in  , & ! if .true., use level-ice pond tracer
           tr_pond_topo_in , & ! if .true., use explicit topography-based ponds
           tr_snow_in      , & ! if .true., use snow trcrs (smice, smliq, rhos_cmp)
           tr_rsnw_in      , & ! if .true., use snow grain radius tracer
           tr_aero_in      , & ! if .true., use aerosol tracers
           tr_brine_in     , & ! if .true., brine height differs from ice thickness
           tr_bgc_S_in     , & ! if .true., use zsalinity
           tr_zaero_in     , & ! if .true., black carbon is tracers  (n_zaero)
           tr_bgc_Nit_in   , & ! if .true., Nitrate tracer in ice 
           tr_bgc_N_in     , & ! if .true., algal nitrogen tracers  (n_algae)
           tr_bgc_DON_in   , & ! if .true., DON pools are tracers  (n_don)
           tr_bgc_C_in     , & ! if .true., algal carbon tracers + DOC and DIC 
           tr_bgc_chl_in   , & ! if .true., algal chlorophyll tracers 
           tr_bgc_Am_in    , & ! if .true., ammonia/um as nutrient tracer 
           tr_bgc_Sil_in   , & ! if .true., silicon as nutrient tracer 
           tr_bgc_DMS_in   , & ! if .true., DMS as product tracer 
           tr_bgc_Fe_in    , & ! if .true., Fe as product tracer 
           tr_bgc_hum_in   , & ! if .true., hum as tracer 
           tr_bgc_PON_in)      ! if .true., PON as product tracer 


        use ice_colpkg_tracers, only: &
             tr_iage      , & ! if .true., use age tracer
             tr_FY        , & ! if .true., use first-year area tracer
             tr_lvl       , & ! if .true., use level ice tracer
             tr_pond      , & ! if .true., use melt pond tracer
             tr_pond_cesm , & ! if .true., use cesm pond tracer
             tr_pond_lvl  , & ! if .true., use level-ice pond tracer
             tr_pond_topo , & ! if .true., use explicit topography-based ponds
             tr_snow      , & ! if .true., use snow trcrs (smice, smliq, rhos_cmp)
             tr_rsnw      , & ! if .true., use snow grain radius tracer
             tr_aero      , & ! if .true., use aerosol tracers
             tr_brine     , & ! if .true., brine height differs from ice thickness
             tr_bgc_S     , & ! if .true., use zsalinity
             tr_zaero     , & ! if .true., black carbon is tracers  (n_zaero)
             tr_bgc_Nit   , & ! if .true., Nitrate tracer in ice 
             tr_bgc_N     , & ! if .true., algal nitrogen tracers  (n_algae)
             tr_bgc_DON   , & ! if .true., DON pools are tracers  (n_don)
             tr_bgc_C     , & ! if .true., algal carbon tracers + DOC and DIC 
             tr_bgc_chl   , & ! if .true., algal chlorophyll tracers 
             tr_bgc_Am    , & ! if .true., ammonia/um as nutrient tracer 
             tr_bgc_Sil   , & ! if .true., silicon as nutrient tracer 
             tr_bgc_DMS   , & ! if .true., DMS as product tracer 
             tr_bgc_Fe    , & ! if .true., Fe as product tracer 
             tr_bgc_hum   , & ! if .true., hum as product tracer 
             tr_bgc_PON       ! if .true., PON as product tracer 


        logical, intent(in) :: &
             tr_iage_in      , & ! if .true., use age tracer
             tr_FY_in        , & ! if .true., use first-year area tracer
             tr_lvl_in       , & ! if .true., use level ice tracer
             tr_pond_in      , & ! if .true., use melt pond tracer
             tr_pond_cesm_in , & ! if .true., use cesm pond tracer
             tr_pond_lvl_in  , & ! if .true., use level-ice pond tracer
             tr_pond_topo_in , & ! if .true., use explicit topography-based ponds
             tr_snow_in      , & ! if .true., use snow trcrs (smice, smliq, rhos_cmp)
             tr_rsnw_in      , & ! if .true., use snow grain radius tracer
             tr_aero_in      , & ! if .true., use aerosol tracers
             tr_brine_in     , & ! if .true., brine height differs from ice thickness
             tr_bgc_S_in     , & ! if .true., use zsalinity
             tr_zaero_in     , & ! if .true., black carbon is tracers  (n_zaero)
             tr_bgc_Nit_in   , & ! if .true., Nitrate tracer in ice 
             tr_bgc_N_in     , & ! if .true., algal nitrogen tracers  (n_algae)
             tr_bgc_DON_in   , & ! if .true., DON pools are tracers  (n_don)
             tr_bgc_C_in     , & ! if .true., algal carbon tracers + DOC and DIC 
             tr_bgc_chl_in   , & ! if .true., algal chlorophyll tracers 
             tr_bgc_Am_in    , & ! if .true., ammonia/um as nutrient tracer 
             tr_bgc_Sil_in   , & ! if .true., silicon as nutrient tracer 
             tr_bgc_DMS_in   , & ! if .true., DMS as product tracer 
             tr_bgc_Fe_in    , & ! if .true., Fe as product tracer 
             tr_bgc_hum_in   , & ! if .true., hum as product tracer 
             tr_bgc_PON_in       ! if .true., PON as product tracer 

        tr_iage      = tr_iage_in
        tr_FY        = tr_FY_in
        tr_lvl       = tr_lvl_in
        tr_pond      = tr_pond_in
        tr_pond_cesm = tr_pond_cesm_in
        tr_pond_lvl  = tr_pond_lvl_in
        tr_pond_topo = tr_pond_topo_in
        tr_snow      = tr_snow_in
        tr_rsnw      = tr_rsnw_in
        tr_aero      = tr_aero_in
        tr_brine     = tr_brine_in
        tr_bgc_S     = tr_bgc_S_in
        tr_zaero     = tr_zaero_in 
        tr_bgc_Nit   = tr_bgc_Nit_in
        tr_bgc_N     = tr_bgc_N_in 
        tr_bgc_DON   = tr_bgc_DON_in
        tr_bgc_C     = tr_bgc_C_in 
        tr_bgc_chl   = tr_bgc_chl_in
        tr_bgc_Am    = tr_bgc_Am_in
        tr_bgc_Sil   = tr_bgc_Sil_in
        tr_bgc_DMS   = tr_bgc_DMS_in
        tr_bgc_Fe    = tr_bgc_Fe_in 
        tr_bgc_hum   = tr_bgc_hum_in
        tr_bgc_PON   = tr_bgc_PON_in 

      end subroutine colpkg_init_tracer_flags

!=======================================================================

      subroutine colpkg_init_tracer_indices(&
           nt_Tsfc_in, & ! ice/snow temperature
           nt_qice_in, & ! volume-weighted ice enthalpy (in layers)
           nt_qsno_in, & ! volume-weighted snow enthalpy (in layers)
           nt_sice_in, & ! volume-weighted ice bulk salinity (CICE grid layers)
           nt_fbri_in, & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
           nt_iage_in, & ! volume-weighted ice age
           nt_FY_in, & ! area-weighted first-year ice area
           nt_alvl_in, & ! level ice area fraction
           nt_vlvl_in, & ! level ice volume fraction
           nt_apnd_in, & ! melt pond area fraction
           nt_hpnd_in, & ! melt pond depth
           nt_ipnd_in, & ! melt pond refrozen lid thickness
           nt_aero_in, & ! starting index for aerosols in ice 
           nt_smice_in, & ! snow ice mass
           nt_smliq_in, & ! snow liquid mass
           nt_rsnw_in, & ! snow grain radius
           nt_rhos_in, & ! snow density
           nt_zaero_in,   & !  black carbon and other aerosols
           nt_bgc_N_in ,  & ! diatoms, phaeocystis, pico/small   
           nt_bgc_C_in ,  & ! diatoms, phaeocystis, pico/small   
           nt_bgc_chl_in, & ! diatoms, phaeocystis, pico/small 
           nt_bgc_DOC_in, & !  dissolved organic carbon
           nt_bgc_DON_in, & !  dissolved organic nitrogen
           nt_bgc_DIC_in, & !  dissolved inorganic carbon
           nt_bgc_Fed_in, & !  dissolved iron
           nt_bgc_Fep_in, & !  particulate iron
           nt_bgc_Nit_in, & ! nutrients  
           nt_bgc_Am_in,  & ! 
           nt_bgc_Sil_in, & !
           nt_bgc_DMSPp_in,&! trace gases (skeletal layer)
           nt_bgc_DMSPd_in,&! 
           nt_bgc_DMS_in, & ! 
           nt_bgc_hum_in, & ! 
           nt_bgc_PON_in, & ! zooplankton and detritus  
           nlt_zaero_in,  & !  black carbon and other aerosols
           nlt_bgc_N_in , & ! diatoms, phaeocystis, pico/small   
           nlt_bgc_C_in , & ! diatoms, phaeocystis, pico/small   
           nlt_bgc_chl_in,& ! diatoms, phaeocystis, pico/small 
           nlt_bgc_DOC_in,& !  dissolved organic carbon
           nlt_bgc_DON_in,& !  dissolved organic nitrogen
           nlt_bgc_DIC_in,& !  dissolved inorganic carbon
           nlt_bgc_Fed_in,& !  dissolved iron
           nlt_bgc_Fep_in,& !  particulate iron
           nlt_bgc_Nit_in,& ! nutrients  
           nlt_bgc_Am_in, & ! 
           nlt_bgc_Sil_in,& !
           nlt_bgc_DMSPp_in,&! trace gases (skeletal layer)
           nlt_bgc_DMSPd_in,&! 
           nlt_bgc_DMS_in,& ! 
           nlt_bgc_hum_in,& ! 
           nlt_bgc_PON_in,& ! zooplankton and detritus  
           nt_zbgc_frac_in,&! fraction of tracer in the mobile phase
           nt_bgc_S_in,   & ! Bulk salinity in fraction ice with dynamic salinity (Bio grid))
           nlt_chl_sw_in, & ! points to total chla in trcrn_sw
           nlt_zaero_sw_in,&! black carbon and dust in trcrn_sw
                            ! Index Dimensions: 
           n_algae, n_algalC, & !
           n_algalchl, n_DOC, & !
           n_DON,n_DIC,n_dFe, & !
           n_pFe, n_aerosols, & !
           bio_index_o_in,    & ! nlt index to fixed data array
           bio_index_in,      & ! nlt index to nt index
           nbtrcr)

        use ice_colpkg_tracers, only: &
             nt_Tsfc, & ! ice/snow temperature
             nt_qice, & ! volume-weighted ice enthalpy (in layers)
             nt_qsno, & ! volume-weighted snow enthalpy (in layers)
             nt_sice, & ! volume-weighted ice bulk salinity (CICE grid layers)
             nt_fbri, & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
             nt_iage, & ! volume-weighted ice age
             nt_FY, & ! area-weighted first-year ice area
             nt_alvl, & ! level ice area fraction
             nt_vlvl, & ! level ice volume fraction
             nt_apnd, & ! melt pond area fraction
             nt_hpnd, & ! melt pond depth
             nt_ipnd, & ! melt pond refrozen lid thickness
             nt_aero, & ! starting index for aerosols in ice
             nt_smice, & ! snow ice mass
             nt_smliq, & ! snow liquid mass
             nt_rsnw, & ! snow grain radius
             nt_rhos, & ! snow density
             nt_zaero,   & !  black carbon and other aerosols
             nt_bgc_N ,  & ! diatoms, phaeocystis, pico/small   
             nt_bgc_C ,  & ! diatoms, phaeocystis, pico/small   
             nt_bgc_chl, & ! diatoms, phaeocystis, pico/small 
             nt_bgc_DOC, & !  dissolved organic carbon
             nt_bgc_DON, & !  dissolved organic nitrogen
             nt_bgc_DIC, & !  dissolved inorganic carbon
             nt_bgc_Fed, & !  dissolved iron
             nt_bgc_Fep, & !  particulate iron
             nt_bgc_Nit, & ! nutrients  
             nt_bgc_Am,  & ! 
             nt_bgc_Sil, & !
             nt_bgc_DMSPp,&! trace gases (skeletal layer)
             nt_bgc_DMSPd,&! 
             nt_bgc_DMS, & ! 
             nt_bgc_hum, & ! 
             nt_bgc_PON, & ! zooplankton and detritus 
             nlt_zaero,  & !  black carbon and other aerosols
             nlt_bgc_N , & ! diatoms, phaeocystis, pico/small   
             nlt_bgc_C , & ! diatoms, phaeocystis, pico/small   
             nlt_bgc_chl,& ! diatoms, phaeocystis, pico/small 
             nlt_bgc_DOC,& !  dissolved organic carbon
             nlt_bgc_DON,& !  dissolved organic nitrogen
             nlt_bgc_DIC,& !  dissolved inorganic carbon
             nlt_bgc_Fed,& !  dissolved iron
             nlt_bgc_Fep,& !  particulate iron
             nlt_bgc_Nit,& ! nutrients  
             nlt_bgc_Am, & ! 
             nlt_bgc_Sil,& !
             nlt_bgc_DMSPp,&! trace gases (skeletal layer)
             nlt_bgc_DMSPd,&! 
             nlt_bgc_DMS,& ! 
             nlt_bgc_hum,& ! 
             nlt_bgc_PON,& ! zooplankton and detritus   
             nt_zbgc_frac,&! fraction of tracer in the mobile phase
             nt_bgc_S,   & ! Bulk salinity in fraction ice with dynamic salinity (Bio grid))
             nlt_chl_sw, & ! points to total chla in trcrn_sw
             nlt_zaero_sw,&! black carbon and dust in trcrn_sw
             bio_index_o,& !
             bio_index 
        
        integer, intent(in) :: &
             nt_Tsfc_in, & ! ice/snow temperature
             nt_qice_in, & ! volume-weighted ice enthalpy (in layers)
             nt_qsno_in, & ! volume-weighted snow enthalpy (in layers)
             nt_sice_in, & ! volume-weighted ice bulk salinity (CICE grid layers)
             nt_fbri_in, & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
             nt_iage_in, & ! volume-weighted ice age
             nt_FY_in, & ! area-weighted first-year ice area
             nt_alvl_in, & ! level ice area fraction
             nt_vlvl_in, & ! level ice volume fraction
             nt_apnd_in, & ! melt pond area fraction
             nt_hpnd_in, & ! melt pond depth
             nt_ipnd_in, & ! melt pond refrozen lid thickness
             nt_aero_in, & ! starting index for aerosols in ice
             nt_smice_in, & ! snow ice mass
             nt_smliq_in, & ! snow liquid mass
             nt_rsnw_in, & ! snow grain radius
             nt_rhos_in, & ! snow density
             nt_bgc_Nit_in, & ! nutrients  
             nt_bgc_Am_in,  & ! 
             nt_bgc_Sil_in, & !
             nt_bgc_DMSPp_in,&! trace gases (skeletal layer)
             nt_bgc_DMSPd_in,&! 
             nt_bgc_DMS_in, & ! 
             nt_bgc_hum_in, & ! 
             nt_bgc_PON_in, & ! zooplankton and detritus   
             nlt_bgc_Nit_in,& ! nutrients  
             nlt_bgc_Am_in, & ! 
             nlt_bgc_Sil_in,& !
             nlt_bgc_DMSPp_in,&! trace gases (skeletal layer)
             nlt_bgc_DMSPd_in,&! 
             nlt_bgc_DMS_in,& ! 
             nlt_bgc_hum_in,& ! 
             nlt_bgc_PON_in,& ! zooplankton and detritus  
             nt_zbgc_frac_in,&! fraction of tracer in the mobile phase
             nt_bgc_S_in,   & ! Bulk salinity in fraction ice with dynamic salinity (Bio grid))
             nlt_chl_sw_in    ! points to total chla in trcrn_sw

       integer, intent(in) :: &
             n_algae,    & !  Dimensions
             n_algalC,   & !
             n_algalchl, & !
             n_DOC,      & !
             n_DON,      & !
             n_DIC,      & !
             n_dFe,      & !
             n_pFe,      & ! 
             n_aerosols, & !
             nbtrcr

        integer (kind=int_kind), dimension(:), intent(in) :: &
             bio_index_o_in, & 
             bio_index_in  

        integer (kind=int_kind), dimension(:), intent(in) :: &
             nt_bgc_N_in ,  & ! diatoms, phaeocystis, pico/small   
             nt_bgc_C_in ,  & ! diatoms, phaeocystis, pico/small   
             nt_bgc_chl_in, & ! diatoms, phaeocystis, pico/small 
             nlt_bgc_N_in , & ! diatoms, phaeocystis, pico/small   
             nlt_bgc_C_in , & ! diatoms, phaeocystis, pico/small   
             nlt_bgc_chl_in   ! diatoms, phaeocystis, pico/small 

        integer (kind=int_kind), dimension(:), intent(in) :: &
             nt_bgc_DOC_in, & !  dissolved organic carbon
             nlt_bgc_DOC_in   !  dissolved organic carbon

        integer (kind=int_kind), dimension(:), intent(in) :: &
             nt_bgc_DON_in, & !  dissolved organic nitrogen
             nlt_bgc_DON_in   !  dissolved organic nitrogen

        integer (kind=int_kind), dimension(:), intent(in) :: &
             nt_bgc_DIC_in, & ! dissolved inorganic carbon
             nlt_bgc_DIC_in   !  dissolved inorganic carbon

        integer (kind=int_kind), dimension(:), intent(in) :: &
             nt_bgc_Fed_in, & !  dissolved iron
             nt_bgc_Fep_in, & !  particulate iron
             nlt_bgc_Fed_in,& !  dissolved iron
             nlt_bgc_Fep_in   !  particulate iron

        integer (kind=int_kind), dimension(:), intent(in) :: &
             nt_zaero_in,   & !  black carbon and other aerosols
             nlt_zaero_in,  & !  black carbon and other aerosols
             nlt_zaero_sw_in  ! black carbon and dust in trcrn_sw

        ! local
        integer (kind=int_kind) :: k

        nt_Tsfc = nt_Tsfc_in
        nt_qice = nt_qice_in
        nt_qsno = nt_qsno_in
        nt_sice = nt_sice_in
        nt_fbri = nt_fbri_in
        nt_iage = nt_iage_in
        nt_FY = nt_FY_in
        nt_alvl = nt_alvl_in
        nt_vlvl = nt_vlvl_in
        nt_apnd = nt_apnd_in
        nt_hpnd = nt_hpnd_in
        nt_ipnd = nt_ipnd_in
        nt_aero = nt_aero_in
        nt_smice = nt_smice_in
        nt_smliq = nt_smliq_in
        nt_rsnw = nt_rsnw_in
        nt_rhos = nt_rhos_in
        nt_bgc_Nit = nt_bgc_Nit_in
        nt_bgc_Am  = nt_bgc_Am_in
        nt_bgc_Sil = nt_bgc_Sil_in
        nt_bgc_DMSPp=nt_bgc_DMSPp_in
        nt_bgc_DMSPd=nt_bgc_DMSPd_in
        nt_bgc_DMS = nt_bgc_DMS_in
        nt_bgc_hum = nt_bgc_hum_in
        nt_bgc_PON = nt_bgc_PON_in
        nlt_bgc_Nit = nlt_bgc_Nit_in
        nlt_bgc_Am  = nlt_bgc_Am_in
        nlt_bgc_Sil = nlt_bgc_Sil_in
        nlt_bgc_DMSPp=nlt_bgc_DMSPp_in
        nlt_bgc_DMSPd=nlt_bgc_DMSPd_in
        nlt_bgc_DMS = nlt_bgc_DMS_in
        nlt_bgc_hum = nlt_bgc_hum_in
        nlt_bgc_PON = nlt_bgc_PON_in
        nlt_chl_sw  = nlt_chl_sw_in
        nt_zbgc_frac=nt_zbgc_frac_in
        nt_bgc_S   = nt_bgc_S_in

        nt_bgc_N(:)    = 0
        nt_bgc_C(:)    = 0
        nt_bgc_chl(:)  = 0
        nlt_bgc_N(:)   = 0
        nlt_bgc_C(:)   = 0
        nlt_bgc_chl(:) = 0
        nt_bgc_DOC(:)  = 0
        nlt_bgc_DOC(:) = 0
        nt_bgc_DIC(:)  = 0
        nlt_bgc_DIC(:) = 0
        nt_bgc_DON(:)  = 0
        nlt_bgc_DON(:) = 0
        nt_bgc_Fed(:)  = 0
        nt_bgc_Fep(:)  = 0
        nlt_bgc_Fed(:) = 0
        nlt_bgc_Fep(:) = 0
        nt_zaero(:)    = 0
        nlt_zaero(:)   = 0
        nlt_zaero_sw(:)= 0
        bio_index(:)   = 0
        bio_index_o(:) = 0

        do k = 1, nbtrcr
           bio_index_o(k)= bio_index_o_in(k)
           bio_index(k)  = bio_index_in(k)
        enddo
        do k = 1, n_algae
           nt_bgc_N(k) = nt_bgc_N_in(k) 
           nlt_bgc_N(k)= nlt_bgc_N_in(k) 
        enddo
        do k = 1, n_algalC
           nt_bgc_C(k) = nt_bgc_C_in(k) 
           nlt_bgc_C(k)= nlt_bgc_C_in(k) 
        enddo
        do k = 1, n_algalchl
           nt_bgc_chl(k) = nt_bgc_chl_in(k) 
           nlt_bgc_chl(k)= nlt_bgc_chl_in(k) 
        enddo
        do k = 1, n_DOC
           nt_bgc_DOC(k) = nt_bgc_DOC_in(k) 
           nlt_bgc_DOC(k)= nlt_bgc_DOC_in(k) 
        enddo
        do k = 1, n_DON
           nt_bgc_DON(k) = nt_bgc_DON_in(k) 
           nlt_bgc_DON(k)= nlt_bgc_DON_in(k) 
        enddo
        do k = 1, n_DIC
           nt_bgc_DIC(k) = nt_bgc_DIC_in(k) 
           nlt_bgc_DIC(k)= nlt_bgc_DIC_in(k) 
        enddo
        do k = 1, n_dFe  
           nt_bgc_Fed(k) = nt_bgc_Fed_in(k) 
           nlt_bgc_Fed(k)= nlt_bgc_Fed_in(k) 
        enddo
        do k = 1, n_pFe  
           nt_bgc_Fep(k) = nt_bgc_Fep_in(k) 
           nlt_bgc_Fep(k)= nlt_bgc_Fep_in(k) 
        enddo
        do k = 1, n_aerosols
           nt_zaero(k)    = nt_zaero_in(k)   
           nlt_zaero(k)   = nlt_zaero_in(k)   
           nlt_zaero_sw(k)= nlt_zaero_sw_in(k)   
        enddo

      end subroutine colpkg_init_tracer_indices

!=======================================================================
! set the number of column tracers

      subroutine colpkg_init_tracer_numbers(&
         ntrcr_in, nbtrcr_in, nbtrcr_sw_in)

      use ice_colpkg_tracers, only: &
         ntrcr, nbtrcr, nbtrcr_sw

      integer (kind=int_kind), intent(in) :: &
         ntrcr_in  , &! number of tracers in use
         nbtrcr_in , &! number of bio tracers in use
         nbtrcr_sw_in ! number of shortwave bio tracers in use
        
         ntrcr     = ntrcr_in
         nbtrcr    = nbtrcr_in
         nbtrcr_sw = nbtrcr_sw_in

      end subroutine colpkg_init_tracer_numbers

!=======================================================================

      subroutine colpkg_biogeochemistry(dt_R4,   &
                           ntrcr,  nbtrcr,   &
                           upNO_R4,  upNH_R4,  iDi_R4,  iki_R4,  zfswin_R4,  &
                           zsal_tot_R4,  darcy_V_R4,  grow_net_R4,   &
                           PP_net_R4,  hbri_R4, dhbr_bot_R4,  dhbr_top_R4,  Zoo_R4, &
                           fbio_snoice_R4,  fbio_atmice_R4,  ocean_bio_R4,  &
                           first_ice,  fswpenln_R4,  bphi_R4,  bTiz_R4,  ice_bio_net_R4,   &
                           snow_bio_net_R4,  totalChla_R4,  fswthrun_R4,  Rayleigh_criteria,  &
                           sice_rho_R4,  fzsal_R4,  fzsal_g_R4,  &
                           bgrid_R4,  igrid_R4,  icgrid_R4,  cgrid_R4,   &
                           nblyr,  nilyr,  nslyr,  n_algae,  n_zaero,  ncat,  &
                           n_doc,  n_dic,   n_don,  n_fed,  n_fep,   &
                           meltbn_R4,  melttn_R4,  congeln_R4,  snoicen_R4,  &
                           sst_R4,  sss_R4,  Tf_R4,  fsnow_R4,  meltsn_R4,  hmix_R4,  salinz_R4,  &
                           hin_old_R4,  flux_bio_R4,  flux_bio_atm_R4,  &
                           aicen_init_R4,  vicen_init_R4,  aicen_R4,  vicen_R4,  vsnon_R4,  &
                           aice0_R4,  trcrn_R4,  vsnon_init_R4,  skl_bgc,  &
                           max_algae,  max_nbtrcr,  &
                           flux_bion_R4,  &
                           l_stop,  stop_label)

      use ice_algae, only: zbio, sklbio
      use ice_brine, only: preflushing_changes, compute_microS_mushy, &
                           update_hbrine, compute_microS
      use ice_colpkg_shared, only: solve_zsal, z_tracers, phi_snow
      use ice_colpkg_tracers, only: nt_fbri, tr_brine, &
          nt_bgc_S, nt_qice, nt_sice, nt_zbgc_frac, bio_index, bio_index_o
      use ice_constants_colpkg, only: c0, c1, puny, p5
      use ice_zsalinity, only: zsalinity
      use ice_zbgc_shared, only:  zbgc_frac_init

      real (kind=real_kind), intent(in) :: &
         dt_R4      ! time step

      integer (kind=int_kind), intent(in) :: &
         ncat, &
         nilyr, &
         nslyr, &
         nblyr, &
         ntrcr, &
         nbtrcr, &
         n_algae, n_zaero, &
         n_doc, n_dic,  n_don, n_fed, n_fep, &
         max_algae, max_nbtrcr

      real (kind=real_kind), dimension (:), intent(inout) :: &
         bgrid_R4         , &  ! biology nondimensional vertical grid points_R4
         igrid_R4         , &  ! biology vertical interface points
         cgrid_R4         , &  ! CICE vertical coordinate
         icgrid_R4        , &  ! interface grid for CICE (shortwave variable_R4)
         ocean_bio_R4     , &  ! contains all the ocean bgc tracer concentrations_R4
         fbio_snoice_R4   , &  ! fluxes from snow to ice
         fbio_atmice_R4   , &  ! fluxes from atm to ice
         dhbr_top_R4      , &  ! brine top change
         dhbr_bot_R4      , &  ! brine bottom change
         darcy_V_R4       , &  ! darcy velocity positive up (m/s)
         hin_old_R4       , &  ! old ice thickness
         sice_rho_R4      , &  ! avg sea ice density  (kg/m^3)
         ice_bio_net_R4   , &  ! depth integrated tracer (mmol/m^2)
         snow_bio_net_R4  , &  ! depth integrated snow tracer (mmol/m^2)
         flux_bio_R4           ! all bio fluxes to ocean

      logical (kind=log_kind), dimension (:), intent(inout) :: &
         first_ice      ! distinguishes ice that disappears (e.g. melts)
                        ! and reappears (e.g. transport) in a grid cell
                        ! during a single time step from ice that was
                        ! there the entire time step (true until ice forms)

      real (kind=real_kind), dimension (:,:), intent(out) :: &
         flux_bion_R4      ! per categeory ice to ocean biogeochemistry flux (mmol/m2/s)

      real (kind=real_kind), dimension (:,:), intent(inout) :: &
         Zoo_R4            , & ! N losses accumulated in timestep (ie. zooplankton/bacteria)
                            ! mmol/m^3
         bphi_R4           , & ! porosity of layers
         bTiz_R4           , & ! layer temperatures interpolated on bio grid (C)
         zfswin_R4         , & ! Shortwave flux into layers interpolated on bio grid  (W/m^2)
         iDi_R4            , & ! igrid Diffusivity (m^2/s)
         iki_R4            , & ! Ice permeability (m^2)
         trcrn_R4     ! tracers

      real (kind=real_kind), intent(inout) :: &
         grow_net_R4       , & ! Specific growth rate (/s) per grid cell
         PP_net_R4         , & ! Total production (mg C/m^2/s) per grid cell
         hbri_R4           , & ! brine height, area-averaged for comparison with hi (m)
         zsal_tot_R4       , & ! Total ice salinity in per grid cell (g/m^2)
         fzsal_R4          , & ! Total flux  of salt to ocean at time step for conservation
         fzsal_g_R4        , & ! Total gravity drainage flux
         upNO_R4           , & ! nitrate uptake rate (mmol/m^2/d) times aice
         upNH_R4           , & ! ammonium uptake rate (mmol/m^2/d) times aice
         totalChla_R4          ! ice integrated chla and summed over all algal groups (mg/m^2)

      logical (kind=log_kind), intent(inout) :: &
         Rayleigh_criteria    ! .true. means Ra_c was reached

      real (kind=real_kind), dimension (:,:), intent(in) :: &
         fswpenln_R4        ! visible SW entering ice layers (W m-2)

      real (kind=real_kind), dimension (:), intent(in) :: &
         fswthrun_R4    , & ! SW through ice to ocean            (W/m^2)
         meltsn_R4      , & ! snow melt in category n (m)
         melttn_R4      , & ! top melt in category n (m)
         meltbn_R4      , & ! bottom melt in category n (m)
         congeln_R4     , & ! congelation ice formation in category n (m)
         snoicen_R4     , & ! snow-ice formation in category n (m)
         salinz_R4      , & ! initial salinity  profile (ppt)
         flux_bio_atm_R4, & ! all bio fluxes to ice from atmosphere
         aicen_init_R4  , & ! initial ice concentration, for linear ITD
         vicen_init_R4  , & ! initial ice volume (m), for linear ITD
         vsnon_init_R4  , & ! initial snow volume (m), for aerosol
         aicen_R4 , & ! concentration of ice
         vicen_R4 , & ! volume per unit area of ice          (m)
         vsnon_R4     ! volume per unit area of snow         (m)

      real (kind=real_kind), intent(in) :: &
         aice0_R4   , & ! open water area fraction
         sss_R4     , & ! sea surface salinity (ppt)
         sst_R4     , & ! sea surface temperature (C)
         hmix_R4    , & ! mixed layer depth (m)
         Tf_R4      , & ! basal freezing temperature (C)
         fsnow_R4       ! snowfall rate (kg/m^2 s)

      logical (kind=log_kind), intent(in) :: &
         skl_bgc       ! if true, solve skeletal biochemistry

      logical (kind=log_kind), intent(inout) :: &
         l_stop          ! if true, abort the model

      character (len=*), intent(inout) :: stop_label


      ! local variables

      integer (kind=int_kind) :: &
         k              , & ! vertical index
         n, mm              ! thickness category index

      real (kind=dbl_kind) :: &
         hin         , & ! new ice thickness
         hsn         , & ! snow thickness  (m)
         hbr_old     , & ! old brine thickness before growh/melt
         dhice       , & ! change due to sublimation/condensation (m)
         kavg        , & ! average ice permeability (m^2)
         bphi_o      , & ! surface ice porosity
         hbrin       , & ! brine height
         dh_direct       ! surface flooding or runoff

      real (kind=dbl_kind), dimension (nblyr+2) :: &
      ! Defined on Bio Grid points
         bSin        , & ! salinity on the bio grid  (ppt)
         brine_sal   , & ! brine salinity (ppt)
         brine_rho       ! brine_density (kg/m^3)

      real (kind=dbl_kind), dimension (nblyr+1) :: &
      ! Defined on Bio Grid interfaces
         iphin       , & ! porosity
         ibrine_sal  , & ! brine salinity  (ppt)
         ibrine_rho  , & ! brine_density (kg/m^3)
         iTin            ! Temperature on the interface grid (oC)

      real (kind=dbl_kind) :: &
         sloss            ! brine flux contribution from surface runoff (g/m^2)

      real (kind=dbl_kind), dimension (ncat) :: &
         hbrnInitial, & ! inital brine height
         hbrnFinal      ! category initial and final brine heights

      ! for bgc sk
      real (kind=dbl_kind) :: &
         dh_bot_chl  , & ! Chlorophyll may or may not flush
         dh_top_chl  , & ! Chlorophyll may or may not flush
         darcy_V_chl

      real (kind=dbl_kind), dimension (nblyr+1) :: &
         zspace    ! vertical grid spacing


      real (kind=dbl_kind) :: &
         dt      ! time step

      real (kind=dbl_kind), dimension (:), allocatable :: &
         bgrid         , &  ! biology nondimensional vertical grid points
         igrid         , &  ! biology vertical interface points
         cgrid         , &  ! CICE vertical coordinate
         icgrid        , &  ! interface grid for CICE (shortwave variable)
         ocean_bio     , &  ! contains all the ocean bgc tracer concentrations
         fbio_snoice   , &  ! fluxes from snow to ice
         fbio_atmice   , &  ! fluxes from atm to ice
         dhbr_top      , &  ! brine top change
         dhbr_bot      , &  ! brine bottom change
         darcy_V       , &  ! darcy velocity positive up (m/s)
         hin_old       , &  ! old ice thickness
         sice_rho      , &  ! avg sea ice density  (kg/m^3)
         ice_bio_net   , &  ! depth integrated tracer (mmol/m^2)
         snow_bio_net  , &  ! depth integrated snow tracer (mmol/m^2)
         flux_bio           ! all bio fluxes to ocean

      real (kind=dbl_kind), dimension (:,:), allocatable :: &
         flux_bion      ! per categeory ice to ocean biogeochemistry flux (mmol/m2/s)

      real (kind=dbl_kind), dimension (:,:), allocatable :: &
         Zoo            , & ! N losses accumulated in timestep (ie. zooplankton/bacteria)
                            ! mmol/m^3
         bphi           , & ! porosity of layers
         bTiz           , & ! layer temperatures interpolated on bio grid (C)
         zfswin         , & ! Shortwave flux into layers interpolated on bio grid  (W/m^2)
         iDi            , & ! igrid Diffusivity (m^2/s)
         iki            , & ! Ice permeability (m^2)
         trcrn     ! tracers

      real (kind=dbl_kind) :: &
         grow_net       , & ! Specific growth rate (/s) per grid cell
         PP_net         , & ! Total production (mg C/m^2/s) per grid cell
         hbri           , & ! brine height, area-averaged for comparison with hi (m)
         zsal_tot       , & ! Total ice salinity in per grid cell (g/m^2)
         fzsal          , & ! Total flux  of salt to ocean at time step for conservation
         fzsal_g        , & ! Total gravity drainage flux
         upNO           , & ! nitrate uptake rate (mmol/m^2/d) times aice
         upNH           , & ! ammonium uptake rate (mmol/m^2/d) times aice
         totalChla          ! ice integrated chla and summed over all algal groups (mg/m^2)


      real (kind=dbl_kind), dimension (:,:), allocatable :: &
         fswpenln        ! visible SW entering ice layers (W m-2)

      real (kind=dbl_kind), dimension (:), allocatable :: &
         fswthrun    , & ! SW through ice to ocean            (W/m^2)
         meltsn      , & ! snow melt in category n (m)
         melttn      , & ! top melt in category n (m)
         meltbn      , & ! bottom melt in category n (m)
         congeln     , & ! congelation ice formation in category n (m)
         snoicen     , & ! snow-ice formation in category n (m)
         salinz      , & ! initial salinity  profile (ppt)
         flux_bio_atm, & ! all bio fluxes to ice from atmosphere
         aicen_init  , & ! initial ice concentration, for linear ITD
         vicen_init  , & ! initial ice volume (m), for linear ITD
         vsnon_init  , & ! initial snow volume (m), for aerosol
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind) :: &
         aice0   , & ! open water area fraction
         sss     , & ! sea surface salinity (ppt)
         sst     , & ! sea surface temperature (C)
         hmix    , & ! mixed layer depth (m)
         Tf      , & ! basal freezing temperature (C)
         fsnow       ! snowfall rate (kg/m^2 s)

      


      
         
      dt = real(dt_R4, kind = dbl_kind)         

     bgrid =real(bgrid_R4, kind = dbl_kind)         
     igrid =real(igrid_R4, kind = dbl_kind)         
     cgrid =real(cgrid_R4, kind = dbl_kind)         
     icgrid =real(icgrid_R4, kind = dbl_kind)        
     ocean_bio =real(ocean_bio_R4, kind = dbl_kind)     
     fbio_snoice =real(fbio_snoice_R4, kind = dbl_kind)   
     fbio_atmice =real(fbio_atmice_R4, kind = dbl_kind)   
     dhbr_top =real(dhbr_top_R4, kind = dbl_kind)      
     dhbr_bot =real(dhbr_bot_R4, kind = dbl_kind)      
     darcy_V =real(darcy_V_R4, kind = dbl_kind)       
     hin_old =real(hin_old_R4, kind = dbl_kind)       
     sice_rho =real(sice_rho_R4, kind = dbl_kind)      
     ice_bio_net =real(ice_bio_net_R4, kind = dbl_kind)   
     snow_bio_net =real(snow_bio_net_R4, kind = dbl_kind)  
     flux_bio =real(flux_bio_R4, kind = dbl_kind)           

     flux_bion  = real(flux_bion_R4, kind = dbl_kind)      

     Zoo  = real(Zoo_R4, kind = dbl_kind)            
     bphi  = real(bphi_R4, kind = dbl_kind)           
     bTiz  = real(bTiz_R4, kind = dbl_kind)           
     zfswin  = real(zfswin_R4, kind = dbl_kind)         
     iDi  = real(iDi_R4, kind = dbl_kind)            
     iki  = real(iki_R4, kind = dbl_kind)            
     trcrn  = real(trcrn_R4, kind = dbl_kind)     
     fswpenln  = real(fswpenln_R4, kind = dbl_kind)        
   
      grow_net = real(grow_net_R4, kind = dbl_kind)       
      PP_net = real(PP_net_R4, kind = dbl_kind)         
      hbri = real(hbri_R4, kind = dbl_kind)           
      zsal_tot = real(zsal_tot_R4, kind = dbl_kind)       
      fzsal = real(fzsal_R4, kind = dbl_kind)          
      fzsal_g = real(fzsal_g_R4, kind = dbl_kind)        
      upNO = real(upNO_R4, kind = dbl_kind)           
      upNH = real(upNH_R4, kind = dbl_kind)           
      totalChla = real(totalChla_R4, kind = dbl_kind) 

     fswthrun  = real(fswthrun_R4, kind = dbl_kind)    
     meltsn  = real(meltsn_R4, kind = dbl_kind)      
     melttn  = real(melttn_R4, kind = dbl_kind)      
     meltbn  = real(meltbn_R4, kind = dbl_kind)      
     congeln  = real(congeln_R4, kind = dbl_kind)     
     snoicen  = real(snoicen_R4, kind = dbl_kind)     
     salinz  = real(salinz_R4, kind = dbl_kind)      
     flux_bio_atm  = real(flux_bio_atm_R4, kind = dbl_kind)
     aicen_init  = real(aicen_init_R4, kind = dbl_kind)  
     vicen_init  = real(vicen_init_R4, kind = dbl_kind)  
     vsnon_init  = real(vsnon_init_R4, kind = dbl_kind)  
     aicen  = real(aicen_R4, kind = dbl_kind) 
     vicen  = real(vicen_R4, kind = dbl_kind) 
     vsnon  = real(vsnon_R4, kind = dbl_kind)     

      aice0 = real(aice0_R4, kind = dbl_kind)   
      sss = real(sss_R4, kind = dbl_kind)     
      sst = real(sst_R4, kind = dbl_kind)     
      hmix = real(hmix_R4, kind = dbl_kind)    
      Tf = real(Tf_R4, kind = dbl_kind)      
      fsnow = real(fsnow_R4, kind = dbl_kind)       

      zspace(:)       = c1/real(nblyr,kind=dbl_kind)
      zspace(1)       = p5*zspace(1)
      zspace(nblyr+1) = p5*zspace(nblyr+1)

      l_stop = .false.

      do n = 1, ncat

      !-----------------------------------------------------------------
      ! initialize
      !-----------------------------------------------------------------
         flux_bion(:,n) = c0
         hin_old(n) = c0
         hbrnFinal(n) = c0
         hbrnInitial(n) = c0

         if (aicen_init(n) > puny) then
            hin_old(n) = vicen_init(n) &
                                / aicen_init(n)
         else

            first_ice(n) = .true.
            if (tr_brine) trcrn(nt_fbri,n) = c1
            do mm = 1,nbtrcr
               trcrn(nt_zbgc_frac-1+mm,n) = zbgc_frac_init(mm)
            enddo
            if (n == 1) Rayleigh_criteria = .false.
            if (solve_zsal) trcrn(nt_bgc_S:nt_bgc_S+nblyr-1,n) = c0
         endif

         if (aicen(n) > puny) then

            dh_top_chl = c0
            dh_bot_chl = c0
            darcy_V_chl= c0
            bSin(:)    = c0
            hsn        = c0
            hin        = c0
            hbrin      = c0
            kavg       = c0
            bphi_o     = c0
            sloss      = c0

      !-----------------------------------------------------------------
      ! brine dynamics
      !-----------------------------------------------------------------

            dhbr_top(n) = c0
            dhbr_bot(n) = c0

            if (tr_brine) then

               dhice = c0
               call preflushing_changes  (n,  aicen  (n),   &
                                 vicen   (n), vsnon  (n),   &
                                 meltbn  (n), melttn (n),   &
                                 congeln (n), snoicen(n),   &
                                 hin_old (n), dhice,        &
                                 trcrn(nt_fbri,n),          &
                                 dhbr_top(n), dhbr_bot(n),  &
                                 hbr_old,     hin,          &
                                 hsn,         first_ice(n), &
                                 l_stop,      stop_label)

               hbrnInitial(n) = hbr_old

               if (l_stop) return

               if (solve_zsal)  then

                  call compute_microS (n,         nilyr,       nblyr,             &
                                bgrid,            cgrid,       igrid,             &
                                trcrn(1:ntrcr,n), hin_old(n),  hbr_old,           &
                                sss,              sst,         bTiz(:,n),         &
                                iTin,             bphi(:,n),   kavg,              &
                                bphi_o,           phi_snow,    Rayleigh_criteria, &
                                first_ice(n),     bSin,        brine_sal,         &
                                brine_rho,        iphin,       ibrine_rho,        &
                                ibrine_sal,       sice_rho(n), sloss,             &
                                salinz(1:nilyr),  l_stop,      stop_label)

                  if (l_stop) return
               else

                 ! Requires the average ice permeability = kavg(:)
                 ! and the surface ice porosity = zphi_o(:)
                 ! computed in "compute_microS" or from "thermosaline_vertical"

                  iDi(:,n) = c0

                  call compute_microS_mushy (n,   nilyr,         nblyr,       &
                                   bgrid,         cgrid,         igrid,       &
                                   trcrn(:,n),    hin_old(n),    hbr_old,     &
                                   sss,           sst,           bTiz(:,n),   &
                                   iTin(:),       bphi(:,n),     kavg,        &
                                   bphi_o,        phi_snow,      bSin(:),     &
                                   brine_sal(:),  brine_rho(:),  iphin(:),    &
                                   ibrine_rho(:), ibrine_sal(:), sice_rho(n), &
                                   iDi(:,n),      l_stop,        stop_label)

               endif ! solve_zsal

               call update_hbrine (meltbn  (n), melttn(n),   &
                                   meltsn  (n), dt,          &
                                   hin,         hsn,         &
                                   hin_old (n), hbrin,       &

                                   hbr_old,     phi_snow,    &
                                   trcrn(nt_fbri,n),         &
                                   snoicen(n),               &
                                   dhbr_top(n), dhbr_bot(n), &
                                   dh_top_chl,  dh_bot_chl,  &
                                   kavg,        bphi_o,      &
                                   darcy_V (n), darcy_V_chl, &
                                   bphi(2,n),   aice0,       &
                                   dh_direct)

               hbri = hbri + hbrin * aicen(n)
               hbrnFinal(n) = hbrin

               if (solve_zsal) then

                  call zsalinity (n,             dt,                  &
                                  nilyr,         bgrid,               &
                                  cgrid,         igrid,               &
                                  trcrn(nt_bgc_S:nt_bgc_S+nblyr-1,n), &
                                  trcrn(nt_qice:nt_qice+nilyr-1,n),   &
                                  trcrn(nt_sice:nt_sice+nilyr-1,n),   &
                                  ntrcr,         trcrn(nt_fbri,n),    &
                                  bSin,          bTiz(:,n),           &
                                  bphi(:,n),     iphin,               &
                                  iki(:,n),      hbr_old,             &
                                  hbrin,         hin,                 &
                                  hin_old(n),    iDi(:,n),            &
                                  darcy_V(n),    brine_sal,           &
                                  brine_rho,     ibrine_sal,          &
                                  ibrine_rho,    dh_direct,           &
                                  Rayleigh_criteria,                  &
                                  first_ice(n),  sss,                 &
                                  sst,           dhbr_top(n),         &
                                  dhbr_bot(n),                        &
                                  l_stop,        stop_label,          &
                                  fzsal,         fzsal_g,             &
                                  bphi_o,        nblyr,               &
                                  vicen(n),      aicen_init(n),       &
                                  zsal_tot)

                  if (l_stop) return

               endif  ! solve_zsal

            endif ! tr_brine

      !-----------------------------------------------------------------
      ! biogeochemistry
      !-----------------------------------------------------------------

            if (z_tracers) then

               call zbio (dt,                    nblyr,                  &
                          nslyr,                 nilyr,                  &
                          melttn(n),                                     &
                          meltsn(n),             meltbn  (n),            &
                          congeln(n),            snoicen(n),             &
                          nbtrcr,                fsnow,                  &
                          ntrcr,                 trcrn(1:ntrcr,n),       &
                          bio_index(1:nbtrcr),   bio_index_o(:),         &
                          aicen_init(n),                                 &
                          vicen_init(n),         vsnon_init(n),          &
                          vicen(n),              vsnon(n),               &
                          aicen(n),              flux_bio_atm(:), &
                          n,                     n_algae,                &
                          n_doc,                 n_dic,                  &
                          n_don,                                         &
                          n_fed,                 n_fep,                  &
                          n_zaero,               first_ice(n),           &
                          hin_old(n),            ocean_bio(:),    &
                          bphi(:,n),             iphin,                  &
                          iDi(:,n),              sss,                    &
                          fswpenln(:,n),                                 &
                          dhbr_top(n),           dhbr_bot(n),            &
                          dh_top_chl,            dh_bot_chl,             &
                          zfswin(:,n),                                   &
                          hbrin,                 hbr_old,                &
                          darcy_V(n),            darcy_V_chl,            &
                          bgrid,                 cgrid,                  &
                          igrid,                 icgrid,                 &
                          bphi_o,                                        &
                          dhice,                 iTin,                   &
                          Zoo(:,n),                                      &
                          flux_bio(:),           dh_direct,              &
                          upNO,                  upNH,                   &
                          fbio_snoice,           fbio_atmice,            &
                          PP_net,                ice_bio_net (:),        &
                          snow_bio_net(:),       grow_net,               &
                          totalChla,                                     &
                          flux_bion(:,n),                                &
                          l_stop,                stop_label)

               if (l_stop) return

            elseif (skl_bgc) then

               call sklbio (dt,                      Tf,                  &
                            ntrcr,                   nilyr,               &
                            nbtrcr,                  n_algae,             &
                            n_zaero,                 n_doc,               &
                            n_dic,                   n_don,               &
                            n_fed,                   n_fep,               &
                            flux_bio (1:nbtrcr),     ocean_bio(:),        &
                            hmix,                    aicen    (n),        &
                            meltbn   (n),            congeln  (n),        &
                            fswthrun (n),            first_ice(n),        &
                            trcrn    (1:ntrcr,n),    hin,                 &
                            PP_net,                  upNO,                &
                            upNH,                    grow_net,            &
                            totalChla,                                    &
                            l_stop,                  stop_label)

               if (l_stop) return

            endif  ! skl_bgc

            first_ice(n) = .false.
         else
            do mm = 1, nbtrcr
               do k  = 1, nblyr+1
                  flux_bion(mm,n) = flux_bion(mm,n) + trcrn(bio_index(mm) + k-1,n) *  &
                     hin_old(n) * zspace(k)/dt * trcrn(nt_fbri,n)
                  flux_bio(mm) = flux_bio(mm) + trcrn(bio_index(mm) + k-1,n) * &
                     vicen_init(n) * zspace(k)/dt * trcrn(nt_fbri,n)
                  trcrn(bio_index(mm) + k-1,n) = c0
                enddo
            enddo
         endif             ! aicen > puny
      enddo                ! ncat

      
      ! dt_R4 = real(dt, kind = real_kind)      
      
      bgrid_R4 = real(bgrid, kind = real_kind)         
      igrid_R4 = real(igrid, kind = real_kind)         
      cgrid_R4 = real(cgrid, kind = real_kind)         
      icgrid_R4 = real(icgrid, kind = real_kind)        
      ocean_bio_R4 = real(ocean_bio, kind = real_kind)     
      fbio_snoice_R4 = real(fbio_snoice, kind = real_kind)   
      fbio_atmice_R4 = real(fbio_atmice, kind = real_kind)   
      dhbr_top_R4 = real(dhbr_top, kind = real_kind)      
      dhbr_bot_R4 = real(dhbr_bot, kind = real_kind)      
      darcy_V_R4 = real(darcy_V, kind = real_kind)       
      hin_old_R4 = real(hin_old, kind = real_kind)       
      sice_rho_R4 = real(sice_rho, kind = real_kind)      
      ice_bio_net_R4 = real(ice_bio_net, kind = real_kind)   
      snow_bio_net_R4 = real(snow_bio_net, kind = real_kind)  
      flux_bio_R4 = real(flux_bio, kind = real_kind)           

   
      flux_bion_R4 = real(flux_bion, kind = real_kind)      

   
      Zoo_R4 = real(Zoo, kind = real_kind)            
      bphi_R4 = real(bphi, kind = real_kind)           
      bTiz_R4 = real(bTiz, kind = real_kind)           
      zfswin_R4 = real(zfswin, kind = real_kind)         
      iDi_R4 = real(iDi, kind = real_kind)            
      iki_R4 = real(iki, kind = real_kind)            
      trcrn_R4 = real(trcrn, kind = real_kind)     

   
      grow_net_R4 = real(grow_net, kind = real_kind)       
      PP_net_R4 = real(PP_net, kind = real_kind)         
      hbri_R4 = real(hbri, kind = real_kind)           
      zsal_tot_R4 = real(zsal_tot, kind = real_kind)       
      fzsal_R4 = real(fzsal, kind = real_kind)          
      fzsal_g_R4 = real(fzsal_g, kind = real_kind)        
      upNO_R4 = real(upNO, kind = real_kind)           
      upNH_R4 = real(upNH, kind = real_kind)           
      totalChla_R4 = real(totalChla, kind = real_kind)          



   
      ! fswpenln_R4 = real(fswpenln, kind = real_kind)        

   
      ! fswthrun_R4 = real(fswthrun, kind = real_kind)    
      ! meltsn_R4 = real(meltsn, kind = real_kind)      
      ! melttn_R4 = real(melttn, kind = real_kind)      
      ! meltbn_R4 = real(meltbn, kind = real_kind)      
      ! congeln_R4 = real(congeln, kind = real_kind)     
      ! snoicen_R4 = real(snoicen, kind = real_kind)     
      ! salinz_R4 = real(salinz, kind = real_kind)      
      ! flux_bio_atm_R4 = real(flux_bio_atm, kind = real_kind)
      ! aicen_init_R4 = real(aicen_init, kind = real_kind)  
      ! vicen_init_R4 = real(vicen_init, kind = real_kind)  
      ! vsnon_init_R4 = real(vsnon_init, kind = real_kind)  
      ! aicen_R4 = real(aicen, kind = real_kind) 
      ! vicen_R4 = real(vicen, kind = real_kind) 
      ! vsnon_R4 = real(vsnon, kind = real_kind)     

   
      ! aice0_R4 = real(aice0, kind = real_kind)   
      ! sss_R4 = real(sss, kind = real_kind)     
      ! sst_R4 = real(sst, kind = real_kind)     
      ! hmix_R4 = real(hmix, kind = real_kind)    
      ! Tf_R4 = real(Tf, kind = real_kind)      
      ! fsnow_R4 = real(fsnow, kind = real_kind)       

      deallocate(bgrid)
      deallocate(igrid)
      deallocate(cgrid)
      deallocate(icgrid)
      deallocate(ocean_bio)
      deallocate(fbio_snoice)
      deallocate(fbio_atmice)
      deallocate(dhbr_top)
      deallocate(dhbr_bot)
      deallocate(darcy_V)
      deallocate(hin_old)
      deallocate(sice_rho)
      deallocate(ice_bio_net)
      deallocate(snow_bio_net)
      deallocate(flux_bio)

      deallocate(flux_bion)
      deallocate(Zoo)
      deallocate(bphi)
      deallocate(bTiz)
      deallocate(zfswin)
      deallocate(iDi)
      deallocate(iki)
      deallocate(trcrn)
      deallocate(fswpenln)
      deallocate(fswthrun)
      deallocate(meltsn)
      deallocate(melttn)
      deallocate(meltbn)
      deallocate(congeln)
      deallocate(snoicen)
      deallocate(salinz)
      deallocate(flux_bio_atm)
      deallocate(aicen_init)
      deallocate(vicen_init)
      deallocate(vsnon_init)
      deallocate(aicen)
      deallocate(vicen)
      deallocate(vsnon)
      end subroutine colpkg_biogeochemistry

!=======================================================================

!  Initialize brine height tracer

      subroutine colpkg_init_hbrine(bgrid_R4, igrid_R4, cgrid_R4, &
          icgrid_R4, swgrid_R4, nblyr, nilyr, phi_snow_R4)

      use ice_constants_colpkg, only: c1, c1p5, c2, p5, c0, rhoi, rhos, p25

      integer (kind=int_kind), intent(in) :: &
         nilyr, & ! number of ice layers
         nblyr    ! number of bio layers

      real (kind=real_kind), intent(inout) :: &
         phi_snow_R4           !porosity at the ice-snow interface

      real (kind=real_kind), dimension (nblyr+2), intent(out) :: &
         bgrid_R4              ! biology nondimensional vertical grid points

      real (kind=real_kind), dimension (nblyr+1), intent(out) :: &
         igrid_R4              ! biology vertical interface points

      real (kind=real_kind), dimension (nilyr+1), intent(out) :: &
         cgrid_R4            , &  ! CICE vertical coordinate
         icgrid_R4           , &  ! interface grid for CICE (shortwave variable)
         swgrid_R4                ! grid for ice tracers used in dEdd scheme


      integer (kind=int_kind) :: &
         k           , & ! vertical index
         n               ! thickness category index

      real (kind=dbl_kind) :: &
         zspace            ! grid spacing for CICE vertical grid

      real (kind=dbl_kind) :: &
         phi_snow           !porosity at the ice-snow interface

      real (kind=dbl_kind), dimension (nblyr+2) :: &
         bgrid              ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1) :: &
         igrid              ! biology vertical interface points

      real (kind=dbl_kind), dimension (nilyr+1) :: &
         cgrid            , &  ! CICE vertical coordinate
         icgrid           , &  ! interface grid for CICE (shortwave variable)
         swgrid                ! grid for ice tracers used in dEdd scheme

      phi_snow = real(phi_snow_R4, kind = dbl_kind)           
      bgrid = real(bgrid_R4, kind = dbl_kind)              
      igrid = real(igrid_R4, kind = dbl_kind)              
      cgrid = real(cgrid_R4, kind = dbl_kind)            
      icgrid = real(icgrid_R4, kind = dbl_kind)           
      swgrid = real(swgrid_R4, kind = dbl_kind)                



      if (phi_snow .le. c0) phi_snow = c1-rhos/rhoi

      !-----------------------------------------------------------------
      ! Calculate bio gridn: 0 to 1 corresponds to ice top to bottom
      !-----------------------------------------------------------------

      bgrid(:)       = c0 ! zsalinity grid points
      bgrid(nblyr+2) = c1 ! bottom value
      igrid(:)       = c0 ! bgc interface grid points
      igrid(1)       = c0 ! ice top
      igrid(nblyr+1) = c1 ! ice bottom

      zspace = c1/max(c1,(real(nblyr,kind=dbl_kind)))
      do k = 2, nblyr+1
         bgrid(k) = zspace*(real(k,kind=dbl_kind) - c1p5)
      enddo

      do k = 2, nblyr
         igrid(k) = p5*(bgrid(k+1)+bgrid(k))
      enddo

      !-----------------------------------------------------------------
      ! Calculate CICE cgrid for interpolation ice top (0) to bottom (1)
      !-----------------------------------------------------------------

      cgrid(1) = c0                           ! CICE vertical grid top point
      zspace = c1/(real(nilyr,kind=dbl_kind)) ! CICE grid spacing

      do k = 2, nilyr+1
         cgrid(k) = zspace * (real(k,kind=dbl_kind) - c1p5)
      enddo

      !-----------------------------------------------------------------
      ! Calculate CICE icgrid for ishortwave interpolation top(0) , bottom (1)
      !-----------------------------------------------------------------

      icgrid(1) = c0
      zspace = c1/(real(nilyr,kind=dbl_kind)) ! CICE grid spacing

      do k = 2, nilyr+1
         icgrid(k) = zspace * (real(k,kind=dbl_kind)-c1)
      enddo

      !------------------------------------------------------------------------
      ! Calculate CICE swgrid for dEdd ice: top of ice (0) , bottom of ice (1)
      ! Does not include snow
      ! see ice_shortwave.F90
      ! swgrid represents the layer index of the delta-eddington ice layer index
      !------------------------------------------------------------------------
      zspace = c1/(real(nilyr,kind=dbl_kind)) ! CICE grid spacing
      swgrid(1) = min(c1/60.0_dbl_kind, zspace*p25)
      swgrid(2) = zspace/c2                   !+ swgrid(1)
      do k = 3, nilyr+1
         swgrid(k) = zspace * (real(k,kind=dbl_kind)-c1p5)
      enddo

      phi_snow_R4 = real(phi_snow, kind = real_kind)           
      bgrid_R4 = real(bgrid, kind = real_kind)              
      igrid_R4 = real(igrid, kind = real_kind)              
      cgrid_R4 = real(cgrid, kind = real_kind)            
      icgrid_R4 = real(icgrid, kind = real_kind)           
      swgrid_R4 = real(swgrid, kind = real_kind)                



      end subroutine colpkg_init_hbrine

!=======================================================================

!  Initialize ocean concentration

      subroutine colpkg_init_ocean_conc (amm, dmsp, dms, algalN, doc, dic, don, &
             fed, fep, hum, nit, sil, zaeros, max_dic, max_don, max_fe, max_aero,&
             CToN, CToN_DON)

      use ice_constants_colpkg, only: c1,  c2, p5, c0, p1
      use ice_colpkg_shared, only: R_C2N, R_C2N_DON
 
      integer (kind=int_kind), intent(in) :: &
        max_dic, &
        max_don, &
        max_fe, &
        max_aero

      real (kind=real_kind), intent(out):: &
       amm      , & ! ammonium
       dmsp     , & ! DMSPp
       dms      , & ! DMS
       hum      , & ! humic material
       nit      , & ! nitrate
       sil          ! silicate

      real (kind=real_kind), dimension(:), intent(out):: &
       algalN   , & ! algae
       doc      , & ! DOC
       dic      , & ! DIC
       don      , & ! DON
       fed      , & ! Dissolved Iron
       fep      , & ! Particulate Iron
       zaeros       ! BC and dust

      real (kind=real_kind), dimension(:), intent(inout), optional :: &
       CToN     , & ! carbon to nitrogen ratio for algae
       CToN_DON     ! nitrogen to carbon ratio for proteins

      integer (kind=int_kind) :: &
        k 

       if (present(CToN)) then
         CToN(1) = real(R_C2N(1), kind = real_kind)
         CToN(2) = real(R_C2N(2), kind = real_kind)    
         CToN(3) = real(R_C2N(3), kind = real_kind)     
       endif

       if (present(CToN_DON)) then
         CToN_DON(1) = real(R_C2N_DON(1), kind=real_kind)
       endif

       amm  = real(c1, kind=real_kind) ! ISPOL < 1 mmol/m^3 
       dmsp = real(p1, kind=real_kind)  
       dms  = real(p1, kind=real_kind)    
       algalN(1) = c1  !0.0026_real_kind ! ISPOL, Lannuzel 2013(pennate) 
       algalN(2) = 0.0057_real_kind ! ISPOL, Lannuzel 2013(small plankton)
       algalN(3) = 0.0027_real_kind ! ISPOL, Lannuzel 2013(Phaeocystis)
                                     ! 0.024_real_kind ! 5% of 1 mgchl/m^3 
       doc(1) = 16.2_real_kind ! 18% saccharides
       doc(2) = 9.0_real_kind  ! lipids
       doc(3) = real(c1, kind=real_kind) ! 
       do k = 1, max_dic
            dic(k) = 1950.0_real_kind  ! 1950-2260 mmol C/m3 (Tynan et al. 2015)
       enddo  
       do k = 1, max_don
            don(k) = 12.9_real_kind              
            ! 64.3_real_kind ! 72% Total DOC~90 mmolC/m^3  ISPOL with N:C of 0.2
       enddo  
       !ki = 1
       !if (trim(fe_data_type) == 'clim') ki = 2
       do k = 1, max_fe ! ki, max_fe
            fed(k) = 0.4_real_kind ! c1 (nM) Lannuzel2007 DFe, 
                                  ! range 0.14-2.6 (nM) van der Merwe 2011
                                  ! Tagliabue 2012 (0.4 nM)
            fep(k) = real(c2, kind=real_kind)  ! (nM) van der Merwe 2011
                        ! (0.6 to 2.9 nM ocean)
       enddo 
       hum  = real(c1, kind=real_kind)         ! mmol C/m^3
       nit  = 12.0_real_kind
       sil  = 25.0_real_kind
       do k = 1, max_aero
         zaeros(k) = real(c0, kind=real_kind) 
       enddo
 

      end subroutine colpkg_init_ocean_conc

      subroutine colpkg_init_ocean_conc_double (amm, dmsp, dms, algalN, doc, dic, don, &
         fed, fep, hum, nit, sil, zaeros, max_dic, max_don, max_fe, max_aero,&
         CToN, CToN_DON)

      use ice_constants_colpkg, only: c1,  c2, p5, c0, p1
      use ice_colpkg_shared, only: R_C2N, R_C2N_DON

      integer (kind=int_kind), intent(in) :: &
         max_dic, &
         max_don, &
         max_fe, &
         max_aero

      real (kind=dbl_kind), intent(out):: &
         amm      , & ! ammonium
         dmsp     , & ! DMSPp
         dms      , & ! DMS
         hum      , & ! humic material
         nit      , & ! nitrate
         sil          ! silicate

      real (kind=dbl_kind), dimension(:), intent(out):: &
         algalN   , & ! algae
         doc      , & ! DOC
         dic      , & ! DIC
         don      , & ! DON
         fed      , & ! Dissolved Iron
         fep      , & ! Particulate Iron
         zaeros       ! BC and dust

      real (kind=dbl_kind), dimension(:), intent(inout), optional :: &
         CToN     , & ! carbon to nitrogen ratio for algae
         CToN_DON     ! nitrogen to carbon ratio for proteins

      integer (kind=int_kind) :: &
         k 

         if (present(CToN)) then
         CToN(1) = R_C2N(1)
         CToN(2) = R_C2N(2)    
         CToN(3) = R_C2N(3)     
         endif

         if (present(CToN_DON)) then
         CToN_DON(1) = R_C2N_DON(1)
         endif

         amm  = c1 ! ISPOL < 1 mmol/m^3 
         dmsp = p1  
         dms  = p1    
         algalN(1) = c1  !0.0026_dbl_kind ! ISPOL, Lannuzel 2013(pennate) 
         algalN(2) = 0.0057_dbl_kind ! ISPOL, Lannuzel 2013(small plankton)
         algalN(3) = 0.0027_dbl_kind ! ISPOL, Lannuzel 2013(Phaeocystis)
                                       ! 0.024_dbl_kind ! 5% of 1 mgchl/m^3 
         doc(1) = 16.2_dbl_kind ! 18% saccharides
         doc(2) = 9.0_dbl_kind  ! lipids
         doc(3) = c1 ! 
         do k = 1, max_dic
            dic(k) = 1950.0_dbl_kind  ! 1950-2260 mmol C/m3 (Tynan et al. 2015)
         enddo  
         do k = 1, max_don
            don(k) = 12.9_dbl_kind              
            ! 64.3_dbl_kind ! 72% Total DOC~90 mmolC/m^3  ISPOL with N:C of 0.2
         enddo  
         !ki = 1
         !if (trim(fe_data_type) == 'clim') ki = 2
         do k = 1, max_fe ! ki, max_fe
            fed(k) = 0.4_dbl_kind ! c1 (nM) Lannuzel2007 DFe, 
                                    ! range 0.14-2.6 (nM) van der Merwe 2011
                                    ! Tagliabue 2012 (0.4 nM)
            fep(k) = c2 ! (nM) van der Merwe 2011
                        ! (0.6 to 2.9 nM ocean)
         enddo 
         hum  = c1        ! mmol C/m^3
         nit  = 12.0_dbl_kind
         sil  = 25.0_dbl_kind
         do k = 1, max_aero
         zaeros(k) = c0
         enddo


      end subroutine colpkg_init_ocean_conc_double
!=======================================================================

!  Initialize zSalinity

      subroutine colpkg_init_zsalinity(nblyr,ntrcr_o, restart_zsal,  Rayleigh_criteria, &
               Rayleigh_real_R4, trcrn_R4, nt_bgc_S, ncat, sss_R4)

      use ice_constants_colpkg, only: c1,  c2, p5, c0, p1
      use ice_colpkg_shared, only: dts_b, salt_loss
 
      integer (kind=int_kind), intent(in) :: &
       nblyr, & ! number of biolayers
       ntrcr_o, & ! number of non bio tracers
       ncat , & ! number of categories
       nt_bgc_S ! zsalinity index

      logical (kind=log_kind), intent(in) :: &
       restart_zsal

      logical (kind=log_kind), intent(inout) :: &
       Rayleigh_criteria

      real (kind=real_kind), intent(inout):: &
       Rayleigh_real_R4

      real (kind=real_kind), intent(in):: &
       sss_R4

      real (kind=real_kind), dimension(:,:), intent(inout):: &
       trcrn_R4 ! bgc subset of trcrn

      real (kind=dbl_kind), dimension(:,:), allocatable:: &
       trcrn ! bgc subset of trcrn

      integer (kind=int_kind) :: &
        k , n

      real (kind=dbl_kind):: &
       Rayleigh_real

      real (kind=dbl_kind):: &
       sss

      Rayleigh_real = real(Rayleigh_real_R4, kind = dbl_kind)
      sss = real(sss_R4, kind = dbl_kind)
      trcrn  = real(trcrn_R4, kind = dbl_kind)

      if (nblyr .LE. 7) then
          dts_b = 300.0_dbl_kind
      else
          dts_b = 50.0_dbl_kind 
      endif

      if (.not. restart_zsal) then
         Rayleigh_criteria = .false.    ! no ice initial condition 
         Rayleigh_real     = c0
         do n = 1,ncat
             do k = 1,nblyr
                trcrn(nt_bgc_S+k-1-ntrcr_o,n) = sss*salt_loss
             enddo   ! k
         enddo      ! n
      endif

      trcrn_R4 = real(trcrn, kind = real_kind)
      deallocate(trcrn)
      end subroutine colpkg_init_zsalinity

!=======================================================================

! basic initialization for ocean_bio_all

      subroutine colpkg_init_OceanConcArray(max_nbtrcr, &
          max_algae, max_don, max_doc, max_dic, max_aero, max_fe, &
          nit, amm, sil, dmsp, dms, algalN, &
          doc, don, dic, fed, fep, zaeros, ocean_bio_all_R4, hum)

      use ice_constants_colpkg, only: c0
      use ice_colpkg_shared, only:  R_C2N, R_chl2N
      use ice_zbgc_shared, only: R_S2N

      integer (kind=int_kind), intent(in) :: &
         max_algae   , & ! maximum number of algal types 
         max_dic     , & ! maximum number of dissolved inorganic carbon types 
         max_doc     , & ! maximum number of dissolved organic carbon types
         max_don     , & ! maximum number of dissolved organic nitrogen types
         max_fe      , & ! maximum number of iron types
         max_aero    , & ! maximum number of aerosols 
         max_nbtrcr      ! maximum number of bio tracers

      real (kind=real_kind), intent(in) :: &
         nit         , & ! ocean nitrate (mmol/m^3)          
         amm         , & ! ammonia/um (mmol/m^3)
         sil         , & ! silicate (mmol/m^3)
         dmsp        , & ! dmsp (mmol/m^3)
         dms         , & ! dms (mmol/m^3)
         hum             ! humic material (mmol/m^3)

      real (kind=real_kind), dimension (max_algae), intent(in) :: &
         algalN          ! ocean algal nitrogen (mmol/m^3) (diatoms, phaeo, pico)

      real (kind=real_kind), dimension (max_doc), intent(in) :: &
         doc             ! ocean doc (mmol/m^3)  (proteins, EPS, lipid)

      real (kind=real_kind), dimension (max_don), intent(in) :: &
         don             ! ocean don (mmol/m^3) 

      real (kind=real_kind), dimension (max_dic), intent(in) :: &
         dic             ! ocean dic (mmol/m^3) 

      real (kind=real_kind), dimension (max_fe), intent(in) :: &
         fed, fep        ! ocean disolved and particulate fe (nM) 

      real (kind=real_kind), dimension (max_aero), intent(in) :: &
         zaeros          ! ocean aerosols (mmol/m^3) 

      real (kind=real_kind), dimension (max_nbtrcr), intent(inout) :: &
         ocean_bio_all_R4   ! fixed order, all values even for tracers false

      real (kind=dbl_kind), dimension (max_nbtrcr) :: &
         ocean_bio_all   ! fixed order, all values even for tracers false
      ! local variables

      integer (kind=int_kind) :: &
         k, ks           ! tracer indices

      ocean_bio_all(:) = c0

      ocean_bio_all = real(ocean_bio_all_R4, kind=dbl_kind)

      do k = 1, max_algae           
         ocean_bio_all(k)      = algalN(k)           ! N
         ks = max_algae + max_doc + max_dic + 1
         ocean_bio_all(ks + k) = R_chl2N(k)*algalN(k)!chl
      enddo   

      ks = max_algae + 1
      do k = 1, max_doc
         ocean_bio_all(ks + k) = doc(k)              ! doc
      enddo  
      ks = ks + max_doc
      do k = 1, max_dic
         ocean_bio_all(ks + k) = dic(k)              ! dic
      enddo 

      ks = 2*max_algae + max_doc + max_dic + 7
      do k = 1, max_don
         ocean_bio_all(ks + k) = don(k)              ! don
      enddo  

      ks = max_algae + 1
      ocean_bio_all(ks) = nit                        ! nit

      ks = 2*max_algae + max_doc + 2 + max_dic
      ocean_bio_all(ks) = amm                        ! Am
      ks = ks + 1
      ocean_bio_all(ks) = sil                        ! Sil
      ks = ks + 1
      ocean_bio_all(ks) =  R_S2N(1)*algalN(1) &      ! DMSPp
                        +  R_S2N(2)*algalN(2) &
                        +  R_S2N(3)*algalN(3) 
      ks = ks + 1
      ocean_bio_all(ks) = dmsp                       ! DMSPd
      ks = ks + 1
      ocean_bio_all(ks) = dms                        ! DMS
      ks = ks + 1
      ocean_bio_all(ks) = nit                        ! PON
      ks = 2*max_algae + max_doc + 7 + max_dic + max_don
      do k = 1, max_fe
         ocean_bio_all(ks + k) = fed(k)              ! fed
      enddo  
      ks = ks + max_fe
      do k = 1, max_fe
         ocean_bio_all(ks + k) = fep(k)              ! fep
      enddo  
      ks = ks + max_fe
      do k = 1, max_aero
         ocean_bio_all(ks+k) = zaeros(k)             ! zaero
      enddo
      ks = ks + max_aero + 1 
      ocean_bio_all(ks)  = hum                       ! humics
      ocean_bio_all_R4 = real(ocean_bio_all, kind=real_kind)

      end subroutine colpkg_init_OceanConcArray

      subroutine colpkg_init_OceanConcArray_double(max_nbtrcr, &
         max_algae, max_don, max_doc, max_dic, max_aero, max_fe, &
         nit, amm, sil, dmsp, dms, algalN, &
         doc, don, dic, fed, fep, zaeros, ocean_bio_all, hum)

      use ice_constants_colpkg, only: c0
      use ice_colpkg_shared, only:  R_C2N, R_chl2N
      use ice_zbgc_shared, only: R_S2N

      integer (kind=int_kind), intent(in) :: &
         max_algae   , & ! maximum number of algal types 
         max_dic     , & ! maximum number of dissolved inorganic carbon types 
         max_doc     , & ! maximum number of dissolved organic carbon types
         max_don     , & ! maximum number of dissolved organic nitrogen types
         max_fe      , & ! maximum number of iron types
         max_aero    , & ! maximum number of aerosols 
         max_nbtrcr      ! maximum number of bio tracers

      real (kind=dbl_kind), intent(in) :: &
         nit         , & ! ocean nitrate (mmol/m^3)          
         amm         , & ! ammonia/um (mmol/m^3)
         sil         , & ! silicate (mmol/m^3)
         dmsp        , & ! dmsp (mmol/m^3)
         dms         , & ! dms (mmol/m^3)
         hum             ! humic material (mmol/m^3)

      real (kind=dbl_kind), dimension (max_algae), intent(in) :: &
         algalN          ! ocean algal nitrogen (mmol/m^3) (diatoms, phaeo, pico)

      real (kind=dbl_kind), dimension (max_doc), intent(in) :: &
         doc             ! ocean doc (mmol/m^3)  (proteins, EPS, lipid)

      real (kind=dbl_kind), dimension (max_don), intent(in) :: &
         don             ! ocean don (mmol/m^3) 

      real (kind=dbl_kind), dimension (max_dic), intent(in) :: &
         dic             ! ocean dic (mmol/m^3) 

      real (kind=dbl_kind), dimension (max_fe), intent(in) :: &
         fed, fep        ! ocean disolved and particulate fe (nM) 

      real (kind=dbl_kind), dimension (max_aero), intent(in) :: &
         zaeros          ! ocean aerosols (mmol/m^3) 

      real (kind=dbl_kind), dimension (max_nbtrcr), intent(inout) :: &
         ocean_bio_all   ! fixed order, all values even for tracers false

      ! local variables

      integer (kind=int_kind) :: &
         k, ks           ! tracer indices

      ocean_bio_all(:) = c0

      do k = 1, max_algae           
         ocean_bio_all(k)      = algalN(k)           ! N
         ks = max_algae + max_doc + max_dic + 1
         ocean_bio_all(ks + k) = R_chl2N(k)*algalN(k)!chl
      enddo   

      ks = max_algae + 1
      do k = 1, max_doc
         ocean_bio_all(ks + k) = doc(k)              ! doc
      enddo  
      ks = ks + max_doc
      do k = 1, max_dic
         ocean_bio_all(ks + k) = dic(k)              ! dic
      enddo 

      ks = 2*max_algae + max_doc + max_dic + 7
      do k = 1, max_don
         ocean_bio_all(ks + k) = don(k)              ! don
      enddo  

      ks = max_algae + 1
      ocean_bio_all(ks) = nit                        ! nit

      ks = 2*max_algae + max_doc + 2 + max_dic
      ocean_bio_all(ks) = amm                        ! Am
      ks = ks + 1
      ocean_bio_all(ks) = sil                        ! Sil
      ks = ks + 1
      ocean_bio_all(ks) =  R_S2N(1)*algalN(1) &      ! DMSPp
                        +  R_S2N(2)*algalN(2) &
                        +  R_S2N(3)*algalN(3) 
      ks = ks + 1
      ocean_bio_all(ks) = dmsp                       ! DMSPd
      ks = ks + 1
      ocean_bio_all(ks) = dms                        ! DMS
      ks = ks + 1
      ocean_bio_all(ks) = nit                        ! PON
      ks = 2*max_algae + max_doc + 7 + max_dic + max_don
      do k = 1, max_fe
         ocean_bio_all(ks + k) = fed(k)              ! fed
      enddo  
      ks = ks + max_fe
      do k = 1, max_fe
         ocean_bio_all(ks + k) = fep(k)              ! fep
      enddo  
      ks = ks + max_fe
      do k = 1, max_aero
         ocean_bio_all(ks+k) = zaeros(k)             ! zaero
      enddo
      ks = ks + max_aero + 1 
      ocean_bio_all(ks)  = hum                       ! humics

      end subroutine colpkg_init_OceanConcArray_double 
!=======================================================================
! Warning messages
!=======================================================================

      subroutine colpkg_clear_warnings()

        use ice_warnings, only: reset_warnings

        call reset_warnings()

      end subroutine colpkg_clear_warnings

!=======================================================================
      
      subroutine colpkg_get_warnings(warningsOut)

        use ice_warnings, only: &
             get_number_warnings, &
             get_warning

        character(len=char_len_long), dimension(:), allocatable, intent(out) :: &
             warningsOut

        integer :: &
             iWarning, &
             nWarnings

        nWarnings = get_number_warnings()

        if (allocated(warningsOut)) deallocate(warningsOut)
        allocate(warningsOut(nWarnings))

        do iWarning = 1, nWarnings
           warningsOut(iWarning) = trim(get_warning(iWarning))
        enddo

      end subroutine colpkg_get_warnings

!=======================================================================

      subroutine colpkg_print_warnings(nu_diag)

        use ice_warnings, only: &
             get_number_warnings, &
             get_warning

        integer, intent(in) :: nu_diag

        integer :: &
             iWarning

        do iWarning = 1, get_number_warnings()
           write(nu_diag,*) trim(get_warning(iWarning))
        enddo

      end subroutine colpkg_print_warnings

!=======================================================================

      end module ice_colpkg

!=======================================================================

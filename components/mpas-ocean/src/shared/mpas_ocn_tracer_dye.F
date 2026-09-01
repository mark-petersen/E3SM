! copyright (c) 2013,  los alamos national security, llc (lans)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_dye
!
!> \brief MPAS ocean dye
!> \author Mathew Maltrud
!> \date   10/27/2021
!> \details
!>  This module contains routines for computing evolution of dyes
!
!-----------------------------------------------------------------------

module ocn_tracer_dye

   use mpas_timer
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_forcing
   use ocn_framework_forcing
   use mpas_timekeeping
   use ocn_constants
   use ocn_config
   use mpas_stream_manager

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_tracer_dye_compute, &
             ocn_tracer_dye_surface_flux_compute,  &
             ocn_tracer_sed_surface_flux_compute,  &
             ocn_tracer_dye_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_dye_compute
!
!> \brief   computes a tracer tendency due to dye
!> \author  Mathew Maltrud
!> \date    12/22/2016
!> \details
!>  This routine computes a tracer tendency for dye
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_dye_compute!{{{

! no dye sources/sinks

   end subroutine ocn_tracer_dye_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_dye_surface_flux_compute
!
!> \brief   computes a tracer tendency due to dye surface fluxes
!> \author  Mathew Maltrud
!> \date    12/22/2016
!> \details
!>  This routine computes a tracer tendency due to dye surface fluxes
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_dye_surface_flux_compute(dyeTracers, forcingPool, meshPool, nTracers,   &
      nCellsSolve, minLevelCell, maxLevelCell, nVertLevels, layerThickness, &
      dyeSurfaceFlux, dyeSurfaceFluxRunoff, &
      indexdyeVolMackenzie, indexdyeVolOb, indexdyeVolYenesei, indexdyeVolLena, indexdyeVolKolyma, &
      indexdyeVolPrecipEvap, indexdyeVolIceMeltForm, &
      indexdyeVolMackenzieNS, indexdyeVolObNS, indexdyeVolYeneseiNS, indexdyeVolLenaNS, &
      indexdyeVolYukonNS, indexdyeVolKolymaNS, &
      indexdyeVolPrecipEvapNS, indexdyeVolIceMeltFormNS, &
      err)!{{{
      

 use ocn_config

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         dyeTracers

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness

      integer, dimension(:), intent(in) :: &
         minLevelCell, maxLevelCell

      ! scalars
      integer, intent(in) :: nTracers, nCellsSolve, nVertLevels, &
         indexdyeVolMackenzie,  &
         indexdyeVolOb,  &
         indexdyeVolYenesei,  &
         indexdyeVolLena,  &
         indexdyeVolKolyma,  &
         indexdyeVolPrecipEvap,  &
         indexdyeVolIceMeltForm,  &
         indexdyeVolMackenzieNS,  &
         indexdyeVolObNS,  &
         indexdyeVolYeneseiNS,  &
         indexdyeVolLenaNS,  &
         indexdyeVolYukonNS,  &
         indexdyeVolKolymaNS,  &
         indexdyeVolPrecipEvapNS,  &
         indexdyeVolIceMeltFormNS

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      ! two dimensional arrays
      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         dyeSurfaceFlux, dyeSurfaceFluxRunoff

      type (mpas_pool_type), intent(inout) :: forcingPool
      type (mpas_pool_type), intent(in) :: meshPool

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, iTracer, k

      real (kind=RKIND), dimension(:), pointer :: seaIceFreshWaterFlux
      real (kind=RKIND), dimension(:), pointer :: riverRunoffFlux, iceRunoffFlux,  &
         rainFlux, snowFlux, evaporationFlux

      type (mpas_pool_type), pointer :: dyeDiagFields

      real (kind=RKIND), dimension(:), pointer :: &
         dyeVolMackenzie_zint,  &
         dyeMackenzieMask,      &
         dyeVolOb_zint,  &
         dyeObMask,      &
         dyeVolYenesei_zint,  &
         dyeYeneseiMask,      &
         dyeVolLena_zint,  &
         dyeLenaMask,  &
         dyeVolKolyma_zint,  &
         dyeKolymaMask,  &
         dyeVolPrecipEvap_zint,  &
         dyeVolIceMeltForm_zint,  &

         dyeVolMackenzieNS_zint,  &
         dyeVolObNS_zint,  &
         dyeVolYeneseiNS_zint,  &
         dyeVolLenaNS_zint,  &
         dyeVolYukonNS_zint,  &
         dyeYukonMask,  &
         dyeVolKolymaNS_zint,  &
         dyeVolPrecipEvapNS_zint,  &
         dyeVolIceMeltFormNS_zint,  &
         dyeArcticMask

      real (kind=RKIND) :: totalRunoff
      real (kind=RKIND) :: totalNonRunoff

      err = 0

      call mpas_timer_start("dye surface flux")

      if (.not. config_use_dyeTracers) return

! get runoff fluxes
      call mpas_pool_get_array(forcingPool, 'riverRunoffFlux', riverRunoffFlux)
      call mpas_pool_get_array(forcingPool, 'iceRunoffFlux', iceRunoffFlux)
! get other freshwater fluxes
      call mpas_pool_get_array(forcingPool, 'rainFlux', rainFlux)
      call mpas_pool_get_array(forcingPool, 'snowFlux', snowFlux)
      call mpas_pool_get_array(forcingPool, 'evaporationFlux', evaporationFlux)
      call mpas_pool_get_array(forcingPool, 'seaIceFreshWaterFlux', seaIceFreshWaterFlux)

! get masks and diagnostic arrays
      call mpas_pool_get_subpool(forcingPool, 'dyeDiagFields', dyeDiagFields)
      call mpas_pool_get_array(dyeDiagFields, 'dyeVolMackenzie_zint', dyeVolMackenzie_zint)
      call mpas_pool_get_array(meshPool, 'dyeMackenzieMask', dyeMackenzieMask)
      call mpas_pool_get_array(dyeDiagFields, 'dyeVolOb_zint', dyeVolOb_zint)
      call mpas_pool_get_array(meshPool, 'dyeObMask', dyeObMask)
      call mpas_pool_get_array(dyeDiagFields, 'dyeVolYenesei_zint', dyeVolYenesei_zint)
      call mpas_pool_get_array(meshPool, 'dyeYeneseiMask', dyeYeneseiMask)
      call mpas_pool_get_array(dyeDiagFields, 'dyeVolLena_zint', dyeVolLena_zint)
      call mpas_pool_get_array(meshPool, 'dyeLenaMask', dyeLenaMask)
      call mpas_pool_get_array(dyeDiagFields, 'dyeVolKolyma_zint', dyeVolKolyma_zint)
      call mpas_pool_get_array(meshPool, 'dyeKolymaMask', dyeKolymaMask)
      call mpas_pool_get_array(dyeDiagFields, 'dyeVolPrecipEvap_zint', dyeVolPrecipEvap_zint)
      call mpas_pool_get_array(dyeDiagFields, 'dyeVolIceMeltForm_zint', dyeVolIceMeltForm_zint)
      call mpas_pool_get_array(meshPool, 'dyeArcticMask', dyeArcticMask)

      call mpas_pool_get_array(dyeDiagFields, 'dyeVolMackenzieNS_zint', dyeVolMackenzieNS_zint)
      call mpas_pool_get_array(dyeDiagFields, 'dyeVolObNS_zint', dyeVolObNS_zint)
      call mpas_pool_get_array(dyeDiagFields, 'dyeVolYeneseiNS_zint', dyeVolYeneseiNS_zint)
      call mpas_pool_get_array(dyeDiagFields, 'dyeVolLenaNS_zint', dyeVolLenaNS_zint)
      call mpas_pool_get_array(dyeDiagFields, 'dyeVolYukonNS_zint', dyeVolYukonNS_zint)
      call mpas_pool_get_array(meshPool, 'dyeYukonMask', dyeYukonMask)
      call mpas_pool_get_array(dyeDiagFields, 'dyeVolKolymaNS_zint', dyeVolKolymaNS_zint)
      call mpas_pool_get_array(dyeDiagFields, 'dyeVolPrecipEvapNS_zint', dyeVolPrecipEvapNS_zint)
      call mpas_pool_get_array(dyeDiagFields, 'dyeVolIceMeltFormNS_zint', dyeVolIceMeltFormNS_zint)

      dyeVolMackenzie_zint = 0.0_RKIND
      dyeVolOb_zint = 0.0_RKIND
      dyeVolYenesei_zint = 0.0_RKIND
      dyeVolLena_zint = 0.0_RKIND
      dyeVolKolyma_zint = 0.0_RKIND
      dyeVolPrecipEvap_zint = 0.0_RKIND
      dyeVolIceMeltForm_zint = 0.0_RKIND

      dyeVolMackenzieNS_zint = 0.0_RKIND
      dyeVolObNS_zint = 0.0_RKIND
      dyeVolYeneseiNS_zint = 0.0_RKIND
      dyeVolLenaNS_zint = 0.0_RKIND
      dyeVolYukonNS_zint = 0.0_RKIND
      dyeVolKolymaNS_zint = 0.0_RKIND
      dyeVolPrecipEvapNS_zint = 0.0_RKIND
      dyeVolIceMeltFormNS_zint = 0.0_RKIND

      !DWJ 08/05/2016: This loop needs OpenMP added to it.
      do iCell=1,nCellsSolve

! only use positive values for flux
      totalRunoff = max( (riverRunoffFlux(iCell) + iceRunoffFlux(iCell)), 0.0_RKIND)

! convert to mg/m2/s since using kg results in huge numbers
        if (dyeMackenzieMask(iCell) > 0.0_RKIND)   &
          dyeSurfaceFluxRunoff(indexdyeVolMackenzie, iCell) = totalRunoff*1.e-6_RKIND
        if (dyeObMask(iCell) > 0.0_RKIND)   &
          dyeSurfaceFluxRunoff(indexdyeVolOb, iCell) = totalRunoff*1.e-6_RKIND
        if (dyeYeneseiMask(iCell) > 0.0_RKIND)   &
          dyeSurfaceFluxRunoff(indexdyeVolYenesei, iCell) = totalRunoff*1.e-6_RKIND
        if (dyeLenaMask(iCell) > 0.0_RKIND)   &
          dyeSurfaceFluxRunoff(indexdyeVolLena, iCell) = totalRunoff*1.e-6_RKIND
        if (dyeKolymaMask(iCell) > 0.0_RKIND)   &
          dyeSurfaceFluxRunoff(indexdyeVolKolyma, iCell) = totalRunoff*1.e-6_RKIND

        if (dyeMackenzieMask(iCell) > 0.0_RKIND)   &
          dyeSurfaceFluxRunoff(indexdyeVolMackenzieNS, iCell) = totalRunoff*1.e-6_RKIND
        if (dyeObMask(iCell) > 0.0_RKIND)   &
          dyeSurfaceFluxRunoff(indexdyeVolObNS, iCell) = totalRunoff*1.e-6_RKIND
        if (dyeYeneseiMask(iCell) > 0.0_RKIND)   &
          dyeSurfaceFluxRunoff(indexdyeVolYeneseiNS, iCell) = totalRunoff*1.e-6_RKIND
        if (dyeLenaMask(iCell) > 0.0_RKIND)   &
          dyeSurfaceFluxRunoff(indexdyeVolLenaNS, iCell) = totalRunoff*1.e-6_RKIND
        if (dyeYukonMask(iCell) > 0.0_RKIND)   &
          dyeSurfaceFluxRunoff(indexdyeVolYukonNS, iCell) = totalRunoff*1.e-6_RKIND
        if (dyeKolymaMask(iCell) > 0.0_RKIND)   &
          dyeSurfaceFluxRunoff(indexdyeVolKolymaNS, iCell) = totalRunoff*1.e-6_RKIND

!maltrud debug
!if (iCell == 1) write(*,*)' DEBUG02: ',iCell, dyeArcticMask(iCell), indexdyeVolPrecipEvap,indexdyeVolIceMeltForm
        if (dyeArcticMask(iCell) > 0.0_RKIND) then
          dyeSurfaceFluxRunoff(indexdyeVolPrecipEvap, iCell) =  &
             rainFlux(iCell) + snowFlux(iCell) + evaporationFlux(iCell)
          dyeSurfaceFluxRunoff(indexdyeVolPrecipEvapNS, iCell) =  &
             rainFlux(iCell) + snowFlux(iCell) + evaporationFlux(iCell)
          dyeSurfaceFluxRunoff(indexdyeVolIceMeltForm, iCell) = seaIceFreshWaterFlux(iCell)
        endif
! entire globe source for dyeVolIceMeltFormNS
! COMMENT OUT AFTER SOURCES ARE TURNED OFF
!        dyeSurfaceFluxRunoff(indexdyeVolIceMeltFormNS, iCell) = seaIceFreshWaterFlux(iCell)

! do vertical integrals
        do k = minLevelCell(iCell), maxLevelCell(iCell)
          dyeVolMackenzie_zint(iCell) = dyeVolMackenzie_zint(iCell) &
                            + layerThickness(k,iCell)*dyeTracers(indexdyeVolMackenzie, k, iCell)
          dyeVolOb_zint(iCell) = dyeVolOb_zint(iCell) &
                            + layerThickness(k,iCell)*dyeTracers(indexdyeVolOb, k, iCell)
          dyeVolYenesei_zint(iCell) = dyeVolYenesei_zint(iCell) &
                            + layerThickness(k,iCell)*dyeTracers(indexdyeVolYenesei, k, iCell)
          dyeVolLena_zint(iCell) = dyeVolLena_zint(iCell) &
                            + layerThickness(k,iCell)*dyeTracers(indexdyeVolLena, k, iCell)
          dyeVolKolyma_zint(iCell) = dyeVolKolyma_zint(iCell) &
                            + layerThickness(k,iCell)*dyeTracers(indexdyeVolKolyma, k, iCell)
          dyeVolPrecipEvap_zint(iCell) = dyeVolPrecipEvap_zint(iCell) &
                            + layerThickness(k,iCell)*dyeTracers(indexdyeVolPrecipEvap, k, iCell)
          dyeVolIceMeltForm_zint(iCell) = dyeVolIceMeltForm_zint(iCell) &
                            + layerThickness(k,iCell)*dyeTracers(indexdyeVolIceMeltForm, k, iCell)
          dyeVolMackenzieNS_zint(iCell) = dyeVolMackenzieNS_zint(iCell) &
                            + layerThickness(k,iCell)*dyeTracers(indexdyeVolMackenzieNS, k, iCell)
          dyeVolObNS_zint(iCell) = dyeVolObNS_zint(iCell) &
                            + layerThickness(k,iCell)*dyeTracers(indexdyeVolObNS, k, iCell)
          dyeVolYeneseiNS_zint(iCell) = dyeVolYeneseiNS_zint(iCell) &
                            + layerThickness(k,iCell)*dyeTracers(indexdyeVolYeneseiNS, k, iCell)
          dyeVolLenaNS_zint(iCell) = dyeVolLenaNS_zint(iCell) &
                            + layerThickness(k,iCell)*dyeTracers(indexdyeVolLenaNS, k, iCell)
          dyeVolYukonNS_zint(iCell) = dyeVolYukonNS_zint(iCell) &
                            + layerThickness(k,iCell)*dyeTracers(indexdyeVolYukonNS, k, iCell)
          dyeVolKolymaNS_zint(iCell) = dyeVolKolymaNS_zint(iCell) &
                            + layerThickness(k,iCell)*dyeTracers(indexdyeVolKolymaNS, k, iCell)
          dyeVolPrecipEvapNS_zint(iCell) = dyeVolPrecipEvapNS_zint(iCell) &
                            + layerThickness(k,iCell)*dyeTracers(indexdyeVolPrecipEvapNS, k, iCell)
          dyeVolIceMeltFormNS_zint(iCell) = dyeVolIceMeltFormNS_zint(iCell) &
                            + layerThickness(k,iCell)*dyeTracers(indexdyeVolIceMeltFormNS, k, iCell)
        enddo  !  k

! maltrud now do virtual flux losses
! ASSUME all fluxes in top layer (usually ok for non-runoff)
        k = 1

      totalNonRunoff = rainFlux(iCell) + snowFlux(iCell) + &
                       evaporationFlux(iCell) + seaIceFreshWaterFlux(iCell)

! only use negative values for flux
      dyeSurfaceFlux(:, iCell) = 0.0_RKIND
      if (totalNonRunoff < 0.0_RKIND) dyeSurfaceFlux(:, iCell) =   &
!       dyeSurfaceFlux(:, iCell) + totalNonRunoff*dyeTracers(:, k, iCell)/rho_fw
        totalNonRunoff*dyeTracers(:, k, iCell)/rho_fw

      enddo  !  iCell

      call mpas_timer_stop("dye surface flux")

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_dye_surface_flux_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_sed_surface_flux_compute
!
! \brief   computes a tracer tendency due to sediment surface fluxes
!> \author  Mark Petersen
!> \date    08/22/2026
!> \details
!>  This routine computes a tracer tendency due to sediment surface fluxes
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_sed_surface_flux_compute(sedTracers, forcingPool, meshPool, nTracers,   &
      nCellsSolve, minLevelCell, maxLevelCell, nVertLevels, layerThickness, &
      sedSurfaceFlux, sedSurfaceFluxRunoff, &
      indexSedMackenzie, indexSedOb, indexSedYenesei, indexSedLena, &
      indexSedYukon, indexSedKolyma, &
      err)!{{{
      

 use ocn_config

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         sedTracers

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness

      integer, dimension(:), intent(in) :: &
         minLevelCell, maxLevelCell

      ! scalars
      integer, intent(in) :: nTracers, nCellsSolve, nVertLevels, &
         indexSedMackenzie,  &
         indexSedOb,  &
         indexSedYenesei,  &
         indexSedLena,  &
         indexSedYukon,  &
         indexSedKolyma

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      ! two dimensional arrays
      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         sedSurfaceFlux, sedSurfaceFluxRunoff

      type (mpas_pool_type), intent(inout) :: forcingPool
      type (mpas_pool_type), intent(in) :: meshPool

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, iTracer, k

      real (kind=RKIND), dimension(:), pointer :: seaIceFreshWaterFlux
      real (kind=RKIND), dimension(:), pointer :: riverRunoffFlux, iceRunoffFlux,  &
         rainFlux, snowFlux, evaporationFlux


      real (kind=RKIND), dimension(:), pointer :: &
         dyeMackenzieMask,      &
         dyeObMask,      &
         dyeYeneseiMask,      &
         dyeLenaMask,  &
         dyeKolymaMask,  &
         dyeYukonMask

      real (kind=RKIND) :: totalRunoff
      real (kind=RKIND) :: totalNonRunoff

      err = 0

      call mpas_timer_start("sed surface flux")

! get runoff fluxes
      call mpas_pool_get_array(forcingPool, 'riverRunoffFlux', riverRunoffFlux)
      call mpas_pool_get_array(forcingPool, 'iceRunoffFlux', iceRunoffFlux)
! get other freshwater fluxes
      call mpas_pool_get_array(forcingPool, 'rainFlux', rainFlux)
      call mpas_pool_get_array(forcingPool, 'snowFlux', snowFlux)
      call mpas_pool_get_array(forcingPool, 'evaporationFlux', evaporationFlux)
      call mpas_pool_get_array(forcingPool, 'seaIceFreshWaterFlux', seaIceFreshWaterFlux)

! get masks and diagnostic arrays
      call mpas_pool_get_array(meshPool, 'dyeMackenzieMask', dyeMackenzieMask)
      call mpas_pool_get_array(meshPool, 'dyeObMask', dyeObMask)
      call mpas_pool_get_array(meshPool, 'dyeYeneseiMask', dyeYeneseiMask)
      call mpas_pool_get_array(meshPool, 'dyeLenaMask', dyeLenaMask)
      call mpas_pool_get_array(meshPool, 'dyeYukonMask', dyeYukonMask)
      call mpas_pool_get_array(meshPool, 'dyeKolymaMask', dyeKolymaMask)

      do iCell=1,nCellsSolve

! only use positive values for flux
      totalRunoff = max( (riverRunoffFlux(iCell) + iceRunoffFlux(iCell)), 0.0_RKIND)

! convert to mg/m2/s since using kg results in huge numbers
        if (dyeMackenzieMask(iCell) > 0.0_RKIND)   &
          sedSurfaceFluxRunoff(indexSedMackenzie, iCell) = totalRunoff*1.e-6_RKIND
        if (dyeObMask(iCell) > 0.0_RKIND)   &
          sedSurfaceFluxRunoff(indexSedOb, iCell) = totalRunoff*1.e-6_RKIND
        if (dyeYeneseiMask(iCell) > 0.0_RKIND)   &
          sedSurfaceFluxRunoff(indexSedYenesei, iCell) = totalRunoff*1.e-6_RKIND
        if (dyeLenaMask(iCell) > 0.0_RKIND)   &
          sedSurfaceFluxRunoff(indexSedLena, iCell) = totalRunoff*1.e-6_RKIND
        if (dyeYukonMask(iCell) > 0.0_RKIND)   &
          sedSurfaceFluxRunoff(indexSedYukon, iCell) = totalRunoff*1.e-6_RKIND
        if (dyeKolymaMask(iCell) > 0.0_RKIND)   &
          sedSurfaceFluxRunoff(indexSedKolyma, iCell) = totalRunoff*1.e-6_RKIND

      enddo  !  iCell

      call mpas_timer_stop("sed surface flux")

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_sed_surface_flux_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_dye_init
!
!> \brief   Initializes ocean dyes
!> \author  Mathew Maltrud
!> \date    12/28/2016
!> \details
!>  This routine initializes fields required for dye tracers
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_dye_init(domain,err)!{{{

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information

      integer, intent(out) :: err !< Output: error flag

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      err = 0

      if (.not. config_use_dyeTracers) return

!maltrud no init necessary

!     call mpas_stream_mgr_read(domain % streamManager, streamID='dyeMasks', ierr=err)
!maltrud debug
!if(err /= 0) write(*,*)' DEBUG03: ', err

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_dye_init!}}}

!***********************************************************************

end module ocn_tracer_dye

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker

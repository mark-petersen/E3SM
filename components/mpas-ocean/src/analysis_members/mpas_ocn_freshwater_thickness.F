! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_freshwater_thickness
!
!> \brief MPAS ocean analysis mode member: Computes freshwater thickness
!>          standard layers
!> \author Dylan Schlichting
!> \date   February 06, 2025
!> \details
!>   Computes freshwater thickness by depth-integrating river tracers
!>   Details described in Zhang et al. 2012 JGR Oceans
!>   https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2012JC008108
!>   Returns a nCells array for each river.
!>
!-----------------------------------------------------------------------

module ocn_freshwater_thickness

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use ocn_constants
   use ocn_diagnostics_variables

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

   public :: ocn_init_freshwater_thickness, &
             ocn_compute_freshwater_thickness, &
             ocn_restart_freshwater_thickness, &
             ocn_finalize_freshwater_thickness

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_freshwater_thickness
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author Dylan Schlichting
!> \date    February 06, 2025
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_freshwater_thickness(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type(domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_init_freshwater_thickness!}}}

!***********************************************************************
!
!  routine ocn_compute_freshwater_thickness
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Dylan Schlichting
!> \date    February 06, 2025
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_freshwater_thickness(domain, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type(domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type(mpas_pool_type), pointer :: freshwaterThicknessAMPool
      type(dm_info) :: dminfo
      type(block_type), pointer :: block
      type(mpas_pool_type), pointer :: statePool
      type(mpas_pool_type), pointer :: meshPool
      type(mpas_pool_type), pointer :: tracersPool

      ! Here are some example variables which may be needed for your analysis member
      integer, pointer :: nVertLevels, nCellsSolve, nEdgesSolve, nVerticesSolve
      integer :: iTracer, k, iCell
      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop, maxLevelVertexBot

      real(kind=RKIND), dimension(:, :, :), pointer :: tracers, riverTagTracers
      real(kind=RKIND), dimension(:, :), pointer :: riverTagTracer1RiverRunoff, riverTagTracer2RiverRunoff
      real(kind=RKIND), dimension(:, :), pointer :: layerThickness
      real(kind=RKIND), dimension(:), pointer ::  riverTagTracer1FwThickness, riverTagTracer2FwThickness

      real(kind=RKIND) :: freshwaterLayer_river1, freshwaterLayer_river2
      err = 0

      dminfo = domain%dminfo

      block => domain%blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block%structs, 'state', statePool)
         call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block%structs, 'freshwaterThickness', freshwaterThicknessAMPool)
         call mpas_pool_get_dimension(block%dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block%dimensions, 'nCellsSolve', nCellsSolve)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
         call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)
         call mpas_pool_get_array(tracersPool, 'riverTagTracers', riverTagTracers, timeLevel)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
         call mpas_pool_get_array(freshwaterThicknessAMPool, 'riverTagTracer1FwThickness', riverTagTracer1FwThickness) 
         call mpas_pool_get_array(freshwaterThicknessAMPool, 'riverTagTracer2FwThickness', riverTagTracer2FwThickness) 

         ! Depth integrate the river tracer concentration to get freshwater thickness. 
         do iCell = 1, nCellsSolve
            k = 1
            riverTagTracer1FwThickness(iCell) = 0.0_RKIND
            riverTagTracer2FwThickness(iCell) = 0.0_RKIND
            do while (k .le. maxLevelCell(iCell))
               freshwaterLayer_river1 = layerThickness(k, iCell)*riverTagTracers(1, k, iCell)
               freshwaterLayer_river2 = layerThickness(k, iCell)*riverTagTracers(2, k, iCell)
               riverTagTracer1FwThickness(iCell) = riverTagTracer1FwThickness(iCell) + freshwaterLayer_river1
               riverTagTracer2FwThickness(iCell) = riverTagTracer2FwThickness(iCell) + freshwaterLayer_river2
               k = k + 1
            enddo
         end do

         block => block%next
      end do

   end subroutine ocn_compute_freshwater_thickness!}}}

!***********************************************************************
!
!  routine ocn_restart_freshwater_thickness
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Dylan Schlichting
!> \date    February 06, 2025
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_freshwater_thickness(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type(domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_restart_freshwater_thickness!}}}

!***********************************************************************
!
!  routine ocn_finalize_freshwater_thickness
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Dylan Schlichting
!> \date    February 06, 2025
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_freshwater_thickness(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type(domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_finalize_freshwater_thickness!}}}

end module ocn_freshwater_thickness

! vim: foldmethod=marker

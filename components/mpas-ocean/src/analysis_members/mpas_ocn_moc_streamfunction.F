! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_moc_streamfunction
!
!> \brief MPAS ocean analysis mode member: moc_streamfunction
!> \author  Nils H. Feige, Mark R. Petersen
!> \date    2016-04-08
!> \brief   Computes Meridional Overturning Circulation streamfunction.
!>
!-----------------------------------------------------------------------

module ocn_moc_streamfunction

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager
   use mpas_spline_interpolation

   use ocn_constants
   use ocn_config
   use ocn_diagnostics_variables
   use ocn_mesh

   !use ppr_1d

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------
   !type(rmap_opts) :: opts
   !integer :: bc_upper, bc_lower

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_init_moc_streamfunction, &
             ocn_compute_moc_streamfunction, &
             ocn_restart_moc_streamfunction, &
             ocn_finalize_moc_streamfunction

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer :: nMocStreamfunctionBinsUsed

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_moc_streamfunction
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Nils H. Feige, Mark R. Petersen
!> \date    2016-04-08
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_moc_streamfunction(domain, err)!{{{

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

      type (domain_type), intent(inout) :: domain

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

      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: mocStreamfunctionAMPool

      integer ::  iBin, k

      real (kind=RKIND) :: binWidth
      ! These are array size 1 because mpas_dmpar_min_real_array calls require arrays.
      real (kind=RKIND), dimension(1) :: minBin, maxBin, minBinDomain, maxBinDomain
      ! the variable used to discriminate cells into Bins (either the y-value or the latitude)
      real (kind=RKIND), dimension(:), pointer :: binBoundaryMocStreamfunction
      real (kind=RKIND), dimension(:), pointer :: mocInterpBottomDepth, mocInterpZMid, mocInterpLayerThickness

      !determines if the simulation was run on a sphere or on a plane
      logical, pointer :: on_a_sphere

      !!!! REGION STUFF
      !! region moc calculation variables
      integer :: currentRegion, i, iCell

      !! region arrays/variables
      character (len=StrKIND), dimension(:), pointer :: regionGroupNames
      integer, dimension(:, :), pointer :: regionCellMasks, regionsInGroup
      integer, dimension(:), pointer ::  nRegionsInGroup
      integer, pointer :: nRegions, nRegionGroups, maxRegionsInGroup
      real (kind=RKIND), dimension(:,:), pointer :: minMaxLatRegion
      real (kind=RKIND), dimension(:)  , pointer :: minLatRegionLocal, maxLatRegionLocal, tminLatRegionLocal, tmaxLatRegionLocal

      !! region preliminary variables
      integer :: regionGroupNumber, regionsInAddGroup

      !!region pool
      type (mpas_pool_type), pointer :: regionPool

      !! region dimensions
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nRegions', nRegions)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nRegionGroups', nRegionGroups)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'maxRegionsInGroup', maxRegionsInGroup)

      !! get region values
      call mpas_pool_get_subpool(domain % blocklist % structs, 'regions', regionPool)
      call mpas_pool_get_array(regionPool, 'regionsInGroup', regionsInGroup)
      call mpas_pool_get_array(regionPool, 'nRegionsInGroup', nRegionsInGroup)
      call mpas_pool_get_array(regionPool, 'regionGroupNames', regionGroupNames)

      !!! region preliminaries
      regionGroupNumber=0
      do i = 1, nRegionGroups
         if (regionGroupNames(i) .eq. config_AM_mocStreamfunction_region_group) then
            regionGroupNumber = i
         end if
      end do
      if (regionGroupNumber.eq.0) then
         call mpas_log_write &
            ('ocn_moc_streamfunction AM: config_AM_mocStreamfunction_region_group not found in regionGroupNames')
         call mpas_log_write &
            ('ocn_moc_streamfunction AM: setting regionGroupNumber=1')
         regionGroupNumber = 1
      end if

      regionsInAddGroup = MIN(nRegionsInGroup(regionGroupNumber), maxRegionsInGroup)
      !!!! END REGION STUFF

      allocate(minLatRegionLocal(maxRegionsInGroup))
      allocate(maxLatRegionLocal(maxRegionsInGroup))
      allocate(tminLatRegionLocal(maxRegionsInGroup))
      allocate(tmaxLatRegionLocal(maxRegionsInGroup))

      dminfo = domain % dminfo

      err = 0

      minBin =  1.0e34_RKIND
      maxBin = -1.0e34_RKIND

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mocStreamfunctionAM', mocStreamfunctionAMPool)

      call mpas_pool_get_array(mocStreamfunctionAMPool, 'minMaxLatRegion', minMaxLatRegion)

      minLatRegionLocal(:) =  4.0_RKIND
      maxLatRegionLocal(:) = -4.0_RKIND

      nMocStreamfunctionBinsUsed = config_AM_mocStreamfunction_num_bins

      call mpas_pool_get_array(mocStreamfunctionAMPool, 'binBoundaryMocStreamfunction', binBoundaryMocStreamfunction)

      ! Find min and max values of binning variable. For the whole domain as well as for each region
      ! in the current region group.
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_array(regionPool, 'regionCellMasks', regionCellMasks)

         minBin = min(minBin, minval(latCell) )
         maxBin = max(maxBin, maxval(latCell) )

         do i = 1, regionsInAddGroup
            currentRegion = regionsInGroup(i, regionGroupNumber)
            do iCell = 1, nCellsOwned
               if (regionCellMasks(currentRegion, iCell) .eq. 1) then
                  minLatRegionLocal(i) = min(minLatRegionLocal(i), latCell(iCell))
                  maxLatRegionLocal(i) = max(maxLatRegionLocal(i), latCell(iCell))
               end if
            end do
         end do

         block => block % next
      end do

      call mpas_dmpar_min_real_array(dminfo, 1, minBin, minBinDomain)
      call mpas_dmpar_max_real_array(dminfo, 1, maxBin, maxBinDomain)

      call mpas_dmpar_min_real_array(dminfo, maxRegionsInGroup, minLatRegionLocal(:), tminLatRegionLocal(:))
      call mpas_dmpar_max_real_array(dminfo, maxRegionsInGroup, maxLatRegionLocal(:), tmaxLatRegionLocal(:))
      minMaxLatRegion(1, :) = tminLatRegionLocal(:)
      minMaxLatRegion(2, :) = tmaxLatRegionLocal(:)

      deallocate(minLatRegionLocal)
      deallocate(maxLatRegionLocal)
      deallocate(tminLatRegionLocal)
      deallocate(tmaxLatRegionLocal)

      ! Set up bins.
      binBoundaryMocStreamfunction = -1.0e34_RKIND

      ! Change min and max bin bounds to configuration settings, if applicable.
      if (config_AM_mocStreamfunction_min_bin > -1.0e33_RKIND) then
         minBinDomain(1) = config_AM_mocStreamfunction_min_bin
      else
         ! use measured min value, but decrease slightly to include least value.
         minBinDomain(1) = minBinDomain(1) - 1.0e-10_RKIND * abs(minBinDomain(1))
      end if

      if (config_AM_mocStreamfunction_max_bin > -1.0e33_RKIND) then
         maxBinDomain(1) = config_AM_mocStreamfunction_max_bin
      else
         ! use measured max value, but increase slightly to include max value.
         maxBinDomain(1) = maxBinDomain(1) + 1.0e-10_RKIND * abs(maxBinDomain(1))
      end if

      binBoundaryMocStreamfunction(1) = minBinDomain(1)
      binWidth = (maxBinDomain(1) - minBinDomain(1)) / nMocStreamfunctionBinsUsed

      ! Use the same bin boundaries for the regions and the global MOC.
      do iBin = 2, nMocStreamfunctionBinsUsed + 1
         binBoundaryMocStreamfunction(iBin) = binBoundaryMocStreamfunction(iBin-1) + binWidth
      end do

      if (config_AM_mocStreamfunction_vert_interp) then
         ! Options for remapping
         !if ( config_vert_remap_order == 1 ) then
         !   opts%cell_meth = pcm_method ! PCM method
         !elseif ( config_vert_remap_order == 2 ) then
         !   opts%cell_meth = plm_method ! PLM method
         !elseif ( config_vert_remap_order == 3 ) then
         !   opts%cell_meth = ppm_method ! PPM method
         !   opts%edge_meth = p3e_method ! 3rd-order edge interp.
         !elseif ( config_vert_remap_order == 5 ) then
         !   opts%cell_meth = pqm_method ! PPM method
         !   opts%edge_meth = p5e_method ! 3rd-order edge interp.
         !endif
         !if ( config_limiter == 'monotonic' ) then
         !   opts%cell_lims = mono_limit ! monotone slope limits
         !elseif ( config_limiter == 'weno' ) then
         !   opts%cell_lims = weno_limit ! WENO slope limits
         !   opts%wall_lims = mono_limit ! WENO slope limits
         !else
         !   opts%cell_lims = null_limit ! no slope limits
         !endif
         !bc_upper = bcon_loose
         !bc_lower = bcon_loose
      end if

   end subroutine ocn_init_moc_streamfunction!}}}

!***********************************************************************
!
!  routine ocn_compute_moc_streamfunction
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Nils H. Feige, Mark R. Petersen
!> \date    2016-04-08
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_moc_streamfunction(domain, timeLevel, err)!{{{
      implicit none
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

      type (domain_type), intent(inout) :: domain

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

      type (mpas_pool_type), pointer :: mocStreamfunctionAMPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: regionPool
      type (mpas_pool_type), pointer :: transectPool

      integer :: iCell, iBin, genericCounter, k, kGC, kZL
      real (kind=RKIND) :: binWidth
      real (kind=RKIND), dimension(:,:), pointer :: mocStreamValLatAndDepthLocal
      real (kind=RKIND), dimension(:), pointer ::  binBoundaryMocStreamfunction
      real (kind=RKIND), dimension(:,:), pointer :: mocStreamvalLatAndDepth, mocStreamValLatAndDepthTotal
      real (kind=RKIND), dimension(:,:), pointer :: mocStreamvalLatAndDepthGM
      real (kind=RKIND), dimension(:,:), pointer :: sumVertBinVelocity

      !!!! TRANSECT VARIABLES !!!!
      integer, pointer :: num_tracers
      integer :: iEdge, iTransect, c1, c2, currentTransect
      integer, dimension(:,:), pointer :: transectEdgeMasks, transectEdgeMaskSigns, &
            transectsInGroup

      real (kind=RKIND) :: m3ps_to_Sv
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, normalVelocity
      real (kind=RKIND), dimension(:,:), allocatable ::  sumTransport, totalSumTransport
      character (len=STRKIND), dimension(:), pointer :: transectNames, transectGroupNames
      integer, dimension(:), pointer ::  nTransectsInGroup
      integer, pointer :: nTransects, nTransectGroups, maxTransectsInGroup
      integer :: transectGroupNumber, transectsInAddGroup
      !!!! END TRANSECT VARIABLES !!!!

      !!!! REGION VARIABLES
      real (kind=RKIND) :: maskFactor
      integer :: currentRegion, i
      real (kind=RKIND), dimension(:,:,:), pointer :: mocStreamValLatAndDepthRegionLocal, &
                         mocStreamvalLatAndDepthRegion, mocStreamValLatAndDepthRegionTotal, &
                         sumVertBinVelocityRegion, mocStreamvalLatAndDepthRegionGM
      integer, dimension(:, :), pointer :: regionCellMasks, regionVertexMasks, regionsInGroup
      character (len=STRKIND), dimension(:), pointer :: regionNames, regionGroupNames
      integer, dimension(:), pointer ::  nRegionsInGroup
      integer, pointer :: nRegions, nRegionGroups, maxRegionsInGroup
      integer :: regionGroupNumber, regionsInAddGroup
      !!!! END REGION VARIABLES
    
      ! interpolation variables
      real (kind=RKIND) :: bottomDepthEdge, dzUpper, dzLower
      real (kind=RKIND), dimension(:), allocatable :: zMidEdge, zTopEdge, refZMid, transportZL, layerThicknessEdge
      real (kind=RKIND), dimension(:), allocatable :: refZTop, vertVelocityTopZL, xIn, xOut, yIn
      real (kind=RKIND), dimension(:,:), allocatable :: normalVelocityHZLAll
      integer :: nIn, nOut
      real(kind=RKIND) :: invAreaCell1, edgeNorm
      real (kind=RKIND), dimension(:), allocatable:: div_hu,div_huTransport,div_huGMBolus

      !!!! REGION INITIALIZATION
      !! region dimensions
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nRegions', nRegions)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nRegionGroups', nRegionGroups)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'maxRegionsInGroup', maxRegionsInGroup)

      !! get region values
      call mpas_pool_get_subpool(domain % blocklist % structs, 'regions', regionPool)
      call mpas_pool_get_array(regionPool, 'regionsInGroup', regionsInGroup)
      call mpas_pool_get_array(regionPool, 'nRegionsInGroup', nRegionsInGroup)
      call mpas_pool_get_array(regionPool, 'regionNames', regionNames)
      call mpas_pool_get_array(regionPool, 'regionGroupNames', regionGroupNames)

      regionGroupNumber=0
      do i = 1, nRegionGroups
         if (regionGroupNames(i) .eq. config_AM_mocStreamfunction_region_group) then
            regionGroupNumber = i
         end if
      end do
      if (regionGroupNumber.eq.0) then
         call mpas_log_write &
            ('ocn_moc_streamfunction AM: config_AM_mocStreamfunction_region_group not found in regionGroupNames')
         call mpas_log_write &
            ('ocn_moc_streamfunction AM: setting regionGroupNumber=1')
         regionGroupNumber = 1
      end if

      regionsInAddGroup = MIN(nRegionsInGroup(regionGroupNumber), maxRegionsInGroup)

      !! allocate regional moc calculation arrays
      allocate(mocStreamValLatAndDepthRegionLocal(nMocStreamfunctionBinsUsed + 1, nVertLevels, maxRegionsInGroup))
      allocate(sumVertBinVelocityRegion(nMocStreamfunctionBinsUsed + 1, nVertLevels, maxRegionsInGroup))
      allocate(mocStreamValLatAndDepthRegionTotal(nMocStreamfunctionBinsUsed + 1, nVertLevels, maxRegionsInGroup))
      allocate(zMidEdge(nVertLevels), zTopEdge(nVertLevelsP1), refZMid(nVertLevels), transportZL(nVertLevels), layerThicknessEdge(nVertLevels))
      allocate(normalVelocityHZLAll(nVertLevels, nEdgesAll))
      allocate(refZTop(nVertLevelsP1), vertVelocityTopZL(nVertLevelsP1))
      allocate(xIn(nVertLevels+2), xOut(nVertLevelsP1), yIn(nVertLevels+2))
      allocate(div_hu(nVertLevels),div_huTransport(nVertLevels),div_huGMBolus(nVertLevels))

      ! compute refZTop. This could be moved to the mesh later.
      refZTop(1) = 0.0_RKIND
      do k = 1,nVertLevels
         refZTop(k+1) = refZTop(k) - refLayerThickness(k)
      end do
      xOut = -refZTop
      nOut = nVertLevels + 1

      mocStreamValLatAndDepthRegionLocal = 0.0_RKIND
      mocStreamValLatAndDepthRegionTotal = 0.0_RKIND
      sumVertBinVelocityRegion = 0.0_RKIND
      !!!! END REGION INITIALIZATION

      !!!! TRANSECT INITIALIZATION
      !! transect dimensions
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nTransects', nTransects)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nTransectGroups', nTransectGroups)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'maxTransectsInGroup', maxTransectsInGroup)

      !! get transect values
      call mpas_pool_get_subpool(domain % blocklist % structs, 'transects', transectPool)
      call mpas_pool_get_array(transectPool, 'transectsInGroup', transectsInGroup)
      call mpas_pool_get_array(transectPool, 'nTransectsInGroup', nTransectsInGroup)
      call mpas_pool_get_array(transectPool, 'transectNames', transectNames)
      call mpas_pool_get_array(transectPool, 'transectGroupNames', transectGroupNames)

      transectGroupNumber=0
      do i = 1, nTransectGroups
         if (transectGroupNames(i) .eq. config_AM_mocStreamfunction_transect_group) then
            transectGroupNumber = i
         end if
      end do
      if (transectGroupNumber.eq.0) then
         call mpas_log_write &
            ('ocn_moc_streamfunction AM: config_AM_mocStreamfunction_transect_group not found in transectGroupNames')
         call mpas_log_write &
            ('ocn_moc_streamfunction AM: setting transectGroupNumber=1')
         transectGroupNumber = 1
      end if

      transectsInAddGroup = nTransectsInGroup(transectGroupNumber)

      allocate(sumTransport(nVertLevels,maxTransectsInGroup))
      allocate(totalSumTransport(nVertLevels,maxTransectsInGroup))

      m3ps_to_Sv = 1e-6
      !!!! END TRANSECT INITIALIZATION

      if (transectsInAddGroup .ne. regionsInAddGroup) then
         ! This writes output every step, and is too verbose.
         !call mpas_log_write ('ocn_moc_streamfunction AM: transectsInGroup count does not ' &
         !   // 'match regionsInGroup count: $i, $i', intArgs = (/ transectsInAddGroup, regionsInAddGroup /) )
         i = min(transectsInAddGroup, regionsInAddGroup)
         transectsInAddGroup = i
         regionsInAddGroup = i
         !call mpas_log_write ('Setting both to min: $i, $i', intArgs = (/ transectsInAddGroup, regionsInAddGroup /) )
      end if

      err = 0

      dminfo = domain % dminfo

      allocate(mocStreamValLatAndDepthLocal(nMocStreamfunctionBinsUsed + 1, nVertLevels))
      allocate(sumVertBinVelocity(nMocStreamfunctionBinsUsed + 1, nVertLevels))
      allocate(mocStreamValLatAndDepthTotal(nMocStreamfunctionBinsUsed + 1, nVertLevels))

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mocStreamfunctionAM', mocStreamfunctionAMPool)

         call mpas_pool_get_array(mocStreamfunctionAMPool, 'binBoundaryMocStreamfunction', binBoundaryMocStreamfunction)

         binWidth = (binBoundaryMocStreamfunction(nMocStreamfunctionBinsUsed + 1) - binBoundaryMocStreamfunction(1)) &
         / nMocStreamfunctionBinsUsed

         call mpas_pool_get_array(regionPool, 'regionCellMasks', regionCellMasks)

         !!!! TRANSECT DOMAINSPLIT VARIABLES
         call mpas_pool_get_array(transectPool,'transectEdgeMaskSigns',transectEdgeMaskSigns)
         call mpas_pool_get_array(transectPool,'transectEdgeMasks',transectEdgeMasks)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)

! example of calling rmap1d
!         if (config_AM_mocStreamfunction_vert_interp) then
!            !$omp parallel
!            !$omp do schedule(runtime) private(iEdge, kTop, kBot, nLayers, nLevels, dataNow, dataNew)
!            do iEdge = 1, nEdges
!
!               kTop = minLevelEdgeBot(iEdge)
!               kBot = maxLevelEdgeTop(iEdge)
!
!               nLayers = kBot - kTop + 1
!               nLevels = nLayers + 1
!
!               if (nLevels .lt. 2) cycle
!
!               dataNow(1, 1, 1:nLayers) = normalVelocity(kTop:kBot, iEdge)
!
!               call rmap1d(nLevels, nLevels, nVars ,nDoFs, &
!                  heightEdgeNow(kTop:kBot+1, iEdge), &
!                  heightEdgeNew(kTop:kBot+1, iEdge), &
!                  dataNow, dataNew, &
!                  bcUpper, bcLower, work, opts)

!    !subroutine rmap1d(npos,nnew,nvar,ndof,xpos, &
!    !    &             xnew,fdat,fnew,bclo,bcup, &
!    !    &             work,opts,tCPU)
!    !
!    ! NPOS  no. edges in old grid.
!    ! NNEW  no. edges in new grid.
!    ! NVAR  no. discrete variables to remap.
!    ! NDOF  no. degrees-of-freedom per cell.
!    ! XPOS  old grid edge positions. XPOS is a length NPOS
!    !       array .
!    ! XNEW  new grid edge positions. XNEW is a length NNEW
!    !       array .
!    ! FDAT  grid-cell moments on old grid. FNEW has SIZE = 
!    !       NDOF-by-NVAR-by-NNEW-1 .
!    ! FNEW  grid-cell moments on new grid. FNEW has SIZE = 
!    !       NDOF-by-NVAR-by-NNEW-1 .
!    ! BCLO  boundary condition at lower endpoint .
!    ! BCHI  boundary condition at upper endpoint . 
!    ! WORK  method work-space. See RCON-WORK for details .
!    ! OPTS  method parameters. See RCON-OPTS for details .
!    ! TCPU  method tcpu-timer.
!    !
!
!               normalVelocity   (kTop:kBot, iEdge) = dataNew(1, 1, 1:nLayers)
!            
!            end do
!            !$omp end do
!            !$omp end parallel
!         end if


            !  Compute normalVelocityHZLAll for all edges
            if (config_AM_mocStreamfunction_vert_interp) then
               ! vertical interpolation, for general vertical coordinates
               do iEdge = 1,nEdgesAll
                  normalVelocityHZLAll(:,iEdge) = 0.0_RKIND
                  if (maxLevelEdgeTop(iEdge).le.1) cycle
                  c1 = cellsOnEdge(1,iEdge)
                  c2 = cellsOnEdge(2,iEdge)
                  bottomDepthEdge = 0.5_RKIND*(bottomDepth(c1) + bottomDepth(c2))
                  do k = 1, maxLevelEdgeTop(iEdge)
                      layerThicknessEdge(k) = 0.5_RKIND*(layerThickness(k,c1) + layerThickness(k,c2))
                  end do
                  zTopEdge(1) = 0.0_RKIND
                  zMidEdge(1) = -0.5_RKIND*layerThicknessEdge(1)
                  do k = 2, maxLevelEdgeTop(iEdge)
                     zTopEdge(k) = zTopEdge(k-1) - layerThicknessEdge(k-1)
                     zMidEdge(k) = zMidEdge(k-1) - 0.5_RKIND*(layerThicknessEdge(k-1) + layerThicknessEdge(k))
                  end do
                  k = maxLevelEdgeTop(iEdge) + 1
                  zTopEdge(k) = zTopEdge(k-1) - layerThicknessEdge(k-1)
                  ! Integrate normalVelocity:
                  ! (zTopEdge, normalVelocity) -> (refZTop, normalVelocityHZLAll)
                  kZL = 1
                  ! loop through general coordinate levels
                  do kGC = 1,maxLevelEdgeTop(iEdge)
                     if (zTopEdge(kGC+1) > refZTop(kZL+1) - 1.0e-8) then
                        ! if GC layer is full within z-level, give the whole thing
                        normalVelocityHZLAll(kZL,iEdge) = normalVelocityHZLAll(kZL,iEdge) + normalVelocity(kGC,iEdge) * layerThicknessEdge(kGC)
                     else
                        ! if GC layer is partially within z-level, give just part.
                        normalVelocityHZLAll(kZL,iEdge) = normalVelocityHZLAll(kZL,iEdge) + normalVelocity(kGC,iEdge) * (zTopEdge(kGC) - refZTop(kZL+1))
                        ! and give the rest to the next level
                        kZL = min(nVertLevels, kZL + 1)
                        normalVelocityHZLAll(kZL,iEdge) = normalVelocityHZLAll(kZL,iEdge) + normalVelocity(kGC,iEdge) * (refZTop(kZL) - zTopEdge(kGC+1))
                     endif
                  end do
               end do
            endif

         !!!! TRANSECT CALCULATION
         sumTransport = 0.0_RKIND
         mocStreamValLatAndDepthLocal = 0.0_RKIND
         mocStreamValLatAndDepthTotal = 0.0_RKIND
         sumVertBinVelocity = 0.0_RKIND
         mocStreamValLatAndDepthRegionLocal = 0.0_RKIND
         mocStreamValLatAndDepthRegionTotal = 0.0_RKIND
         sumVertBinVelocityRegion = 0.0_RKIND

         do iTransect = 1,transectsInAddGroup
            currentTransect = transectsInGroup(iTransect, transectGroupNumber)
            if (config_AM_mocStreamfunction_vert_interp) then
               ! vertical interpolation, for general vertical coordinates
               do iEdge = 1,nEdgesOwned
                if (transectEdgeMasks(currentTransect,iEdge).eq.1.and.maxLevelEdgeTop(iEdge).gt.0) then
                  c1 = cellsOnEdge(1,iEdge)
                  c2 = cellsOnEdge(2,iEdge)
                  bottomDepthEdge = 0.5_RKIND*(bottomDepth(c1) + bottomDepth(c2))
                  do k = 1, maxLevelEdgeTop(iEdge)
                      layerThicknessEdge(k) = 0.5_RKIND*(layerThickness(k,c1) + layerThickness(k,c2))
                  end do
                  zTopEdge(1) = 0.0_RKIND
                  zMidEdge(1) = -0.5_RKIND*layerThicknessEdge(1)
                  do k = 2, maxLevelEdgeTop(iEdge)
                     zTopEdge(k) = zTopEdge(k-1) - layerThicknessEdge(k-1)
                     zMidEdge(k) = zMidEdge(k-1) - 0.5_RKIND*(layerThicknessEdge(k-1) + layerThicknessEdge(k))
                  end do
                  k = maxLevelEdgeTop(iEdge) + 1
                  zTopEdge(k) = zTopEdge(k-1) - layerThicknessEdge(k-1)
                  ! Integrate normalVelocity:
                  ! (zTopEdge, normalVelocity) -> (refZTop, transportZL)
                  transportZL(:) = 0.0_RKIND
                  kZL = 1
                  ! loop through general coordinate levels
                  do kGC = 1,maxLevelEdgeTop(iEdge)
                     if (zTopEdge(kGC+1) > refZTop(kZL+1) - 1.0e-8) then
                        ! if GC layer is full within z-level, give the whole thing
                        transportZL(kZL) = transportZL(kZL) + normalVelocity(kGC,iEdge) * layerThicknessEdge(kGC)
                     else
                        ! if GC layer is partially within z-level, give just part.
                        transportZL(kZL) = transportZL(kZL) + normalVelocity(kGC,iEdge) * (zTopEdge(kGC) - refZTop(kZL+1))
                        ! and give the rest to the next level
                        kZL = min(nVertLevels, kZL + 1)
                        transportZL(kZL) = transportZL(kZL) + normalVelocity(kGC,iEdge) * (refZTop(kZL) - zTopEdge(kGC+1))
                     endif
                  end do
                  do k = 1, nVertLevels
                     sumTransport(k,iTransect) = sumTransport(k,iTransect) + &
                        transectEdgeMaskSigns(currentTransect,iEdge) &
                        * transportZL(k)*dvEdge(iEdge)
                  end do
                 endif !(transectEdgeMasks(currentTransect, iEdge) .eq. 1)
               end do
            else
               ! no vertical interpolation, for z-level runs
               do iEdge = 1,nEdgesOwned
                  c1 = cellsOnEdge(1,iEdge)
                  c2 = cellsOnEdge(2,iEdge)
                  do k = 1, maxLevelEdgeTop(iEdge)
                     sumTransport(k,iTransect) = sumTransport(k,iTransect) + &
                        transectEdgeMaskSigns(currentTransect,iEdge) &
                        * transectEdgeMasks(currentTransect, iEdge) &
                        * normalVelocity(k,iEdge)*dvEdge(iEdge) &
                        * 0.5_RKIND*(layerThickness(k,c1) + layerThickness(k,c2))
                  end do
               end do
            end if
            do k = 2, nVertLevels
               mocStreamValLatAndDepthRegionLocal(1, k, iTransect) = &
                  mocStreamValLatAndDepthRegionLocal(1, k - 1, iTransect) &
                  + sumTransport(k - 1, iTransect)
            end do
         end do

         !!!! END TRANSECT CALCULATION

          if (config_AM_mocStreamfunction_vert_interp) then
            ! vertical interpolation, for general vertical coordinates
            do iCell = 1,nCellsOwned
               iBin = max(int((latCell(iCell) - binBoundaryMocStreamfunction(1)) / binWidth), 2)

!!!!!!!! This is linear interpolation of vertVelocityTop
!               ! Interpolate vertVelocityTop:
!               ! (zTop(:,iCell), vertVelocityTop(:,iCell)) -> (refZTop(:), vertVelocityTopZL(:))
!               !  xIn            yIn size: nIn=maxLevelCell+2  xOut        yOut  nVertLevels
!
!               nIn = maxLevelCell(iCell) + 2
!               xIn(1:nIn-2) = -zTop(1:nIn-2,iCell)
!               yIn(1:nIn-2) = vertVelocityTop(1:nIn-2,iCell)
!               ! Add two extra points so all velocities are zero below bottom
!               xIn(nIn-1) = bottomDepth(iCell)
!               yIn(nIn-1) = 0.0_RKIND
!               xIn(nIn) = 6000 ! deepest ocean cell in m
!               yIn(nIn) = 0.0_RKIND
!               ! xOut = -refZTop (set above)
!               ! nOut = nVertLevels + 1 (set above)
!               xIn(1) = 0.0_RKIND ! set SSH to match xOut
!               call mpas_interpolate_linear(xIn,yIn,nIn, &
!                  xOut,vertVelocityTopZL,nOut)
               !if (iCell < 10000) then
               !    print *, 'iCell, nIn, nOut',iCell, nIn, nOut
               !print *, 'xIn',xIn
               !print *, 'yIn',yIn
               !print *, 'xOut',xOut
               !print *, 'yOut',vertVelocityTopZL
               !endif
!!!!!!!! End linear interpolation of vertVelocityTop

               div_hu(:) = 0.0_RKIND
               !div_huTransport(:) = 0.0_RKIND
               !div_huGMBolus(:) = 0.0_RKIND
               invAreaCell1 = invAreaCell(iCell)
               do i = 1, nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(i, iCell)
                  edgeNorm = dvEdge(iEdge) * invAreaCell1 * edgeSignOnCell(i, iCell)
                  do k = 1,nVertLevels
                     div_hu(k) = div_hu(k) - normalVelocityHZLAll(k, iEdge) * edgeNorm
                     !div_huTransport(k) = div_huTransport(k) - normalTransportVelocity(k, iEdge) * edgeNorm * refLayerThickness(k)
                     !div_huGMBolus(k) = div_huGMBolus(k) - normalGMBolusVelocity(k, iEdge) * edgeNorm * refLayerThickness(k)
                  end do
               end do
               ! Vertical velocity at bottom (maxLevelCell(iCell)+1) is zero, initialized above.
               vertVelocityTopZL(1) = 0.0_RKIND
               vertVelocityTopZL(nVertLevels+1) = 0.0_RKIND
               do k = nVertLevels, 2, -1
                  vertVelocityTopZL(k) = vertVelocityTopZL(k+1) - div_hu(k)
                  !vertTransportVelocityTopZL(k) = vertTransportVelocityTopZL(k+1) - div_huTransport(k)
                  !vertGMBolusVelocityTopZL(k) = vertGMBolusVelocityTopZL(k+1) - div_huGMBolus(k)
               end do
               do k = 1,nVertLevels
                  do i = 1, regionsInAddGroup
                     currentRegion = regionsInGroup(i, regionGroupNumber)
                     sumVertBinVelocityRegion(iBin, k, i) = sumVertBinVelocityRegion(iBin, k, i) + vertVelocityTopZL(k) * &
                           areaCell(iCell) * regionCellMasks(currentRegion, iCell)
                  end do
                  sumVertBinVelocity(iBin, k) = sumVertBinVelocity(iBin, k) + (vertVelocityTopZL(k) * areaCell(iCell))
               end do
            end do

         else ! no vertical interpolation, for z-level runs
          do iCell = 1,nCellsOwned
            iBin = MAX(int((latCell(iCell) - binBoundaryMocStreamfunction(1)) / binWidth), 2)
            do k = 1,maxLevelCell(iCell)
               do i = 1, regionsInAddGroup
                  currentRegion = regionsInGroup(i, regionGroupNumber)
                  sumVertBinVelocityRegion(iBin, k, i) = sumVertBinVelocityRegion(iBin, k, i) + (vertVelocityTop(k, iCell) * &
                        areaCell(iCell) * regionCellMasks(currentRegion, iCell))
               end do
               sumVertBinVelocity(iBin, k) = sumVertBinVelocity(iBin, k) + (vertVelocityTop(k, iCell) * areaCell(iCell))
            end do
          end do
         endif

         do i = 1, regionsInAddGroup
            do k = 1,nVertLevels
               do iBin = 2, nMocStreamfunctionBinsUsed + 1
                  mocStreamValLatAndDepthLocal(iBin, k) = mocStreamValLatAndDepthLocal(iBin-1, k) &
                     + sumVertBinVelocity(iBin, k)
                  mocStreamValLatAndDepthRegionLocal(iBin, k, i) = mocStreamValLatAndDepthRegionLocal(iBin-1, k, i) &
                     + sumVertBinVelocityRegion(iBin, k, i)
               end do
            end do
         end do

         block => block % next
     end do

     call mpas_dmpar_sum_real_array(dminfo, nVertLevels * (nMocStreamfunctionBinsUsed + 1), mocStreamValLatAndDepthLocal, &
         mocStreamvalLatAndDepthTotal)

     call mpas_dmpar_sum_real_array(dminfo, nVertLevels * (nMocStreamfunctionBinsUsed + 1) * maxRegionsInGroup, &
         mocStreamValLatAndDepthRegionLocal, mocStreamvalLatAndDepthRegionTotal)

     call mpas_pool_get_subpool(domain % blocklist % structs, 'mocStreamfunctionAM', mocStreamfunctionAMPool)
     call mpas_pool_get_array(mocStreamfunctionAMPool, 'mocStreamvalLatAndDepth', mocStreamvalLatAndDepth)
     mocStreamvalLatAndDepth = mocStreamvalLatAndDepthTotal * m3ps_to_Sv

     call mpas_pool_get_array(mocStreamfunctionAMPool, 'mocStreamvalLatAndDepthRegion', mocStreamvalLatAndDepthRegion)
     mocStreamvalLatAndDepthRegion = mocStreamvalLatAndDepthRegionTotal * m3ps_to_Sv

     if(config_use_GM) THEN !compute GM bolus contribution to the streamfunction
       block => domain % blocklist
       do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'mocStreamfunctionAM', mocStreamfunctionAMPool)

           call mpas_pool_get_array(mocStreamfunctionAMPool, 'binBoundaryMocStreamfunction', binBoundaryMocStreamfunction)

           binWidth = (binBoundaryMocStreamfunction(nMocStreamfunctionBinsUsed + 1) - binBoundaryMocStreamfunction(1)) &
           / nMocStreamfunctionBinsUsed

           call mpas_pool_get_array(regionPool, 'regionCellMasks', regionCellMasks)

           !!!! TRANSECT DOMAINSPLIT VARIABLES
           call mpas_pool_get_array(transectPool,'transectEdgeMaskSigns',transectEdgeMaskSigns)
           call mpas_pool_get_array(transectPool,'transectEdgeMasks',transectEdgeMasks)
           normalVelocity => normalGMBolusVelocity
           call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)

           !!!! TRANSECT CALCULATION
           sumTransport = 0.0_RKIND
           mocStreamValLatAndDepthLocal = 0.0_RKIND
           mocStreamValLatAndDepthTotal = 0.0_RKIND
           sumVertBinVelocity = 0.0_RKIND
           mocStreamValLatAndDepthRegionLocal = 0.0_RKIND
           mocStreamValLatAndDepthRegionTotal = 0.0_RKIND
           sumVertBinVelocityRegion = 0.0_RKIND

           do iTransect = 1,transectsInAddGroup
              currentTransect = transectsInGroup(iTransect, transectGroupNumber)
              do iEdge = 1,nEdgesOwned
                 c1 = cellsOnEdge(1,iEdge)
                 c2 = cellsOnEdge(2,iEdge)
                 do k = 1, maxLevelEdgeTop(iEdge)
                    sumTransport(k,iTransect) = sumTransport(k,iTransect) + &
                       transectEdgeMaskSigns(currentTransect,iEdge) &
                       * transectEdgeMasks(currentTransect, iEdge) &
                       * normalVelocity(k,iEdge)*dvEdge(iEdge) &
                       * 0.5_RKIND*(layerThickness(k,c1) + layerThickness(k,c2))
                 end do
              end do
              do k = 2, nVertLevels
                 mocStreamValLatAndDepthRegionLocal(1, k, iTransect) = &
                    mocStreamValLatAndDepthRegionLocal(1, k - 1, iTransect) &
                    + sumTransport(k - 1, iTransect)
              end do
           end do

           !!!! END TRANSECT CALCULATION

           do iCell = 1,nCellsOwned
              iBin = max(int((latCell(iCell) - binBoundaryMocStreamfunction(1)) / binWidth), 2)
              do k = 1,maxLevelCell(iCell)
                 do i = 1, regionsInAddGroup
                    currentRegion = regionsInGroup(i, regionGroupNumber)
                    sumVertBinVelocityRegion(iBin, k, i) = sumVertBinVelocityRegion(iBin, k, i) + (vertGMBolusVelocityTop(k, iCell) * &
                          areaCell(iCell) * regionCellMasks(currentRegion, iCell))
                 end do
                 sumVertBinVelocity(iBin, k) = sumVertBinVelocity(iBin, k) + (vertGMBolusVelocityTop(k, iCell) * areaCell(iCell))
              end do
           end do

           do i = 1, regionsInAddGroup
              do k = 1,nVertLevels
                 do iBin = 2, nMocStreamfunctionBinsUsed + 1
                    mocStreamValLatAndDepthLocal(iBin, k) = mocStreamValLatAndDepthLocal(iBin-1, k) &
                       + sumVertBinVelocity(iBin, k)
                    mocStreamValLatAndDepthRegionLocal(iBin, k, i) = mocStreamValLatAndDepthRegionLocal(iBin-1, k, i) &
                       + sumVertBinVelocityRegion(iBin, k, i)
                 end do
              end do
           end do

           block => block % next
       end do

       call mpas_dmpar_sum_real_array(dminfo, nVertLevels * (nMocStreamfunctionBinsUsed + 1), mocStreamValLatAndDepthLocal, &
           mocStreamvalLatAndDepthTotal)

       call mpas_dmpar_sum_real_array(dminfo, nVertLevels * (nMocStreamfunctionBinsUsed + 1) * maxRegionsInGroup, &
           mocStreamValLatAndDepthRegionLocal, mocStreamvalLatAndDepthRegionTotal)

       call mpas_pool_get_subpool(domain % blocklist % structs, 'mocStreamfunctionAM', mocStreamfunctionAMPool)
       call mpas_pool_get_array(mocStreamfunctionAMPool, 'mocStreamvalLatAndDepthGM', mocStreamvalLatAndDepthGM)
       mocStreamvalLatAndDepthGM = mocStreamvalLatAndDepthTotal * m3ps_to_Sv

       call mpas_pool_get_array(mocStreamfunctionAMPool, 'mocStreamvalLatAndDepthRegionGM', mocStreamvalLatAndDepthRegionGM)
       mocStreamvalLatAndDepthRegionGM = mocStreamvalLatAndDepthRegionTotal * m3ps_to_Sv

       !Add GM bolus contribution to resolved streamfunction to create total streamfunction
       mocStreamvalLatAndDepthRegion = mocStreamvalLatAndDepthRegion + mocStreamvalLatAndDepthRegionGM
       mocStreamvalLatAndDepth = mocStreamvalLatAndDepth + mocStreamvalLatAndDepthGM

     endif !config_use_GM

     deallocate(mocStreamvalLatAndDepthTotal)
     deallocate(mocStreamvalLatAndDepthLocal)
     deallocate(sumVertBinVelocity)

     deallocate(mocStreamvalLatAndDepthRegionTotal)
     deallocate(mocStreamvalLatAndDepthRegionLocal)
     deallocate(sumVertBinVelocityRegion)
     deallocate(zMidEdge, refZMid, transportZL, layerThicknessEdge)
     deallocate(refZTop, vertVelocityTopZL)
     deallocate(normalVelocityHZLAll)
     deallocate(div_hu,div_huTransport,div_huGMBolus)

     !!!! TRANSECT CELANUP
     deallocate(sumTransport)
     deallocate(totalSumTransport)

   end subroutine ocn_compute_moc_streamfunction!}}}

!***********************************************************************
!
!  routine ocn_restart_moc_streamfunction
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Nils H. Feige, Mark R. Petersen
!> \date    2016-04-08
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_moc_streamfunction(domain, err)!{{{

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

      type (domain_type), intent(inout) :: domain

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

   end subroutine ocn_restart_moc_streamfunction!}}}

!***********************************************************************
!
!  routine ocn_finalize_moc_streamfunction
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Nils H. Feige, Mark R. Petersen
!> \date    2016-04-08
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_moc_streamfunction(domain, err)!{{{

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

      type (domain_type), intent(inout) :: domain

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

   end subroutine ocn_finalize_moc_streamfunction!}}}

end module ocn_moc_streamfunction

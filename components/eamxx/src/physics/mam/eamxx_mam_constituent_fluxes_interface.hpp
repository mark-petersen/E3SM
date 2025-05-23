#ifndef EAMXX_MAM_CONSTITUENT_FLUXES_HPP
#define EAMXX_MAM_CONSTITUENT_FLUXES_HPP

// For declaring contituent fluxes class derived from atm process
// class
#include <physics/mam/eamxx_mam_generic_process_interface.hpp>

// For MAM4 aerosol configuration
#include <physics/mam/mam_coupling.hpp>
#include <string>

namespace scream {

// The process responsible for applying MAM4 constituent fluxes. The
// AD stores exactly ONE instance of this class in its list of subcomponents.
class MAMConstituentFluxes final : public MAMGenericInterface {
 public:
  using KT            = ekat::KokkosTypes<DefaultDevice>;
  using const_view_2d = Field::view_dev_t<const Real **>;

 private:
  const_view_2d constituent_fluxes_;
  // aerosol state variables
  mam_coupling::AerosolState wet_aero_, dry_aero_;
  // wet mixing ratios (water species)
  mam_coupling::WetAtmosphere wet_atm_;
  // dry mixing ratios (water species)
  mam_coupling::DryAtmosphere dry_atm_;
  // workspace manager for internal local variables
  mam_coupling::Buffer buffer_;

 public:
  // Constructor
  MAMConstituentFluxes(const ekat::Comm &comm,
                       const ekat::ParameterList &params);

  // --------------------------------------------------------------------------
  // AtmosphereProcess overrides (see share/atm_process/atmosphere_process.hpp)
  // --------------------------------------------------------------------------

  // The name of the subcomponent
  std::string name() const { return "mam_constituent_fluxes"; }

  // grid
  void set_grids(
      const std::shared_ptr<const GridsManager> grids_manager) override;

  // management of common atm process memory
  size_t requested_buffer_size_in_bytes() const override;
  void init_buffers(const ATMBufferManager &buffer_manager) override;

  // Initialize variables
  void initialize_impl(const RunType run_type) override;

  // Run the process by one time step
  void run_impl(const double dt) override;

  // Finalize
  void finalize_impl(){/*Do nothing*/};

};  // MAMConstituentFluxes

}  // namespace scream

#endif  // EAMXX_MAM_CONSTITUENT_FLUXES_HPP

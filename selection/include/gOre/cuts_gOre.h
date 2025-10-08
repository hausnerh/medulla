/**
 * @file cuts_gOre.h
 * @brief Header file for defining cuts useful for NC analyses
 * specficially for saving photon/electron discrimination for the end
 * @author hhausner@fnal.gov
 **/

#ifndef CUTS_GORE_H
#define CUTS_GORE_H

#include "include/particle_variables.h"

#include "include/gOre/core_gOre.h"

/**
 * @namespace cuts::gOre
 * @brief Cut specific to NC single photon analyses
 **/
namespace cuts::gOre
{
  /**
   * @brief Is there a single gOre in the interaction?
   * @tparam T the type of interaction (true or reco).
   * @param obj the interaction in question
   * @param params a vector whose first element is the energy threshold
   * @return true if there is either a single photon xor single electron
   * @return false if there is not exactly one photon or exactly one electron
   **/
  template<class T>
    bool single_gOre(const T& obj, std::vector<double> params={GORE_MIN_GORE_ENERGY})
    {
      // 0 if none, 1 if either one photon or one electron, 2+ if multiple
      size_t gOre_multiplicity = particle_multiplicity(obj, 1, pvars::kPhoton,   params)
                               + particle_multiplicity(obj, 1, pvars::kElectron, params);
      size_t    leading_idx = selectors::gOre::   leading_primary_gOre(obj);
      size_t subleading_idx = selectors::gOre::subleading_primary_gOre(obj);
      bool    leading_abv_thresh = (   leading_idx != kNoMatch) && (pvars::ke(obj.particles.at(   leading_idx)) > params.at(0));
      bool subleading_abv_thresh = (subleading_idx != kNoMatch) && (pvars::ke(obj.particles.at(subleading_idx)) > params.at(0));
      if ((gOre_multiplicity == 1) && subleading_abv_thresh)
      {
        std::cout << "Found single shower event with subleading shower above threshold:\n"
                  << "  Photon Multiplicity: " << particle_multiplicity(obj, 1, pvars::kPhoton,   params) << '\n'
                  << "  Electron Multipicity: " << particle_multiplicity(obj, 1, pvars::kElectron, params) << '\n'
                  << "  Leading gOre:\n"
                  << "    pid: " << obj.particles.at(   leading_idx).pid << '\n'
                  << "    KE: " <<  pvars::ke(obj.particles.at(   leading_idx)) << '\n'
                  << "  Subeading gOre:\n"
                  << "    pid: " << obj.particles.at(subleading_idx).pid << '\n'
                  << "    KE: " <<  pvars::ke(obj.particles.at(subleading_idx)) << '\n'
                  << std::endl;
      }
      return (gOre_multiplicity == 1);
    }
  REGISTER_CUT_SCOPE(RegistrationScope::Both, single_gOre, single_gOre);

  /**
   * @brief Does the interaction have a single photon or single electron topology?
   * @details This cut makes use of the thresholds defined in core::gOre::Interaction.
   * Require there are no muons or pions above threshold, either a single photon or a single electron,
   * but no requirement on the protons. The photon and electrons are to be distinguished later with an
   * optimized PID cut.
   * @tparam T the type of interaction (true or reco).
   * @param obj the interaction in question
   * @return true if there is a single shower-like particle and no muons or pions, false otherwise.
   **/
  template<class T>
    bool gOre_topology(const T& obj, std::vector<double> params={GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      return cuts::gOre::single_gOre(obj, {params.at(0)}) &&
             cuts::no_muons         (obj, {params.at(1)}) &&
             cuts::no_charged_pions (obj, {params.at(2)}) ;
    }
  REGISTER_CUT_SCOPE(RegistrationScope::Both, gOre_topology, gOre_topology); 

  /**
   * @brief Does the gOre shower fall in our energy range of interest?
   * @details ∆(1232) resonance peak means the photon's energy in an Nγ decay
   * should be (M∆^2 - MN^2)/(2M∆) in the rest frame, approx 258 MeV (neutron vs proton matters little).
   * There is some spread around this value due to, eg, nuclear motion but it's a good target. Give a +/- 150 MeV spread
   * @tparam T the type of interaction (true or reco)
   * @param obj the interaction in question 
   * @return if the kinetic energy of the shower is in the range of 258 +/- 150 MeV
   **/
  template<class T>
    bool gOre_target_KE(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      double gOre_KE = interaction.primary_gOre()->ke;
      return (108 < gOre_KE) && (gOre_KE < 408);
    }
  REGISTER_CUT_SCOPE(RegistrationScope::Both, gOre_target_KE, gOre_target_KE);

  /**
   * @brief Does the interaction fall in our stricter fiducial requirements
   * @details The allow distance between the vertex and each detector wall.
   * These limits are defiend in core_gOre.h
   * @tparam T the type of interaction (true or reco)
   * @param obj the interaction in question 
   * @return true if the vertex is inside our fiducial volume
   **/
  template<class T>
    bool gOre_fiducial_cut(const T& obj)
    {
      auto [vtx_x, vtx_y, vtx_z] = utilities::to_three_vector(obj.vertex);
      return (GORE_WALL_X_POS - vtx_x > GORE_FID_THRESH_X_POS) &&
             (vtx_x - GORE_WALL_X_NEG > GORE_FID_THRESH_X_NEG) &&
             (GORE_WALL_Y_POS - vtx_y > GORE_FID_THRESH_Y_POS) &&
             (vtx_y - GORE_WALL_Y_NEG > GORE_FID_THRESH_Y_NEG) &&
             (GORE_WALL_Z_POS - vtx_z > GORE_FID_THRESH_Z_POS) &&
             (vtx_z - GORE_WALL_Z_NEG > GORE_FID_THRESH_Z_NEG) ;
    }
  REGISTER_CUT_SCOPE(RegistrationScope::Both, gOre_fiducial_cut, gOre_fiducial_cut);

  /** @brief gOre_topology & no protons **/
  template<class T>
    bool gOre_0p(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      return (interaction.is_valid) && (interaction.nProtons() == 0);
    }
  REGISTER_CUT_SCOPE(RegistrationScope::Both, gOre_0p, gOre_0p); 

  /** @brief gOre_topology & 1 proton **/
  template<class T>
    bool gOre_1p(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      return (interaction.is_valid) && (interaction.nProtons() == 1);
    }
  REGISTER_CUT_SCOPE(RegistrationScope::Both, gOre_1p, gOre_1p);

  /** @brief gOre_topology & at least 1 proton **/
  template<class T>
    bool gOre_Np(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      return (interaction.is_valid) && (interaction.nProtons() > 0);
    }
  REGISTER_CUT_SCOPE(RegistrationScope::Both, gOre_Np, gOre_Np); 

  /** @brief fiducial/contained gOre_0p with flash cut **/
  template<class T>
    bool fid_con_flash_gOre_0p(const T& obj)
    {
      return cuts::fiducial_cut(obj) && cuts::containment_cut(obj) && cuts::flash_cut(obj) && gOre_0p(obj);
    }
  REGISTER_CUT_SCOPE(RegistrationScope::Both, fid_con_flash_gOre_0p, fid_con_flash_gOre_0p);

  /** @brief fiducial/contained gOre_1p with flash cut **/
  template<class T>
    bool fid_con_flash_gOre_1p(const T& obj)
    {
      return cuts::fiducial_cut(obj) && cuts::containment_cut(obj) && cuts::flash_cut(obj) && gOre_1p(obj);
    }
  REGISTER_CUT_SCOPE(RegistrationScope::Both, fid_con_flash_gOre_1p, fid_con_flash_gOre_1p); 

  /** @brief fiducial/contained gOre_Np with flash cut **/
  template<class T>
    bool fid_con_flash_gOre_Np(const T& obj)
    {
      return cuts::fiducial_cut(obj) && cuts::containment_cut(obj) && cuts::flash_cut(obj) && gOre_Np(obj);
    }
  REGISTER_CUT_SCOPE(RegistrationScope::Both, fid_con_flash_gOre_Np, fid_con_flash_gOre_Np); 

  /**
   * @brief does the interaction contain a sub-threshold photon candidate?
   * @tparam T the type of interaction (true or reco).
   * @param obj the interaction in question
   * @return true if there is a valid sub-leading photon candidate
   **/
  template<class T>
    bool has_subthreshold_gOre(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      return (interaction.subleading_gore_ke > 0);
    }
  REGISTER_CUT_SCOPE(RegistrationScope::Both, has_subthreshold_gOre, has_subthreshold_gOre); 

  // *** TRUTH CUTS ***//

  /** @brief true NC interaction **/
  template<class T>
  bool isnc(const T& obj) { return obj.current_type == 1; }
  REGISTER_CUT_SCOPE(RegistrationScope::True, isnc, isnc);

  /** @brief true fiducial/contained neutrino **/
  template<class T>
  bool is_fid_con_nu(const T& obj)
  {
    return cuts::fiducial_cut(obj) && cuts::containment_cut(obj) && cuts::neutrino(obj);
  }
  REGISTER_CUT_SCOPE(RegistrationScope::True, is_fid_con_nu, is_fid_con_nu);

  /** @brief true fiducia/contained NC resonant event **/
  template<class T>
  bool is_fid_con_nc_res(const T& obj)
  {
    return is_fid_con_nu(obj) && isnc(obj) && (obj.interaction_mode == caf::kRes);
  }
  REGISTER_CUT_SCOPE(RegistrationScope::True, is_fid_con_nc_res, is_fid_con_nc_res);

  /** @brief true single photon event **/
  template<class T>
  bool gOre_is_photon(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
  {
    core::gOre::True_Interaction interaction(obj, params);
    return interaction.is_valid && (pvars::pid(*interaction.primary_gOre()) == pvars::kPhoton);
  }
  REGISTER_CUT_SCOPE(RegistrationScope::True, gOre_is_photon, gOre_is_photon);

  /** @brief true single electron event **/
  template<class T>
  bool gOre_is_electron(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
  {
    core::gOre::True_Interaction interaction(obj, params);
    return interaction.is_valid && (pvars::pid(*interaction.primary_gOre()) == pvars::kElectron);
  }
  REGISTER_CUT_SCOPE(RegistrationScope::True, gOre_is_electron, gOre_is_electron);

  /** @brief has muon above threshold **/
  template<class T>
  bool muon_abv_thresh(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
  {
    return (core::gOre::count_primaries(obj, params.at(0), params.at(1), params.at(2), params.at(3)).at(pvars::kMuon) > 0);
  }
  REGISTER_CUT_SCOPE(RegistrationScope::True, muon_abv_thresh, muon_abv_thresh);

  /** @brief has pion above threshold **/
  template<class T>
  bool pion_abv_thresh(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
  {
    return (core::gOre::count_primaries(obj, params.at(0), params.at(1), params.at(2), params.at(3)).at(pvars::kMuon) > 0);
  }
  REGISTER_CUT_SCOPE(RegistrationScope::True, pion_abv_thresh, pion_abv_thresh);

  /** @brief has electron above threshold **/
  template<class T>
  bool electron_abv_thresh(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
  {
    return (core::gOre::count_primaries(obj,  params.at(0), params.at(1), params.at(2), params.at(3)).at(pvars::kElectron) > 0);
  }
  REGISTER_CUT_SCOPE(RegistrationScope::True, electron_abv_thresh, electron_abv_thresh);

  /** @brief more than one shower above threshold**/
  template<class T>
  bool more_than_one_gOre(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
  {
    core::gOre::True_Interaction interaction(obj, params);
    return (interaction.ngOres() > 1);
  }
  REGISTER_CUT_SCOPE(RegistrationScope::True, more_than_one_gOre, more_than_one_gOre);

  /** @brief true fiducial contained single shower event **/
  template<class T>
  bool is_fid_con_gOre(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
  {
    core::gOre::True_Interaction interaction(obj, params);
    return interaction.is_valid && is_fid_con_nu(obj);
  }
  REGISTER_CUT_SCOPE(RegistrationScope::True, is_fid_con_gOre, is_fid_con_gOre);

  /**
   * @brief Check if MC interaction is fiducial contained NC Res event
   * @details This is a logical and of `is_fid_con_nc_nu` and `obj.interaction_mode == caf::kRes`
   * @param obj the interaction (MC only)
   * @return true if obj is a true fiducial contained NC neutrino event, false otherwise
   **/
  template<class T>
  bool is_fid_con_nc_nu_res(const T& obj)
  {
    return is_fid_con_nu(obj) && isnc(obj) && (obj.interaction_mode == caf::kRes);
  } 
  REGISTER_CUT_SCOPE(RegistrationScope::True, is_fid_con_nc_nu_res, is_fid_con_nc_nu_res);

  /**
   * @brief Check if the gOre is a photon produced from a decay
   * @details In the case the particle is primary the ancestor pdg code is the same as that of the particle in question
   * @param obj the interaction (MC only)
   * @return true if the object has a single photon and the photon was a decay product
   **/
  template<class T>
  bool decay_gOre_photon(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
  {
    core::gOre::True_Interaction interaction(obj, params);
    return interaction.is_valid &&
           interaction.primary_gOre()->pdg_code == 22 &&
           interaction.primary_gOre()->ancestor_pdg_code != 22;
  }
  REGISTER_CUT_SCOPE(RegistrationScope::True, decay_gOre_photon, decay_gOre_photon);

  /**
   * @brief Check if the gOre is a primary photon (eg not from a π0 decay)
   * @details In the case the particle is primary the ancestor pdg code is the same as that of the particle in question
   * @param obj the interaction (MC only)
   * @return true if the object has a single photon and the photon was a primary
   **/
  template<class T>
  bool primary_gOre_photon(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
  {
    core::gOre::True_Interaction interaction(obj, params);
    return interaction.is_valid &&
           interaction.primary_gOre()->pdg_code == 22 &&
           interaction.primary_gOre()->ancestor_pdg_code == 22;
  }
  REGISTER_CUT_SCOPE(RegistrationScope::True, primary_gOre_photon, primary_gOre_photon);

  template<class T>
  bool debug(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
  {
    core::gOre::True_Interaction interaction(obj, params);
    if (not interaction.is_valid)
      return false;
    int64_t          pdg_code = interaction.primary_gOre()->pdg_code; 
    int64_t ancestor_pdg_code = interaction.primary_gOre()->ancestor_pdg_code;
    int64_t interaction_mode = obj.interaction_mode;
    if (ancestor_pdg_code == 22 && pdg_code == 22 && interaction_mode == caf::kRes)
    {
      std::cout << "Primary photon! SRInteractionTruthDLP Info:"                << '\n'
                << "  bjorken_x: "               << obj.bjorken_x               << '\n'
                << "  cathode_offset: "          << obj.cathode_offset          << '\n'
                << "  current_type: "            << obj.current_type            << '\n'
                << "  depositions_adapt_q_sum: " << obj.depositions_adapt_q_sum << '\n'
                << "  depositions_adapt_sum: "   << obj.depositions_adapt_sum   << '\n'
                << "  depositions_g4_sum: "      << obj.depositions_g4_sum      << '\n'
                << "  depositions_q_sum: "       << obj.depositions_q_sum       << '\n'
                << "  depositions_sum: "         << obj.depositions_sum         << '\n'
                << "  distance_travel: "         << obj.distance_travel         << '\n'
                << "  energy_init: "             << obj.energy_init             << '\n'
                << "  energy_transfer: "         << obj.energy_transfer         << '\n'
                << "  flash_hypo_pe: "           << obj.flash_hypo_pe           << '\n'
                << "  flash_total_pe: "          << obj.flash_total_pe          << '\n'
                << "  hadronic_invariant_mass: " << obj.hadronic_invariant_mass << '\n'
                << "  id: "                      << obj.id                      << '\n'
                << "  inelasticity: "            << obj.inelasticity            << '\n'
                << "  interaction_mode: "        << obj.interaction_mode        << '\n'
                << "  interaction_type: "        << obj.interaction_type        << '\n'
                << "  is_cathode_crosser: "      << obj.is_cathode_crosser      << '\n'
                << "  is_contained: "            << obj.is_contained            << '\n'
                << "  is_fiducial: "             << obj.is_fiducial             << '\n'
                << "  is_flash_matched: "        << obj.is_flash_matched        << '\n'
                << "  is_matched: "              << obj.is_matched              << '\n'
                << "  is_truth: "                << obj.is_truth                << '\n'
                << "  lepton_p: "                << obj.lepton_p                << '\n'
                << "  lepton_pdg_code: "         << obj.lepton_pdg_code         << '\n'
                << "  lepton_track_id: "         << obj.lepton_track_id         << '\n'
                << "  mct_index: "               << obj.mct_index               << '\n'
                << "  momentum_transfer: "       << obj.momentum_transfer       << '\n'
                << "  momentum_transfer_mag: "   << obj.momentum_transfer_mag   << '\n'
                << "  nu_id: "                   << obj.nu_id                   << '\n'
                << "  nucleon: "                 << obj.nucleon                 << '\n'
                << "  num_particles: "           << obj.num_particles           << '\n'
                << "  num_primary_particles: "   << obj.num_primary_particles   << '\n'
                << "  orig_id: "                 << obj.orig_id                 << '\n'
                << "  pdg_code: "                << obj.pdg_code                << '\n'
                << "  quark: "                   << obj.quark                   << '\n'
                << "  size: "                    << obj.size                    << '\n'
                << "  size_adapt: "              << obj.size_adapt              << '\n'
                << "  size_g4: "                 << obj.size_g4                 << '\n'
                << "  target: "                  << obj.target                  << '\n'
                << "  theta: "                   << obj.theta                   << '\n'
                << "  track_id: "                << obj.track_id                << std::endl;
    }
    return true;
  }
  REGISTER_CUT_SCOPE(RegistrationScope::True, debug, debug);

  template<class T>
  bool mc_debug(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY,
                                                            GORE_FID_THRESH_X_POS, GORE_FID_THRESH_X_NEG, GORE_FID_THRESH_Y_POS, GORE_FID_THRESH_Y_NEG, GORE_FID_THRESH_Z_POS, GORE_FID_THRESH_Z_NEG})
  {
    // check only NC ∆ Res events
    bool isnc = obj.isnc;
    int resnum = obj.resnum;
    caf::genie_interaction_type_ genie_inttype = obj.genie_inttype;
    bool non_pion = (genie_inttype != caf::kResNCNuProtonPi0)       &&
                    (genie_inttype != caf::kResNCNuProtonPiPlus)    &&
                    (genie_inttype != caf::kResNCNuNeutronPi0)      &&
                    (genie_inttype != caf::kResNCNuNeutronPiMinus)  &&
                    (genie_inttype != caf::kResNCNuBarProtonPi0)    &&
                    (genie_inttype != caf::kResNCNuBarProtonPiPlus) &&
                    (genie_inttype != caf::kResNCNuBarNeutronPi0)   &&
                    (genie_inttype != caf::kResNCNuBarNeutronPiMinus);
    bool nc_delta_res = isnc && (resnum == 0);
    if (not nc_delta_res || not non_pion)
      return false;
    std::cout << "GTruth INFO: \n"
              << "  ischarm: "    << obj.ischarm
              << "  isseaquark: " << obj.isseaquark
              << "  xsec: "       << obj.xsec << std::endl;
    std::cout << "PRE-FSI COUNTS: \n"
              << "  π+: " << obj.npiplus
              << "  π-: " << obj.npiminus
              << "  π0: " << obj.npizero
              << "   p: " << obj.nproton
              << "   n: " << obj.nneutron << std::endl;
    core::gOre::mc_topology topology(obj.prim, params);
    topology.report();
    return true;
  }
  REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, mc_debug, mc_debug);

  //*** MC CUTS***//
  /**
   * @brief is the interaction a true NC Δ res event?
   * @param obj The MC Truth of the interaction
   **/
  template<class T> 
  bool nc_delta_res(const T& obj)
  {
    return (obj.isnc && (obj.resnum == 0));
  }
  REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, nc_delta_res, nc_delta_res);

  /**
   * @brief is the interaction a true NC Δ res event which isn't catagorized as producing pions?
   * @details check the GENIE interaction type against the different NC Res modes enumerated
   * Would prefer to have the photon producing resonancy to check instead, but w/e
   * @param The MC Truth of the interaction
   **/
  template<class T>
  bool nc_delta_res_no_pion(const T& obj)
  {
    return nc_delta_res(obj) &&
           (obj.genie_inttype != caf::kResNCNuProtonPi0)       &&
           (obj.genie_inttype != caf::kResNCNuProtonPiPlus)    &&
           (obj.genie_inttype != caf::kResNCNuNeutronPi0)      &&
           (obj.genie_inttype != caf::kResNCNuNeutronPiMinus)  &&
           (obj.genie_inttype != caf::kResNCNuBarProtonPi0)    &&
           (obj.genie_inttype != caf::kResNCNuBarProtonPiPlus) &&
           (obj.genie_inttype != caf::kResNCNuBarNeutronPi0)   &&
           (obj.genie_inttype != caf::kResNCNuBarNeutronPiMinus);
  }
  REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, nc_delta_res_no_pion, nc_delta_res_no_pion);

  /**
   * @brief is the interaction a true NC Δ res event which produces a pi0?
   * @details check the GENIE interaction type against the different NC Res modes enumerated
   * Would prefer to have the photon producing resonancy to check instead, but w/e
   * @param The MC Truth of the interaction
   **/
  template <class T>
  bool nc_delta_res_pi0(const T& obj)
  {
    return nc_delta_res(obj) &&
           ((obj.genie_inttype == caf::kResNCNuProtonPi0)       ||
            (obj.genie_inttype == caf::kResNCNuNeutronPi0)      ||
            (obj.genie_inttype == caf::kResNCNuBarProtonPi0)    ||
            (obj.genie_inttype == caf::kResNCNuBarNeutronPi0)    );
  }
  REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, nc_delta_res_pi0, nc_delta_res_pi0);
} // end cuts::gOre namespace

 #endif

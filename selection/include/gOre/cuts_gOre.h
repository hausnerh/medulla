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
 * @namespace cuts::nc::gOre
 * @brief Cut specific to NC single photon analyses
 **/
namespace cuts::nc::gOre
{
  /**
   * @brief Does the interaction have a single photon or single electron topology?
   * @details This cut makes use of the thresholds defined in core::nc::gOre::Interaction.
   * Require there are no muons or pions above threshold, either a single photon or a single electron,
   * but no requirement on the protons. The photon and electrons are to be distinguished later with an
   * optimized PID cut.
   * @tparam T the type of interaction (true or reco).
   * @param obj the interaction in question
   * @return true if there is a single shower-like particle and no muons or pions, false otherwise.
   **/
  template<class T>
    bool gOre_topology(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      return interaction.is_valid;
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
    bool gOre_target_KE(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
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
    bool gOre_0p(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      return (interaction.is_valid) && (interaction.nProtons() == 0);
    }
  REGISTER_CUT_SCOPE(RegistrationScope::Both, gOre_0p, gOre_0p); 

  /** @brief gOre_topology & 1 proton **/
  template<class T>
    bool gOre_1p(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      return (interaction.is_valid) && (interaction.nProtons() == 1);
    }
  REGISTER_CUT_SCOPE(RegistrationScope::Both, gOre_1p, gOre_1p);

  /** @brief gOre_topology & at least 1 proton **/
  template<class T>
    bool gOre_Np(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
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
    bool has_subthreshold_gOre(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
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
  bool gOre_is_photon(const T& obj)
  {
    core::nc::gOre::True_Interaction interaction(obj);
    return interaction.is_valid && (pvars::pid(*interaction.primary_gOre()) == pvars::kPhoton);
  }
  REGISTER_CUT_SCOPE(RegistrationScope::True, gOre_is_photon, gOre_is_photon);

  /** @brief true single electron event **/
  template<class T>
  bool gOre_is_electron(const T& obj)
  {
    core::nc::gOre::True_Interaction interaction(obj);
    return interaction.is_valid && (pvars::pid(*interaction.primary_gOre()) == pvars::kElectron);
  }
  REGISTER_CUT_SCOPE(RegistrationScope::True, gOre_is_electron, gOre_is_electron);

  /** @brief has muon above threshold **/
  template<class T>
  bool muon_abv_thresh(const T& obj)
  {
    return (core::nc::gOre::count_primaries(obj).at(pvars::kMuon) > 0);
  }
  REGISTER_CUT_SCOPE(RegistrationScope::True, muon_abv_thresh, muon_abv_thresh);

  /** @brief has pion above threshold **/
  template<class T>
  bool pion_abv_thresh(const T& obj)
  {
    return (core::nc::gOre::count_primaries(obj).at(pvars::kMuon) > 0);
  }
  REGISTER_CUT_SCOPE(RegistrationScope::True, pion_abv_thresh, pion_abv_thresh);

  /** @brief has electron above threshold **/
  template<class T>
  bool electron_abv_thresh(const T& obj)
  {
    return (core::nc::gOre::count_primaries(obj).at(pvars::kElectron) > 0);
  }
  REGISTER_CUT_SCOPE(RegistrationScope::True, electron_abv_thresh, electron_abv_thresh);

  /** @brief more than one shower above threshold**/
  template<class T>
  bool more_than_one_gOre(const T& obj)
  {
    core::nc::gOre::True_Interaction interaction(obj);
    return (interaction.ngOres() > 1);
  }
  REGISTER_CUT_SCOPE(RegistrationScope::True, more_than_one_gOre, more_than_one_gOre);

  /** @brief true fiducial contained single shower event **/
  template<class T>
  bool is_fid_con_gOre(const T& obj)
  {
    core::nc::gOre::True_Interaction interaction(obj);
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
  bool decay_gOre_photon(const T& obj)
  {
    core::nc::gOre::True_Interaction interaction(obj);
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
  bool primary_gOre_photon(const T& obj)
  {
    core::nc::gOre::True_Interaction interaction(obj);
    return interaction.is_valid &&
           interaction.primary_gOre()->pdg_code == 22 &&
           interaction.primary_gOre()->ancestor_pdg_code == 22;
  }
  REGISTER_CUT_SCOPE(RegistrationScope::True, primary_gOre_photon, primary_gOre_photon);

  template<class T>
  bool debug(const T& obj)
  {
    core::nc::gOre::True_Interaction interaction(obj);
    if (not interaction.is_valid)
      return false;
    int64_t ancestor_pdg_code = interaction.primary_gOre()->ancestor_pdg_code;
    int64_t parent_pdg_code = interaction.primary_gOre()->parent_pdg_code;
    std::cout << "***\n"
              << "Ancestor " << ancestor_pdg_code << '\n'
              << "Parent " << ancestor_pdg_code << std::endl;
    return true;
  }
  REGISTER_CUT_SCOPE(RegistrationScope::True, debug, debug);

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
} // end cuts::nc::gOre namespace

 #endif

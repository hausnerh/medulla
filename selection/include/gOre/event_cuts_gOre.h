/**
 * @file event_cuts_gOre.h
 * @brief Header file for defining event cuts useful for NC analyses
 * specficially for saving photon/electron discrimination for the end
 * @author hhausner@fnal.gov
 **/

#ifndef EVENT_CUTS_GORE_H
#define EVENT_CUTS_GORE_H

#include "include/event_variables.h"

#include "include/gOre/cuts_gOre.h"

/**
 * @namespace ecuts::gOre
 * @brief Event cut specific to NC single photon analyses
 **/
namespace ecuts::gOre
{
  /**
   * @brief Does the event contain a true NC ∆->Nγ?
   * @details This is an event cut so it can be supplied for our signal trees
   * @tparam T the top-level record.
   * @param sr the StandardRecord to apply the cut on.
   * @return true if the event has a true NC ∆->Nγ 
   **/
  template <typename T>
    bool has_nc_delta_ng_interaction(const T& sr, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY,
                                                                                GORE_FID_THRESH_X_POS, GORE_FID_THRESH_X_NEG, GORE_FID_THRESH_Y_POS, GORE_FID_THRESH_Y_NEG, GORE_FID_THRESH_Z_POS, GORE_FID_THRESH_Z_NEG})
    {
      // loop over the sr truth info and check if there's a hit
      for (auto const& interaction : sr.mc.nu)
      {
        bool isnc = interaction.isnc;
        bool delta_res = (interaction.resnum == 0);
        core::gOre::mc_topology topology(interaction.prim, params);
        bool single_photon = topology.single_photon() && topology.only_photons_and_nucleons();
        if (isnc && delta_res && single_photon)
          return true;
      }
      return false;
    }
  REGISTER_CUT_SCOPE(RegistrationScope::Event, has_nc_delta_ng_interaction, has_nc_delta_ng_interaction);

  /**
   * @brief Does the event contain only a true NC ∆->Nγ?
   * @details This is an event cut so it can be supplied for our signal trees
   * @tparam T the top-level record.
   * @param sr the StandardRecord to apply the cut on.
   * @return true if the event has only a true NC ∆->Nγ 
   **/
  template <typename T>
    bool only_nc_delta_ng_interaction(const T& sr)
    {
      return (sr.mc.nnu == 1) && has_nc_delta_ng_interaction(sr);
    }
  REGISTER_CUT_SCOPE(RegistrationScope::Event, only_nc_delta_ng_interaction, only_nc_delta_ng_interaction);

  /**
   * @brief veto gOre interactions where there is another interaction in the event with a shower over threshold
   * @details look at the recos. if there is no gOre, reject. if there is a gOre but no other photon producing interaction, accept.
   * if there is a gOre and another photon producing interaction, reject.
   * @tparam T the top-level record.
   * @param sr the StandardRecord to apply the cut on.
   * @return true only if there is a gOre interaction and no other shower producing interaction
   **/
  template <typename T>
    bool gOre_veto(const T& sr, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      size_t nreco = sr.ndlp;
      if (nreco < 2)
        return false;
      bool gOre(false);
      bool otherPhoton(false);
      for (auto const& reco : sr.dlp)
      {
        core::gOre::Reco_Interaction interaction(reco, params);
        bool this_gOre = interaction.is_valid;
        bool this_photon(false);
        if (not this_gOre)
          for (auto const& gOre : interaction.gOres)
            this_photon = (this_photon || (gOre->pdg_code == 22));
        gOre = (gOre || this_gOre);
        otherPhoton = (otherPhoton || this_photon);
      }
      return gOre && not otherPhoton;
    }
  REGISTER_CUT_SCOPE(RegistrationScope::Event, gOre_veto, gOre_veto);

  /**
   * @brief veto events which have two photons in different reco interactions which get close to the pion mass peak and a
   * "vertex" for the decay that can be located to within 5 cm
   * @details look over the reco events for what looks like a split π0. check if the primary photons can reconstruct the
   * π0 mass to within 10% (configurable with the params). aslo check that projecting the photons backwards we find a separation of
   * less than 5 cm. veto those events.
   * @tparam T the top-level record.
   * @param sr the StandardRecord to apply the cut on.
   * @param params the fraction of how close to the π0 mass needed to be vetoed (default 10%)
   * @return true if there is not what looks like a split π0
   **/
  template <typename T>
    bool pion_mass_veto(const T& sr, std::vector<double> params={0.1, 5.0, GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      double frac_diff_limit = (params.size() > 0) ? params.at(0) : 0.1;
      double vtx_resolution  = (params.size() > 1) ? params.at(1) : 5.0;
      double gOre_theshold   = (params.size() > 2) ? params.at(2) : GORE_MIN_GORE_ENERGY;
      double muon_theshold   = (params.size() > 3) ? params.at(3) : GORE_MIN_MUON_ENERGY;
      double pion_theshold   = (params.size() > 4) ? params.at(4) : GORE_MIN_PION_ENERGY;
      if (not (frac_diff_limit > 0))
        std::__throw_domain_error(("π0 mass veto acceptance set to "+std::to_string(frac_diff_limit)+", but must be greater than 0").c_str());

      // ignore if there's not at least 2 recos
      size_t nreco = sr.ndlp;
      if (nreco < 2)
        return true;
     
      // make a list of the primary photons, their directions, and vertices
      std::vector<std::tuple<double, utilities::three_vector, utilities::three_vector>> gOre_ke_dir_vtx;
      std::vector<std::tuple<double, utilities::three_vector, utilities::three_vector>> photon_ke_dir_vtx;
      for (auto const& reco : sr.dlp)
      {
        bool is_sel_gOre = (cuts::fiducial_cut     (reco)                  &&
                            cuts::containment_cut  (reco)                  &&
                            cuts::valid_flashmatch (reco)                  &&
                            cuts::gOre::single_gOre(reco, {gOre_theshold}) &&
                            cuts::no_muons         (reco, {muon_theshold}) &&
                            cuts::no_charged_pions (reco, {pion_theshold}) );
        if (is_sel_gOre)
        {
          size_t photon_idx   = selectors::leading_photon(reco);
          size_t electron_idx = selectors::leading_electron(reco);
          if (photon_idx == kNoMatch && electron_idx == kNoMatch)
          {
            std::cout << "ecuts::gOre::pion_mass_veto - Could not find valid indices for gOre in selected gOre event." << std::endl;
            continue;
          }
          caf::SRParticleDLPProxy& gOre = (photon_idx != kNoMatch) ? reco.particles.at(photon_idx) : reco.particles.at(electron_idx);
          gOre_ke_dir_vtx.emplace_back(pvars::ke(gOre), utilities::to_three_vector(gOre.start_dir), utilities::to_three_vector(gOre.start_point));
        } else
        {
          size_t photon_idx = selectors::leading_photon(reco);
          if (photon_idx == kNoMatch)
            continue;
          caf::SRParticleDLPProxy& photon = reco.particles.at(photon_idx);
          photon_ke_dir_vtx.emplace_back(pvars::ke(photon), utilities::to_three_vector(photon.start_dir), utilities::to_three_vector(photon.start_point));
        }
      }

      // if there's not exactly one gOre return false
      if (gOre_ke_dir_vtx.size() != 1)
        return false;
      
      // if there's not at least one non-gOre photon event there's nothing to veto
      if (photon_ke_dir_vtx.size() == 0)
        return true;

      // loop over the photon pairs and compare the reco π0 mass to the true value
      for (size_t idx1 = 0; idx1 < gOre_ke_dir_vtx.size(); ++idx1)
      {
        auto const& [ke1, dir1, vtx1] = gOre_ke_dir_vtx.at(idx1);
        for (size_t idx2 = 0; idx2 < photon_ke_dir_vtx.size(); ++idx2)
        {
          auto const& [ke2, dir2, vtx2] = photon_ke_dir_vtx.at(idx2);
          double reco_pion_mass = std::sqrt(ke1*ke2*(1. - utilities::dot_product(dir1, dir2)));
          double frac_diff = std::abs(1. - (reco_pion_mass / 134.98));
          bool near_pion_mass = (frac_diff < frac_diff_limit);
          bool not_parallel = (utilities::dot_product(dir1, dir2) != 1);
          utilities::three_vector cross = utilities::cross_product(dir1, dir2);
          double dist = (not_parallel) ? std::abs(utilities::dot_product(cross, utilities::subtract(vtx1, vtx2))) / utilities::magnitude(cross)
                                       : utilities::magnitude(utilities::cross_product(utilities::subtract(vtx1, vtx2), dir1));
          if (near_pion_mass && (dist < vtx_resolution))
            return false;
        }
      }

      // if you make it here you're in the clear
      return true;
    }
  REGISTER_CUT_SCOPE(RegistrationScope::Event, pion_mass_veto, pion_mass_veto);

  /**
   * @brief veto any event with more than N reco interaction
   * @tparam T the top-level record.
   * @param sr the StandardRecord to apply the cut on.
   * @param params the number of reco interactions over which we veto the event (default to N = 1)
   * @return true if there N or fewer reco interactions
   **/
  template <typename T>
    bool n_reco_veto(const T& sr, std::vector<double> params={1})
    {
      return (sr.ndlp <= params.at(0));
    }
  REGISTER_CUT_SCOPE(RegistrationScope::Event, n_reco_veto, n_reco_veto);

  /**
   * @brief Does the event contain a selected reconstructed interaction?
   * @tparam T the top-level record
   * @param sr the StandardRecord to apply the cut on
   * @param params the threshold energies for the cuts
   * @return true if there is at least one reco interaction which passes the selection
   **/
  template <typename T>
    bool has_selected_gOre(const T& sr, std::vector<double> params={GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      for (auto const& reco : sr.dlp)
      {
        if (cuts::fiducial_cut     (reco)                 &&
            cuts::containment_cut  (reco)                 &&
            cuts::valid_flashmatch (reco)                 &&
            cuts::gOre::single_gOre(reco, {params.at(0)}) &&
            cuts::no_muons         (reco, {params.at(1)}) &&
            cuts::no_charged_pions (reco, {params.at(2)}))
          return true;
      }
      return false;
    }
  REGISTER_CUT_SCOPE(RegistrationScope::Event, has_selected_gOre, has_selected_gOre);
  
} // end ecuts::gOre

#endif
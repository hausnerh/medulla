/**
 * @file event_cuts_gOre.h
 * @brief Header file for defining event cuts useful for NC analyses
 * specficially for saving photon/electron discrimination for the end
 * @author hhausner@fnal.gov
 **/

#ifndef EVENT_CUTS_GORE_H
#define EVENT_CUTS_GORE_H

#include "include/event_variables.h"

#include "include/gOre/core_gOre.h"

/**
 * @namespace ecuts::nc::gOre
 * @brief Event cut specific to NC single photon analyses
 **/
namespace ecuts::nc::gOre
{
  /**
   * @brief Does the event contain a true NC ∆->Nγ?
   * @details This is an event cut so it can be supplied for our signal trees
   * @tparam T the top-level record.
   * @param sr the StandardRecord to apply the cut on.
   * @return true if the event has a true NC ∆->Nγ 
   **/
  template <typename T>
    bool has_nc_delta_ng_interaction(const T& sr)
    {
      // loop over the sr truth info and check if there's a hit
      for (auto const& interaction : sr.mc.nu)
      {
        bool isnc = interaction.isnc;
        bool delta_res = (interaction.resnum == 0);
        core::nc::gOre::mc_topology topology(interaction.prim);
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
    bool gOre_veto(const T& sr)
    {
      size_t nreco = sr.ndlp;
      if (nreco < 2)
        return false;
      bool gOre(false);
      bool otherPhoton(false);
      for (auto const& reco : sr.dlp)
      {
        core::nc::gOre::Reco_Interaction interaction(reco);
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
   * @brief veto events which have two photons in different reco interactions which get close to the pion mass peak
   * @details look over the reco events for what looks like a split π0. check if the primary photons can reconstruct the
   * π0 mass to within 10% (configurable with the params). veto those events.
   * @tparam T the top-level record.
   * @param sr the StandardRecord to apply the cut on.
   * @param params the fraction of how close to the π0 mass needed to be vetoed (default 10%)
   * @return true if there is not what looks like a split π0
   **/
  template <typename T>
    bool pion_mass_veto(const T& sr, std::vector<double> params={})
    {
      double frac_diff_limit = (params.size() == 1) ? params.at(0) : 0.1;
      if (frac_diff_limit == std::numeric_limits<double>::min())
        std::cout << "frac_diff_limit set to  std::numeric_limits<double>::min()" << std::endl;
      if (not (frac_diff_limit > 0))
        std::__throw_domain_error(("π0 mass veto acceptance set to "+std::to_string(frac_diff_limit)+", but must be greater than 0").c_str());

      // ignore if there's not at least 2 recos
      size_t nreco = sr.ndlp;
      if (nreco < 2)
        return true;
     
      // make a list of the primary photons and their directions
      std::vector<std::pair<double, utilities::three_vector>> photon_ke_dir;  
      for (auto const& reco : sr.dlp)
      {
        size_t photon_idx = selectors::leading_photon(reco);
        if (photon_idx == kNoMatch)
          continue;
        caf::SRParticleDLPProxy& photon = reco.particles.at(photon_idx);
        photon_ke_dir.emplace_back(pvars::ke(photon), utilities::to_three_vector(photon.start_dir));
      }

      // ignore if there aren't at least 2 photons
      if (photon_ke_dir.size() < 2)
        return true;
      
      // loop over the photon pairs and compare the reco π0 mass to the true value
      for (size_t idx1 = 1; idx1 < photon_ke_dir.size(); ++idx1)
      {
        auto const& [ke1, dir1] = photon_ke_dir.at(idx1);
        for (size_t idx2 = 0; idx2 < idx1; ++idx2)
        {
          auto const& [ke2, dir2] = photon_ke_dir.at(idx2);
          double reco_pion_mass = std::sqrt(ke1*ke2*(1. - utilities::dot_product(dir1, dir2)));
          double frac_diff = std::abs(1. - (reco_pion_mass / 134.98));
          if (frac_diff < frac_diff_limit)
            return false;
        }
      }

      // if you make it here you're in the clear
      return true;
    }
  REGISTER_CUT_SCOPE(RegistrationScope::Event, pion_mass_veto, pion_mass_veto);
  
} // end ecuts::nc::gOre

#endif
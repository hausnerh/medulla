
/**
 * @file event_vars_gOre.h
 * @brief Header file for defining event variables useful for NC analyses
 * specficially for saving photon/electron discrimination for the end
 * @author hhausner@fnal.gov
 **/

#ifndef EVENT_VARS_GORE_H
#define EVENT_VARS_GORE_H

#include "include/gOre/event_cuts_gOre.h"

/**
 * @namespace evars::gOre
 * @brief Event variables specific to NC single photon analyses
 **/
namespace evars::gOre
{
  /**
   * @brief How many interactions in the event have a photon above threshold
   * @tparam T the top-level record.
   * @param sr the StandardRecord to apply the cut on.
   * @return the number of interactions reconstructed with a photon above threshold
   **/
  template <typename T>
    double n_photon_dlp(const T& sr)
    {
      double counts(0);
      for (auto const& reco : sr.dlp)
      {
        size_t photon_idx = selectors::leading_photon(reco);
        if (photon_idx != kNoMatch)
          counts += 1;
      }
      return counts;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Event, n_photon_dlp, n_photon_dlp);

  /**
   * @brief If there are more than two photon producing events, how close does it look the the photons were from eachother
   * if we project their directions backwards
   * @tparam T the top-level record.
   * @param sr the StandardRecord to apply the cut on.
   * @return distance of closest approach of the primary photons in two reco interactions
   **/
  template <typename T>
    double photons_closest_dist(const T& sr)
    {
      // make a list of the primary photons, their directions, and vertices
      std::vector<std::tuple<utilities::three_vector, utilities::three_vector>> photon_dir_vtx;
      for (auto const& reco : sr.dlp)
      {
        size_t photon_idx = selectors::leading_photon(reco);
        if (photon_idx == kNoMatch)
          continue;
        caf::SRParticleDLPProxy& photon = reco.particles.at(photon_idx);
        photon_dir_vtx.emplace_back(utilities::to_three_vector(photon.start_dir), utilities::to_three_vector(photon.start_point));
      }

      // ignore if there aren't at least 2 photons
      if (photon_dir_vtx.size() < 2)
        return std::numeric_limits<double>::quiet_NaN();;
     
      // loop over the photon pairs and compare the reco π0 mass to the true value
      double min_dist(std::numeric_limits<double>::max());
      for (size_t idx1 = 1; idx1 < photon_dir_vtx.size(); ++idx1)
      {
        auto const& [dir1, vtx1] = photon_dir_vtx.at(idx1);
        for (size_t idx2 = 0; idx2 < idx1; ++idx2)
        {
          auto const& [dir2, vtx2] = photon_dir_vtx.at(idx2);
          bool not_parallel = (utilities::dot_product(dir1, dir2) != 1);
          utilities::three_vector cross = utilities::cross_product(dir1, dir2);
          double dist = (not_parallel) ? std::abs(utilities::dot_product(cross, utilities::subtract(vtx1, vtx2))) / utilities::magnitude(cross)
                                       : utilities::magnitude(utilities::cross_product(utilities::subtract(vtx1, vtx2), dir1));
          if (dist < min_dist)
            min_dist = dist;
        }
      }
      return min_dist;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Event, photons_closest_dist, photons_closest_dist);

  /**
   * @brief If there are more than two photon producing events, how close does the event get to the
   * π0 mass peak assuming the photons were split?
   * @tparam T the top-level record.
   * @param sr the StandardRecord to apply the cut on.
   * @return 1 - (reco m_π0 / true m_π0), closest to zero if there's more than 2 photonic events
   **/
  template <typename T>
    double split_pion_mass(const T& sr)
    {
      // make a list of the primary photons, their directions, and vertices
      std::vector<std::tuple<double, utilities::three_vector>> photon_ke_dir;
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
        return std::numeric_limits<double>::quiet_NaN();;
      
      // loop over the photon pairs and compare the reco π0 mass to the true value
      double min_frac_diff(std::numeric_limits<double>::max());
      for (size_t idx1 = 1; idx1 < photon_ke_dir.size(); ++idx1)
      {
        auto const& [ke1, dir1] = photon_ke_dir.at(idx1);
        for (size_t idx2 = 0; idx2 < idx1; ++idx2)
        {
          auto const& [ke2, dir2] = photon_ke_dir.at(idx2);
          double reco_pion_mass = std::sqrt(ke1*ke2*(1. - utilities::dot_product(dir1, dir2)));
          double frac_diff = std::abs(1. - (reco_pion_mass / 134.98));
          if (frac_diff < min_frac_diff)
            min_frac_diff = frac_diff;
        }
      }
      return min_frac_diff;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Event, split_pion_mass, split_pion_mass);

  /**
   * @brief Count the number of MCTruths in the event of vars::gOre::mc_category cat
   * @tparam T the top-level record.
   * @param sr the StandardRecord to apply the cut on.
   * @param cat the category (see vars::gOre::mc_category)
   * @return double how many MCTruths are of the specified category?
   **/
  template <typename T>
    double n_mc_category(const T& sr, double cat)
    {
      double counts(0);
      for (auto const& interaction : sr.mc.nu)
        if (vars::gOre::mc_category(interaction) == cat)
          counts += 1;
      return counts;
    }

  /**
   * @brief Count the number of MCTruths in the event of vars::gOre::mc_category 0
   * @tparam T the top-level record.
   * @param sr the StandardRecord to apply the cut on.
   * @return double how many MCTruths are of category 0?
   **/
  template <typename T>
    double n_mc_category_0(const T& sr)
    {
      return n_mc_category(sr, 0);
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Event, n_mc_category_0, n_mc_category_0);

  /**
   * @brief Count the number of MCTruths in the event of vars::gOre::mc_category 1
   * @tparam T the top-level record.
   * @param sr the StandardRecord to apply the cut on.
   * @return double how many MCTruths are of category 1?
   **/
  template <typename T>
    double n_mc_category_1(const T& sr)
    {
      return n_mc_category(sr, 1);
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Event, n_mc_category_1, n_mc_category_1);

  /**
   * @brief Count the number of MCTruths in the event of vars::gOre::mc_category 2
   * @tparam T the top-level record.
   * @param sr the StandardRecord to apply the cut on.
   * @return double how many MCTruths are of category 2?
   **/
  template <typename T>
    double n_mc_category_2(const T& sr)
    {
      return n_mc_category(sr, 2);
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Event, n_mc_category_2, n_mc_category_2);

  /**
   * @brief Count the number of MCTruths in the event of vars::gOre::mc_category 3
   * @tparam T the top-level record.
   * @param sr the StandardRecord to apply the cut on.
   * @return double how many MCTruths are of category 3?
   **/
  template <typename T>
    double n_mc_category_3(const T& sr)
    {
      return n_mc_category(sr, 3);
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Event, n_mc_category_3, n_mc_category_3);

  /**
   * @brief Count the number of MCTruths in the event of vars::gOre::mc_category 4
   * @tparam T the top-level record.
   * @param sr the StandardRecord to apply the cut on.
   * @return double how many MCTruths are of category 4?
   **/
  template <typename T>
    double n_mc_category_4(const T& sr)
    {
      return n_mc_category(sr, 4);
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Event, n_mc_category_4, n_mc_category_4);

  /**
   * @brief Count the number of MCTruths in the event of vars::gOre::mc_category 5
   * @tparam T the top-level record.
   * @param sr the StandardRecord to apply the cut on.
   * @return double how many MCTruths are of category 5?
   **/
  template <typename T>
    double n_mc_category_5(const T& sr)
    {
      return n_mc_category(sr, 5);
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Event, n_mc_category_5, n_mc_category_5);

  /**
   * @brief Count the number of MCTruths in the event of vars::gOre::mc_category 6
   * @tparam T the top-level record.
   * @param sr the StandardRecord to apply the cut on.
   * @return double how many MCTruths are of category 6?
   **/
  template <typename T>
    double n_mc_category_6(const T& sr)
    {
      return n_mc_category(sr, 6);
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Event, n_mc_category_6, n_mc_category_6);
  
} // end namespace evars::gOre

#endif
/**
 * @file vars_gOre.h
 * @brief Header file for defining variables useful for NC analyses
 * specficially for saving photon/electron discrimination for the end
 * @author hhausner@fnal.gov
 **/

#ifndef VARS_GORE_H
#define VARS_GORE_H

#include "include/gOre/cuts_gOre.h"


/**
 * @namespace vars::gOre
 * @brief Variables specific to NC single photon analyses
 **/
namespace vars::gOre
{
  /**
   * @brief return if it passes the gOre topology cut
   * @param obj the interaction of interest (reco or MC)
   * @return 1 for passing, 0 for failing
   **/
  template <class T>
    double is_gOre(const T& obj, std::vector<double> params={GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      return cuts::gOre::gOre_topology(obj, params);
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, is_gOre, is_gOre);

  /**
   * @brief categorize the gOre into primary photon, secondary photon, or electron
   * @details check the ancestor pdg to determine if it is a primary (identical pdg to the gOre),
   * from some decay (the ancestor pdg will be the origin), or is an electron
   * @param obj the interaction of interest (MC only)
   * @return 0 for primary photon, 1 for secondary photon, 2 for electron, -1 for invalid
   **/
  template <class T>
    double gOre_category(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      if (not interaction.is_valid)
        return -1.;
      if (interaction.primary_gOre()->pdg_code != 22)
        return 2;
      if (interaction.primary_gOre()->ancestor_pdg_code != 22)
        return 1;
      return 0;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::True, gOre_category, gOre_category);

  /**
   * @brief the ancestor pdg code for the gOre
   * @details what particle from the original FS produced this shower? If it is the same as the particle in question it was a primary
   * @param obj the interaction of interest (MC only)
   * @return the pdg code of the ancestor, (NaN if invalid)
   **/
  template <class T>
    double gOre_ancestor(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return interaction.primary_gOre()->ancestor_pdg_code;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::True, gOre_ancestor, gOre_ancestor);

  /**
   * @brief which volumes have flash matches?
   * @return 0 for East only, 1 for West only, 2 for both, -1 for neither
   **/
  template <class T>
    double flash_volumes(const T& obj)
    {
      bool west(false);
      bool east(false);
      for (auto const& flsh : obj.flash_volume_ids)
      {
        if (flsh == 0)
          east = true;
        if (flsh == 1)
          west = true;
      }
      // interpreting west & east as bits, construct out return value
      int vol = ((west<<1) | east) - 1;
      return vol;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, flash_volumes, flash_volumes);

  /**
   * @brief how many protons are in the interaction?
   * @param obj the interaction of interest (data or MC)
   * @return the number of protons as a double
  **/
  template <class T>
    double n_protons(const T& obj, std::vector<double> params = {GORE_MIN_PROTON_ENERGY})
    {
      //core::gOre::Interaction<T> interaction(obj);
      //if (not interaction.is_valid)
      //  return std::numeric_limits<double>::quiet_NaN();
      //return static_cast<double>(interaction.nProtons());
      return cuts::particle_multiplicity(obj, std::numeric_limits<size_t>::max(), pvars::kProton, params);
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, n_protons, n_protons);

  /**
   * @brief the total KE of protons above threshold
   * @param obj the interaction of interest (data or MC)
   **/
  template <class T>
    double total_proton_ke(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      double ke = 0;
      for (auto const& proton : interaction.protons)
        ke += proton->ke;
      return ke;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, total_proton_ke, total_proton_ke);

  /**
   * @brief the minimum KE of subthreshold muons
   * @details if there are no subthreshold muons the energy is 0;
   **/
  template <class T>
    double min_muon_ke(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return interaction.min_muon_ke;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, min_muon_ke, min_muon_ke);

  /**
   * @brief the minimum KE of subthreshold pions
   * @details if there are no subthreshold pions the energy is 0;
   **/
  template <class T>
    double min_pion_ke(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return interaction.min_pion_ke;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, min_pion_ke, min_pion_ke);

  /**
   * @brief what is the kinetic energy of the candidate photon?
   * @param the interaction of interest (data or MC)
   * @return the kinetic energy of the photon/electron
  **/
  template <class T>
    double gOre_ke(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return interaction.primary_gOre()->ke;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, gOre_ke, gOre_ke);

  /**
   * @brief what is the kinetic energy of the subleading candidate photon?
   * @details if there is a candidate with a lower KE than threshold, what is it?
   * If there is no such candidate the return is std::numeric_limits<double>::lowest()
   * @param the interaction of interest (data or MC)
   * @return the kinetic energy of the photon/electron
  **/
  template <class T>
    double subleading_gOre_ke(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return interaction.subleading_gore_ke;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, subleading_gOre_ke, subleading_gOre_ke);

  /**
   * @brief axial spread of the candidate photon
   * @param the interaction of interest (data or MC)
   **/
  template <class T>
    double gOre_axial_spread(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return static_cast<double>(interaction.primary_gOre()->axial_spread);
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Reco, gOre_axial_spread, gOre_axial_spread);

  /**
   * @brief directional spread of the candidate photon
   * @param the interaction of interest (data or MC)
   **/
  template <class T>
    double gOre_directional_spread(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return static_cast<double>(interaction.primary_gOre()->directional_spread);
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Reco, gOre_directional_spread, gOre_directional_spread);

  /**
   * @brief length of the candidate photon
   * @param the interaction of interest (data or MC)
   **/
  template <class T>
    double gOre_length(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return static_cast<double>(interaction.primary_gOre()->length);
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, gOre_length, gOre_length);

  /**
   * @brief the momentum magnitude of the candicate photon 
   **/
  template <class T>
    double gOre_momentum(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return static_cast<double>(interaction.primary_gOre()->p);
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, gOre_momentum, gOre_momentum);

  /**
   * @brief the transerverse momentum of the candidate photon
   **/
  template <class T>
    double gOre_dpT(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return pvars::dpT(*interaction.primary_gOre());
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, gOre_dpT, gOre_dpT);

  /**
   * @brief the start dedx of the candidate photon
   **/
  template <class T>
    double gOre_start_dedx(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return interaction.primary_gOre()->start_dedx;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Reco, gOre_start_dedx, gOre_start_dedx);

  /**
   * @brief the straightness of the candidate photon
   **/
  template <class T>
    double gOre_straightness(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return interaction.primary_gOre()->start_straightness;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Reco, gOre_straightness, gOre_straightness);

  /**
   * @brief the vertex gap of the candidate photon
   **/
  template <class T>
    double gOre_gap(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      utilities::three_vector interaction_vtx = utilities::to_three_vector(obj.vertex);
      utilities::three_vector gOre_vtx = utilities::to_three_vector(interaction.primary_gOre()->start_point);
      return utilities::magnitude(utilities::subtract(gOre_vtx, interaction_vtx));
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, gOre_gap, gOre_gap);

  /**
   * @brief the photon/electron pid of the candidate photon
   * @details 1 is photon-like, -1 is electron-like
   **/
  template <class T>
    double gOre_score(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      std::vector<double> softmax_vec = pvars::particle_softmax_vec(*interaction.primary_gOre());
      return softmax_vec[pvars::kPhoton] - softmax_vec[pvars::kElectron];
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Reco, gOre_score, gOre_score);

  /**
   * @brief the photon/electron polar angle w.r.t the beam
   **/
  template <class T>
    double gOre_polar_angle(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return pvars::polar_angle(*interaction.primary_gOre());
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, gOre_polar_angle, gOre_polar_angle);

  /**
   * @brief the photon/electron azimuthal angle w.r.t the beam
   **/
  template <class T>
    double gOre_azimuthal_angle(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return pvars::azimuthal_angle(*interaction.primary_gOre());
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, gOre_azimuthal_angle, gOre_azimuthal_angle);

  /**
   * @brief distance between the vertex and the x-/y-side active volume boundaries
   * @details getting geom info from https://sbn-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=21693&filename=ICARUS_geometry_update_26Apr21_v3.pdf&version=3
   **/
  template <class T>
    double xy_wall_dist(const T& obj)
    {
      double min_dist = std::numeric_limits<double>::quiet_NaN();
      double vtx_x = vars::vertex_x(obj);
      double vtx_y = vars::vertex_y(obj);
      // taxi cab metric, so just min dist in each coordinate
      // vtx should be between, but if any are negative that's worth noting
      double pos_x_dist = (GORE_WALL_X_POS - vtx_x);
      double neg_x_dist = (vtx_x - GORE_WALL_X_NEG);
      double pos_y_dist = (GORE_WALL_Y_POS - vtx_y);
      double neg_y_dist = (vtx_y - GORE_WALL_Y_NEG);

      min_dist = std::min({pos_x_dist, neg_x_dist, pos_y_dist, neg_y_dist});

      return min_dist;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, xy_wall_dist, xy_wall_dist);

  /**
   * @brief distance between the vertex and the z-side active volume boundaries
   * @details getting geom info from https://sbn-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=21693&filename=ICARUS_geometry_update_26Apr21_v3.pdf&version=3
   **/
  template <class T>
    double z_wall_dist(const T& obj)
    {
      double min_dist = std::numeric_limits<double>::quiet_NaN();
      double vtx_z = vars::vertex_z(obj);
      // taxi cab metric, so just min dist in each coordinate
      // vtx should be between, but if any are negative that's worth noting
      double pos_z_dist = (GORE_WALL_Z_POS - vtx_z);
      double neg_z_dist = (vtx_z - GORE_WALL_Z_NEG);

      min_dist = std::min({pos_z_dist, neg_z_dist});

      return min_dist;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, z_wall_dist, z_wall_dist);

  /**
   * @brief sum of gOre shower KE (including sub-thresholds) 
   **/
  template <class T>
    double total_gOre_KE(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return interaction.total_gore_ke;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, total_gOre_KE, total_gOre_KE);

  /**
   * @brief how many subthreshold gOre showers are there?
   **/
  template <class T>
    double n_gOre_showers(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY})
    {
      core::gOre::Interaction<T> interaction(obj, params);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return interaction.nShowers;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, n_gOre_showers, n_gOre_showers);

  /**
   * @brief how many primary showers (above threshold) are in the event?
   **/
  template <class T>
    double n_primary_showers(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY})
    {
      size_t count(0);
      for(const auto & p : obj.particles)
      {
        bool is_shower = (pvars::pid(p) == pvars::kPhoton) || (pvars::pid(p) == pvars::kElectron); 
        if(is_shower && pvars::primary_classification(p) && pvars::ke(p) >= params[0])
          ++count;
      }
      return count;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, n_primary_showers, n_primary_showers);

  /**
   * @brief how many secondary showers (above threshold) are in the event?
   **/
  template <class T>
    double n_secondary_showers(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY})
    {
      size_t count(0);
      for(const auto & p : obj.particles)
      {
        bool is_shower = (pvars::pid(p) == pvars::kPhoton) || (pvars::pid(p) == pvars::kElectron); 
        if(is_shower && not pvars::primary_classification(p) && pvars::ke(p) >= params[0])
          ++count;
      }
      return count;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, n_secondary_showers, n_secondary_showers);

  /**
   * @brief approximate the pion_mass peak
   * @details m = sqrt(2*E_1*E_2*(1 - CosTh))
   **/
  template <class T>
    double pion_mass(const T& obj)
    {
      double pmass = -999.9;
      size_t idx_1 = selectors::gOre::leading_primary_gOre(obj);
      if (idx_1 == kNoMatch)
        return pmass;
      size_t idx_2 = selectors::gOre::subleading_primary_gOre(obj);
      if (idx_2 == kNoMatch)
        return pmass;
      auto const& gOre_1 = obj.particles.at(idx_1);
      double ke_1 = pvars::gOre::shower_ke(gOre_1);
      utilities::three_vector dir_1 = utilities::to_three_vector(gOre_1.start_dir);
      auto const& gOre_2 = obj.particles.at(idx_2);
      double ke_2 = pvars::gOre::shower_ke(gOre_2);
      utilities::three_vector dir_2 = utilities::to_three_vector(gOre_2.start_dir);
      double cosTh = utilities::dot_product(dir_1, dir_2);
      pmass = std::sqrt(2.*ke_1*ke_2*(1.-cosTh));
      return pmass;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Reco, pion_mass, pion_mass);

  /**
   * @brief approximate the Delta baryon/nucleon mass square splitting
   * @details if there is a proton, this is possible
   **/
  template <class T>
    double baryon_mass_splitting(const T& obj)
    {
      double massSqSplit = std::numeric_limits<double>::quiet_NaN();
      size_t idx_gamma = selectors::gOre::leading_primary_gOre(obj);
      size_t idx_proton = selectors::gOre::leading_primary_proton(obj);
      if (idx_gamma == kNoMatch || idx_proton == kNoMatch)
        return massSqSplit;
      auto const& gamma  = obj.particles.at(idx_gamma);
      auto const& proton = obj.particles.at(idx_proton);
      double p_gamma  = pvars::gOre::shower_p(gamma);
      double p_proton = pvars::p(proton);
      utilities::three_vector dir_gamma  = utilities::to_three_vector(gamma.start_dir);
      utilities::three_vector dir_proton = utilities::to_three_vector(proton.start_dir);
      double cosTh = utilities::dot_product(dir_gamma, dir_proton);
      massSqSplit = 2.*p_gamma*(std::sqrt(p_proton*p_proton + PROTON_MASS*PROTON_MASS) - p_proton*cosTh);
      return massSqSplit;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, baryon_mass_splitting, baryon_mass_splitting);

  /**
   * @brief approximate the Delta baryon mass
   * @details if there is a proton, this is possible
   **/
  template <class T>
    double delta_mass(const T& obj)
    {
      double mass = -999.9;
      size_t idx_gamma = selectors::gOre::leading_primary_gOre(obj);
      size_t idx_proton = selectors::gOre::leading_primary_proton(obj);
      if (idx_gamma == kNoMatch || idx_proton == kNoMatch)
        return mass;
      auto const& gamma  = obj.particles.at(idx_gamma);
      auto const& proton = obj.particles.at(idx_proton);
      double p_gamma  = pvars::gOre::shower_p(gamma);
      double p_proton = pvars::p(proton);
      utilities::three_vector dir_gamma  = utilities::to_three_vector(gamma.start_dir);
      utilities::three_vector dir_proton = utilities::to_three_vector(proton.start_dir);
      double cosTh = utilities::dot_product(dir_gamma, dir_proton);
      mass = std::sqrt(2.*p_gamma*(std::sqrt(p_proton*p_proton + PROTON_MASS*PROTON_MASS) - p_proton*cosTh) + PROTON_MASS*PROTON_MASS);
      return mass;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, delta_mass, delta_mass);

  /**
   * @brief approximate the Delta baryon mass (1g1p)
   * @details Use the proton
   **/
  template <class T>
    double delta_mass_P(const T& obj)
    {
      double mass = -999.9;
      size_t idx_gamma = selectors::gOre::leading_primary_gOre(obj);
      size_t idx_proton = selectors::gOre::leading_primary_proton(obj);
      if (idx_gamma == kNoMatch || idx_proton == kNoMatch)
        return mass;
      auto const& gamma  = obj.particles.at(idx_gamma);
      auto const& proton = obj.particles.at(idx_proton);
      double p_gamma  = pvars::gOre::shower_p(gamma);
      double p_proton = pvars::p(proton);
      utilities::three_vector dir_gamma  = utilities::to_three_vector(gamma.start_dir);
      utilities::three_vector dir_proton = utilities::to_three_vector(proton.start_dir);
      double cosTh = utilities::dot_product(dir_gamma, dir_proton);
      mass = std::sqrt(2.*p_gamma*(std::sqrt(p_proton*p_proton + PROTON_MASS*PROTON_MASS) - p_proton*cosTh) + PROTON_MASS*PROTON_MASS);
      return mass;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, delta_mass_P, delta_mass_P);

  /**
   * @brief approximate the Delta baryon mass (1g0p)
   * @details Make approximations for the neutron. Params are first the scaling of the PE to neutron momentum,
   * and second the scaling factor for transverse direction cosine to cosTh
   **/
  template <class T>
    double delta_mass_N(const T& obj, std::vector<double> params = {0.01, 1.0})
    {
      double mass = -999.9;
      size_t idx_gamma = selectors::gOre::leading_primary_gOre(obj);
      if (idx_gamma == kNoMatch)
        return mass;
      auto const& gamma  = obj.particles.at(idx_gamma);
      double p_gamma   = pvars::gOre::shower_p(gamma);
      double p_neutron = params[0]*vars::flash_total_pe(obj);
      double cosTh = params[1]*std::sqrt(1 - std::pow(pvars::gOre::shower_dpT(gamma) / p_gamma, 2.0));
      mass = std::sqrt(2.*p_gamma*(std::sqrt(p_neutron*p_neutron + NEUTRON_MASS*NEUTRON_MASS) - p_neutron*cosTh) + NEUTRON_MASS*NEUTRON_MASS);
      return mass;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Both, delta_mass_N, delta_mass_N);

  //*** TRUTH ONLY VARS ***//

  /**
   * @brief The Delta baryon mass square splitting for 1g0p
   * @tparam T The type of the interaction (MC Truth only)
   * @param obj the interaction
   * @return double the Delta baryon mass
   **/
  template <class T>
    double delta_mass_N_MC(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY,
                                                                        GORE_FID_THRESH_X_POS, GORE_FID_THRESH_X_NEG, GORE_FID_THRESH_Y_POS, GORE_FID_THRESH_Y_NEG, GORE_FID_THRESH_Z_POS, GORE_FID_THRESH_Z_NEG})
    {
      double mass = kNoMatchValue;
      bool isnc = obj.isnc;
      caf::genie_interaction_mode_ genie_mode = obj.genie_mode;
      caf::genie_interaction_type_ genie_inttype = obj.genie_inttype;
      int resnum = obj.resnum;
      // is NC ∆ res
      bool is_nc_delta_res = isnc && (resnum == 0);
      if (not is_nc_delta_res)
        return mass;
      // post-FSI primary particles
      core::gOre::mc_topology topology(obj.prim, params);
      // single photon topology (1γ and maybe some nucleons)
      // here want only 1γ0p
      bool is_single_photon_topology = topology.single_photon() && topology.only_photons_and_nucleons();
      bool no_protons = (topology.count_with_antiparticles(2212) == 0);
      if ((not is_single_photon_topology) || (not no_protons))
        return mass;

      // take the first neutron, which is probably the right one
      utilities::three_vector p_gamma_vec   = topology.get(  22, 0).momentum();
      utilities::three_vector p_neutron_vec = topology.get(2112, 0).momentum();
      double p_gamma   = utilities::magnitude(p_gamma_vec);
      double p_neutron = utilities::magnitude(p_neutron_vec);
      double cosTh = utilities::dot_product(p_gamma_vec, p_neutron_vec) / (p_gamma * p_neutron);
     
      // what is the delta mass?
      mass = std::sqrt(2.*p_gamma*(std::sqrt(p_neutron*p_neutron + NEUTRON_MASS*NEUTRON_MASS) - p_neutron*cosTh) + NEUTRON_MASS*NEUTRON_MASS);
      return mass; 
    }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, delta_mass_N_MC, delta_mass_N_MC);

  /**
   * @brief The neutron momentum in MC truth
   * @tparam T The type of the interaction (MC Truth only)
   * @param obj the interaction
   * @return double the Delta neutron momentum
   **/
  template <class T>
    double neutron_momentum(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY,
                                                                        GORE_FID_THRESH_X_POS, GORE_FID_THRESH_X_NEG, GORE_FID_THRESH_Y_POS, GORE_FID_THRESH_Y_NEG, GORE_FID_THRESH_Z_POS, GORE_FID_THRESH_Z_NEG})
    {
      double momentum = kNoMatchValue;
      bool isnc = obj.isnc;
      caf::genie_interaction_mode_ genie_mode = obj.genie_mode;
      caf::genie_interaction_type_ genie_inttype = obj.genie_inttype;
      int resnum = obj.resnum;
      // is NC ∆ res
      bool is_nc_delta_res = isnc && (resnum == 0);
      if (not is_nc_delta_res)
        return momentum;
      // post-FSI primary particles
      core::gOre::mc_topology topology(obj.prim, params);
      // single photon topology (1γ and maybe some nucleons)
      // here want only 1γ0p
      bool is_single_photon_topology = topology.single_photon() && topology.only_photons_and_nucleons();
      bool no_protons = (topology.count_with_antiparticles(2212) == 0);
      if ((not is_single_photon_topology) || (not no_protons))
        return momentum;

      // take the first neutron, which is probably the right one
      utilities::three_vector p_neutron_vec = topology.get(2112, 0).momentum();
      momentum = utilities::magnitude(p_neutron_vec);
      return momentum;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, neutron_momentum, neutron_momentum);
    
  /**
   * @brief get the resonance number from the MC Truth
   * @details from genie::EResonance enum (I think)
   *  kNoResonance = -1,
   *  kP33_1232    =  0,
   *  kS11_1535    =  1,
   *  kD13_1520    =  2,
   *  kS11_1650    =  3,
   *  kD13_1700    =  4,
   *  kD15_1675    =  5,
   *  kS31_1620    =  6,
   *  kD33_1700    =  7,
   *  kP11_1440    =  8,
   *  kP33_1600    =  9,
   *  kP13_1720    = 10,
   *  kF15_1680    = 11,
   *  kP31_1910    = 12,
   *  kP33_1920    = 13,
   *  kF35_1905    = 14,
   *  kF37_1950    = 15,
   *  kP11_1710    = 16,
   *  kF17_1970    = 17
   * @param obj the interaction of interest (MCTruth only)
   * @return the resonance number (straight from GENIE)
   **/
  template <class T>
    double baryon_res_code(const T& obj)
    {
      return obj.resnum;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, baryon_res_code, baryon_res_code);


  /**
   * @brief categorize for MC Truth
   * @details The fiducialization is being a pain, so focus on FSI topology and mode/interaction tyoe
   * 0: NC ∆->Nγ
   * 1: NC Other Single Photon
   * 2: NC π0 (∆ Res)
   * 3: NC π0 (Other)
   * 4: NC Charged π
   * 5: Other NC
   * 6: CC-Electron
   * 7: Other CC
   * 8: Non-neutrino
   **/
  template <class T>
    double mc_category(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY,
                                                                   GORE_FID_THRESH_X_POS, GORE_FID_THRESH_X_NEG, GORE_FID_THRESH_Y_POS, GORE_FID_THRESH_Y_NEG, GORE_FID_THRESH_Z_POS, GORE_FID_THRESH_Z_NEG})
    {
      double cat(8);
      bool isnc = obj.isnc;
      caf::genie_interaction_mode_ genie_mode = obj.genie_mode;
      caf::genie_interaction_type_ genie_inttype = obj.genie_inttype;
      int resnum = obj.resnum;
      // post-FSI primary particles
      core::gOre::mc_topology topology(obj.prim, params);
      // is NC ∆ res
      bool is_nc_delta_res = isnc && (resnum == 0);
      // single photon topology (1γ and maybe some nucleons)
      bool is_single_photon_topology = topology.single_photon() && topology.only_photons_and_nucleons();
      // is a neutrino interaction
      bool is_nu = obj.index != -1;
      if (is_nc_delta_res && is_single_photon_topology)
        cat = 0;
      else if (isnc && is_single_photon_topology)
        cat = 1;
      else if (is_nc_delta_res && topology.has_pi0())
        cat = 2;
      else if (isnc && topology.has_pi0())
        cat = 3;
      else if (isnc && topology.has_pi_pm())
        cat = 4;
      else if (isnc)
        cat = 5;
      else if (not isnc && topology.count_with_antiparticles(11) > 0)
        cat = 6;
      else if (is_nu)
        cat = 7;
      return cat;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, mc_category, mc_category);
} //end vars::gOre namespace

#endif

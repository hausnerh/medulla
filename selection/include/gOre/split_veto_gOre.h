/**
 * @file split_veto_gOre.h
 * @brief Split-π0 diagnostics and veto for the single photon selection
 * @details A π0 whose second photon is reconstructed as a SEPARATE interaction
 * looks like a single photon candidate. These per-candidate variables look at
 * the event's OTHER reco interactions (via context::current_event) for a shower
 * that could be that missing photon: how close it starts to the candidate
 * vertex, whether it points back to it, and the π0 mass it makes with the
 * candidate photon. Other interactions are included whether or not they are
 * flash-matched (a split-off shower rarely gets a clean flash match). In ICARUS
 * only showers in the candidate's cryostat (same sign of x) are considered.
 * @author hhausner@fnal.gov
 **/

#ifndef SPLIT_VETO_GORE_H
#define SPLIT_VETO_GORE_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "include/gOre/vars_gOre.h"

#define GORE_SPLIT_DIST_MAX   200.0 // cm, search radius around the candidate vertex
#define GORE_SPLIT_IMPACT_MAX  10.0 // cm, placeholder veto threshold (tune on N-1)
#define GORE_PI0_MASS         134.98 // MeV

/**
 * @namespace vars::gOre
 * @brief Variables specific to NC single photon analyses
 **/
namespace vars::gOre
{
  /**
   * @brief Summary of the other-interaction showers around a candidate
   **/
  struct split_shower_summary
  {
    double dist          = std::numeric_limits<double>::quiet_NaN(); // closest other-shower start to the vertex [cm]
    double p             = std::numeric_limits<double>::quiet_NaN(); // momentum of that closest shower [MeV]
    double impact        = std::numeric_limits<double>::quiet_NaN(); // min impact parameter to the vertex, showers within dist_max pointing away [cm]
    double pi0_mass_diff = std::numeric_limits<double>::quiet_NaN(); // min |1 - m(γ_cand, γ_other) / m_π0|, showers within dist_max
  };

  /**
   * @brief Scan the event's other reco interactions for showers near the candidate
   * @param obj the candidate reco interaction
   * @param params [0] minimum shower momentum (MeV), [1] search radius (cm)
   **/
  template <class T>
    split_shower_summary split_shower_scan(const T& obj, const std::vector<double>& params)
    {
      split_shower_summary out;
      const EventType * sr = context::current_event;
      if (sr == nullptr)
        return out;
      double p_min = (params.size() > 0) ? params.at(0) : GORE_MIN_GORE_ENERGY;
      double d_max = (params.size() > 1) ? params.at(1) : GORE_SPLIT_DIST_MAX;
      bool   icarus = (context::current_detector == caf::Det_t::kICARUS);

      utilities::three_vector vtx = utilities::to_three_vector(obj.vertex);

      // candidate photon, for the π0 mass with the other shower
      size_t idx_gamma = selectors::gOre::leading_primary_gOre(obj);
      double p_gamma = std::numeric_limits<double>::quiet_NaN();
      utilities::three_vector dir_gamma;
      if (idx_gamma != kNoMatch)
      {
        auto const& gamma = obj.particles.at(idx_gamma);
        p_gamma   = pvars::gOre::shower_p(gamma);
        dir_gamma = utilities::normalize(utilities::to_three_vector(gamma.start_dir));
      }

      double best_dist(std::numeric_limits<double>::max());
      double best_impact(std::numeric_limits<double>::max());
      double best_mdiff(std::numeric_limits<double>::max());
      for (auto const& other : sr->dlp)
      {
        if ((int64_t)other.id == (int64_t)obj.id)
          continue;
        for (auto const& p : other.particles)
        {
          if (p.pid != pvars::kPhoton && p.pid != pvars::kElectron)
            continue;
          double p_other = pvars::gOre::shower_p(p);
          if (!(p_other >= p_min))
            continue;
          utilities::three_vector start = utilities::to_three_vector(p.start_point);
          if (icarus && (std::get<0>(start) * std::get<0>(vtx) < 0))
            continue; // other cryostat
          utilities::three_vector sep = utilities::subtract(start, vtx);
          double dist = utilities::magnitude(sep);
          if (dist < best_dist)
          {
            best_dist = dist;
            out.p = p_other;
          }
          if (dist > d_max)
            continue;

          // impact parameter of the shower's backward-projected axis to the
          // vertex, for showers travelling away from the vertex
          utilities::three_vector dir = utilities::normalize(utilities::to_three_vector(p.start_dir));
          if (utilities::dot_product(dir, sep) > 0)
          {
            double impact = utilities::magnitude(utilities::cross_product(sep, dir));
            if (impact < best_impact)
              best_impact = impact;
          }

          // π0 mass, taking the other photon to come from the candidate vertex
          if (idx_gamma != kNoMatch && dist > 0)
          {
            double cos_th = utilities::dot_product(dir_gamma, utilities::normalize(sep));
            double mass = std::sqrt(std::max(0., 2. * p_gamma * p_other * (1. - cos_th)));
            double mdiff = std::abs(1. - mass / GORE_PI0_MASS);
            if (mdiff < best_mdiff)
              best_mdiff = mdiff;
          }
        }
      }
      if (best_dist   < std::numeric_limits<double>::max()) out.dist          = best_dist;
      if (best_impact < std::numeric_limits<double>::max()) out.impact        = best_impact;
      if (best_mdiff  < std::numeric_limits<double>::max()) out.pi0_mass_diff = best_mdiff;
      return out;
    }

  /**
   * @brief Distance (cm) from the candidate vertex to the closest shower start in another interaction
   * @param params [0] minimum shower momentum (MeV)
   * @return NaN if there is no such shower
   **/
  template <class T>
    double split_shower_dist(const T& obj, const std::vector<double>& params)
    {
      return split_shower_scan(obj, params).dist;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Reco, gOre_split_shower_dist, split_shower_dist);

  /**
   * @brief Momentum (MeV) of the closest shower in another interaction
   * @param params [0] minimum shower momentum (MeV)
   * @return NaN if there is no such shower
   **/
  template <class T>
    double split_shower_p(const T& obj, const std::vector<double>& params)
    {
      return split_shower_scan(obj, params).p;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Reco, gOre_split_shower_p, split_shower_p);

  /**
   * @brief Smallest impact parameter (cm) to the candidate vertex of another
   * interaction's shower starting within the search radius and pointing away
   * @param params [0] minimum shower momentum (MeV), [1] search radius (cm)
   * @return NaN if there is no such shower
   **/
  template <class T>
    double split_shower_impact(const T& obj, const std::vector<double>& params)
    {
      return split_shower_scan(obj, params).impact;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Reco, gOre_split_shower_impact, split_shower_impact);

  /**
   * @brief Best |1 - m/m_π0| of the candidate photon with another interaction's
   * shower within the search radius (that shower's direction taken from the vertex)
   * @param params [0] minimum shower momentum (MeV), [1] search radius (cm)
   * @return NaN if there is no candidate photon or no such shower
   **/
  template <class T>
    double split_pi0_mass_diff(const T& obj, const std::vector<double>& params)
    {
      return split_shower_scan(obj, params).pi0_mass_diff;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Reco, gOre_split_pi0_mass_diff, split_pi0_mass_diff);

  /**
   * @brief Truth diagnostic: how many OTHER reco interactions are truth-matched
   * to the same neutrino as the candidate (i.e. split fragments of it)
   * @return NaN for data, or if the candidate is not matched to a neutrino
   **/
  template <class T>
    double split_same_nu_fragments(const T& obj)
    {
      const EventType * sr = context::current_event;
      double nan = std::numeric_limits<double>::quiet_NaN();
      if (sr == nullptr || !sr->hdr.ismc)
        return nan;
      auto nu_of = [sr](const auto& reco) -> int64_t {
        if (reco.match_ids.size() == 0)
          return -1;
        size_t idx = (size_t)reco.match_ids[0];
        return (idx < sr->dlp_true.size()) ? (int64_t)sr->dlp_true[idx].nu_id : -1;
      };
      int64_t nu = nu_of(obj);
      if (nu < 0)
        return nan;
      double count(0);
      for (auto const& other : sr->dlp)
        if ((int64_t)other.id != (int64_t)obj.id && nu_of(other) == nu)
          count += 1;
      return count;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::Reco, gOre_split_same_nu_fragments, split_same_nu_fragments);
} // end vars::gOre namespace

/**
 * @namespace cuts::gOre
 * @brief Cuts specific to NC single photon analyses
 **/
namespace cuts::gOre
{
  /**
   * @brief Split-π0 veto: reject the candidate if another interaction has a
   * shower near the vertex that points back to it (and, optionally, makes a
   * π0 mass with the candidate photon)
   * @param params [0] minimum shower momentum (MeV), [1] search radius (cm),
   * [2] max impact parameter (cm), [3] optional: also veto if |1 - m/m_π0| is
   * below this (<= 0 disables)
   * @return true if the candidate passes (no split-π0-like shower)
   **/
  template <class T>
    bool split_shower_veto(const T& obj, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_SPLIT_DIST_MAX, GORE_SPLIT_IMPACT_MAX, 0.})
    {
      double impact_max = (params.size() > 2) ? params.at(2) : GORE_SPLIT_IMPACT_MAX;
      double mass_frac  = (params.size() > 3) ? params.at(3) : 0.;
      vars::gOre::split_shower_summary s = vars::gOre::split_shower_scan(obj, params);
      if (s.impact < impact_max)
        return false;
      if (mass_frac > 0 && s.pi0_mass_diff < mass_frac)
        return false;
      return true;
    }
  REGISTER_CUT_SCOPE(RegistrationScope::Reco, gOre_split_shower_veto, split_shower_veto);
} // end cuts::gOre namespace

#endif

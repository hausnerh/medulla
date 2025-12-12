/**
 * @file particle_vars_gOre.h
 * @brief particle variables for gOre analyses
 * @author Harry Hausner (hhausner@fnal.gov)
 **/

#ifndef PARTICLE_VARS_GORE_H
#define PARTICLE_VARS_GORE_H

#include "include/particle_variables.h"

/**
 * @namespace pvars::gOre
 * @brief Particle variables specific for gOre analyses
 **/
namespace pvars::gOre
{
  /**
   * @brief Return the ancestor of the selected particle
   * @tparam T the type of particle (True Particle only)
   * @param p the selected particle
   * @return double the pdg code of the ancestor, (NaN if invalid)
   **/
  template <class T>
    double ancestor(const T& p)
    {
      return static_cast<double>(p.ancestor_pdg_code);
    }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, ancestor, ancestor);

  /**
   * @brief Correction factor for shower KE
   * @details There are three factors which go into this correction,
   * 1. the pion mass peak correction
   * 2. removal of an outdated shower correction factor
   * 3. a new shower correction factor to replace the old
   * The pion mass correction is different for data & MC, but just do MC for now
   * @tparam T the type of particle (true or reco)
   * @param p the selected particle
   * @return double the correction factor for the shower kinetic energy
   **/
  template <class T>
    double corr_ke(const T& p)
    {
      double corr = 134.9768;
      corr /= 130.68;
      corr /= 1.2359;
      corr /= 0.82;
      return corr;
    }

  /**
   * @brief Correction factor for shower momentum
   * @details Like the kinetic energy correction, but with a mass factor
   * @tparam T the type of particle (true or reco)
   * @param p the selected particle
   * @return double the correction factor for the shower kinetic energy
   **/
  template <class T>
    double corr_p(const T& p)
    {
      double corr = corr_ke(p);
      double ke = p.calo_ke;
      double m = pvars::mass(p);
      return corr*std::sqrt((ke + 2.*corr*m) / (ke + 2.*corr*m));
    }

  /**
   * @brief Corrected shower KE
   * @tparam T the type of particle (true or reco)
   * @param p the selected particle
   * @return the corrected shower kinetic energy
   **/
  template <class T>
    double shower_ke(const T& p)
    {
      return corr_ke(p) * p.calo_ke;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, shower_ke, shower_ke);

  /**
   * @brief Corrected shower momentum
   * @details use the cottection factors from the KE to scale the momentum
   * In the case where it's a photon the factors are applied directly, but there is also a factor for the mass
   * @tparam T the type of particle (true or reco)
   * @param p the selected particle
   * @return the corrected shower momentum
   **/
  template <class T>
    double shower_p(const T& p)
    {
      return corr_p(p) * pvars::p(p);
    }
  REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, shower_p, shower_p);

  /**
   * @brief Corrected shower dpT
   * @details use the cottection factors from the KE to scale the momentum
   * In the case where it's a photon the factors are applied directly, but there is also a factor for the mass
   * @tparam T the type of particle (true or reco)
   * @param p the selected particle
   * @return the corrected shower transverse momentum
   **/
  template <class T>
    double shower_dpT(const T& p)
    {
      return corr_p(p) * pvars::dpT(p);
    }
  REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, shower_dpT, shower_dpT);
}// end pvars::gOre namespace

#endif
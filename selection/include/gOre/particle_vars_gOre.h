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
   * @brief Return the axial spread of the particle
   * @tparam T the type of particle (True or Reco) 
   * @param p the selected particle
   * @return double the axial spread of the particle
   **/
  template <class T>
    double axial_spread(const T& p)
    {
      return p.axial_spread;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, axial_spread, axial_spread);

  /**
   * @brief Return the directional spread of the particle
   * @tparam T the type of particle (True or Reco) 
   * @param p the selected particle
   * @return double the directional spread of the particle
   **/
  template <class T>
    double directional_spread(const T& p)
    {
      return p.directional_spread;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, directional_spread, directional_spread);
}// end pvars::gOre namespace

#endif
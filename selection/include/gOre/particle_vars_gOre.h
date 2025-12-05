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
}// end pvars::gOre namespace

#endif
/**
 * @file scorers_gOre.h
 * @brief Particle Scores for the gOre analysis 
 * @author Harry Hausner (hhausner@fnal.gov)
 **/

#ifndef SCORERS_GORE_H
#define SCORERS_GORE_H

#include "include/scorers.h"

/**
 * @namespace scorers::gOre
 * @brief gOre particle selectors
 **/
namespace scorers::gOre
{
  /**
   * @brief Variable for the particle's PID. Looser photon score for gOre.
   * @param p the particle to apply the variable on.
   * @return the PID of the particle.
   */
  template<class T>
    double gOre_pid(const T & p)
    {
      double pid = std::numeric_limits<double>::quiet_NaN();
      if (p.pid == pvars::kPhoton || p.pid == pvars::kElectron)
      {
        pid = (p.pid_scores[pvars::kPhoton] > 0.077) ? pvars::kPhoton : pvars::kElectron;
      }
      else
      {
        pid = pvars::default_pid(p);
      }
      return pid;
    }
  REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, gOre_pid, gOre_pid);
}

#endif
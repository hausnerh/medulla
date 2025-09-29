/**
 * @file selectors_gOre.h
 * @brief Particle Selectors for the gOre analysis 
 * @author Harry Hausner (hhausner@fnal.gov)
 **/

#ifndef SELECTORS_GORE_H
#define SELECTORS_GORE_H

#include "include/selectors.h"

/**
 * @namespace selectors::gOre
 * @brief gOre particle selectors
 **/
namespace selectors::gOre
{
  /**
   * @brief What is the index of the leading gOre?
   * @tparam T the type of interaction (true or reco).
   * @param obj the interaction to operate on.
   * @return size_t the index of the gOre with the highest kinetic energy
   **/
  template <class T>
    size_t leading_gOre(const T & obj)
    {
      double leading_ke(0);
      size_t index(kNoMatch);
      for(size_t i(0); i < obj.particles.size(); ++i)
      {
        const auto & p = obj.particles.at(i);
        double energy(pvars::ke(p));
        bool is_gOre = (pvars::pid(p) == pvars::kPhoton) || (pvars::pid(p) == pvars::kElectron);
        if (is_gOre && energy > leading_ke)
        {
          leading_ke = energy;
          index = i;
        }
      }
      return index;
    }
  REGISTER_SELECTOR(leading_gOre, leading_gOre);
  
  /**
   * @brief What is the index of the subleading gOre?
   * @tparam T the type of interaction (true or reco).
   * @param obj the interaction to operate on.
   * @return size_t the index of the gOre with the second highest kinetic energy
   **/
  template <class T>
    size_t subleading_gOre(const T & obj)
    {
      double leading_ke(0);
      size_t leading_index(kNoMatch);
      double subleading_ke(0);
      size_t subleading_index(kNoMatch);
      for(size_t i(0); i < obj.particles.size(); ++i)
      {
        const auto & p = obj.particles.at(i);
        double energy(pvars::ke(p));
        bool is_gOre = (pvars::pid(p) == pvars::kPhoton) || (pvars::pid(p) == pvars::kElectron);
        if (is_gOre && energy > leading_ke)
        {
          subleading_ke = leading_ke;
          subleading_index = leading_index;
          leading_ke = energy;
          leading_index = i;
        } else if (energy > subleading_ke)
        {
          subleading_ke = energy;
          subleading_index = i;
        }
      }
      return subleading_index;
    }
  REGISTER_SELECTOR(subleading_gOre, subleading_gOre);
}// end namespace selectors::gOre

#endif

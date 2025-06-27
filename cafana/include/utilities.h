/**
 * @file utilities.h
 * @brief Header file for definitions of utility functions acting on
 * interactions.
 * @details This file contains definitions of utility functions which are used
 * to support the implementation of analysis variables and cuts. These functions
 * are intended to be used to simplify the implementation of variables and cuts
 * by providing common functionality which can be reused across multiple
 * variables and cuts.
 * @author mueller@fnal.gov
 */
#ifndef UTILITIES_H
#define UTILITIES_H

#include <vector>

#include "include/particle_variables.h"
#include "include/particle_cuts.h"

/**
 * @namespace utilities
 * @brief Namespace for organizing utility functions for supporting analysis
 * variables and cuts.
 * @details This namespace is intended to be used for organizing utility
 * functions which are used to support the implementation of analysis variables
 * and cuts. These functions are intended to be used to simplify the
 * implementation of variables and cuts by providing common functionality which
 * can be reused across multiple variables and cuts.
 * @note The namespace is intended to be used in conjunction with the
 * vars and cuts namespaces, which are used for organizing variables and cuts
 * which act on interactions.
 */
namespace utilities
{
    /**
     * @brief Count the primaries of the interaction with cuts applied to each particle.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to find the topology of.
     * @return the count of primaries of each particle type within the
     * interaction.
     */
    template<class T>
        std::vector<uint32_t> count_primaries(const T & obj)
        {
            std::vector<uint32_t> counts(5, 0);
            for(auto &p : obj.particles)
            {
                if(pcuts::final_state_signal(p))
                    ++counts[PIDFUNC(p)];
            }
            return counts;
        }

    /**
     * @brief Finds the index corresponding to the leading particle of the specifed
     * particle type.
     * @details The leading particle is defined as the particle with the highest
     * kinetic energy. If the interaction is a true interaction, the initial kinetic
     * energy is used instead of the CSDA kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @param pid of the particle type.
     * @return the index of the leading particle (highest KE). 
     */
    template <class T>
        size_t leading_particle_index(const T & obj, uint16_t pid)
        {
            double leading_ke(0);
            size_t index(0);
            for(size_t i(0); i < obj.particles.size(); ++i)
            {
                const auto & p = obj.particles[i];
                double energy(p.csda_ke);
                if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                    energy = pvars::ke(p);
                if(PIDFUNC(p) == pid && energy > leading_ke)
                {
                    leading_ke = energy;
                    index = i;
                }
            }
            return index;
        }

    /**
     * @brief Finds the index corresponding to the leading muon.
     * @details The leading muon is defined as the muon with the highest
     * kinetic energy. If the interaction is a true interaction, the initial
     * kinetic energy is used instead of the CSDA kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading muon (highest KE).
     */
    template<class T>
        size_t leading_muon_index(const T & obj)
        {
            return leading_particle_index(obj, 2);
        }
    
    /**
     * @brief Finds the index corresponding to the leading proton.
     * @details The leading proton is defined as the proton with the highest
     * kinetic energy. If the interaction is a true interaction, the initial
     * kinetic energy is used instead of the CSDA kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading proton (highest KE).
     */
    template<class T>
        size_t leading_proton_index(const T & obj)
        {
            return leading_particle_index(obj, 4);
        }

  /**
   * @brief get the particle vector for an interaction
   * @tparam TP the type of particles (true or reco).
   * @tparam T the type of interaction (true or reco).
   * @param obj the interaction to operate on.
   * @return the vector of primary particles above threshold in the interaction
   **/
  template<class TP, class T>
    std::vector<const TP*> get_all_particles(const T& obj)
    {
      std::vector<const TP*> retParticles;
      for (auto const& p : obj.particles)
        retParticles.push_back(&p);
      return retParticles;
    }

  /**
   * @brief get the indices of the particles for the specified type
   * @tparam T the type of interaction (true or reco).
   * @param obj the interaction to operate on.
   * @param id the type of particle you want (see inclue/particle_utilities.h)
   * @return a vector of indices coresponding to the primary particles of the given type
   **/
  template<class T>
    std::vector<unsigned int> get_particle_indices(const T& obj, const pvars::Particle_t id)
    {
      std::vector<unsigned int> indices;
      if      constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
      {
        std::vector<const caf::SRParticleTruthDLPProxy*> particles = get_all_particles<caf::SRParticleTruthDLPProxy>(obj);
        for (unsigned int idx = 0; idx < particles.size(); ++idx)
          if (pcuts::final_state_signal(*particles[idx]) && (PIDFUNC(*particles[idx]) == id))
            indices.push_back(idx);
      }
      else if constexpr (std::is_same_v<T, caf::SRInteractionDLPProxy>)
      {
        std::vector<const caf::SRParticleDLPProxy*> particles = get_all_particles<caf::SRParticleDLPProxy>(obj);
        for (unsigned int idx = 0; idx < particles.size(); ++idx)
          if (pcuts::final_state_signal(*particles[idx]) && (PIDFUNC(*particles[idx]) == id))
            indices.push_back(idx);
      }
      return indices;
    }

  /**
   * @brief get the specified particles from the interaction
   * @tparam TP the type of particles (true or reco).
   * @tparam T the type of interaction (true or reco).
   * @param obj the interaction to operate on.
   * @param id the type of particle you want (see inclue/particle_utilities.h)
   * @return the vector of primary particles of specified type above threshold in the interaction
   **/
  template<class TP, class T>
    std::vector<const TP*> get_specified_particles(const T& obj, const pvars::Particle_t id)
    {
      std::vector<const TP*> targetParticles;
      if      constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
      {
        std::vector<const caf::SRParticleTruthDLPProxy*> particles = get_all_particles<caf::SRParticleTruthDLPProxy>(obj);
        for (unsigned int idx = 0; idx < particles.size(); ++idx)
          if (pcuts::final_state_signal(*particles[idx]) && (PIDFUNC(*particles[idx]) == id))
            targetParticles.push_back(particles[idx]);
      }
      else if constexpr (std::is_same_v<T, caf::SRInteractionDLPProxy>)
      {
        std::vector<const caf::SRParticleDLPProxy*> particles = get_all_particles<caf::SRParticleDLPProxy>(obj);
        for (unsigned int idx = 0; idx < particles.size(); ++idx)
          if (pcuts::final_state_signal(*particles[idx]) && (PIDFUNC(*particles[idx]) == id))
            targetParticles.push_back(particles[idx]);
      }
      return targetParticles;
    }
}
#endif // UTILITIES_H

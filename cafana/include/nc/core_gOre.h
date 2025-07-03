/**
 * @file core_gOre.h
 * @brief Header file for defining variables useful for NC analyses
 * specficially for saving photon/electron discrimination for the end
 * @author hhausner@fnal.gov
 **/

#ifndef CORE_GORE_H
#define CORE_GORE_H

/**
 * Threshold levels for particle KE
 * Nominal is 143.425 MeV for muons
 * 50 MeV for protons
 * 25 for everything else
 **/
#define GORE_MIN_GORE_ENERGY   25.0    // MeV
#define GORE_MIN_MUON_ENERGY  143.425  // MeV
#define GORE_MIN_PROTON_ENERGY 50.0    // MeV
#define GORE_MIN_PION_ENERGY   25.0    // MeV

#include "include/utilities.h"

namespace core::nc::gOre
{
  template<class T>
    bool final_state_signal(const T& p,
                            const double gore_min   = GORE_MIN_GORE_ENERGY,
                            const double muon_min   = GORE_MIN_MUON_ENERGY,
                            const double proton_min = GORE_MIN_PROTON_ENERGY,
                            const double pion_min   = GORE_MIN_PION_ENERGY)
    {
      bool is_signal(false);
      if (pcuts::is_primary(p))
      {
        pvars::Particle_t pid = static_cast<pvars::Particle_t>(PIDFUNC(p));
        switch (pid)
        {
          case pvars::kMuon:
            is_signal = (pvars::ke(p) > muon_min);
            break;
          case pvars::kProton:
            is_signal = (pvars::ke(p) > proton_min);
            break;
          case pvars::kPion:
            is_signal = (pvars::ke(p) > pion_min);
            break;
          case pvars::kPhoton:
            is_signal = (pvars::ke(p) > gore_min);
            break;
          case pvars::kElectron:
            is_signal = (pvars::ke(p) > gore_min);
            break;
          default:
            break;
        }
      }
      return is_signal;
    }

  template<class T>
    std::vector<size_t> count_primaries(const T& obj,
                                        const double gore_min   = GORE_MIN_GORE_ENERGY,
                                        const double muon_min   = GORE_MIN_MUON_ENERGY,
                                        const double proton_min = GORE_MIN_PROTON_ENERGY,
                                        const double pion_min   = GORE_MIN_PION_ENERGY)
    {
      std::vector<size_t> primaries(5, 0); // five valid pvars::Paricle_t
      for (auto const& particle : obj.particles)
      {
        if (final_state_signal(particle, gore_min, muon_min, proton_min, pion_min))
          ++primaries.at(PIDFUNC(particle));
      }
      return primaries;
    }

  template<class>
  struct ParticleType;
  template<>
    struct ParticleType<caf::SRInteractionDLPProxy>
    {
      typedef caf::SRParticleDLPProxy type;
    };
  template<>
    struct ParticleType<caf::SRInteractionTruthDLPProxy>
    {
      typedef caf::SRParticleTruthDLPProxy type;
    };

  /**
   * @brief hold your interaction (MC or data) and the single photon/electron
   **/
  template <class T>
    struct Interaction
    {
      typedef typename ParticleType<T>::type PT;
      Interaction(const T& object) : obj(object)
      {
        bool found_gOre(false);
        bool found_Other(false);
        min_muon_ke = std::numeric_limits<double>::max();
        min_pion_ke = std::numeric_limits<double>::max();
        subleading_gore_ke = std::numeric_limits<double>::lowest();
        double leading_gore_ke = std::numeric_limits<double>::lowest();
        nShowers = 0;
        total_gore_ke = 0;
        utilities::three_vector    leading_p;
        utilities::three_vector subleading_p;
        for (auto& particle : obj.particles)
        {
          pvars::Particle_t pid = static_cast<pvars::Particle_t>(PIDFUNC(particle));
          if (pcuts::is_primary(particle))
          {
            double pKE = pvars::ke(particle);
            switch (pid)
            {
              case pvars::kMuon:
                min_muon_ke = (min_muon_ke > pKE) ? pKE : min_muon_ke;
                break;
              case pvars::kPion:
                min_pion_ke = (min_pion_ke > pKE) ? pKE : min_pion_ke;
                break;
              case pvars::kProton:
                nProtons_subthresh += (pKE < GORE_MIN_PROTON_ENERGY);
              default:
                if (pid == pvars::kPhoton || pid == pvars::kElectron)
                {
                  ++nShowers;
                  total_gore_ke += pKE;
                  if (pKE > leading_gore_ke)
                  {
                    subleading_gore_ke = leading_gore_ke;
                    subleading_p = leading_p;
                    leading_gore_ke = pKE;
                    leading_p = utilities::to_three_vector(particle.momentum);
                  }
                  else if (pKE > subleading_gore_ke)
                  {
                    subleading_gore_ke = pKE;
                    subleading_p = utilities::to_three_vector(particle.momentum);
                  }
                }
                break;
            }
          }
          if (final_state_signal(particle))
          {
            switch (pid)
            {
              case pvars::kProton:
                protons.push_back(&particle);
                break;
              case pvars::kPhoton:
                if (found_gOre)
                {
                  found_Other = true;
                }
                else
                {
                  photon_or_electron = &particle;
                  found_gOre = true;
                }
                break;
              case pvars::kElectron:
                if (found_gOre)
                {
                  found_Other = true;
                }
                else
                {
                  photon_or_electron = &particle;
                  found_gOre = true;
                }
                break;
              default:
                found_Other = true;
                break;
            }
          }
        }
        is_valid = (found_gOre) && (not found_Other);
        nProtons = protons.size();
        min_muon_ke = (min_muon_ke == std::numeric_limits<double>::max()) ? 0. : min_muon_ke;
        min_pion_ke = (min_pion_ke == std::numeric_limits<double>::max()) ? 0. : min_pion_ke;
        if (is_valid && pvars::ke(*photon_or_electron) < subleading_gore_ke)
        {
          std::string msg = "Primary gOre KE is less than subleading gOre KE!!!\n  Primary gOre KE: " + std::to_string(pvars::ke(*photon_or_electron)) + "\n  Leading gOre KE: " + std::to_string(leading_gore_ke) + " (should be the same as above)\n  subleading gOre KE: " + std::to_string(subleading_gore_ke);
          throw std::runtime_error(msg);
        }
        pion_costh = std::numeric_limits<double>::max();
        if (nShowers > 1)
        {
          utilities::three_vector gOre_1_p = utilities::to_three_vector(photon_or_electron->momentum);
          utilities::three_vector gOre_2_p = subleading_p;
          pion_costh = utilities::dot_product(gOre_1_p, gOre_2_p) / (utilities::magnitude(gOre_1_p) * utilities::magnitude(gOre_2_p));
        }
      }
      const T& obj;
      PT* photon_or_electron;
      std::vector<PT*> protons;
      double min_muon_ke; // lowest muon KE below threshold
      double min_pion_ke; // lowest pion KE below threshold
      double subleading_gore_ke; // what is the KE of the subleading gOre candidate?
      double pion_costh; // if there is a subleading gOre candidate, what is the cosine of the angle between it and the primary?
      double total_gore_ke; // what is the sum of KE for all gOre showers (above and below theshold)?
      bool is_valid; // are there no pions or muons above threshold?
      size_t nProtons; // how many protons are there?
      size_t nProtons_subthresh; // how many protons are there below threshold?
      size_t nShowers; // how many gOre showers are there (should be at least 1 above theshold)
    };

  typedef Interaction<caf::SRInteractionDLPProxy>      Reco_Interaction;
  typedef Interaction<caf::SRInteractionTruthDLPProxy> True_Interaction;
} // end core::nc::gOre namespace

#endif

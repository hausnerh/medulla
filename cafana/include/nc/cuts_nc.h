/**
 * @file cuts_nc.h
 * @brief Header file for defining cuts useful for NC analyses
 * @details This file defines true and reco selection criteria (cuts)
 * for use in analyses using NC interactions. These cuts scaffold off
 * those in include/cuts.h
 * @author hhausner@fnal.gov
 **/

#ifndef CUTS_NC_H
#define CUTS_NC_H
#include "include/utilities.h"
#include "include/nc/vars_nc.h"

/**
 * @namespace cuts::nc
 * @brief Namespace containing cuts useful for NC analyses
 * @details This namespace is intended to be used for organizing cuts which act
 * on interactions specific to NC analyses. Each cut is implemented as
 * a function which takes an interaction object as an argument and returns a
 * boolean. The function should be templated on the type of interaction object if
 * the cut is intended to be used on both true and reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the cuts
 * namespace, which is used for organizing generic cuts which act on interactions.
 **/

namespace cuts::nc
{
  // ****** PARTICLE CUTS ******

    /**
     * @brief Check if the particle meets final state signal requirements (NC Edition).
     * @details must be primary and have an energy above threshold.
     * Muons must have a length of at least 50 cm (143.425 MeV), protons
     * must have an energy above 50 MeV, and all other particles must have
     * an energy above 25 MeV.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to check.
     * @return true if the particle is a final state signal particle.
     */
    template<class T>
        bool final_state_signal(const T & p)
        {
            bool passes(false);
            if(is_primary(p))
            {
                double energy(pvars::ke(p));
                switch (p.pid)
                {
                  case pvars::kPhoton:
                    passes = (energy > 25.);
                    break;
                  case pvars::kProton:
                    passes = (energy > 50.);
                    break;
                  case pvars::kElectron:
                    passes = (energy > 25.);
                    break;
                  case pvars::kMuon:
                    passes = (energy > 143.425);
                    break;
                  case pvars::kPion:
                    passes = (energy > 25.);
                    break;
                }
            }
            return passes;
        }

  // ****** SELECTION CUTS ******

  /**
   * @brief Check if the interaction has charged lepton primaries
   * @details Iterate throught the daughter particles of the interaction, and
   * check if any of the primary daughters are either an election or a muon.
   * @param obj the interaction to examine (MC or data)
   * @return true if there are no charged lepton primaries, false otherwise
   **/
   template<class T>
     bool no_charged_leptons(const T& obj)
     {
       for (auto const& daughter : obj.particles)
         if (is_primary(daughter) && ((daughter.pid == pvars::kElectron) || 
                                      (daughter.pid == pvars::kMuon)     )) return false;
       return true;
     } // end no_charged_leptons

 /**
  * @brief Is the final state topology a signle photon?
  * @details Checks if the final state is exactly one photon above threshold
  * @param obj the interaction to examine (MC or data)
  * @return true if there is exactly 1 photon above threshold in the final state, false otherwise
  **/
 template<class T>
   bool topology_single_photon(const T& obj)
   {
     std::vector<uint32_t> topology = utilities::count_primaries(obj);
     return (topology[pvars::kPhoton]   == 1) &&
            (topology[pvars::kElectron] == 0) &&
            (topology[pvars::kMuon]     == 0) &&
            (topology[pvars::kPion]     == 0) &&
            (topology[pvars::kProton]   == 0) ;
   } // end topology_single_photon

 /**
  * @brief Is the final state topology inclusive of a signle photon?
  * @details Checks if the final state is one photon + other particles
  * @param obj the interaction to examine (MC or data)
  * @return true if there is exactly 1 photon above threshold in the final state + other particles, false otherwise
  **/
 template<class T>
   bool topology_single_photon_inclusive(const T& obj)
   {
     std::vector<uint32_t> topology = utilities::count_primaries(obj);
     return (topology[pvars::kPhoton] == 1);
   } // end topology_single_photon

 /**
  * @brief Is the final state topology a signle photon + a single proton?
  * @details Checks if the final state is exactly one photon & one proton above threshold
  * @param obj the interaction to examine (MC or data)
  * @return true if there is exactly 1 photon & 1 proton in the final state, false otherwise
  **/
 template<class T>
   bool topology_1photon_1proton(const T& obj)
   {
     std::vector<uint32_t> topology = utilities::count_primaries(obj);
     return (topology[pvars::kPhoton]   == 1) &&
            (topology[pvars::kElectron] == 0) &&
            (topology[pvars::kMuon]     == 0) &&
            (topology[pvars::kPion]     == 0) &&
            (topology[pvars::kProton]   == 1) ;
   } // end topology_1photon_1proton

  /**
   * @brief apply a no lepton topology cut based on the topology string
   * @details instead of counting primaries by hand, check the topology string
   * in the event. Do some regex magic.
   * @param obj the interaction of interest (MC or data)
   * @return true if no leptons exist in the topology string
   **/
  template<class T>
    bool no_charged_leptons_string(const T& obj)
    {
      return (vars::nc::count_electrons(obj) == 0) &&
             (vars::nc::count_muons    (obj) == 0) ;
    }

  /**
   * @brief Apply fiducial, containment, & flash cuts to single photon topology
   * @details Logical and of `fiducial_containment_flash_cut` and `topology_single_photon`
   * @param obj the interaction of interest (MC or data)
   * @return true if the interaction is fiducial, contained, passes flash cuts, and has a single photon topology
   **/
  template<class T>
    bool fiducial_containment_flash_cut_single_photon(const T& obj)
    {
      return cuts::fiducial_containment_flash_cut(obj) && topology_single_photon(obj);
    }

  /**
   * @brief Apply fiducial, containment, & flash cuts to inclusive single photon topology
   * @details Logical and of `fiducial_containment_flash_cut` and `topology_single_photon_inclusive`
   * @param obj the interaction of interest (MC or data)
   * @return true if the interaction is fiducial, contained, passes flash cuts, and has
   * an inclusive single photon topology
   **/
  template<class T>
    bool fiducial_containment_flash_cut_single_photon_inclusive(const T& obj)
    {
      return cuts::fiducial_containment_flash_cut(obj) && topology_single_photon_inclusive(obj);
    }

  /**
   * @brief Apply fiducial, containment, & flash cuts to 1photon1proton topology
   * @details Logical and of `fiducial_containment_flash_cut` and `topology_1photon_1proton`
   * @param obj the interaction of interest (MC or data)
   * @return true if the interaction is fiducial, contained, passes flash cuts, and has a 1photon1proton topology
   **/
  template<class T>
    bool fiducial_containment_flash_cut_1photon_1proton(const T& obj)
    {
      return cuts::fiducial_containment_flash_cut(obj) && topology_1photon_1proton(obj);
    }

  // ****** TRUTH CUTS ******
  
  /**
   * @brief Check if the MC interaction is NC
   * @details This checks if the true interaction is NC,
   * which is not exactly not CC due to internal enum issues,
   * ie obj.current_type could be -1 (unknown type)
   * @param obj the interaction (MC only)
   * @return true if the interaction is NC, otherwise false
   **/
  bool isnc(const caf::SRInteractionTruthDLPProxy& obj) { return obj.current_type == 1; }

  /**
   * @brief Check if MC interaction is NC neutrino
   * @details This simply combines the cuts::nc::isnc and cuts::neutrino selections
   * @param obj the interaction (MC only)
   * @return true if obj is true NC neutrino intraction, otherwise false
   **/
  bool isnc_neutrino(const caf::SRInteractionTruthDLPProxy& obj) { return isnc(obj) && cuts::neutrino(obj); }

  /**
   * @brief Check if MC interaction is fiducial contained neutrino
   * @details This is a logical and of `fiducial_cut`, `containment_cut`, & `neutrino`
   * @param obj the interaction (MC only)
   * @return true if obj is a true fiducial contained neutrino event, false otherwise
   **/
  bool is_fid_con_nu(const caf::SRInteractionTruthDLPProxy& obj)
  {
    return cuts::fiducial_cut   (obj) &&
           cuts::containment_cut(obj) &&
           cuts::neutrino       (obj);
  }

  /**
   * @brief Check if MC interaction is fiducial contained NC neutrino
   * @details This is a logical and of `is_fid_con_nu` and `isnc`
   * @param obj the interaction (MC only)
   * @return true if obj is a true fiducial contained NC neutrino event, false otherwise
   **/
  bool is_fid_con_nc_nu(const caf::SRInteractionTruthDLPProxy& obj)
  {
    return is_fid_con_nu(obj) && isnc(obj);
  }

  /**
   * @brief Check if MC interaction is fiducial contained CC neutrino
   * @details This is a logical and of `is_fid_con_nu` and `iscc`
   * @param obj the interaction (MC only)
   * @return true if obj is a true fiducial contained CC neutrino event, false otherwise
   **/
  bool is_fid_con_cc_nu(const caf::SRInteractionTruthDLPProxy& obj)
  {
    return is_fid_con_nu(obj) && cuts::iscc(obj);
  }

  /**
   * @brief Is the true topology is a single photon above threshold?
   * @details We must make sure that there aren't other particles below
   * threshold we are missing in utilities::count_primaries
   * @param obj the interaction to examine (MC only)
   * @return true if there is exactly 1 photon in the final state, false otherwise
   **/
  bool true_topology_single_photon(const caf::SRInteractionTruthDLPProxy& obj)
  {
     unsigned int countPhotons = 0;
     unsigned int countOther = 0;
     for(auto& particle : obj.particles)
     {
       if(not pcuts::is_primary(particle))
         continue;
       switch (particle.pid)
       {
         case pvars::kPhoton:
           if (pcuts::final_state_signal(particle))
             ++countPhotons;
           break;
         default:
           if (pcuts::final_state_signal(particle))
             ++countOther;
           break;
       }
     }
     return (countOther == 0) && (countPhotons == 1);
  }

  /**
   * @brief Is the true topology is a single photon and a signle proton above threshold?
   * @details We must make sure that there aren't other particles below
   * threshold we are missing in utilities::count_primaries
   * @param obj the interaction to examine (MC only)
   * @return true if there is exactly 1 photon & 1 proton in the final state, false otherwise
   **/
   bool true_topology_1photon_1proton(const caf::SRInteractionTruthDLPProxy& obj)
   {
     unsigned int countPhotons = 0;
     unsigned int countProtons = 0;
     unsigned int countOther = 0;
     for(auto& particle : obj.particles)
     {
       if(not pcuts::is_primary(particle))
         continue;
       switch (particle.pid)
       {
         case pvars::kPhoton:
           if (pcuts::final_state_signal(particle))
             ++countPhotons;
           break;
         case pvars::kProton:
           if (pcuts::final_state_signal(particle))
             ++countProtons;
           break;
         default:
           if (pcuts::final_state_signal(particle))
             ++countOther;
           break;
       }
     }
     return (countOther == 0) && (countPhotons == 1) && (countProtons == 1);
   }

  /**
   * @brief Check if MC interaction is fiducial contained single photon event
   * @details This is a logical and of `is_fid_con_nu` & `topology_single_photon`
   * @param obj the interaction (MC only)
   * @return true if obj is a true fiducial contained single photon event
   **/
  bool is_fid_con_single_photon(const caf::SRInteractionTruthDLPProxy& obj)
  {
    return is_fid_con_nu(obj) && topology_single_photon(obj);
  }

  /**
   * @brief Check if MC interaction is fiducial contained NC single photon event
   * @details This is a logical and of `is_fid_con_nc_nu` & `topology_single_photon`
   * @param obj the interaction (MC only)
   * @return true if obj is a true fiducial contained NC single photon event
   **/
  bool is_fid_con_nc_single_photon(const caf::SRInteractionTruthDLPProxy& obj)
  {
    return is_fid_con_nc_nu(obj) && topology_single_photon(obj);
  }

  /**
   * @brief Check if MC interaction is inclusive of fiducial contained single photon event
   * @details This is a logical and of `is_fid_con_nc_nu` & `topology_single_photon`
   * @param obj the interaction (MC only)
   * @return true if obj is a true fiducial contained NC single photon event
   **/
  bool is_fid_con_single_photon_inclusive(const caf::SRInteractionTruthDLPProxy& obj)
  {
    return cuts::fiducial_cut   (obj) &&
           cuts::containment_cut(obj) &&
           cuts::neutrino(obj)        &&
           topology_single_photon_inclusive(obj);
  }

  /**
   * @brief Check if MC interaction is fiducial contained 1photon1proton
   * @details This is a logical and of `is_fid_con_nu` & `topology_1photon_1proton`
   * @param obj the interaction (MC only)
   * @return true if obj is a true fiducial contained 1photon1proton
   **/
  bool is_fid_con_1photon_1proton(const caf::SRInteractionTruthDLPProxy& obj)
  {
    return is_fid_con_nu(obj) && true_topology_1photon_1proton(obj);
  }

  /**
   * @brief Check if MC interaction is fiducial contained NC 1photon1proton
   * @details This is a logical and of `is_fid_con_nc_nu` & `topology_1photon_1proton`
   * @param obj the interaction (MC only)
   * @return true if obj is a true fiducial contained NC 1photon1proton
   **/
  bool is_fid_con_nc_1photon_1proton(const caf::SRInteractionTruthDLPProxy& obj)
  {
    return is_fid_con_nc_nu(obj) && true_topology_1photon_1proton(obj);
  }
 
  /**
   * @brief Check if the MC interaction looks like it contains a single pi0
   * @details We look for two photons (the signal of a pi0 decay)
   * This is a background for NC 1photon serches
   * @param obj the interaction (MC only)
   * @return true if obj is a 1pi0 inclusive topology
   **/
  bool true_topology_pi0(const caf::SRInteractionTruthDLPProxy & obj)
  {
    std::vector<uint32_t> topology = utilities::count_primaries(obj);
    return (topology[pvars::kPhoton] == 2);
  }

  /**
   * @brief Check if the interaction is a true NC pi0
   * @details Is the interaction a NC neutrino with 2 photons?
   * @param obj the interaction (MC only)
   * @return true if obj is a 1pi0 inclusive NC neutrino
   **/
  bool nc_pi0_background(const caf::SRInteractionTruthDLPProxy & obj)
  {
    return is_fid_con_nc_nu(obj) && true_topology_pi0(obj);
  }

  /**
   * @brief Check if the interaction is a true CC pi0
   * @details Is the interaction a CC neutrino with 2 photons?
   * @param obj the interaction (MC only)
   * @return true if obj is a 1pi0 inclusive CC neutrino
   **/
  bool cc_pi0_background(const caf::SRInteractionTruthDLPProxy & obj)
  {
    return cuts::neutrino(obj) && cuts::iscc(obj) && cuts::fiducial_cut(obj) && cuts::containment_cut(obj) && true_topology_pi0(obj);
  }

  /**
   * @brief Check if the interaction is a non-signal 1photon NC interaction
   * @details Event should match our signal apart from the fiducial or containment cuts
   * @param obj the interaction (MC only)
   * @return true if obj is non-fiducial or uncontained NC-1photon event
   **/
  bool nonsignal_single_photon(const caf::SRInteractionTruthDLPProxy& obj)
  {
    return cuts::neutrino(obj) && isnc(obj) && not (cuts::fiducial_cut(obj) && cuts::containment_cut(obj)) && topology_single_photon(obj);
  }

  /**
   * @brief Check if the interaction is a non-signal 1photon1proton NC interaction
   * @details Event should match our signal apart from the fiducial or containment cuts
   * @param obj the interaction (MC only)
   * @return true if obj is non-fiducial or uncontained NC-1photon1proton event
   **/
  bool nonsignal_1photon_1proton(const caf::SRInteractionTruthDLPProxy& obj)
  {
    return cuts::neutrino(obj) && isnc(obj) && not (cuts::fiducial_cut(obj) && cuts::containment_cut(obj)) && topology_1photon_1proton(obj);
  }
} // end namespace cuts::nc
 
// This is for variables that need the cuts
// they belong in the vars::nc namespace but can't go in that file b/c circular references
namespace vars::nc
{ 
  /**
   * @brief Variable for enumerating interaction categories for 1photon selection.
   * @details This variable provides a basic categorization of interactions
   * 0: 1g (contained and fiducial)
   * 1: NC 2g (ie pi0-like) (contained and fiducial)
   * 2: CC 2g (ie pi0-like) (contained and fiducial)
   * 3: 1g (not contained or not fiducial)
   * 4: Other NC nu
   * 5: Other CC nu
   * 6: Cosmic
   * @param obj The interaction to apply the variable on.
   * @return the enumerated category of the interaction.
  */
  double category_single_photon(const caf::SRInteractionTruthDLPProxy& obj)
  {
    if      (cuts::nc::is_fid_con_nc_single_photon(obj)) return 0;
    else if (cuts::nc::nc_pi0_background(obj))           return 1;
    else if (cuts::nc::cc_pi0_background(obj))           return 2;
    else if (cuts::nc::nonsignal_single_photon(obj))     return 3;
    else if (cuts::neutrino(obj) && cuts::nc::isnc(obj)) return 4;
    else if (cuts::neutrino(obj) && cuts::iscc(obj))     return 5;
    else                                                 return 6;
  }

  /**
   * @brief Variable for enumerating interaction categories for inclusive 1photon selection.
   * @details This variable provides a basic categorization of interactions
   * 0: 1g + stuff (contained and fiducial)
   * 1: NC 2g (ie pi0-like) (contained and fiducial)
   * 2: CC 2g (ie pi0-like) (contained and fiducial)
   * 3: 1g (not contained or not fiducial)
   * 4: Other NC nu
   * 5: Other CC nu
   * 6: Cosmic
   * @param obj The interaction to apply the variable on.
   * @return the enumerated category of the interaction.
  */
  double category_single_photon_inclusive(const caf::SRInteractionTruthDLPProxy& obj)
  {
    if      (cuts::nc::is_fid_con_single_photon_inclusive(obj)) return 0;
    else if (cuts::nc::nc_pi0_background(obj))                  return 1;
    else if (cuts::nc::cc_pi0_background(obj))                  return 2;
    else if (cuts::nc::nonsignal_single_photon(obj))            return 3;
    else if (cuts::neutrino(obj) && cuts::nc::isnc(obj))        return 4;
    else if (cuts::neutrino(obj) && cuts::iscc(obj))            return 5;
    else                                                        return 6;
  }

  /**
   * @brief Variable for enumerating interaction categories for 1photon1proton selection.
   * @details This variable provides a basic categorization of interactions
   * 0: 1g1p (contained and fiducial)
   * 1: NC 2g (ie pi0-like) (contained and fiducial)
   * 2: CC 2g (ie pi0-like) (contained and fiducial)
   * 3: 1g1p (not contained or not fiducial)
   * 4: Other NC nu
   * 5: Other CC nu
   * 6: Cosmic
   * @param obj The interaction to apply the variable on.
   * @return the enumerated category of the interaction.
  */
  double category_1photon_1proton(const caf::SRInteractionTruthDLPProxy& obj)
  {
    if      (cuts::nc::is_fid_con_nc_1photon_1proton(obj)) return 0;
    else if (cuts::nc::nc_pi0_background(obj))             return 1;
    else if (cuts::nc::cc_pi0_background(obj))             return 2;
    else if (cuts::nc::nonsignal_1photon_1proton(obj))     return 3;
    else if (cuts::neutrino(obj) && cuts::nc::isnc(obj))   return 4;
    else if (cuts::neutrino(obj) && cuts::iscc(obj))       return 5;
    else                                                   return 6;
  }

  /**
   * @brief Alternative enumeration of interaction categories for 1photon selection.
   * @details This variable provides a basic categorization of interactions
   * 0: NC 1g (contained and fiducial)
   * 1: CC 1g
   * 2: NC 2g (ie pi0-like) (contained and fiducial)
   * 3: CC 2g (ie pi0-like) (contained and fiducial)
   * 4: 1g (not contained or not fiducial)
   * 5: Other NC nu
   * 6: Other CC nu
   * 7: Cosmic
   * @param obj The interaction to apply the variable on.
   * @return the enumerated category of the interaction.
  */
  double alt_category_single_photon(const caf::SRInteractionTruthDLPProxy& obj)
  {
    if      (cuts::nc::is_fid_con_nc_single_photon(obj))                               return 0;
    else if (cuts::nc::is_fid_con_cc_nu(obj) && cuts::nc::topology_single_photon(obj)) return 1;
    else if (cuts::nc::nc_pi0_background(obj))                                         return 2;
    else if (cuts::nc::cc_pi0_background(obj))                                         return 3;
    else if (cuts::nc::nonsignal_single_photon(obj))                                   return 4;
    else if (cuts::neutrino(obj) && cuts::nc::isnc(obj))                               return 5;
    else if (cuts::neutrino(obj) && cuts::iscc(obj))                                   return 6;
    else                                                                               return 7;
  }

  /**
   * @brief Alternative enumeration of interaction categories for inclusive 1photon selection.
   * @details This variable provides a basic categorization of interactions
   * 0: NC 1g + stuff (contained and fiducial)
   * 1: CC 1g + stuff (contained and fiducial)
   * 2: NC 2g (ie pi0-like) (contained and fiducial)
   * 3: CC 2g (ie pi0-like) (contained and fiducial)
   * 4: 1g (not contained or not fiducial)
   * 5: Other NC nu
   * 6: Other CC nu
   * 7: Cosmic
   * @param obj The interaction to apply the variable on.
   * @return the enumerated category of the interaction.
  */
  double alt_category_single_photon_inclusive(const caf::SRInteractionTruthDLPProxy& obj)
  {
    if      (cuts::nc::is_fid_con_single_photon_inclusive(obj))                                  return 0;
    else if (cuts::nc::is_fid_con_cc_nu(obj) && cuts::nc::topology_single_photon_inclusive(obj)) return 1;
    else if (cuts::nc::nc_pi0_background(obj))                                                   return 2;
    else if (cuts::nc::cc_pi0_background(obj))                                                   return 3;
    else if (cuts::nc::nonsignal_single_photon(obj))                                             return 4;
    else if (cuts::neutrino(obj) && cuts::nc::isnc(obj))                                         return 5;
    else if (cuts::neutrino(obj) && cuts::iscc(obj))                                             return 6;
    else                                                                                         return 7;
  }

  /**
   * @brief Alternative enumeration of interaction categories for 1photon1proton selection.
   * @details This variable provides a basic categorization of interactions
   * 0: NC 1g1p (contained and fiducial)
   * 1: CC 1g1p (contained and fiducial)
   * 2: NC 2g (ie pi0-like) (contained and fiducial)
   * 3: CC 2g (ie pi0-like) (contained and fiducial)
   * 4: 1g1p (not contained or not fiducial)
   * 5: Other NC nu
   * 6: Other CC nu
   * 7: Cosmic
   * @param obj The interaction to apply the variable on.
   * @return the enumerated category of the interaction.
  */
  double alt_category_1photon_1proton(const caf::SRInteractionTruthDLPProxy& obj)
  {
    if      (cuts::nc::is_fid_con_nc_1photon_1proton(obj))                               return 0;
    else if (cuts::nc::is_fid_con_cc_nu(obj) && cuts::nc::topology_1photon_1proton(obj)) return 1;
    else if (cuts::nc::nc_pi0_background(obj))                                           return 2;
    else if (cuts::nc::cc_pi0_background(obj))                                           return 3;
    else if (cuts::nc::nonsignal_1photon_1proton(obj))                                   return 4;
    else if (cuts::neutrino(obj) && cuts::nc::isnc(obj))                                 return 5;
    else if (cuts::neutrino(obj) && cuts::iscc(obj))                                     return 6;
    else                                                                                 return 7;
  }
} // end namespace vars::nc
#endif

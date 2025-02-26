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

/**
 * @namespace cuts::nc
 * @brief Namespace containing cuts useful for NC analyses
 * @details This namespace is intended to be used for organizing cuts which act
 * on interactions specific to the muon2024 analysis. Each cut is implemented as
 * a function which takes an interaction object as an argument and returns a
 * boolean. The function should be templated on the type of interaction object if
 * the cut is intended to be used on both true and reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the cuts
 * namespace, which is used for organizing generic cuts which act on interactions.
 **/

namespace cuts::nc
{
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
  * @details Loops over the daughter particles of the interaction, if there is a non-photon
  * primary, or there is more than one photon primary, return false;
  * @param obj the interaction to examine (MC or data)
  * @return true if there is exactly 1 photon in the final state, false otherwise
  **/
 template<class T>
   bool topology_single_photon(const T& obj)
   {
     bool foundPhoton = false;
     for (auto const& daughter : obj.particles)
       if (pcuts::is_primary(daughter))
       {
         if ((daughter.pid != pvars::kPhoton) || (foundPhoton == true))
           return false;
         else
           foundPhoton = true;
       }
     return foundPhoton;
   } // end topology_single_photon

 /**
  * @brief Is the final state topology a signle photon + a single proton?
  * @details Loops over the daughter particles of the interaction
  * @param obj the interaction to examine (MC or data)
  * @return true if there is exactly 1 photon & 1 proton in the final state, false otherwise
  **/
 template<class T>
   bool topology_1photon_1proton(const T& obj)
   {
     bool foundPhoton = false;
     bool foundProton = false;
     for (auto const& daughter : obj.particles)
       switch (daughter.pid)
       {
         case pvars::kPhoton:
           if (foundPhoton) return false;
           foundPhoton = true;
           break;
         case pvars::kProton:
           if (foundProton) return false;
           foundProton = true;
           break;
         default:
           return false;
       }
     return foundPhoton && foundProton;
   } // end topology_1photon_1proton

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
   * @brief Check if MC interaction is fiducial contained NC neutrino
   * @details This is a logical and of `fiducial_cut`, `containment_cut`, & `isnc_neutrino`
   * @param obj the interaction (MC only)
   * @return true if obj is a true fiducial contained NC neutrino event, false otherwise
   **/
  bool is_fid_con_nc_nu(const caf::SRInteractionTruthDLPProxy& obj)
  {
    return cuts::fiducial_cut   (obj) &&
           cuts::containment_cut(obj) &&
           isnc_neutrino        (obj);
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
   * @brief Check if MC interaction is fiducial contained NC 1photon1proton
   * @details This is a logical and of `is_fid_con_nc_nu` & `topology_1photon_1proton`
   * @param obj the interaction (MC only)
   * @return true if obj is a true fiducial contained NC 1photon1proton
   **/
  bool is_fid_con_nc_1photon_1proton(const caf::SRInteractionTruthDLPProxy& obj)
  {
    return is_fid_con_nc_nu(obj) && topology_1photon_1proton(obj);
  }
  
} // end namespace cuts::nc
#endif

/**
 * @file cuts_gOre.h
 * @brief Header file for defining cuts useful for NC analyses
 * specficially for saving photon/electron discrimination for the end
 * @author hhausner@fnal.gov
 **/

#ifndef CUTS_GORE_H
#define CUTS_GORE_H

#include "include/nc/core_gOre.h"

namespace cuts::nc::gOre
{
  template<class T>
    bool gOre_topology(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      return interaction.is_valid;
    }

  template<class T>
    bool gOre_0p(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      return (interaction.is_valid) && (interaction.nProtons == 0);
    }

  template<class T>
    bool gOre_1p(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      return (interaction.is_valid) && (interaction.nProtons == 1);
    }

  template<class T>
    bool gOre_Np(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      return (interaction.is_valid) && (interaction.nProtons > 0);
    }

  template<class T>
    bool fid_con_flash_gOre_0p(const T& obj)
    {
      return cuts::fiducial_containment_flash_cut(obj) && gOre_0p(obj);
    }

  template<class T>
    bool fid_con_flash_gOre_1p(const T& obj)
    {
      return cuts::fiducial_containment_flash_cut(obj) && gOre_1p(obj);
    }

  template<class T>
    bool fid_con_flash_gOre_Np(const T& obj)
    {
      return cuts::fiducial_containment_flash_cut(obj) && gOre_Np(obj);
    }

  template<class T>
    bool has_subthreshold_gOre(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      return (interaction.subleading_gore_ke > 0);
    }

  // *** TRUTH CUTS ***//
  bool is_cc(const caf::SRInteractionTruthDLPProxy& obj) { return obj.current_type == 0; }
  bool is_nc(const caf::SRInteractionTruthDLPProxy& obj) { return obj.current_type == 1; }
  bool is_fid_con_nu(const caf::SRInteractionTruthDLPProxy& obj)
  {
    return cuts::fiducial_cut(obj) && cuts::containment_cut(obj) && cuts::neutrino(obj);
  }
  bool is_fid_con_nc_res(const caf::SRInteractionTruthDLPProxy& obj)
  {
    return is_fid_con_nu(obj) && is_nc(obj) && (obj.interaction_mode == caf::kRes);
  }
  bool gOre_is_photon(const caf::SRInteractionTruthDLPProxy& obj)
  {
    core::nc::gOre::True_Interaction interaction(obj);
    return interaction.is_valid && (PIDFUNC(*interaction.photon_or_electron) == pvars::kPhoton);
  }
  bool gOre_is_electron(const caf::SRInteractionTruthDLPProxy& obj)
  {
    core::nc::gOre::True_Interaction interaction(obj);
    return interaction.is_valid && (PIDFUNC(*interaction.photon_or_electron) == pvars::kElectron);
  }
  bool is_fid_con_gOre(const caf::SRInteractionTruthDLPProxy& obj)
  {
    core::nc::gOre::True_Interaction interaction(obj);
    return interaction.is_valid && is_fid_con_nu(obj);
  }
  /**
   * @brief Check if MC interaction is fiducial contained NC Res event
   * @details This is a logical and of `is_fid_con_nc_nu` and `obj.interaction_mode == caf::kRes`
   * @param obj the interaction (MC only)
   * @return true if obj is a true fiducial contained NC neutrino event, false otherwise
   **/
  bool is_fid_con_nc_nu_res(const caf::SRInteractionTruthDLPProxy& obj)
  {
    return is_fid_con_nu(obj) && is_nc(obj) && (obj.interaction_mode == caf::kRes);
  } 

  //*** MC CUTS***//
  /**
   * @brief is the interaction a true NC Δ res event?
   * @param The MC Truth of the interaction
   **/
  bool nc_delta_res(const caf::Proxy<caf::SRTrueInteraction>& obj)
  {
    return (obj.isnc && (obj.resnum == 0));
  }

  /**
   * @brief is the interaction a true NC Δ res event which isn't catagorized as producing pions?
   * @details check the GENIE interaction type against the different NC Res modes enumerated
   * Would prefer to have the photon producing resonancy to check instead, but w/e
   * @param The MC Truth of the interaction
   **/
  bool nc_delta_res_no_pion(const caf::Proxy<caf::SRTrueInteraction>& obj)
  {
    return nc_delta_res(obj) &&
           (obj.genie_inttype != caf::kResNCNuProtonPi0)       &&
           (obj.genie_inttype != caf::kResNCNuProtonPiPlus)    &&
           (obj.genie_inttype != caf::kResNCNuNeutronPi0)      &&
           (obj.genie_inttype != caf::kResNCNuNeutronPiMinus)  &&
           (obj.genie_inttype != caf::kResNCNuBarProtonPi0)    &&
           (obj.genie_inttype != caf::kResNCNuBarProtonPiPlus) &&
           (obj.genie_inttype != caf::kResNCNuBarNeutronPi0)   &&
           (obj.genie_inttype != caf::kResNCNuBarNeutronPiMinus);
  }

  /**
   * @brief is the interaction a true NC Δ res event which produces a pi0?
   * @details check the GENIE interaction type against the different NC Res modes enumerated
   * Would prefer to have the photon producing resonancy to check instead, but w/e
   * @param The MC Truth of the interaction
   **/
  bool nc_delta_res_pi0(const caf::Proxy<caf::SRTrueInteraction>& obj)
  {
    return nc_delta_res(obj) &&
           ((obj.genie_inttype == caf::kResNCNuProtonPi0)       ||
            (obj.genie_inttype == caf::kResNCNuNeutronPi0)      ||
            (obj.genie_inttype == caf::kResNCNuBarProtonPi0)    ||
            (obj.genie_inttype == caf::kResNCNuBarNeutronPi0)    );
  }
} // end cuts::nc::gOre namespace

 #endif

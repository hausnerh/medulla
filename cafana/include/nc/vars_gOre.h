/**
 * @file vars_gOre.h
 * @brief Header file for defining variables useful for NC analyses
 * specficially for saving photon/electron discrimination for the end
 * @author hhausner@fnal.gov
 **/

#ifndef VARS_GORE_H
#define VARS_GORE_H

#include "include/nc/cuts_gOre.h"

namespace vars::nc::gOre
{
  /**
   * @brief how many protons are in the interaction?
   * @param the interaction of interest (data or MC)
   * @return the number of protons as a double
  **/
  template <class T>
    double n_protons(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return static_cast<double>(interaction.nProtons);
    }

  /**
   * @brief how many protons are in the interaction below threshold?
   * @param the interaction of interest (data or MC)
   * @return the number of protons as a double
  **/
  template <class T>
    double n_protons_subthreshold(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return static_cast<double>(interaction.nProtons_subthresh);
    }

  /**
   * @brief the total KE of protons above threshold
   * @param the interaction of interest (data or MC)
   **/
  template <class T>
    double total_proton_ke(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      double ke = 0;
      for (auto const& proton : interaction.protons)
        ke += proton->ke;
      return ke;
    }

  /**
   * @brief the minimum KE of subthreshold muons
   * @details if there are no subthreshold muons the energy is 0;
   **/
  template <class T>
    double min_muon_ke(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return interaction.min_muon_ke;
    }

  /**
   * @brief the minimum KE of subthreshold pions
   * @details if there are no subthreshold pions the energy is 0;
   **/
  template <class T>
    double min_pion_ke(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return interaction.min_pion_ke;
    }

  /**
   * @brief what is the kinetic energy of the candidate photon?
   * @param the interaction of interest (data or MC)
   * @return the kinetic energy of the photon/electron
  **/
  template <class T>
    double gOre_ke(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return interaction.photon_or_electron->ke;
    }

  /**
   * @brief what is the kinetic energy of the subleading candidate photon?
   * @details if there is a candidate with a lower KE than threshold, what is it?
   * If there is no such candidate the return is std::numeric_limits<double>::lowest()
   * @param the interaction of interest (data or MC)
   * @return the kinetic energy of the photon/electron
  **/
  template <class T>
    double subleading_gOre_ke(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return interaction.subleading_gore_ke;
    }

  /**
   * @brief axial spread of the candidate photon
   * @param the interaction of interest (data or MC)
   **/
  template <class T>
    double gOre_axial_spread(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return static_cast<double>(interaction.photon_or_electron->axial_spread);
    }

  /**
   * @brief directional spread of the candidate photon
   * @param the interaction of interest (data or MC)
   **/
  template <class T>
    double gOre_directional_spread(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return static_cast<double>(interaction.photon_or_electron->directional_spread);
    }

  /**
   * @brief length of the candidate photon
   * @param the interaction of interest (data or MC)
   **/
  template <class T>
    double gOre_length(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return static_cast<double>(interaction.photon_or_electron->length);
    }

  /**
   * @brief the momentum magnitude of the candicate photon 
   **/
  template <class T>
    double gOre_momentum(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return static_cast<double>(interaction.photon_or_electron->p);
    }

  /**
   * @brief the transerverse momentum of the candidate photon
   **/
  template <class T>
    double gOre_dpT(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return pvars::dpT(*interaction.photon_or_electron);
    }

  /**
   * @brief the start dedx of the candidate photon
   **/
  template <class T>
    double gOre_start_dedx(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return interaction.photon_or_electron->start_dedx;
    }

  /**
   * @brief the straightness of the candidate photon
   **/
  template <class T>
    double gOre_straightness(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return interaction.photon_or_electron->start_straightness;
    }

  /**
   * @brief the vertex gap of the candidate photon
   **/
  template <class T>
    double gOre_gap(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return interaction.photon_or_electron->vertex_distance;
    }

  /**
   * @brief the photon/electron pid of the candidate photon
   * @details 1 is photon-like, -1 is electron-like
   **/
  template <class T>
    double gOre_score(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      std::vector<double> softmax_vec = pvars::particle_softmax_vec(*interaction.photon_or_electron);
      return softmax_vec[pvars::kPhoton] - softmax_vec[pvars::kElectron];
    }

  /**
   * @brief the photon/electron polar angle w.r.t the beam
   **/
  template <class T>
    double gOre_polar_angle(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return pvars::polar_angle(*interaction.photon_or_electron);
    }

  /**
   * @brief the photon/electron azimuthal angle w.r.t the beam
   **/
  template <class T>
    double gOre_azimuthal_angle(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return pvars::azimuthal_angle(*interaction.photon_or_electron);
    }

  /**
   * @brief distance between the vertex and the x-/y-side active volume boundaries
   * @details getting geom info from https://sbn-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=21693&filename=ICARUS_geometry_update_26Apr21_v3.pdf&version=3
   **/
  template <class T>
    double xy_wall_dist(const T& obj)
    {
      double min_dist = std::numeric_limits<double>::quiet_NaN();
      double vtx_x = vars::vertex_x(obj);
      double vtx_y = vars::vertex_y(obj);
      // taxi cab metric, so just min dist in each coordinate
      // vtx should be between, but if any are negative that's worth noting
      double pos_x_dist = (GORE_WALL_X_POS - vtx_x);
      double neg_x_dist = (vtx_x - GORE_WALL_X_NEG);
      double pos_y_dist = (GORE_WALL_Y_POS - vtx_y);
      double neg_y_dist = (vtx_y - GORE_WALL_Y_NEG);

      min_dist = std::min({pos_x_dist, neg_x_dist, pos_y_dist, neg_y_dist});

      return min_dist;
    }

  /**
   * @brief distance between the vertex and the z-side active volume boundaries
   * @details getting geom info from https://sbn-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=21693&filename=ICARUS_geometry_update_26Apr21_v3.pdf&version=3
   **/
  template <class T>
    double z_wall_dist(const T& obj)
    {
      double min_dist = std::numeric_limits<double>::quiet_NaN();
      double vtx_z = vars::vertex_z(obj);
      // taxi cab metric, so just min dist in each coordinate
      // vtx should be between, but if any are negative that's worth noting
      double pos_z_dist = (GORE_WALL_Z_POS - vtx_z);
      double neg_z_dist = (vtx_z - GORE_WALL_Z_NEG);

      min_dist = std::min({pos_z_dist, neg_z_dist});

      return min_dist;
    }

  /**
   * @brief sum of gOre shower KE (including sub-thresholds) 
   **/
  template <class T>
    double total_gOre_KE(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return interaction.total_gore_ke;
    }

  /**
   * @brief how many subthreshold gOre showers are there?
   **/
  template <class T>
    double n_gOre_showers(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      return interaction.nShowers;
    }

  /**
   * @brief if there is a subthreshold shower, approximate the pion_mass peak
   * @details m = sqrt(2*E_1*E_2*(1 - CosTh))
   **/
  template <class T>
    double pion_mass(const T& obj)
    {
      core::nc::gOre::Interaction<T> interaction(obj);
      if (not interaction.is_valid)
        return std::numeric_limits<double>::quiet_NaN();
      // if there is only one shower reconstructed, assume it is actually 2 with (1 - CosTh) ~ directional spread
      // ie sqrt(2*spread)*E, The directional spread may need a scale factor, but this is a first pass
      if (interaction.nShowers == 1)
      {
        // truth particles don't have directional spread, so return a quiet NaN (for now at least)
        if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
        {
          return std::numeric_limits<double>::quiet_NaN();
        }
        else
        {
          return std::sqrt(2.*interaction.photon_or_electron->directional_spread)*interaction.photon_or_electron->ke;
        }
      }
      // in this case there must be a subleading shower. Use it for the pion mass estimate.
      return std::sqrt(2.*interaction.photon_or_electron->ke*interaction.subleading_gore_ke*(1.-interaction.pion_costh));
    }

  //*** TRUTH ONLY VARS ***//
  /**
   * @brief get the resonance number from the MC Truth
   * @details from genie::EResonance enum (I think)
   *  kNoResonance = -1,
   *  kP33_1232    =  0,
   *  kS11_1535    =  1,
   *  kD13_1520    =  2,
   *  kS11_1650    =  3,
   *  kD13_1700    =  4,
   *  kD15_1675    =  5,
   *  kS31_1620    =  6,
   *  kD33_1700    =  7,
   *  kP11_1440    =  8,
   *  kP33_1600    =  9,
   *  kP13_1720    = 10,
   *  kF15_1680    = 11,
   *  kP31_1910    = 12,
   *  kP33_1920    = 13,
   *  kF35_1905    = 14,
   *  kF37_1950    = 15,
   *  kP11_1710    = 16,
   *  kF17_1970    = 17
   * @param obj the interaction of interest (MCTruth only)
   * @return the resonance number (straight from GENIE)
   **/
  double baryon_res_code(const caf::Proxy<caf::SRTrueInteraction>& obj)
  {
    return obj.resnum;
  }

  //*** CATEGORIES ***//
  /**
   * @brief categorize the signals
   * @details 0: signal 1gXp
   * 1: background 1eXp
   * 2: background 2gXp
   * 3: background non-fiducial/contained 1gXp
   * 4: background Other NC
   * 5: background Other CC
   * 6: background non-neutrino
   **/
  double gOre_category_DEPRECIATED(const caf::SRInteractionTruthDLPProxy& obj)
  {
    double cat(6);
    std::vector<size_t> topology = core::nc::gOre::count_primaries(obj);
    bool no_muon_no_pion = (topology.at(pvars::kMuon) + topology.at(pvars::kPion) == 0);
    if      (cuts::neutrino(obj)                 &&
             cuts::fiducial_containment_cut(obj) &&
             cuts::nc::gOre::gOre_is_photon(obj)  )
      cat = 0;
    else if (cuts::neutrino(obj)                 &&
             cuts::fiducial_containment_cut(obj) &&
             cuts::nc::gOre::gOre_is_electron(obj))
      cat = 1;
    else if (cuts::neutrino(obj)                 &&
             cuts::fiducial_containment_cut(obj) &&
             (topology.at(pvars::kPhoton) == 2)  &&
             no_muon_no_pion                      )
      cat = 2;
    else if (cuts::neutrino(obj)                     &&
             not cuts::fiducial_containment_cut(obj) &&
             cuts::nc::gOre::gOre_is_photon(obj)      )
      cat = 3;
    else if (cuts::neutrino(obj)      &&
             cuts::nc::gOre::is_nc(obj))
      cat = 4;
    else if (cuts::neutrino(obj)      &&
             cuts::nc::gOre::is_cc(obj))
      cat = 5;
    return cat;
  }

  /**
   * @brief categorize by interaction type/mode
   * @details All categories before 6 are fiducial/contained
   * 0: signal 1gXp
   * 1: background 1eXp
   * 2: Res Pi0
   * 3: Res Pi+/Pi-
   * 4: QE
   * 5: DIS
   * 6: Non-fiducial/uncontained Neutrino
   * 7: Other Neutrino
   * 8: Non-Neutrino
   **/
  double gOre_category(const caf::SRInteractionTruthDLPProxy& obj)
  {
    double cat(8);
    int64_t type = obj.interaction_type;
    int64_t mode = obj.interaction_mode;
    bool is_ResPi0 = (type == caf::kResCCNuNeutronPi0)   ||
                     (type == caf::kResNCNuProtonPi0)    ||
                     (type == caf::kResNCNuNeutronPi0)   ||
                     (type == caf::kResCCNuBarProtonPi0) ||
                     (type == caf::kResNCNuBarProtonPi0) ||
                     (type == caf::kResNCNuBarNeutronPi0) ;
    bool is_ResPiPlus = (type == caf::kResCCNuProtonPiPlus)      ||
                        (type == caf::kResCCNuNeutronPiPlus)     ||
                        (type == caf::kResNCNuProtonPiPlus)      ||
                        (type == caf::kResNCNuBarProtonPiPlus)   ||
                        (type == caf::kResCCNuDeltaPlusPiPlus)   ||
                        (type == caf::kResCCNuBarDeltaMinusPiPlus);
    bool is_ResPiMinus = (type == caf::kResNCNuNeutronPiMinus)    ||
                         (type == caf::kResCCNuBarNeutronPiMinus) ||
                         (type == caf::kResCCNuBarProtonPiMinus)  ||
                         (type == caf::kResNCNuBarNeutronPiMinus) ||
                         (type == caf::kResCCNuDelta2PlusPiMinus) ||
                         (type == caf::kResCCNuBarDelta0PiMinus)   ;
    if      (cuts::neutrino(obj)                   &&
             cuts::fiducial_containment_cut(obj)   &&
             cuts::nc::gOre::gOre_is_photon(obj)    )
      cat = 0;
    else if (cuts::neutrino(obj)                   &&
             cuts::fiducial_containment_cut(obj)   &&
             cuts::nc::gOre::gOre_is_electron(obj  ))
      cat = 1;
    else if (cuts::neutrino(obj)                   &&
             cuts::fiducial_containment_cut(obj)   &&
             is_ResPi0                              )
      cat = 2;
    else if (cuts::neutrino(obj)                   &&
             cuts::fiducial_containment_cut(obj)   &&
             (is_ResPiPlus || is_ResPiMinus)        )
      cat = 3;
    else if (cuts::neutrino(obj)                   &&
             cuts::fiducial_containment_cut(obj)   &&
             (mode == caf::kQE)                     )
      cat = 4;
    else if (cuts::neutrino(obj)                   &&
             cuts::fiducial_containment_cut(obj)   &&
             (mode == caf::kDIS)                    )
      cat = 5;
    else if (cuts::neutrino(obj)                   &&
             not cuts::fiducial_containment_cut(obj))
      cat = 6;
    else if (cuts::neutrino(obj))
      cat = 7;
    return cat;
  }

  /**
   * @brief categorize for MC Truth
   * @details The fiducialization is being a pain, so focus on FSI topology and mode/interaction tyoe
   * 0: NC ∆->Nγ
   * 1: NC Other Single Photon
   * 2: NC π0 No Photon
   * 3: Other NC
   * 4: CC
   * 5: Non-neutrino
   **/
  double gOre_mc_category(const caf::Proxy<caf::SRTrueInteraction>& obj)
  {
    double cat(5);
    bool isnc = obj.isnc;
    caf::genie_interaction_mode_ genie_mode = obj.genie_mode;
    caf::genie_interaction_type_ genie_inttype = obj.genie_inttype;
    int resnum = obj.resnum;
    // post-FSI primary particles
    core::nc::gOre::mc_topology topology(obj.prim);
    // is NC ∆ res
    bool is_nc_delta_res = isnc && (resnum == 0) && (genie_inttype == 1000);
    // single photon topology (1γ and maybe some nucleons)
    bool is_single_photon_topology = topology.single_photon() && topology.only_photons_and_nucleons();
    // has π0 no photons
    bool has_pi0_0g = topology.has_pi0() && not topology.has(22);
    // is a neutrino interaction
    bool is_nu = obj.index != -1;
    if (is_nc_delta_res && is_single_photon_topology)
      cat = 0;
    else if (isnc && is_single_photon_topology)
      cat = 1;
    else if (isnc && has_pi0_0g)
      cat = 2;
    else if (isnc)
      cat = 3;
    else if (is_nu)
      cat = 4;
    return cat;
  }
} //end vars::nc::gOre namespace

#endif

/**
 * @file utils_gOre.h
 * @brief utilities useful for the gOre analysis
 * @details define analysis functions, constants, lists, etc here
 **/

#ifndef UTIL_GORE_H
#define UTIL_GORE_H

// preprocessor consts
#define PLACEHOLDERVALUE std::numeric_limits<double>::quiet_NaN()
#define PROTON_BINDING_ENERGY 30.9 // MeV
#define BEAM_IS_NUMI false
#define PIDFUNC pvars::pid

// SPINE includes
#include "include/mctruth.h"
#include "include/variables.h"
#include "include/muon2024/variables_muon2024.h"
#include "include/cuts.h"
#include "include/muon2024/cuts_muon2024.h"
#include "include/spinevar.h"
#include "include/analysis.h"
#include "include/nc/vars_gOre.h"
#include "include/nc/cuts_gOre.h"
#include "include/yell_try_die.h"

// sbncode includes
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "TDirectory.h"
#include "TFile.h"

namespace gOreVar = vars::nc::gOre;
namespace gOreCut = cuts::nc::gOre;

namespace gOre
{
  /**
   * Make maps of the vars we want to plot, keyed by the name we want to use in the TTrees
   * Define some aliases for the function ptrs b/c I think it's an ugly signature
   **/
  using  TVAR = double (*)(const TTYPE&);
  using  RVAR = double (*)(const RTYPE&);
  using MCVAR = double (*)(const MCTRUTH&);
  std::unordered_map<std::string, TVAR> ttypeVars =
  {
    {"neutrino_id",              &vars::neutrino_id                       },
    {"true_edep",                &vars::visible_energy                    },
    {"true_edep_calosub",        &vars::visible_energy_calosub            },
    {"is_true_gOre",             WRAP_BOOL(gOreCut::gOre_topology)        },
    {"has_true_subthresh_gOre",  WRAP_BOOL(gOreCut::has_subthreshold_gOre)},
    {"gOre_is_photon",           WRAP_BOOL(gOreCut::gOre_is_photon)       },
    {"gOre_is_electron",         WRAP_BOOL(gOreCut::gOre_is_electron)     },
    {"category_gOre",            &gOreVar::gOre_category                  },
    {"n_true_protons",           &gOreVar::n_protons                      },
    {"n_true_subthresh_protons", &gOreVar::n_protons_subthreshold         },
    {"n_true_gOre_showers",      &gOreVar::n_gOre_showers                 },
    {"total_true_proton_ke",     &gOreVar::total_proton_ke                },
    {"gOre_true_ke",             &gOreVar::gOre_ke                        },
    {"true_total_gOre_ke",       &gOreVar::total_gOre_KE                  },
    {"true_pion_mass",           &gOreVar::pion_mass                      },
    {"true_min_muon_ke",         &gOreVar::min_muon_ke                    },
    {"true_min_pion_ke",         &gOreVar::min_pion_ke                    },
    {"true_vertex_x",            &vars::vertex_x                          },
    {"true_vertex_y",            &vars::vertex_y                          },
    {"true_vertex_z",            &vars::vertex_z                          },
    {"true_wall_xy",             &gOreVar::xy_wall_dist                   },
    {"true_wall_z",              &gOreVar::z_wall_dist                    }
  };
  std::unordered_map<std::string, RVAR> rtypeVars =
  {
    {"reco_edep",                &vars::visible_energy             },
    {"reco_edep_calosub",        &vars::visible_energy_calosub     },
    {"flash_time",               &vars::flash_time                 },
    {"flash_total_PE",           &vars::flash_total_pe             },
    {"flash_hypothesis",         &vars::flash_hypothesis           },
    {"vertex_x",                 &vars::vertex_x                   },
    {"vertex_y",                 &vars::vertex_y                   },
    {"vertex_z",                 &vars::vertex_z                   },
    {"wall_xy",                  &gOreVar::xy_wall_dist            },
    {"wall_z",                   &gOreVar::z_wall_dist             },
    {"fiducial_cut",             WRAP_BOOL(cuts::fiducial_cut)     },
    {"containment_cut",          WRAP_BOOL(cuts::containment_cut)  },
    {"flash_cut",                WRAP_BOOL(cuts::flash_cut)        },
    {"is_gOre",                  WRAP_BOOL(gOreCut::gOre_topology) },
    {"n_protons",                &gOreVar::n_protons               },
    {"n_subthresh_protons",      &gOreVar::n_protons_subthreshold  },
    {"n_gOre_showers",           &gOreVar::n_gOre_showers          },
    {"total_proton_ke",          &gOreVar::total_proton_ke         },
    {"gOre_ke",                  &gOreVar::gOre_ke                 },
    {"total_gOre_ke",            &gOreVar::total_gOre_KE           },
    {"subleading_gOre_ke",       &gOreVar::subleading_gOre_ke      },
    {"pion_mass",                &gOreVar::pion_mass               },
    {"gOre_axial_spread",        &gOreVar::gOre_axial_spread       },
    {"gOre_directional_spread",  &gOreVar::gOre_directional_spread },
    {"gOre_azimuthal_angle",     &gOreVar::gOre_azimuthal_angle    },
    {"gOre_polar_angle",         &gOreVar::gOre_polar_angle        },
    {"gOre_length",              &gOreVar::gOre_length             },
    {"gOre_momentum",            &gOreVar::gOre_momentum           },
    {"gOre_dpT",                 &gOreVar::gOre_dpT                },
    {"gOre_start_dedx",          &gOreVar::gOre_start_dedx         },
    {"gOre_straightness",        &gOreVar::gOre_straightness       },
    {"gOre_gap",                 &gOreVar::gOre_gap                },
    {"gOre_score",               &gOreVar::gOre_score              },
    {"min_muon_ke",              &gOreVar::min_muon_ke             },
    {"min_pion_ke",              &gOreVar::min_pion_ke             }
  };
  std::unordered_map<std::string, MCVAR> mctruthVars =
  {
    {"baseline",                 &mctruth::true_neutrino_baseline        },
    {"pdg_code",                 &mctruth::true_neutrino_pdg             },
    {"is_cc",                    &mctruth::true_neutrino_cc              },
    {"interaction_mode",         &mctruth::interaction_mode              },
    {"interaction_type",         &mctruth::interaction_type              },
    {"true_energy",              &mctruth::true_neutrino_energy          },
    {"res_code",                 &gOreVar::baryon_res_code               },
    {"nc_delta_res",             WRAP_BOOL(gOreCut::nc_delta_res)        },
    {"nc_delta_res_no_pion",     WRAP_BOOL(gOreCut::nc_delta_res_no_pion)},
    {"nc_delta_res_pi0",         WRAP_BOOL(gOreCut::nc_delta_res_pi0)    }
  };
  
  /**
   * @brief Fill a SpillMultiVarMap based on the provided cuts and variables
   * @details Each of the variable types, (TVAR, RVAR, & MCVAR) take slightly different signatures to make the SpineVar,
   * so those maps are handles separeately.
   * @tparam CUT_TYPE What sort of type is being applied (typically TTYPE or RTYPE)
   * @param varMap The map to be filled
   * @param cut What is the selection cut?
   * @param cat What true category are you applying your selection to?
   **/
  template <class CUT_TYPE>
    void fillMultiVarMap(std::map<std::string, ana::SpillMultiVar>& varMap, bool (*cut)(const CUT_TYPE&), bool (*cat)(const TTYPE&))
    {
      for (const auto varPair : ttypeVars)
        varMap.insert({varPair.first, SpineVar<  TTYPE, CUT_TYPE>(varPair.second, cut, cat)});
      for (const auto varPair : rtypeVars)
        varMap.insert({varPair.first, SpineVar<  RTYPE, CUT_TYPE>(varPair.second, cut, cat)});
      for (const auto varPair : mctruthVars)
        varMap.insert({varPair.first, SpineVar<MCTRUTH, CUT_TYPE>(varPair.second, cut, cat)});
    }
  } // end gOre namespace

#endif
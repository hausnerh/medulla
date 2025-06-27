/**
 * @file nc_1gamma.C
 * @brief A basic analysis macro with a loose NC selection
 * @details This macro takes liberal use of Justin Mueller's example analysis macro
 * @author hhausner@fnal.gov
 **/

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
#include "include/nc/cuts_nc.h"

// sbncode includes
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "TDirectory.h"
#include "TFile.h"

/**
 * Which SPINE files to load, updates occationally, so double check
 * The extended CV sample is much larger, so test on the smaller 'nominal' sample
 **/
ana::SpectrumLoader var00("/pnfs/icarus/persistent/users/mueller/mixed/simulation/cvext/input*.flat.root");
//ana::SpectrumLoader var00("/pnfs/icarus/persistent/users/mueller/mixed/simulation/nominal/input*.flat.root");

/**
 * Make maps of the vars we want to plot, keyed by the name we want to use in the TTrees
 * Define some aliases for the function ptrs b/c I think it's an ugly signature
 **/
using  TVAR = double (*)(const TTYPE&);
using  RVAR = double (*)(const RTYPE&);
using MCVAR = double (*)(const MCTRUTH&);
std::unordered_map<std::string, TVAR> ttypeVars =
{
  {"neutrino_id",              &vars::neutrino_id                               },
  {"true_edep",                &vars::visible_energy                            },
  {"true_photon_E",            &vars::nc::photon_energy                         },
  {"true_proton_E",            &vars::nc::proton_energy                         },
  {"true_photon_E_frac",       &vars::nc::photon_energy_frac                    },
  {"true_proton_E_frac",       &vars::nc::proton_energy_frac                    },
  {"category_1g",              &vars::nc::category_single_photon                },
  {"category_1g1p",            &vars::nc::category_1photon_1proton              },
  {"alt_category_1g",          &vars::nc::alt_category_single_photon            },
  {"alt_category_1g1p",        &vars::nc::alt_category_1photon_1proton          },
  {"category_1gOr1e",          &vars::nc::category_single_photon_or_electron    },
  {"true_quality_photon",      WRAP_BOOL(cuts::nc::photon_quality)              },
  {"true_photons_abvThrsh",    &vars::nc::count_photons_above_threshold_double  },
  {"true_protons_abvThrsh",    &vars::nc::count_protons_above_threshold_double  },
  {"true_pions_abvThrsh",      &vars::nc::count_pions_above_threshold_double    },
  {"true_muons_abvThrsh",      &vars::nc::count_muons_above_threshold_double    },
  {"true_electrons_abvThrsh",  &vars::nc::count_electrons_above_threshold_double},
  {"true_photons_noThrsh",     &vars::nc::count_photons_no_threshold_double     },
  {"true_protons_noThrsh",     &vars::nc::count_protons_no_threshold_double     },
  {"true_pions_noThrsh",       &vars::nc::count_pions_no_threshold_double       },
  {"true_muons_noThrsh",       &vars::nc::count_muons_no_threshold_double       },
  {"true_electrons_noThrsh",   &vars::nc::count_electrons_no_threshold_double   },
  {"true_photon_dpT",          &vars::nc::photon_dpT                            },
  {"true_proton_dpT",          &vars::nc::proton_dpT                            },
  {"true_photon_dpT_frac",     &vars::nc::photon_dpT_frac                       },
  {"true_proton_dpT_frac",     &vars::nc::proton_dpT_frac                       },
  {"true_proton_photon_gap",   &vars::nc::photon_proton_gap                     },
  {"true_proton_photon_cosTh", &vars::nc::photon_proton_cosTh                   },
  {"true_photon_polar",        &vars::nc::photon_polar_angle                    },
  {"true_photon_azimuthal",    &vars::nc::photon_azimuthal_angle                },
  {"true_electron_vtx_dist",   &vars::nc::electron_vtx_dist                     },
  {"true_photon_vtx_dist",     &vars::nc::photon_vtx_dist                       }
};
std::unordered_map<std::string, RVAR> rtypeVars =
{
  {"reco_edep",                &vars::visible_energy                                                        },
  {"reco_photon_E",            &vars::nc::photon_energy                                                     },
  {"reco_proton_E",            &vars::nc::proton_energy                                                     },
  {"reco_photon_E_frac",       &vars::nc::photon_energy_frac                                                },
  {"reco_proton_E_frac",       &vars::nc::proton_energy_frac                                                },
  {"flash_time",               &vars::flash_time                                                            },
  {"flash_total_PE",           &vars::flash_total_pe                                                        },
  {"flash_hypothesis",         &vars::flash_hypothesis                                                      },
  {"fiducial_cut",             WRAP_BOOL(cuts::fiducial_cut)                                                },
  {"containment_cut",          WRAP_BOOL(cuts::containment_cut)                                             },
  {"flash_cut",                WRAP_BOOL(cuts::flash_cut)                                                   },
  {"vtx_cut",                  WRAP_BOOL(cuts::nc::photon_not_at_vertex)                                    },
  {"selected_1g",              WRAP_BOOL(cuts::nc::fiducial_containment_flash_cut_single_photon)            },
  {"selected_1gOr1e",          WRAP_BOOL(cuts::nc::fiducial_containment_flash_cut_single_photon_or_electron)},
  {"selected_1g1p",            WRAP_BOOL(cuts::nc::fiducial_containment_flash_cut_1photon_1proton)          },
  {"quality_photon",           WRAP_BOOL(cuts::nc::photon_quality)                                          },
  {"reco_photons_abvThrsh",    &vars::nc::count_photons_above_threshold_double                              },
  {"reco_protons_abvThrsh",    &vars::nc::count_protons_above_threshold_double                              },
  {"reco_pions_abvThrsh",      &vars::nc::count_pions_above_threshold_double                                },
  {"reco_muons_abvThrsh",      &vars::nc::count_muons_above_threshold_double                                },
  {"reco_electrons_abvThrsh",  &vars::nc::count_electrons_above_threshold_double                            },
  {"reco_photon_dpT",          &vars::nc::photon_dpT                                                        },
  {"reco_proton_dpT",          &vars::nc::proton_dpT                                                        },
  {"reco_photon_dpT_frac",     &vars::nc::photon_dpT_frac                                                   },
  {"reco_proton_dpT_frac",     &vars::nc::proton_dpT_frac                                                   },
  {"reco_proton_photon_gap",   &vars::nc::photon_proton_gap                                                 },
  {"reco_proton_photon_cosTh", &vars::nc::photon_proton_cosTh                                               },
  {"reco_photon_softmax",      &vars::nc::primary_photon_photon_score                                       },
  {"reco_electron_softmax",    &vars::nc::primary_photon_electron_score                                     },
  {"reco_photon_secondary_id", &vars::nc::photon_secondary_classification                                   },
  {"reco_photon_polar",        &vars::nc::photon_polar_angle                                                },
  {"reco_photon_azimuthal",    &vars::nc::photon_azimuthal_angle                                            },
  {"electron_electron_id",     &vars::nc::electron_electron_id                                              },
  {"electron_photon_id",       &vars::nc::electron_photon_id                                                },
  {"electron_photon_skew",     &vars::nc::electron_photon_skew                                              },
  {"photon_alignment",         &vars::nc::photon_alignment                                                  },
  {"electron_vtx_dist",        &vars::nc::electron_vtx_dist                                                 },
  {"electron_start_dedx",      &vars::nc::electron_start_dedx                                               },
  {"electron_straightness",    &vars::nc::electron_straightness                                             },
  {"photon_vtx_dist",          &vars::nc::photon_vtx_dist                                                   },
  {"photon_start_dedx",        &vars::nc::photon_start_dedx                                                 },
  {"photon_straightness",      &vars::nc::photon_straightness                                               },
  {"photon_directional_spread",&vars::nc::photon_directional_spread                                         },
  {"photon_axial_spread",      &vars::nc::photon_axial_spread                                               }
};
std::unordered_map<std::string, MCVAR> mctruthVars =
{
  {"baseline",                 &mctruth::true_neutrino_baseline              },
  {"pdg_code",                 &mctruth::true_neutrino_pdg                   },
  {"is_cc",                    &mctruth::true_neutrino_cc                    },
  {"interaction_mode",         &mctruth::interaction_mode                    },
  {"interaction_type",         &mctruth::interaction_type                    },
  {"true_energy",              &mctruth::true_neutrino_energy                },
  {"res_code",                 &vars::nc::baryon_res_code                    },
  {"mc_nprim",                 &vars::nc::mc_nprim                           },
  {"nc_delta_res",             WRAP_BOOL(cuts::nc::true_nc_delta_res)        },
  {"nc_delta_res_no_pion",     WRAP_BOOL(cuts::nc::true_nc_delta_res_no_pion)}
};

/**
 * @brief Fill a SpillMultiVarMap based on the provided cuts and variables
 * @details Each of the variable types, (TVAR, RVAR, & MCVAR) take slightly different signatures to make the SpineVar,
 * so those maps are handles separeately.
 * @@tparam CUT_TYPE What sort of type is being applied (typically TTYPE or RTYPE)
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

void nc_1gamma()
{
    /**
     * @brief Create an instance of the Analysis class for the 1gamma analysis.
     * @details This creates an instance of the Analysis class, which is used
     * to run the analysis on the specified samples. The name of the analysis,
     * and therefore the name of the output file, is specified as an argument
     * to the constructor.
     **/
    ana::Analysis analysis("nc_1gamma");

    /**
     * @brief Add a sample to the analysis.
     * @details This adds a sample to the analysis by creating a SpectrumLoader
     * object and adding it to the Analysis class. The SpectrumLoader object
     * represents the sample in the analysis, and is used to load the data from
     * the ROOT file and apply the cuts and variables. The name passed to the
     * AddLoader function is used to create a directory in the output ROOT file
     * to store the results of the analysis.
     **/
    //ana::SpectrumLoader var00("/pnfs/icarus/persistent/users/mueller/fall2024/nominal/flat/*.root");
    analysis.AddLoader("nominal", &var00, true);

    /**
     * @brief Add a set of variables for selected interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     **/
    #define CUT cuts::nc::topology_single_photon
    #define TCUT cuts::neutrino
    std::map<std::string, ana::SpillMultiVar> vars_sel_nu_1g;
    fillMultiVarMap<RTYPE>(vars_sel_nu_1g, &CUT, &TCUT);
    analysis.AddTree("Nu_Topology_1photon", vars_sel_nu_1g, true);

    #undef TCUT
    #define TCUT cuts::cosmic
    std::map<std::string, ana::SpillMultiVar> vars_sel_cos_1g;
    fillMultiVarMap<RTYPE>(vars_sel_cos_1g, &CUT, &TCUT);
    analysis.AddTree("Cos_Topology_1photon", vars_sel_cos_1g, true);
   
    #undef CUT
    #define CUT cuts::nc::topology_single_photon_or_electron
    #undef TCUT
    #define TCUT cuts::neutrino
    std::map<std::string, ana::SpillMultiVar> vars_sel_nu_1gOr1e;
    fillMultiVarMap<RTYPE>(vars_sel_nu_1gOr1e, &CUT, &TCUT);
    analysis.AddTree("Nu_Topology_1photonOr1electron", vars_sel_nu_1gOr1e, true);

    #undef TCUT
    #define TCUT cuts::cosmic
    std::map<std::string, ana::SpillMultiVar> vars_sel_cos_1gOr1e;
    fillMultiVarMap<RTYPE>(vars_sel_cos_1gOr1e, &CUT, &TCUT);
    analysis.AddTree("Cos_Topology_1photonOr1electron", vars_sel_cos_1gOr1e, true);

    #undef CUT
    #define CUT cuts::nc::topology_1photon_1proton
    #undef TCUT
    #define TCUT cuts::neutrino
    std::map<std::string, ana::SpillMultiVar> vars_sel_nu_1g1p;
    fillMultiVarMap<RTYPE>(vars_sel_nu_1g1p, &CUT, &TCUT);
    analysis.AddTree("Nu_Topology_1photon1proton", vars_sel_nu_1g1p, true);

    #undef TCUT
    #define TCUT cuts::cosmic
    std::map<std::string, ana::SpillMultiVar> vars_sel_cos_1g1p;
    fillMultiVarMap<RTYPE>(vars_sel_cos_1g1p, &CUT, &TCUT);
    analysis.AddTree("Cos_Topology_1photon1proton", vars_sel_cos_1g1p, true);
    
    /**
     * @brief Add a set of variables for signal interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     **/
    #define SIGCUT cuts::nc::is_fid_con_single_photon
    std::map<std::string, ana::SpillMultiVar> vars_sig_1g;
    fillMultiVarMap<TTYPE>(vars_sig_1g, &SIGCUT, &SIGCUT);
    analysis.AddTree("signalNu_Topology_1photon", vars_sig_1g, true);

    #undef SIGCUT
    #define SIGCUT cuts::nc::is_fid_con_1photon_1proton
    std::map<std::string, ana::SpillMultiVar> vars_sig_1g1p;
    fillMultiVarMap<TTYPE>(vars_sig_1g1p, &SIGCUT, &SIGCUT);
    analysis.AddTree("signalNu_Topology_1photon1proton", vars_sig_1g1p, true);

    // NB: must also apply (nc_delta_res_no_pion) to get the desired signal
    #undef SIGCUT
    #define SIGCUT cuts::nc::is_fid_con_nc_nu_res
    std::map<std::string, ana::SpillMultiVar> vars_fid_con_nc_res;
    fillMultiVarMap<TTYPE>(vars_fid_con_nc_res, &SIGCUT, &SIGCUT);
    analysis.AddTree("signalNu_FidCon_NCRes", vars_fid_con_nc_res, true);

    /**
     * @brief Run the analysis.
     * @details This runs the analysis on the samples specified by the
     * SpectrumLoaders and variables added to the Analysis class. It loops over
     * each sample (here only one), applies the cuts and variables to the data,
     * and stores the results in a TFile.
     **/
    analysis.Go();
}

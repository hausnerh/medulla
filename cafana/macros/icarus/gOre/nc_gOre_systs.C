/**
 * @file nc_gOre_train.C
 * @brief A basic analysis macro with a loose NC selection
 * Make different trees for training samples
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

/**
 * Make maps of the vars we want to plot, keyed by the name we want to use in the TTrees
 * Define some aliases for the function ptrs b/c I think it's an ugly signature
 **/
using  TVAR = double (*)(const TTYPE&);
using  RVAR = double (*)(const RTYPE&);
using MCVAR = double (*)(const MCTRUTH&);
std::unordered_map<std::string, TVAR> ttypeVars =
{
  {"neutrino_id",              &vars::neutrino_id                  },
  {"true_edep",                &vars::visible_energy               },
  {"true_edep_calosub",        &vars::visible_energy_calosub       },
  {"is_true_gOre",             WRAP_BOOL(gOreCut::gOre_topology)   },
  {"gOre_is_photon",           WRAP_BOOL(gOreCut::gOre_is_photon)  },
  {"gOre_is_electron",         WRAP_BOOL(gOreCut::gOre_is_electron)},
  {"category_gOre",            &gOreVar::gOre_category             },
  {"n_true_protons",           &gOreVar::n_protons                 },
  {"total_true_proton_ke",     &gOreVar::total_proton_ke           },
  {"gOre_true_ke",             &gOreVar::gOre_ke                   }
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
  {"fiducial_cut",             WRAP_BOOL(cuts::fiducial_cut)     },
  {"containment_cut",          WRAP_BOOL(cuts::containment_cut)  },
  {"flash_cut",                WRAP_BOOL(cuts::flash_cut)        },
  {"is_gOre",                  WRAP_BOOL(gOreCut::gOre_topology) },
  {"n_protons",                &gOreVar::n_protons               },
  {"total_proton_ke",          &gOreVar::total_proton_ke         },
  {"gOre_ke",                  &gOreVar::gOre_ke                 },
  {"subleading_gOre_ke",       &gOreVar::subleading_gOre_ke      },
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
  {"gOre_score",               &gOreVar::gOre_score              }
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

int nc_gOre_systs()
{
    /**
     * @brief Create an instance of the Analysis class for the 1gamma analysis.
     * @details This creates an instance of the Analysis class, which is used
     * to run the analysis on the specified samples. The name of the analysis,
     * and therefore the name of the output file, is specified as an argument
     * to the constructor.
     **/
    ana::Analysis analysis("nc_gOre_systs");

    /**
     * @brief Samples to the analysis.
     **/
    std::string inputs_var00 = "/pnfs/icarus/persistent/users/mueller/mixed/simulation/var00/input*.flat.root";
    ana::SpectrumLoader var_var00(inputs_var00);
    analysis.AddLoader("var00", &var_var00, true);

    std::string inputs_var01 = "/pnfs/icarus/persistent/users/mueller/mixed/simulation/var01/input*.flat.root";
    ana::SpectrumLoader var_var01(inputs_var01);
    analysis.AddLoader("var01", &var_var01, true);

    std::string inputs_var02 = "/pnfs/icarus/persistent/users/mueller/mixed/simulation/var02/input*.flat.root";
    ana::SpectrumLoader var_var02(inputs_var02);
    analysis.AddLoader("var02", &var_var02, true);

    std::string inputs_var03m = "/pnfs/icarus/persistent/users/mueller/mixed/simulation/var03m/input*.flat.root";
    ana::SpectrumLoader var_var03m(inputs_var03m);
    analysis.AddLoader("var03m", &var_var03m, true);

    std::string inputs_var03p = "/pnfs/icarus/persistent/users/mueller/mixed/simulation/var03p/input*.flat.root";
    ana::SpectrumLoader var_var03p(inputs_var03p);
    analysis.AddLoader("var03p", &var_var03p, true);

    std::string inputs_var04 = "/pnfs/icarus/persistent/users/mueller/mixed/simulation/var04/input*.flat.root";
    ana::SpectrumLoader var_var04(inputs_var04);
    analysis.AddLoader("var04", &var_var04, true);

    std::string inputs_var05 = "/pnfs/icarus/persistent/users/mueller/mixed/simulation/var05/input*.flat.root";
    ana::SpectrumLoader var_var05(inputs_var05);
    analysis.AddLoader("var05", &var_var05, true);

    std::string inputs_var06 = "/pnfs/icarus/persistent/users/mueller/mixed/simulation/var06/input*.flat.root";
    ana::SpectrumLoader var_var06(inputs_var06);
    analysis.AddLoader("var06", &var_var06, true);

    std::string inputs_var07 = "/pnfs/icarus/persistent/users/mueller/mixed/simulation/var07/input*.flat.root";
    ana::SpectrumLoader var_var07(inputs_var07);
    analysis.AddLoader("var07", &var_var07, true);

    std::string inputs_var08 = "/pnfs/icarus/persistent/users/mueller/mixed/simulation/var08/input*.flat.root";
    ana::SpectrumLoader var_var08(inputs_var08);
    analysis.AddLoader("var08", &var_var08, true);

    std::string inputs_var09m = "/pnfs/icarus/persistent/users/mueller/mixed/simulation/var09m/input*.flat.root";
    ana::SpectrumLoader var_var09m(inputs_var09m);
    analysis.AddLoader("var09m", &var_var09m, true);

    std::string inputs_var09p = "/pnfs/icarus/persistent/users/mueller/mixed/simulation/var09p/input*.flat.root";
    ana::SpectrumLoader var_var09p(inputs_var09p);
    analysis.AddLoader("var09p", &var_var09p, true);

    std::string inputs_var10m = "/pnfs/icarus/persistent/users/mueller/mixed/simulation/var10m/input*.flat.root";
    ana::SpectrumLoader var_var10m(inputs_var10m);
    analysis.AddLoader("var10m", &var_var10m, true);

    std::string inputs_var10p = "/pnfs/icarus/persistent/users/mueller/mixed/simulation/var10p/input*.flat.root";
    ana::SpectrumLoader var_var10p(inputs_var10p);
    analysis.AddLoader("var10p", &var_var10p, true);

    std::string inputs_nominal = "/pnfs/icarus/persistent/users/mueller/mixed/simulation/nominal/input*.flat.root";
    ana::SpectrumLoader var_nominal(inputs_nominal);
    analysis.AddLoader("nominal", &var_nominal, true);

    /**
     * @brief Add a set of variables for selected interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     **/
    #define CUT gOreCut::gOre_topology 
    #define TCUT cuts::neutrino
    std::map<std::string, ana::SpillMultiVar> vars_sel_nu_gOre;
    fillMultiVarMap<RTYPE>(vars_sel_nu_gOre, &CUT, &TCUT);
    analysis.AddTree("Nu_Topology_gOre", vars_sel_nu_gOre, true);

    #undef TCUT
    #define TCUT cuts::cosmic
    std::map<std::string, ana::SpillMultiVar> vars_sel_cos_gOre;
    fillMultiVarMap<RTYPE>(vars_sel_cos_gOre, &CUT, &TCUT);
    analysis.AddTree("Cos_Topology_gOre", vars_sel_cos_gOre, true);
    
    /**
     * @brief Add a set of variables for signal interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     **/
    #define SIGCUT gOreCut::is_fid_con_gOre
    std::map<std::string, ana::SpillMultiVar> vars_sig_gOre;
    fillMultiVarMap<TTYPE>(vars_sig_gOre, &SIGCUT, &SIGCUT);
    analysis.AddTree("signalNu_Topology_gOre", vars_sig_gOre, true);

    // NB: must also apply (nc_delta_res_no_pion) to get the desired signal
    #undef SIGCUT
    #define SIGCUT gOreCut::is_fid_con_nc_nu_res
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
    return 0;
}

int main(int argc, char* argv[])
{
  gErrorIgnoreLevel=3000;
  std::string sample(argv[1]);
  int ret = try_call("Make gOre Trees with Systematics", nc_gOre_systs);
  return ret;
}

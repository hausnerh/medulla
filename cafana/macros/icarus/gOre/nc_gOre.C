/**
 * @file nc_gOre.C
 * @brief A basic analysis macro with a loose NC selection
 * @details This macro takes liberal use of Justin Mueller's example analysis macro
 * @author hhausner@fnal.gov
 **/

#include "include/nc/utils_gOre.h"

int nc_gOre(const std::string& sample)
{
    /**
     * @brief Create an instance of the Analysis class for the 1gamma analysis.
     * @details This creates an instance of the Analysis class, which is used
     * to run the analysis on the specified samples. The name of the analysis,
     * and therefore the name of the output file, is specified as an argument
     * to the constructor.
     **/
    ana::Analysis analysis("nc_gOre_"+sample);

    /**
     * @brief Add a sample to the analysis.
     * @details This adds a sample to the analysis by creating a SpectrumLoader
     * object and adding it to the Analysis class. The SpectrumLoader object
     * represents the sample in the analysis, and is used to load the data from
     * the ROOT file and apply the cuts and variables. The name passed to the
     * AddLoader function is used to create a directory in the output ROOT file
     * to store the results of the analysis.
     **/
    std::string inputs = "/pnfs/icarus/persistent/users/mueller/mixed/simulation/" + sample + "/input*.flat.root";
    ana::SpectrumLoader var00(inputs);
    analysis.AddLoader("nominal", &var00, true);

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
    gOre::fillMultiVarMap<RTYPE>(vars_sel_nu_gOre, &CUT, &TCUT);
    analysis.AddTree("Nu_Topology_gOre", vars_sel_nu_gOre, true);

    #undef TCUT
    #define TCUT cuts::cosmic
    std::map<std::string, ana::SpillMultiVar> vars_sel_cos_gOre;
    gOre::fillMultiVarMap<RTYPE>(vars_sel_cos_gOre, &CUT, &TCUT);
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
    gOre::fillMultiVarMap<TTYPE>(vars_sig_gOre, &SIGCUT, &SIGCUT);
    analysis.AddTree("signalNu_Topology_gOre", vars_sig_gOre, true);

    // NB: must also apply (nc_delta_res_no_pion) to get the desired signal
    #undef SIGCUT
    #define SIGCUT gOreCut::is_fid_con_nc_nu_res
    std::map<std::string, ana::SpillMultiVar> vars_fid_con_nc_res;
    gOre::fillMultiVarMap<TTYPE>(vars_fid_con_nc_res, &SIGCUT, &SIGCUT);
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
  int ret = try_call("Make gOre Trees", nc_gOre, sample);
  return ret;
}

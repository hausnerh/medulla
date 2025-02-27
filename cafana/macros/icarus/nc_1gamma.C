/**
 * @file nc_selection.C
 * @brief A basic analysis macro with a loose NC selection
 * @details This macro takes liberal use of Justin Mueller's example analysis macro
 * @author hhausner@fnal.gov
*/

// preprocessor consts
#define PLACEHOLDERVALUE std::numeric_limits<double>::quiet_NaN()
#define PROTON_BINDING_ENERGY 30.9 // MeV
#define BEAM_IS_NUMI false
#define WRITE_PURITY_TREES false

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

void nc_1gamma()
{
    /**
     * @brief Create an instance of the Analysis class for the 1gamma analysis.
     * @details This creates an instance of the Analysis class, which is used
     * to run the analysis on the specified samples. The name of the analysis,
     * and therefore the name of the output file, is specified as an argument
     * to the constructor.
     */
    ana::Analysis analysis("nc_1gamma");

    /**
     * @brief Add a sample to the analysis.
     * @details This adds a sample to the analysis by creating a SpectrumLoader
     * object and adding it to the Analysis class. The SpectrumLoader object
     * represents the sample in the analysis, and is used to load the data from
     * the ROOT file and apply the cuts and variables. The name passed to the
     * AddLoader function is used to create a directory in the output ROOT file
     * to store the results of the analysis.
     */
    ana::SpectrumLoader var00("/pnfs/icarus/persistent/users/mueller/fall2024/nominal/flat/*.root");
    analysis.AddLoader("nominal", &var00, true);

    /**
     * @brief Add a set of variables for selected interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     */
    #define CUT cuts::nc::fiducial_containment_flash_cut_single_photon
    #define TCUT cuts::neutrino
    std::map<std::string, ana::SpillMultiVar> vars_sel_nu_1g;
    vars_sel_nu_1g.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_sel_nu_1g.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_sel_nu_1g.insert({"pdg", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_pdg, &CUT, &TCUT)});
    vars_sel_nu_1g.insert({"cc", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_cc, &CUT, &TCUT)});
    vars_sel_nu_1g.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_sel_nu_1g.insert({"interaction_type", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_type, &CUT, &TCUT)});
    vars_sel_nu_1g.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::muon2024::category, &CUT, &TCUT)});
    vars_sel_nu_1g.insert({"true_energy", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_energy, &CUT, &TCUT)});
    vars_sel_nu_1g.insert({"reco_edep", SpineVar<RTYPE,RTYPE>(&vars::visible_energy, &CUT, &TCUT)});
    vars_sel_nu_1g.insert({"flash_time", SpineVar<RTYPE,RTYPE>(&vars::flash_time, &CUT, &TCUT)});
    vars_sel_nu_1g.insert({"flash_total", SpineVar<RTYPE,RTYPE>(&vars::flash_total_pe, &CUT, &TCUT)});
    vars_sel_nu_1g.insert({"flash_hypothesis", SpineVar<RTYPE,RTYPE>(&vars::flash_hypothesis, &CUT, &TCUT)});
    vars_sel_nu_1g.insert({"reco_primary_photons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_photons_double, &CUT, &TCUT)});
    vars_sel_nu_1g.insert({"true_primary_photons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_photons_double, &CUT, &TCUT)});
    vars_sel_nu_1g.insert({"reco_primary_electrons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_electrons_double, &CUT, &TCUT)});
    vars_sel_nu_1g.insert({"true_primary_electrons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_electrons_double, &CUT, &TCUT)});
    vars_sel_nu_1g.insert({"reco_primary_muons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_muons_double, &CUT, &TCUT)});
    vars_sel_nu_1g.insert({"true_primary_muons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_muons_double, &CUT, &TCUT)});
    vars_sel_nu_1g.insert({"reco_primary_pions", SpineVar<RTYPE,RTYPE>(&vars::nc::count_pions_double, &CUT, &TCUT)});
    vars_sel_nu_1g.insert({"true_primary_pions", SpineVar<TTYPE,RTYPE>(&vars::nc::count_pions_double, &CUT, &TCUT)});
    vars_sel_nu_1g.insert({"reco_primary_protons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_protons_double, &CUT, &TCUT)});
    vars_sel_nu_1g.insert({"true_primary_protons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_protons_double, &CUT, &TCUT)});

    analysis.AddTree("selectedNu_1photon", vars_sel_nu_1g, true);

    #undef TCUT
    #define TCUT cuts::cosmic
    std::map<std::string, ana::SpillMultiVar> vars_sel_cos_1g;
    vars_sel_cos_1g.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"pdg", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_pdg, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"cc", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_cc, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"interaction_type", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_type, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::muon2024::category, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"true_energy", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_energy, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"reco_edep", SpineVar<RTYPE,RTYPE>(&vars::visible_energy, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"flash_time", SpineVar<RTYPE,RTYPE>(&vars::flash_time, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"flash_total", SpineVar<RTYPE,RTYPE>(&vars::flash_total_pe, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"flash_hypothesis", SpineVar<RTYPE,RTYPE>(&vars::flash_hypothesis, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"reco_primary_photons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_photons_double, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"true_primary_photons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_photons_double, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"reco_primary_electrons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_electrons_double, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"true_primary_electrons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_electrons_double, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"reco_primary_muons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_muons_double, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"true_primary_muons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_muons_double, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"reco_primary_pions", SpineVar<RTYPE,RTYPE>(&vars::nc::count_pions_double, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"true_primary_pions", SpineVar<TTYPE,RTYPE>(&vars::nc::count_pions_double, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"reco_primary_protons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_protons_double, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"true_primary_protons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_protons_double, &CUT, &TCUT)});
    
    analysis.AddTree("selectedCos_1photon", vars_sel_cos_1g, true);
    
    #undef CUT
    #define CUT cuts::nc::fiducial_containment_flash_cut_1photon_1proton
    #undef TCUT
    #define TCUT cuts::neutrino
    std::map<std::string, ana::SpillMultiVar> vars_sel_nu_1g1p;
    vars_sel_nu_1g1p.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_sel_nu_1g1p.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_sel_nu_1g1p.insert({"pdg", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_pdg, &CUT, &TCUT)});
    vars_sel_nu_1g1p.insert({"cc", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_cc, &CUT, &TCUT)});
    vars_sel_nu_1g1p.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_sel_nu_1g1p.insert({"interaction_type", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_type, &CUT, &TCUT)});
    vars_sel_nu_1g1p.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::muon2024::category, &CUT, &TCUT)});
    vars_sel_nu_1g1p.insert({"true_energy", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_energy, &CUT, &TCUT)});
    vars_sel_nu_1g1p.insert({"reco_edep", SpineVar<RTYPE,RTYPE>(&vars::visible_energy, &CUT, &TCUT)});
    vars_sel_nu_1g1p.insert({"flash_time", SpineVar<RTYPE,RTYPE>(&vars::flash_time, &CUT, &TCUT)});
    vars_sel_nu_1g1p.insert({"flash_total", SpineVar<RTYPE,RTYPE>(&vars::flash_total_pe, &CUT, &TCUT)});
    vars_sel_nu_1g1p.insert({"flash_hypothesis", SpineVar<RTYPE,RTYPE>(&vars::flash_hypothesis, &CUT, &TCUT)});
    vars_sel_nu_1g1p.insert({"reco_primary_photons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_photons_double, &CUT, &TCUT)});
    vars_sel_nu_1g1p.insert({"true_primary_photons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_photons_double, &CUT, &TCUT)});
    vars_sel_nu_1g1p.insert({"reco_primary_electrons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_electrons_double, &CUT, &TCUT)});
    vars_sel_nu_1g1p.insert({"true_primary_electrons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_electrons_double, &CUT, &TCUT)});
    vars_sel_nu_1g1p.insert({"reco_primary_muons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_muons_double, &CUT, &TCUT)});
    vars_sel_nu_1g1p.insert({"true_primary_muons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_muons_double, &CUT, &TCUT)});
    vars_sel_nu_1g1p.insert({"reco_primary_pions", SpineVar<RTYPE,RTYPE>(&vars::nc::count_pions_double, &CUT, &TCUT)});
    vars_sel_nu_1g1p.insert({"true_primary_pions", SpineVar<TTYPE,RTYPE>(&vars::nc::count_pions_double, &CUT, &TCUT)});
    vars_sel_nu_1g1p.insert({"reco_primary_protons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_protons_double, &CUT, &TCUT)});
    vars_sel_nu_1g1p.insert({"true_primary_protons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_protons_double, &CUT, &TCUT)});

    analysis.AddTree("selectedNu_1photon1proton", vars_sel_nu_1g1p, true);

    #undef TCUT
    #define TCUT cuts::cosmic
    std::map<std::string, ana::SpillMultiVar> vars_sel_cos_1g1p;
    vars_sel_cos_1g1p.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_sel_cos_1g1p.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_sel_cos_1g1p.insert({"pdg", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_pdg, &CUT, &TCUT)});
    vars_sel_cos_1g1p.insert({"cc", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_cc, &CUT, &TCUT)});
    vars_sel_cos_1g1p.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_sel_cos_1g1p.insert({"interaction_type", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_type, &CUT, &TCUT)});
    vars_sel_cos_1g1p.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::muon2024::category, &CUT, &TCUT)});
    vars_sel_cos_1g1p.insert({"true_energy", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_energy, &CUT, &TCUT)});
    vars_sel_cos_1g1p.insert({"reco_edep", SpineVar<RTYPE,RTYPE>(&vars::visible_energy, &CUT, &TCUT)});
    vars_sel_cos_1g1p.insert({"flash_time", SpineVar<RTYPE,RTYPE>(&vars::flash_time, &CUT, &TCUT)});
    vars_sel_cos_1g1p.insert({"flash_total", SpineVar<RTYPE,RTYPE>(&vars::flash_total_pe, &CUT, &TCUT)});
    vars_sel_cos_1g1p.insert({"flash_hypothesis", SpineVar<RTYPE,RTYPE>(&vars::flash_hypothesis, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"reco_primary_photons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_photons_double, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"true_primary_photons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_photons_double, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"reco_primary_electrons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_electrons_double, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"true_primary_electrons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_electrons_double, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"reco_primary_muons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_muons_double, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"true_primary_muons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_muons_double, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"reco_primary_pions", SpineVar<RTYPE,RTYPE>(&vars::nc::count_pions_double, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"true_primary_pions", SpineVar<TTYPE,RTYPE>(&vars::nc::count_pions_double, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"reco_primary_protons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_protons_double, &CUT, &TCUT)});
    vars_sel_cos_1g.insert({"true_primary_protons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_protons_double, &CUT, &TCUT)});
    
    analysis.AddTree("selectedCos_1photon1proton", vars_sel_cos_1g1p, true);
    
    /**
     * @brief Add a set of variables for signal interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     */
    #define SIGCUT cuts::nc::is_fid_con_nc_single_photon
    std::map<std::string, ana::SpillMultiVar> vars_sig_1g;
    vars_sig_1g.insert({"nu_id", SpineVar<TTYPE,TTYPE>(&vars::neutrino_id, &SIGCUT, &SIGCUT)});
    vars_sig_1g.insert({"baseline", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_baseline, &SIGCUT, &SIGCUT)});
    vars_sig_1g.insert({"pdg", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_pdg, &SIGCUT, &SIGCUT)});
    vars_sig_1g.insert({"cc", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_cc, &SIGCUT, &SIGCUT)});
    vars_sig_1g.insert({"interaction_mode", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_mode, &SIGCUT, &SIGCUT)});
    vars_sig_1g.insert({"interaction_type", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_type, &SIGCUT, &SIGCUT)});
    vars_sig_1g.insert({"category", SpineVar<TTYPE,TTYPE>(&vars::muon2024::category, &SIGCUT, &SIGCUT)});
    vars_sig_1g.insert({"fiducial_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::fiducial_cut), &SIGCUT, &SIGCUT)});
    vars_sig_1g.insert({"contain_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::containment_cut), &SIGCUT, &SIGCUT)});
    vars_sig_1g.insert({"flash_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::flash_cut), &SIGCUT, &SIGCUT)});
    vars_sig_1g.insert({"reco_primary_photons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_photons_double, &CUT, &TCUT)});
    vars_sig_1g.insert({"true_primary_photons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_photons_double, &CUT, &TCUT)});
    vars_sig_1g.insert({"reco_primary_electrons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_electrons_double, &CUT, &TCUT)});
    vars_sig_1g.insert({"true_primary_electrons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_electrons_double, &CUT, &TCUT)});
    vars_sig_1g.insert({"reco_primary_muons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_muons_double, &CUT, &TCUT)});
    vars_sig_1g.insert({"true_primary_muons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_muons_double, &CUT, &TCUT)});
    vars_sig_1g.insert({"reco_primary_pions", SpineVar<RTYPE,RTYPE>(&vars::nc::count_pions_double, &CUT, &TCUT)});
    vars_sig_1g.insert({"true_primary_pions", SpineVar<TTYPE,RTYPE>(&vars::nc::count_pions_double, &CUT, &TCUT)});
    vars_sig_1g.insert({"reco_primary_protons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_protons_double, &CUT, &TCUT)});
    vars_sig_1g.insert({"true_primary_protons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_protons_double, &CUT, &TCUT)});

    analysis.AddTree("signalNu_1photon", vars_sig_1g, true);

    #undef SIGCUT
    #define SIGCUT cuts::nc::is_fid_con_nc_1photon_1proton
    std::map<std::string, ana::SpillMultiVar> vars_sig_1g1p;
    vars_sig_1g1p.insert({"nu_id", SpineVar<TTYPE,TTYPE>(&vars::neutrino_id, &SIGCUT, &SIGCUT)});
    vars_sig_1g1p.insert({"baseline", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_baseline, &SIGCUT, &SIGCUT)});
    vars_sig_1g1p.insert({"pdg", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_pdg, &SIGCUT, &SIGCUT)});
    vars_sig_1g1p.insert({"cc", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_cc, &SIGCUT, &SIGCUT)});
    vars_sig_1g1p.insert({"interaction_mode", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_mode, &SIGCUT, &SIGCUT)});
    vars_sig_1g1p.insert({"interaction_type", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_type, &SIGCUT, &SIGCUT)});
    vars_sig_1g1p.insert({"category", SpineVar<TTYPE,TTYPE>(&vars::muon2024::category, &SIGCUT, &SIGCUT)});
    vars_sig_1g1p.insert({"fiducial_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::fiducial_cut), &SIGCUT, &SIGCUT)});
    vars_sig_1g1p.insert({"contain_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::containment_cut), &SIGCUT, &SIGCUT)});
    vars_sig_1g1p.insert({"flash_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::flash_cut), &SIGCUT, &SIGCUT)});
    vars_sig_1g1p.insert({"reco_primary_photons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_photons_double, &CUT, &TCUT)});
    vars_sig_1g1p.insert({"true_primary_photons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_photons_double, &CUT, &TCUT)});
    vars_sig_1g1p.insert({"reco_primary_electrons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_electrons_double, &CUT, &TCUT)});
    vars_sig_1g1p.insert({"true_primary_electrons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_electrons_double, &CUT, &TCUT)});
    vars_sig_1g1p.insert({"reco_primary_muons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_muons_double, &CUT, &TCUT)});
    vars_sig_1g1p.insert({"true_primary_muons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_muons_double, &CUT, &TCUT)});
    vars_sig_1g1p.insert({"reco_primary_pions", SpineVar<RTYPE,RTYPE>(&vars::nc::count_pions_double, &CUT, &TCUT)});
    vars_sig_1g1p.insert({"true_primary_pions", SpineVar<TTYPE,RTYPE>(&vars::nc::count_pions_double, &CUT, &TCUT)});
    vars_sig_1g1p.insert({"reco_primary_protons", SpineVar<RTYPE,RTYPE>(&vars::nc::count_protons_double, &CUT, &TCUT)});
    vars_sig_1g1p.insert({"true_primary_protons", SpineVar<TTYPE,RTYPE>(&vars::nc::count_protons_double, &CUT, &TCUT)});

    analysis.AddTree("signalNu_1photon1proton", vars_sig_1g1p, true);

    /**
     * @brief Run the analysis.
     * @details This runs the analysis on the samples specified by the
     * SpectrumLoaders and variables added to the Analysis class. It loops over
     * each sample (here only one), applies the cuts and variables to the data,
     * and stores the results in a TFile.
     */
    analysis.Go();
}

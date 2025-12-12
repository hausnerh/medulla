/**
 * @file run_analusis_gOre.cc
 * @brief Analyze the selection produced by gOre.toml
 * @details Make plots of the different variables of interest in terms of the categories
 * @author hhausner@fnal.gov
 **/

#include "TError.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include <filesystem>

#include "include/analysis_tree.h"
#include "include/yell_try_die.h"

inline const double pion_mass = 134.9768;

double landauPlusGauss(double* x, double* par)
{
  return par[0] * TMath::Landau(x[0], par[1], par[2], true) + par[3] * TMath::Gaus(x[0], par[4], par[5], true);
}

inline std::vector<std::pair<std::string, ana::tools::cut_sequence>> sel_cats =
{
  {"1#gammaXp (Fid/Con)",     "true_category == 0"},
  {"1#gammaXp (Not Fid/Con)", "true_category == 1"},
  {"1eXp (Fid/Con)",          "true_category == 2"},
  {"1eXp (Not Fid/Con)",      "true_category == 3"},
  {"2#gamma (Fid/Con)",       "true_category == 4"},
  {"2#gamma (Not Fid/Con)",   "true_category == 5"},
  {"Other NC",                "true_category == 6"},
  {"Other CC",                "true_category == 7"},
  {"Cosmic",                  "true_category == 8"}
};

//inline std::vector<std::pair<std::string, ana::tools::cut_sequence>> sig_cats =
//{
//  {"NC #Delta#rightarrowN#gamma", "true_mc_category == 0"},
//  {"Other NC 1#gammaXp Post-FSI", "true_mc_category == 1"},
//  {"NC #pi^{0}_{}",               "true_mc_category == 2"},
//  {"NC #pi^{+/-}_{}",             "true_mc_category == 3"},
//  {"Other NC",                    "true_mc_category == 4"},
//  {"CC e",                        "true_mc_category == 5"},
//  {"Other CC",                    "true_mc_category == 6"}
//  //{"Cosmic",                      "true_mc_category == 7"}
//};

inline std::map<int, std::string> genie_modes =
{
  { -1, "Unknown"},
  {  0, "Quasi-Elastic"},
  {  1, "Resonance"},
  {  2, "DIS"},
  {  3, "Coherent"},
  {  4, "Coherenet Elastic"},
  {  5, "Electron Scattering"},
  {  6, "IMD Annihilation"},
  {  7, "IBD"},
  {  8, "Glashow Resonance"},
  {  9, "AM #nu-#gamma"},
  { 10, "MEC"},
  { 11, "Diffractive"},
  { 12, "EM"},
  { 13, "Weak Mixing"}
};

inline std::map<int, std::string> genie_inttypes =
{
  {  -1, "Unknown"},
  {1000, "Nuance Offset"},
  {1001, "CC QE"},
  {1002, "NC QE"},
  {1003, "CC Res #nu p #pi^{+}"},
  {1004, "CC Res #nu p #pi^{0}"},
  {1005, "CC Res #nu n #pi^{+}"},
  {1006, "NC Res #nu p #pi^{0}"},
  {1007, "NC Res #nu n #pi^{+}"},
  {1008, "NC Res #nu n #pi^{0}"},
  {1009, "NC Res #nu p #pi^{-}"},
  {1010, "CC Res anti-#nu n #pi^{-}"},
  {1011, "CC Res anti-#nu n #pi^{0}"},
  {1012, "CC Res anti-#nu p #pi^{-}"},
  {1013, "NC Res anti-#nu p #pi^{0}"},
  {1014, "NC Res anti-#nu n #pi^{+}"},
  {1015, "NC Res anti-#nu n #pi^{0}"},
  {1016, "NC Res anti-#nu p #pi^{-}"},
  {1017, "CC Res #nu #Delta^{-}"},
  {1021, "CC Res #nu #Delta^{++} #pi^{-}"},
  {1028, "CC Res anti-#nu #Delta^{0} #pi^{-}"},
  {1032, "CC Res anti-#nu #Delta^{-} #pi^{+}"},
  {1039, "CC Res #nu p #rho^{+}"},
  {1041, "CC Res #nu n #rho^{+}"},
  {1046, "CC Res anti-#nu n #rho^{-}"},
  {1048, "CC Res anti-#nu n #rho^{0}"},
  {1053, "CC Res #nu #Sigma^{+} #Kappa^{+}"},
  {1055, "CC Res #nu #Sigma^{+} #Kappa^{0}"},
  {1060, "CC Res anti-#nu #Sigma^{-} #Kappa^{0}"},
  {1062, "CC Res anti-#nu #Sigma^{0} #Kappa^{0}"},
  {1067, "CC Res #nu p #eta"},
  {1070, "CC Res anti-#nu n #eta"},
  {1073, "CC Res #nu #Kappa^{+} #Lambda^{0}"},
  {1076, "CC Res anti-#nu #Kappa^{0} #Lambda^{0}"},
  {1079, "CC Res #nu p #pi^{+} #pi^{-}"},
  {1080, "CC Res #nu p #pi^{0} #pi^{0}"},
  {1085, "CC Res anti-#nu n #pi^{+} #pi^{-}"},
  {1086, "CC Res anti-#nu n #pi^{0} #pi^{0}"},
  {1090, "CC Res anti-#nu p #pi^{0} #pi^{0}"},
  {1091, "CC DIS"},
  {1092, "NC DIS"},
  {1093, "Unused 1"},
  {1094, "Unused 2"},
  {1095, "CC QE Hyperon"},
  {1096, "NC Coherent"},
  {1097, "CC Coherent"},
  {1098, "Electron Scattering"},
  {1099, "IMD"},
  {1100, "MEC 2p2h"}
};

inline std::map<int, std::string> res_codes =
{
  {-1, "No Resonance"},
  { 0, "#Delta(1232)"},
  { 1, "#Nu(1535)"},
  { 2, "#Nu(1520)"},
  { 3, "#Nu(1650)"},
  { 4, "#Delta(1700)"},
  { 5, "#Nu(1675)"},
  { 6, "#Delta(1620)"},
  { 7, "#Delta(1700)"},
  { 8, "#Nu(1440)"},
  { 9, "#Delta(1600)"},
  {10, "#Nu(1720)"},
  {11, "#Nu(1680)"},
  {12, "#Delta(1910)"},
  {13, "#Delta(1920)"},
  {14, "#Delta(1905)"},
  {15, "#Delta(1950)"},
  {16, "#Nu(1710)"},
  {17, "#Nu(1970)"}
};

//inline std::map<int, std::string> mc_cats =
//{
//  {0, "NC #Delta#rightarrowN#gamma (1#gammaXp Post-FSI Topology)"},
//  {1, "NC #Delta#rightarrowN#gamma (Other Post-FSI Topology)"},
//  {2, "Other NC 1#gammaXp"},
//  {3, "NC #pi^{0}_{} No Primary #gamma"},
//  {4, "Other NC"},
//  {5, "CC"},
//  {6, "Cosmic"}
//};
//inline std::map<int, std::string> mc_cats =
//{
//  {0, sig_cats.at(0).first},
//  {1, sig_cats.at(1).first},
//  {2, sig_cats.at(2).first},
//  {3, sig_cats.at(3).first},
//  {4, sig_cats.at(4).first},
//  {5, sig_cats.at(5).first},
//  {6, sig_cats.at(6).first}
//  //{7, sig_cats.at(7).first}
//};

inline std::vector<std::tuple<std::string, bool, double, double>> vars_to_optimize =
{
  //{"reco_flash_total_pe",                          false,  0, 20000},
  //{"reco_leading_primary_gOre_dpT",                true,   0,  1000},
  {"reco_leading_primary_gOre_primary_softmax",    false,  0,     1},
  {"reco_leading_primary_gOre_start_dedx",         false,  0,     5},
  {"reco_leading_primary_gOre_axial_spread",       false, -0.5,   1},
  {"reco_leading_primary_gOre_directional_spread", true,   0,     1},
  {"reco_leading_primary_gOre_photon_softmax",     false,  0,     1}
  //{"reco_gOre_gap",                                false,  0,   100}
};

inline std::vector<std::tuple<std::string, std::string, bool, double, double>> vars_to_optimize_with_condition =
{
  //{"reco_n_protons > 1",  "reco_gOre_score", false, -1, 1},
  //{"reco_n_protons == 0", "reco_gOre_score", false, -1, 1}
  //{"reco_n_protons == 0", "reco_flash_total_pe", false, 0, 20000},
  {"reco_pion_mass > 0", "reco_pion_mass", true, 0, 300},
  {"reco_subleading_primary_gOre_shower_ke > 0", "reco_subleading_primary_gOre_primary_softmax", true, 0, 1},
  {"reco_subleading_primary_gOre_shower_ke > 0", "reco_subleading_primary_gOre_shower_ke", true, 0, 25},
  {"reco_n_protons >= 1", "reco_gOre_gap", false, 0, 100}
};

std::tuple<ana::tools::cut_sequence, ana::tools::cut_sequence> optimize_threshold(ana::tools::analysis_tree& the_analysis_tree, // NOT CONST TO UPDATE SIGNAL DEFS
                                                                                  const ana::tools::cut_sequence& old_cut,
                                                                                  const ana::tools::cut_sequence& old_truth_cut,
                                                                                  const std::string& var,
                                                                                  const std::string& truth_var,
                                                                                  const double& lw_end, const double up_end,
                                                                                  const std::string& sample, const std::string& pdf_suffix)
{
  std::string report = "Optiminze " + var + " Threshold";
  std::cout << report << std::endl;
  auto [cut_opt, cut_plot]
    = try_call(report,
        [&the_analysis_tree, &old_cut, &old_truth_cut, &var, &truth_var, &lw_end, &up_end]
        {
          return the_analysis_tree.optimize_threshold(var, truth_var, lw_end, up_end, old_cut, old_truth_cut);
        });
  ana::tools::cut_sequence cut = cut_opt.cut;
  ana::tools::cut_sequence truth_cut = cut_opt.truth_cut;
  try_call(cut.string(), [&the_analysis_tree, &cut]{ the_analysis_tree.report_on_cut(cut); });
  cut_plot.sc.Print("plots/"+sample+"/"+var+"_optimized"+pdf_suffix);
  cut_opt.Print("plots/"+sample+"/"+var+"_FOM"+pdf_suffix);
  return std::tie(cut, truth_cut);
}

ana::tools::cut_sequence optimize_cut(const ana::tools::analysis_tree& the_analysis_tree,
                                      const ana::tools::cut_sequence& old_cut,
                                      const std::string& var, const bool& upper,
                                      const double& lw_end, const double up_end,
                                      const std::string& sample, const std::string& pdf_suffix)
{
  std::string report = (upper) ? "Optimize Upper Bound on " + var :
                                 "Optimize Lower Bound on " + var;
  std::cout << report << std::endl;
  auto [cut_opt, cut_plot]
    = try_call(report,
        [&the_analysis_tree, &old_cut, &var, &upper, &lw_end, &up_end]
        {
          return the_analysis_tree.optimize_bound(var, lw_end, up_end, old_cut, upper);
        });
  ana::tools::cut_sequence cut = cut_opt.cut;
  try_call(cut.string(), [&the_analysis_tree, &cut]{ the_analysis_tree.report_on_cut(cut); });
  cut_plot.sc.Print("plots/"+sample+"/"+var+"_optimized"+pdf_suffix);
  cut_opt.Print("plots/"+sample+"/"+var+"_FOM"+pdf_suffix);
  return cut;
}

ana::tools::cut_sequence optimize_conditional_cut(const ana::tools::analysis_tree& the_analysis_tree,
                                                  const ana::tools::cut_sequence& old_cut,
                                                  const std::string& condition, const std::string& var,
                                                  const bool& upper, const double& lw_end, const double up_end,
                                                  const std::string& sample, const std::string& pdf_suffix)
{
  std::string report = (upper) ? "Optimize Upper Bound on " + var + " Under Condition " + condition:
                                 "Optimize Lower Bound on " + var + " Under Condition " + condition;
  std::cout << report << std::endl;
  auto [cut_opt, cut_plot]
    = try_call(report,
        [&the_analysis_tree, &old_cut, &condition, &var, &upper, &lw_end, &up_end]
        {
          return the_analysis_tree.optimize_conditional_bound(condition, var, lw_end, up_end, old_cut, upper);
        });
  ana::tools::cut_sequence cut = cut_opt.cut;
  try_call(cut.string(), [&the_analysis_tree, &cut]{ the_analysis_tree.report_on_cut(cut); });
  cut_plot.sc.Print("plots/"+sample+"/"+var+"_optimized_with_"+condition+pdf_suffix);
  cut_opt.Print("plots/"+sample+"/"+var+"_FOM_with_"+condition+pdf_suffix);
  return cut;
}

int run_analysis(const std::string& fileName,
                 const std::string& sampleBase,
                 const std::string& sampleName,
                 const std::string& topology,
                 const std::string& sel,
                 const std::string& sig,
                 const bool& optimizeCuts)
{
  // sample string
  const std::string sample = sampleBase + "/" + sampleName + "/" + topology;

  // setup output dir
  if (not std::filesystem::is_directory("plots") || not std::filesystem::exists("plots"))
    std::filesystem::create_directory("plots");
  if (not std::filesystem::is_directory("plots/"+sampleBase) || not std::filesystem::exists("plots/"+sampleBase))
    std::filesystem::create_directory("plots/"+sampleBase);
  if (not std::filesystem::is_directory("plots/"+sampleBase+"/"+sampleName) || not std::filesystem::exists("plots/"+sampleBase+"/"+sampleName))
    std::filesystem::create_directory("plots/"+sampleBase+"/"+sampleName);
  if (not std::filesystem::is_directory("plots/"+sample) || not std::filesystem::exists("plots/"+sample))
    std::filesystem::create_directory("plots/"+sample);
 
  // setup up categories
  const std::vector<std::pair<std::string, ana::tools::cut_sequence>> sig_cats = (topology == "1g0p") ?
    std::vector<std::pair<std::string, ana::tools::cut_sequence>>({
      {"NC #Delta#rightarrowN#gamma (1#gamma0p)",      "(true_mc_category == 0) && (true_n_protons == 0)"},
      {"NC #Delta#rightarrowN#gamma (Other Topology)", "(true_mc_category == 0) && (true_n_protons != 0)"},
      {"Other NC 1#gammaXp Post-FSI",                  "true_mc_category == 1"},
      {"NC #pi^{0}_{} (#Delta Res)",                   "true_mc_category == 2"},
      {"NC #pi^{0}_{} (Other)",                        "true_mc_category == 3"},
      {"NC #pi^{+/-}_{}",                              "true_mc_category == 4"},
      {"Other NC",                                     "true_mc_category == 5"},
      {"CC e",                                         "true_mc_category == 6"},
      {"Other CC",                                     "true_mc_category == 7"}
    }) :
    std::vector<std::pair<std::string, ana::tools::cut_sequence>>({
      {"NC #Delta#rightarrowN#gamma (1#gamma1p)",      "(true_mc_category == 0) && (true_n_protons == 1)"},
      {"NC #Delta#rightarrowN#gamma (Other Topology)", "(true_mc_category == 0) && (true_n_protons != 1)"},
      {"Other NC 1#gammaXp Post-FSI",                  "true_mc_category == 1"},
      {"NC #pi^{0}_{} (#Delta Res)",                   "true_mc_category == 2"},
      {"NC #pi^{0}_{} (Other)",                        "true_mc_category == 3"},
      {"NC #pi^{+/-}_{}",                              "true_mc_category == 4"},
      {"Other NC",                                     "true_mc_category == 5"},
      {"CC e",                                         "true_mc_category == 6"},
      {"Other CC",                                     "true_mc_category == 7"}
    });
  const std::map<int, std::string> mc_cats =
  {
    {0, sig_cats.at(0).first},
    {1, sig_cats.at(1).first},
    {2, sig_cats.at(2).first},
    {3, sig_cats.at(3).first},
    {4, sig_cats.at(4).first},
    {5, sig_cats.at(5).first},
    {6, sig_cats.at(6).first},
    {7, sig_cats.at(7).first},
    {8, sig_cats.at(8).first}
  };

  // setup analysis tree
  std::string fileLocation = std::filesystem::current_path().string();
  std::string directoryName = "events/" + sampleName;
  std::string selTreeName = "selected_" + topology;
  std::string sigTreeName = "signal_" + topology;
  ana::tools::cut_sequence selection_cut = "reco_gOre_topology";
  ana::tools::cut_sequence signal_cut = "true_mc_category == 0";
  ana::tools::cut_sequence reco_topology_cut = (topology == "1g0p") ? "reco_n_protons == 0" : "reco_n_protons == 1";
  ana::tools::cut_sequence true_topology_cut = (topology == "1g0p") ? "true_n_protons == 0" : "true_n_protons == 1";
  selection_cut += reco_topology_cut;
  signal_cut += true_topology_cut;
  ana::tools::analysis_tree my_analysis_tree(fileLocation+"/"+fileName, directoryName,
                                             selection_cut, signal_cut,
                                             selTreeName, sigTreeName,
                                             sig_cats, sig_cats);
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~VARIABLE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~BINS~~~~MIN~~~~~MAX~~~~~TITLE
  //my_analysis_tree.add_variable("true_gOre_topology",                                2,    -0.5,    1.5,  "True Single Shower Topology");
  //my_analysis_tree.add_variable("reco_gOre_topology",                                2,    -0.5,    1.5,  "Reconstucted Single Shower Topology");
  my_analysis_tree.add_variable("reco_flash_total_pe",                              50,     0, 100000,    "Total Flash PE");
  //my_analysis_tree.add_variable("true_n_primary_showers",                            6,     0,      6,    "True Primary Showers above Threshold");
  //my_analysis_tree.add_variable("reco_n_primary_showers",                            6,     0,      6,    "Primary Showers above Threshold");
  //my_analysis_tree.add_variable("true_n_secondary_showers",                          6,     0,      6,    "True Secondary Showers above Threshold");
  //my_analysis_tree.add_variable("reco_n_secondary_showers",                          6,     0,      6,    "Secondary Showers above Threshold");
  my_analysis_tree.add_variable("true_n_protons",                                    6,     0,      6,    "True Protons above Threshold");
  my_analysis_tree.add_variable("reco_n_protons",                                    6,     0,      6,    "Protons above Threshold");
  my_analysis_tree.add_variable("true_vertex_x",                                    50,  -400,    400,    "True Vertex X (cm)");
  my_analysis_tree.add_variable("true_vertex_y",                                    50,  -200,    200,    "True Vertex Y (cm)");
  my_analysis_tree.add_variable("true_vertex_z",                                    50, -1000,   1000,    "True Vertex Z (cm)");
  my_analysis_tree.add_variable("reco_vertex_x",                                    50,  -400,    400,    "Vertex X (cm)");
  my_analysis_tree.add_variable("reco_vertex_y",                                    50,  -200,    200,    "Vertex Y (cm)");
  my_analysis_tree.add_variable("reco_vertex_z",                                    50, -1000,   1000,    "Vertex Z (cm)");
  //my_analysis_tree.add_variable("true_xy_wall_dist",                                50,     0,    200,    "True Minimum Distance from Vertex to X-/Y-side Detector Wall");
  //my_analysis_tree.add_variable("true_z_wall_dist",                                 50,     0,   1000,    "True Minimum Distance from Vertex to Z-side Detector Wall");
  //my_analysis_tree.add_variable("reco_xy_wall_dist",                                50,     0,    200,    "Minimum Distance from Vertex to X-/Y-side Detector Wall");
  //my_analysis_tree.add_variable("reco_z_wall_dist",                                 50,     0,   1000,    "Minimum Distance from Vertex to Z-side Detector Wall");
  my_analysis_tree.add_variable("reco_gOre_gap",                                    50,     0,    100,    "#gamma-candiate Distance from Vertex (cm)");
  my_analysis_tree.add_variable("reco_pion_mass",                                   50,     0,    250,    "Reconstructed Neutral Pion Mass Peak (MeV/c^{2}_{})");
  my_analysis_tree.add_variable("reco_delta_mass_P",                                25,   950,   1600,    "Reconstructed M_{#Delta} Resonance Peak (MeV/c^{2}_{})");
  my_analysis_tree.add_variable("reco_delta_mass_N",                                25,   950,   1600,    "Reconstructed M_{#Delta} Resonance Peak (MeV/c^{2}_{})");
  my_analysis_tree.add_variable("true_delta_mass_P",                                25,   950,   1600,    "True M_{#Delta} Resonance Peak (MeV/c^{2}_{})");
  my_analysis_tree.add_variable("true_delta_mass_N",                                25,   950,   1600,    "True M_{#Delta} Resonance Peak (MeV/c^{2}_{})");
  //my_analysis_tree.add_variable("true_detla_mass",                                  50,   800,   1800,    "True M_{#Delta} Resonance Peak (MeV/c^{2}_{})");
  my_analysis_tree.add_variable("true_interaction_mode",                            15,    -1.5,   13.5,  "GENIE Interaction Mode");
  my_analysis_tree.add_variable("true_interaction_type",                           101,   999.5, 1100.5,  "GENIE Interaction Type");
  my_analysis_tree.add_variable("true_baryon_res_code",                             19,    -1.5,   17.5,  "Resonance Number");
  my_analysis_tree.add_variable("true_mc_category",                                  8,    -0.5,    7.5,  "MC Truth Category"); 
  my_analysis_tree.add_variable("true_category",                                     9,    -0.5,    8.5,  "Topological Truth Category"); 
  my_analysis_tree.add_variable("reco_leading_primary_gOre_primary_softmax",       100,     0,      1,    "Shower Primary Softmax");
  my_analysis_tree.add_variable("reco_leading_primary_gOre_photon_softmax",         50,     0,      1,    "Shower Photon Softmax");
  //my_analysis_tree.add_variable("reco_leading_primary_gOre_electron_softmax",       50,     0,      1,    "Shower Electron Softmax");
  my_analysis_tree.add_variable("reco_leading_primary_gOre_start_dedx",             50,     0,     25,    "Shower Start dE/dx (MeV/cm)");
  my_analysis_tree.add_variable("true_leading_primary_gOre_azimuthal_angle",        50,    -3.14,   3.14, "True Shower Azimuthal Angle (rad)");
  my_analysis_tree.add_variable("reco_leading_primary_gOre_azimuthal_angle",        50,    -3.14,   3.14, "Shower Azimuthal Angle (rad)");
  my_analysis_tree.add_variable("true_leading_primary_gOre_polar_angle",            50,     0,      3.14, "True Shower Polar Angle (rad)");
  my_analysis_tree.add_variable("reco_leading_primary_gOre_polar_angle",            50,     0,      3.14, "Shower Polar Angle (rad)");
  my_analysis_tree.add_variable("true_leading_primary_gOre_shower_ke",              50,     0,   1000,    "True Shower KE (MeV)");
  my_analysis_tree.add_variable("reco_leading_primary_gOre_shower_ke",              50,     0,   1000,    "Shower KE (MeV)");
  my_analysis_tree.add_variable("true_leading_primary_gOre_shower_dpT",            100,     0,   1000,    "True Shower Transverse Momentum (MeV/c)");
  my_analysis_tree.add_variable("reco_leading_primary_gOre_shower_dpT",            100,     0,   1000,    "Shower Transverse Momentum (MeV/c)");
  my_analysis_tree.add_variable("reco_leading_primary_gOre_axial_spread",           75,    -0.5,    1,    "Shower Axial Spread");
  my_analysis_tree.add_variable("reco_leading_primary_gOre_directional_spread",    100,     0,      1,    "Shower Directional Spread");
  my_analysis_tree.add_variable("true_leading_primary_gOre_iou",                   100,     0,      1,    "True Shower IoU");
  my_analysis_tree.add_variable("reco_leading_primary_gOre_iou",                   100,     0,      1,    "Shower IoU");
  my_analysis_tree.add_variable("reco_subleading_primary_gOre_primary_softmax",    100,     0,      1,    "Subleading Shower Primary Softmax");
  my_analysis_tree.add_variable("reco_subleading_primary_gOre_photon_softmax",      50,     0,      1,    "Subleading Shower Photon Softmax");
  //my_analysis_tree.add_variable("reco_subleading_primary_gOre_electron_softmax",    50,     0,      1,    "Subleading Shower Electron Softmax");
  my_analysis_tree.add_variable("reco_subleading_primary_gOre_start_dedx",          50,     0,     25,    "Subleading Shower Start dE/dx (MeV/cm)");
  my_analysis_tree.add_variable("true_subleading_primary_gOre_azimuthal_angle",     50,    -3.14,   3.14, "True Subleading Shower Azimuthal Angle (rad)");
  my_analysis_tree.add_variable("reco_subleading_primary_gOre_azimuthal_angle",     50,    -3.14,   3.14, "Subleading Shower Azimuthal Angle (rad)");
  my_analysis_tree.add_variable("true_subleading_primary_gOre_polar_angle",         50,     0,      3.14, "True Subleading Shower Polar Angle (rad)");
  my_analysis_tree.add_variable("reco_subleading_primary_gOre_polar_angle",         50,     0,      3.14, "Subleading Shower Polar Angle (rad)");
  my_analysis_tree.add_variable("true_subleading_primary_gOre_shower_ke",           25,     0,     25,    "True Subleading Shower KE (MeV)");
  my_analysis_tree.add_variable("reco_subleading_primary_gOre_shower_ke",           25,     0,     25,    "Subleading Shower KE (MeV)");
  my_analysis_tree.add_variable("true_subleading_primary_gOre_shower_dpT",          25,     0,     25,    "True Subleading Shower Transverse Momentum (MeV/c)");
  my_analysis_tree.add_variable("reco_subleading_primary_gOre_shower_dpT",          25,     0,     25,    "Subleading Shower Transverse Momentum (MeV/c)");
  my_analysis_tree.add_variable("reco_subleading_primary_gOre_axial_spread",        75,    -0.5,    1,    "Subleading Shower Axial Spread");
  my_analysis_tree.add_variable("reco_subleading_primary_gOre_directional_spread", 100,     0,      1,    "Subleading Shower Directional Spread");
  my_analysis_tree.add_variable("true_subleading_primary_gOre_iou",                100,     0,      1,    "True Subleading Shower IoU");
  my_analysis_tree.add_variable("reco_subleading_primary_gOre_iou",                100,     0,      1,    "Subleading Shower IoU");
  my_analysis_tree.add_variable("reco_leading_primary_proton_primary_softmax",     100,     0,      1,    "Leading Proton Primary Softmax");
  my_analysis_tree.add_variable("reco_leading_primary_proton_proton_softmax",       50,     0,      1,    "Leading Proton Photon Softmax");
  my_analysis_tree.add_variable("true_leading_primary_proton_azimuthal_angle",      50,    -3.14,   3.14, "True Proton Azimuthal Angle (rad)");
  my_analysis_tree.add_variable("reco_leading_primary_proton_azimuthal_angle",      50,    -3.14,   3.14, "Proton Azimuthal Angle (rad)");
  my_analysis_tree.add_variable("true_leading_primary_proton_polar_angle",          50,     0,      3.14, "True Proton Polar Angle (rad)");
  my_analysis_tree.add_variable("reco_leading_primary_proton_polar_angle",          50,     0,      3.14, "Proton Polar Angle (rad)");
  my_analysis_tree.add_variable("true_leading_primary_proton_ke",                   50,     0,   1000,    "True Proton KE (MeV)");
  my_analysis_tree.add_variable("reco_leading_primary_proton_ke",                   50,     0,   1000,    "Proton KE (MeV)");
  my_analysis_tree.add_variable("true_leading_primary_proton_dpT",                 100,     0,   1000,    "True Proton Transverse Momentum (MeV/c)");
  my_analysis_tree.add_variable("reco_leading_primary_proton_dpT",                 100,     0,   1000,    "Proton Transverse Momentum (MeV/c)");
  my_analysis_tree.add_variable("true_leading_primary_proton_length",              100,     0,    200,    "True Proton Length (cm)");
  my_analysis_tree.add_variable("reco_leading_primary_proton_length",              100,     0,    200,    "Proton Length (cm)");

  std::string pdf_suffix = ".pdf";
  ana::tools::cut_sequence cut;
  ana::tools::cut_sequence truth_cut;
  std::cout
    << "************************\n"
    << "* "    << sel <<     " *\n"
    << "************************" << std::endl;
  //*** TOPOLOGY ***//
  // should alread be implemented as part of making the selected TTree,
  // but we do this for completeness
  std::cout << "//*** TOPOLOGY ***//" << std::endl;
  cut += "reco_gOre_topology == 1";
  cut += reco_topology_cut;
  try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
  //*** ENERGY THRESH ***/
  //std::tie(cut, truth_cut) = optimize_threshold(my_analysis_tree, cut, truth_cut, "reco_leading_primary_gOre_ke", "true_leading_primary_gOre_ke", 0, 500, sample, pdf_suffix);
  //*** FLASH CUT ***//
  std::cout << "//*** FLASH CUT ***//" << std::endl;
  cut += "reco_flash == 1";
  try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
  //*** FIDUCIAL CUT ***//
  std::cout << "//*** FIDUCIAL CUT ***//" << std::endl;
  cut += "reco_fiducial == 1";
  try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
  //*** CONTAINMENT CUT***//
  std::cout << "//*** CONTAINMENT CUT ***//" << std::endl;
  cut += "reco_containment == 1";
  try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });

  if (optimizeCuts)
  {
    for (auto const& [var, upper, lw_end, up_end] : vars_to_optimize)
      cut = optimize_cut(my_analysis_tree, cut, var, upper, lw_end, up_end, sample, pdf_suffix);
    auto plot_gOre_shower_ke_electron_cuts = 
      try_call("plot_gOre_shower_ke_electron_cuts",
        [&my_analysis_tree, &cut]{ return my_analysis_tree.plot_var_sel("reco_leading_primary_gOre_shower_ke", cut); });
    plot_gOre_shower_ke_electron_cuts.PrintPreliminary("plots/"+sample+"/electron_cuts_reco_leading_primary_gOre_shower_ke"+pdf_suffix);
    auto plot_n_protons_electron_cuts = 
      try_call("plot_n_protons_electron_cuts",
        [&my_analysis_tree, &cut]{ return my_analysis_tree.plot_var_sel("reco_n_protons", cut); });
    plot_n_protons_electron_cuts.PrintPreliminary("plots/"+sample+"/electron_cuts_reco_n_protons"+pdf_suffix);
    for (auto const& [condition, var, upper, lw_end, up_end] : vars_to_optimize_with_condition)
      cut = optimize_conditional_cut(my_analysis_tree, cut, condition, var, upper, lw_end, up_end, sample, pdf_suffix);
    auto plot_gOre_shower_ke_pion_cuts = 
      try_call("plot_gOre_shower_ke_pion_cuts",
        [&my_analysis_tree, &cut]{ return my_analysis_tree.plot_var_sel("reco_leading_primary_gOre_shower_ke", cut); });
    plot_gOre_shower_ke_pion_cuts.PrintPreliminary("plots/"+sample+"/pion_cuts_reco_leading_primary_gOre_shower_ke"+pdf_suffix);
    auto plot_n_protons_pion_cuts = 
      try_call("plot_n_protons_pion_cuts",
        [&my_analysis_tree, &cut]{ return my_analysis_tree.plot_var_sel("reco_n_protons", cut); });
    plot_n_protons_pion_cuts.PrintPreliminary("plots/"+sample+"/pion_cuts_reco_n_protons"+pdf_suffix);
  }
  else
  {
    auto precut_leading_primary_gOre_primary_softmax =
      try_call("leading_primary_gOre_primary_softmax",
         [&my_analysis_tree, &cut]{ return my_analysis_tree.plot_var_sel("reco_leading_primary_gOre_primary_softmax", cut); });
    precut_leading_primary_gOre_primary_softmax.PrintPreliminary("plots/"+sample+"/precut_leading_primary_gOre_primary_softmax"+pdf_suffix);

    auto precut_leading_primary_gOre_start_dedx =
      try_call("leading_primary_gOre_start_dedx",
         [&my_analysis_tree, &cut]{ return my_analysis_tree.plot_var_sel("reco_leading_primary_gOre_start_dedx", cut); });
    precut_leading_primary_gOre_start_dedx.PrintPreliminary("plots/"+sample+"/precut_leading_primary_gOre_start_dedx"+pdf_suffix);

    auto precut_leading_primary_gOre_axial_spread =
      try_call("leading_primary_gOre_axial_spread",
         [&my_analysis_tree, &cut]{ return my_analysis_tree.plot_var_sel("reco_leading_primary_gOre_axial_spread", cut); });
    precut_leading_primary_gOre_axial_spread.PrintPreliminary("plots/"+sample+"/precut_leading_primary_gOre_axial_spread"+pdf_suffix);

    auto precut_leading_primary_gOre_directional_spread =
      try_call("leading_primary_gOre_directional_spread",
         [&my_analysis_tree, &cut]{ return my_analysis_tree.plot_var_sel("reco_leading_primary_gOre_directional_spread", cut); });
    precut_leading_primary_gOre_directional_spread.PrintPreliminary("plots/"+sample+"/precut_leading_primary_gOre_directional_spread"+pdf_suffix);

    auto precut_leading_primary_gOre_photon_softmax =
      try_call("leading_primary_gOre_photon_softmax",
         [&my_analysis_tree, &cut]{ return my_analysis_tree.plot_var_sel("reco_leading_primary_gOre_photon_softmax", cut); });
    precut_leading_primary_gOre_photon_softmax.PrintPreliminary("plots/"+sample+"/precut_leading_primary_gOre_photon_softmax"+pdf_suffix);

    std::cout << "//*** Electron/Photon Seperation ***//" << std::endl;
    cut += "reco_leading_primary_gOre_primary_softmax > 0.992";
    cut += "reco_leading_primary_gOre_start_dedx > 3.495";
    cut += "reco_leading_primary_gOre_axial_spread > -0.1895";
    cut += "reco_leading_primary_gOre_directional_spread < 0.113";
    cut += "reco_leading_primary_gOre_photon_softmax > 0.077";
    //cut += "reco_leading_primary_gOre_start_dedx > 3.645";
    //cut += "reco_leading_primary_gOre_axial_spread > -0.005";
    //cut += "reco_leading_primary_gOre_directional_spread < 0.068";
    //cut += "reco_leading_primary_gOre_photon_softmax > 0.092";
    try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });

    auto postcut_leading_primary_gOre_primary_softmax =
      try_call("leading_primary_gOre_primary_softmax",
         [&my_analysis_tree, &cut]{ return my_analysis_tree.plot_var_sel("reco_leading_primary_gOre_primary_softmax", cut); });
    postcut_leading_primary_gOre_primary_softmax.PrintPreliminary("plots/"+sample+"/postcut_leading_primary_gOre_primary_softmax"+pdf_suffix);

    auto postcut_leading_primary_gOre_start_dedx =
      try_call("leading_primary_gOre_start_dedx",
         [&my_analysis_tree, &cut]{ return my_analysis_tree.plot_var_sel("reco_leading_primary_gOre_start_dedx", cut); });
    postcut_leading_primary_gOre_start_dedx.PrintPreliminary("plots/"+sample+"/postcut_leading_primary_gOre_start_dedx"+pdf_suffix);

    auto postcut_leading_primary_gOre_axial_spread =
      try_call("leading_primary_gOre_axial_spread",
         [&my_analysis_tree, &cut]{ return my_analysis_tree.plot_var_sel("reco_leading_primary_gOre_axial_spread", cut); });
    postcut_leading_primary_gOre_axial_spread.PrintPreliminary("plots/"+sample+"/postcut_leading_primary_gOre_axial_spread"+pdf_suffix);

    auto postcut_leading_primary_gOre_directional_spread =
      try_call("leading_primary_gOre_directional_spread",
         [&my_analysis_tree, &cut]{ return my_analysis_tree.plot_var_sel("reco_leading_primary_gOre_directional_spread", cut); });
    postcut_leading_primary_gOre_directional_spread.PrintPreliminary("plots/"+sample+"/postcut_leading_primary_gOre_directional_spread"+pdf_suffix);

    auto postcut_leading_primary_gOre_photon_softmax =
      try_call("leading_primary_gOre_photon_softmax",
         [&my_analysis_tree, &cut]{ return my_analysis_tree.plot_var_sel("reco_leading_primary_gOre_photon_softmax", cut); });
    postcut_leading_primary_gOre_photon_softmax.PrintPreliminary("plots/"+sample+"/postcut_leading_primary_gOre_photon_softmax"+pdf_suffix);

    // *** Fit Pion Mass Peak ***

    auto precut_subleading_primary_gOre_ke =
      try_call("subleading_primary_gOre_ke",
         [&my_analysis_tree, &cut]{ return my_analysis_tree.plot_var_sel("reco_subleading_primary_gOre_shower_ke", cut); });
    precut_subleading_primary_gOre_ke.PrintPreliminary("plots/"+sample+"/precut_subleading_primary_gOre_shower_ke"+pdf_suffix);

    auto precut_pion_mass =
      try_call("pion_mass",
         [&my_analysis_tree, &cut]{ return my_analysis_tree.plot_var_sel("reco_pion_mass", cut); });
    precut_pion_mass.PrintPreliminary("plots/"+sample+"/precut_pion_mass"+pdf_suffix);
    precut_pion_mass.SumHist();
    TF1* pion_fit = new TF1("selected", landauPlusGauss, 0, 250, 6);
    pion_fit->SetParNames("LandauNorm", "LandauMean", "LandauSigma", "GausNorm", "GausMean", "GausSigma");
    pion_fit->SetParameters(precut_pion_mass.sum->GetMaximum(), 21, 7.8, 0.25*precut_pion_mass.sum->GetMaximum(), pion_mass, 10);
    pion_fit->SetParLimits(1, 0, 50);
    pion_fit->SetParLimits(2, 0, 250);
    pion_fit->SetParLimits(4, 50, 250);
    pion_fit->SetParLimits(5, 0, 250);
    TFitResultPtr res = precut_pion_mass.sum->Fit(pion_fit, "LS", "0 R");
    std::cout << "Chi^2: " << res->Chi2() << ", NDoF: " << res->Ndf() << std::endl;
    precut_pion_mass.canvas->cd();
    precut_pion_mass.stack->Draw("hist");
    pion_fit->Draw("L,same");
    precut_pion_mass.canvas->Update();
    precut_pion_mass.PrintPreliminary("plots/"+sample+"/precut_pion_mass_FIT"+pdf_suffix);

    std::cout << "//*** Pion Rejection ***//" << std::endl;
    cut.add_conditional_cut("reco_pion_mass > 0", "reco_pion_mass < 94.2");
    cut.add_conditional_cut("reco_subleading_primary_gOre_shower_ke > 0", "reco_subleading_primary_gOre_shower_ke < 6.8");
    try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });

    auto postcut_subleading_primary_gOre_ke =
      try_call("subleading_primary_gOre_shower_ke",
         [&my_analysis_tree, &cut]{ return my_analysis_tree.plot_var_sel("reco_subleading_primary_gOre_shower_ke", cut); });
    postcut_subleading_primary_gOre_ke.PrintPreliminary("plots/"+sample+"/postcut_subleading_primary_gOre_shower_ke"+pdf_suffix);

    auto postcut_pion_mass =
      try_call("pion_mass",
         [&my_analysis_tree, &cut]{ return my_analysis_tree.plot_var_sel("reco_pion_mass", cut); });
    postcut_pion_mass.PrintPreliminary("plots/"+sample+"/postcut_pion_mass"+pdf_suffix);

  }

  //*** Plot Vars ***//
  for (auto const& var : my_analysis_tree.variables())
  {
    auto var_plot =
      try_call("plot "+var,
        [&my_analysis_tree, &var, &cut]{ return my_analysis_tree.plot_var_sel(var, cut); });
    auto var_plot_sig =
     try_call("plot "+var+" signal",
       [&my_analysis_tree, &var, &cut]{ return my_analysis_tree.plot_var_sig(var, cut); });
    std::string pdfName = "plots/"+sample+"/"+var+pdf_suffix;
    std::string pdfName_sig = "plots/"+sample+"/signal_"+var+pdf_suffix;
    //if (var == var_delta_mass_N)
    //{
    //  pdfName = "plots/"+sample+"/reco_delta_mass_neutron_approx"+pdf_suffix;
    //  pdfName_sig = "plots/"+sample+"/signal_reco_delta_mass_neutron_approx"+pdf_suffix;
    //}
    // some vars use alphanumeric labels
    if (var == "true_interaction_type" ||
        var == "true_interaction_mode" ||
        var == "true_baryon_res_code"  ||
        var == "true_mc_category"      ||
        var == "true_category"          )
    {
      using StrFromDouble = std::function<std::string(double)>;
      StrFromDouble get_label = (var == "true_interaction_type") ? StrFromDouble([](const int& code){ return (genie_inttypes.count(code) == 1) ? genie_inttypes.at(code) : ""; })
                              : (var == "true_interaction_mode") ? StrFromDouble([](const int& code){ return genie_modes.at(code); })
                              : (var == "true_baryon_res_code")  ? StrFromDouble([](const int& code){ return res_codes.at(code); })
                              : (var == "true_mc_category")      ? StrFromDouble([&mc_cats](const int& code){
                                                                                                              if (code == 0)
                                                                                                              {
                                                                                                                return std::string("NC #Delta#rightarrowN#gamma");
                                                                                                              }
                                                                                                              else
                                                                                                                return mc_cats.at(code + 1);
                                                                                                            })
                              : (var == "true_category")         ? StrFromDouble([](const int& code){ return sel_cats.at(code).first; })
                              :                                    StrFromDouble([](const int& code){ return std::to_string(code); });
      size_t nBins = var_plot.hists.front()->GetXaxis()->GetNbins();
      for (size_t bin = 1; bin < nBins + 1; ++bin)
      {
        int bin_code = var_plot.hists.front()->GetXaxis()->GetBinCenter(bin);
        std::string bin_label = get_label(bin_code);
        for (auto& hist : var_plot.hists)
          hist->GetXaxis()->SetBinLabel(bin, bin_label.c_str());
        for (auto& hist : var_plot_sig.hists)
          hist->GetXaxis()->SetBinLabel(bin, bin_label.c_str());
      }
      std::string title_str = (var == "true_interaction_type") ? "GENIE Interaction Type;;Events"
                            : (var == "true_interaction_mode") ? "GENIE Interaction Mode;;Events"
                            : (var == "true_baryon_res_code")  ?              "Resonance;;Events"
                            : (var == "true_categoty")         ?   "Topological Category;;Events"
                            : (var == "true_mc_categoty")      ?            "MC Category;;Events"
                            :                                                          ";;Events";
      var_plot.stack->SetTitle(title_str.c_str());
      var_plot_sig.stack->SetTitle(title_str.c_str());
      if (var == "true_interaction_type")
      {
        var_plot.canvas->cd();
        var_plot.stack->GetXaxis()->SetLabelSize(0.015);
        var_plot.stack->GetXaxis()->LabelsOption("v");
        var_plot_sig.canvas->cd();
        var_plot_sig.stack->GetXaxis()->SetLabelSize(0.015);
        var_plot_sig.stack->GetXaxis()->LabelsOption("v");
      }
      var_plot.canvas->Update();
      var_plot_sig.canvas->Update();
    }
    // if Delta mass plot draw a line for the mass value of 1232 MeV/c^2
    if (var.find("delta") != std::string::npos)
    {
      var_plot.canvas->cd();
      var_plot.canvas->Update();
      double max = std::ceil(var_plot.stack->GetMaximum() + std::sqrt(var_plot.stack->GetMaximum()));
      TLine massLine(1232, 0, 1232, max);
      massLine.SetLineColor(kBlack);
      massLine.SetLineStyle(kDashed);
      massLine.Draw();
      var_plot.PrintPreliminary(pdfName);
      var_plot_sig.canvas->cd();
      var_plot_sig.canvas->Update();
      double max_sig = std::ceil(var_plot_sig.stack->GetMaximum() + std::sqrt(var_plot_sig.stack->GetMaximum()));
      TLine massLine_sig(1232, 0, 1232, max_sig);
      massLine_sig.SetLineColor(kBlack);
      massLine_sig.SetLineStyle(kDashed);
      massLine_sig.Draw();
      var_plot_sig.PrintPreliminary(pdfName_sig);
    }
    else
    {
      var_plot.PrintPreliminary(pdfName);
      var_plot_sig.PrintPreliminary(pdfName_sig);
    }
  }

  // try res peak plots
  //auto delta_mass_P =
  //  try_call("delta_mass_P",
  //     [&my_analysis_tree, &cut]{ return my_analysis_tree.plot_var_sel("reco_delta_mass", cut.with_addition("reco_n_protons == 1")); });
  //delta_mass_P.stack->SetTitle("1#gamma1p");
  //delta_mass_P.PrintPreliminary("plots/"+sample+"/delta_mass_P.pdf");
  //auto delta_mass_N =
  //  try_call("delta_mass_N",
  //     [&my_analysis_tree, &cut, &var_delta_mass_N]{ return my_analysis_tree.plot_var_sel(var_delta_mass_N, cut.with_addition("reco_n_protons == 0")); });
  //delta_mass_N.stack->SetTitle("1#gamma0p");
  //delta_mass_N.PrintPreliminary("plots/"+sample+"/delta_mass_N.pdf");

  return 0;
}

int main(int argc, char* argv[])
{
  gErrorIgnoreLevel=3000;
  //gDebug=3;
  int ret = 0;
  bool optimize_cuts = true;
  //bool optimize_cuts = false;
  if (not (argc > 1))
    die("Must pass in the name of the ROOT file to analyze. Bail.");
  if (not (argc > 2))
    die("Must give name of sample in file to process");
  if (not (argc > 3))
    die("Must give the name of the topology used");
  const std::string suffix = ".root";
  const std::string fileName(argv[1]);
  const std::string sampleName(argv[2]);
  const std::string topology(argv[3]);
  if (fileName.compare(fileName.size() - suffix.size(), fileName.size(), suffix) != 0)
    die("File must be a ROOT file (.root extension). Bail.");
  std::string sampleBase = fileName.substr(0, fileName.size() - suffix.size());
  ret = try_call("Run gOre Analysis", run_analysis, fileName, sampleBase, sampleName, topology, "gOre", "gOre", optimize_cuts);
  return ret;
}
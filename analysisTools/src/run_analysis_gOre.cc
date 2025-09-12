/**
 * @file run_analusis_gOre.cc
 * @brief Analyze the selection produced by gOre.toml
 * @details Make plots of the different variables of interest in terms of the categories
 * @author hhausner@fnal.gov
 **/

#include "TError.h"

#include <filesystem>

#include "include/analysis_tree.h"
#include "include/yell_try_die.h"

inline std::vector<std::pair<std::string, ana::tools::cut_sequence>> sel_cats =
{
  {"1#gammaXp (primary photon)", "true_category == 0"},
  {"1#gammaXp (decay photon)",   "true_category == 1"},
  {"1eXp",                       "true_category == 2"},
  {"Mulitple Showers",           "true_category == 3"},
  {"Muon Above Threshold",       "true_category == 4"},
  {"Pion Above Threshold",       "true_category == 5"},
  {"Other Topology Failure",     "true_category == 6"},
  {"Uncontained",                "true_category == 7"},
  {"Non-Fiducial",               "true_category == 8"},
  {"Cosmic",                     "true_category == 9"}
};

inline std::vector<std::pair<std::string, ana::tools::cut_sequence>> sig_cats =
{
  {"NC #Delta#rightarrowN#gamma", "true_mc_category == 0"},
  {"Other NC Single Photon",      "true_mc_category == 1"},
  {"NC #pi^{0}_{}",               "true_mc_category == 2"},
  {"Other NC",                    "true_mc_category == 3"},
  {"CC",                          "true_mc_category == 4"}
};

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

inline std::map<int, std::string> mc_cats =
{
  {0, "NC #Delta#rightarrowN#gamma"},
  {1, "NC Other post-FSI 1#gammaXp"},
  {2, "NC #pi^{0}_{}, post-FSI 0#gamma"},
  {3, "Other NC"},
  {4, "CC"},
  {5, "Cosmic"}
};

inline std::vector<std::tuple<std::string, bool, double, double>> vars_to_optimize =
{
  //{"reco_xy_wall_dist",            false,  0,   200},
  //{"reco_z_wall_dist",             false,  0,   200},
  {"reco_gOre_start_dedx",         false,  0,     5},
  {"reco_flash_total_pe",          false,  0, 20000},
  {"reco_gOre_directional_spread", true,   0,     1}
};

inline std::vector<std::tuple<std::string, std::string, bool, double, double>> vars_to_optimize_with_condition =
{
  {"reco_n_protons > 1",  "reco_gOre_score", false, -1, 1},
  {"reco_n_protons == 0", "reco_gOre_score", false, -1, 1}
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
                 const std::string& sample,
                 const std::string& sel,
                 const std::string& sig,
                 const bool& optimizeCuts)
{
  // setup output dir
  if (not std::filesystem::is_directory("plots") || not std::filesystem::exists("plots"))
    std::filesystem::create_directory("plots");
  if (not std::filesystem::is_directory("plots/"+sample) || not std::filesystem::exists("plots/"+sample))
    std::filesystem::create_directory("plots/"+sample);
  
  // setup analysis tree
  std::string fileLocation = std::filesystem::current_path().string();
  std::string directoryName = "events/simulation";
  std::string selTreeName = "selected";
  std::string sigTreeName = "signal";
  ana::tools::cut_sequence selection_cut = "reco_is_gOre";
  selection_cut += "108 < reco_gOre_ke";
  selection_cut += "reco_gOre_ke < 408";
  ana::tools::analysis_tree my_analysis_tree(fileLocation+"/"+fileName, directoryName,
                                             selection_cut,
                                             selTreeName, sigTreeName,
                                             sel_cats, sig_cats);
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~VARIABLE~~~~~~~~~~~~~~~~~~~BINS~~~~MIN~~~~~MAX~~~~~TITLE
  my_analysis_tree.add_variable("reco_flash_total_pe",         100,     0, 100000,    "Total Flash PE");
  my_analysis_tree.add_variable("reco_gOre_score",              50,    -1,      1,    "#gamma-candiate PID Score");
  my_analysis_tree.add_variable("reco_n_protons",                6,     0,      6,    "Protons above Threshold");
  my_analysis_tree.add_variable("true_vertex_x",                50,  -400,    400,    "True Vertex X (cm)");
  my_analysis_tree.add_variable("true_vertex_y",                50,  -200,    200,    "True Vertex Y (cm)");
  my_analysis_tree.add_variable("true_vertex_z",                50, -1000,   1000,    "True Vertex Z (cm)");
  my_analysis_tree.add_variable("reco_vertex_x",                50,  -400,    400,    "Vertex X (cm)");
  my_analysis_tree.add_variable("reco_vertex_y",                50,  -200,    200,    "Vertex Y (cm)");
  my_analysis_tree.add_variable("reco_vertex_z",                50, -1000,   1000,    "Vertex Z (cm)");
  my_analysis_tree.add_variable("true_xy_wall_dist",            50,     0,    200,    "True Minimum Distance from Vertex to X-/Y-side Detector Wall");
  my_analysis_tree.add_variable("true_z_wall_dist",             50,     0,   1000,    "True Minimum Distance from Vertex to Z-side Detector Wall");
  my_analysis_tree.add_variable("reco_xy_wall_dist",            50,     0,    200,    "Minimum Distance from Vertex to X-/Y-side Detector Wall");
  my_analysis_tree.add_variable("reco_z_wall_dist",             50,     0,   1000,    "Minimum Distance from Vertex to Z-side Detector Wall");
  my_analysis_tree.add_variable("reco_gOre_start_dedx",         50,     0,     10,    "#gamma-candiate Start dE/dx (MeV/cm)");
  my_analysis_tree.add_variable("reco_gOre_azimuthal_angle",    50,     0,      3.14, "#gamma-candiate Azimuthal Angle (rad)");
  my_analysis_tree.add_variable("reco_gOre_polar_angle",        50,     0,      3.14, "#gamma-candiate Polar Angle (rad)");
  my_analysis_tree.add_variable("true_gOre_ke",                100,     0,   1000,    "#gamma-candiate KE (MeV)");
  my_analysis_tree.add_variable("reco_gOre_ke",                100,     0,   1000,    "#gamma-candiate KE (MeV)");
  my_analysis_tree.add_variable("true_subleading_gOre_ke",      25,     0,     25,    "Highest True Subthreshold #gamma-candidate KE (MeV)");
  my_analysis_tree.add_variable("reco_subleading_gOre_ke",      25,     0,     25,    "Highest Subthreshold #gamma-candidate KE (MeV)");
  my_analysis_tree.add_variable("true_min_muon_ke",             25,     0,    200,    "Lowest Muon KE (MeV)");
  my_analysis_tree.add_variable("reco_min_muon_ke",             25,     0,    200,    "Lowest Muon KE (MeV)");
  my_analysis_tree.add_variable("true_min_pion_ke",             25,     0,     25,    "Lowest Pion KE (MeV)");
  my_analysis_tree.add_variable("reco_min_pion_ke",             25,     0,     25,    "Lowest Pion KE (MeV)");
  my_analysis_tree.add_variable("reco_gOre_straightness",       75,     0,      1,    "#gamma-candiate Straightness");
  my_analysis_tree.add_variable("reco_gOre_axial_spread",       75,     0,      1,    "#gamma-candiate Axial Spread");
  my_analysis_tree.add_variable("reco_gOre_directional_spread", 75,     0,      1,    "#gamma-candiate Directional Spread");
  my_analysis_tree.add_variable("reco_gOre_gap",               100,     0,    100,    "#gamma-candiate Distance from Vertex (cm)");
  my_analysis_tree.add_variable("reco_total_gOre_KE",          150,     0,    150,    "Total KE in Showers (incl. subthreshold)");
  my_analysis_tree.add_variable("reco_pion_mass",               50,     0,    200,    "Reconstructed Neutral Pion Mass Peak (MeV/c^{2}_{})");
  my_analysis_tree.add_variable("true_interaction_mode",        15,    -1.5,   13.5,  "GENIE Interaction Mode");
  my_analysis_tree.add_variable("true_interaction_type",       101,   999.5, 1100.5,  "GENIE Interaction Type");
  my_analysis_tree.add_variable("true_baryon_res_code",         19,    -1.5,   17.5,  "Resonance Number");
  my_analysis_tree.add_variable("true_mc_category",              6,    -0.5,   5.5,  "MC Truth Category"); 

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
  cut += "reco_is_gOre";
  try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
  //*** GORE KE ***/
  std::cout << "//** GORE KE ***//" << std::endl;
  cut += "108 < reco_gOre_ke";
  cut += "reco_gOre_ke < 408";
  try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
  //*** FLASH CUT ***//
  std::cout << "//*** FLASH CUT ***//" << std::endl;
  cut += "reco_flash";
  try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
  //*** FIDUCIAL CUT ***//
  std::cout << "//*** FIDUCIAL CUT ***//" << std::endl;
  //cut += "reco_fiducial";
  cut += "reco_xy_wall_dist > 50";
  cut += "reco_z_wall_dist > 50";
  try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
  //*** CONTAINMENT CUT***//
  std::cout << "//*** CONTAINMENT CUT ***//" << std::endl;
  cut += "reco_containment";
  try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
  ////*** PROTON VTX GAP ***//
  //std::cout << "//*** PROTON VTX GAP ***//" << std::endl;
  //cut.add_conditional_cut("reco_n_protons > 0", "reco_gOre_gap > 5");
  //try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); }); 

  if (optimizeCuts)
  {
    for (auto const& [var, upper, lw_end, up_end] : vars_to_optimize)
      cut = optimize_cut(my_analysis_tree, cut, var, upper, lw_end, up_end, sample, pdf_suffix);
    for (auto const& [condition, var, upper, lw_end, up_end] : vars_to_optimize_with_condition)
      cut = optimize_conditional_cut(my_analysis_tree, cut, condition, var, upper, lw_end, up_end, sample, pdf_suffix);
  }
  else
  {
    std::cout << "//*** OPTIMIZED CUTS ***//" << std::endl;
    cut += "reco_flash_total_pe > 2720";
    cut += "reco_gOre_directional_spread < 0.112";
    cut += "reco_gOre_score > -0.506";
    try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
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
    // some vars use alphanumeric labels
    if (var == "true_interaction_type" ||
        var == "true_interaction_mode" ||
        var == "true_baryon_res_code"  ||
        var == "true_mc_category"       )
    {
      std::function<std::string(double)> get_label = (var == "true_interaction_type") ? [](const int& code){ return (genie_inttypes.count(code) == 1) ? genie_inttypes.at(code) : ""; }
                                                   : (var == "true_interaction_mode") ? [](const int& code){ return genie_modes.at(code); }
                                                   : (var == "true_baryon_res_code")  ? [](const int& code){ return res_codes.at(code); }
                                                   : (var == "true_mc_category")      ? [](const int& code){ return mc_cats.at(code); }
                                                   :                                    [](const int& code){ return std::to_string(code); };
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
    var_plot.canvas->SaveAs(pdfName.c_str());
    var_plot_sig.canvas->SaveAs(pdfName_sig.c_str());
  }

  return 0;
}

int main(int argc, char* argv[])
{
  gErrorIgnoreLevel=3000;
  //gDebug=3;
  int ret = 0;
  bool optimize_cuts = true;
  if (not (argc > 1))
    die("Must pass in the name of the ROOT file to analyze. Bail.");
  const std::string suffix = ".root";
  const std::string fileName(argv[1]);
  if (fileName.compare(fileName.size() - suffix.size(), fileName.size(), suffix) != 0)
    die("File must be a ROOT file (.root extension). Bail.");
  std::string sample = fileName.substr(0, fileName.size() - suffix.size());
  ret = try_call("Run gOre Analysis", run_analysis, fileName, sample, "gOre", "gOre", optimize_cuts);
  return ret;
}
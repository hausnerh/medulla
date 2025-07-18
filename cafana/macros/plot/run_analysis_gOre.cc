#include "include/sig_bkg.h"
#include "include/yell_try_die.h"
#include <filesystem>

inline std::vector<std::pair<std::string, nc::cut_sequence>> sel_cats_delta =
{
  {"Fiducial/Contained NC #Delta#rightarrowN#gamma",     {"nc_delta_res_no_pion == 1", "category_gOre == 0"}},
  {"Non-Fiducial/Contained NC #Delta#rightarrowN#gamma", {"nc_delta_res_no_pion == 1", "category_gOre != 0"}},
  {"Other Fiducial/Contained True 1#gammaXp",            {"nc_delta_res_no_pion == 0", "category_gOre == 0"}},
  {"1eXp",                                               {"nc_delta_res_no_pion == 0", "category_gOre == 1"}},
  {"Res #pi^{0}_{}",                                     {"nc_delta_res_no_pion == 0", "category_gOre == 2"}},
  {"Res #pi^{+}_{}/#pi^{-}_{}",                          {"nc_delta_res_no_pion == 0", "category_gOre == 3"}},
  {"QE",                                                 {"nc_delta_res_no_pion == 0", "category_gOre == 4"}},
  {"DIS",                                                {"nc_delta_res_no_pion == 0", "category_gOre == 5"}},
  {"Non-fiducial or Uncontained #nu",                    {"nc_delta_res_no_pion == 0", "category_gOre == 6"}},
  {"Other #nu",                                          {"nc_delta_res_no_pion == 0", "category_gOre == 7"}},
  {"Cosmic",                                             {"nc_delta_res_no_pion == 0", "category_gOre == 8"}}
};

inline std::vector<std::pair<std::string, nc::cut_sequence>> sel_cats_gore =
{
  {"1#gammaXp",                       "category_gOre == 0"},
  {"1eXp",                            "category_gOre == 1"},
  {"Res #pi^{0}_{}",                  "category_gOre == 2"},
  {"Res #pi^{+}_{}/#pi^{-}_{}",       "category_gOre == 3"},
  {"QE",                              "category_gOre == 4"},
  {"DIS",                             "category_gOre == 5"},
  {"Non-fiducial or Uncontained #nu", "category_gOre == 6"},
  {"Other #nu",                       "category_gOre == 7"},
  {"Cosmic",                          "category_gOre == 8"}
};

inline std::vector<std::pair<std::string, nc::cut_sequence>> sig_cats =
{
  {"True 1#gammaXp Topology",   "gOre_is_photon == 1"}
  //{"True 1eXp Topology",        "gOre_is_electron == 1"}
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
  {1, "Other Single Photon"},
  {2, "#pi^{0} No Photons"},
  {3, "Other NC"},
  {4, "Other CC"},
  {5, "Non-Neutrino"}
};

int run_analysis(const std::string& sample,
                 const std::string& sel,
                 const std::string& sig,
                 const bool& deltaResSwitch,
                 const bool& optimizeCuts)
{
  // setup output dir
  if (not std::filesystem::is_directory("plots") || not std::filesystem::exists("plots"))
    std::filesystem::create_directory("plots");
  if (not std::filesystem::is_directory("plots/"+sample) || not std::filesystem::exists("plots/"+sample))
    std::filesystem::create_directory("plots/"+sample);
  
  // setup analysis tree
  std::string fileName = "nc_gOre_"+sample+".root";
  std::string directoryName = "events/nominal";
  std::string selTreeName = "Nu_Topology_gOre";
  std::string sigTreeName = (deltaResSwitch) ? "signalNu_FidCon_NCRes" : "signalNu_Topology_gOre";
  std::string cosTreeName = "Cos_Topology_gOre";
  std::vector<std::pair<std::string, nc::cut_sequence>> sel_cats = (deltaResSwitch) ? sel_cats_delta
                                                                               : sel_cats_gore;
  nc::analysis_tree my_analysis_tree(fileName, directoryName, "(is_gOre == 1)",
                                     selTreeName, sigTreeName, cosTreeName,
                                     sel_cats, sig_cats);
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~VARIABLE~~~~~~~~~~~~~~BINS~~~~MIN~~~~~MAX~~~~~TITLE
  my_analysis_tree.add_variable("flash_total_PE",         100,     0, 100000,    "Total Flash PE");
  my_analysis_tree.add_variable("reco_edep",               15,     0,      1.5,  "Reconstructed E^{ }_{Dep} (GeV)");
  my_analysis_tree.add_variable("min_muon_ke",            150,     0,    150,    "Minimum Subthreshold Muon KE (MeV)");
  my_analysis_tree.add_variable("min_pion_ke",            250,     0,     25,    "Minimum Subthreshold Pion KE (MeV)");
  my_analysis_tree.add_variable("gOre_score",              50,    -1,      1,    "#gamma-candiate PID Score");
  my_analysis_tree.add_variable("n_protons",                6,     0,      6,    "Protons above Threshold");
  my_analysis_tree.add_variable("true_vertex_x",           50,  -400,    400,    "True Vertex X (cm)");
  my_analysis_tree.add_variable("true_vertex_y",           50,  -200,    200,    "True Vertex Y (cm)");
  my_analysis_tree.add_variable("true_vertex_z",           50, -1000,   1000,    "True Vertex Z (cm)");
  my_analysis_tree.add_variable("vertex_x",                50,  -400,    400,    "Vertex X (cm)");
  my_analysis_tree.add_variable("vertex_y",                50,  -200,    200,    "Vertex Y (cm)");
  my_analysis_tree.add_variable("vertex_z",                50, -1000,   1000,    "Vertex Z (cm)");
  my_analysis_tree.add_variable("true_wall_xy",            50,     0,    200,    "True Minimum Distance from Vertex to X-/Y-side Detector Wall");
  my_analysis_tree.add_variable("true_wall_z",             50,     0,   1000,    "True Minimum Distance from Vertex to Z-side Detector Wall");
  my_analysis_tree.add_variable("wall_xy",                 50,     0,    200,    "Minimum Distance from Vertex to X-/Y-side Detector Wall");
  my_analysis_tree.add_variable("wall_z",                  50,     0,   1000,    "Minimum Distance from Vertex to Z-side Detector Wall");
  my_analysis_tree.add_variable("gOre_start_dedx",         50,     0,     10,    "#gamma-candiate Start dE/dx (MeV/cm)");
  my_analysis_tree.add_variable("gOre_azimuthal_angle",    50,     0,      3.14, "#gamma-candiate Azimuthal Angle (rad)");
  my_analysis_tree.add_variable("gOre_polar_angle",        50,     0,      3.14, "#gamma-candiate Polar Angle (rad)");
  my_analysis_tree.add_variable("gOre_ke",                 20,     0,   2000,    "#gamma-candiate KE (MeV)");
  my_analysis_tree.add_variable("subleading_gOre_ke",      25,     0,     25,    "Highest Subthreshold #gamma-candidate KE (MeV)");
  my_analysis_tree.add_variable("true_subleading_gOre_ke", 25,     0,     25,    "Highest True Subthreshold #gamma-candidate KE (MeV)");
  my_analysis_tree.add_variable("gOre_straightness",       75,     0,      1,    "#gamma-candiate Straightness");
  my_analysis_tree.add_variable("gOre_axial_spread",       75,     0,      1,    "#gamma-candiate Axial Spread");
  my_analysis_tree.add_variable("gOre_directional_spread", 75,     0,      1,    "#gamma-candiate Directional Spread");
  my_analysis_tree.add_variable("gOre_gap",               100,     0,    100,    "#gamma-candiate Distance from Vertex (cm)");
  my_analysis_tree.add_variable("total_gOre_ke",          150,     0,    150,    "Total KE in Showers (incl. subthreshold)");
  my_analysis_tree.add_variable("pion_mass",               50,     0,    200,    "Reconstructed Neutral Pion Mass Peak (MeV/c^{2}_{}");
  my_analysis_tree.add_variable("interaction_mode",        15,    -1.5,   13.5,  "GENIE Interaction Mode");
  my_analysis_tree.add_variable("interaction_type",       101,   999.5, 1100.5,  "GENIE Interaction Type");
  my_analysis_tree.add_variable("res_code",                19,    -1.5,   17.5,  "Resonance Number");
  my_analysis_tree.add_variable("gOre_mc_category",         6,    -0.5,    5.5,  "MC Category");

  std::string signal_def = (deltaResSwitch) ? "NC Delta Res No Pions " + sel : "TOPOLOGICAL "+sel;
  std::string pdf_suffix = "_"+sample+"_sig-"+sig+"_sel-"+sel;
  pdf_suffix += (deltaResSwitch) ? "_NCDeltaRes.pdf" : ".pdf";
  nc::cut_sequence cut;
  std::cout
    << "************************\n"
    << "* " << signal_def << " *\n"
    << "************************" << std::endl;
  //*** TOPOLOGY ***//
  // should alread be implemented as part of making the selected TTree,
  // but we do this for completeness
  std::cout << "//*** TOPOLOGY ***//" << std::endl;
  cut += "is_gOre";
  try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
  //*** FLASH CUT ***//
  std::cout << "//*** FLASH CUT ***//" << std::endl;
  cut += "flash_cut";
  try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
  //*** FIDUCIAL CUT ***//
  std::cout << "//*** FIDUCIAL CUT ***//" << std::endl;
  cut += "fiducial_cut";
  try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
  //*** CONTAINMENT CUT***//
  std::cout << "//*** CONTAINMENT CUT ***//" << std::endl;
  cut += "containment_cut";
  try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });

  if (optimizeCuts)
  {
    //*** OPTIMIZE WALL XY ***//
    std::cout << "//*** OPTIMIZE WALL XY ***//" << std::endl;
    auto [wall_xy_cut_opt, wall_xy_cut_plot]
      = try_call("optimize fidutial cut (XY)",
          [&my_analysis_tree, &cut]{ return my_analysis_tree.optimize_lower_bound("wall_xy", 0, 200, cut); });
    cut = wall_xy_cut_opt.cut;
    try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
    wall_xy_cut_plot.sc.canvas->Print(("plots/"+sample+"/wall_xy_optimized"+pdf_suffix).c_str());
    wall_xy_cut_opt.canvas    ->Print(("plots/"+sample+"/wall_xy_FOM"+pdf_suffix).c_str());

    //*** OPTIMIZE WALL Z ***//
    std::cout << "//*** OPTIMIZE WALL Z ***//" << std::endl;
    auto [wall_z_cut_opt, wall_z_cut_plot]
      = try_call("optimize fidutial cut (Z)",
          [&my_analysis_tree, &cut]{ return my_analysis_tree.optimize_lower_bound("wall_z", 0, 200, cut); });
    cut = wall_z_cut_opt.cut;
    try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
    wall_z_cut_plot.sc.canvas->Print(("plots/"+sample+"/wall_z_optimized"+pdf_suffix).c_str());
    wall_z_cut_opt.canvas    ->Print(("plots/"+sample+"/wall_z_FOM"+pdf_suffix).c_str());

    //*** OPTIMIZE PE CUT ***//
    std::cout << "//*** OPTIMIZE PE CUT ***//" << std::endl;
    auto [pe_cut_opt, pe_cut_plot]
      = try_call("optimize photon/electron score",
          [&my_analysis_tree, &cut]{ return my_analysis_tree.optimize_lower_bound("flash_total_PE", 0, 20000, cut); });
    cut = pe_cut_opt.cut;
    try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
    pe_cut_plot.sc.canvas->Print(("plots/"+sample+"/gOre_pe_optimized"+pdf_suffix).c_str());
    pe_cut_opt.canvas    ->Print(("plots/"+sample+"/gOre_pe_FOM"+pdf_suffix).c_str());

    //*** OPTIMIZE DIRECTIONAL SPREAD CUT  ***//
    std::cout << "//*** OPTIMIZE DIRECTIONAL SPREAD CUT ***//" << std::endl;
    auto [dir_cut_opt, dir_cut_plot]
      = try_call("optimize straightness",
          [&my_analysis_tree, &cut]{ return my_analysis_tree.optimize_upper_bound("gOre_directional_spread", 0, 1, cut); });
    cut = dir_cut_opt.cut;
    try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
    dir_cut_plot.sc.canvas->Print(("plots/"+sample+"/gOre_directional_spread_optimized"+pdf_suffix).c_str());
    dir_cut_opt.canvas    ->Print(("plots/"+sample+"/gOre_directional_spread_FOM"+pdf_suffix).c_str());

    //*** OPTIMIZE PION MASS CUT ***//
    //std::cout << "//*** OPTIMIZE PION MASS CUT ***//" << std::endl;
    auto [pi0_cut_opt, pi0_cut_plot]
      = try_call("optimize pion mass",
          [&my_analysis_tree, &cut]{ return my_analysis_tree.optimize_lower_bound("pion_mass", 0, 900, cut); });
    cut = pi0_cut_opt.cut;
    try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
    pi0_cut_plot.sc.canvas->Print(("plots/"+sample+"/pion_mass_optimized"+pdf_suffix).c_str());
    pi0_cut_opt.canvas    ->Print(("plots/"+sample+"/pion_mass_FOM"+pdf_suffix).c_str());

    //*** OPTIMIZE PID CUT ***//
    std::cout << "//*** OPTIMIZE PID CUT ***//" << std::endl;
    auto [pid_cut_opt, pid_cut_plot]
      = try_call("optimize photon/electron score",
          [&my_analysis_tree, &cut]{ return my_analysis_tree.optimize_lower_bound("gOre_score", -1, 1, cut); });
    cut = pid_cut_opt.cut;
    try_call(cut.string(), [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
    pid_cut_plot.sc.canvas->Print(("plots/"+sample+"/gOre_pid_optimized"+pdf_suffix).c_str());
    pid_cut_opt.canvas    ->Print(("plots/"+sample+"/gOre_pid_FOM"+pdf_suffix).c_str());
  }
  else
  {
    std::cout << "//*** OPTIMIZED CUTS ***//" << std::endl;
    cut += "flash_total_PE > 2720";
    cut += "gOre_directional_spread < 0.112";
    cut += "gOre_score > -0.506";
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
    if (var == "interaction_type" ||
        var == "interaction_mode" ||
        var == "res_code"         ||
        var == "gOre_mc_category"  )
    {
      std::function<std::string(double)> get_label = (var == "interaction_type") ? [](const int& code){ return (genie_inttypes.count(code) == 1) ? genie_inttypes.at(code) : ""; }
                                                   : (var == "interaction_mode") ? [](const int& code){ return genie_modes.at(code); }
                                                   : (var == "res_code")         ? [](const int& code){ return res_codes.at(code); }
                                                   : (var == "gOre_mc_category") ? [](const int& code){ return mc_cats.at(code); }
                                                   :                               [](const int& code){ return std::to_string(code); };
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
      std::string title_str = (var == "interaction_type") ? "GENIE Interaction Type;;Events"
                            : (var == "interaction_mode") ? "GENIE Interaction Mode;;Events"
                            : (var == "res_code")         ?              "Resonance;;Events"
                            : (var == "gOre_mc_category") ?        "Post-FSI State;;Evenets"
                            :                                                     ";;Events";
      var_plot.stack->SetTitle(title_str.c_str());
      var_plot_sig.stack->SetTitle(title_str.c_str());
      if (var == "interaction_type")
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

  // specifically check the pion mass in the 1 shower case and the multi-shower case
  nc::cut_sequence cut_1_shower = cut.with_addition("n_gOre_showers == 1");
  nc::cut_sequence cut_n_shower = cut.with_addition("n_gOre_showers > 1");
  auto plot_1_shower =
    try_call("plot pion mass for 1 shower",
      [&my_analysis_tree, &cut_1_shower]{ return my_analysis_tree.plot_var_sel("pion_mass", cut_1_shower); }); 
  auto plot_n_shower =
    try_call("plot pion mass for n>1 showers",
      [&my_analysis_tree, &cut_n_shower]{ return my_analysis_tree.plot_var_sel("pion_mass", cut_n_shower); });
  plot_1_shower.stack->SetTitle("1 Shower");
  plot_n_shower.stack->SetTitle("N>1 Showers");
  std::string pdfName_1_shower = "plots/"+sample+"/pion_mass_1_shower"+pdf_suffix;
  std::string pdfName_n_shower = "plots/"+sample+"/pion_mass_n_shower"+pdf_suffix;
  plot_1_shower.canvas->SaveAs(pdfName_1_shower.c_str());
  plot_n_shower.canvas->SaveAs(pdfName_n_shower.c_str());

  return 0;
}

int main(int argc, char* argv[])
{
  gErrorIgnoreLevel=3000;
  //gDebug=3;
  int ret = 0;
  if (argc == 0)
    die("must supply a sample to select which file is used.");
  std::string sample (argv[1]);
  bool optimize_cuts = false;
  ret = try_call("Run gOre Analysis", run_analysis, sample, "gOre", "gOre", false, optimize_cuts);
  return ret;
}

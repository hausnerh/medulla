#include "include/sig_bkg.h"
#include "include/yell_try_die.h"
#include <filesystem>

inline std::vector<std::pair<std::string, std::string>> sel_cats_delta =
{
  {"Fiducial/Contained NC #Delta#rightarrowN#gamma",     "(nc_delta_res_no_pion == 1) && (category_gOre == 0)"},
  {"Non-Fiducial/Contained NC #Delta#rightarrowN#gamma", "(nc_delta_res_no_pion == 1) && (category_gOre != 0)"},
  {"Other True 1#gammaXp",                               "(nc_delta_res_no_pion == 0) && (category_gOre == 0)"},
  {"1eXp",                                               "(nc_delta_res_no_pion == 0) && (category_gOre == 1)"},
  {"Res #pi^{0}_{}",                                     "(nc_delta_res_no_pion == 0) && (category_gOre == 2)"},
  {"Res #pi^{+}_{}/#pi^{-}_{}",                          "(nc_delta_res_no_pion == 0) && (category_gOre == 3)"},
  {"QE",                                                 "(nc_delta_res_no_pion == 0) && (category_gOre == 4)"},
  {"DIS",                                                "(nc_delta_res_no_pion == 0) && (category_gOre == 5)"},
  {"Other #nu",                                          "(nc_delta_res_no_pion == 0) && (category_gOre == 6)"},
  {"Cosmic",                                             "(nc_delta_res_no_pion == 0) && (category_gOre == 7)"}
};

inline std::vector<std::pair<std::string, std::string>> sel_cats_gore =
{
  {"1#gammaXp",                 "(category_gOre == 0)"},
  {"1eXp",                      "(category_gOre == 1)"},
  {"Res #pi^{0}_{}",            "(category_gOre == 2)"},
  {"Res #pi^{+}_{}/#pi^{-}_{}", "(category_gOre == 3)"},
  {"QE",                        "(category_gOre == 4)"},
  {"DIS",                       "(category_gOre == 5)"},
  {"Other #nu",                 "(category_gOre == 6)"},
  {"Cosmic",                    "(category_gOre == 7)"}
};

inline std::vector<std::pair<std::string, std::string>> sig_cats =
{
  {"True 1#gammaXp Topology",   "(gOre_is_photon == 1)"}
  //{"True 1eXp Topology",        "(gOre_is_electron == 1)"}
};

int run_analysis(const std::string& treeDir,
                 const std::string& sample, 
                 const std::string& sel,
                 const std::string& sig,
                 const bool& deltaResSwitch,
                 const bool& optimizeCuts)
{
  // setup output dir
  if (not std::filesystem::is_directory("plots") || not std::filesystem::exists("plots"))
    std::filesystem::create_directory("plots");
  if (not std::filesystem::is_directory("plots/"+treeDir) || not std::filesystem::exists("plots/"+treeDir))
    std::filesystem::create_directory("plots/"+treeDir);
  if (not std::filesystem::is_directory("plots/"+treeDir+"/"+sample) || not std::filesystem::exists("plots/"+treeDir+"/"+sample))
    std::filesystem::create_directory("plots/"+treeDir+"/"+sample);
  
  // setup analysis tree
  //std::string fileName = "trees/"+treeDir+"/nc_gOre_"+sample+".root";
  std::string fileName = "nc_gOre_"+sample+".root";
  std::string directoryName = "events/nominal";
  std::string selTreeName = "Nu_Topology_gOre";
  std::string sigTreeName = (deltaResSwitch) ? "signalNu_FidCon_NCRes" : "signalNu_Topology_gOre";
  std::string cosTreeName = "Cos_Topology_gOre";
  std::vector<std::pair<std::string, std::string>> sel_cats = (deltaResSwitch) ? sel_cats_delta
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
  my_analysis_tree.add_variable("interaction_mode",        16,    -1.5,   14.5,  "GENIE Interaction Mode");
  my_analysis_tree.add_variable("interaction_type",       100,   999.5, 1100.5,  "GENIE Interaction Type");
  my_analysis_tree.add_variable("res_code",                18,    -1.5,   17.5,  "Resonance Number");
  my_analysis_tree.add_variable("gOre_mc_category",         6,    -0.5,    5.5,  "MC Category");

  std::string signal_def = (deltaResSwitch) ? "NC Delta Res No Pions " + sel : "TOPOLOGICAL "+sel;
  std::string pdf_suffix = "_"+sample+"_sig-"+sig+"_sel-"+sel;
  pdf_suffix += (deltaResSwitch) ? "_NCDeltaRes.pdf" : ".pdf";
  std::cout
    << "************************\n"
    << "* " << signal_def << " *\n"
    << "************************" << std::endl;
  std::string cut = "";
  //*** TOPOLOGY ***//
  // should alread be implemented as part of making the selected TTree,
  // but we do this for completeness
  std::cout << "//*** TOPOLOGY ***//" << std::endl;
  cut += "(is_gOre)";
  try_call(cut, [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
  //*** FLASH CUT ***//
  std::cout << "//*** FLASH CUT ***//" << std::endl;
  cut += " && (flash_cut)";
  try_call(cut, [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
  //*** FIDUCIAL CUT ***//
  std::cout << "//*** FIDUCIAL CUT ***//" << std::endl;
  cut += " && (fiducial_cut)";
  try_call(cut, [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
  //*** CONTAINMENT CUT***//
  std::cout << "//*** CONTAINMENT CUT ***//" << std::endl;
  cut += " && (containment_cut)";
  try_call(cut, [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });

  if (optimizeCuts)
  {
    //*** OPTIMIZE WALL XY ***//
    std::cout << "//*** OPTIMIZE WALL XY ***//" << std::endl;
    auto [wall_xy_cut_opt, wall_xy_cut_plot]
      = try_call("optimize fidutial cut (XY)",
          [&my_analysis_tree, &cut]{ return my_analysis_tree.optimize_lower_bound("wall_xy", 0, 200, cut); });
    cut = wall_xy_cut_opt.cut;
    try_call(cut, [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
    wall_xy_cut_plot.sc.canvas->Print(("plots/"+treeDir+"/"+sample+"/wall_xy_optimized"+pdf_suffix).c_str());
    wall_xy_cut_opt.canvas    ->Print(("plots/"+treeDir+"/"+sample+"/wall_xy_FOM"+pdf_suffix).c_str());

    //*** OPTIMIZE WALL Z ***//
    std::cout << "//*** OPTIMIZE WALL Z ***//" << std::endl;
    auto [wall_z_cut_opt, wall_z_cut_plot]
      = try_call("optimize fidutial cut (Z)",
          [&my_analysis_tree, &cut]{ return my_analysis_tree.optimize_lower_bound("wall_z", 0, 200, cut); });
    cut = wall_z_cut_opt.cut;
    try_call(cut, [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
    wall_z_cut_plot.sc.canvas->Print(("plots/"+treeDir+"/"+sample+"/wall_z_optimized"+pdf_suffix).c_str());
    wall_z_cut_opt.canvas    ->Print(("plots/"+treeDir+"/"+sample+"/wall_z_FOM"+pdf_suffix).c_str());

    //*** OPTIMIZE PE CUT ***//
    std::cout << "//*** OPTIMIZE PE CUT ***//" << std::endl;
    auto [pe_cut_opt, pe_cut_plot]
      = try_call("optimize photon/electron score",
          [&my_analysis_tree, &cut]{ return my_analysis_tree.optimize_lower_bound("flash_total_PE", 0, 20000, cut); });
    cut = pe_cut_opt.cut;
    try_call(cut, [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
    pe_cut_plot.sc.canvas->Print(("plots/"+treeDir+"/"+sample+"/gOre_pe_optimized"+pdf_suffix).c_str());
    pe_cut_opt.canvas    ->Print(("plots/"+treeDir+"/"+sample+"/gOre_pe_FOM"+pdf_suffix).c_str());

    //*** OPTIMIZE DIRECTIONAL SPREAD CUT  ***//
    std::cout << "//*** OPTIMIZE DIRECTIONAL SPREAD CUT ***//" << std::endl;
    auto [dir_cut_opt, dir_cut_plot]
      = try_call("optimize straightness",
          [&my_analysis_tree, &cut]{ return my_analysis_tree.optimize_upper_bound("gOre_directional_spread", 0, 1, cut); });
    cut = dir_cut_opt.cut;
    try_call(cut, [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
    dir_cut_plot.sc.canvas->Print(("plots/"+treeDir+"/"+sample+"/gOre_directional_spread_optimized"+pdf_suffix).c_str());
    dir_cut_opt.canvas    ->Print(("plots/"+treeDir+"/"+sample+"/gOre_directional_spread_FOM"+pdf_suffix).c_str());

    //*** OPTIMIZE PION MASS CUT ***//
    //std::cout << "//*** OPTIMIZE PION MASS CUT ***//" << std::endl;
    auto [pi0_cut_opt, pi0_cut_plot]
      = try_call("optimize pion mass",
          [&my_analysis_tree, &cut]{ return my_analysis_tree.optimize_lower_bound("pion_mass", 0, 900, cut); });
    cut = pi0_cut_opt.cut;
    try_call(cut, [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
    pi0_cut_plot.sc.canvas->Print(("plots/"+treeDir+"/"+sample+"/pion_mass_optimized"+pdf_suffix).c_str());
    pi0_cut_opt.canvas    ->Print(("plots/"+treeDir+"/"+sample+"/pion_mass_FOM"+pdf_suffix).c_str());

    //*** OPTIMIZE PID CUT ***//
    std::cout << "//*** OPTIMIZE PID CUT ***//" << std::endl;
    auto [pid_cut_opt, pid_cut_plot]
      = try_call("optimize photon/electron score",
          [&my_analysis_tree, &cut]{ return my_analysis_tree.optimize_lower_bound("gOre_score", -1, 1, cut); });
    cut = pid_cut_opt.cut;
    try_call(cut, [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
    pid_cut_plot.sc.canvas->Print(("plots/"+treeDir+"/"+sample+"/gOre_pid_optimized"+pdf_suffix).c_str());
    pid_cut_opt.canvas    ->Print(("plots/"+treeDir+"/"+sample+"/gOre_pid_FOM"+pdf_suffix).c_str());
  }
  else
  {
    std::cout << "//*** PID CUT ***//" << std::endl;
    cut += " && (gOre_score > -0.388)";
    try_call(cut, [&my_analysis_tree, &cut]{ my_analysis_tree.report_on_cut(cut); });
  }

  //*** Plot Vars ***//
  for (auto const& var : my_analysis_tree.variables())
  {
    auto var_plot =
      try_call("plot "+var,
        [&my_analysis_tree, &var, &cut]{ return my_analysis_tree.plot_var_sel(var, cut); });
    std::string pdfName = "plots/"+treeDir+"/"+sample+"/"+var+pdf_suffix;
    var_plot.canvas->SaveAs(pdfName.c_str());
    auto var_plot_sig =
     try_call("plot "+var+" signal",
       [&my_analysis_tree, &var, &cut]{ return my_analysis_tree.plot_var_sig(var, cut); });
    std::string pdfName_sig = "plots/"+treeDir+"/"+sample+"/signal_"+var+pdf_suffix;
    var_plot_sig.canvas->SaveAs(pdfName_sig.c_str());
  }

  // specifically check the pion mass in the 1 shower case and the multi-shower case
  std::string cut_1_shower = cut + " && (n_gOre_showers == 1)";
  std::string cut_n_shower = cut + " && (n_gOre_showers > 1)";
  auto plot_1_shower =
    try_call("plot pion mass for 1 shower",
      [&my_analysis_tree, &cut_1_shower]{ return my_analysis_tree.plot_var_sel("pion_mass", cut_1_shower); }); 
  auto plot_n_shower =
    try_call("plot pion mass for n>1 showers",
      [&my_analysis_tree, &cut_n_shower]{ return my_analysis_tree.plot_var_sel("pion_mass", cut_n_shower); });
  plot_1_shower.stack->SetTitle("1 Shower");
  plot_n_shower.stack->SetTitle("N>1 Showers");
  std::string pdfName_1_shower = "plots/"+treeDir+"/"+sample+"/pion_mass_1_shower"+pdf_suffix;
  std::string pdfName_n_shower = "plots/"+treeDir+"/"+sample+"/pion_mass_n_shower"+pdf_suffix;
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
  std::string treeDir(argv[1]);
  std::string sample (argv[2]);
  bool optimize_cuts = true;
  ret += try_call("Run gOre Analysis", run_analysis, treeDir, sample, "gOre", "gOre", false, optimize_cuts);
  ret += try_call("Run gOre Analysis", run_analysis, treeDir, sample, "gOre", "gOre", true,  optimize_cuts);
  return ret;
}

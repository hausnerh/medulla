/**
 * @file run_analysis_gOre.cc
 * @brief Analyze the 1g1p selection produced by gOre.toml
 * @details Plots selection variables at each stage of the cut sequence:
 *
 *   Stage 1 — Preselection
 *     reco_gOre_topology, reco_n_protons, reco_flash, reco_fiducial, reco_containment
 *
 *   Stage 2 — pi0 Rejection                                  [KEY: subleading_ke]
 *     reco_subleading_primary_gOre_shower_ke   <-- primary pi0 handle
 *     reco_leading_primary_gOre_primary_softmax
 *     reco_leading_primary_proton_primary_softmax
 *     abs(reco_delta_mass_P - 1232)
 *
 *   Stage 3 — e/gamma Separation                             [KEY: start_dedx, photon_softmax]
 *     reco_leading_primary_gOre_start_dedx     <-- single most powerful e/gamma discriminant
 *     reco_leading_primary_gOre_photon_softmax <-- NN score trained directly for gamma/e sep
 *     reco_leading_primary_gOre_axial_spread
 *     reco_leading_primary_gOre_directional_spread
 *
 * @author hhausner@fnal.gov
 **/

#include "TError.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include <filesystem>

#include "include/analysis_tree.h"
#include "include/yell_try_die.h"

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------
inline const double pion_mass = 134.9768; // MeV/c^2

// ---------------------------------------------------------------------------
// Fit: Landau (combinatorial background) + Gaussian (pi0 peak)
// ---------------------------------------------------------------------------
double landauPlusGauss(double* x, double* par)
{
  return par[0] * TMath::Landau(x[0], par[1], par[2], true)
       + par[3] * TMath::Gaus(x[0], par[4], par[5], true);
}

// ---------------------------------------------------------------------------
// Topological truth categories (used to label axes in the variable loop)
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// GENIE labels (used to label axes in the variable loop)
// ---------------------------------------------------------------------------
inline std::map<int, std::string> genie_modes =
{
  { -1, "Unknown"},
  {  0, "Quasi-Elastic"},
  {  1, "Resonance"},
  {  2, "DIS"},
  {  3, "Coherent"},
  {  4, "Coherent Elastic"},
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

// ---------------------------------------------------------------------------
// Variables and scan ranges used when optimizing cut thresholds
// ---------------------------------------------------------------------------
inline std::vector<std::tuple<std::string, bool, double, double>> vars_to_optimize =
{
  {"reco_leading_primary_gOre_primary_softmax",    false, 0.9,  1  },
  {"reco_leading_primary_proton_primary_softmax",  false, 0.9,  1  },
  {"abs(reco_delta_mass_P - 1232)",                true,  0,    300},
  {"reco_leading_primary_gOre_start_dedx",         false, 0,    5  },
  {"reco_leading_primary_gOre_axial_spread",       false, -0.5, 1  },
  {"reco_leading_primary_gOre_directional_spread", true,  0,    1  },
  {"reco_leading_primary_gOre_photon_softmax",     false, 0.9,  1  }
};

inline std::vector<std::tuple<std::string, std::string, bool, double, double>> vars_to_optimize_with_condition =
{
  {"reco_n_protons >= 1", "reco_gOre_gap", false, 0, 100}
};

// ---------------------------------------------------------------------------
// Helper functions: threshold/bound optimization
// ---------------------------------------------------------------------------
std::tuple<ana::tools::cut_sequence, ana::tools::cut_sequence>
optimize_threshold(ana::tools::analysis_tree& the_analysis_tree,
                   const ana::tools::cut_sequence& old_cut,
                   const ana::tools::cut_sequence& old_truth_cut,
                   const std::string& var, const std::string& truth_var,
                   const double& lw_end, const double up_end,
                   const std::string& sample, const std::string& pdf_suffix)
{
  const std::string report = "Optimize " + var + " Threshold";
  std::cout << report << std::endl;
  auto [cut_opt, cut_plot] =
    try_call(report,
      [&]{ return the_analysis_tree.optimize_threshold(var, truth_var, lw_end, up_end,
                                                       old_cut, old_truth_cut); });
  ana::tools::cut_sequence cut       = cut_opt.cut;
  ana::tools::cut_sequence truth_cut = cut_opt.truth_cut;
  try_call(cut.string(), [&]{ the_analysis_tree.report_on_cut(cut); });
  cut_plot.sc.PrintWIP("plots/"+sample+"/"+var+"_optimized"+pdf_suffix);
  cut_opt.PrintWIP("plots/"+sample+"/"+var+"_FOM"+pdf_suffix);
  return std::tie(cut, truth_cut);
}

ana::tools::cut_sequence
optimize_cut(const ana::tools::analysis_tree& the_analysis_tree,
             const ana::tools::cut_sequence& old_cut,
             const std::string& var, const bool& upper,
             const double& lw_end, const double up_end,
             const std::string& sample, const std::string& pdf_suffix)
{
  const std::string report = (upper) ? "Optimize Upper Bound on " + var
                                     : "Optimize Lower Bound on " + var;
  std::cout << report << std::endl;
  auto [cut_opt, cut_plot] =
    try_call(report,
      [&]{ return the_analysis_tree.optimize_bound(var, lw_end, up_end, old_cut, upper); });
  ana::tools::cut_sequence cut = cut_opt.cut;
  try_call(cut.string(), [&]{ the_analysis_tree.report_on_cut(cut); });
  cut_plot.sc.PrintWIP("plots/"+sample+"/"+var+"_optimized"+pdf_suffix);
  cut_opt.PrintWIP("plots/"+sample+"/"+var+"_FOM"+pdf_suffix);
  return cut;
}

ana::tools::cut_sequence
optimize_conditional_cut(const ana::tools::analysis_tree& the_analysis_tree,
                         const ana::tools::cut_sequence& old_cut,
                         const std::string& condition, const std::string& var,
                         const bool& upper, const double& lw_end, const double up_end,
                         const std::string& sample, const std::string& pdf_suffix)
{
  const std::string report = (upper) ? "Optimize Upper Bound on " + var + " Under Condition " + condition
                                     : "Optimize Lower Bound on " + var + " Under Condition " + condition;
  std::cout << report << std::endl;
  auto [cut_opt, cut_plot] =
    try_call(report,
      [&]{ return the_analysis_tree.optimize_conditional_bound(condition, var, lw_end, up_end,
                                                               old_cut, upper); });
  ana::tools::cut_sequence cut = cut_opt.cut;
  try_call(cut.string(), [&]{ the_analysis_tree.report_on_cut(cut); });
  cut_plot.sc.PrintWIP("plots/"+sample+"/"+var+"_optimized_with_"+condition+pdf_suffix);
  cut_opt.PrintWIP("plots/"+sample+"/"+var+"_FOM_with_"+condition+pdf_suffix);
  return cut;
}

// ===========================================================================
// run_analysis
// ===========================================================================
int run_analysis(const std::string& fileName,
                 const std::string& sampleBase,
                 const std::string& sampleName,
                 const std::string& topology,
                 const std::string& sel,
                 const std::string& sig,
                 const bool& optimizeCuts)
{
  // -------------------------------------------------------------------------
  // 1g1p final cut values
  //   Edit these to change the selection without hunting through the logic.
  // -------------------------------------------------------------------------
  // Stage 2 — pi0 rejection
  const double cut_subleading_ke          = 0.0;    // veto: !(ke > 0)
  const double cut_primary_softmax        = 0.9991; // shower GNN primary score
  const double cut_proton_primary_softmax = 0.9513; // proton GNN primary score
  const double cut_delta_mass_window      = 59.7;   // |M_reco - 1232| < X MeV/c^2
  // Stage 3 — e/gamma separation
  const double cut_start_dedx             = 4.07;   // [KEY] gamma ~4-5 MeV/cm, e ~2 MeV/cm
  const double cut_axial_spread_lo        = -0.5;   // NOTE: equals histogram lower edge; effectively no-op
  const double cut_directional_spread_hi  = 0.037;  // tight; shower collimation
  const double cut_photon_softmax         = 0.9225; // [KEY] NN gamma/e score

  // -------------------------------------------------------------------------
  // Output directory setup
  // -------------------------------------------------------------------------
  const std::string sample = sampleBase + "/" + sampleName + "/" + topology;
  for (const auto& dir : {std::string("plots"),
                          std::string("plots/"+sampleBase),
                          std::string("plots/"+sampleBase+"/"+sampleName),
                          std::string("plots/"+sample)})
    if (not std::filesystem::exists(dir))
      std::filesystem::create_directory(dir);

  // -------------------------------------------------------------------------
  // Signal/MC truth categories (topology-dependent)
  // -------------------------------------------------------------------------
  const std::vector<std::pair<std::string, ana::tools::cut_sequence>> sig_cats =
    (topology == "1g0p") ?
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
    }) : (topology == "1g1p") ?
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
    }) :
    std::vector<std::pair<std::string, ana::tools::cut_sequence>>({
      {"NC #Delta#rightarrowN#gamma (1#gammaNp)",      "(true_mc_category == 0) && (true_n_protons >= 1)"},
      {"NC #Delta#rightarrowN#gamma (1#gamma0p)",      "(true_mc_category == 0) && (true_n_protons == 0)"},
      {"Other NC 1#gammaXp Post-FSI",                  "true_mc_category == 1"},
      {"NC #pi^{0}_{} (#Delta Res)",                   "true_mc_category == 2"},
      {"NC #pi^{0}_{} (Other)",                        "true_mc_category == 3"},
      {"NC #pi^{+/-}_{}",                              "true_mc_category == 4"},
      {"Other NC",                                     "true_mc_category == 5"},
      {"CC e",                                         "true_mc_category == 6"},
      {"Other CC",                                     "true_mc_category == 7"}
    });

  const std::map<int, std::string> mc_cats = []
    (const std::vector<std::pair<std::string, ana::tools::cut_sequence>>& cats)
  {
    std::map<int, std::string> m;
    for (int i = 0; i < static_cast<int>(cats.size()); ++i) m[i] = cats[i].first;
    return m;
  }(sig_cats);

  // -------------------------------------------------------------------------
  // Analysis tree: file, directory, initial selection/signal cuts
  // -------------------------------------------------------------------------
  const ana::tools::cut_sequence reco_topology_cut = (topology == "1g0p") ? "reco_n_protons == 0"
                                                    : (topology == "1g1p") ? "reco_n_protons == 1"
                                                                           : "reco_n_protons >= 1";
  const ana::tools::cut_sequence true_topology_cut = (topology == "1g0p") ? "true_n_protons == 0"
                                                    : (topology == "1g1p") ? "true_n_protons == 1"
                                                                           : "reco_n_protons >= 1";
  ana::tools::cut_sequence selection_cut = "reco_gOre_topology";
  ana::tools::cut_sequence signal_cut    = "true_mc_category == 0";
  selection_cut += reco_topology_cut;
  signal_cut    += true_topology_cut;

  ana::tools::analysis_tree my_analysis_tree(
    std::filesystem::current_path().string() + "/" + fileName,
    "events/" + sampleName,
    selection_cut, signal_cut,
    "selected_" + topology,
    (topology != "1gNp") ? "signal_" + topology : "signal_1gMp",
    sig_cats, sig_cats);

  // -------------------------------------------------------------------------
  // Variable registration
  //   Format: add_variable(name, nBins, min, max, title)
  // -------------------------------------------------------------------------
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~BINS~~~~MIN~~~~~~MAX~~~~~~TITLE
  my_analysis_tree.add_variable("reco_gOre_topology",                             2,  -0.5,     1.5,   "Reconstructed gOre Topology");
  my_analysis_tree.add_variable("reco_flash",                                     2,  -0.5,     1.5,   "Reconstructed Flash");
  my_analysis_tree.add_variable("reco_fiducial",                                  2,  -0.5,     1.5,   "Reconstructed Fiducial");
  my_analysis_tree.add_variable("reco_containment",                               2,  -0.5,     1.5,   "Reconstructed Containment");
  my_analysis_tree.add_variable("reco_flash_total_pe",                           50,     0,  100000,   "Total Flash PE");
  my_analysis_tree.add_variable("true_opening_costh",                            50,    -1,       1,   "True Opening Angle");
  my_analysis_tree.add_variable("reco_opening_costh",                            50,    -1,       1,   "Reconstructed Opening Angle");
  my_analysis_tree.add_variable("true_n_protons",                                 6,     0,       6,   "True Protons above Threshold");
  my_analysis_tree.add_variable("reco_n_protons",                                 6,     0,       6,   "Protons above Threshold");
  my_analysis_tree.add_variable("true_vertex_x",                                 50,  -400,     400,   "True Vertex X (cm)");
  my_analysis_tree.add_variable("true_vertex_y",                                 50,  -200,     200,   "True Vertex Y (cm)");
  my_analysis_tree.add_variable("true_vertex_z",                                 50, -1000,    1000,   "True Vertex Z (cm)");
  my_analysis_tree.add_variable("reco_vertex_x",                                 50,  -400,     400,   "Vertex X (cm)");
  my_analysis_tree.add_variable("reco_vertex_y",                                 50,  -200,     200,   "Vertex Y (cm)");
  my_analysis_tree.add_variable("reco_vertex_z",                                 50, -1000,    1000,   "Vertex Z (cm)");
  my_analysis_tree.add_variable("reco_gOre_gap",                                 50,     0,     100,   "#gamma-candidate Distance from Vertex (cm)");
  my_analysis_tree.add_variable("reco_pion_mass",                                50,     0,     250,   "Reconstructed Neutral Pion Mass Peak (MeV/c^{2}_{})");
  my_analysis_tree.add_variable("reco_delta_mass_P",                             25,   950,    1600,   "Reconstructed M_{#Delta} Resonance Peak (MeV/c^{2}_{})");
  my_analysis_tree.add_variable("reco_delta_mass_N",                             25,   950,    1600,   "Reconstructed M_{#Delta} Resonance Peak (MeV/c^{2}_{})");
  my_analysis_tree.add_variable("true_delta_mass_P",                             25,   950,    1600,   "True M_{#Delta} Resonance Peak (MeV/c^{2}_{})");
  my_analysis_tree.add_variable("true_delta_mass_N",                             25,   950,    1600,   "True M_{#Delta} Resonance Peak (MeV/c^{2}_{})");
  my_analysis_tree.add_variable("true_delta_mass_N_MC",                          25,   950,    1600,   "True M_{#Delta} Resonance Peak (MeV/c^{2}_{})");
  my_analysis_tree.add_variable("true_interaction_mode",                         15,    -1.5,    13.5, "GENIE Interaction Mode");
  my_analysis_tree.add_variable("true_interaction_type",                        101,   999.5, 1100.5,  "GENIE Interaction Type");
  my_analysis_tree.add_variable("true_baryon_res_code",                          19,    -1.5,    17.5, "Resonance Number");
  my_analysis_tree.add_variable("true_mc_category",                               8,    -0.5,     7.5, "MC Truth Category");
  my_analysis_tree.add_variable("true_category",                                  9,    -0.5,     8.5, "Topological Truth Category");
  // Leading shower
  my_analysis_tree.add_variable("reco_leading_primary_gOre_primary_softmax",    100,    0.99,     1,   "Shower Primary Softmax");
  my_analysis_tree.add_variable("reco_leading_primary_gOre_photon_softmax",      50,    0,        1,   "Shower Photon Softmax");
  my_analysis_tree.add_variable("reco_leading_primary_gOre_start_dedx",          50,    0,       25,   "Shower Start dE/dx (MeV/cm)");
  my_analysis_tree.add_variable("reco_leading_primary_gOre_axial_spread",        75,   -0.5,      1,   "Shower Axial Spread");
  my_analysis_tree.add_variable("reco_leading_primary_gOre_directional_spread", 100,    0,        1,   "Shower Directional Spread");
  my_analysis_tree.add_variable("reco_leading_primary_gOre_iou",                100,    0,        1,   "Shower IoU");
  my_analysis_tree.add_variable("true_leading_primary_gOre_iou",                100,    0,        1,   "True Shower IoU");
  my_analysis_tree.add_variable("reco_leading_primary_gOre_shower_ke",           50,    0,     1000,   "Shower KE (MeV)");
  my_analysis_tree.add_variable("true_leading_primary_gOre_shower_ke",           50,    0,     1000,   "True Shower KE (MeV)");
  my_analysis_tree.add_variable("reco_leading_primary_gOre_shower_dpT",         100,    0,     1000,   "Shower Transverse Momentum (MeV/c)");
  my_analysis_tree.add_variable("true_leading_primary_gOre_shower_dpT",         100,    0,     1000,   "True Shower Transverse Momentum (MeV/c)");
  my_analysis_tree.add_variable("reco_leading_primary_gOre_azimuthal_angle",     50,   -3.14,    3.14, "Shower Azimuthal Angle (rad)");
  my_analysis_tree.add_variable("true_leading_primary_gOre_azimuthal_angle",     50,   -3.14,    3.14, "True Shower Azimuthal Angle (rad)");
  my_analysis_tree.add_variable("reco_leading_primary_gOre_polar_angle",         50,    0,       3.14, "Shower Polar Angle (rad)");
  my_analysis_tree.add_variable("true_leading_primary_gOre_polar_angle",         50,    0,       3.14, "True Shower Polar Angle (rad)");
  // Subleading shower
  my_analysis_tree.add_variable("reco_subleading_primary_gOre_shower_ke",        25,    0,       25,   "Subleading Shower KE (MeV)");
  my_analysis_tree.add_variable("true_subleading_primary_gOre_shower_ke",        25,    0,       25,   "True Subleading Shower KE (MeV)");
  my_analysis_tree.add_variable("reco_subleading_primary_gOre_primary_softmax", 100,    0,        1,   "Subleading Shower Primary Softmax");
  my_analysis_tree.add_variable("reco_subleading_primary_gOre_photon_softmax",   50,    0,        1,   "Subleading Shower Photon Softmax");
  my_analysis_tree.add_variable("reco_subleading_primary_gOre_start_dedx",       50,    0,       25,   "Subleading Shower Start dE/dx (MeV/cm)");
  my_analysis_tree.add_variable("reco_subleading_primary_gOre_axial_spread",     75,   -0.5,      1,   "Subleading Shower Axial Spread");
  my_analysis_tree.add_variable("reco_subleading_primary_gOre_directional_spread",100,  0,        1,   "Subleading Shower Directional Spread");
  my_analysis_tree.add_variable("reco_subleading_primary_gOre_iou",             100,    0,        1,   "Subleading Shower IoU");
  my_analysis_tree.add_variable("true_subleading_primary_gOre_iou",             100,    0,        1,   "True Subleading Shower IoU");
  my_analysis_tree.add_variable("reco_subleading_primary_gOre_shower_dpT",       25,    0,       25,   "Subleading Shower Transverse Momentum (MeV/c)");
  my_analysis_tree.add_variable("true_subleading_primary_gOre_shower_dpT",       25,    0,       25,   "True Subleading Shower Transverse Momentum (MeV/c)");
  my_analysis_tree.add_variable("reco_subleading_primary_gOre_azimuthal_angle",  50,   -3.14,    3.14, "Subleading Shower Azimuthal Angle (rad)");
  my_analysis_tree.add_variable("true_subleading_primary_gOre_azimuthal_angle",  50,   -3.14,    3.14, "True Subleading Shower Azimuthal Angle (rad)");
  my_analysis_tree.add_variable("reco_subleading_primary_gOre_polar_angle",      50,    0,       3.14, "Subleading Shower Polar Angle (rad)");
  my_analysis_tree.add_variable("true_subleading_primary_gOre_polar_angle",      50,    0,       3.14, "True Subleading Shower Polar Angle (rad)");
  // Leading proton
  my_analysis_tree.add_variable("reco_leading_primary_proton_primary_softmax",  100,    0,        1,   "Leading Proton Primary Softmax");
  my_analysis_tree.add_variable("reco_leading_primary_proton_proton_softmax",    50,    0,        1,   "Leading Proton Softmax");
  my_analysis_tree.add_variable("reco_leading_primary_proton_ke",                50,    0,     1000,   "Proton KE (MeV)");
  my_analysis_tree.add_variable("true_leading_primary_proton_ke",                50,    0,     1000,   "True Proton KE (MeV)");
  my_analysis_tree.add_variable("reco_leading_primary_proton_dpT",              100,    0,     1000,   "Proton Transverse Momentum (MeV/c)");
  my_analysis_tree.add_variable("true_leading_primary_proton_dpT",              100,    0,     1000,   "True Proton Transverse Momentum (MeV/c)");
  my_analysis_tree.add_variable("reco_leading_primary_proton_length",           100,    0,      200,   "Proton Length (cm)");
  my_analysis_tree.add_variable("true_leading_primary_proton_length",           100,    0,      200,   "True Proton Length (cm)");
  my_analysis_tree.add_variable("reco_leading_primary_proton_azimuthal_angle",   50,   -3.14,    3.14, "Proton Azimuthal Angle (rad)");
  my_analysis_tree.add_variable("true_leading_primary_proton_azimuthal_angle",   50,   -3.14,    3.14, "True Proton Azimuthal Angle (rad)");
  my_analysis_tree.add_variable("reco_leading_primary_proton_polar_angle",       50,    0,       3.14, "Proton Polar Angle (rad)");
  my_analysis_tree.add_variable("true_leading_primary_proton_polar_angle",       50,    0,       3.14, "True Proton Polar Angle (rad)");
  // Composite kinematics
  my_analysis_tree.add_variable("abs(reco_delta_mass_P - 1232)",                 50,    0,      300,   "#vertReco M_{#Delta} - 1232 MeV/c^{2}_{}#vert");
  my_analysis_tree.add_variable("true_neutron_momentum",                          50,    0,     1000,   "True Neutron Momentum (MeV/c)");
  my_analysis_tree.add_variable("true_neutron_cosTh",                             50,    0,        1,   "True Photon-Neutron cos#Theta");

  // -------------------------------------------------------------------------
  // Helper lambda: plot var at the current cut state and save with stage prefix
  // -------------------------------------------------------------------------
  const std::string pdf_suffix = ".pdf";
  ana::tools::cut_sequence cut;
  ana::tools::cut_sequence truth_cut;

  auto plot_at_stage = [&](const std::string& var, const std::string& stage)
  {
    auto p = try_call(var, [&]{ return my_analysis_tree.plot_var_sel(var, cut); });
    p.PrintWIP("plots/"+sample+"/"+stage+"_"+var+pdf_suffix);
    return p;
  };

  std::cout << "\n*** " << sel << " ***\n" << std::endl;

  // =========================================================================
  // STAGE 1: PRESELECTION
  //   Apply: topology, n_protons, flash, fiducial, containment
  //
  //   [KEY - look at after preselection]
  //     reco_subleading_primary_gOre_shower_ke  — tells you how much pi0 background
  //       is still present; a long tail here means pi0 rejection has work to do
  //     reco_leading_primary_gOre_start_dedx    — early look at e/gamma mix entering
  //       the next two stages
  // =========================================================================
  std::cout << "//*** STAGE 1: PRESELECTION ***//" << std::endl;

  cut += "reco_gOre_topology == 1";
  cut += reco_topology_cut;
  cut += "reco_flash == 1";
  cut += "reco_fiducial == 1";
  cut += "reco_containment == 1";
  try_call(cut.string(), [&]{ my_analysis_tree.report_on_cut(cut); });

  // Preselection variables post-preselection (confirm selection is working)
  plot_at_stage("reco_n_protons",       "presel");
  plot_at_stage("reco_flash_total_pe",  "presel");
  plot_at_stage("reco_vertex_x",        "presel");
  plot_at_stage("reco_vertex_y",        "presel");
  plot_at_stage("reco_vertex_z",        "presel");

  // [KEY] Preview pi0 rejection variables — look for excess subleading shower events
  //       and how much of the delta mass window is filled by background
  plot_at_stage("reco_subleading_primary_gOre_shower_ke",          "presel"); // [KEY]
  plot_at_stage("reco_leading_primary_gOre_primary_softmax",       "presel");
  plot_at_stage("reco_leading_primary_proton_primary_softmax",     "presel");
  plot_at_stage("abs(reco_delta_mass_P - 1232)",                   "presel");
  plot_at_stage("reco_pion_mass",                                  "presel");

  // [KEY] Preview e/gamma separation — is signal already separating from CC-e background?
  plot_at_stage("reco_leading_primary_gOre_start_dedx",            "presel"); // [KEY]
  plot_at_stage("reco_leading_primary_gOre_photon_softmax",        "presel"); // [KEY]
  plot_at_stage("reco_leading_primary_gOre_axial_spread",          "presel");
  plot_at_stage("reco_leading_primary_gOre_directional_spread",    "presel");
  plot_at_stage("reco_leading_primary_gOre_primary_softmax",       "presel");
  plot_at_stage("reco_gOre_gap",                                   "presel");

  if (optimizeCuts)
  {
    // -----------------------------------------------------------------------
    // OPTIMIZE: scan cuts on all variables sequentially
    // -----------------------------------------------------------------------
    std::cout << "//*** OPTIMIZE: Subleading Shower Veto ***//" << std::endl;
    cut += "!(reco_subleading_primary_gOre_shower_ke > 0)";
    try_call(cut.string(), [&]{ my_analysis_tree.report_on_cut(cut); });

    for (auto const& [var, upper, lw_end, up_end] : vars_to_optimize)
      cut = optimize_cut(my_analysis_tree, cut, var, upper, lw_end, up_end, sample, pdf_suffix);

    for (auto const& [condition, var, upper, lw_end, up_end] : vars_to_optimize_with_condition)
      cut = optimize_conditional_cut(my_analysis_tree, cut, condition, var, upper, lw_end, up_end, sample, pdf_suffix);
  }
  else
  {
    // =========================================================================
    // STAGE 2: pi0 REJECTION
    //   Apply: subleading shower veto, shower primary softmax,
    //          proton primary softmax, delta mass window
    //
    //   [KEY - look at before applying these cuts]
    //     reco_subleading_primary_gOre_shower_ke  — the dominant pi0 handle;
    //       events with a second shower above threshold are mostly pi0→gg
    //     reco_leading_primary_gOre_primary_softmax — secondary pi0 handle;
    //       look at the tail below 0.9991 to understand what's being cut
    //
    //   [KEY - look at after these cuts]
    //     reco_leading_primary_gOre_start_dedx    — how much CC-e background
    //       remains at this point; this is the primary handle for Stage 3
    //     reco_leading_primary_gOre_photon_softmax — ditto
    // =========================================================================
    std::cout << "//*** STAGE 2: pi0 REJECTION ***//" << std::endl;

    // Plot pi0 rejection variables PRE-cut (entering this stage)
    plot_at_stage("reco_subleading_primary_gOre_shower_ke",          "pre_pi0rej"); // [KEY]
    plot_at_stage("reco_leading_primary_gOre_primary_softmax",       "pre_pi0rej");
    plot_at_stage("reco_leading_primary_proton_primary_softmax",     "pre_pi0rej");
    plot_at_stage("abs(reco_delta_mass_P - 1232)",                   "pre_pi0rej");
    {
      // Fit the pi0 mass peak before applying pion rejection
      auto pion_mass_plot = try_call("pion_mass_pre_pi0rej",
        [&]{ return my_analysis_tree.plot_var_sel("reco_pion_mass", cut); });
      pion_mass_plot.SumHist();
      TF1* pion_fit = new TF1("pion_fit_pre", landauPlusGauss, 0, 250, 6);
      pion_fit->SetParNames("LandauNorm", "LandauMean", "LandauSigma", "GausNorm", "GausMean", "GausSigma");
      pion_fit->SetParameters(pion_mass_plot.sum->GetMaximum(), 21, 7.8,
                              0.25*pion_mass_plot.sum->GetMaximum(), pion_mass, 10);
      pion_fit->SetParLimits(1, 0,   50);
      pion_fit->SetParLimits(2, 0,  250);
      pion_fit->SetParLimits(4, 50, 250);
      pion_fit->SetParLimits(5, 0,  250);
      TFitResultPtr res = pion_mass_plot.sum->Fit(pion_fit, "LS", "0 R");
      std::cout << "Pion mass fit — Chi^2: " << res->Chi2() << ", NDoF: " << res->Ndf() << std::endl;
      pion_mass_plot.canvas->cd();
      pion_mass_plot.stack->Draw("hist");
      pion_fit->Draw("L,same");
      pion_mass_plot.canvas->Update();
      pion_mass_plot.PrintWIP("plots/"+sample+"/pre_pi0rej_pion_mass_FIT"+pdf_suffix);
    }

    // Apply pi0 rejection cuts
    cut += "!(reco_subleading_primary_gOre_shower_ke > 0)";
    cut += "reco_leading_primary_gOre_primary_softmax > "    + std::to_string(cut_primary_softmax);
    cut += "reco_leading_primary_proton_primary_softmax > "  + std::to_string(cut_proton_primary_softmax);
    cut += "abs(reco_delta_mass_P - 1232) < "               + std::to_string(cut_delta_mass_window);
    try_call(cut.string(), [&]{ my_analysis_tree.report_on_cut(cut); });

    // Plot pi0 rejection variables POST-cut
    plot_at_stage("reco_subleading_primary_gOre_shower_ke",          "post_pi0rej");
    plot_at_stage("reco_leading_primary_gOre_primary_softmax",       "post_pi0rej");
    plot_at_stage("reco_leading_primary_proton_primary_softmax",     "post_pi0rej");
    plot_at_stage("abs(reco_delta_mass_P - 1232)",                   "post_pi0rej");
    plot_at_stage("reco_pion_mass",                                  "post_pi0rej");

    // [KEY] Preview e/gamma stage — look for the remaining CC-e contamination
    plot_at_stage("reco_leading_primary_gOre_start_dedx",            "post_pi0rej"); // [KEY]
    plot_at_stage("reco_leading_primary_gOre_photon_softmax",        "post_pi0rej"); // [KEY]
    plot_at_stage("reco_leading_primary_gOre_axial_spread",          "post_pi0rej");
    plot_at_stage("reco_leading_primary_gOre_directional_spread",    "post_pi0rej");
    plot_at_stage("reco_leading_primary_gOre_primary_softmax",       "post_pi0rej");
    plot_at_stage("reco_gOre_gap",                                   "post_pi0rej");

    // =========================================================================
    // STAGE 3: e/gamma SEPARATION
    //   Apply: start_dedx, axial_spread (no-op), directional_spread, photon_softmax
    //
    //   [KEY - look at before applying these cuts]
    //     reco_leading_primary_gOre_start_dedx    — gamma showers peak ~4-5 MeV/cm
    //       (Compton-driven), electrons peak ~2 MeV/cm (MIP-like). This single
    //       variable carries most of the e/gamma discrimination power.
    //     reco_leading_primary_gOre_photon_softmax — NN score trained directly on
    //       gamma vs electron; strongly correlated with start_dedx but independent
    //       information from shower shape features
    //
    //   Note: reco_leading_primary_gOre_axial_spread > -0.5 is equal to the
    //         histogram lower edge and is effectively a no-op in this selection.
    //
    //   [KEY - look at after these cuts]
    //     reco_delta_mass_P  — check the Delta(1232) peak is clean post-selection
    //     reco_n_protons     — confirm 1p topology is preserved
    // =========================================================================
    std::cout << "//*** STAGE 3: e/gamma SEPARATION ***//" << std::endl;

    // Plot e/gamma separation variables PRE-cut (entering this stage)
    plot_at_stage("reco_leading_primary_gOre_start_dedx",            "pre_egam"); // [KEY]
    plot_at_stage("reco_leading_primary_gOre_photon_softmax",        "pre_egam"); // [KEY]
    plot_at_stage("reco_leading_primary_gOre_axial_spread",          "pre_egam");
    plot_at_stage("reco_leading_primary_gOre_directional_spread",    "pre_egam");

    // Apply e/gamma separation cuts
    cut += "reco_leading_primary_gOre_start_dedx > "          + std::to_string(cut_start_dedx);
    cut += "reco_leading_primary_gOre_axial_spread > "        + std::to_string(cut_axial_spread_lo);
    cut += "reco_leading_primary_gOre_directional_spread < "  + std::to_string(cut_directional_spread_hi);
    cut += "reco_leading_primary_gOre_photon_softmax > "      + std::to_string(cut_photon_softmax);
    try_call(cut.string(), [&]{ my_analysis_tree.report_on_cut(cut); });

    // Plot e/gamma separation variables POST-cut
    plot_at_stage("reco_leading_primary_gOre_start_dedx",            "post_egam");
    plot_at_stage("reco_leading_primary_gOre_photon_softmax",        "post_egam");
    plot_at_stage("reco_leading_primary_gOre_axial_spread",          "post_egam");
    plot_at_stage("reco_leading_primary_gOre_directional_spread",    "post_egam");
    plot_at_stage("reco_leading_primary_gOre_primary_softmax",       "post_egam");
    plot_at_stage("reco_gOre_gap",                                   "post_egam");
    
    // [KEY] Final sanity checks — delta mass peak shape and proton multiplicity
    plot_at_stage("reco_delta_mass_P",                               "post_egam"); // [KEY]
    plot_at_stage("abs(reco_delta_mass_P - 1232)",                   "post_egam"); // [KEY]
    plot_at_stage("reco_n_protons",                                  "post_egam");
    plot_at_stage("reco_leading_primary_gOre_shower_ke",             "post_egam");
    plot_at_stage("reco_leading_primary_proton_ke",                  "post_egam");

  } // end !optimizeCuts

  // =========================================================================
  // Final variable loop — plot all registered variables at the final cut state
  // =========================================================================
  for (auto const& var : my_analysis_tree.variables())
  {
    auto var_plot =
      try_call("plot "+var,
        [&]{ return my_analysis_tree.plot_var_sel(var, cut); });
    auto var_plot_sig =
      try_call("plot "+var+" signal",
        [&]{ return my_analysis_tree.plot_var_sig(var, cut); });

    const std::string pdfName     = "plots/"+sample+"/"+var+pdf_suffix;
    const std::string pdfName_sig = "plots/"+sample+"/signal_"+var+pdf_suffix;

    // Apply alphanumeric axis labels for categorical truth variables
    if (var == "true_interaction_type" ||
        var == "true_interaction_mode" ||
        var == "true_baryon_res_code"  ||
        var == "true_mc_category"      ||
        var == "true_category"          )
    {
      using StrFromInt = std::function<std::string(int)>;
      StrFromInt get_label =
        (var == "true_interaction_type") ? StrFromInt([](int c){ return genie_inttypes.count(c) ? genie_inttypes.at(c) : ""; })
      : (var == "true_interaction_mode") ? StrFromInt([](int c){ return genie_modes.at(c); })
      : (var == "true_baryon_res_code")  ? StrFromInt([](int c){ return res_codes.at(c); })
      : (var == "true_mc_category")      ? StrFromInt([&mc_cats](int c){
                                               return (c == 0) ? std::string("NC #Delta#rightarrowN#gamma")
                                                               : mc_cats.at(c + 1);
                                           })
      : (var == "true_category")         ? StrFromInt([](int c){ return sel_cats.at(c).first; })
      :                                    StrFromInt([](int c){ return std::to_string(c); });

      const size_t nBins = var_plot.hists.front()->GetXaxis()->GetNbins();
      for (size_t bin = 1; bin <= nBins; ++bin)
      {
        const int bin_code = var_plot.hists.front()->GetXaxis()->GetBinCenter(bin);
        const std::string bin_label = get_label(bin_code);
        for (auto& h : var_plot.hists)     h->GetXaxis()->SetBinLabel(bin, bin_label.c_str());
        for (auto& h : var_plot_sig.hists) h->GetXaxis()->SetBinLabel(bin, bin_label.c_str());
      }

      const std::string title_str =
        (var == "true_interaction_type") ? "GENIE Interaction Type;;Events"
      : (var == "true_interaction_mode") ? "GENIE Interaction Mode;;Events"
      : (var == "true_baryon_res_code")  ? "Resonance;;Events"
      : (var == "true_category")         ? "Topological Category;;Events"
      : (var == "true_mc_category")      ? "MC Category;;Events"
      :                                    ";;Events";
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

    // Draw a vertical line at M_Delta = 1232 MeV/c^2 on delta mass plots
    if (var.find("delta") != std::string::npos)
    {
      var_plot.canvas->cd();
      var_plot.canvas->Update();
      const double ymax = std::ceil(var_plot.stack->GetMaximum() + std::sqrt(var_plot.stack->GetMaximum()));
      TLine massLine(1232, 0, 1232, ymax);
      massLine.SetLineColor(kBlack);
      massLine.SetLineStyle(kDashed);
      massLine.Draw();

      var_plot_sig.canvas->cd();
      var_plot_sig.canvas->Update();
      const double ymax_sig = std::ceil(var_plot_sig.stack->GetMaximum() + std::sqrt(var_plot_sig.stack->GetMaximum()));
      TLine massLine_sig(1232, 0, 1232, ymax_sig);
      massLine_sig.SetLineColor(kBlack);
      massLine_sig.SetLineStyle(kDashed);
      massLine_sig.Draw();

      var_plot.PrintWIP(pdfName);
      var_plot_sig.PrintWIP(pdfName_sig);
    }
    else
    {
      var_plot.PrintWIP(pdfName);
      var_plot_sig.PrintWIP(pdfName_sig);
    }
  }

  return 0;
}

// ===========================================================================
// main
// ===========================================================================
int main(int argc, char* argv[])
{
  gErrorIgnoreLevel = 3000;

  if (argc < 2) die("Must pass in the name of the ROOT file to analyze. Bail.");
  if (argc < 3) die("Must give name of sample in file to process.");
  if (argc < 4) die("Must give the name of the topology used.");

  const std::string suffix(".root");
  const std::string fileName(argv[1]);
  const std::string sampleName(argv[2]);
  const std::string topology(argv[3]);

  if (fileName.compare(fileName.size() - suffix.size(), suffix.size(), suffix) != 0)
    die("File must be a ROOT file (.root extension). Bail.");

  const std::string sampleBase = fileName.substr(0, fileName.size() - suffix.size());
  const bool optimizeCuts = false;

  return try_call("Run gOre Analysis",
                  run_analysis,
                  fileName, sampleBase, sampleName, topology,
                  "gOre", "gOre", optimizeCuts);
}

/**
 * @file analysis_tree.cc
 * @brief Implement the analysis_tree methods
 * @author hhausner@fnal.gov
 */

#include <iomanip>
#include <iostream>
#include <sstream>

#include "include/analysis_tree.h"

/**
 * @namespace ana::tools
 * @brief Namespace for tools helpful to running analyses off SPINE outputs
 */
namespace ana::tools
{
  /**
   * @brief Construct an analysis_tree
   * @details Pass in an analysis file (see macros/sbn/example.C) and wrap the TTrees in a way that makes developing an analysis easier.
   * You need to tell the analisys_tree where the signal/selected TTrees live, supply an initial selection cut, and
   * the signal/selection categories you are interested in analyzing.
   */
  analysis_tree::analysis_tree(const std::string& inFileName, 
                               const std::string& directory,
                               const ana::tools::cut_sequence& selection_cut,
                               const std::string& sel_tree,
                               const std::string& sig_tree,
                               const std::vector<std::pair<std::string, ana::tools::cut_sequence>> sel_cats,
                               const std::vector<std::pair<std::string, ana::tools::cut_sequence>> sig_cats) : inFile(TFile::Open(inFileName.c_str(), "READ"))
  {
    // Error Checks
    if ((inFile.get() == nullptr) || inFile->IsZombie())
      throw std::runtime_error("Could not open input file.");
    if (sel_cats.size() == 0)
      throw std::runtime_error("Must supply at least the signal category.");
    if (sig_cats.size() == 0)
      throw std::runtime_error("Must supply at least the signal category.");

    // Initialize Color Pallet
    // see https://jfly.uni-koeln.de/color/ for details on the pallet
    colors.push_back(new TColor(0.00, 0.00, 0.00)); // Black
    colors.push_back(new TColor(0.90, 0.60, 0.00)); // Orange
    colors.push_back(new TColor(0.35, 0.70, 0.90)); // Sky Blue
    colors.push_back(new TColor(0.00, 0.60, 0.50)); // Bluish Green
    colors.push_back(new TColor(0.95, 0.90, 0.25)); // Yellow
    colors.push_back(new TColor(0.00, 0.45, 0.70)); // Blue
    colors.push_back(new TColor(0.80, 0.40, 0.00)); // Vermillion
    colors.push_back(new TColor(0.80, 0.60, 0.70)); // Reddish Purple
    colors.push_back(new TColor(0.60, 0.50, 0.80)); // An Extra Purple (not from OG pallet)
    // pad colors with grays
    for (double gry = 0.10; gry < 0.55; gry += 0.10)
    {
      // alternate to maximize contrast
      colors.push_back(new TColor(0.05 + gry, 0.05 + gry, 0.05 + gry));
      colors.push_back(new TColor(1.05 - gry, 1.05 - gry, 1.05 - gry));
    }
    //colors.push_back(new TColor(0.25, 0.25, 0.25)); // Gray (light gray)
    //colors.push_back(new TColor(0.45, 0.45, 0.45)); // Gray (just gray)
    //colors.push_back(new TColor(0.85, 0.85, 0.85)); // Gray (but darker)

    // Memory Handling
    std::string tmpName = "temp_TTrees-"+sel_tree+"-"+sig_tree+"_Cuts-"+sel_cats.front().second.string()+"-"+sig_cats.front().second.string()+".root";
    memFile = std::make_unique<TMemFile>(tmpName.c_str(), "RECREATE");
    memFile->cd();

    // Read In Trees
    std::cout << "Reading from " << inFile->GetName() << std::endl;
    std::string selection_tree_name = directory + "/" + sel_tree;
    std::cout << "Looking for selection in " << selection_tree_name << std::endl;
    std::string signal_tree_name    = directory + "/" + sig_tree;
    std::cout << "Looking for signal in " << signal_tree_name << std::endl;
    ana::tools::cut_sequence signal_cut    = sel_cats.front().second;
    ana::tools::cut_sequence bkg_cut       = ana::tools::negate_cut(signal_cut); 
    TTree* tmpSignalTree = inFile->Get<TTree>(signal_tree_name.c_str());
    signalTree = std::unique_ptr<TTree>(tmpSignalTree->CopyTree(signal_cut.c_str()));
    signalTree->SetDirectory(memFile.get());
    TTree* selectedTree = inFile->Get<TTree>(selection_tree_name.c_str());
    selectedSignalTree = std::unique_ptr<TTree>(signalTree->CopyTree(selection_cut.c_str()));
    selectedSignalTree->SetDirectory(memFile.get());
    backgroundTree   = std::unique_ptr<TTree>(selectedTree->CopyTree(bkg_cut.c_str()));
    backgroundTree->SetDirectory(memFile.get());
    total_signal = signalTree->GetEntries();

    // Store Selection/Signal categories
    for (auto const& [sel_cat_name, sel_cat_cut] : sel_cats)
    {
      sel_cat_labels.push_back(sel_cat_name);
      sel_cat_cuts  .push_back(sel_cat_cut);
    }
    for (auto const& [sig_cat_name, sig_cat_cut] : sig_cats)
    {
      sig_cat_labels.push_back(sig_cat_name);
      sig_cat_cuts  .push_back(sig_cat_cut);
    }
  }

  /**
   * @brief Add a vairable to the analysis
   * @details If you want to plot a variable it needs a name, a binning scheme, and a title for plotting
   */
  void analysis_tree::add_variable(const std::string& name,
                                   const unsigned int& bins,
                                   const double& lw, const double& up,
                                   const std::string& title)
  {
    if (var_to_range.count(name) != 0 ||
        var_to_title.count(name) != 0)
      throw std::runtime_error("analysis_tree::add_variable -- Variable " + name + " already added to tree");
    vars.push_back(name);
    var_to_range.emplace(name, std::make_tuple(bins, lw, up));
    var_to_title.emplace(name, title);
  }

  /**
   * @brief Optimize a threshod for a signal definition
   * @details This modifies the signal cut, and so changes the number of signal events and migrates events between the
   * selected signal and background neutrino trees. It also assumes you want to update the selection to match the signal change
   */
  std::tuple<ana::tools::opt_canvas, ana::tools::limit_canvas> analysis_tree::optimize_threshold(const std::string& var,
                                                                                                 const std::string& true_var,
                                                                                                 const double& lw, const double& up,
                                                                                                 const ana::tools::cut_sequence& old_cut,
                                                                                                 const ana::tools::cut_sequence& old_truth_cut)
  {
    ana::tools::opt_canvas opt;
    ana::tools::limit_canvas lc = plot_var_sel(var, old_cut);
    opt.canvas = std::make_shared<TCanvas>((var+"_threshold").c_str(),
                                           (var+"_threshold").c_str(), 1618, 1000);
    opt.canvas->cd();
    opt.canvas->SetBit(kCanDelete, false);
    opt.purity_graph = std::make_shared<TGraph>();
    opt.purity_graph->SetBit(kCanDelete, false);
    opt.efficiency_graph = std::make_shared<TGraph>();
    opt.efficiency_graph->SetBit(kCanDelete, false);
    opt.fom_graph = std::make_shared<TGraph>();
    opt.fom_graph->SetBit(kCanDelete, false);
    opt.purity_efficiency_fom = std::make_shared<TMultiGraph>((var+"_threshold_epf").c_str(),
                                                              (var+"_threshold_epf").c_str());
    opt.purity_efficiency_fom->SetBit(kCanDelete, false);
    opt.is_upper_bound = false;
    unsigned int nSteps = 1000;
    double step = (up - lw) / nSteps;
    opt.fom = std::numeric_limits<double>::lowest();
    opt.limit = std::numeric_limits<double>::lowest();
    for (double limit = lw; limit < up; limit += step)
    {
      std::ostringstream limStrm;
      limStrm << std::setprecision(10) << limit;
      ana::tools::cut_sequence cut = old_cut.with_addition(var + " < " + limStrm.str());
      ana::tools::cut_sequence truth_cut = old_truth_cut.with_addition(true_var + " < " + limStrm.str());
      auto [new_total_signal, selected_signal, selected_background, efficiency, purity] = signal_def_e_p(cut, truth_cut);
      double fom = get_fom(selected_signal, selected_background, efficiency, purity);
      opt.purity_graph->AddPoint(limit, purity);
      opt.efficiency_graph->AddPoint(limit, efficiency);
      opt.fom_graph->AddPoint(limit, fom);
      if (opt.fom < fom)
      {
        opt.fom = fom;
        opt.cut = cut;
        opt.truth_cut = truth_cut;
        opt.limit = limit;
      }
    }

    // update signal/selection trees
    // ensure the events which are selected but no longer signal are migrated properly
    std::unique_ptr<TTree> newSignalTree = std::unique_ptr<TTree>(signalTree->CopyTree(opt.truth_cut.c_str()));
    newSignalTree->SetDirectory(memFile.get());
    signalTree = std::move(newSignalTree);
    std::unique_ptr<TTree> newSelectedSignalTree = std::unique_ptr<TTree>(selectedSignalTree->CopyTree(opt.cut.with_addition(opt.truth_cut).c_str()));
    newSelectedSignalTree->SetDirectory(memFile.get());
    selectedSignalTree = std::move(newSelectedSignalTree);
    TTree* bkgFromBkg = backgroundTree->CopyTree(opt.cut.c_str());
    bkgFromBkg->SetDirectory(memFile.get());
    TTree* bkgFromSig = selectedSignalTree->CopyTree(opt.cut.with_addition(ana::tools::negate_cut(opt.truth_cut)).c_str());
    bkgFromSig->SetDirectory(memFile.get());
    TList bkgList;
    bkgList.Add(bkgFromBkg);
    bkgList.Add(bkgFromSig);
    std::unique_ptr<TTree> newNuBackgroundTree = std::unique_ptr<TTree>(TTree::MergeTrees(&bkgList));
    newNuBackgroundTree->SetDirectory(memFile.get());
    backgroundTree = std::move(newNuBackgroundTree);
    total_signal = signalTree->GetEntries();
    
    for (size_t idx = 0; idx < opt.fom_graph->GetN(); ++idx)
    {
      double fom, bound;
      opt.fom_graph->GetPoint(idx, bound, fom);
      opt.fom_graph->SetPoint(idx, bound, fom / std::ceil(opt.fom));
    }
    opt.purity_graph->SetLineColor(colors.at(1)->GetNumber());
    opt.efficiency_graph->SetLineColor(colors.at(2)->GetNumber());
    opt.fom_graph->SetLineColor(colors.at(3)->GetNumber());
    opt.legend = std::make_shared<TLegend>(0.6, 0.6, 0.8, 0.8);
    opt.legend->SetBit(kCanDelete, false);
    opt.legend->SetFillStyle(0);
    opt.legend->AddEntry(opt.purity_graph.get(), "Purity", "l");
    opt.legend->AddEntry(opt.efficiency_graph.get(), "Efficiency", "l");
    opt.legend->AddEntry(opt.fom_graph.get(), "Figure of Merit", "l");
    opt.purity_efficiency_fom->Add(opt.purity_graph.get());
    opt.purity_efficiency_fom->Add(opt.efficiency_graph.get());
    opt.purity_efficiency_fom->Add(opt.fom_graph.get());
    opt.purity_efficiency_fom->GetXaxis()->SetTitle((var_to_title.at(var) + "Threshold").c_str());
    opt.purity_efficiency_fom->GetYaxis()->SetTitle("Purity, Efficiency");
    opt.purity_efficiency_fom->Draw("al");
    opt.purity_efficiency_fom->GetHistogram()->GetXaxis()->SetRangeUser(lw, up);
    opt.purity_efficiency_fom->SetMinimum(0);
    opt.purity_efficiency_fom->SetMaximum(1);
    opt.canvas->Update();
    opt.line = std::make_shared<TLine>(opt.limit, 0, opt.limit, 1);
    opt.line->SetBit(kCanDelete, false);
    opt.line->SetLineColor(kRed);
    opt.line->SetLineStyle(kDashed);
    opt.legend->Draw("same");
    opt.line->Draw("same");
    double bin_height = 0;
    double arrow_length = 0;
    for (TObject* obj : *static_cast<TList*>(lc.stack->GetHists()))
    {
      TH1* hist = static_cast<TH1*>(obj);
      int bin = hist->FindBin(opt.limit);
      bin_height += hist->GetBinContent(bin);
      if (arrow_length == 0)
        arrow_length = 5*(hist->GetXaxis()->GetTickLength());
    }
    lc.arrow = std::make_shared<TArrow>(opt.limit, bin_height,
                                        opt.limit + arrow_length, bin_height, 0.05, "|->");
    lc.arrow->SetBit(kCanDelete, false);
    lc.arrow->SetLineColor(kRed);
    lc.arrow->SetLineWidth(2);
    lc.arrow->SetLineStyle(kDashed);
    opt.gaxis = std::make_shared<TGaxis>(up, 0, up, 1, 0, std::ceil(opt.fom), 510, "+L");
    opt.gaxis->SetBit(kCanDelete, false);
    opt.gaxis->SetTitle("Figure of Merit");
    opt.gaxis->SetTitleSize(opt.purity_efficiency_fom->GetHistogram()->GetYaxis()->GetTitleSize());
    opt.gaxis->SetTitleOffset(0.8);
    opt.gaxis->SetTitleFont(opt.purity_efficiency_fom->GetHistogram()->GetYaxis()->GetTitleFont());
    opt.gaxis->SetLabelSize(opt.purity_efficiency_fom->GetHistogram()->GetYaxis()->GetLabelSize());
    opt.gaxis->SetLabelOffset(opt.purity_efficiency_fom->GetHistogram()->GetYaxis()->GetLabelOffset());
    opt.gaxis->SetLabelFont(opt.purity_efficiency_fom->GetHistogram()->GetYaxis()->GetLabelFont());
    opt.gaxis->SetTickLength(opt.purity_efficiency_fom->GetHistogram()->GetYaxis()->GetTickLength());
    opt.gaxis->SetNdivisions(opt.purity_efficiency_fom->GetHistogram()->GetYaxis()->GetNdivisions());
    opt.gaxis->Draw();
    lc.canvas->cd();
    lc.arrow->Draw();
    return std::make_tuple(opt, lc);
  }

  /**
   * @brief Return the variables currently in the analysis
   * @details Useful for looping over variables when making plots
   */
  std::vector<std::string> analysis_tree::variables() const
  {
    return vars;
  }

  /**
   * @brief Calculate the purity and efficiency for a new signal definition
   * @details The new signal definition is a subset of the original (no getting back thrown away signal events)
   * @return A tuple of the new total signal, the new selected signal, the new selected background, the efficiency, and the purity
   */
  std::tuple<unsigned int, unsigned int, unsigned int, double, double> analysis_tree::signal_def_e_p(const ana::tools::cut_sequence& cut,
                                                                                                     const ana::tools::cut_sequence& truth_cut) const
  {
    unsigned int new_total_signal     = signalTree        ->GetEntries(truth_cut.c_str());
    unsigned int selected_signal      = selectedSignalTree->GetEntries(cut.with_addition(truth_cut).c_str());
    unsigned int selected_background  = signalTree        ->GetEntries(cut.with_addition(ana::tools::negate_cut(truth_cut)).c_str());
                 selected_background += backgroundTree  ->GetEntries(cut.c_str());
    double efficiency = static_cast<double>(selected_signal) / new_total_signal;
    double purity     = static_cast<double>(selected_signal) / (selected_signal + selected_background);
    return std::make_tuple(new_total_signal, selected_signal, selected_background, efficiency, purity);
  }

  /**
   * @brief Calculate the purity and efficiency for a given selection
   * @return A tuple of the selected signal, the selected background, the efficiency, and the purity
   */
  std::tuple<unsigned int, unsigned int, double, double> analysis_tree::cut_e_p(const ana::tools::cut_sequence& cut) const
  {
    unsigned int selected_signal      = selectedSignalTree->GetEntries(cut.c_str());
    unsigned int selected_background  = backgroundTree  ->GetEntries(cut.c_str());
    double efficiency = static_cast<double>(selected_signal) / total_signal;
    double purity     = static_cast<double>(selected_signal) / (selected_signal + selected_background);
    return std::make_tuple(selected_signal, selected_background, efficiency, purity);
  }

  /**
   * @brief Report on how the cut performs
   * @details Pipe the cut string, the selected signal/background events, the background categorization breakdown, the efficiency,
   * and the purity to std::cout. Note that that the signal category (Cat 0) is included in the background breakdown, but it should
   * be empty.
   */
  void analysis_tree::report_on_cut(const ana::tools::cut_sequence& cut) const
  {
    auto [selected_signal, selected_background, efficiency, purity] = cut_e_p(cut);
    std::vector<unsigned int> bkg_cat_counts;
    size_t nCats = sel_cat_labels.size();
    std::ostringstream cat_report;
    for (size_t cat = 0; cat < nCats; ++cat)
    {
      ana::tools::cut_sequence cat_cut = cut.with_addition(sel_cat_cuts.at(cat));
      bkg_cat_counts.push_back(backgroundTree->GetEntries(cat_cut.c_str()));
      cat_report << "    Cat " << cat << " (" << sel_cat_labels.at(cat) << "): " << bkg_cat_counts.at(cat);
      if (cat != nCats - 1)
        cat_report << '\n';
    }
    std::cout
      << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << '\n' 
      << cut.string()                                            << '\n'
      << "  Total Signal: "               << total_signal        << '\n'
      << "  Selected Signal Events: "     << selected_signal     << '\n'
      << "  Selected Background Events: " << selected_background << '\n'
      <<    cat_report.str()                                     << '\n'
      << "  efficiency: "                 << efficiency          << '\n'
      << "  purity: "                     << purity              << '\n'
      << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  }

  /**
   * @brief Plot the distribution of selected events in a given variable for a given selection and return a canvas which can be saved
   * @details The variable must have already been added to the variable list with analysis_tree::add_variable
   */
  ana::tools::stack_canvas analysis_tree::plot_var_sel(const std::string& var, const ana::tools::cut_sequence& cut) const
  {
    ana::tools::stack_canvas sc;
    sc.canvas = std::make_shared<TCanvas>(var.c_str(), var.c_str(), 1618, 1000);
    sc.canvas->SetBit(kCanDelete, false);
    sc.canvas->cd();

    auto [bins, lw, up] = var_to_range.at(var); 
    std::shared_ptr<TH1F> sel_sig = std::make_shared<TH1F>("sel_sig", "", bins, lw, up);
    selectedSignalTree->Draw((var+">>sel_sig").c_str(), cut.c_str(), "goff");
    sel_sig->SetDirectory(nullptr);
    sc.hists.push_back(std::move(sel_sig));
    std::vector<std::shared_ptr<TH1F>> bkg_hists;
    size_t nCats = sel_cat_labels.size();
    for (size_t cat = 0; cat < nCats; ++cat)
    {
      std::string histName = "bkg_" + std::to_string(cat);
      ana::tools::cut_sequence cut_cat = cut.with_addition(sel_cat_cuts.at(cat));
      std::shared_ptr<TH1F> tmp_hist = std::make_shared<TH1F>(histName.c_str(), "", bins, lw, up);
      backgroundTree->Draw((var+">>"+histName).c_str(), cut_cat.c_str(), "goff");
      tmp_hist->SetDirectory(nullptr);
      bkg_hists.push_back(tmp_hist);
      if (bkg_hists.at(cat) == nullptr)
        throw std::runtime_error("bkg_hists for cat " + std::to_string(cat) + " is null");
    }
    for (size_t cat = 1; cat < nCats; ++cat)
      sc.hists.push_back(std::move(bkg_hists.at(cat)));
    sc.stack = std::make_shared<THStack>((var+"_stack").c_str(), (var+"_stack").c_str());
    sc.stack->SetBit(kCanDelete, false);
    sc.legend = std::make_shared<TLegend>(0.6, 0.6, 0.8, 0.8);
    sc.legend->SetBit(kCanDelete, false);
    sc.legend->SetFillStyle(0);
    sc.stack->SetTitle((";"+var_to_title.at(var)+";Events").c_str());
    for (size_t cat = 0; cat < nCats; ++cat)
    {
      sc.hists.at(nCats - cat - 1)->SetFillColor(colors.at(nCats - cat)->GetNumber());
      sc.hists.at(nCats - cat - 1)->SetLineColor(kBlack);
      sc.stack->Add(sc.hists.at(nCats - cat - 1).get());
      sc.legend->AddEntry(sc.hists.at(cat).get(), sel_cat_labels.at(cat).c_str(), "f");
    }
    sc.stack->Draw("hist");
    sc.legend->Draw("same");
    return sc;
  }

  /**
   * @brief Plot the distribution of signal events in a given variable for a given selection and return a canvas which can be saved
   * @details The variable must have already been added to the variable list with analysis_tree::add_variable
   */
  ana::tools::stack_canvas analysis_tree::plot_var_sig(const std::string& var, const ana::tools::cut_sequence& cut) const
  {
    ana::tools::stack_canvas sc;
    sc.canvas = std::make_shared<TCanvas>(var.c_str(), var.c_str(), 1618, 1000);
    sc.canvas->SetBit(kCanDelete, false);
    sc.canvas->cd();

    auto [bins, lw, up] = var_to_range.at(var); 
    size_t nCats = sig_cat_labels.size();
    for (size_t cat = 0; cat < nCats; ++cat)
    {
      std::string histName = "sel_sig_" + std::to_string(cat);
      ana::tools::cut_sequence cut_cat = cut.with_addition(sig_cat_cuts.at(cat));
      std::shared_ptr<TH1F> tmp_hist = std::make_shared<TH1F>(histName.c_str(), "", bins, lw, up);
      signalTree->Draw((var+">>"+histName).c_str(), cut_cat.c_str(), "goff");
      tmp_hist->SetFillColor(colors.at(nCats + cat)->GetNumber());
      tmp_hist->SetLineColor(kBlack);
      tmp_hist->SetDirectory(nullptr);
      sc.hists.push_back(tmp_hist);
      if (sc.hists.at(cat) == nullptr)
        throw std::runtime_error("sel_sig_hists for cat " + std::to_string(cat) + " is null");
    }
    for (size_t cat = 0; cat < nCats; ++cat)
    {
      std::string histName = "unsel_sig_" + std::to_string(cat);
      ana::tools::cut_sequence cut_cat = ana::tools::negate_cut(cut).with_addition(sig_cat_cuts.at(cat));
      std::shared_ptr<TH1F> tmp_hist = std::make_shared<TH1F>(histName.c_str(), "", bins, lw, up);
      signalTree->Draw((var+">>"+histName).c_str(), cut_cat.c_str(), "goff");
      tmp_hist->SetFillColor(colors.at(2*nCats + cat)->GetNumber());
      tmp_hist->SetLineColor(kBlack);
      tmp_hist->SetDirectory(nullptr);
      sc.hists.push_back(tmp_hist);
      if (sc.hists.at(nCats + cat) == nullptr)
        throw std::runtime_error("unsel_hists for cat " + std::to_string(cat) + " is null");
    }

    sc.stack = std::make_shared<THStack>((var+"_stack").c_str(), (var+"_stack").c_str());
    sc.stack->SetBit(kCanDelete, false);
    sc.legend = std::make_shared<TLegend>(0.6, 0.6, 0.8, 0.8);
    sc.legend->SetBit(kCanDelete, false);
    sc.legend->SetFillStyle(0);
    sc.stack->SetTitle((";"+var_to_title.at(var)+";Events").c_str());
    for (size_t cat = 0; cat < 2*nCats; ++cat)
    {
      sc.stack->Add(sc.hists.at(2*nCats - cat - 1).get());
      std::string label = (cat < nCats) ? "Selected " + sig_cat_labels.at(cat) : "Unselected " + sig_cat_labels.at(cat - nCats);
      sc.legend->AddEntry(sc.hists.at(cat).get(), label.c_str(), "f");
    }
    sc.stack->Draw("hist");
    sc.legend->Draw("same");
    return sc;
  }

  /**
   * @brief Get the figure of merit for cut optimization based on the selected signal, selected background, purity, & efficiency
   * @details Optimizes for signal over uncertainty
   * TODO: allow for alternative figures of merit
   */
  double analysis_tree::get_fom(const unsigned int& selected_signal,
                                const unsigned int& selected_background,
                                const double& efficiency, 
                                const double& purity) const
  {
    double fom(0.0);
    // by default to optimizing for signal over uncertainty
    if (selected_signal + selected_background > 0)
      fom = static_cast<double>(selected_signal) / std::sqrt(selected_signal + selected_background);
    //if (purity + efficiency > 0)
    //  fom = 2. * purity * efficiency / (purity + efficiency);
    return fom;
  }

  /**
   * @brief Optimize a cut which is a lower bound on a given variable
   * @details Requires the range over which to optimize as well as the selection to build on
   */
  std::tuple<ana::tools::opt_canvas, ana::tools::limit_canvas> analysis_tree::optimize_lower_bound(const std::string& var,
                                                                                                   const double& lw, const double& up,
                                                                                                   const ana::tools::cut_sequence& old_cut) const
  {
    ana::tools::opt_canvas opt;
    ana::tools::limit_canvas lc = plot_var_sel(var, old_cut);
    opt.canvas = std::make_shared<TCanvas>((var+"_lower_bound").c_str(),
                                           (var+"_lower_bound").c_str(), 1618, 1000);
    opt.canvas->cd();
    opt.canvas->SetBit(kCanDelete, false);
    opt.purity_graph = std::make_shared<TGraph>();
    opt.purity_graph->SetBit(kCanDelete, false);
    opt.efficiency_graph = std::make_shared<TGraph>();
    opt.efficiency_graph->SetBit(kCanDelete, false);
    opt.fom_graph = std::make_shared<TGraph>();
    opt.fom_graph->SetBit(kCanDelete, false);
    opt.purity_efficiency_fom = std::make_shared<TMultiGraph>((var+"_lower_bound_epf").c_str(),
                                                              (var+"_lower_bound_epf").c_str());
    opt.purity_efficiency_fom->SetBit(kCanDelete, false);
    opt.is_upper_bound = false;
    unsigned int nSteps = 1000;
    double step = (up - lw) / nSteps;
    opt.fom = std::numeric_limits<double>::lowest();
    opt.limit = std::numeric_limits<double>::lowest();
    for (double limit = lw; limit < up; limit += step)
    {
      std::ostringstream limStrm;
      limStrm << std::setprecision(10) << limit;
      ana::tools::cut_sequence cut = old_cut.with_addition(var + " > " + limStrm.str());
      auto [selected_signal, selected_background, efficiency, purity] = cut_e_p(cut);
      double fom = get_fom(selected_signal, selected_background, efficiency, purity);
      opt.purity_graph->AddPoint(limit, purity);
      opt.efficiency_graph->AddPoint(limit, efficiency);
      opt.fom_graph->AddPoint(limit, fom);
      if (opt.fom < fom)
      {
        opt.fom = fom;
        opt.cut = cut;
        opt.limit = limit;
      }
    }
    for (size_t idx = 0; idx < opt.fom_graph->GetN(); ++idx)
    {
      double fom, bound;
      opt.fom_graph->GetPoint(idx, bound, fom);
      opt.fom_graph->SetPoint(idx, bound, fom / std::ceil(opt.fom));
    }
    opt.purity_graph->SetLineColor(colors.at(1)->GetNumber());
    opt.efficiency_graph->SetLineColor(colors.at(2)->GetNumber());
    opt.fom_graph->SetLineColor(colors.at(3)->GetNumber());
    opt.legend = std::make_shared<TLegend>(0.6, 0.6, 0.8, 0.8);
    opt.legend->SetBit(kCanDelete, false);
    opt.legend->SetFillStyle(0);
    opt.legend->AddEntry(opt.purity_graph.get(), "Purity", "l");
    opt.legend->AddEntry(opt.efficiency_graph.get(), "Efficiency", "l");
    opt.legend->AddEntry(opt.fom_graph.get(), "Figure of Merit", "l");
    opt.purity_efficiency_fom->Add(opt.purity_graph.get());
    opt.purity_efficiency_fom->Add(opt.efficiency_graph.get());
    opt.purity_efficiency_fom->Add(opt.fom_graph.get());
    opt.purity_efficiency_fom->GetXaxis()->SetTitle(("Lower Bound on " + var_to_title.at(var)).c_str());
    opt.purity_efficiency_fom->GetYaxis()->SetTitle("Purity, Efficiency");
    opt.purity_efficiency_fom->Draw("al");
    opt.purity_efficiency_fom->GetHistogram()->GetXaxis()->SetRangeUser(lw, up);
    opt.purity_efficiency_fom->SetMinimum(0);
    opt.purity_efficiency_fom->SetMaximum(1);
    opt.canvas->Update();
    opt.line = std::make_shared<TLine>(opt.limit, 0, opt.limit, 1);
    opt.line->SetBit(kCanDelete, false);
    opt.line->SetLineColor(kRed);
    opt.line->SetLineStyle(kDashed);
    opt.legend->Draw("same");
    opt.line->Draw("same");
    double bin_height = 0;
    double arrow_length = 0;
    for (TObject* obj : *static_cast<TList*>(lc.stack->GetHists()))
    {
      TH1* hist = static_cast<TH1*>(obj);
      int bin = hist->FindBin(opt.limit);
      bin_height += hist->GetBinContent(bin);
      if (arrow_length == 0)
        arrow_length = 5*(hist->GetXaxis()->GetTickLength());
    }
    lc.arrow = std::make_shared<TArrow>(opt.limit, bin_height,
                                        opt.limit + arrow_length, bin_height, 0.05, "|->");
    lc.arrow->SetBit(kCanDelete, false);
    lc.arrow->SetLineColor(kRed);
    lc.arrow->SetLineWidth(2);
    lc.arrow->SetLineStyle(kDashed);
    opt.gaxis = std::make_shared<TGaxis>(up, 0, up, 1, 0, std::ceil(opt.fom), 510, "+L");
    opt.gaxis->SetBit(kCanDelete, false);
    opt.gaxis->SetTitle("Figure of Merit");
    opt.gaxis->SetTitleSize(opt.purity_efficiency_fom->GetHistogram()->GetYaxis()->GetTitleSize());
    opt.gaxis->SetTitleOffset(0.8);
    opt.gaxis->SetTitleFont(opt.purity_efficiency_fom->GetHistogram()->GetYaxis()->GetTitleFont());
    opt.gaxis->SetLabelSize(opt.purity_efficiency_fom->GetHistogram()->GetYaxis()->GetLabelSize());
    opt.gaxis->SetLabelOffset(opt.purity_efficiency_fom->GetHistogram()->GetYaxis()->GetLabelOffset());
    opt.gaxis->SetLabelFont(opt.purity_efficiency_fom->GetHistogram()->GetYaxis()->GetLabelFont());
    opt.gaxis->SetTickLength(opt.purity_efficiency_fom->GetHistogram()->GetYaxis()->GetTickLength());
    opt.gaxis->SetNdivisions(opt.purity_efficiency_fom->GetHistogram()->GetYaxis()->GetNdivisions());
    opt.gaxis->Draw();
    lc.canvas->cd();
    lc.arrow->Draw();
    return std::make_tuple(opt, lc);
  }

  /**
   * @brief Optimize a cut which is an upper bound on a given variable
   * @details Requires the range over which to optimize as well as the selection to build on
   */
  std::tuple<ana::tools::opt_canvas, ana::tools::limit_canvas> analysis_tree::optimize_upper_bound(const std::string& var,
                                                                                                   const double& lw, const double& up,
                                                                                                   const ana::tools::cut_sequence& old_cut) const
  {
    ana::tools::opt_canvas opt;
    ana::tools::limit_canvas lc = plot_var_sel(var, old_cut);
    opt.canvas = std::make_shared<TCanvas>((var+"_upper_bound").c_str(),
                                           (var+"_upper_bound").c_str(), 1618, 1000);
    opt.canvas->cd();
    opt.canvas->SetBit(kCanDelete, false);
    opt.purity_graph = std::make_shared<TGraph>();
    opt.purity_graph->SetBit(kCanDelete, false);
    opt.efficiency_graph = std::make_shared<TGraph>();
    opt.efficiency_graph->SetBit(kCanDelete, false);
    opt.fom_graph = std::make_shared<TGraph>();
    opt.fom_graph->SetBit(kCanDelete, false);
    opt.purity_efficiency_fom = std::make_shared<TMultiGraph>((var+"_upper_bound_epf").c_str(),
                                                              (var+"_upper_bound_epf").c_str());
    opt.purity_efficiency_fom->SetBit(kCanDelete, false);
    opt.is_upper_bound = true;
    unsigned int nSteps = 1000;
    double step = (up - lw) / nSteps;
    opt.fom = std::numeric_limits<double>::lowest();
    opt.limit = std::numeric_limits<double>::lowest();
    for (double limit = lw; limit < up; limit += step)
    {
      std::ostringstream limStrm;
      limStrm << std::setprecision(10) << limit;
      ana::tools::cut_sequence cut = old_cut.with_addition(var + " < " + limStrm.str());
      auto [selected_signal, selected_background, efficiency, purity] = cut_e_p(cut);
      double fom = get_fom(selected_signal, selected_background, efficiency, purity);
      opt.purity_graph->AddPoint(limit, purity);
      opt.efficiency_graph->AddPoint(limit, efficiency);
      opt.fom_graph->AddPoint(limit, fom);
      if (opt.fom < fom)
      {
        opt.fom = fom;
        opt.cut = cut;
        opt.limit = limit;
      }
    }
    for (size_t idx = 0; idx < opt.fom_graph->GetN(); ++idx)
    {
      double fom, bound;
      opt.fom_graph->GetPoint(idx, bound, fom);
      opt.fom_graph->SetPoint(idx, bound, fom / std::ceil(opt.fom));
    }
    opt.purity_graph->SetLineColor(colors.at(1)->GetNumber());
    opt.efficiency_graph->SetLineColor(colors.at(2)->GetNumber());
    opt.fom_graph->SetLineColor(colors.at(3)->GetNumber());
    opt.legend = std::make_shared<TLegend>(0.6, 0.6, 0.8, 0.8);
    opt.legend->SetBit(kCanDelete, false);
    opt.legend->SetFillStyle(0);
    opt.legend->AddEntry(opt.purity_graph.get(), "Purity", "l");
    opt.legend->AddEntry(opt.efficiency_graph.get(), "Efficiency", "l");
    opt.legend->AddEntry(opt.fom_graph.get(), "Figure of Merit", "l");
    opt.purity_efficiency_fom->Add(opt.purity_graph.get());
    opt.purity_efficiency_fom->Add(opt.efficiency_graph.get());
    opt.purity_efficiency_fom->Add(opt.fom_graph.get());
    opt.purity_efficiency_fom->GetXaxis()->SetTitle(("Upper Bound on " + var_to_title.at(var)).c_str());
    opt.purity_efficiency_fom->GetYaxis()->SetTitle("Purity, Efficiency");
    opt.purity_efficiency_fom->Draw("al");
    opt.purity_efficiency_fom->GetHistogram()->GetXaxis()->SetRangeUser(lw, up);
    opt.purity_efficiency_fom->SetMinimum(0);
    opt.purity_efficiency_fom->SetMaximum(1);
    opt.canvas->Update();
    opt.line = std::make_shared<TLine>(opt.limit, 0, opt.limit, 1);
    opt.line->SetBit(kCanDelete, false);
    opt.line->SetLineColor(kRed);
    opt.line->SetLineStyle(kDashed);
    opt.legend->Draw("same");
    opt.line->Draw("same");
    double bin_height = 0;
    double arrow_length = 0;
    for (TObject* obj : *static_cast<TList*>(lc.stack->GetHists()))
    {
      TH1* hist = static_cast<TH1*>(obj);
      int bin = hist->FindBin(opt.limit);
      bin_height += hist->GetBinContent(bin);
      if (arrow_length == 0)
        arrow_length = 5*(hist->GetXaxis()->GetTickLength());
    }
    lc.arrow = std::make_shared<TArrow>(opt.limit, bin_height,
                                        opt.limit - arrow_length, bin_height, 0.05, "|->");
    lc.arrow->SetBit(kCanDelete, false);
    lc.arrow->SetLineColor(kRed);
    lc.arrow->SetLineWidth(2);
    lc.arrow->SetLineStyle(kDashed);
    opt.gaxis = std::make_shared<TGaxis>(up, 0, up, 1, 0, std::ceil(opt.fom), 510, "+L");
    opt.gaxis->SetBit(kCanDelete, false);
    opt.gaxis->SetTitle("Figure of Merit");
    opt.gaxis->SetTitleSize(opt.purity_efficiency_fom->GetHistogram()->GetYaxis()->GetTitleSize());
    opt.gaxis->SetTitleOffset(0.8);
    opt.gaxis->SetTitleFont(opt.purity_efficiency_fom->GetHistogram()->GetYaxis()->GetTitleFont());
    opt.gaxis->SetLabelSize(opt.purity_efficiency_fom->GetHistogram()->GetYaxis()->GetLabelSize());
    opt.gaxis->SetLabelOffset(opt.purity_efficiency_fom->GetHistogram()->GetYaxis()->GetLabelOffset());
    opt.gaxis->SetLabelFont(opt.purity_efficiency_fom->GetHistogram()->GetYaxis()->GetLabelFont());
    opt.gaxis->SetTickLength(opt.purity_efficiency_fom->GetHistogram()->GetYaxis()->GetTickLength());
    opt.gaxis->SetNdivisions(opt.purity_efficiency_fom->GetHistogram()->GetYaxis()->GetNdivisions());
    opt.gaxis->Draw();
    lc.canvas->cd();
    lc.arrow->Draw();
    return std::make_tuple(opt, lc);
  }

  /**
   * @brief Optimize either an upper or lower bound
   */
  std::tuple<ana::tools::opt_canvas, ana::tools::limit_canvas> analysis_tree::optimize_bound(const std::string& var,
                                                                                             const double& lw, const double& up,
                                                                                             const ana::tools::cut_sequence& old_cut,
                                                                                             const bool& is_upper_bound) const
  {
    return (is_upper_bound) ? optimize_upper_bound(var, lw, up, old_cut)
                            : optimize_lower_bound(var, lw, up, old_cut);
  }
  
}// end ana::tools namespace
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TColor.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TError.h"
#include "TMarker.h"
#include "TLine.h"
#include "TMemFile.h"
#include "TEntryList.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGaxis.h"
#include "TArrow.h"

#include <memory>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <functional>
#include <unordered_map>
#include <tuple>

#ifndef SIG_BKG_H
#define SIG_BKG_H

namespace nc
{
  struct stack_canvas
  {
    std::shared_ptr<TCanvas> canvas;
    std::shared_ptr<TLegend> legend;
    std::shared_ptr<THStack> stack;
    std::vector<std::shared_ptr<TH1F>> hists;
  };

  struct limit_canvas : stack_canvas
  {
    limit_canvas(stack_canvas&& base)
    {
      sc = base;
      canvas = sc.canvas;
      legend = sc.legend;
      stack = sc.stack;
      hists = sc.hists;
    }
    stack_canvas sc;
    std::shared_ptr<TCanvas> canvas;
    std::shared_ptr<TLegend> legend;
    std::shared_ptr<THStack> stack;
    std::vector<std::shared_ptr<TH1F>> hists;
    std::shared_ptr<TArrow> arrow;
  };

  struct optimization
  {
    std::shared_ptr<TCanvas> canvas;
    std::shared_ptr<TGaxis> gaxis;
    std::shared_ptr<TMultiGraph> purity_efficiency_fom;
    std::shared_ptr<TLegend> legend;
    std::shared_ptr<TGraph> purity_graph;
    std::shared_ptr<TGraph> efficiency_graph;
    std::shared_ptr<TGraph> fom_graph;
    std::shared_ptr<TLine> line;
    std::string cut;
    double limit;
    double fom;
    bool is_upper_bound; // false -> lower bound
  };

  class analysis_tree
  {
    public:
      analysis_tree(const std::string& inFileName,
                    const std::string& directory,
                    const std::string& selection_cut,
                    const std::string& sel_tree,
                    const std::string& sig_tree,
                    const std::string& cos_tree,
                    const std::vector<std::pair<std::string, std::string>> sel_cats,
                    const std::vector<std::pair<std::string, std::string>> sig_cats) : inFile(TFile::Open(inFileName.c_str(), "READ"))
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
        colors.push_back(new TColor(0.45, 0.45, 0.45)); // Gray (just gray)
        colors.push_back(new TColor(0.85, 0.85, 0.85)); // Gray (but darker)

        // Memory Handling
        std::string tmpName = "temp_TTrees-"+sel_tree+"-"+sig_tree+"-"+cos_tree+"_Cuts-"+sel_cats.front().second+"-"+sig_cats.front().second+".root";
        memFile = std::make_unique<TMemFile>(tmpName.c_str(), "RECREATE");
        memFile->cd();

        // Read In Trees
        std::cout << "Reading from " << inFile->GetName() << std::endl;
        std::string selection_tree_name = directory + "/" + sel_tree;
        std::string signal_tree_name    = directory + "/" + sig_tree;
        std::string cosmic_tree_name    = directory + "/" + cos_tree;
        std::string signal_cut    = sel_cats.front().second;
        std::string bkg_cut       = get_bkg_cut(signal_cut); 
        TTree* tmpSignalTree = inFile->Get<TTree>(signal_tree_name.c_str());
        signalTree = std::unique_ptr<TTree>(tmpSignalTree->CopyTree(signal_cut.c_str()));
        signalTree->SetDirectory(memFile.get());
        TTree* selectedTree = inFile->Get<TTree>(selection_tree_name.c_str());
        selectedSignalTree = std::unique_ptr<TTree>(signalTree->CopyTree(selection_cut.c_str()));
        selectedSignalTree->SetDirectory(memFile.get());
        nuBackgroundTree   = std::unique_ptr<TTree>(selectedTree->CopyTree(bkg_cut.c_str()));
        nuBackgroundTree->SetDirectory(memFile.get());
        TTree* tmpCosmicTree = inFile->Get<TTree>(cosmic_tree_name.c_str());
        tmpCosmicTree->LoadBaskets();
        tmpCosmicTree->GetEntry(1);
        cosmicTree = std::unique_ptr<TTree>(tmpCosmicTree->CopyTree("1"));
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

        //memFile->Write();
      }

      void add_variable(const std::string& name,
                        const unsigned int& bins,
                        const double& lw, const double& up,
                        const std::string& title)
       {
         if (var_to_range.count(name) != 0 ||
             var_to_title.count(name) != 0)
           throw std::runtime_error("analysis_tree::add_variable -- Variable " + name + " already added to tree");
         var_to_range.emplace(name, std::make_tuple(bins, lw, up));
         var_to_title.emplace(name, title);
       }

      std::string get_bkg_cut(const std::string& signal_cut) const
      {
        std::string bkg_cut = "(! ("+signal_cut+"))";
        return bkg_cut;
      }

      std::tuple<unsigned int, unsigned int, double, double> cut_e_p(const std::string& cut) const
      {
        unsigned int selected_signal      = selectedSignalTree->GetEntries(cut.c_str());
        unsigned int selected_background  = nuBackgroundTree  ->GetEntries(cut.c_str());
                     selected_background += cosmicTree        ->GetEntries(cut.c_str());
        double efficiency = static_cast<double>(selected_signal) / total_signal;
        double purity     = static_cast<double>(selected_signal) / (selected_signal + selected_background);
        return std::make_tuple(selected_signal, selected_background, efficiency, purity);
      }

      void report_on_cut(const std::string& cut) const
      {
        auto [selected_signal, selected_background, efficiency, purity] = cut_e_p(cut);
        std::vector<unsigned int> bkg_cat_counts;
        size_t nCats = sel_cat_labels.size();
        std::ostringstream cat_report;
        for (size_t cat = 0; cat < nCats; ++cat)
        {
          std::string cat_cut = "("+cut+") && ("+sel_cat_cuts.at(cat)+")";
          bkg_cat_counts.push_back(nuBackgroundTree->GetEntries(cat_cut.c_str()));
          if (cat == nCats - 1)
            bkg_cat_counts.back() += cosmicTree->GetEntries(cut.c_str());
          cat_report << "    Cat " << cat << " (" << sel_cat_labels.at(cat) << "): " << bkg_cat_counts.at(cat);
          if (cat != nCats - 1)
            cat_report << '\n';
        }
        std::cout
          << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << '\n' 
          << cut                                                     << '\n'
          << "  Total Signal: "               << total_signal        << '\n'
          << "  Selected Signal Events: "     << selected_signal     << '\n'
          << "  Selected Background Events: " << selected_background << '\n'
          <<    cat_report.str()                                     << '\n'
          << "  efficiency: "                 << efficiency          << '\n'
          << "  purity: "                     << purity              << '\n'
          << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      }

      stack_canvas plot_var_sel(const std::string& var, const std::string& cut) const
      {
        // check the var is in each of the TTrees
        auto check_branch = [](const std::unique_ptr<TTree>& tree, const std::string var)
        {
          if (not tree->GetBranch(var.c_str()) &&
              not tree->GetLeaf  (var.c_str())  )
              throw std::runtime_error("Tree "+std::string(tree->GetName())+" has no "+var+" branch. Cannot plot it.");
        };
        check_branch(selectedSignalTree, var);
        check_branch(nuBackgroundTree,   var);
        check_branch(cosmicTree,         var);

        stack_canvas sc;
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
          std::string histName = "nu_bkg_" + std::to_string(cat);
          std::string cut_cat = cut + " && " + sel_cat_cuts.at(cat);
          std::shared_ptr<TH1F> tmp_hist = std::make_shared<TH1F>(histName.c_str(), "", bins, lw, up);
          nuBackgroundTree->Draw((var+">>"+histName).c_str(), cut_cat.c_str(), "goff");
          tmp_hist->SetDirectory(nullptr);
          bkg_hists.push_back(tmp_hist);
          if (bkg_hists.at(cat) == nullptr)
            throw std::runtime_error("bkg_hists for cat " + std::to_string(cat) + " is null");
        }
        std::shared_ptr<TH1F> cos_bkg = std::make_shared<TH1F>("cos_bkg", "", bins, lw, up);
        cosmicTree->Draw((var+">>cos_bkg").c_str(), cut.c_str(), "goff");
        cos_bkg->SetDirectory(nullptr);
        cos_bkg->Add(bkg_hists.at(nCats - 1).get());
        for (size_t cat = 1; cat < nCats - 1; ++cat)
          sc.hists.push_back(std::move(bkg_hists.at(cat)));
        sc.hists.push_back(std::move(cos_bkg));
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

      stack_canvas plot_var_sig(const std::string& var, const std::string& cut) const
      {
        // check the var is in each of the TTrees
        auto check_branch = [](const std::unique_ptr<TTree>& tree, const std::string var)
        {
          if (not tree->GetBranch(var.c_str()) &&
              not tree->GetLeaf  (var.c_str())  )
              throw std::runtime_error("Tree "+std::string(tree->GetName())+" has no "+var+" branch. Cannot plot it.");
        };
        check_branch(signalTree, var);

        stack_canvas sc;
        sc.canvas = std::make_shared<TCanvas>(var.c_str(), var.c_str(), 1618, 1000);
        sc.canvas->SetBit(kCanDelete, false);
        sc.canvas->cd();

        auto [bins, lw, up] = var_to_range.at(var); 
        size_t nCats = sig_cat_labels.size();
        for (size_t cat = 0; cat < nCats; ++cat)
        {
          std::string histName = "sel_sig_" + std::to_string(cat);
          std::string cut_cat = cut + " && " + sig_cat_cuts.at(cat);
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
          std::string cut_cat = "(! ("+cut+")) && " + sig_cat_cuts.at(cat);
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

      double get_fom(const unsigned int& selected_signal,
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

      std::tuple<optimization, limit_canvas> optimize_lower_bound(const std::string& var,
                                                                  const double& lw, const double& up,
                                                                  const std::string& old_cut) const
      {
        optimization opt;
        limit_canvas lc = plot_var_sel(var, old_cut);
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
        opt.cut = "";
        opt.limit = std::numeric_limits<double>::lowest();
        for (double limit = lw; limit < up; limit += step)
        {
          std::ostringstream limStrm;
          limStrm << std::setprecision(10) << limit;
          std::string cut = old_cut + " && (" + var + " > " + limStrm.str() + ")";
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

      std::tuple<optimization, limit_canvas> optimize_upper_bound(const std::string& var,
                                                                  const double& lw, const double& up,
                                                                  const std::string& old_cut) const
      {
        optimization opt;
        limit_canvas lc = plot_var_sel(var, old_cut);
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
        opt.is_upper_bound = true;
        unsigned int nSteps = 1000;
        double step = (up - lw) / nSteps;
        opt.fom = std::numeric_limits<double>::lowest();
        opt.cut = "";
        opt.limit = std::numeric_limits<double>::lowest();
        for (double limit = lw; limit < up; limit += step)
        {
          std::ostringstream limStrm;
          limStrm << std::setprecision(10) << limit;
          std::string cut = old_cut + " && (" + var + " < " + limStrm.str() + ")";
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
    private:
      std::unique_ptr<TMemFile> memFile;
      std::unique_ptr<TFile> inFile;
      std::vector<TColor*> colors;
      std::unique_ptr<TTree> signalTree;
      std::unique_ptr<TTree> selectedSignalTree;
      std::unique_ptr<TTree> nuBackgroundTree;
      std::unique_ptr<TTree> cosmicTree;
      unsigned int total_signal;
      std::vector<std::string> sel_cat_labels;
      std::vector<std::string> sel_cat_cuts;
      std::vector<std::string> sig_cat_labels;
      std::vector<std::string> sig_cat_cuts;
      std::unordered_map<std::string, std::tuple<unsigned int, double, double>> var_to_range;
      std::unordered_map<std::string, std::string> var_to_title;
  };
}

#endif // SIG_BKG_H

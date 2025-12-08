/**
 * @file plot_structs.cc
 * @brief Functions for the plot structs
 * @author hhausner@fnal.gov
 */

#include "include/plot_structs.h"

#include <iostream>

/**
 * @namespace ana::tools
 * @brief Namespace for tools helpful to running analyses off SPINE outputs
 */
namespace ana::tools
{
  /**
   * @brief make a TH1F with the sum of the sections
   **/
  void stack_canvas::SumHist()
  {
    sum = std::make_shared<TH1F>("sum", "Sum", hists.front()->GetNbinsX(), stack->GetXaxis()->GetXmin(), stack->GetXaxis()->GetXmax());
    for (auto const& hist : hists)
      sum->Add(hist.get());
  }

  /**
   * @brief Find where to put the legend
   **/
  std::tuple<double, double, double, double> stack_canvas::LegendPlace() const
  {
    canvas->cd();
    stack->Draw("hist");

    canvas->GetBottomMargin();

    double x1 = 0.6;
    double y1 = 0.45;
    double x2 = 0.85;
    double y2 = 0.875;

    double xMin = stack->GetXaxis()->GetXmin();
    double yMin = stack->GetMinimum();
    double xMax = stack->GetXaxis()->GetXmax();
    double yMax = stack->GetMaximum();

    double yaxis_low = (y1 - canvas->GetBottomMargin()) / (1 - canvas->GetTopMargin() - canvas->GetBottomMargin()) * (yMax - yMin) + yMin;

    size_t nBins = hists.front()->GetNbinsX();
    std::vector<double> bin_heights(nBins, 0.0);
    for (size_t bin = 1; bin < nBins + 1; ++bin)
    {
      for (const auto& hist : hists)
        bin_heights.at(bin - 1) += hist->GetBinContent(bin);
    }

    double x_step = 0.01;

    do
    {
      double xaxis_low  = (x1 - canvas->GetLeftMargin()) / (1 - canvas->GetRightMargin() - canvas->GetLeftMargin()) * (xMax - xMin) + xMin;
      double xaxis_high = (x2 - canvas->GetLeftMargin()) / (1 - canvas->GetRightMargin() - canvas->GetLeftMargin()) * (xMax - xMin) + xMin;
      size_t first_bin = hists.front()->GetXaxis()->FindBin(xaxis_low);
      size_t last_bin  = hists.front()->GetXaxis()->FindBin(xaxis_high);
      double bins_max = std::numeric_limits<double>::lowest();
      for (size_t bin = first_bin; bin < last_bin + 1; ++bin)
        bins_max = (bins_max < bin_heights.at(bin - 1)) ? bin_heights.at(bin - 1) : bins_max; 
      if (bins_max < yaxis_low - 0.1)
        break;
      x1 -= x_step;
      x2 -= x_step;
    }
    while (x1 > canvas->GetLeftMargin());

    return std::tie(x1, y1, x2, y2);
  }

  /**
   * @brief Add Preliminary to canvas
   **/
  void stack_canvas::Preliminary()
  {
    double margin = canvas->GetLeftMargin();
    TLatex* prelimText = new TLatex();
    prelimText->SetNDC();
    prelimText->SetTextColor(kRed);
    prelimText->SetTextAlign(12);
    canvas->cd();
    prelimText->DrawLatex(margin, 0.95, "ICARUS Preliminary");
  }
 
  /**
   * @brief Draw statistical uncertainty
   **/
  void stack_canvas::SetStatUncertainty()
  {
    err = std::make_shared<TH1F>(*static_cast<TH1F*>(hists.front()->Clone("ErrHist")));

    for (size_t bin = 1; bin < err->GetNbinsX() + 1; ++bin)
    {
      double binContent = 0;
      for (const auto& hist : hists)
        binContent += hist->GetBinContent(bin);
      err->SetBinContent(bin, binContent);
    }

    err->SetFillStyle(1001);
    err->SetFillColorAlpha(12, 0.55);
    err->SetLineStyle(0);
    err->SetLineWidth(0);
    err->SetMarkerStyle(0);
  }

  /**
   * @brief Construct the legend
   **/
  void stack_canvas::ConstructLegend()
  {
    auto const [x1, y1, x2, y2] = LegendPlace();
    legend = std::make_shared<TLegend>(x1, y1, x2, y2); 
    legend->SetBit(kCanDelete, false);
    legend->SetFillStyle(0);
  }
   
  /**
   * @brief Print the canvas
   * @details ROOT doesn't play nicely with std::strings, so this wraps the TCanvas::Print method
   * so we can print our canvases without having to call c_str all the dang time
   */
  void stack_canvas::Print(const std::string& file_name) const
  {
    canvas->Update();
    canvas->Print(file_name.c_str());
  }

  /**
   * @brief Print with Preliminary label
   **/
  void stack_canvas::PrintPreliminary(const std::string& file_name)
  {
    Preliminary();
    SetStatUncertainty();
    legend->AddEntry(err.get(), "Statistical Uncertainty", "f");
    err->Draw("E2,same");
    Print(file_name);
  }

  /**
   * @brief Find where to put the legend
   **/
  std::tuple<double, double, double, double> opt_canvas::LegendPlace() const
  {
    double x1 = 0.6;
    double y1 = 0.4;
    double x2 = 0.85;
    double y2 = 0.85;
    return std::tie(x1, y1, x2, y2);
  }
  /**
   * @brief Add Preliminary to canvas
   **/
  void opt_canvas::Preliminary()
  {
    double margin = canvas->GetLeftMargin();
    TLatex* prelimText = new TLatex();
    prelimText->SetNDC();
    prelimText->SetTextColor(kRed);
    prelimText->SetTextAlign(12);
    canvas->cd();
    prelimText->DrawLatex(margin, 0.95, "ICARUS Preliminary");
  }

  /**
   * @brief Print the canvas
   * @details ROOT doesn't play nicely with std::strings, so this wraps the TCanvas::Print method
   * so we can print our canvases without having to call c_str all the dang time
   */
  void opt_canvas::Print(const std::string& file_name) const
  {
    canvas->Print(file_name.c_str());
  }

  /**
   * @brief Print with Preliminary label
   **/
  void opt_canvas::PrintPreliminary(const std::string& file_name)
  {
    canvas->Update();
    Preliminary();
    Print(file_name);
  }

  /**
   * @brief Construct a limit_canvas from a stack_canvas
   * @details The ordering is important
   */
  limit_canvas::limit_canvas(stack_canvas&& base)
  {
    sc = base;
    canvas = sc.canvas;
    legend = sc.legend;
    stack = sc.stack;
    hists = sc.hists;
  }

}// end ana::tools namespace
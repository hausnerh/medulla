/**
 * @file plot_structs.h
 * @brief Header file for structs useful for plotting analyses
 * @details Defines three structs, stack_canvas, limit_canvas, and opt_canvas, used different plotting ideas.
 * The stack_canvas holds a canvas and a stack histogram for plotting a selection with different categories.
 * The limit_canvas allows you to also show a cut on a specific value. The opt_canvas shows
 * an optimization procedure done to fine a specific cut value.
 * @author hhausner@fnal.gov
 */

#ifndef PLOT_STRUCTS_H
#define PLOT_STRUCTS_H

#include "cut_sequence.h"

#include "TArrow.h"
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMarker.h"
#include "TMultiGraph.h"

#include <memory>
#include <string>

/**
 * @namespace ana::tools
 * @brief Namespace for tools helpful to running analyses off SPINE outputs
 */
namespace ana::tools
{
  /**
   * @struct stack_canvas
   * @brief A structure for holding and plotting a list of histograms we want to draw in a stack
   * @details We need to carry around the TH1Fs in order for the THStack to work right.
   */
  struct stack_canvas
  {
    std::shared_ptr<TCanvas> canvas;
    std::shared_ptr<TLegend> legend;
    std::shared_ptr<THStack> stack;
    std::vector<std::shared_ptr<TH1F>> hists;
    void Print(const std::string& file_name) const;
  };

  /**
   * @struct limit_canvas
   * @brief A structure for holding and plotting a list of histograms to draw in a stack plus an arrow to demark a cut
   * @details Inherit from stack_canvas
   */
  struct limit_canvas : stack_canvas
  {
    limit_canvas(stack_canvas&& base);
    stack_canvas sc;
    std::shared_ptr<TCanvas> canvas;
    std::shared_ptr<TLegend> legend;
    std::shared_ptr<THStack> stack;
    std::vector<std::shared_ptr<TH1F>> hists;
    std::shared_ptr<TArrow> arrow;
  };

  /**
   * @struct opt_canvas
   * @brief A structure for plotting an optimization, with a purity, efficiency, and figure of merit
   */
  struct opt_canvas
  {
    std::shared_ptr<TCanvas> canvas;
    std::shared_ptr<TGaxis> gaxis;
    std::shared_ptr<TMultiGraph> purity_efficiency_fom;
    std::shared_ptr<TLegend> legend;
    std::shared_ptr<TGraph> purity_graph;
    std::shared_ptr<TGraph> efficiency_graph;
    std::shared_ptr<TGraph> fom_graph;
    std::shared_ptr<TLine> line;
    ana::tools::cut_sequence cut;
    ana::tools::cut_sequence truth_cut; ///< unused in non-threshold optimizations
    double limit;
    double fom;
    bool is_upper_bound; ///< false -> lower bound
    void Print(const std::string& file_name) const;
  };

}// end ana::tools namespace

#endif
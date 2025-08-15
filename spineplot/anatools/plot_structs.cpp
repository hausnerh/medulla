/**
 * @file plot_structs.cpp
 * @brief Functions for the plot structs
 * @author hhausner@fnal.gov
 */

#include "plot_structs.h"

/**
 * @namespace ana::tools
 * @brief Namespace for tools helpful to running analyses off SPINE outputs
 */
namespace ana::tools
{
  /**
   * @brief Print the canvas
   * @details ROOT doesn't play nicely with std::strings, so this wraps the TCanvas::Print method
   * so we can print our canvases without having to call c_str all the dang time
   */
  void stack_canvas::Print(const std::string& file_name) const
  {
    canvas->Print(file_name.c_str());
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
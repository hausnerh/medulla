
/**
 * @file analysis_tree.h
 * @brief Header file for the analysis_tree class, which wrap the TTrees produced in SPINE analyses
 * @details The idea is to put all of the ROOT TTree interactions under the hood so analyzers can focus on plotting
 * and applying cuts.
 * @author hhausner@fnal.gov
 */

#ifndef ANALYSIS_TREE_H
#define ANALYSIS_TREE_H

#include "TColor.h"
#include "TFile.h"
#include "TMemFile.h"
#include "TTree.h"

#include <tuple>
#include <unordered_map>

#include "cut_sequence.h"
#include "plot_structs.h"

/**
 * @namespace ana::tools
 * @brief Namespace for tools helpful to running analyses off SPINE outputs
 */
namespace ana::tools
{
  /**
   * @class analysis_tree
   * @brief This class holds the signal and background events in your selection and allows you to apply cuts and plot variables
   * @details Specify the ROOT file which holds the TTrees, as well as where to find the signal/selected nu/selected cosmic TTrees,
   * which should be the standard outputs of a SPINE analysis. Also supply the signal and selection catagories for plotting.
   */
  class analysis_tree
  {
    public:
      // constructor
      analysis_tree(const std::string& inFileName,
                    const std::string& directory,
                    const ana::tools::cut_sequence& selection_cut,
                    const std::string& sel_tree,
                    const std::string& sig_tree,
                    const std::string& cos_tree,
                    const std::vector<std::pair<std::string, ana::tools::cut_sequence>> sel_cats,
                    const std::vector<std::pair<std::string, ana::tools::cut_sequence>> sig_cats);

      // Modifiers
      void add_variable(const std::string& name,
                        const unsigned int& bins,
                        const double& lw, const double& up,
                        const std::string& title);
      std::tuple<ana::tools::opt_canvas, ana::tools::limit_canvas> optimize_threshold(const std::string& var,
                                                                                      const std::string& true_var,
                                                                                      const double& lw, const double& up,
                                                                                      const ana::tools::cut_sequence& old_cut,
                                                                                      const ana::tools::cut_sequence& old_truth_cut);
      // Const Methods
      std::vector<std::string> variables() const;
      std::tuple<unsigned int, unsigned int, unsigned int, double, double> signal_def_e_p(const ana::tools::cut_sequence& cut,
                                                                                          const ana::tools::cut_sequence& truth_cut) const;
      std::tuple<unsigned int, unsigned int, double, double> cut_e_p(const ana::tools::cut_sequence& cut) const;
      void report_on_cut(const ana::tools::cut_sequence& cut) const;
      ana::tools::stack_canvas plot_var_sel(const std::string& var, const ana::tools::cut_sequence& cut) const;
      ana::tools::stack_canvas plot_var_sig(const std::string& var, const ana::tools::cut_sequence& cut) const;
      double get_fom(const unsigned int& selected_signal,
                     const unsigned int& selected_background,
                     const double& efficiency, 
                     const double& purity) const;
      std::tuple<ana::tools::opt_canvas, ana::tools::limit_canvas> optimize_lower_bound(const std::string& var,
                                                                                        const double& lw, const double& up,
                                                                                        const ana::tools::cut_sequence& old_cut) const;
      std::tuple<ana::tools::opt_canvas, ana::tools::limit_canvas> optimize_upper_bound(const std::string& var,
                                                                                        const double& lw, const double& up,
                                                                                        const ana::tools::cut_sequence& old_cut) const;
      std::tuple<ana::tools::opt_canvas, ana::tools::limit_canvas> optimize_bound(const std::string& var,
                                                                                  const double& lw, const double& up,
                                                                                  const ana::tools::cut_sequence& old_cut,
                                                                                  const bool& is_upper_bound) const; 

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
      std::vector<std::string> sig_cat_labels;
      std::vector<ana::tools::cut_sequence> sel_cat_cuts;
      std::vector<ana::tools::cut_sequence> sig_cat_cuts;
      std::vector<std::string> vars;
      std::unordered_map<std::string, std::tuple<unsigned int, double, double>> var_to_range;
      std::unordered_map<std::string, std::string> var_to_title;
  };

}// end ana::tools

#endif
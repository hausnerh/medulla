/**
 * @file cut_sequence.h
 * @brief Header file for Cut Sequences to streamline analysis with SPINE analysis trees
 * @details The cut_sequence makes handling cuts in ways ROOT TTrees like easier without while letting
 * the user write simple cuts that can be chained together. It also handles conditional cuts and negating cuts.
 * @author hhausner@fnal.gov
 */

#ifndef CUT_SEQ_H
#define CUT_SEQ_H

#include <initializer_list>
#include <string>
#include <vector>

/**
 * @namespace ana::tools
 * @brief Namespace for tools helpful to running analyses off SPINE outputs
 */
namespace ana::tools
{
  /**
   * @class cut_sequence
   * @brief Oraganizes analysis cuts for human-readable inputs and ROOT TTree-interpretable outputs
   * @details Store each stage of the cut as a std::string and chain them together in ways a TTree
   * can understand. Be able to add and modify your cuts in logical ways and not have to manually
   * handle the strings yourself
   */
  class cut_sequence
  {
    public:
      // constructors
      cut_sequence();
      cut_sequence(const std::vector<std::string>& cut_list);
      cut_sequence(const std::initializer_list<std::string>& init);
      cut_sequence(const std::string& initial_cut);
      cut_sequence(const char* cstr);

      // non-const methods
      void add_cut(const std::string& cut);
      void add_conditional_cut(const std::string& cut_if, const std::string& cut_then);
      void add_sequence(const cut_sequence& new_seq);

      // const methods
      size_t size() const;
      std::string get_cut(const size_t& index) const;
      std::string string() const;
      const char* c_str() const;
      cut_sequence with_addition(const cut_sequence& additional_cut) const;

      // friends
      friend cut_sequence& operator+=(cut_sequence& lhs, const cut_sequence& rhs);
  
    private:
      std::vector<std::string> cuts;
      std::string cut_as_string;
  };

  cut_sequence& operator+=(cut_sequence& lhs, const cut_sequence& rhs);
  cut_sequence negate_cut(const cut_sequence& old_cut);
}// end ana::tools namespace

#endif
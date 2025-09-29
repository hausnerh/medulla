/**
 * @file cut_sequence.cc
 * @brief Functions for the cut_sequence class
 * @author hhausner@fnal.gov
 */

#include "include/cut_sequence.h"

/**
 * @namespace ana::tools
 * @brief Namespace for tools helpful to running analyses off SPINE outputs
 */
namespace ana::tools
{
  /**
   * @brief Default constructor for cut_sequence
   */
  cut_sequence::cut_sequence() {}

  /**
   * @brief Copy constructor
   */
  cut_sequence::cut_sequence(const cut_sequence& other) : cuts(other.cuts)
  {
    cut_as_string = string();
  }

  /**
   * @brief Construct a cut sequence from a vector of individual cut strings
   * @details The strings should be of a form such that a TTree could interpret them, eg "vertex_x < 7"
   */
  cut_sequence::cut_sequence(const std::vector<std::string>& cut_list) : cuts(cut_list)
  {
    cut_as_string = string();
  }

  /**
   * @brief Construct a cut sequence from an initializer list
   * @details The strings should be of a form such that a TTree could interpret them, eg "vertex_x < 7"
   */
  cut_sequence::cut_sequence(const std::initializer_list<std::string>& init) : cuts(init)
  {
    cut_as_string = string();
  }

  /**
   * @brief Construct a cut sequence from a single initial cut
   * @details The string should be of the form that a TTree could interpret, eg "vertex_x < 7"
   */
  cut_sequence::cut_sequence(const std::string& initial_cut) : cuts({initial_cut})
  {
    cut_as_string = string();
  }
 
  /**
   * @brief Construct a cut sequence from a single initial cut
   * @details The c-string should be of the form that a TTree could interpret, eg "vertex_x < 7"
   */
  cut_sequence::cut_sequence(const char* cstr) : cuts({cstr})
  {
    cut_as_string = string();
  }

  /**
   * @brief Add a cut to a cut sequence
   * @details The cut should be a string of the form that a TTree could interpret, eg "vertex_x < 7"
   */
  void cut_sequence::add_cut(const std::string& cut)
  {
    cuts.push_back(cut);
    cut_as_string = string();
  }

  /**
   * @brief Add a conditional cut to the sequence
   * @details Conditionals (P -> Q) are equivalent to (!P or Q), which is easier for the TTree to parse
   * The cut and condition should be a strings of a form that a TTree could interpret, eg "vertex_x < 7"
   */
  void cut_sequence::add_conditional_cut(const std::string& cut_if, const std::string& cut_then)
  {
    std::string cut = "(!(" + cut_if +")) || (" + cut_then+ ")";
    add_cut(cut);
  }

  /**
   * @brief Append another cut sequence to the current one
   */
  void cut_sequence::add_sequence(const cut_sequence& new_seq)
  {
    if (new_seq.size() != 0)
      for (size_t idx = 0; idx < new_seq.size(); ++idx)
        add_cut(new_seq.get_cut(idx));
  }

  /**
   * @brief Assignment operator
   **/
  cut_sequence& cut_sequence::operator=(const cut_sequence& other)
  {
    if (this != &other)
    {
      cuts = other.cuts;
      cut_as_string = string();
    }
    return *this;
  }

  /**
   * @brief How many cuts are in the sequence?
   */
  size_t cut_sequence::size() const
  {
    return cuts.size();
  }

  /**
   * @brief Return the string for the Nth cut
   */
  std::string cut_sequence::get_cut(const size_t& index) const
  {
    return cuts.at(index);
  }

  /**
   * @brief Return a string corresponding to the whole cut sequence
   * @details This is mostly for debugging/print statements. To apply your cuts to a TTree
   * see cut_sequence::c_str
   */
  std::string cut_sequence::string() const
  {
    std::string full_cut = "";
    for (size_t idx = 0; idx < cuts.size(); ++idx)
    {
      full_cut += "(" + cuts.at(idx) + ")";
      if (idx != cuts.size() - 1)
        full_cut += " && ";
    }
    return full_cut;
  }

  /**
   * @brief Return a c-string correstponding to the whole cut sequence
   * @details TTrees only like c-strings, not std::strings, so this conversion method exists to
   * make it easier to apply the cuts in, for example TTree::GetEntries
   */
  const char* cut_sequence::c_str() const
  {
    static thread_local std::string tmp_str;
    tmp_str = string();
    return tmp_str.c_str();
  }

  /**
   * @brief Return a cut sequence with with an additional cut sequence appended
   * @details This is mostly useful for chaining cuts together without having to declare a new one
   */
  cut_sequence cut_sequence::with_addition(const cut_sequence& additional_cut) const
  {
    cut_sequence new_cut(cuts);
    new_cut.add_sequence(additional_cut);
    return new_cut;
  }

  /**
   * @brief An update opperator to make appending a new cut easier in scripts
   */
  cut_sequence& operator+=(cut_sequence& lhs, const cut_sequence& rhs)
  {
    lhs.add_sequence(rhs);
    return lhs;
  }

  /**
   * @brief Negate a cut sequence
   * @details This cut will return false when the input cut was true. Note that the returned cut is of size 1,
   * that is to say, the negated cut is one string. This negated cut can then be appended to other cut sequences.
   * Negation is a simple application of de Morgan's Law
   */
  cut_sequence negate_cut(const cut_sequence& old_cut)
  {
    std::string neg_str = "";
    for (size_t idx = 0; idx < old_cut.size(); ++idx)
    {
      neg_str += "(!(" + old_cut.get_cut(idx) + "))";
      if (idx != old_cut.size() - 1)
        neg_str += " || ";
    }
    cut_sequence neg_seq(neg_str);
    return neg_seq;
  }

}// end namespace ana::tools
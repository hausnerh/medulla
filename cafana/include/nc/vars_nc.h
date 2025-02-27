/**
 * @file vars_nc.h
 * @brief Header file for defining variables useful for NC analyses
 * @details This file defines some helpful counting functions
 * @author hhausner@fnal.gov
 **/

#ifndef VARS_NC_H
#define VARS_NC_H
#include "include/utilities.h"
#include <regex>

/**
 * @namespace vars::nc
 * @brief This namespace is for organizing variables and functions
 * useful for NC analyses.
 **/

namespace vars::nc
{
  // ****** UTILITY FUNCTIONS ******

  /**
   * @brief Count the number of particles in the topology string
   * @details search topology string for the particle specified by the Particle_t enum
   * and count the number of primaries. We construct the regex string based on the pvars::Particle_t
   * (and throw an error if we cannot match it). The regex looks for a number followed by the particle string
   * followed by another number or the end of the string. Pass the string to the golfer and see if there is a match.
   * If there is no match assume there are no primaries of that type.
   * @param topology the string containing an interaction topology
   * @param pid the pvars::Particle_t particle enumeration we are trying to count
   * @return the number of primaries of type pid in the topology string
   **/
  unsigned int count_topology_particle(std::string topology, int pid)
  {
    std::string particleStr = (pid == pvars::kPhoton)   ? "g"  :
                              (pid == pvars::kElectron) ? "e"  :
                              (pid == pvars::kMuon)     ? "mu" :
                              (pid == pvars::kPion)     ? "pi" :
                              (pid == pvars::kProton)   ? "p"  :
                                                           ""  ;
    if (particleStr == "")
      throw std::invalid_argument("in vars::nc::count_particles — cannot parse particle "+std::to_string(pid));
    std::string searchStr = R"((\d+))" + particleStr + R"((\d+|$))";
    std::regex golfer(searchStr);
    std::smatch matches;
    if (std::regex_search(topology, matches, golfer))
    {
      return std::stoul(matches[1].str());
    }
    else
    {
      return 0;
    }
  }

  /**
   * @brief count photons in the topology string for an interaction
   * @details call count_topology_particle with pid = pvars::kPhoton
   * @param obj the interaction of interest (MC or data)
   * @return the number of photons counted in the topology
   **/
  template<class T>
    unsigned int count_photons(const T& obj)
    {
      return count_topology_particle(obj.topology, pvars::kPhoton);
    }

  /**
   * @brief count electrons in the topology string for an interaction
   * @details call count_topology_particle with pid = pvars::kElectron
   * @param obj the interaction of interest (MC or data)
   * @return the number of electrons counted in the topology
   **/
  template<class T>
    unsigned int count_electrons(const T& obj)
    {
      return count_topology_particle(obj.topology, pvars::kElectron);
    }

  /**
   * @brief count muons in the topology string for an interaction
   * @details call count_topology_particle with pid = pvars::kMuon
   * @param obj the interaction of interest (MC or data)
   * @return the number of muons counted in the topology
   **/
  template<class T>
    unsigned int count_muons(const T& obj)
    {
      return count_topology_particle(obj.topology, pvars::kMuon);
    }

  /**
   * @brief count pions in the topology string for an interaction
   * @details call count_topology_particle with pid = pvars::kPion
   * @param obj the interaction of interest (MC or data)
   * @return the number of pions counted in the topology
   **/
  template<class T>
    unsigned int count_pions(const T& obj)
    {
      return count_topology_particle(obj.topology, pvars::kPion);
    }

  /**
   * @brief count protons in the topology string for an interaction
   * @details call count_topology_particle with pid = pvars::kProton
   * @param obj the interaction of interest (MC or data)
   * @return the number of protons counted in the topology
   **/
  template<class T>
    unsigned int count_protons(const T& obj)
    {
      return count_topology_particle(obj.topology, pvars::kProton);
    }
    
  /**
   * @brief count photons in the topology string for an interaction (double version)
   * @details call count_photons and cast the return as a double
   * @param obj the interaction of interest (MC or data)
   * @return the number of photons counted in the topology as a double
   **/
  template<class T>
    double count_photons_double(const T& obj)
    {
      return count_photons(obj);
    }
    
  /**
   * @brief count electrons in the topology string for an interaction (double version)
   * @details call count_electrons and cast the return as a double
   * @param obj the interaction of interest (MC or data)
   * @return the number of electrons counted in the topology as a double
   **/
  template<class T>
    double count_electrons_double(const T& obj)
    {
      return count_electrons(obj);
    }
    
  /**
   * @brief count muons in the topology string for an interaction (double version)
   * @details call count_muons and cast the return as a double
   * @param obj the interaction of interest (MC or data)
   * @return the number of muons counted in the topology as a double
   **/
  template<class T>
    double count_muons_double(const T& obj)
    {
      return count_muons(obj);
    }
    
  /**
   * @brief count pions in the topology string for an interaction (double version)
   * @details call count_pions and cast the return as a double
   * @param obj the interaction of interest (MC or data)
   * @return the number of pions counted in the topology as a double
   **/
  template<class T>
    double count_pions_double(const T& obj)
    {
      return count_pions(obj);
    }
    
  /**
   * @brief count protons in the topology string for an interaction (double version)
   * @details call count_protons and cast the return as a double
   * @param obj the interaction of interest (MC or data)
   * @return the number of protons counted in the topology as a double
   **/
  template<class T>
    double count_protons_double(const T& obj)
    {
      return count_protons(obj);
    }

} // end vars::nc namespace
#endif

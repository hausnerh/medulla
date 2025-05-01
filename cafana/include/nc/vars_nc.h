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
  
  /**
   * @brief Count the number of particles in the interaction which are above threshold
   * @details this method is similar to that used by count_primaries defined in utilities.h,
   * However here we only count the particles of the specified pid
   * @param obj the interaction of interest (MC or data)
   * @param pid the pvars::Particle_t particle enumeration we are trying to count
   **/
  template<class T>
    unsigned int count_particle_above_threshold(const T& obj, int pid)
    {
      unsigned int count = 0;
      for (auto& particle : obj.particles)
      {
        if ((particle.pid == pid) && pcuts::final_state_signal(particle))
         ++count;
      }
      return count;
    }

  /**
   * @brief count photons in the interaction which are above threshold
   * @details call count_particle_above_threshold with pid = pvars::kPhoton
   * @param obj the interaction of interest (MC or data)
   * @return the number of photons with energy above 25 MeV 
   **/
  template<class T>
    unsigned int count_photons_above_threshold(const T& obj)
    {
      return count_particle_above_threshold(obj, pvars::kPhoton);
    }

  /**
   * @brief count electrons in the interaction which are above threshold
   * @details call count_particle_above_threshold with pid = pvars::kElectron
   * @param obj the interaction of interest (MC or data)
   * @return the number of electrons with energy above 25 MeV 
   **/
  template<class T>
    unsigned int count_electrons_above_threshold(const T& obj)
    {
      return count_particle_above_threshold(obj, pvars::kElectron);
    }

  /**
   * @brief count muons in the interaction which are above threshold
   * @details call count_particle_above_threshold with pid = pvars::kMuon
   * @param obj the interaction of interest (MC or data)
   * @return the number of muons with energy above 143.425 MeV 
   **/
  template<class T>
    unsigned int count_muons_above_threshold(const T& obj)
    {
      return count_particle_above_threshold(obj, pvars::kMuon);
    }

  /**
   * @brief count pions in the interaction which are above threshold
   * @details call count_particle_above_threshold with pid = pvars::kPion
   * @param obj the interaction of interest (MC or data)
   * @return the number of pions with energy above 25 MeV 
   **/
  template<class T>
    unsigned int count_pions_above_threshold(const T& obj)
    {
      return count_particle_above_threshold(obj, pvars::kPion);
    }

  /**
   * @brief count protons in the interaction which are above threshold
   * @details call count_particle_above_threshold with pid = pvars::kProton
   * @param obj the interaction of interest (MC or data)
   * @return the number of protons with energy above 50 MeV 
   **/
  template<class T>
    unsigned int count_protons_above_threshold(const T& obj)
    {
      return count_particle_above_threshold(obj, pvars::kProton);
    }
    
  /**
   * @brief count photons in the interaction which are above threshold (double version)
   * @details call count_photons_above_threshold and cast the return as a double
   * @param obj the interaction of interest (MC or data)
   * @return the number of photons with energy above 25 MeV
   **/
  template<class T>
    double count_photons_above_threshold_double(const T& obj)
    {
      return count_photons_above_threshold(obj);
    }
    
  /**
   * @brief count electrons in the interaction which are above threshold (double version)
   * @details call count_electrons_above_threshold and cast the return as a double
   * @param obj the interaction of interest (MC or data)
   * @return the number of electrons with energy above 25 MeV
   **/
  template<class T>
    double count_electrons_above_threshold_double(const T& obj)
    {
      return count_electrons_above_threshold(obj);
    }
    
  /**
   * @brief count muons in the interaction which are above threshold (double version)
   * @details call count_muons_above_threshold and cast the return as a double
   * @param obj the interaction of interest (MC or data)
   * @return the number of muons with energy above 143.425 MeV
   **/
  template<class T>
    double count_muons_above_threshold_double(const T& obj)
    {
      return count_muons_above_threshold(obj);
    }
    
  /**
   * @brief count pions in the interaction which are above threshold (double version)
   * @details call count_pions_above_threshold and cast the return as a double
   * @param obj the interaction of interest (MC or data)
   * @return the number of pions with energy above 25 MeV
   **/
  template<class T>
    double count_pions_above_threshold_double(const T& obj)
    {
      return count_pions_above_threshold(obj);
    }
    
  /**
   * @brief count protons in the interaction which are above threshold (double version)
   * @details call count_protons_above_threshold and cast the return as a double
   * @param obj the interaction of interest (MC or data)
   * @return the number of protons with energy above 25 MeV
   **/
  template<class T>
    double count_protons_above_threshold_double(const T& obj)
    {
      return count_protons_above_threshold(obj);
    }

  /**
   * @brief What is the max energy of the specified particle
   * @details iterate through the interaction's primary daughters
   * and return the max energy of the specified particle
   * @param obj the interaction of interest (MC or data)
   * @param particleType which particle from the pvars::Particle_t enum
   * @return the max energy of all primary particleType
   **/
  template<class T>
    double max_particle_energy(const T& obj, const pvars::Particle_t& particleType)
    {
      double pE = std::numeric_limits<double>::lowest();
      for (auto& particle : obj.particles)
        if ((particle.pid == particleType) && pcuts::is_primary(particle))
          pE = (pE > (pvars::energy(particle) / 1000.)) ? pE : pvars::energy(particle) / 1000.;
      return pE;
    }

  /**
   * @brief Max energy of primary photons
   * @details Apply max_particle_energy with particleType = pvars::kPhoton
   * @param obj the interaction of interest (MC or data)
   * @return the max energy of all primary photons
   **/
  template<class T>
    double max_photon_energy(const T& obj)
    {
      return max_particle_energy(obj, pvars::kPhoton);
    }

  /**
   * @brief Max energy of primary protons
   * @details Apply max_particle_energy with particleType = pvars::kProton
   * @param obj the interaction of interest (MC or data)
   * @return the max energy of all primary protons
   **/
  template<class T>
    double max_proton_energy(const T& obj)
    {
      return max_particle_energy(obj, pvars::kProton);
    }

  /**
   * @brief The transverse momentum fraction of the particle (of specified type) with the most energy
   * @details Find the primary particleType particle with the largest energy,
   * and return that particle's transverse momentum fraction. If no such particle exists return nonsense
   * @param obj the interaction of interest (MC or data)
   * @param particleType the pvars::Particle_t particle enumeration we are trying to identify
   * @return the transverse momentum fraction
   **/
  template<class T>
    double primary_particle_dpT_frac(const T& obj, const pvars::Particle_t& particleType)
    {
      double pE = std::numeric_limits<double>::lowest();
      double dpT_frac = std::numeric_limits<double>::lowest();
      for (auto& particle : obj.particles)
        if ((particle.pid == particleType) && pcuts::is_primary(particle) && (pE < (pvars::energy(particle) / 1000.)))
        {
          pE = pvars::energy(particle) / 1000.;
          dpT_frac = pvars::dpT(particle) / utilities::magnitude(utilities::to_three_vector(particle.momentum));
        }
      return dpT_frac;
    }

  /**
   * @brief The transverse momentum of the particle (of specified type) with the most energy
   * @details Find the primary particleType particle with the largest energy,
   * and return that particle's transverse momentum. If no such particle exists return nonsense
   * @param obj the interaction of interest (MC or data)
   * @param particleType the pvars::Particle_t particle enumeration we are trying to identify
   * @return the transverse momentum
   **/
  template<class T>
    double primary_particle_dpT(const T& obj, const pvars::Particle_t& particleType)
    {
      double pE = std::numeric_limits<double>::lowest();
      double dpT = std::numeric_limits<double>::lowest();
      for (auto& particle : obj.particles)
        if ((particle.pid == particleType) && pcuts::is_primary(particle) && (pE < (pvars::energy(particle) / 1000.)))
        {
          pE = pvars::energy(particle) / 1000.;
          dpT = pvars::dpT(particle);
        }
      return dpT;
    }

  /**
   * @brief The length of the particle (of specified type) with the most energy
   * @details Find the primary particleType particle with the largest energy,
   * and return that particle's length. If no such particle exists return nonsense
   * @param obj the interaction of interest (MC or data)
   * @param particleType the pvars::Particle_t particle enumeration we are trying to identify
   * @return the length of the particle track/shower
   **/
  template<class T>
    double primary_particle_length(const T& obj, const pvars::Particle_t& particleType)
    {
      double pE = std::numeric_limits<double>::lowest();
      double length = std::numeric_limits<double>::lowest();
      for (auto& particle : obj.particles)
        if ((particle.pid == particleType) && pcuts::is_primary(particle) && (pE < (pvars::energy(particle) / 1000.)))
        {
          pE = pvars::energy(particle) / 1000.;
          length = pvars::length(particle);
        }
      return length;
    } 

  /**
   * @brief Primary photon transverse momentum fraction
   * @details Apply primary_particle_dpT_frac with particleType == pvars::kPhoton
   * @param obj the interaction of interest (MC or data)
   * @return the transverse momentum fraction of the photon
   **/
  template<class T>
    double primary_photon_dpT_frac(const T& obj)
    {
      return primary_particle_dpT_frac(obj, pvars::kPhoton);
    }

  /**
   * @brief Primary photon transverse momentum
   * @details Apply primary_particle_dpT with particleType == pvars::kPhoton
   * @param obj the interaction of interest (MC or data)
   * @return the transverse momentum of the photon
   **/
  template<class T>
    double primary_photon_dpT(const T& obj)
    {
      return primary_particle_dpT(obj, pvars::kPhoton);
    }

  /**
   * @brief Primary photon length
   * @details Apply primary_particle_length with particleType == pvars::kPhoton
   * @param obj the interaction of interest (MC or data)
   * @return the length of the photon
   **/
  template<class T>
    double primary_photon_length(const T& obj)
    {
      return primary_particle_length(obj, pvars::kPhoton);
    }

  /**
   * @brief Primary proton transverse momentum fraction
   * @details Apply primary_particle_dpT_frac with particleType == pvars::kProton
   * @param obj the interaction of interest (MC or data)
   * @return the transverse momentum fraction of the proton
   **/
  template<class T>
    double primary_proton_dpT_frac(const T& obj)
    {
      return primary_particle_dpT_frac(obj, pvars::kProton);
    }

  /**
   * @brief Primary proton transverse momentum
   * @details Apply primary_particle_dpT with particleType == pvars::kProton
   * @param obj the interaction of interest (MC or data)
   * @return the transverse momentum of the proton
   **/
  template<class T>
    double primary_proton_dpT(const T& obj)
    {
      return primary_particle_dpT(obj, pvars::kProton);
    }

  /**
   * @brief Primary proton length
   * @details Apply primary_particle_length with particleType == pvars::kProton
   * @param obj the interaction of interest (MC or data)
   * @return the length of the proton
   **/
  template<class T>
    double primary_proton_length(const T& obj)
    {
      return primary_particle_length(obj, pvars::kProton);
    }

  /**
   * @brief Calculate the distance between the primary photon's and primary proton's vertices
   * @details Intended for the 1g1p selection. Find the photon/proton with the most energy
   * and caluclate the distance between their start positions
   * @param obj the interaction of interest (MC or data)
   * @return the distance between the starts of the photon and proton
   **/
  template<class T>
    double photon_proton_gap(const T& obj)
    {
      utilities::three_vector phStart;
      utilities::three_vector prStart;
      double phE = std::numeric_limits<double>::lowest();
      double prE = std::numeric_limits<double>::lowest();
      for (auto& particle : obj.particles)
      {
        if ((particle.pid == pvars::kPhoton) && pcuts::is_primary(particle) && (phE < (pvars::energy(particle) / 1000.)))
        {
          phE = pvars::energy(particle) / 1000.;
          phStart = utilities::to_three_vector(particle.start_point);
        }
        if ((particle.pid == pvars::kProton) && pcuts::is_primary(particle) && (prE < (pvars::energy(particle) / 1000.)))
        {
          prE = pvars::energy(particle) / 1000.;
          prStart = utilities::to_three_vector(particle.start_point);
        }
      }
      double gap = std::numeric_limits<double>::lowest();
      if ((phE != std::numeric_limits<double>::lowest()) &&
          (prE != std::numeric_limits<double>::lowest())  )
      {
        utilities::three_vector gap_vec = utilities::subtract(phStart, prStart);
        gap = utilities::magnitude(gap_vec);
      }
      return gap;
    }

  /**
   * @brief Calculate the cosine the angle between the primary photon's and primary proton's momenta
   * @details Intended for the 1g1p selection. Find the photon/proton with the most energy
   * and caluclate the cosine of the angle between their momenta
   * @param obj the interaction of interest (MC or data)
   * @return the cosine of the angle between the starts of the photon and proton momenta
   **/
  template<class T>
    double photon_proton_cosTh(const T& obj)
    {
      utilities::three_vector phP;
      utilities::three_vector prP;
      double phE = std::numeric_limits<double>::lowest();
      double prE = std::numeric_limits<double>::lowest();
      for (auto& particle : obj.particles)
      {
        if ((particle.pid == pvars::kPhoton) && pcuts::is_primary(particle) && (phE < (pvars::energy(particle) / 1000.)))
        {
          phE = pvars::energy(particle) / 1000.;
          phP = utilities::to_three_vector(particle.momentum);
        }
        if ((particle.pid == pvars::kProton) && pcuts::is_primary(particle) && (prE < (pvars::energy(particle) / 1000.)))
        {
          prE = pvars::energy(particle) / 1000.;
          prP = utilities::to_three_vector(particle.momentum);
        }
      }
      double cosTh = std::numeric_limits<double>::lowest();
      if ((phE != std::numeric_limits<double>::lowest()) &&
          (prE != std::numeric_limits<double>::lowest())  )
      {
        cosTh = utilities::dot_product(phP, prP) / (utilities::magnitude(phP) * utilities::magnitude(prP));
      }
      return cosTh;
    }

  /**
   * @brief The polar angle of the primary photon
   * @param obj the interaction of interest (MC or data)
   * @return The polar angle of the primary photon
   **/
  template<class T>
    double photon_polar_angle(const T& obj)
    {
      double pE = std::numeric_limits<double>::lowest();
      double angle = std::numeric_limits<double>::lowest();
      for (auto& particle : obj.particles)
        if ((particle.pid == pvars::kPhoton) && pcuts::is_primary(particle) && (pE < (pvars::energy(particle) / 1000.)))
        {
          pE = pvars::energy(particle) / 1000.;
          angle = pvars::polar_angle(particle);
        }
      return angle;
    }

  /**
   * @brief The polar angle of the primary photon
   * @param obj the interaction of interest (MC or data)
   * @return The polar angle of the primary photon
   **/
  template<class T>
    double photon_azimuthal_angle(const T& obj)
    {
      double pE = std::numeric_limits<double>::lowest();
      double angle = std::numeric_limits<double>::lowest();
      for (auto& particle : obj.particles)
        if ((particle.pid == pvars::kPhoton) && pcuts::is_primary(particle) && (pE < (pvars::energy(particle) / 1000.)))
        {
          pE = pvars::energy(particle) / 1000.;
          angle = pvars::azimuthal_angle(particle);
        }
      return angle;
    }

  // RECO ONLY VARS

  /**
   * @brief If the event has a single reco electron, what is the electron softmax score?
   * @param obj the interaction of interest (data only)
   * @return the electron's electron softmax score
   **/
  double electron_electron_id(const caf::SRInteractionDLPProxy& obj)
  {
    double softmax = std::numeric_limits<double>::lowest();
    std::vector<uint32_t> topology = utilities::count_primaries(obj);
    if (topology[pvars::kElectron] == 1)
    {
      for (auto& particle : obj.particles)
      {
        if ((particle.pid == pvars::kElectron) && pcuts::final_state_signal(particle))
        {
          std::vector<double> softmax_vec = pvars::particle_softmax_vec(particle);
          softmax = softmax_vec[pvars::kElectron];
          break;
        }
      }
    }
    return softmax;
  }

  /**
   * @brief If the event has a single reco electron, what is the electron's photon softmax score?
   * @param obj the interaction of interest (data only)
   * @return the electron's photon softmax score
   **/
  double electron_photon_id(const caf::SRInteractionDLPProxy& obj)
  {
    double softmax = std::numeric_limits<double>::lowest();
    std::vector<uint32_t> topology = utilities::count_primaries(obj);
    if (topology[pvars::kElectron] == 1)
    {
      for (auto& particle : obj.particles)
      {
        if ((particle.pid == pvars::kElectron) && pcuts::final_state_signal(particle))
        {
          std::vector<double> softmax_vec = pvars::particle_softmax_vec(particle);
          softmax = softmax_vec[pvars::kPhoton];
          break;
        }
      }
    }
    return softmax;
  }

  /**
   * @brief If the event has a single reco electron, what is the difference
   * between the electron's photon softmax and electron softmax
   * @details since the softmax goes from 0 to 1 for both scorres,
   * 1 is photon-like, -1 is electron-like, 0 is indistiguishable
   * @param obj the interaction of interest (data only)
   * @return the electron's electron softmax score
   **/
  double electron_photon_skew(const caf::SRInteractionDLPProxy& obj)
  {
    double skew = std::numeric_limits<double>::lowest();
    std::vector<uint32_t> topology = utilities::count_primaries(obj);
    if (topology[pvars::kElectron] == 1)
    {
      for (auto& particle : obj.particles)
      {
        if ((particle.pid == pvars::kElectron) && pcuts::final_state_signal(particle))
        {
          std::vector<double> softmax_vec = pvars::particle_softmax_vec(particle);
          skew = softmax_vec[pvars::kPhoton] - softmax_vec[pvars::kElectron];
          break;
        }
      }
    }
    return skew;
  }

  /**
   * @brief If the interaction has 2 photons, how aligned are they
   * @details -1 is anti-aligned, 0 is unaligned, 1 is aligned
   * @param obj the interaction of interest (data only)
   * @return photons' alignment
   **/
  double photon_alignment(const caf::SRInteractionDLPProxy& obj)
  {
    double alignment = std::numeric_limits<double>::lowest();
    std::vector<uint32_t> topology = utilities::count_primaries(obj);
    if (topology[pvars::kPhoton] == 2)
    {
      std::vector<utilities::three_vector> momenta;
      for (auto& particle : obj.particles)
      {
        if ((particle.pid == pvars::kPhoton) && pcuts::final_state_signal(particle))
          momenta.push_back(utilities::to_three_vector(particle.momentum));
      }
      if (momenta.size() != 2)
        return alignment;
      alignment = utilities::dot_product(momenta[0], momenta[1]);
      alignment /= (utilities::magnitude(momenta[0]) * utilities::magnitude(momenta[1]));
    }
    return alignment;
  }

  /**
   * @brief If the interaction has 2 photons, what is the gap, if any, from the end of the first
   * to the start of the second
   * @details Check that the separation could be explained by the momentum, and return -5 if it is not
   * @param obj the interaction of interest (data only)
   * @return photons' separation
   **/
  double photon_separation(const caf::SRInteractionDLPProxy& obj)
  {
    double separation = std::numeric_limits<double>::lowest();
    std::vector<uint32_t> topology = utilities::count_primaries(obj);
    if (topology[pvars::kPhoton] == 2)
    {
      std::vector<utilities::three_vector> starts;
      std::vector<utilities::three_vector> stops;
      std::vector<utilities::three_vector> momenta;
      for (auto& particle : obj.particles)
      {
        if ((particle.pid == pvars::kPhoton) && pcuts::final_state_signal(particle))
        {
          starts.push_back(utilities::to_three_vector(particle.start_point));
          stops.push_back(utilities::to_three_vector(particle.end_point));
          momenta.push_back(utilities::to_three_vector(particle.momentum));
        }
      }
      if (starts.size() != 2)
        return separation;
      utilities::three_vector oneToTwo = utilities::subtract(starts[1], stops[0]);
      utilities::three_vector twoToOne = utilities::subtract(starts[0], stops[1]);
      if      (utilities::dot_product(oneToTwo, momenta[0]) > 0)
      {
        separation = utilities::magnitude(oneToTwo);
      }
      else if (utilities::dot_product(twoToOne, momenta[1]) > 0)
      {
        separation = utilities::magnitude(twoToOne);
      }
      else
      {
        separation = -5.;
      }
    }
    return separation;
  }

  /**
   * @brief The pid_type softmax scores paricle with the highest energy of type target_type
   * @details Loop over the particles in the event to select the particle of the specified type
   * and return a vector of pid scores
   * @param obj the interaction of interest (data only)
   * @param target_type the type of particle in the interaction you want the pid scores for
   * @return the softmax scores for the particle
   **/
  std::vector<double> primary_particle_ids(const caf::SRInteractionDLPProxy& obj, const pvars::Particle_t& target_type)
  { 
    double pE = std::numeric_limits<double>::lowest();
    std::vector<double> pids(5, std::numeric_limits<double>::lowest());
    for (auto& particle : obj.particles)
      if ((particle.pid == target_type) && pcuts::is_primary(particle) && (pE < (pvars::energy(particle) / 1000.)))
      {
        pE = pvars::energy(particle) / 1000.;
        pids = pvars::particle_softmax_vec(particle);
      }
    return pids;
  }

  /**
   * @brief The pid_type softmax score for the target_type paricle with the highest energy
   * @details Loop over the particles in the event to select the particle of the specified type
   * and extract how much it looks like the second specified type
   * @param obj the interaction of interest (data only)
   * @param target_type the type of particle in the interaction you want the pid score for
   * @param pid_type the type of pid score you want
   * @return the softmax score for the particle
   **/
  double primary_particle_particle_id(const caf::SRInteractionDLPProxy& obj, const pvars::Particle_t& target_type,
                                                    const pvars::Particle_t& pid_type)
  {
    if (pid_type == pvars::kUnknown)
      return std::numeric_limits<double>::lowest();
    return primary_particle_ids(obj, target_type)[pid_type];
  }

  /**
   * @brief How much does the photon look like a photon?
   * @details call primary_particle_particle_id on (pvars::kPhoton, pvars::kPhoton)
   * @param obj the interaction of interest (data only)
   * @return the primary photon's photon softmax score
   **/
  double primary_photon_photon_score(const caf::SRInteractionDLPProxy& obj)
  {
    return primary_particle_particle_id(obj, pvars::kPhoton, pvars::kPhoton);
  }

  /**
   * @brief How much does the photon look like an electron?
   * @details call primary_particle_particle_id on (pvars::kPhoton, pvars::kElectron)
   * @param obj the interaction of interest (data only)
   * @return the primary photon's electron softmax score
   **/
  double primary_photon_electron_score(const caf::SRInteractionDLPProxy& obj)
  {
    return primary_particle_particle_id(obj, pvars::kPhoton, pvars::kElectron);
  }

  /**
   * @brief What does the primary photon second most look like?
   * @details Return the particle the primary photon is second closest to.
   * @param obj the interaction of interest (data only)
   * @return the particle type, as a double, which the photon second most looks like
   **/
  double photon_secondary_classification(const caf::SRInteractionDLPProxy& obj)
  {
    std::vector<double> pids = primary_particle_ids(obj, pvars::kPhoton);
    double maxPid = std::numeric_limits<double>::lowest();
    double sndPid = std::numeric_limits<double>::lowest();
    pvars::Particle_t frstClass = pvars::kUnknown;
    pvars::Particle_t scndClass = pvars::kUnknown;
    std::vector<pvars::Particle_t> types = {pvars::kPhoton,  
                                            pvars::kElectron,
                                            pvars::kMuon,
                                            pvars::kPion,
                                            pvars::kProton};
    for (auto const& idx : types)
    {
      double pid = pids[idx];
      if (pid > maxPid)
      {
        sndPid = maxPid;
        scndClass = frstClass;
        maxPid = pid;
        frstClass = idx;
      }
      else if (pid > sndPid)
      {
        sndPid = pid;
        scndClass = idx;
      }
    }
    return scndClass;
  }
} // end vars::nc namespace
#endif

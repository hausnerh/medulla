/**
 * @file core_gOre.h
 * @brief Header file for defining variables useful for NC analyses
 * specficially for saving photon/electron discrimination for the end
 * @author hhausner@fnal.gov
 **/

#ifndef CORE_GORE_H
#define CORE_GORE_H

/**
 * Threshold levels for particle KE
 * Nominal is 143.425 MeV for muons
 * 50 MeV for protons
 * 25 for everything else
 **/
#define GORE_MIN_GORE_ENERGY   10.0 // MeV
#define GORE_MIN_MUON_ENERGY   25.0 // MeV
#define GORE_MIN_PROTON_ENERGY 25.0 // MeV
#define GORE_MIN_PION_ENERGY   25.0 // MeV

/**
 * Wall Locations in cm
 * getting geom info from https://sbn-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=21693&filename=ICARUS_geometry_update_26Apr21_v3.pdf&version=3
 **/
#define GORE_WALL_X_POS  358.49
#define GORE_WALL_X_NEG -358.49
#define GORE_WALL_Y_POS  134.96
#define GORE_WALL_Y_NEG -181.89
#define GORE_WALL_Z_POS  894.95
#define GORE_WALL_Z_NEG -894.95

/**
 * Fiducial Thresholds in cm
 * Nominally 25 cm for X/Y walls, 50 cm from +Z and 30 cm from -Z
 **/
#define GORE_FID_THRESH_X_POS 50.0
#define GORE_FID_THRESH_X_NEG 50.0
#define GORE_FID_THRESH_Y_POS 50.0
#define GORE_FID_THRESH_Y_NEG 50.0
#define GORE_FID_THRESH_Z_POS 50.0
#define GORE_FID_THRESH_Z_NEG 50.0

#include <optional>

#include "include/utilities.h"

#include "include/gOre/selectors_gOre.h"

/**
 * @namespace core::gOre
 * @brief Base classes and definitions for NC single photon analyses
 * @details These classes in priciple could be used for other analyses, so do not hesitate to used them
 **/
namespace core::gOre
{
  /**
   * @class mc_topo_particle
   * @brief Lightweight class to store/access essential info from the SRTrueParticles
   **/
  class mc_topo_particle
  {
    public:
      /**
       * @brief basic constructor
       * @details get the PDG code, energy, particle start, and containment flag
       * @note the energy is converted into MeV from GeV
       * @todo would like to have a bounding box around the particle so we can check containment on the fly
       **/
      mc_topo_particle(const caf::Proxy<caf::SRTrueParticle>& particle)
      {
        pdg_code = particle.pdg;
        gen_energy = particle.genE * 1000.; // convert from GeV to MeV
        start = {particle.start.x, particle.start.y, particle.start.z};
        contained = particle.contained;
      }
      /** @brief the PDG code identifying the particle **/
      int pdg() const
      {
        return pdg_code;
      }
      /** @brief the particle's energy (at generation) in MeV **/
      double energy() const
      {
        return gen_energy;
      }
      /** @brief the starting location of the particle in detector coordinates (cm) **/
      utilities::three_vector get_start() const
      {
        return start;
      }
      /** @brief is the particle origin contained in the detector's fiducial volume? **/
      bool is_fiducial(std::vector<double> thresholds = {GORE_FID_THRESH_X_POS, GORE_FID_THRESH_X_NEG, GORE_FID_THRESH_Y_POS, GORE_FID_THRESH_Y_NEG, GORE_FID_THRESH_Z_POS, GORE_FID_THRESH_Z_NEG}) const
      {
        double vtx_x = std::get<0>(start);
        double vtx_y = std::get<1>(start);
        double vtx_z = std::get<2>(start);
        return (GORE_WALL_X_POS - vtx_x > thresholds.at(0)) &&
               (vtx_x - GORE_WALL_X_NEG > thresholds.at(1)) &&
               (GORE_WALL_Y_POS - vtx_y > thresholds.at(2)) &&
               (vtx_y - GORE_WALL_Y_NEG > thresholds.at(3)) &&
               (GORE_WALL_Z_POS - vtx_z > thresholds.at(4)) &&
               (vtx_z - GORE_WALL_Z_NEG > thresholds.at(5)) ;
      }
      /** @brief is the particle contained **/
      bool is_contained() const
      {
        return contained;
      }
      /** @brief is the particle a photon? **/
      bool is_photon() const
      {
        return (std::abs(pdg_code) == 22);
      }
      /** @brief is the particle an electron? **/
      bool is_electron() const
      {
        return (std::abs(pdg_code) == 11);
      }
      /** @brief is the particle a proton? **/
      bool is_proton() const
      {
        return (std::abs(pdg_code) == 2212);
      }
      /** @brief is the particle a neutron? **/
      bool is_neutron() const
      {
        return (std::abs(pdg_code) == 2112);
      }
      /** @brief is the particle a muon? **/
      bool is_muon() const
      {
        return (std::abs(pdg_code) == 13);
      }
      /** @brief is the particle a charged pion? **/
      bool is_charged_pion() const
      {
        return (std::abs(pdg_code) == 211);
      }
      /** @brief is the particle a neutral pion? **/
      bool is_neutral_pion() const
      {
        return (std::abs(pdg_code) == 111);
      }
      /** @brief is the particle a neutrino? **/
      bool is_neutrino() const
      {
        // flavor doesn't matter here
        return (std::abs(pdg_code) == 12 ||
                std::abs(pdg_code) == 14 ||
                std::abs(pdg_code) == 16 );
      }
      /** @brief is the particle above threshold? **/
      bool is_above_threshold(std::vector<double> thresholds = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY}) const
      {
        // only check photons, protons, electrons, muons, and pions (since those are the ones with thresholds)
        // since the thresholds are for kinetic energy, subtract off the mass
        // neutrinos are are subthreshold (they're outgoing and basicaly invisible)
        // everything else say it's above
        if (is_neutrino())
        {
          return false;
        }
        else if (is_photon())
        {
          return (gen_energy > thresholds.at(0));
        }
        else if (is_electron())
        {
          return (gen_energy - ELECTRON_MASS > thresholds.at(0));
        }
        else if (is_muon())
        {
          return (gen_energy - MUON_MASS > thresholds.at(1));
        }
        else if (is_proton())
        {
          return (gen_energy - PROTON_MASS > thresholds.at(2));
        }
        else if (is_charged_pion())
        {
          return (gen_energy - PION_MASS > thresholds.at(3));
        }
        return true;
      }
    private:
      int pdg_code;                   ///< particle's PDG code 
      double gen_energy;              ///< particle's energy at generation
      utilities::three_vector start;  ///< where the particle was generated
      bool contained;                 ///< was the particle contained?
  };

  /**
   * @class mc_topology
   * @brief A collection of mc_topo_particles
   **/
  class mc_topology
  {
    /** 
     * @brief Basic mc_topology constructor
     * @note Assume the first particle in the vector is the primary lepton
     **/
    public:
      mc_topology(const caf::Proxy<std::vector<caf::SRTrueParticle>>& particle_vec, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY,
                                                                                                                  GORE_FID_THRESH_X_POS, GORE_FID_THRESH_X_NEG, GORE_FID_THRESH_Y_POS, GORE_FID_THRESH_Y_NEG, GORE_FID_THRESH_Z_POS, GORE_FID_THRESH_Z_NEG})
      {
        // init thresholds
        energy_thresholds = (params.size() < 4) ? std::vector({GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY}) : std::vector({params.at(0), params.at(1), params.at(2), params.at(3)});
        fiducial_thresholds = (params.size() < 10) ? std::vector({GORE_FID_THRESH_X_POS, GORE_FID_THRESH_X_NEG, GORE_FID_THRESH_Y_POS, GORE_FID_THRESH_Y_NEG, GORE_FID_THRESH_Z_POS, GORE_FID_THRESH_Z_NEG}) : std::vector({params.at(4), params.at(5), params.at(6), params.at(7), params.at(8), params.at(9)});
        // verify non-empty vector
        if (particle_vec.size() == 0)
          throw std::runtime_error("mc_topology: cannot initialize from empty vector.");
        // first particle is the lepton. Use for fiducialization
        // caf::Proxy is annoying so we can't use first
        mc_topo_particle lepton(particle_vec[0]);
        fiducial = lepton.is_fiducial(fiducial_thresholds);
        vertex = lepton.get_start();
        // default assume event is contained. if anything over threshold is not, then the event isn't
        contained = true;
        for (auto const& particle : particle_vec)
        {
          mc_topo_particle cur_topo_particle(particle);
          if (cur_topo_particle.is_above_threshold(energy_thresholds))
          {
            contained = contained && cur_topo_particle.is_contained();
            add(cur_topo_particle);
          }
        }
      }
      /**
       * @brief Does the topology contian any particle of the given PDG code
       * @details This checks that the map has been initialized, the actual number
       * of particles would be given by the vector returned by the map.
       * @param pdg The target PDG code
       * @return true if there is at least 1 of the specified particle
       **/
      bool has(const int& pdg) const
      {
        return (particles_by_pdg.count(pdg) == 1);
      }
      /** @brief add a particle to the topology **/
      void add(const mc_topo_particle& topo_particle)
      {
        int pdg = topo_particle.pdg();
        if(has(pdg))
        {
          particles_by_pdg.at(pdg).push_back(topo_particle);
        }
        else
        {
          particles_by_pdg[pdg] = {topo_particle};
        }
      }
      /** 
       * @brief count the number of particles of PDG code in the topology
       * @param pdg The target PDG code
       * @return the number of specified particles in the topology
       **/
      unsigned int count(const int& pdg) const
      {
        return (has(pdg)) ? particles_by_pdg.at(pdg).size() : 0;
      }
      /** 
       * @brief count the number of (anti-)particles of PDG code in the topology
       * @details treat the particles and anti-particles the same (eg count both positive and negative pions)
       * @param pdg The target PDG code
       * @return the number of specified particles in the topology
       **/
      unsigned int count_with_antiparticles(const int& pdg) const
      {
        unsigned int counts = 0;
        if(has( pdg)) counts += particles_by_pdg.at( pdg).size();
        if(has(-pdg)) counts += particles_by_pdg.at(-pdg).size();
        return counts;
      }
      /** @brief pipe a basic overview of the topology to std::cout **/
      void report() const
      {
        double vtx_x = std::get<0>(vertex);
        double vtx_y = std::get<1>(vertex);
        double vtx_z = std::get<2>(vertex);
        std::cout << "~~ MC TOPOLOGY ~~\n"
                  << "~~ Vertex X: " << vtx_x << " ~~\n"
                  << "~~ Vertex Y: " << vtx_y << " ~~\n"
                  << "~~ Vertex Z: " << vtx_z << " ~~" << std::endl;
        for (auto const& [pdg, particles] : particles_by_pdg)
          std::cout << "  " << pdg << ": " << particles.size() << std::endl;
        std::cout << "~~ ~~ ~~ ~~ ~~ ~~" << std::endl;
      }
      /** @brief the total number of particles in the topology **/
      unsigned int total_particles() const
      {
        unsigned int counts = 0;
        for (auto const& [pdg, particles] : particles_by_pdg)
          counts += particles.size();
        return counts;
      }
      /** @brief does the topology contain exactly 1 photon? **/
      bool single_photon() const
      {
        return count(22) == 1;
      }
      /** @brief does the topology contain only photons and nucleons? **/
      bool only_photons_and_nucleons() const
      {
        bool is_it = true;
        for (auto const& [pdg, particles] : particles_by_pdg)
          if (pdg != 22 && pdg != 2112 && pdg != 2212)
          {
            is_it = false;
            break;
          }
        return is_it;
      }
      /** @brief does the toplogy have a neutral pion? **/
      bool has_pi0() const
      {
        return has(111);
      }
      /** @brief does the topology have a charged pion? **/
      bool has_pi_pm() const
      {
        return count_with_antiparticles(211) > 0;
      }
      /** @brief is the interaction vertex contained in the fiducial volume? **/
      bool is_fiducial() const
      {
        return fiducial;
      }
      /** @brief are all of the particles contained? **/
      bool is_contained() const
      {
        return contained;
      }
    private:
      std::unordered_map<int, std::vector<mc_topo_particle>> particles_by_pdg; ///< A map of the PDG codes to the particles
      bool contained;                                                          ///< Are all of the particles contained?
      bool fiducial;                                                           ///< Is the vertex inside the fiducial volume?
      utilities::three_vector vertex;                                          ///< The location of the vertex in the detector (cm)
      std::vector<double> energy_thresholds;                                   ///< The particle energy thresholds
      std::vector<double> fiducial_thresholds;                                 ///< The detector fiducial thresholds
  };

  /**
   * @brief is the partcle final state signal based on the thresholds set?
   * @note all particles which are not electrons, photons, muons, protons, or pions fail
   * @tparam T data or MC particle
   * @param p the particle
   * @param gore_min the energy threshold of signal photons/electrons (MeV)
   * @param muon_min the energy threshold of signal muons (MeV)
   * @param proton_min the energy threshold of signal protons (MeV)
   * @param pion_min the energy threshold of signal pions (MeV)
   * @return true if the particle is one we care about for signal definitions
   **/
  template<class T>
    bool final_state_signal(const T& p,
                            const double gore_min   = GORE_MIN_GORE_ENERGY,
                            const double muon_min   = GORE_MIN_MUON_ENERGY,
                            const double proton_min = GORE_MIN_PROTON_ENERGY,
                            const double pion_min   = GORE_MIN_PION_ENERGY)
    {
      bool is_signal(false);
      if (pcuts::is_primary(p))
      {
        pvars::Particle_t pid = static_cast<pvars::Particle_t>(pvars::pid(p));
        switch (pid)
        {
          case pvars::kMuon:
            is_signal = (pvars::ke(p) > muon_min);
            break;
          case pvars::kProton:
            is_signal = (pvars::ke(p) > proton_min);
            break;
          case pvars::kPion:
            is_signal = (pvars::ke(p) > pion_min);
            break;
          case pvars::kPhoton:
            is_signal = (pvars::ke(p) > gore_min);
            break;
          case pvars::kElectron:
            is_signal = (pvars::ke(p) > gore_min);
            break;
          default:
            break;
        }
      }
      return is_signal;
    }

  /**
   * @brief How many final state signal particles of each pvars::Paricle_t are there?
   * @tparam T data or MC interaction
   * @param obj the interaction
   * @param gore_min the energy threshold of signal photons/electrons (MeV)
   * @param muon_min the energy threshold of signal muons (MeV)
   * @param proton_min the energy threshold of signal protons (MeV)
   * @param pion_min the energy threshold of signal pions (MeV)
   * @return The vector of final state signal particle counts
   **/
  template<class T>
    std::vector<size_t> count_primaries(const T& obj,
                                        const double gore_min   = GORE_MIN_GORE_ENERGY,
                                        const double muon_min   = GORE_MIN_MUON_ENERGY,
                                        const double proton_min = GORE_MIN_PROTON_ENERGY,
                                        const double pion_min   = GORE_MIN_PION_ENERGY)
    {
      std::vector<size_t> primaries(5, 0); // five valid pvars::Paricle_t
      for (auto const& particle : obj.particles)
      {
        if (final_state_signal(particle, gore_min, muon_min, proton_min, pion_min))
          ++primaries.at(pvars::pid(particle));
      }
      return primaries;
    }

  /**
   * @struct ParticleType
   * @brief allow us to call data and MC particles with a templetated name
   * @details for an interaction of type T, ParticleType<T>::type is the type of particle in that interaction
   **/
  template<class>
  struct ParticleType;
  template<>
    struct ParticleType<caf::SRInteractionDLPProxy>
    {
      typedef caf::SRParticleDLPProxy type;
    };
  template<>
    struct ParticleType<caf::SRInteractionTruthDLPProxy>
    {
      typedef caf::SRParticleTruthDLPProxy type;
    };
  
  template <class T>
  using ParticleProxy_t = typename ParticleType<T>::type;

  /**
   * @brief return the primary gOre, if it exists
   * @tparam T the type of interaction (data or MC)
   * @param obj the interction
   * @return std::optional<ParticleProxy_t<T>> the primary gOre, if it exists
   **/
  template <class T>
    std::optional<ParticleProxy_t<T>> primary_gOre(const T& obj)
    {
      size_t idx = selectors::gOre::leading_gOre(obj);
      if (idx == kNoMatch)
        return std::nullopt;
      return obj.particles.at(idx);
    }

  /**
   * @brief return the secondary gOre, if it exists
   * @tparam T the type of interaction (data or MC)
   * @param obj the interction
   * @return std::optional<ParticleProxy_t<T>> the secondary gOre, if it exists
   **/
  template <class T>
    std::optional<ParticleProxy_t<T>> secondary_gOre(const T& obj)
    {
      size_t idx = selectors::gOre::subleading_gOre(obj);
      if (idx == kNoMatch)
        return std::nullopt;
      return obj.particles.at(idx);
    }

  /**
   * @struct Interaction
   * @brief hold your interaction (MC or data) and the single photon/electron
   * @details wraps around a data or MC interaction and allows easy access to parameters of interest
   * @tparam T data or MC interaction class
   **/
  template <class T>
    struct Interaction
    {
      typedef typename ParticleType<T>::type PT; ///< The type of particle stored in the interaction
      /** @brief basic Interaction constructor **/
      Interaction(const T& object, std::vector<double> params = {GORE_MIN_GORE_ENERGY, GORE_MIN_MUON_ENERGY, GORE_MIN_PROTON_ENERGY, GORE_MIN_PION_ENERGY}) : obj(object), thresholds(params)
      {
        bool found_Other(false);
        min_muon_ke = std::numeric_limits<double>::max();
        min_pion_ke = std::numeric_limits<double>::max();
        subleading_gore_ke = std::numeric_limits<double>::lowest();
        double leading_gore_ke = std::numeric_limits<double>::lowest();
        utilities::three_vector    leading_p;
        utilities::three_vector subleading_p;
        for (auto& particle : obj.particles)
        {
          pvars::Particle_t pid = static_cast<pvars::Particle_t>(pvars::pid(particle));
          if (pcuts::is_primary(particle))
          {
            double pKE = pvars::ke(particle);
            switch (pid)
            {
              case pvars::kMuon:
                min_muon_ke = (min_muon_ke > pKE) ? pKE : min_muon_ke;
                break;
              case pvars::kPion:
                min_pion_ke = (min_pion_ke > pKE) ? pKE : min_pion_ke;
                break;
              default:
                if (pid == pvars::kPhoton || pid == pvars::kElectron)
                {
                  ++nShowers;
                  total_gore_ke += pKE;
                  if (pKE > leading_gore_ke)
                  {
                    subleading_gore_ke = leading_gore_ke;
                    subleading_p = leading_p;
                    leading_gore_ke = pKE;
                    leading_p = utilities::to_three_vector(particle.momentum);
                  }
                  else if (pKE > subleading_gore_ke)
                  {
                    subleading_gore_ke = pKE;
                    subleading_p = utilities::to_three_vector(particle.momentum);
                  }
                }
                break;
            }
          }
          if (final_state_signal(particle, thresholds.at(0), thresholds.at(1), thresholds.at(2), thresholds.at(3)))
          {
            switch (pid)
            {
              case pvars::kProton:
                protons.push_back(&particle);
                break;
              case pvars::kPhoton:
                gOres.push_back(&particle);
                break;
              case pvars::kElectron:
                gOres.push_back(&particle);
                break;
              default:
                found_Other = true;
                break;
            }
          }
        }
        std::sort(protons.begin(), protons.end(), [](PT* a, PT* b){ return pvars::ke(*a) > pvars::ke(*b); });
        std::sort(gOres.begin(), gOres.end(), [](PT* a, PT* b){ return pvars::ke(*a) > pvars::ke(*b); });
        is_valid = (gOres.size() == 1) && (not found_Other);
        min_muon_ke = (min_muon_ke == std::numeric_limits<double>::max()) ? 0. : min_muon_ke;
        min_pion_ke = (min_pion_ke == std::numeric_limits<double>::max()) ? 0. : min_pion_ke;
        if (is_valid && pvars::ke(*primary_gOre()) < subleading_gore_ke)
        {
          std::string msg = "Primary gOre KE is less than subleading gOre KE!!!\n  Primary gOre KE: " + std::to_string(pvars::ke(*primary_gOre())) + "\n  Leading gOre KE: " + std::to_string(leading_gore_ke) + " (should be the same as above)\n  subleading gOre KE: " + std::to_string(subleading_gore_ke);
          throw std::runtime_error(msg);
        }
        pion_costh = std::numeric_limits<double>::max();
        if (gOres.size() > 0 && nShowers > 1)
        {
          utilities::three_vector gOre_1_p = utilities::to_three_vector(primary_gOre()->momentum);
          utilities::three_vector gOre_2_p = subleading_p;
          pion_costh = utilities::dot_product(gOre_1_p, gOre_2_p) / (utilities::magnitude(gOre_1_p) * utilities::magnitude(gOre_2_p));
        }
      }
      PT* primary_gOre()
      {
        if (gOres.size() > 0)
          return gOres.at(0);
        return nullptr;
      }
      size_t ngOres()
      {
        return gOres.size();
      }
      size_t nProtons()
      {
        return protons.size();
      }
      const T& obj;                   ///< The wrapped interaction
      std::vector<PT*> gOres;         ///< The photons & electrons in the final state signal, sorted by descending KE
      std::vector<PT*> protons;       ///< The protons in the final state signal, sorted by descending KE
      double min_muon_ke;             ///< lowest muon KE below threshold
      double min_pion_ke;             ///< lowest pion KE below threshold
      double subleading_gore_ke;      ///< what is the KE of the subleading gOre candidate?
      double pion_costh;              ///< if there is a subleading gOre candidate, what is the cosine of the angle between it and the primary?
      double total_gore_ke;           ///< what is the sum of KE for all gOre showers (above and below theshold)?
      bool is_valid;                  ///< are there no pions or muons above threshold?
      size_t nShowers;                ///< How many showers are there (including sub-threshold)
      std::vector<double> thresholds; ///< Particle Energy Thresholds
    };

  typedef Interaction<caf::SRInteractionDLPProxy>      Reco_Interaction; ///< handy typedef for reco interactions
  typedef Interaction<caf::SRInteractionTruthDLPProxy> True_Interaction; ///< handy typedef for MC interactions
} // end core::gOre namespace

#endif

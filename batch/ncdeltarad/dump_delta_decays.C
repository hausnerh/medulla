// Dump every NC Δ+/Δ0 decay in GEN-stage art files to CSV, for pilot_check.py.
//
//   (setup icaruscode first; gallery and nusimdata come with it)
//   root -l -b -q 'dump_delta_decays.C("gen_files.list", "delta_decays.csv")'
//
// GENIEHelper builds each simb::MCParticle with TrackId = GHEP index and
// Mother = GHEP first mother, so the Δ's decay products are the particles
// whose Mother() is the Δ's TrackId. The Δ itself is GHEP status 3
// (kIStDecayedState) after BaryonResonanceDecayer.
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

void dump_delta_decays(const std::string& filelist, const std::string& outcsv,
                       const std::string& label = "generator")
{
  std::vector<std::string> files;
  std::ifstream in(filelist);
  for (std::string f; std::getline(in, f);)
    if (!f.empty()) files.push_back(f);
  if (files.empty()) { std::cerr << "no input files in " << filelist << std::endl; return; }

  std::ofstream out(outcsv);
  out << "run,subrun,event,delta_pdg,w_delta,w_kine,hitnuc,channel,n_daughters,"
         "n_photon_status1_gt30mev,mother_pdg\n";

  art::InputTag const tag(label);
  long nev = 0, ndelta = 0;
  for (gallery::Event ev(files); !ev.atEnd(); ev.next()) {
    ++nev;
    auto const& aux = ev.eventAuxiliary();
    auto const& truths = *ev.getValidHandle<std::vector<simb::MCTruth>>(tag);
    for (auto const& t : truths) {
      if (!t.NeutrinoSet() || t.GetNeutrino().CCNC() != simb::kNC) continue;
      auto const& nu = t.GetNeutrino();
      int nphot = 0;
      for (int i = 0; i < t.NParticles(); ++i) {
        auto const& p = t.GetParticle(i);
        if (p.StatusCode() == 1 && p.PdgCode() == 22 && p.E() > 0.030) ++nphot;
      }
      for (int i = 0; i < t.NParticles(); ++i) {
        auto const& d = t.GetParticle(i);
        if ((d.PdgCode() != 2214 && d.PdgCode() != 2114) || d.StatusCode() != 3) continue;
        ++ndelta;
        std::string channel;
        int nd = 0;
        for (int j = 0; j < t.NParticles(); ++j) {
          auto const& c = t.GetParticle(j);
          if (c.Mother() != d.TrackId()) continue;
          ++nd;
          channel += (channel.empty() ? "" : ":") + std::to_string(c.PdgCode());
        }
        // Mother: the struck nucleon for a primary Δ; a heavier resonance
        // for a secondary one (R -> Δ π). 0 if not found.
        int mother_pdg = 0;
        for (int j = 0; j < t.NParticles(); ++j)
          if (t.GetParticle(j).TrackId() == d.Mother()) { mother_pdg = t.GetParticle(j).PdgCode(); break; }
        out << aux.run() << ',' << aux.subRun() << ',' << aux.event() << ','
            << d.PdgCode() << ',' << d.Momentum().M() << ',' << nu.W() << ','
            << nu.HitNuc() << ',' << channel << ',' << nd << ',' << nphot << ','
            << mother_pdg << '\n';
      }
    }
  }
  std::cout << "dump_delta_decays: " << nev << " events, " << ndelta
            << " NC Δ+/Δ0 decays -> " << outcsv << std::endl;
}

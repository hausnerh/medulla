/**
 * @file mctruth_cuts_gOre.h
 * @brief MCTruth-scope cuts specific to the gOre NC Δ→Nγ analysis.
 * @details These cuts operate directly on SRTrueInteraction and complement
 * the generic ones in `mctruth_cuts.h`. They cover the truth-level
 * conditions required to define the NC Δ→Nγ signal in TOML categories
 * (NC vs CC, Δ resonance, no muons), so the bespoke event-level
 * `has_nc_delta_ng_interaction` filter and the `mc_category` MCTruth
 * variable can be retired in favor of unified `[[category]]` blocks.
 * @author hhausner@fnal.gov
 **/
#ifndef MCTRUTH_CUTS_GORE_H
#define MCTRUTH_CUTS_GORE_H
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"

#include "include/framework.h"
#include "include/variables.h"

namespace mctruth::gOre
{
    /**
     * @brief Cut on the GENIE resonance number.
     * @details `obj.resnum` is the GENIE-internal resonance index.
     * The Δ(1232) corresponds to `resnum == 0`.
     * @tparam T the type of the object to apply the cut on.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @param params single-element vector containing the resonance number
     * to match (cast to int). Defaults to 0 (Δ(1232)).
     * @return true if `obj.resnum` equals the requested resonance number.
     **/
    template<typename T>
    bool is_resonance(const T & obj, std::vector<double> params={0.0,})
    {
        return obj.resnum == static_cast<int>(params.at(0));
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, is_resonance, is_resonance);

    /**
     * @brief Cut for zero true final-state muons above threshold.
     * @details Mirrors the structure of upstream `mctruth::single_muon`
     * but tests for absence rather than count==1. KE computed from `genE`.
     * @tparam T the type of the object to apply the cut on.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @param params KE threshold in MeV. Defaults to the gOre analysis
     * threshold of 25 MeV (note: differs from the SBN standard of 143.425
     * used in `single_muon`).
     * @return true if no muons above threshold.
     **/
    template<typename T>
    bool no_muons(const T & obj, std::vector<double> params={25.0,})
    {
        for(const auto & p : obj.prim)
        {
            if(std::abs(p.pdg) == 13)
            {
                double ke = 1000. * (p.genE - (MUON_MASS/1000.));
                if(ke >= params.at(0))
                    return false;
            }
        }
        return true;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, no_muons, no_muons);

    /**
     * @brief Composite shortcut: NC + Δ resonance.
     * @details Convenience cut equivalent to `!iscc && is_resonance(0)`.
     * Useful when a TOML `[[category]]` block wants to keep the cut list
     * compact.
     * @tparam T the type of the object to apply the cut on.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @return true if the interaction is NC and on the Δ(1232) resonance.
     **/
    template<typename T>
    bool nc_delta_resonance(const T & obj)
    {
        return (!obj.iscc) && (obj.resnum == 0);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, nc_delta_resonance, nc_delta_resonance);

} // namespace mctruth::gOre

#endif

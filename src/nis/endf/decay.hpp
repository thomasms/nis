#ifndef ENDF_DECAY_HPP
#define ENDF_DECAY_HPP

#include <functional>
#include <map>
#include <memory>

#include "nis/endf/decay/product_computer.hpp"


namespace nis {
    namespace endf {

        const std::map<int, std::function<std::unique_ptr<DecayProductComputer>()>> DECAY_MAP = {
            {0, []() { return std::make_unique<GammaRayReactionProducts>(); }},
            {1, []() { return std::make_unique<BetaMinusDecayReactionProducts>(); }},
            {2, []() { return std::make_unique<BetaPlusDecayReactionProducts>(); }},
            {3, []() { return std::make_unique<ITDecayReactionProducts>(); }},
            {4, []() { return std::make_unique<AlphaDecayReactionProducts>(); }},
            {5, []() { return std::make_unique<NeutronEmissionDecayReactionProducts>(); }},
            {6, []() { return std::make_unique<SFDecayReactionProducts>(); }},
            {7, []() { return std::make_unique<ProtonEmissionDecayReactionProducts>(); }},
            {8, []() { return std::make_unique<BetaMinusDecayReactionProducts>(); }},
            {9, []() { return std::make_unique<GammaRayReactionProducts>(); }},
            // unknown origin
            {10, []() { return std::make_unique<ToSinkDecayReactionProducts>(); }}
        };

        std::vector<int> expand_rtypes(double rtyp);

        DecayProductT compute_products(int parent_zai, double rtyp, int final_isomeric_state);
    }
}

#endif //ENDF_DECAY_HPP
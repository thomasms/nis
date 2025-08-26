#ifndef ENDF_DECAY_PRODUCT_COMPUTER_HPP
#define ENDF_DECAY_PRODUCT_COMPUTER_HPP

#include <variant>

#include "nis/zai.hpp"


namespace nis {
    namespace endf {

        // (primary, {secondaries})
        using DecayProductT = std::tuple<int, std::vector<int>>;

        struct GammaRayReactionProducts{

            DecayProductT get_products(int parent_zai) const{
                return DecayProductT(parent_zai, {});
            };

        };

        struct BetaMinusDecayReactionProducts{

            DecayProductT get_products(int parent_zai) const{
                const zai::triple parent_triple = zai::decompose(parent_zai);
                // Z increases by one, as a neutron goes to a proton
                // thus A remains the same.
                // i.e. C14 (60140) -> N14 (70140) + e- + mu_e
                const int daughter_zai = zai::make(std::get<0>(parent_triple)+1,
                                                   std::get<1>(parent_triple),
                                                   std::get<2>(parent_triple));
                // we don't care about the electron or neutrino
                return DecayProductT(daughter_zai, {});
            };
        };

        struct BetaPlusDecayReactionProducts{

            DecayProductT get_products(int parent_zai) const{
                const zai::triple parent_triple = zai::decompose(parent_zai);
                // Z decreases by one, as a proton goes to a neutron
                // thus A remains the same.
                // i.e. A26 (130260) + e- -> Mg26 (120260) + mu_e
                const int daughter_zai = zai::make(std::get<0>(parent_triple)-1,
                                                   std::get<1>(parent_triple),
                                                   std::get<2>(parent_triple));
                // we don't care about the electron or neutrino
                return DecayProductT(daughter_zai, {});
            };
        };

        struct ITDecayReactionProducts{

            DecayProductT get_products(int parent_zai) const{
                const zai::triple parent_triple = zai::decompose(parent_zai);
                // TODO: need to deal with isometric state
                const int daughter_zai = zai::make(std::get<0>(parent_triple),
                                                   std::get<1>(parent_triple),
                                                   std::get<2>(parent_triple));
                return DecayProductT(daughter_zai, {});
            };
        };

        struct AlphaDecayReactionProducts{

            DecayProductT get_products(int parent_zai) const{
                const zai::triple parent_triple = zai::decompose(parent_zai);
                const int daughter_zai = zai::make(std::get<0>(parent_triple)-2,
                                                   std::get<1>(parent_triple)-4,
                                                   std::get<2>(parent_triple));
                return DecayProductT(daughter_zai, {zai::ALPHA});
            };
        };

        struct NeutronEmissionDecayReactionProducts{

            DecayProductT get_products(int parent_zai) const{
                const zai::triple parent_triple = zai::decompose(parent_zai);
                const int daughter_zai = zai::make(std::get<0>(parent_triple),
                                                   std::get<1>(parent_triple)-1,
                                                   std::get<2>(parent_triple));
                return DecayProductT(daughter_zai, {});
            };
        };

        // this one is tricky :) - we need more information/data from MT=454,459
        struct SFDecayReactionProducts{

            DecayProductT get_products(int parent_zai) const{
                // const zai::triple parent_triple = zai::decompose(parent_zai);
                // todo: ....
                return DecayProductT(-1, {});
            };
        };

        struct ProtonEmissionDecayReactionProducts{

            DecayProductT get_products(int parent_zai) const{
                const zai::triple parent_triple = zai::decompose(parent_zai);
                const int daughter_zai = zai::make(std::get<0>(parent_triple)-1,
                                                   std::get<1>(parent_triple),
                                                   std::get<2>(parent_triple));
                return DecayProductT(daughter_zai, {});
            };
        };

        struct ToSinkDecayReactionProducts{

            DecayProductT get_products([[maybe_unused]] int parent_zai) const{
                // const zai::triple parent_triple = zai::decompose(parent_zai);
                return DecayProductT(-1, {});
            };
        };

        using DecayProductComputer = std::variant<GammaRayReactionProducts,
                                                  BetaMinusDecayReactionProducts,
                                                  BetaPlusDecayReactionProducts,
                                                  ITDecayReactionProducts,
                                                  AlphaDecayReactionProducts,
                                                  NeutronEmissionDecayReactionProducts,
                                                  SFDecayReactionProducts,
                                                  ProtonEmissionDecayReactionProducts,
                                                  ToSinkDecayReactionProducts>;

    }
}

#endif //ENDF_DECAY_PRODUCT_COMPUTER_HPP
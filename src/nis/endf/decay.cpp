#include <cmath>
#include <memory>

#include "nis/physics.hpp"
#include "nis/util.hpp"
#include "nis/zai.hpp"

#include "nis/endf/decay.hpp"


namespace nis {
    namespace endf {

        std::vector<int> expand_rtypes(double rtyp){
            // rtype is actually a double which can contain multiple reaction modes
            // i.e. 1.5 is beta minus followed by neutron emission
            // We assume that at most we can have 7 subsequent reactions, hence the 1e6 factor.
            // Since this is stored in the ENDF file, with limited precision, it should hold
            constexpr int _rtyp_factor = 1'000'000;
            const int initial_rtype = static_cast<int>(std::floor(rtyp));

            // if we are zero first then we are gamma
            if(initial_rtype == 0){
                return std::vector<int>{0};
            }

            // if we the initial rtype is not valid
            // then return an empty vector
            // regardless of what else comes next
            if(DECAY_MAP.count(initial_rtype) == 0){
                return std::vector<int>();
            }

            // all other cases we process the subsequent values
            const int rtyp_as_int = static_cast<int>(rtyp*_rtyp_factor);
            const int expanded_rtype = rtyp_as_int - (initial_rtype*_rtyp_factor);
            std::list<int> rtypes = util::get_digits(expanded_rtype);
            // we also need to remove trailing zeros,
            // as this does not mean GammaRay emission
            // it is the end
            std::vector<int> reaction_types{initial_rtype};
            for(int v: rtypes){
                if((v == 0) || (DECAY_MAP.count(v) == 0)){
                    break;
                }
                reaction_types.emplace_back(v);
            }
            return reaction_types;
        }

        DecayProductT compute_products(int parent_zai, double rtyp, int final_isomeric_state){
            const std::vector<int> rtypes = expand_rtypes(rtyp);

            // we need to trim remaining zeros, as there is no way to tell between
            // 1.0 (beta -> gamma) and 1.0 (beta)
            DecayProductComputer dpc = ToSinkDecayReactionProducts();
            int primary_zai = parent_zai;
            std::vector<int> secondaries;
            for(auto reaction_type: rtypes){
                if (DECAY_MAP.count(reaction_type) == 1){
                    dpc = DECAY_MAP.at(reaction_type);
                    // we assume - and because the data tells us this -
                    // that SF (6) doesn't have any modes following it
                    // to make this simpler with the primary
                    // Perhaps a tree is a better data structure here....
                    DecayProductT prods = std::visit([&primary_zai](auto const& obj) {
                        return obj.get_products(primary_zai);
                    }, dpc);
                    primary_zai = std::get<0>(prods);
                    for(auto& secondary: std::get<1>(prods)){
                        secondaries.emplace_back(secondary);
                    }
                }
            }
            // finally check the final isomeric state
            const zai::triple primary_triple = zai::decompose(primary_zai);
            primary_zai = zai::make(std::get<0>(primary_triple),
                                    std::get<1>(primary_triple),
                                    final_isomeric_state);

            return DecayProductT(primary_zai, secondaries);
        }

    }
}

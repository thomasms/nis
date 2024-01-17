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

            // if we are zero first then we are gamma
            const int initial_rtype = static_cast<int>(std::floor(rtyp));
            if(initial_rtype == 0){
                return std::vector<int>{0};
            }

            const int expanded_rtype = static_cast<int>(rtyp*1'000'000);
            std::list<int> rtypes = util::get_digits(expanded_rtype);
            // this is to account for initial procedures bigger than 9
            rtypes.pop_front();
            rtypes.push_front(initial_rtype);
            // we also need to remove trailing zeros, as this does not mean GammaRay emission
            // it is the end
            std::vector<int> reaction_types;
            for(int v: rtypes){
                if((v == 0) || (DECAY_MAP.count(v) == 0)){
                    break;
                }
                reaction_types.push_back(v);
            }
            return reaction_types;
        }

        DecayProductT compute_products(int parent_zai, double rtyp, int final_isomeric_state){
            using Comp_ptr_t = std::unique_ptr<DecayProductComputer>;
            const std::vector<int> rtypes = expand_rtypes(rtyp);

            // we need to trim remaining zeros, as there is no way to tell between
            // 1.0 (beta -> gamma) and 1.0 (beta)
            Comp_ptr_t dpc = std::make_unique<ToSinkDecayReactionProducts>();
            int primary_zai = parent_zai;
            std::vector<int> secondaries;
            for(auto reaction_type: rtypes){
                if (DECAY_MAP.count(reaction_type) == 1){
                    auto create_computer = DECAY_MAP.at(reaction_type);
                    dpc = create_computer();
                    // we assume - and because the data tells us this -
                    // that SF (6) doesn't have any modes following it
                    // to make this simpler with the primary
                    // Perhaps a tree is a better data structure here....
                    DecayProductT prods = dpc->get_products(primary_zai);
                    primary_zai = std::get<0>(prods);
                    for(auto secondary: std::get<1>(prods)){
                        secondaries.push_back(secondary);
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

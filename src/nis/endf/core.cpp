#include <cmath>
#include <memory>

#include "nis/physics.hpp"
#include "nis/nuclides.hpp"
#include "nis/util.hpp"
#include "nis/endf/decay.hpp"

#include "nis/endf/core.hpp"

namespace nis {
    namespace endf {

        int get_parent_zai(const njoy::ENDFtk::tree::Material& material){
            const int mf_number = 1;
            const int mt_number = 451;
            if(material.hasMF(mf_number)){
                njoy::ENDFtk::tree::File file = material.file(mf_number);
                if(file.hasMT(mt_number)){
                    njoy::ENDFtk::tree::Section section = file.MT(mt_number);
                    auto parsed_section = section.parse<mf_number, mt_number>();
                    const int za = parsed_section.ZA();
                    // ENDFtk gives us the isomeric state as an unsigned int
                    const int i = static_cast<int>(parsed_section.LISO());
                    const int zai = za*10 + i;
                    return zai;
                }
            }
        }

    }
}

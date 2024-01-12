#include <cmath>


namespace nis {
    namespace zai {

        int get_z_from_zai(int zai){
            return std::floor(zai / 10000);
        }

        int get_a_from_zai(int zai){
            int za = std::floor(zai / 10);
            return za % 1000;
        }

        int get_i_from_zai(int zai){
            return zai % 10;
        }

    }
}

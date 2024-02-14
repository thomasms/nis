#ifndef ZAI_HPP
#define ZAI_HPP

#include <cmath>
#include <map>
#include <tuple>


namespace nis {
    namespace zai {

        using triple = std::tuple<int, int, int>;

        const int ALPHA = 20040;
        // would be better to make at compile time
        const std::map<int, unsigned char> META_STABLE_MAP = {
            {0, ' '},
            {1, 'm'},
            {2, 'n'},
            {3, 'o'},
            {4, 'p'},
            {5, 'q'},
            {6, 'r'},
            {7, 's'},
            {8, 't'},
            {9, 'u'},
        };

        inline int get_z(int zai){
            return std::floor(zai / 10000);
        }

        inline int get_a(int zai){
            int za = std::floor(zai / 10);
            return za % 1000;
        }

        inline int get_i(int zai){
            return zai % 10;
        }

        inline int make(int z, int a, int i){
            return z*10000 + a*10 + i;
        }

        inline int make(const triple& t){
            return make(std::get<0>(t), std::get<1>(t), std::get<2>(t));
        }

        inline triple decompose(int zai){
            const int z = get_z(zai);
            const int a = get_a(zai);
            const int i = get_i(zai);
            return triple(z, a, i);
        }

        int get_from_name(std::string name);
        std::string get_name_from(int zai);

    }
}

#endif // ZAI_HPP
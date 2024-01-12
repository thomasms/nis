#ifndef UTIL_HPP
#define UTIL_HPP

#include <algorithm>
#include <map>

namespace nis {
    namespace util {

        template<std::size_t N>
        std::map<int, int> make_map_index_lookup(const std::array<int, N>& arr){
            std::map<int, int> hash_map;
            std::transform(arr.begin(), arr.end(), std::inserter(hash_map, hash_map.begin()),
                                [i = 0](int num) mutable { return std::make_pair(num, i++); });
            return hash_map;
        }
    }
}

#endif // UTIL_HPP
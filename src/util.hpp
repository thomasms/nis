#ifndef UTIL_HPP
#define UTIL_HPP

#include <algorithm>
#include <list>
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

        /*
        Sign is removed/ignored.
        */
        template<typename T>
        std::list<int> get_digits(T number){
            std::list<int> digits;
            if(number == 0){
                digits.push_back(0);
                return digits;
            }
            const T sign_removed = (std::is_signed_v<T> && number < 0) ? std::max(number, -number) : number;
            T value = sign_removed;
            while(value > 0){
                int digit = value % 10;
                digits.push_back(digit);
                value /= 10;
            }
            std::reverse(digits.begin(), digits.end());
            return digits;
        }
    }
}

#endif // UTIL_HPP
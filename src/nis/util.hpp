#ifndef UTIL_HPP
#define UTIL_HPP

#include <algorithm>
#include <array>
#include <cctype>
#include <iostream>
#include <list>
#include <map>
#include <string>


namespace nis {
    namespace util {

        template<typename Key, typename Value, std::size_t N>
        struct CompileTimeMap {
            std::array<std::pair<Key, Value>, N> pairs;

            static constexpr Value defaultValue{};
            constexpr const Value& operator[](const Key& key) const {
                auto it = std::find_if(pairs.begin(), pairs.end(),
                    [&key](const auto& pair) { return pair.first == key; });

                if (it != pairs.end()) {
                    return it->second;
                }

                // Return a default value if key is not found
                return defaultValue;
            }
        };

        template<std::size_t N>
        inline std::map<int, int> make_map_index_lookup(const std::array<int, N>& arr){
            std::map<int, int> hash_map;
            std::transform(arr.begin(), arr.end(), std::inserter(hash_map, hash_map.begin()),
                                [i = 0](int num) mutable { return std::make_pair(num, i++); });
            return hash_map;
        }

        // if values are not unique, it will take the first
        template<typename K, typename V>
        inline std::map<V, K> invert_map(const std::map<K, V>& imap){
            std::map<V, K> reversed;
            for (const auto& pair : imap) {
                reversed.emplace(pair.second, pair.first);
            }
            return reversed;
        }

        /*
        Sign is removed/ignored.
        */
        template<typename T>
        inline std::list<int> get_digits(T number){
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

        inline std::string string_lower(const std::string& str) {
            std::string result = str;
            for(char& c : result) {
                c = std::tolower(static_cast<unsigned char>(c));
            }
            return result;
        }

        inline std::string string_rtrim(std::string str){
            size_t end = str.find_last_not_of(" \t\r\n");
            return (end == std::string::npos) ? "" : str.substr(0, end + 1);
        }

    }
}

#endif // UTIL_HPP
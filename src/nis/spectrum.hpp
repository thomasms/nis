#ifndef SPECTRUM_HPP
#define SPECTRUM_HPP

#include <algorithm>
#include <cmath>


namespace nis {
    namespace core {

        template<typename T=double>
        using SpectrumT = std::vector<std::tuple<T, T>>;
        
        template<typename T=double>
        struct Edge{
            explicit Edge(T lower, T upper) : _lower(lower), _upper(upper) {
                assert(lower <= upper);
            }

            bool operator==(const Edge& rhs){
                return (_lower == rhs.lower()) && (_upper == rhs.upper());
            };

            T width() const {
                return this->upper - this->_lower;
            };

            T log_ratio() const {
                if(_lower <= 0){
                    return 1;
                }
                return std::log(_upper/_lower);
            }

            T lower() const { return _lower;};
            T upper() const { return _upper;};

            private:
                T _lower;
                T _upper;
        };

        template <typename T>
        std::tuple<bool, Edge<T>> get_overlap(const Edge<T>& edge1, const Edge<T>& edge2){
            const bool IS_LOWER = true;
            const bool IS_UPPER = false;
            std::vector<std::tuple<T, bool>> points = {
                std::make_pair(edge1.lower(), IS_LOWER),
                std::make_pair(edge1.upper(), IS_UPPER),
                std::make_pair(edge2.lower(), IS_LOWER),
                std::make_pair(edge2.upper(), IS_UPPER),
            };

            std::sort(points.begin(), points.end(), [](std::tuple<T, bool>& a, std::tuple<T, bool>& b)
            {
                // ascending
                return std::get<0>(a) < std::get<0>(b);
            });

            // the second item must be a lower point and the third item must be an upper point
            // for overlap
            const bool overlap = (std::get<1>(points[1]) == IS_LOWER) && (std::get<1>(points[2]) == IS_UPPER);
            const Edge overlap_edge = Edge(std::get<0>(points[1]), std::get<0>(points[2]));

            return std::make_pair(overlap, overlap_edge);
        }


    }
}

#endif // SPECTRUM_HPP
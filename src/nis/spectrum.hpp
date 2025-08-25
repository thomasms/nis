#ifndef SPECTRUM_HPP
#define SPECTRUM_HPP

#include <algorithm>
#include <cmath>
#include <functional>

#include "nis/histogram.hpp"


namespace nis {
    namespace core {

        template<typename T=double>
        struct Edge{
            explicit Edge(T lower, T upper) : _lower(lower), _upper(upper) {
                assert(lower <= upper);
            }

            bool operator==(const Edge& rhs){
                return (_lower == rhs.lower()) && (_upper == rhs.upper());
            };

            T width() const {
                return _upper - _lower;
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

        template<typename T=double>
        std::vector<Edge<T>> _get_edges_from_bounds(const std::vector<T>& bounds){
            std::vector<Edge<T>> edges;
            edges.reserve(bounds.size()-1);
            for(unsigned int i=0;i<bounds.size()-1;++i){
                edges.emplace_back(Edge(bounds[i], bounds[i+1]));
            }
            return edges;
        }

        template<typename T=double>
        using ConvertFuncT = std::function<T(const Edge<T>& input_edge, T input_value, const Edge<T>& overlapping_edge)>;

        template<typename T=double>
        std::tuple<T, unsigned int> _compute_overlap(const std::vector<Edge<T>>& input_edges,
                                                     const std::vector<T>& input_values,
                                                     const Edge<T>& oedge,
                                                     ConvertFuncT<T>& cfunc,
                                                     unsigned int last_overlap_index=0){
            T output_value = 0;
            bool prev_has_overlap = false;
            int new_last_index = last_overlap_index;
            for(unsigned int i=last_overlap_index; i<input_edges.size();++i){
                std::tuple<bool, Edge<T>> overlap = get_overlap(input_edges[i], oedge);
                bool has_overlap = std::get<0>(overlap);
                const auto overlap_edge = std::get<1>(overlap);
                if(has_overlap){
                    new_last_index = i;
                    output_value += cfunc(input_edges[i],
                                          input_values[i],
                                          overlap_edge);
                }

                if(!has_overlap && prev_has_overlap){
                    break;
                }

                prev_has_overlap = has_overlap;
            }

            return std::make_pair(output_value, new_last_index);
        }

        template<typename T=double>
        using Spectrum = Histogram<std::vector, T, T>;

        template<typename T=double>
        Spectrum<T> _convert_imp(const Spectrum<T>& in_spectrum,
                                 const std::vector<T>& output_bounds,
                                 ConvertFuncT<T>& cfunc){
            /*
            Returns the output_values depending on the cfunc.

            output_bounds is a list of energies, units are irrelevant,
            as long as it matches the units of the input_bounds.

            Asserts both input and output bounds are of length greater than 1

            Assumes that input and output bounds are in ascending energy. If not
            then it will go unchecked and will produce odd results
            */
            assert(output_bounds.size() > 1);

            const auto input_edges = _get_edges_from_bounds<T>(in_spectrum.x());
            const auto output_edges = _get_edges_from_bounds<T>(output_bounds);

            unsigned int last_index = 0;
            std::vector<T> output_values;
            output_values.reserve(output_edges.size());
            for(auto& edge: output_edges){
                std::tuple<T, unsigned int> x = _compute_overlap(input_edges, in_spectrum.y(), edge, cfunc, last_index);
                last_index = std::get<1>(x);
                output_values.emplace_back(std::get<0>(x));
            }

            return Spectrum<T>(output_bounds, std::move(output_values));
        }

        template<typename T=double>
        Spectrum<T> convert_by_energy(const Spectrum<T>& in_spectrum,
                                         const std::vector<T>& output_bounds){
            ConvertFuncT<T> cfunc = [](const Edge<T>& input_edge, T input_value, const Edge<T>& overlapping_edge){
                return overlapping_edge.width()*input_value/input_edge.width();
            };
            return _convert_imp(in_spectrum, output_bounds, cfunc);
        }

        template<typename T=double>
        Spectrum<T> convert_by_lethargy(const Spectrum<T>& in_spectrum,
                                           const std::vector<T>& output_bounds){
            ConvertFuncT<T> cfunc = [](const Edge<T>& input_edge, T input_value, const Edge<T>& overlapping_edge){
                return overlapping_edge.log_ratio()*input_value/input_edge.log_ratio();
            };
            return _convert_imp(in_spectrum, output_bounds, cfunc);
        }

    }
}

#endif // SPECTRUM_HPP
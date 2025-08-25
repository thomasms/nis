#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include <algorithm>
#include <cmath>
#include <functional>


namespace nis {
    namespace core {

        // TODO: associate a unit with the energies and values
        // we could use a library for units, but not willing to adopt another dependency
        // boost::units - too heavy
        // mp-units - is nice and modern but c++20 only
        // nholthaus/units - could work maybe
        // aurora-opensource/au - looks good but not very popular and still early days
        template <template <typename> class Container, typename xT=double, typename yT=double>
        struct Histogram{
            explicit Histogram(const Container<xT>& x, const Container<yT>& y) : _x(x), _y(y) {
                assert(_x.size() > 1);
                assert(_y.size() + 1 <= _x.size());
            }

            const Container<xT>& x() const { return _x;};
            const Container<yT>& y() const { return _y;};

            private:
                Container<xT> _x;
                Container<yT> _y;
        };

    }
}

#endif //HISTOGRAM_HPP
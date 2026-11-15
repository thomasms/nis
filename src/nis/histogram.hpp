#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>


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
                // exceptions, not asserts: a spectrum can come from a user, and an assert stops the whole process
                if(_x.size() < 2){
                    throw std::invalid_argument("a histogram needs at least 2 bounds");
                }
                if(_y.size() + 1 > _x.size()){
                    throw std::invalid_argument("a histogram has more values than groups");
                }
            }

            Container<xT> xcopy() const { return _x;};
            Container<yT> ycopy() const { return _y;};

            const Container<xT>& x() const { return _x;};
            const Container<yT>& y() const { return _y;};

            private:
                Container<xT> _x;
                Container<yT> _y;
        };

    }
}

#endif //HISTOGRAM_HPP
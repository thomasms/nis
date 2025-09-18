#ifndef MATHS_HPP
#define MATHS_HPP

#include <cmath>

namespace nis{
    namespace maths {

        // Should we just use Eigen here?
        // TODO: accommodate units
        template<typename T=double>
        T dot_product(const std::vector<T>& vec1, const std::vector<T>& vec2){
            assert(vec1.size() == vec2.size());

            T d_product = 0;
            for(unsigned int i=0;i<vec1.size();++i){
                d_product += vec1[i]*vec2[i];
            }
            return d_product;
        };

    }
}

#endif //MATHS_HPP
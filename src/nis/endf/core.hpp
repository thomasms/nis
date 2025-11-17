#ifndef ENDF_CORE_HPP
#define ENDF_CORE_HPP

#include <ENDFtk.hpp>

#include <Eigen/Core>


namespace nis {
    namespace endf {

        int get_parent_zai(const njoy::ENDFtk::tree::Material& material);

    }
}

#endif //ENDF_CORE_HPP
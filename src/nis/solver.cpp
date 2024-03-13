#include "nis/cram.hpp"
#include "nis/nuclides.hpp"
#include "nis/util.hpp"
#include "nis/zai.hpp"

#include "nis/solver.hpp"


namespace nis {
    namespace solver {

        Eigen::VectorX<double> compute_inventory(const io::input::InputData<double>& input_data, const Eigen::MatrixXd& mat){
            // set initial vector (inventory)
            const int vec_size = static_cast<int>(mat.cols());
            Eigen::VectorX<double> n0 = Eigen::VectorX<double>::Zero(vec_size);
            const auto zai_map = util::make_map_index_lookup(data::ALL_ZAIS);

            for (auto it = input_data.nuclides.begin(); it != input_data.nuclides.end(); ++it) {
                const int nuclide_zai = zai::get_from_name(it->first);
                if(zai_map.count(nuclide_zai) == 1){
                    const auto zai_index = zai_map.at(nuclide_zai);
                    n0(zai_index) = it->second;
                }
            }

            // solve
            const Eigen::VectorX<double> nt = core::cram48::solve<double>(mat, n0, input_data.dt);
            return nt;

        }

    }
}
#include <cmath>


namespace nis{
    namespace physics {

        const double LN_2 = std::log(2);

        double compute_decay_const(double half_life){
            return LN_2/half_life;
        }

    }
}
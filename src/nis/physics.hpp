#ifndef PHYSICS_HPP
#define PHYSICS_HPP

#include <cmath>


namespace nis{
    namespace physics {

        const double LN_2 = std::log(2);

        // exact, from the 2019 SI definition of the elementary charge
        constexpr double EV_TO_J = 1.602176634e-19;

        inline double compute_decay_const(double half_life){
            return LN_2/half_life;
        }

    }
}

#endif //PHYSICS_HPP
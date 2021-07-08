#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "tmp0n6wklye_p.h"
#include "tmp0n6wklye_k.h"
#include "tmp0n6wklye_y.h"
#include "tmp0n6wklye_sigmay.h"
#include "tmp0n6wklye_my.h"

namespace amici {
namespace model_tmp0n6wklye {

void Jy_tmp0n6wklye(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_deaths, 2)) + 0.5*std::pow(-mobservable_deaths + observable_deaths, 2)/std::pow(sigma_observable_deaths, 2);
            break;
    }
}

} // namespace model_tmp0n6wklye
} // namespace amici

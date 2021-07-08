#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "tmp2t76h1f9_p.h"
#include "tmp2t76h1f9_k.h"
#include "tmp2t76h1f9_y.h"
#include "tmp2t76h1f9_sigmay.h"
#include "tmp2t76h1f9_my.h"

namespace amici {
namespace model_tmp2t76h1f9 {

void Jy_tmp2t76h1f9(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_deaths, 2)) + 0.5*std::pow(-mobservable_deaths + observable_deaths, 2)/std::pow(sigma_observable_deaths, 2);
            break;
    }
}

} // namespace model_tmp2t76h1f9
} // namespace amici

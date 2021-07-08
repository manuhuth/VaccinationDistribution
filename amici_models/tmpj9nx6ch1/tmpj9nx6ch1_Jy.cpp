#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

namespace amici {
namespace model_tmpj9nx6ch1 {

void Jy_tmpj9nx6ch1(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(std::pow(sigmaobservable_deaths, 2)) + 0.5*M_LN2 + 0.5*std::log(amici::pi) + 0.5*std::pow(-mobservable_deaths + observable_deaths, 2)/std::pow(sigmaobservable_deaths, 2);
            break;
    }
}

} // namespace amici
} // namespace model_tmpj9nx6ch1
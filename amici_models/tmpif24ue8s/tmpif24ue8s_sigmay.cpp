#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "p.h"
#include "k.h"

namespace amici {
namespace model_tmpif24ue8s {

void sigmay_tmpif24ue8s(realtype *sigmay, const realtype t, const realtype *p, const realtype *k){
    sigmay[0] = 1.0;
}

} // namespace amici
} // namespace model_tmpif24ue8s
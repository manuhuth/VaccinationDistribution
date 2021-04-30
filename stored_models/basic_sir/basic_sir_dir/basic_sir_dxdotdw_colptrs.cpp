#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_basic_sir {

static constexpr std::array<sunindextype, 3> dxdotdw_colptrs_basic_sir_ = {
    0, 2, 4
};

void dxdotdw_colptrs_basic_sir(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_basic_sir_));
}
} // namespace model_basic_sir
} // namespace amici

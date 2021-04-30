#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_basic_sir {

static constexpr std::array<sunindextype, 4> dxdotdw_rowvals_basic_sir_ = {
    0, 1, 1, 2
};

void dxdotdw_rowvals_basic_sir(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_basic_sir_));
}
} // namespace model_basic_sir
} // namespace amici

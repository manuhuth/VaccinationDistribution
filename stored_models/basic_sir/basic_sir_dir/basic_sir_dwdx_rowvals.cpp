#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_basic_sir {

static constexpr std::array<sunindextype, 4> dwdx_rowvals_basic_sir_ = {
    0, 0, 1, 0
};

void dwdx_rowvals_basic_sir(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_basic_sir_));
}
} // namespace model_basic_sir
} // namespace amici

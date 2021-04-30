#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_basic_sir {

static constexpr std::array<sunindextype, 4> dwdx_colptrs_basic_sir_ = {
    0, 1, 3, 4
};

void dwdx_colptrs_basic_sir(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_basic_sir_));
}
} // namespace model_basic_sir
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_basic_sir {

static constexpr std::array<sunindextype, 3> dwdp_colptrs_basic_sir_ = {
    0, 1, 2
};

void dwdp_colptrs_basic_sir(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_basic_sir_));
}
} // namespace model_basic_sir
} // namespace amici

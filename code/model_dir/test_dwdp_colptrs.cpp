#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_test {

static constexpr std::array<sunindextype, 3> dwdp_colptrs_test_ = {
    0, 1, 2
};

void dwdp_colptrs_test(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_test_));
}
} // namespace model_test
} // namespace amici

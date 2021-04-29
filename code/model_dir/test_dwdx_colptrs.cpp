#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_test {

static constexpr std::array<sunindextype, 4> dwdx_colptrs_test_ = {
    0, 1, 3, 4
};

void dwdx_colptrs_test(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_test_));
}
} // namespace model_test
} // namespace amici

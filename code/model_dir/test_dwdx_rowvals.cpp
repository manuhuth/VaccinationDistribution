#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_test {

static constexpr std::array<sunindextype, 4> dwdx_rowvals_test_ = {
    0, 0, 1, 0
};

void dwdx_rowvals_test(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_test_));
}
} // namespace model_test
} // namespace amici

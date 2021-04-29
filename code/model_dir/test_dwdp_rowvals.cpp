#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_test {

static constexpr std::array<sunindextype, 2> dwdp_rowvals_test_ = {
    0, 1
};

void dwdp_rowvals_test(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_test_));
}
} // namespace model_test
} // namespace amici

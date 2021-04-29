#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_test {

static constexpr std::array<sunindextype, 4> dxdotdw_rowvals_test_ = {
    0, 1, 1, 2
};

void dxdotdw_rowvals_test(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_test_));
}
} // namespace model_test
} // namespace amici

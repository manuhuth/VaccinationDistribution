#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_test {

static constexpr std::array<sunindextype, 3> dxdotdw_colptrs_test_ = {
    0, 2, 4
};

void dxdotdw_colptrs_test(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_test_));
}
} // namespace model_test
} // namespace amici

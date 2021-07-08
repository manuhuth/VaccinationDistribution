#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_tmp2t76h1f9 {

static constexpr std::array<sunindextype, 24> dwdw_rowvals_tmp2t76h1f9_ = {
    4, 9, 10, 8, 11, 5, 7, 8, 6, 9, 10, 11, 108, 112, 114, 109, 113, 115, 111, 117, 119, 110, 116, 118
};

void dwdw_rowvals_tmp2t76h1f9(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_tmp2t76h1f9_));
}
} // namespace model_tmp2t76h1f9
} // namespace amici

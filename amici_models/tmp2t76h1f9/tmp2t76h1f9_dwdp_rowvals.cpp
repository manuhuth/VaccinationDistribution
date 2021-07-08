#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_tmp2t76h1f9 {

static constexpr std::array<sunindextype, 22> dwdp_rowvals_tmp2t76h1f9_ = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3
};

void dwdp_rowvals_tmp2t76h1f9(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_tmp2t76h1f9_));
}
} // namespace model_tmp2t76h1f9
} // namespace amici

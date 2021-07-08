#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_tmp2t76h1f9 {

static constexpr std::array<sunindextype, 23> dwdp_colptrs_tmp2t76h1f9_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22
};

void dwdp_colptrs_tmp2t76h1f9(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_tmp2t76h1f9_));
}
} // namespace model_tmp2t76h1f9
} // namespace amici

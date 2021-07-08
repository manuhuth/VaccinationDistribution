#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_tmp0n6wklye {

static constexpr std::array<sunindextype, 22> dwdp_rowvals_tmp0n6wklye_ = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3
};

void dwdp_rowvals_tmp0n6wklye(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_tmp0n6wklye_));
}
} // namespace model_tmp0n6wklye
} // namespace amici

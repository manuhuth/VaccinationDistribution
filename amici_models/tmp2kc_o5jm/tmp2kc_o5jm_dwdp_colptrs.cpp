#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_tmp2kc_o5jm {

static constexpr std::array<int, 23> dwdp_colptrs_tmp2kc_o5jm_ = {
    0, 11, 22, 33, 44, 55, 66, 77, 88, 99, 110, 121, 132, 143, 154, 165, 176, 187, 198, 209, 220, 231, 242
};

void dwdp_colptrs_tmp2kc_o5jm(sunindextype *colptrs){
    std::copy(dwdp_colptrs_tmp2kc_o5jm_.begin(), dwdp_colptrs_tmp2kc_o5jm_.end(), colptrs);
}
} // namespace amici
} // namespace model_tmp2kc_o5jm
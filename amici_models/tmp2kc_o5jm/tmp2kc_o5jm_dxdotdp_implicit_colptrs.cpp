#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_tmp2kc_o5jm {

static constexpr std::array<int, 23> dxdotdp_implicit_colptrs_tmp2kc_o5jm_ = {
    0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440
};

void dxdotdp_implicit_colptrs_tmp2kc_o5jm(sunindextype *colptrs){
    std::copy(dxdotdp_implicit_colptrs_tmp2kc_o5jm_.begin(), dxdotdp_implicit_colptrs_tmp2kc_o5jm_.end(), colptrs);
}
} // namespace amici
} // namespace model_tmp2kc_o5jm
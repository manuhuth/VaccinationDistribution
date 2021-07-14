#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination_piecewise {

static constexpr std::array<int, 111> dxdotdp_implicit_colptrs_vaccination_piecewise_ = {
    0, 36, 68, 86, 104, 108, 112, 116, 120, 124, 128, 132, 136, 148, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 560, 560, 560, 560, 560, 560, 560, 560, 560, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 820, 840, 860, 880, 900, 920, 940, 960, 960, 978, 996, 996
};

void dxdotdp_implicit_colptrs_vaccination_piecewise(sunindextype *colptrs){
    std::copy(dxdotdp_implicit_colptrs_vaccination_piecewise_.begin(), dxdotdp_implicit_colptrs_vaccination_piecewise_.end(), colptrs);
}
} // namespace amici
} // namespace model_vaccination_piecewise
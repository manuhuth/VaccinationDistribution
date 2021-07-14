#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<int, 114> dxdotdp_implicit_colptrs_vaccination_ = {
    0, 36, 68, 86, 104, 108, 112, 116, 120, 124, 128, 132, 136, 148, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 560, 578, 596, 596, 626, 656, 686, 716, 746, 776, 806, 836, 866, 896, 926, 946, 966, 986, 1006, 1026, 1046, 1066, 1086, 1106, 1126, 1146, 1166, 1186, 1206, 1226, 1246, 1266, 1286, 1306, 1326, 1346, 1366
};

void dxdotdp_implicit_colptrs_vaccination(sunindextype *colptrs){
    std::copy(dxdotdp_implicit_colptrs_vaccination_.begin(), dxdotdp_implicit_colptrs_vaccination_.end(), colptrs);
}
} // namespace amici
} // namespace model_vaccination
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<int, 114> dwdp_colptrs_vaccination_ = {
    0, 24, 48, 96, 168, 172, 184, 188, 200, 204, 216, 220, 232, 268, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 314, 324, 334, 344, 354, 364, 374, 384, 394, 404, 414, 424, 434, 444, 454, 464, 474, 484, 494, 504, 504, 540, 576, 576, 598, 620, 642, 664, 686, 708, 730, 752, 774, 796, 818, 829, 840, 851, 862, 873, 884, 895, 906, 917, 928, 939, 950, 961, 972, 983, 994, 1005, 1016, 1027, 1038, 1049, 1060
};

void dwdp_colptrs_vaccination(sunindextype *colptrs){
    std::copy(dwdp_colptrs_vaccination_.begin(), dwdp_colptrs_vaccination_.end(), colptrs);
}
} // namespace amici
} // namespace model_vaccination
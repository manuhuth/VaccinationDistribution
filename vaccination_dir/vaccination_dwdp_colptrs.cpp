#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<int, 117> dwdp_colptrs_vaccination_ = {
    0, 24, 48, 96, 168, 190, 194, 206, 210, 222, 226, 238, 242, 254, 290, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 336, 346, 356, 366, 376, 386, 396, 406, 416, 426, 426, 436, 446, 456, 466, 476, 486, 496, 506, 516, 526, 526, 526, 562, 598, 598, 620, 642, 664, 686, 708, 730, 752, 774, 796, 818, 840, 851, 862, 873, 884, 895, 906, 917, 928, 939, 950, 961, 972, 983, 994, 1005, 1016, 1027, 1038, 1049, 1060, 1071, 1082
};

void dwdp_colptrs_vaccination(sunindextype *colptrs){
    std::copy(dwdp_colptrs_vaccination_.begin(), dwdp_colptrs_vaccination_.end(), colptrs);
}
} // namespace amici
} // namespace model_vaccination
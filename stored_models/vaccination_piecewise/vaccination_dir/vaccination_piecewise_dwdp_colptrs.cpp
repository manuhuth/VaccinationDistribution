#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination_piecewise {

static constexpr std::array<int, 111> dwdp_colptrs_vaccination_piecewise_ = {
    0, 24, 48, 96, 168, 172, 184, 188, 200, 204, 216, 220, 232, 268, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 314, 324, 334, 344, 354, 364, 374, 384, 394, 404, 414, 424, 434, 444, 454, 464, 474, 484, 494, 504, 504, 504, 504, 504, 504, 504, 504, 504, 504, 504, 514, 524, 534, 544, 554, 564, 574, 584, 594, 604, 614, 624, 634, 644, 654, 664, 674, 684, 694, 704, 704, 740, 776, 776
};

void dwdp_colptrs_vaccination_piecewise(sunindextype *colptrs){
    std::copy(dwdp_colptrs_vaccination_piecewise_.begin(), dwdp_colptrs_vaccination_piecewise_.end(), colptrs);
}
} // namespace amici
} // namespace model_vaccination_piecewise
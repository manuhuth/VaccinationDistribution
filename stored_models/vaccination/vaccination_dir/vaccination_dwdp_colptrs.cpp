#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<int, 114> dwdp_colptrs_vaccination_ = {
    0, 24, 48, 96, 168, 172, 184, 188, 200, 204, 216, 220, 232, 268, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 314, 324, 334, 344, 354, 364, 374, 384, 394, 394, 404, 414, 424, 434, 444, 454, 464, 474, 484, 484, 484, 520, 556, 556, 578, 600, 622, 644, 666, 688, 710, 732, 754, 776, 798, 809, 820, 831, 842, 853, 864, 875, 886, 897, 908, 919, 930, 941, 952, 963, 974, 985, 996, 1007, 1018, 1029, 1040
};

void dwdp_colptrs_vaccination(sunindextype *colptrs){
    std::copy(dwdp_colptrs_vaccination_.begin(), dwdp_colptrs_vaccination_.end(), colptrs);
}
} // namespace amici
} // namespace model_vaccination
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<int, 117> dxdotdp_implicit_colptrs_vaccination_ = {
    0, 36, 68, 86, 104, 134, 138, 142, 146, 150, 154, 158, 162, 166, 178, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 210, 230, 250, 270, 290, 310, 330, 350, 370, 390, 390, 410, 430, 450, 470, 490, 510, 530, 550, 570, 590, 590, 590, 608, 626, 626, 656, 686, 716, 746, 776, 806, 836, 866, 896, 926, 956, 976, 996, 1016, 1036, 1056, 1076, 1096, 1116, 1136, 1156, 1176, 1196, 1216, 1236, 1256, 1276, 1296, 1316, 1336, 1356, 1376, 1396
};

void dxdotdp_implicit_colptrs_vaccination(sunindextype *colptrs){
    std::copy(dxdotdp_implicit_colptrs_vaccination_.begin(), dxdotdp_implicit_colptrs_vaccination_.end(), colptrs);
}
} // namespace amici
} // namespace model_vaccination
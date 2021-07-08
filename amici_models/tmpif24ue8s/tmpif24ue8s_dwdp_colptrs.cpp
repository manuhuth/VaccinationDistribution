#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_tmpif24ue8s {

static constexpr std::array<int, 23> dwdp_colptrs_tmpif24ue8s_ = {
    0, 11, 22, 33, 44, 55, 66, 77, 88, 99, 110, 121, 132, 143, 154, 165, 176, 187, 198, 209, 220, 231, 242
};

void dwdp_colptrs_tmpif24ue8s(sunindextype *colptrs){
    std::copy(dwdp_colptrs_tmpif24ue8s_.begin(), dwdp_colptrs_tmpif24ue8s_.end(), colptrs);
}
} // namespace amici
} // namespace model_tmpif24ue8s
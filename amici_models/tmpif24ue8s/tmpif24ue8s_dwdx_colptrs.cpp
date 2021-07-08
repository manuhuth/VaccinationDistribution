#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_tmpif24ue8s {

static constexpr std::array<int, 43> dwdx_colptrs_tmpif24ue8s_ = {
    0, 82, 154, 226, 308, 380, 452, 536, 620, 694, 768, 842, 916, 1000, 1084, 1158, 1232, 1306, 1380, 1462, 1544, 1616, 1688, 1760, 1832, 1914, 1996, 2068, 2140, 2212, 2284, 2284, 2284, 2284, 2284, 2284, 2284, 2284, 2284, 2284, 2284, 2284, 2284
};

void dwdx_colptrs_tmpif24ue8s(sunindextype *colptrs){
    std::copy(dwdx_colptrs_tmpif24ue8s_.begin(), dwdx_colptrs_tmpif24ue8s_.end(), colptrs);
}
} // namespace amici
} // namespace model_tmpif24ue8s
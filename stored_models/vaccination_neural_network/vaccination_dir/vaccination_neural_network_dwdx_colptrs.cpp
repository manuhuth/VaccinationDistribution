#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination_neural_network {

static constexpr std::array<int, 43> dwdx_colptrs_vaccination_neural_network_ = {
    0, 72, 144, 216, 288, 360, 432, 506, 580, 654, 728, 802, 876, 950, 1024, 1098, 1172, 1246, 1320, 1392, 1464, 1536, 1608, 1680, 1752, 1824, 1896, 1968, 2040, 2112, 2184, 2184, 2184, 2184, 2184, 2184, 2184, 2184, 2184, 2184, 2184, 2184, 2184
};

void dwdx_colptrs_vaccination_neural_network(sunindextype *colptrs){
    std::copy(dwdx_colptrs_vaccination_neural_network_.begin(), dwdx_colptrs_vaccination_neural_network_.end(), colptrs);
}
} // namespace amici
} // namespace model_vaccination_neural_network
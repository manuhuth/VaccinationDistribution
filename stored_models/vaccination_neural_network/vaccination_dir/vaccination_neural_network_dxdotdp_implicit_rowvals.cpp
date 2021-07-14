#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination_neural_network {

static constexpr std::array<int, 196> dxdotdp_implicit_rowvals_vaccination_neural_network_ = {
    6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 8, 9, 10, 11, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 20, 26, 32, 38, 1, 4, 8, 14, 22, 28, 34, 40, 2, 5, 10, 16, 21, 27, 33, 39, 1, 4, 9, 15, 23, 29, 35, 41, 2, 5, 11, 17, 0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 0, 1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17
};

void dxdotdp_implicit_rowvals_vaccination_neural_network(sunindextype *rowvals){
    std::copy(dxdotdp_implicit_rowvals_vaccination_neural_network_.begin(), dxdotdp_implicit_rowvals_vaccination_neural_network_.end(), rowvals);
}
} // namespace amici
} // namespace model_vaccination_neural_network
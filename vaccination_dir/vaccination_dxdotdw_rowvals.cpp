#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<int, 232> dxdotdw_rowvals_vaccination_ = {
    6, 18, 6, 30, 7, 19, 7, 31, 8, 20, 8, 32, 9, 21, 9, 33, 10, 22, 10, 34, 11, 23, 11, 35, 12, 24, 12, 36, 13, 25, 13, 37, 14, 26, 14, 38, 15, 27, 15, 39, 16, 28, 16, 40, 17, 29, 17, 41, 0, 6, 1, 8, 2, 10, 3, 12, 4, 14, 5, 16, 0, 7, 1, 9, 2, 11, 3, 13, 4, 15, 5, 17, 0, 6, 1, 8, 2, 10, 3, 12, 4, 14, 5, 16, 0, 7, 1, 9, 2, 11, 3, 13, 4, 15, 5, 17, 0, 6, 1, 8, 2, 10, 3, 12, 4, 14, 5, 16, 0, 7, 1, 9, 2, 11, 3, 13, 4, 15, 5, 17, 0, 6, 1, 8, 2, 10, 3, 12, 4, 14, 5, 16, 0, 7, 1, 9, 2, 11, 3, 13, 4, 15, 5, 17, 0, 6, 1, 8, 2, 10, 3, 12, 4, 14, 5, 16, 0, 7, 1, 9, 2, 11, 3, 13, 4, 15, 5, 17, 0, 6, 1, 8, 2, 10, 3, 12, 4, 14, 5, 16, 0, 7, 1, 9, 2, 11, 3, 13, 4, 15, 5, 17, 0, 1, 0, 2, 3, 4, 3, 5, 6, 8, 6, 10, 7, 9, 7, 11, 12, 14, 12, 16, 13, 15, 13, 17, 18, 20, 18, 22, 19, 21, 19, 23, 24, 26, 24, 28, 25, 27, 25, 29
};

void dxdotdw_rowvals_vaccination(sunindextype *rowvals){
    std::copy(dxdotdw_rowvals_vaccination_.begin(), dxdotdw_rowvals_vaccination_.end(), rowvals);
}
} // namespace amici
} // namespace model_vaccination
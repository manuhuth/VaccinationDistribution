#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<sunindextype, 216> dxdotdw_rowvals_vaccination_ = {
    6, 18, 6, 30, 7, 19, 7, 31, 8, 20, 8, 32, 9, 21, 9, 33, 10, 22, 10, 34, 11, 23, 11, 35, 12, 24, 12, 36, 13, 25, 13, 37, 14, 26, 14, 38, 15, 27, 15, 39, 16, 28, 16, 40, 17, 29, 17, 41, 0, 6, 1, 8, 2, 10, 3, 12, 4, 14, 5, 16, 0, 7, 1, 9, 2, 11, 3, 13, 4, 15, 5, 17, 0, 6, 1, 8, 2, 10, 3, 12, 4, 14, 5, 16, 0, 7, 1, 9, 2, 11, 3, 13, 4, 15, 5, 17, 0, 6, 1, 8, 2, 10, 3, 12, 4, 14, 5, 16, 0, 7, 1, 9, 2, 11, 3, 13, 4, 15, 5, 17, 0, 6, 1, 8, 2, 10, 3, 12, 4, 14, 5, 16, 0, 7, 1, 9, 2, 11, 3, 13, 4, 15, 5, 17, 0, 6, 1, 8, 2, 10, 3, 12, 4, 14, 5, 16, 0, 7, 1, 9, 2, 11, 3, 13, 4, 15, 5, 17, 0, 6, 1, 8, 2, 10, 3, 12, 4, 14, 5, 16, 0, 7, 1, 9, 2, 11, 3, 13, 4, 15, 5, 17, 0, 1, 0, 2, 3, 4, 3, 5, 18, 20, 18, 22, 19, 21, 19, 23, 24, 26, 24, 28, 25, 27, 25, 29
};

void dxdotdw_rowvals_vaccination(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_vaccination_));
}
} // namespace model_vaccination
} // namespace amici

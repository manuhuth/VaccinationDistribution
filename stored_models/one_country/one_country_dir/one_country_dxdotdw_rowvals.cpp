#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_one_country {

static constexpr std::array<sunindextype, 216> dxdotdw_rowvals_one_country_ = {
    6, 18, 6, 30, 7, 19, 7, 31, 8, 20, 8, 32, 9, 21, 9, 33, 10, 22, 10, 34, 11, 23, 11, 35, 12, 24, 12, 36, 13, 25, 13, 37, 14, 26, 14, 38, 15, 27, 15, 39, 16, 28, 16, 40, 17, 29, 17, 41, 0, 6, 1, 6, 2, 6, 3, 6, 4, 6, 5, 6, 0, 7, 1, 7, 2, 7, 3, 7, 4, 7, 5, 7, 0, 8, 1, 8, 2, 8, 3, 8, 4, 8, 5, 8, 0, 9, 1, 9, 2, 9, 3, 9, 4, 9, 5, 9, 0, 10, 1, 10, 2, 10, 3, 10, 4, 10, 5, 10, 0, 11, 1, 11, 2, 11, 3, 11, 4, 11, 5, 11, 0, 12, 1, 12, 2, 12, 3, 12, 4, 12, 5, 12, 0, 13, 1, 13, 2, 13, 3, 13, 4, 13, 5, 13, 0, 14, 1, 14, 2, 14, 3, 14, 4, 14, 5, 14, 0, 15, 1, 15, 2, 15, 3, 15, 4, 15, 5, 15, 0, 16, 1, 16, 2, 16, 3, 16, 4, 16, 5, 16, 0, 17, 1, 17, 2, 17, 3, 17, 4, 17, 5, 17, 0, 1, 0, 2, 3, 4, 3, 5, 18, 20, 18, 22, 19, 21, 19, 23, 24, 26, 24, 28, 25, 27, 25, 29
};

void dxdotdw_rowvals_one_country(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_one_country_));
}
} // namespace model_one_country
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<sunindextype, 72> dxdotdw_rowvals_vaccination_ = {
    3, 9, 3, 15, 4, 10, 4, 16, 5, 11, 5, 17, 6, 12, 6, 18, 7, 13, 7, 19, 8, 14, 8, 20, 0, 3, 1, 3, 2, 3, 0, 4, 1, 4, 2, 4, 0, 5, 1, 5, 2, 5, 0, 6, 1, 6, 2, 6, 0, 7, 1, 7, 2, 7, 0, 8, 1, 8, 2, 8, 0, 1, 0, 2, 9, 11, 9, 13, 10, 12, 10, 14
};

void dxdotdw_rowvals_vaccination(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_vaccination_));
}
} // namespace model_vaccination
} // namespace amici

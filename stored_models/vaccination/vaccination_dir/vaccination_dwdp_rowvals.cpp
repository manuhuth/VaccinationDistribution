#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<sunindextype, 96> dwdp_rowvals_vaccination_ = {
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 6, 7, 15, 21, 27, 10, 11, 16, 22, 28, 8, 9, 18, 24, 30, 12, 13, 19, 25, 31, 14, 15, 16, 20, 21, 22, 26, 27, 28, 17, 18, 19, 23, 24, 25, 29, 30, 31, 0, 1, 0, 1
};

void dwdp_rowvals_vaccination(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_vaccination_));
}
} // namespace model_vaccination
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_one_country {

static constexpr std::array<sunindextype, 98> dwdp_rowvals_one_country_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 4, 5, 13, 19, 25, 8, 9, 14, 20, 26, 6, 7, 16, 22, 28, 10, 11, 17, 23, 29, 12, 13, 14, 18, 19, 20, 24, 25, 26, 15, 16, 17, 21, 22, 23, 27, 28, 29, 30, 32, 34, 31, 33, 35
};

void dwdp_rowvals_one_country(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_one_country_));
}
} // namespace model_one_country
} // namespace amici

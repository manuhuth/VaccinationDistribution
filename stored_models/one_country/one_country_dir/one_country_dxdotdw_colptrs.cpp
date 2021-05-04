#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_one_country {

static constexpr std::array<sunindextype, 37> dxdotdw_colptrs_one_country_ = {
    0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72
};

void dxdotdw_colptrs_one_country(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_one_country_));
}
} // namespace model_one_country
} // namespace amici

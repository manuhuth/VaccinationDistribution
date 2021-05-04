#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_one_country {

static constexpr std::array<sunindextype, 61> dwdp_colptrs_one_country_ = {
    0, 24, 48, 96, 168, 172, 184, 188, 200, 204, 216, 220, 232, 268, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 307, 310, 313, 316
};

void dwdp_colptrs_one_country(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_one_country_));
}
} // namespace model_one_country
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_one_country {

static constexpr std::array<sunindextype, 22> dwdx_colptrs_one_country_ = {
    0, 20, 38, 56, 76, 96, 116, 136, 156, 176, 196, 216, 234, 252, 270, 288, 288, 288, 288, 288, 288, 288
};

void dwdx_colptrs_one_country(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_one_country_));
}
} // namespace model_one_country
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_one_country {

static constexpr std::array<sunindextype, 38> dwdp_colptrs_one_country_ = {
    0, 12, 24, 36, 54, 57, 60, 62, 65, 67, 70, 72, 75, 77, 80, 89, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98
};

void dwdp_colptrs_one_country(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_one_country_));
}
} // namespace model_one_country
} // namespace amici

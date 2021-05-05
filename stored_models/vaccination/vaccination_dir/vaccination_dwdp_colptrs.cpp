#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<sunindextype, 38> dwdp_colptrs_vaccination_ = {
    0, 12, 24, 36, 54, 56, 59, 61, 64, 66, 69, 71, 74, 83, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 95, 98
};

void dwdp_colptrs_vaccination(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_vaccination_));
}
} // namespace model_vaccination
} // namespace amici

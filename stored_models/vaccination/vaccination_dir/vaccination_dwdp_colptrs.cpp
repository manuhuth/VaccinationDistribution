#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<sunindextype, 67> dwdp_colptrs_vaccination_ = {
    0, 24, 48, 96, 168, 172, 184, 188, 200, 204, 216, 220, 232, 268, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 304, 306, 308, 310, 312, 314, 316, 318, 320, 322, 324
};

void dwdp_colptrs_vaccination(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_vaccination_));
}
} // namespace model_vaccination
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<sunindextype, 43> dwdx_colptrs_vaccination_ = {
    0, 76, 148, 220, 296, 368, 440, 518, 596, 674, 752, 830, 908, 986, 1064, 1142, 1220, 1298, 1376, 1452, 1528, 1600, 1672, 1744, 1816, 1892, 1968, 2040, 2112, 2184, 2256, 2256, 2256, 2256, 2256, 2256, 2256, 2256, 2256, 2256, 2256, 2256, 2256
};

void dwdx_colptrs_vaccination(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_vaccination_));
}
} // namespace model_vaccination
} // namespace amici

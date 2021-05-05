#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<sunindextype, 22> dwdx_colptrs_vaccination_ = {
    0, 20, 38, 56, 76, 96, 116, 136, 156, 176, 196, 216, 234, 252, 270, 288, 288, 288, 288, 288, 288, 288
};

void dwdx_colptrs_vaccination(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_vaccination_));
}
} // namespace model_vaccination
} // namespace amici

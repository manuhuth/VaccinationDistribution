#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<sunindextype, 22> dwdx_colptrs_vaccination_ = {
    0, 22, 40, 58, 78, 98, 118, 138, 158, 178, 200, 222, 240, 258, 276, 294, 294, 294, 294, 294, 294, 294
};

void dwdx_colptrs_vaccination(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_vaccination_));
}
} // namespace model_vaccination
} // namespace amici

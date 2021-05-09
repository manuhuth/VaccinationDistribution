#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<sunindextype, 43> dwdx_colptrs_vaccination_ = {
    0, 80, 156, 232, 312, 388, 464, 542, 620, 698, 776, 854, 932, 1010, 1088, 1166, 1244, 1322, 1400, 1480, 1560, 1636, 1712, 1788, 1864, 1944, 2024, 2100, 2176, 2252, 2328, 2328, 2328, 2328, 2328, 2328, 2328, 2328, 2328, 2328, 2328, 2328, 2328
};

void dwdx_colptrs_vaccination(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_vaccination_));
}
} // namespace model_vaccination
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination_partly {

static constexpr std::array<sunindextype, 44> dwdx_colptrs_vaccination_partly_ = {
    0, 76, 148, 220, 296, 368, 440, 514, 588, 662, 736, 810, 884, 958, 1032, 1106, 1180, 1254, 1328, 1404, 1480, 1552, 1624, 1696, 1768, 1844, 1920, 1992, 2064, 2136, 2208, 2208, 2208, 2208, 2208, 2208, 2208, 2208, 2208, 2208, 2208, 2208, 2208, 2210
};

void dwdx_colptrs_vaccination_partly(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_vaccination_partly_));
}
} // namespace model_vaccination_partly
} // namespace amici

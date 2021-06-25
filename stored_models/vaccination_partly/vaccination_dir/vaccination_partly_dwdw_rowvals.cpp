#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination_partly {

static constexpr std::array<sunindextype, 20> dwdw_rowvals_vaccination_partly_ = {
    3, 2, 4, 6, 5, 7, 8, 9, 107, 111, 113, 106, 110, 112, 109, 115, 117, 108, 114, 116
};

void dwdw_rowvals_vaccination_partly(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_vaccination_partly_));
}
} // namespace model_vaccination_partly
} // namespace amici

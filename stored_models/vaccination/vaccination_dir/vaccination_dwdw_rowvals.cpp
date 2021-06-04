#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<sunindextype, 20> dwdw_rowvals_vaccination_ = {
    3, 2, 4, 7, 5, 6, 8, 107, 111, 113, 9, 106, 110, 112, 108, 114, 116, 109, 115, 117
};

void dwdw_rowvals_vaccination(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_vaccination_));
}
} // namespace model_vaccination
} // namespace amici

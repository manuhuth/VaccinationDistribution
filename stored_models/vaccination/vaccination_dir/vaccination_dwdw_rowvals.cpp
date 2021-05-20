#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<sunindextype, 16> dwdw_rowvals_vaccination_ = {
    7, 5, 4, 6, 105, 109, 111, 107, 113, 115, 104, 108, 110, 106, 112, 114
};

void dwdw_rowvals_vaccination(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_vaccination_));
}
} // namespace model_vaccination
} // namespace amici

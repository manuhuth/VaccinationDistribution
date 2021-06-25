#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<sunindextype, 24> dwdw_rowvals_vaccination_ = {
    4, 5, 8, 10, 7, 11, 6, 7, 8, 9, 11, 109, 113, 115, 108, 112, 114, 10, 110, 116, 118, 111, 117, 119
};

void dwdw_rowvals_vaccination(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_vaccination_));
}
} // namespace model_vaccination
} // namespace amici

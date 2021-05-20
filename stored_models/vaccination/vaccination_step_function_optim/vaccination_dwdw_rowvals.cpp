#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<sunindextype, 12> dwdw_rowvals_vaccination_ = {
    101, 105, 107, 100, 104, 106, 102, 108, 110, 103, 109, 111
};

void dwdw_rowvals_vaccination(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_vaccination_));
}
} // namespace model_vaccination
} // namespace amici

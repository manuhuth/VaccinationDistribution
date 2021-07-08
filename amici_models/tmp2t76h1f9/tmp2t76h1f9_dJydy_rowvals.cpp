#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_tmp2t76h1f9 {

static constexpr std::array<std::array<sunindextype, 1>, 1> dJydy_rowvals_tmp2t76h1f9_ = {{
    {0}, 
}};

void dJydy_rowvals_tmp2t76h1f9(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_tmp2t76h1f9_[index]));
}
} // namespace model_tmp2t76h1f9
} // namespace amici

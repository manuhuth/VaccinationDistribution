#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_basic_sir {

static constexpr std::array<std::array<sunindextype, 5>, 4> dJydy_colptrs_basic_sir_ = {{
    {0, 1, 1, 1, 1}, 
    {0, 0, 1, 1, 1}, 
    {0, 0, 0, 1, 1}, 
    {0, 0, 0, 0, 1}, 
}};

void dJydy_colptrs_basic_sir(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_basic_sir_[index]));
}
} // namespace model_basic_sir
} // namespace amici

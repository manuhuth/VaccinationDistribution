#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<std::array<sunindextype, 9>, 8> dJydy_colptrs_vaccination_ = {{
    {0, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 1}, 
}};

void dJydy_colptrs_vaccination(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_vaccination_[index]));
}
} // namespace model_vaccination
} // namespace amici

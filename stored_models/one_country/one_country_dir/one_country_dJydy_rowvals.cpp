#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_one_country {

static constexpr std::array<std::array<sunindextype, 1>, 22> dJydy_rowvals_one_country_ = {{
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_one_country(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_one_country_[index]));
}
} // namespace model_one_country
} // namespace amici

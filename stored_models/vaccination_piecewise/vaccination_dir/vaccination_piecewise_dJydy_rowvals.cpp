#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination_piecewise {

static constexpr std::array<std::array<int, 1>, 12> dJydy_rowvals_vaccination_piecewise_ = {{
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

void dJydy_rowvals_vaccination_piecewise(sunindextype *rowvals, int index){
    std::copy(dJydy_rowvals_vaccination_piecewise_[index].begin(), dJydy_rowvals_vaccination_piecewise_[index].end(), rowvals);
}
} // namespace amici
} // namespace model_vaccination_piecewise
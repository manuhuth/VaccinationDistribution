#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<std::array<int, 1>, 4> dJydy_rowvals_vaccination_ = {{
    {0}, 
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_vaccination(sunindextype *rowvals, int index){
    std::copy(dJydy_rowvals_vaccination_[index].begin(), dJydy_rowvals_vaccination_[index].end(), rowvals);
}
} // namespace amici
} // namespace model_vaccination
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<std::array<int, 5>, 4> dJydy_colptrs_vaccination_ = {{
    {0, 1, 1, 1, 1}, 
    {0, 0, 1, 1, 1}, 
    {0, 0, 0, 1, 1}, 
    {0, 0, 0, 0, 1}, 
}};

void dJydy_colptrs_vaccination(sunindextype *colptrs, int index){
    std::copy(dJydy_colptrs_vaccination_[index].begin(), dJydy_colptrs_vaccination_[index].end(), colptrs);
}
} // namespace amici
} // namespace model_vaccination
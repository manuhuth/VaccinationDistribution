#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_tmp2kc_o5jm {

static constexpr std::array<std::array<int, 1>, 1> dJydy_rowvals_tmp2kc_o5jm_ = {{
    {0}, 
}};

void dJydy_rowvals_tmp2kc_o5jm(sunindextype *rowvals, int index){
    std::copy(dJydy_rowvals_tmp2kc_o5jm_[index].begin(), dJydy_rowvals_tmp2kc_o5jm_[index].end(), rowvals);
}
} // namespace amici
} // namespace model_tmp2kc_o5jm
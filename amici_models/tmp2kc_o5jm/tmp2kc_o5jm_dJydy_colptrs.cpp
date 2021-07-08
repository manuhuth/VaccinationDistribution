#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_tmp2kc_o5jm {

static constexpr std::array<std::array<int, 2>, 1> dJydy_colptrs_tmp2kc_o5jm_ = {{
    {0, 1}, 
}};

void dJydy_colptrs_tmp2kc_o5jm(sunindextype *colptrs, int index){
    std::copy(dJydy_colptrs_tmp2kc_o5jm_[index].begin(), dJydy_colptrs_tmp2kc_o5jm_[index].end(), colptrs);
}
} // namespace amici
} // namespace model_tmp2kc_o5jm
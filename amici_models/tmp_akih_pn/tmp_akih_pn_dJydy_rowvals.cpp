#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_tmp_akih_pn {

static constexpr std::array<std::array<int, 1>, 1> dJydy_rowvals_tmp_akih_pn_ = {{
    {0}, 
}};

void dJydy_rowvals_tmp_akih_pn(sunindextype *rowvals, int index){
    std::copy(dJydy_rowvals_tmp_akih_pn_[index].begin(), dJydy_rowvals_tmp_akih_pn_[index].end(), rowvals);
}
} // namespace amici
} // namespace model_tmp_akih_pn
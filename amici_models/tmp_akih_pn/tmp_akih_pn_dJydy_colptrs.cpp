#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_tmp_akih_pn {

static constexpr std::array<std::array<int, 2>, 1> dJydy_colptrs_tmp_akih_pn_ = {{
    {0, 1}, 
}};

void dJydy_colptrs_tmp_akih_pn(sunindextype *colptrs, int index){
    std::copy(dJydy_colptrs_tmp_akih_pn_[index].begin(), dJydy_colptrs_tmp_akih_pn_[index].end(), colptrs);
}
} // namespace amici
} // namespace model_tmp_akih_pn
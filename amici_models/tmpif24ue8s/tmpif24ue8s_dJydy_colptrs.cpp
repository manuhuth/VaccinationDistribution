#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_tmpif24ue8s {

static constexpr std::array<std::array<int, 2>, 1> dJydy_colptrs_tmpif24ue8s_ = {{
    {0, 1}, 
}};

void dJydy_colptrs_tmpif24ue8s(sunindextype *colptrs, int index){
    std::copy(dJydy_colptrs_tmpif24ue8s_[index].begin(), dJydy_colptrs_tmpif24ue8s_[index].end(), colptrs);
}
} // namespace amici
} // namespace model_tmpif24ue8s
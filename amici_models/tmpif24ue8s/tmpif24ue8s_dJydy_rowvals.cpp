#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_tmpif24ue8s {

static constexpr std::array<std::array<int, 1>, 1> dJydy_rowvals_tmpif24ue8s_ = {{
    {0}, 
}};

void dJydy_rowvals_tmpif24ue8s(sunindextype *rowvals, int index){
    std::copy(dJydy_rowvals_tmpif24ue8s_[index].begin(), dJydy_rowvals_tmpif24ue8s_[index].end(), rowvals);
}
} // namespace amici
} // namespace model_tmpif24ue8s
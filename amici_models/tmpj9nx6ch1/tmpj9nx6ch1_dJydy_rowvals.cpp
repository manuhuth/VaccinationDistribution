#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_tmpj9nx6ch1 {

static constexpr std::array<std::array<int, 1>, 1> dJydy_rowvals_tmpj9nx6ch1_ = {{
    {0}, 
}};

void dJydy_rowvals_tmpj9nx6ch1(sunindextype *rowvals, int index){
    std::copy(dJydy_rowvals_tmpj9nx6ch1_[index].begin(), dJydy_rowvals_tmpj9nx6ch1_[index].end(), rowvals);
}
} // namespace amici
} // namespace model_tmpj9nx6ch1
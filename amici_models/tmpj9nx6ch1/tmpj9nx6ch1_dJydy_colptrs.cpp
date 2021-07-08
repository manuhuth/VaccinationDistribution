#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_tmpj9nx6ch1 {

static constexpr std::array<std::array<int, 2>, 1> dJydy_colptrs_tmpj9nx6ch1_ = {{
    {0, 1}, 
}};

void dJydy_colptrs_tmpj9nx6ch1(sunindextype *colptrs, int index){
    std::copy(dJydy_colptrs_tmpj9nx6ch1_[index].begin(), dJydy_colptrs_tmpj9nx6ch1_[index].end(), colptrs);
}
} // namespace amici
} // namespace model_tmpj9nx6ch1
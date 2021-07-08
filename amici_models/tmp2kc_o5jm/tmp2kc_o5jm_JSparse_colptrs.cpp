#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_tmp2kc_o5jm {

static constexpr std::array<int, 43> JSparse_colptrs_tmp2kc_o5jm_ = {
    0, 24, 42, 60, 84, 102, 120, 145, 170, 190, 210, 230, 250, 275, 300, 320, 340, 360, 380, 404, 428, 446, 464, 482, 500, 524, 548, 566, 584, 602, 620, 620, 620, 620, 620, 620, 620, 620, 620, 620, 620, 620, 620
};

void JSparse_colptrs_tmp2kc_o5jm(sunindextype *colptrs){
    std::copy(JSparse_colptrs_tmp2kc_o5jm_.begin(), JSparse_colptrs_tmp2kc_o5jm_.end(), colptrs);
}
} // namespace amici
} // namespace model_tmp2kc_o5jm
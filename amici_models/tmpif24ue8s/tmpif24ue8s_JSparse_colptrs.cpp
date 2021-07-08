#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_tmpif24ue8s {

static constexpr std::array<int, 43> JSparse_colptrs_tmpif24ue8s_ = {
    0, 24, 42, 60, 84, 102, 120, 145, 170, 190, 210, 230, 250, 275, 300, 320, 340, 360, 380, 404, 428, 446, 464, 482, 500, 524, 548, 566, 584, 602, 620, 620, 620, 620, 620, 620, 620, 620, 620, 620, 620, 620, 620
};

void JSparse_colptrs_tmpif24ue8s(sunindextype *colptrs){
    std::copy(JSparse_colptrs_tmpif24ue8s_.begin(), JSparse_colptrs_tmpif24ue8s_.end(), colptrs);
}
} // namespace amici
} // namespace model_tmpif24ue8s
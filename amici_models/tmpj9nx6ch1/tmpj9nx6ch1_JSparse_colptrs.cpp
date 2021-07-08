#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_tmpj9nx6ch1 {

static constexpr std::array<int, 43> JSparse_colptrs_tmpj9nx6ch1_ = {
    0, 24, 42, 60, 84, 102, 120, 145, 170, 190, 210, 230, 250, 275, 300, 320, 340, 360, 380, 404, 428, 446, 464, 482, 500, 524, 548, 566, 584, 602, 620, 620, 620, 620, 620, 620, 620, 620, 620, 620, 620, 620, 620
};

void JSparse_colptrs_tmpj9nx6ch1(sunindextype *colptrs){
    std::copy(JSparse_colptrs_tmpj9nx6ch1_.begin(), JSparse_colptrs_tmpj9nx6ch1_.end(), colptrs);
}
} // namespace amici
} // namespace model_tmpj9nx6ch1
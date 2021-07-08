#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_tmpj9nx6ch1 {

static constexpr std::array<int, 43> JSparseB_colptrs_tmpj9nx6ch1_ = {
    0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450, 480, 510, 540, 545, 550, 556, 562, 568, 574, 579, 584, 590, 596, 602, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620
};

void JSparseB_colptrs_tmpj9nx6ch1(sunindextype *colptrs){
    std::copy(JSparseB_colptrs_tmpj9nx6ch1_.begin(), JSparseB_colptrs_tmpj9nx6ch1_.end(), colptrs);
}
} // namespace amici
} // namespace model_tmpj9nx6ch1
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination_neural_network {

static constexpr std::array<int, 43> JSparseB_colptrs_vaccination_neural_network_ = {
    0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450, 480, 510, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564
};

void JSparseB_colptrs_vaccination_neural_network(sunindextype *colptrs){
    std::copy(JSparseB_colptrs_vaccination_neural_network_.begin(), JSparseB_colptrs_vaccination_neural_network_.end(), colptrs);
}
} // namespace amici
} // namespace model_vaccination_neural_network
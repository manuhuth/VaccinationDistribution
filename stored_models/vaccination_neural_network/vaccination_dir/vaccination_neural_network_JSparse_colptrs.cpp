#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination_neural_network {

static constexpr std::array<int, 43> JSparse_colptrs_vaccination_neural_network_ = {
    0, 18, 36, 54, 72, 90, 108, 128, 148, 168, 188, 208, 228, 248, 268, 288, 308, 328, 348, 366, 384, 402, 420, 438, 456, 474, 492, 510, 528, 546, 564, 564, 564, 564, 564, 564, 564, 564, 564, 564, 564, 564, 564
};

void JSparse_colptrs_vaccination_neural_network(sunindextype *colptrs){
    std::copy(JSparse_colptrs_vaccination_neural_network_.begin(), JSparse_colptrs_vaccination_neural_network_.end(), colptrs);
}
} // namespace amici
} // namespace model_vaccination_neural_network
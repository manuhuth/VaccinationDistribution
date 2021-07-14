#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination_neural_network {

static constexpr std::array<std::array<int, 10>, 9> dJydy_colptrs_vaccination_neural_network_ = {{
    {0, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, 
}};

void dJydy_colptrs_vaccination_neural_network(sunindextype *colptrs, int index){
    std::copy(dJydy_colptrs_vaccination_neural_network_[index].begin(), dJydy_colptrs_vaccination_neural_network_[index].end(), colptrs);
}
} // namespace amici
} // namespace model_vaccination_neural_network
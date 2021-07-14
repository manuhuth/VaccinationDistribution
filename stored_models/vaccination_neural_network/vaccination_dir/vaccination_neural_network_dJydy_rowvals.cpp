#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination_neural_network {

static constexpr std::array<std::array<int, 1>, 9> dJydy_rowvals_vaccination_neural_network_ = {{
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_vaccination_neural_network(sunindextype *rowvals, int index){
    std::copy(dJydy_rowvals_vaccination_neural_network_[index].begin(), dJydy_rowvals_vaccination_neural_network_[index].end(), rowvals);
}
} // namespace amici
} // namespace model_vaccination_neural_network
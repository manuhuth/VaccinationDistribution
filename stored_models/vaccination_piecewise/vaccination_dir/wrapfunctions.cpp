#include "amici/model.h"
#include "wrapfunctions.h"
#include "vaccination_piecewise.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_vaccination_piecewise::Model_vaccination_piecewise());
}


} // namespace generic_model

} // namespace amici

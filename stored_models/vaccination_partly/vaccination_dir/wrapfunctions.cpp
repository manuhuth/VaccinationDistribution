#include "amici/model.h"
#include "wrapfunctions.h"
#include "vaccination_partly.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_vaccination_partly::Model_vaccination_partly());
}


} // namespace generic_model

} // namespace amici

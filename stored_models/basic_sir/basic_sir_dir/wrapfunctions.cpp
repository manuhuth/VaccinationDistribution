#include "amici/model.h"
#include "wrapfunctions.h"
#include "basic_sir.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_basic_sir::Model_basic_sir());
}


} // namespace generic_model

} // namespace amici

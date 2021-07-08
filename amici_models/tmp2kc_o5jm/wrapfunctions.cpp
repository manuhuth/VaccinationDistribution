#include "amici/model.h"
#include "wrapfunctions.h"
#include "tmp2kc_o5jm.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_tmp2kc_o5jm::Model_tmp2kc_o5jm());
}


} // namespace generic_model

} // namespace amici

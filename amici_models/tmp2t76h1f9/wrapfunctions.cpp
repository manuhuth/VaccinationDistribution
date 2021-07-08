#include "amici/model.h"
#include "wrapfunctions.h"
#include "tmp2t76h1f9.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_tmp2t76h1f9::Model_tmp2t76h1f9());
}


} // namespace generic_model

} // namespace amici

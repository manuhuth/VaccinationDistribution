#include "amici/model.h"
#include "wrapfunctions.h"
#include "tmp0n6wklye.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_tmp0n6wklye::Model_tmp0n6wklye());
}


} // namespace generic_model

} // namespace amici

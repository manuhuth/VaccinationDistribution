#include "amici/model.h"
#include "wrapfunctions.h"
#include "tmpif24ue8s.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_tmpif24ue8s::Model_tmpif24ue8s());
}


} // namespace generic_model

} // namespace amici

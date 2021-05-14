#include "amici/model.h"
#include "wrapfunctions.h"
#include "vaccination.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_vaccination::Model_vaccination());
}


} // namespace generic_model

} // namespace amici

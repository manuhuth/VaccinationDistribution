#include "amici/model.h"
#include "wrapfunctions.h"
#include "one_country.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_one_country::Model_one_country());
}


} // namespace generic_model

} // namespace amici

#include "amici/model.h"
#include "wrapfunctions.h"
#include "test.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_test::Model_test());
}


} // namespace generic_model

} // namespace amici

#include "amici/model.h"
#include "wrapfunctions.h"
#include "tmp_akih_pn.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_tmp_akih_pn::Model_tmp_akih_pn());
}


} // namespace generic_model

} // namespace amici

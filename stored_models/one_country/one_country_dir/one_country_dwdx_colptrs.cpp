#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_one_country {

static constexpr std::array<sunindextype, 43> dwdx_colptrs_one_country_ = {
    0, 74, 146, 218, 292, 364, 436, 510, 584, 658, 732, 806, 880, 954, 1028, 1102, 1176, 1250, 1324, 1398, 1472, 1544, 1616, 1688, 1760, 1834, 1908, 1980, 2052, 2124, 2196, 2196, 2196, 2196, 2196, 2196, 2196, 2196, 2196, 2196, 2196, 2196, 2196
};

void dwdx_colptrs_one_country(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_one_country_));
}
} // namespace model_one_country
} // namespace amici

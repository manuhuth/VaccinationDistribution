#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<sunindextype, 312> dwdp_rowvals_vaccination_ = {
    4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 8, 9, 20, 21, 29, 32, 41, 44, 53, 56, 65, 68, 77, 80, 89, 92, 12, 13, 24, 25, 30, 33, 42, 45, 54, 57, 66, 69, 78, 81, 90, 93, 10, 11, 22, 23, 35, 38, 47, 50, 59, 62, 71, 74, 83, 86, 95, 98, 14, 15, 26, 27, 36, 39, 48, 51, 60, 63, 72, 75, 84, 87, 96, 99, 28, 29, 30, 31, 32, 33, 40, 41, 42, 43, 44, 45, 52, 53, 54, 55, 56, 57, 64, 65, 66, 67, 68, 69, 76, 77, 78, 79, 80, 81, 88, 89, 90, 91, 92, 93, 34, 35, 36, 37, 38, 39, 46, 47, 48, 49, 50, 51, 58, 59, 60, 61, 62, 63, 70, 71, 72, 73, 74, 75, 82, 83, 84, 85, 86, 87, 94, 95, 96, 97, 98, 99, 0, 3, 1, 2, 0, 3, 2, 1
};

void dwdp_rowvals_vaccination(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_vaccination_));
}
} // namespace model_vaccination
} // namespace amici

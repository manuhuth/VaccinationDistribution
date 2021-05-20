#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_vaccination {

static constexpr std::array<sunindextype, 324> dwdp_rowvals_vaccination_ = {
    8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 12, 13, 24, 25, 33, 36, 45, 48, 57, 60, 69, 72, 81, 84, 93, 96, 16, 17, 28, 29, 34, 37, 46, 49, 58, 61, 70, 73, 82, 85, 94, 97, 14, 15, 26, 27, 39, 42, 51, 54, 63, 66, 75, 78, 87, 90, 99, 102, 18, 19, 30, 31, 40, 43, 52, 55, 64, 67, 76, 79, 88, 91, 100, 103, 32, 33, 34, 35, 36, 37, 44, 45, 46, 47, 48, 49, 56, 57, 58, 59, 60, 61, 68, 69, 70, 71, 72, 73, 80, 81, 82, 83, 84, 85, 92, 93, 94, 95, 96, 97, 38, 39, 40, 41, 42, 43, 50, 51, 52, 53, 54, 55, 62, 63, 64, 65, 66, 67, 74, 75, 76, 77, 78, 79, 86, 87, 88, 89, 90, 91, 98, 99, 100, 101, 102, 103, 6, 7, 4, 5, 0, 3, 0, 3, 0, 3, 0, 3, 1, 2, 1, 2, 1, 2, 1, 2
};

void dwdp_rowvals_vaccination(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_vaccination_));
}
} // namespace model_vaccination
} // namespace amici

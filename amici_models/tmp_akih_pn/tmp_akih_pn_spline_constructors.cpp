#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include "amici/splinefunctions.h"
#include <vector>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_tmp_akih_pn {

std::vector<HermiteSpline> spline_constructors_tmp_akih_pn(const realtype *p, const realtype *k){
	std::vector<HermiteSpline> splines;

	std::vector<realtype> nodes0 {xx0, xx1, xx2, xx3, xx4, xx5, xx6, xx7, xx8, xx9, xx10};
	std::vector<realtype> values0 {yy_countryA_vac1_0, yy_countryA_vac1_1, yy_countryA_vac1_2, yy_countryA_vac1_3, yy_countryA_vac1_4, yy_countryA_vac1_5, yy_countryA_vac1_6, yy_countryA_vac1_7, yy_countryA_vac1_8, yy_countryA_vac1_9, yy_countryA_vac1_10};
	std::vector<realtype> slopes0;
	HermiteSpline spline0 = HermiteSpline(nodes0, values0, slopes0, SplineBoundaryCondition::given, SplineBoundaryCondition::given, SplineExtrapolation::noExtrapolation, SplineExtrapolation::noExtrapolation, true, false, false);
	splines.push_back(spline0);

	std::vector<realtype> nodes1 {xx0, xx1, xx2, xx3, xx4, xx5, xx6, xx7, xx8, xx9, xx10};
	std::vector<realtype> values1 {yy_countryA_vac2_0, yy_countryA_vac2_1, yy_countryA_vac2_2, yy_countryA_vac2_3, yy_countryA_vac2_4, yy_countryA_vac2_5, yy_countryA_vac2_6, yy_countryA_vac2_7, yy_countryA_vac2_8, yy_countryA_vac2_9, yy_countryA_vac2_10};
	std::vector<realtype> slopes1;
	HermiteSpline spline1 = HermiteSpline(nodes1, values1, slopes1, SplineBoundaryCondition::given, SplineBoundaryCondition::given, SplineExtrapolation::noExtrapolation, SplineExtrapolation::noExtrapolation, true, false, false);
	splines.push_back(spline1);

return splines;
}

} // namespace amici
} // namespace model_tmp_akih_pn
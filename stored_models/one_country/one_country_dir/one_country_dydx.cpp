#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "one_country_x.h"
#include "one_country_p.h"
#include "one_country_w.h"
#include "one_country_dwdx.h"

namespace amici {
namespace model_one_country {

void dydx_one_country(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[45] = 1;
    dydx[90] = 1;
    dydx[135] = 1;
    dydx[180] = 1;
    dydx[225] = 1;
    dydx[270] = 1;
    dydx[315] = 1;
    dydx[360] = 1;
    dydx[405] = 1;
    dydx[450] = 1;
    dydx[495] = 1;
    dydx[540] = 1;
    dydx[585] = 1;
    dydx[630] = 1;
    dydx[675] = 1;
    dydx[720] = 1;
    dydx[765] = 1;
    dydx[810] = 1;
    dydx[855] = 1;
    dydx[900] = 1;
    dydx[945] = 1;
    dydx[990] = 1;
    dydx[1035] = 1;
    dydx[1080] = 1;
    dydx[1125] = 1;
    dydx[1170] = 1;
    dydx[1215] = 1;
    dydx[1260] = 1;
    dydx[1305] = 1;
    dydx[1350] = 1;
    dydx[1395] = 1;
    dydx[1440] = 1;
    dydx[1485] = 1;
    dydx[1530] = 1;
    dydx[1575] = 1;
    dydx[1620] = 1;
    dydx[1665] = 1;
    dydx[1710] = 1;
    dydx[1755] = 1;
    dydx[1800] = 1;
    dydx[1845] = 1;
}

} // namespace model_one_country
} // namespace amici

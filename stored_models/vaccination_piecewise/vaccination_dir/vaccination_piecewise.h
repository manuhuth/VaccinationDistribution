#ifndef _amici_TPL_MODELNAME_h
#define _amici_TPL_MODELNAME_h
#include <cmath>
#include <memory>

#include "amici/model_ode.h"
#include "amici/solver_cvodes.h"
#include "amici/splinefunctions.h"

#include "sundials/sundials_types.h"

namespace amici {

class Solver;

namespace model_vaccination_piecewise {

extern std::array<const char*, 112> parameterNames;
extern std::array<const char*, 0> fixedParameterNames;
extern std::array<const char*, 42> stateNames;
extern std::array<const char*, 12> observableNames;
extern std::array<const char*, 112> parameterIds;
extern std::array<const char*, 0> fixedParameterIds;
extern std::array<const char*, 42> stateIds;
extern std::array<const char*, 12> observableIds;


extern void J_vaccination_piecewise(realtype *J, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w,
                            const realtype *dwdx);
extern void JB_vaccination_piecewise(realtype *JB, const realtype t, const realtype *x,
                             const realtype *p, const realtype *k,
                             const realtype *h, const realtype *xB,
                             const realtype *w, const realtype *dwdx);
extern void JDiag_vaccination_piecewise(realtype *JDiag, const realtype t,
                                const realtype *x, const realtype *p,
                                const realtype *k, const realtype *h,
                                const realtype *w, const realtype *dwdx);
extern void JSparse_vaccination_piecewise(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);
extern void JSparse_colptrs_vaccination_piecewise(sunindextype *colptrs);
extern void JSparse_rowvals_vaccination_piecewise(sunindextype *rowvals);
extern void JSparseB_vaccination_piecewise(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx);
extern void JSparseB_colptrs_vaccination_piecewise(sunindextype *colptrs);
extern void JSparseB_rowvals_vaccination_piecewise(sunindextype *rowvals);
extern void Jy_vaccination_piecewise(realtype *nllh, const int iy, const realtype *p,
                             const realtype *k, const realtype *y,
                             const realtype *sigmay, const realtype *my);
extern void dJydsigmay_vaccination_piecewise(realtype *dJydsigmay, const int iy,
                                     const realtype *p, const realtype *k,
                                     const realtype *y, const realtype *sigmay,
                                     const realtype *my);
extern void dJydy_vaccination_piecewise(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_colptrs_vaccination_piecewise(sunindextype *colptrs, int index);
extern void dJydy_rowvals_vaccination_piecewise(sunindextype *rowvals, int index);
extern void dwdp_vaccination_piecewise(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp,const realtype *spl, const realtype *sspl);
extern void dwdp_colptrs_vaccination_piecewise(sunindextype *colptrs);
extern void dwdp_rowvals_vaccination_piecewise(sunindextype *rowvals);
extern void dwdx_vaccination_piecewise(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl);
extern void dwdx_colptrs_vaccination_piecewise(sunindextype *colptrs);
extern void dwdx_rowvals_vaccination_piecewise(sunindextype *rowvals);
extern void dxdotdw_vaccination_piecewise(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdw_colptrs_vaccination_piecewise(sunindextype *colptrs);
extern void dxdotdw_rowvals_vaccination_piecewise(sunindextype *rowvals);
extern void dxdotdp_explicit_vaccination_piecewise(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdp_explicit_colptrs_vaccination_piecewise(sunindextype *colptrs);
extern void dxdotdp_explicit_rowvals_vaccination_piecewise(sunindextype *rowvals);
extern void dxdotdp_implicit_colptrs_vaccination_piecewise(sunindextype *colptrs);
extern void dxdotdp_implicit_rowvals_vaccination_piecewise(sunindextype *rowvals);

extern void dydx_vaccination_piecewise(realtype *dydx, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w, const realtype *dwdx);
extern void dydp_vaccination_piecewise(realtype *dydp, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const int ip, const realtype *w,
                               const realtype *dwp);
extern void dsigmaydp_vaccination_piecewise(realtype *dsigmaydp, const realtype t,
                                    const realtype *p, const realtype *k,
                                    const int ip);
extern void sigmay_vaccination_piecewise(realtype *sigmay, const realtype t,
                                 const realtype *p, const realtype *k);
extern std::vector<HermiteSpline> spline_constructors_vaccination_piecewise(const realtype *p,
                                                                    const realtype *k);
extern void dspline_valuesdp_vaccination_piecewise(realtype *dspline_valuesdp,
                                           const realtype *p, const realtype *k);
extern void dspline_slopesdp_vaccination_piecewise(realtype *dspline_slopesdp,
                                           const realtype *p, const realtype *k);

extern void w_vaccination_piecewise(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl);
extern void x0_vaccination_piecewise(realtype *x0, const realtype t, const realtype *p,
                             const realtype *k);
extern void x0_fixedParameters_vaccination_piecewise(realtype *x0, const realtype t,
                                             const realtype *p,
                                             const realtype *k);
extern void sx0_vaccination_piecewise(realtype *sx0, const realtype t,
                              const realtype *x0, const realtype *p,
                              const realtype *k, const int ip);
extern void sx0_fixedParameters_vaccination_piecewise(realtype *sx0, const realtype t,
                                              const realtype *x0,
                                              const realtype *p,
                                              const realtype *k, const int ip);
extern void xdot_vaccination_piecewise(realtype *xdot, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w);
extern void y_vaccination_piecewise(realtype *y, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w);

extern void x_solver_vaccination_piecewise(realtype *x_solver, const realtype *x_rdata);
extern void total_cl_vaccination_piecewise(realtype *total_cl, const realtype *x_rdata);

/**
 * @brief AMICI-generated model subclass.
 */
class Model_vaccination_piecewise : public amici::Model_ODE {
  public:
    /**
     * @brief Default constructor.
     */
    Model_vaccination_piecewise()
        : amici::Model_ODE(
              42,                                // nx_rdata
              42,                            // nxtrue_rdata
              42,                               // nx_solver
              42,                           // nxtrue_solver
              0,                        // nx_solver_reinit
              12,                                      // ny
              12,                                  // nytrue
              0,                                      // nz
              0,                                  // nztrue
              0,                                  // nevent
              1,                              // nobjective
              116,                                      // nw
              2284,                                   // ndwdx
              776,                                   // ndwdp
              232,                                // ndxdotdw
              std::vector<int>{1,1,1,1,1,1,1,1,1,1,1,1},                                  // ndjydy
              620,                                     // nnz
              42,                                     // ubw
              42,                                     // lbw
              amici::SecondOrderMode::none,                                  // o2mode
              std::vector<realtype>{0.01, 0.1, 0.5, 2.0, 0.5, 0.5, 0.6, 0.6, 0.5, 0.5, 0.6, 0.6, 1.0, 1.3, 200000.0, 0.0, 0.0, 200000.0, 0.0, 0.0, 100.0, 100.0, 0.0, 0.0, 0.0, 0.0, 100.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 14.0, 28.0, 42.0, 56.0, 70.0, 84.0, 98.0, 112.0, 126.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.001, 0.001, 1.0},       // dynamic parameters
              std::vector<realtype>{}, // fixedParameters
              std::vector<int>{},                          // plist
              std::vector<realtype>(42, 0.0),   // idlist
              std::vector<int>{},                          // z2event
              true,                                        // pythonGenerated
              0,                       // ndxdotdp_explicit
              996                        // ndxdotdp_implicit
          ) {}

    /**
     * @brief Clone this model instance.
     * @return A deep copy of this instance.
     */
    virtual amici::Model *clone() const override {
        return new Model_vaccination_piecewise(*this);
    }

    /** model specific implementation for fJ
     * @param J Matrix to which the Jacobian will be written
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heaviside vector
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJ(realtype *J, const realtype t, const realtype *x,
                    const realtype *p, const realtype *k, const realtype *h,
                    const realtype *w, const realtype *dwdx) override {
        J_vaccination_piecewise(J, t, x, p, k, h, w, dwdx);
    }

    /** model specific implementation for fJB
     * @param JB Matrix to which the Jacobian will be written
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heaviside vector
     * @param xB Vector with the adjoint states
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJB(realtype *JB, const realtype t, const realtype *x,
                     const realtype *p, const realtype *k, const realtype *h,
                     const realtype *xB, const realtype *w,
                     const realtype *dwdx) override {
        JB_vaccination_piecewise(JB, t, x, p, k, h, xB, w, dwdx);
    }

    /** model specific implementation for fJDiag
     * @param JDiag Matrix to which the Jacobian will be written
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heaviside vector
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJDiag(realtype *JDiag, const realtype t, const realtype *x,
                        const realtype *p, const realtype *k, const realtype *h,
                        const realtype *w, const realtype *dwdx) override {
        JDiag_vaccination_piecewise(JDiag, t, x, p, k, h, w, dwdx);
    }

    virtual void fJSparse(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        JSparse_vaccination_piecewise( JSparse,  t,  x,  p,  k,  h,  w,  dwdx);
    }


    virtual void fJSparse_colptrs(sunindextype *colptrs) override {
        JSparse_colptrs_vaccination_piecewise(colptrs);
    }


    virtual void fJSparse_rowvals(sunindextype *rowvals) override {
        JSparse_rowvals_vaccination_piecewise(rowvals);
    }


    virtual void fJSparseB(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) override {
        JSparseB_vaccination_piecewise( JSparseB,  t,  x,  p,  k,  h,  xB,  w,  dwdx);
    }


    virtual void fJSparseB_colptrs(sunindextype *colptrs) override {
        JSparseB_colptrs_vaccination_piecewise(colptrs);
    }


    virtual void fJSparseB_rowvals(sunindextype *rowvals) override {
        JSparseB_rowvals_vaccination_piecewise(rowvals);
    }


    /** model specific implementation of fJrz
     * @param nllh regularization for event measurements z
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     **/
    virtual void fJrz(realtype *nllh, const int iz, const realtype *p,
                      const realtype *k, const realtype *rz,
                      const realtype *sigmaz) override {}

    /** model specific implementation of fJy
     * @param nllh negative log-likelihood for measurements y
     * @param iy output index
     * @param p parameter vector
     * @param k constant vector
     * @param y model output at timepoint
     * @param sigmay measurement standard deviation at timepoint
     * @param my measurements at timepoint
     **/
    virtual void fJy(realtype *nllh, const int iy, const realtype *p,
                     const realtype *k, const realtype *y,
                     const realtype *sigmay, const realtype *my) override {
        Jy_vaccination_piecewise(nllh, iy, p, k, y, sigmay, my);
    }

    /** model specific implementation of fJz
     * @param nllh negative log-likelihood for event measurements z
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     * @param mz event measurements at timepoint
     **/
    virtual void fJz(realtype *nllh, const int iz, const realtype *p,
                     const realtype *k, const realtype *z,
                     const realtype *sigmaz, const realtype *mz) override {}

    /** model specific implementation of fdJrzdsigma
     * @param dJrzdsigma Sensitivity of event penalization Jrz w.r.t.
     * standard deviation sigmaz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param rz model root output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     **/
    virtual void fdJrzdsigma(realtype *dJrzdsigma, const int iz,
                             const realtype *p, const realtype *k,
                             const realtype *rz,
                             const realtype *sigmaz) override {}

    /** model specific implementation of fdJrzdz
     * @param dJrzdz partial derivative of event penalization Jrz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param rz model root output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     **/
    virtual void fdJrzdz(realtype *dJrzdz, const int iz, const realtype *p,
                         const realtype *k, const realtype *rz,
                         const realtype *sigmaz) override {}

    /** model specific implementation of fdJydsigma
     * @param dJydsigma Sensitivity of time-resolved measurement
     * negative log-likelihood Jy w.r.t. standard deviation sigmay
     * @param iy output index
     * @param p parameter vector
     * @param k constant vector
     * @param y model output at timepoint
     * @param sigmay measurement standard deviation at timepoint
     * @param my measurement at timepoint
     **/
    virtual void fdJydsigma(realtype *dJydsigma, const int iy,
                            const realtype *p, const realtype *k,
                            const realtype *y, const realtype *sigmay,
                            const realtype *my) override {
        dJydsigmay_vaccination_piecewise(dJydsigma, iy, p, k, y, sigmay, my);
    }


    /** model specific implementation of fdJzdsigma
     * @param dJzdsigma Sensitivity of event measurement
     * negative log-likelihood Jz w.r.t. standard deviation sigmaz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     * @param mz event measurement at timepoint
     **/
    virtual void fdJzdsigma(realtype *dJzdsigma, const int iz,
                            const realtype *p, const realtype *k,
                            const realtype *z, const realtype *sigmaz,
                            const realtype *mz) override {}

    /** model specific implementation of fdJzdz
     * @param dJzdz partial derivative of event measurement negative
     *log-likelihood Jz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     * @param mz event measurement at timepoint
     **/
    virtual void fdJzdz(realtype *dJzdz, const int iz, const realtype *p,
                        const realtype *k, const realtype *z,
                        const realtype *sigmaz, const realtype *mz) override {}

    /** model specific implementation of fdeltasx
     * @param deltaqB sensitivity update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     * @param ip sensitivity index
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     * @param xB adjoint state
     **/
    virtual void fdeltaqB(realtype *deltaqB, const realtype t,
                          const realtype *x, const realtype *p,
                          const realtype *k, const realtype *h, const int ip,
                          const int ie, const realtype *xdot,
                          const realtype *xdot_old,
                          const realtype *xB) override {}

    /** model specific implementation of fdeltasx
     * @param deltasx sensitivity update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     * @param w repeating elements vector
     * @param ip sensitivity index
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     * @param sx state sensitivity
     * @param stau event-time sensitivity
     **/
    virtual void fdeltasx(realtype *deltasx, const realtype t,
                          const realtype *x, const realtype *p,
                          const realtype *k, const realtype *h,
                          const realtype *w, const int ip, const int ie,
                          const realtype *xdot, const realtype *xdot_old,
                          const realtype *sx, const realtype *stau) override {}

    /** model specific implementation of fdeltax
     * @param deltax state update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     **/
    virtual void fdeltax(realtype *deltax, const realtype t, const realtype *x,
                         const realtype *p, const realtype *k,
                         const realtype *h, const int ie, const realtype *xdot,
                         const realtype *xdot_old) override {}

    /** model specific implementation of fdeltaxB
     * @param deltaxB adjoint state update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     * @param xB current adjoint state
     **/
    virtual void fdeltaxB(realtype *deltaxB, const realtype t,
                          const realtype *x, const realtype *p,
                          const realtype *k, const realtype *h, const int ie,
                          const realtype *xdot, const realtype *xdot_old,
                          const realtype *xB) override {}

    /** model specific implementation of fdrzdp
     * @param drzdp partial derivative of root output rz w.r.t. model parameters
     *p
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     * @param ip parameter index w.r.t. which the derivative is requested
     **/
    virtual void fdrzdp(realtype *drzdp, const int ie, const realtype t,
                        const realtype *x, const realtype *p, const realtype *k,
                        const realtype *h, const int ip) override {}

    /** model specific implementation of fdrzdx
     * @param drzdx partial derivative of root output rz w.r.t. model states x
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     **/
    virtual void fdrzdx(realtype *drzdx, const int ie, const realtype t,
                        const realtype *x, const realtype *p, const realtype *k,
                        const realtype *h) override {}

    /** model specific implementation of fsigmay
     * @param dsigmaydp partial derivative of standard deviation of measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     * @param ip sensitivity index
     **/
    virtual void fdsigmaydp(realtype *dsigmaydp, const realtype t,
                            const realtype *p, const realtype *k,
                            const int ip) override {
        dsigmaydp_vaccination_piecewise(dsigmaydp, t, p, k, ip);
    }

    /** model specific implementation of fsigmaz
     * @param dsigmazdp partial derivative of standard deviation of event
     *measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     * @param ip sensitivity index
     **/
    virtual void fdsigmazdp(realtype *dsigmazdp, const realtype t,
                            const realtype *p, const realtype *k,
                            const int ip) override {}

    virtual void fdJydy(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        dJydy_vaccination_piecewise( dJydy,  iy,  p,  k,  y,  sigmay,  my);
    }

    virtual void fdJydy_colptrs(sunindextype *colptrs, int index) override {
        dJydy_colptrs_vaccination_piecewise(colptrs, index);
    }

    virtual void fdJydy_rowvals(sunindextype *rowvals, int index) override {
        dJydy_rowvals_vaccination_piecewise(rowvals, index);
    }


    virtual std::vector<HermiteSpline> fspline_constructors(const realtype *p,
                                                            const realtype *k) override {
        return spline_constructors_vaccination_piecewise(p, k);
    }

    virtual void fdspline_valuesdp(realtype *dspline_valuesdp,
                                   const realtype *p,
                                   const realtype *k) override {
        dspline_valuesdp_vaccination_piecewise(dspline_valuesdp, p, k);
    }
    virtual void fdspline_slopesdp(realtype *dspline_slopesdp,
                                   const realtype *p,
                                   const realtype *k) override {
        dspline_valuesdp_vaccination_piecewise(dspline_slopesdp, p, k);
    }

    virtual void fdwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp,const realtype *spl, const realtype *sspl) override {
        dwdp_vaccination_piecewise( dwdp,  t,  x,  p,  k,  h,  w,  tcl,  dtcldp, spl,  sspl);
    }


    virtual void fdwdp_colptrs(sunindextype *colptrs) override {
        dwdp_colptrs_vaccination_piecewise(colptrs);
    }


    virtual void fdwdp_rowvals(sunindextype *rowvals) override {
        dwdp_rowvals_vaccination_piecewise(rowvals);
    }


    virtual void fdwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl) override {
        dwdx_vaccination_piecewise( dwdx,  t,  x,  p,  k,  h,  w,  tcl,  spl);
    }


    virtual void fdwdx_colptrs(sunindextype *colptrs) override {
        dwdx_colptrs_vaccination_piecewise(colptrs);
    }


    virtual void fdwdx_rowvals(sunindextype *rowvals) override {
        dwdx_rowvals_vaccination_piecewise(rowvals);
    }


    virtual void fdxdotdw(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdw_vaccination_piecewise( dxdotdw,  t,  x,  p,  k,  h,  w);
    }

    virtual void fdxdotdw_colptrs(sunindextype *colptrs) override {
        dxdotdw_colptrs_vaccination_piecewise(colptrs);
    }

    virtual void fdxdotdw_rowvals(sunindextype *rowvals) override {
        dxdotdw_rowvals_vaccination_piecewise(rowvals);
    }


    virtual void fdxdotdp_explicit(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdp_explicit_vaccination_piecewise( dxdotdp_explicit,  t,  x,  p,  k,  h,  w);
    }


    virtual void fdxdotdp_explicit_colptrs(sunindextype *colptrs) override {
        dxdotdp_explicit_colptrs_vaccination_piecewise(colptrs);
    }


    virtual void fdxdotdp_explicit_rowvals(sunindextype *rowvals) override {
        dxdotdp_explicit_rowvals_vaccination_piecewise(rowvals);
    }


    virtual void fdxdotdp_implicit_colptrs(sunindextype *colptrs) override {
        dxdotdp_implicit_colptrs_vaccination_piecewise(colptrs);
    }


    virtual void fdxdotdp_implicit_rowvals(sunindextype *rowvals) override {
        dxdotdp_implicit_rowvals_vaccination_piecewise(rowvals);
    }


    /** model specific implementation of fdydx
     * @param dydx partial derivative of observables y w.r.t. model states x
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     **/
    virtual void fdydx(realtype *dydx, const realtype t, const realtype *x,
                       const realtype *p, const realtype *k, const realtype *h,
                       const realtype *w, const realtype *dwdx) override {
        dydx_vaccination_piecewise(dydx, t, x, p, k, h, w, dwdx);
    }

    /** model specific implementation of fdydp
     * @param dydp partial derivative of observables y w.r.t. model parameters p
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     * @param ip parameter index w.r.t. which the derivative is requested
     **/
    virtual void fdydp(realtype *dydp, const realtype t, const realtype *x,
                       const realtype *p, const realtype *k, const realtype *h,
                       const int ip, const realtype *w,
                       const realtype *dwdp) override {
        dydp_vaccination_piecewise(dydp, t, x, p, k, h, ip, w, dwdp);
    }

    /** model specific implementation of fdzdp
     * @param dzdp partial derivative of event-resolved output z w.r.t. model
     *parameters p
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     * @param ip parameter index w.r.t. which the derivative is requested
     **/
    virtual void fdzdp(realtype *dzdp, const int ie, const realtype t,
                       const realtype *x, const realtype *p, const realtype *k,
                       const realtype *h, const int ip) override {}

    /** model specific implementation of fdzdx
     * @param dzdx partial derivative of event-resolved output z w.r.t. model
     *states x
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     **/
    virtual void fdzdx(realtype *dzdx, const int ie, const realtype t,
                       const realtype *x, const realtype *p, const realtype *k,
                       const realtype *h) override {}

    /** model specific implementation for froot
     * @param root values of the trigger function
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heaviside vector
     **/
    virtual void froot(realtype *root, const realtype t, const realtype *x,
                       const realtype *p, const realtype *k,
                       const realtype *h) override {}

    /** model specific implementation of frz
     * @param rz value of root function at current timepoint (non-output events
     *not included)
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     **/
    virtual void frz(realtype *rz, const int ie, const realtype t,
                     const realtype *x, const realtype *p, const realtype *k,
                     const realtype *h) override {}

    /** model specific implementation of fsigmay
     * @param sigmay standard deviation of measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fsigmay(realtype *sigmay, const realtype t, const realtype *p,
                         const realtype *k) override {
        sigmay_vaccination_piecewise(sigmay, t, p, k);
    }

    /** model specific implementation of fsigmaz
     * @param sigmaz standard deviation of event measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fsigmaz(realtype *sigmaz, const realtype t, const realtype *p,
                         const realtype *k) override {}

    /** model specific implementation of fsrz
     * @param srz Sensitivity of rz, total derivative
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param sx current state sensitivity
     * @param h heaviside vector
     * @param ip sensitivity index
     **/
    virtual void fsrz(realtype *srz, const int ie, const realtype t,
                      const realtype *x, const realtype *p, const realtype *k,
                      const realtype *h, const realtype *sx,
                      const int ip) override {}

    /** model specific implementation of fstau
     * @param stau total derivative of event timepoint
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     * @param sx current state sensitivity
     * @param ip sensitivity index
     * @param ie event index
     **/
    virtual void fstau(realtype *stau, const realtype t, const realtype *x,
                       const realtype *p, const realtype *k, const realtype *h,
                       const realtype *sx, const int ip,
                       const int ie) override {}

    /** model specific implementation of fsx0
     * @param sx0 initial state sensitivities
     * @param t initial time
     * @param x0 initial state
     * @param p parameter vector
     * @param k constant vector
     * @param ip sensitivity index
     **/
    virtual void fsx0(realtype *sx0, const realtype t, const realtype *x0,
                      const realtype *p, const realtype *k,
                      const int ip) override {
        sx0_vaccination_piecewise(sx0, t, x0, p, k, ip);
    }

    /** model specific implementation of fsx0_fixedParameters
     * @param sx0 initial state sensitivities
     * @param t initial time
     * @param x0 initial state
     * @param p parameter vector
     * @param k constant vector
     * @param ip sensitivity index
     **/
    virtual void fsx0_fixedParameters(realtype *sx0, const realtype t,
                                      const realtype *x0, const realtype *p,
                                      const realtype *k,
                                      const int ip) override {
        sx0_fixedParameters_vaccination_piecewise(sx0, t, x0, p, k, ip);
    }

    /** model specific implementation of fsz
     * @param sz Sensitivity of rz, total derivative
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     * @param sx current state sensitivity
     * @param ip sensitivity index
     **/
    virtual void fsz(realtype *sz, const int ie, const realtype t,
                     const realtype *x, const realtype *p, const realtype *k,
                     const realtype *h, const realtype *sx,
                     const int ip) override {}

    virtual void fw(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl) override {
        w_vaccination_piecewise( w,  t,  x,  p,  k,  h,  tcl,  spl);
    }


    /** model specific implementation of fx0
     * @param x0 initial state
     * @param t initial time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fx0(realtype *x0, const realtype t, const realtype *p,
                     const realtype *k) override {
        x0_vaccination_piecewise(x0, t, p, k);
    }

    /** model specific implementation of fx0_fixedParameters
     * @param x0 initial state
     * @param t initial time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fx0_fixedParameters(realtype *x0, const realtype t,
                                     const realtype *p,
                                     const realtype *k) override {
        x0_fixedParameters_vaccination_piecewise(x0, t, p, k);
    }

    /** model specific implementation for fxdot
     * @param xdot residual function
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heaviside vector
     * @param w vector with helper variables
     **/
    virtual void fxdot(realtype *xdot, const realtype t, const realtype *x,
                       const realtype *p, const realtype *k, const realtype *h,
                       const realtype *w) override {
        xdot_vaccination_piecewise(xdot, t, x, p, k, h, w);
    }

    /** model specific implementation of fy
     * @param y model output at current timepoint
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     **/
    virtual void fy(realtype *y, const realtype t, const realtype *x,
                    const realtype *p, const realtype *k, const realtype *h,
                    const realtype *w) override {
        y_vaccination_piecewise(y, t, x, p, k, h, w);
    }

    /** model specific implementation of fz
     * @param z value of event output
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     **/
    virtual void fz(realtype *z, const int ie, const realtype t,
                    const realtype *x, const realtype *p, const realtype *k,
                    const realtype *h) override {}

    

    virtual void fx_solver(realtype *x_solver, const realtype *x_rdata) override {
        x_solver_vaccination_piecewise( x_solver,  x_rdata);
    }


    virtual void ftotal_cl(realtype *total_cl, const realtype *x_rdata) override {
        total_cl_vaccination_piecewise( total_cl,  x_rdata);
    }



    std::string getName() const override {
        return "vaccination_piecewise";
    }

    /**
     * @brief Get names of the model parameters
     * @return the names
     */
    virtual std::vector<std::string> getParameterNames() const override {
        return std::vector<std::string>(parameterNames.begin(),
                                        parameterNames.end());
    }

    /**
     * @brief Get names of the model states
     * @return the names
     */
    virtual std::vector<std::string> getStateNames() const override {
        return std::vector<std::string>(stateNames.begin(), stateNames.end());
    }

    /**
     * @brief Get names of the fixed model parameters
     * @return the names
     */
    virtual std::vector<std::string> getFixedParameterNames() const override {
        return std::vector<std::string>(fixedParameterNames.begin(),
                                        fixedParameterNames.end());
    }

    /**
     * @brief Get names of the observables
     * @return the names
     */
    virtual std::vector<std::string> getObservableNames() const override {
        return std::vector<std::string>(observableNames.begin(),
                                        observableNames.end());
    }

    /**
     * @brief Get ids of the model parameters
     * @return the ids
     */
    virtual std::vector<std::string> getParameterIds() const override {
        return std::vector<std::string>(parameterIds.begin(),
                                        parameterIds.end());
    }

    /**
     * @brief Get ids of the model states
     * @return the ids
     */
    virtual std::vector<std::string> getStateIds() const override {
        return std::vector<std::string>(stateIds.begin(), stateIds.end());
    }

    /**
     * @brief Get ids of the fixed model parameters
     * @return the ids
     */
    virtual std::vector<std::string> getFixedParameterIds() const override {
        return std::vector<std::string>(fixedParameterIds.begin(),
                                        fixedParameterIds.end());
    }

    /**
     * @brief Get ids of the observables
     * @return the ids
     */
    virtual std::vector<std::string> getObservableIds() const override {
        return std::vector<std::string>(observableIds.begin(),
                                        observableIds.end());
    }

    /**
     * @brief function indicating whether reinitialization of states depending on
     fixed parameters is permissible
     * @return flag inidication whether reinitialization of states depending on
     fixed parameters is permissible
     */
    virtual bool isFixedParameterStateReinitializationAllowed() const override {
        return true;
    }

    /**
     * @brief returns the amici version that was used to generate the model
     * @return ver amici version string
     */
    virtual const std::string getAmiciVersion() const override {
        return "0.11.6";
    }

    /**
     & @brief returns the amici version that was used to generate the model
     * @return commit amici git commit hash
     */
    virtual const std::string getAmiciCommit() const override {
        return "0ab3af6849844a080d490c9026f074d2c8879dc4";
    }
};


} // namespace model_vaccination_piecewise

} // namespace amici

#endif /* _amici_TPL_MODELNAME_h */

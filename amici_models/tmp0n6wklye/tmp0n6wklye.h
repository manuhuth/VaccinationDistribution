#ifndef _amici_TPL_MODELNAME_h
#define _amici_TPL_MODELNAME_h
#include <cmath>
#include <memory>
#include <gsl/gsl-lite.hpp>

#include "amici/model_ode.h"
#include "amici/solver_cvodes.h"

#include "sundials/sundials_types.h"

namespace amici {

class Solver;

namespace model_tmp0n6wklye {

extern std::array<const char*, 22> parameterNames;
extern std::array<const char*, 93> fixedParameterNames;
extern std::array<const char*, 43> stateNames;
extern std::array<const char*, 1> observableNames;
extern std::array<const char*, 120> expressionNames;
extern std::array<const char*, 22> parameterIds;
extern std::array<const char*, 93> fixedParameterIds;
extern std::array<const char*, 43> stateIds;
extern std::array<const char*, 1> observableIds;
extern std::array<const char*, 120> expressionIds;

extern void Jy_tmp0n6wklye(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydsigma_tmp0n6wklye(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_tmp0n6wklye(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_colptrs_tmp0n6wklye(SUNMatrixWrapper &colptrs, int index);
extern void dJydy_rowvals_tmp0n6wklye(SUNMatrixWrapper &rowvals, int index);
extern void root_tmp0n6wklye(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
extern void dwdp_tmp0n6wklye(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp);
extern void dwdp_colptrs_tmp0n6wklye(SUNMatrixWrapper &colptrs);
extern void dwdp_rowvals_tmp0n6wklye(SUNMatrixWrapper &rowvals);
extern void dwdx_tmp0n6wklye(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl);
extern void dwdx_colptrs_tmp0n6wklye(SUNMatrixWrapper &colptrs);
extern void dwdx_rowvals_tmp0n6wklye(SUNMatrixWrapper &rowvals);
extern void dwdw_tmp0n6wklye(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl);
extern void dwdw_colptrs_tmp0n6wklye(SUNMatrixWrapper &colptrs);
extern void dwdw_rowvals_tmp0n6wklye(SUNMatrixWrapper &rowvals);
extern void dxdotdw_tmp0n6wklye(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdw_colptrs_tmp0n6wklye(SUNMatrixWrapper &colptrs);
extern void dxdotdw_rowvals_tmp0n6wklye(SUNMatrixWrapper &rowvals);






extern void dydx_tmp0n6wklye(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);

extern void sigmay_tmp0n6wklye(realtype *sigmay, const realtype t, const realtype *p, const realtype *k);

extern void w_tmp0n6wklye(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl);
extern void x0_tmp0n6wklye(realtype *x0, const realtype t, const realtype *p, const realtype *k);
extern void x0_fixedParameters_tmp0n6wklye(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs);

extern void sx0_fixedParameters_tmp0n6wklye(realtype *sx0_fixedParameters, const realtype t, const realtype *x0, const realtype *p, const realtype *k, const int ip, gsl::span<const int> reinitialization_state_idxs);
extern void xdot_tmp0n6wklye(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void y_tmp0n6wklye(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void stau_tmp0n6wklye(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip, const int ie);

extern void deltasx_tmp0n6wklye(realtype *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau);

extern void x_solver_tmp0n6wklye(realtype *x_solver, const realtype *x_rdata);


/**
 * @brief AMICI-generated model subclass.
 */
class Model_tmp0n6wklye : public amici::Model_ODE {
  public:
    /**
     * @brief Default constructor.
     */
    Model_tmp0n6wklye()
        : amici::Model_ODE(
              amici::ModelDimensions(
                  43,                            // nx_rdata
                  43,                        // nxtrue_rdata
                  43,                           // nx_solver
                  43,                       // nxtrue_solver
                  0,                    // nx_solver_reinit
                  22,                                  // np
                  93,                                  // nk
                  1,                                  // ny
                  1,                              // nytrue
                  0,                                  // nz
                  0,                              // nztrue
                  22,                              // nevent
                  1,                          // nobjective
                  120,                                  // nw
                  2212,                               // ndwdx
                  22,                               // ndwdp
                  24,                               // ndwdw
                  216,                            // ndxdotdw
                  std::vector<int>{1},                              // ndjydy
                  0,                                       // nnz
                  43,                                 // ubw
                  43                                  // lbw
              ),
              amici::SimulationParameters(
                  std::vector<realtype>{0.01, 0.10000000000000001, 0.5, 2.0, 0.5, 0.5, 0.59999999999999998, 0.59999999999999998, 0.5, 0.5, 0.59999999999999998, 0.59999999999999998, 1.0, 1.3, 200000.0, 0.0, 0.0, 200000.0, 0.0, 0.0, 100.0, 100.0, 0.0, 0.0, 0.0, 0.0, 100.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.001, 0.001, 1.0, 0.0, 14.0, 28.0, 42.0, 56.0, 70.0, 84.0, 98.0, 112.0, 126.0, 140.0}, // fixedParameters
                  std::vector<realtype>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}        // dynamic parameters
              ),
              amici::SecondOrderMode::none,                                  // o2mode
              std::vector<realtype>(43, 0.0),   // idlist
              std::vector<int>{},                          // z2event
              true,                                        // pythonGenerated
              0,                       // ndxdotdp_explicit
              0,                       // ndxdotdx_explicit
              3                        // w_recursion_depth
          ) {}

    /**
     * @brief Clone this model instance.
     * @return A deep copy of this instance.
     */
    virtual amici::Model *clone() const override {
        return new Model_tmp0n6wklye(*this);
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

    virtual void fJy(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        Jy_tmp0n6wklye(Jy, iy, p, k, y, sigmay, my);
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

    virtual void fdJydsigma(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        dJydsigma_tmp0n6wklye(dJydsigma, iy, p, k, y, sigmay, my);
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

    virtual void fdeltasx(realtype *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau) override {
        deltasx_tmp0n6wklye(deltasx, t, x, p, k, h, w, ip, ie, xdot, xdot_old, sx, stau);
    }


    virtual void fdeltax(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old) override {}


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

    virtual void fdsigmaydp(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const int ip) override {}


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
        dJydy_tmp0n6wklye(dJydy, iy, p, k, y, sigmay, my);
    }

    virtual void fdJydy_colptrs(SUNMatrixWrapper &colptrs, int index) override {        dJydy_colptrs_tmp0n6wklye(colptrs, index);
    }

    virtual void fdJydy_rowvals(SUNMatrixWrapper &rowvals, int index) override {        dJydy_rowvals_tmp0n6wklye(rowvals, index);
    }


    virtual void fdwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp) override {
        dwdp_tmp0n6wklye(dwdp, t, x, p, k, h, w, tcl, dtcldp);
    }

    virtual void fdwdp_colptrs(SUNMatrixWrapper &colptrs) override {        dwdp_colptrs_tmp0n6wklye(colptrs);
    }

    virtual void fdwdp_rowvals(SUNMatrixWrapper &rowvals) override {        dwdp_rowvals_tmp0n6wklye(rowvals);
    }


    virtual void fdwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl) override {
        dwdx_tmp0n6wklye(dwdx, t, x, p, k, h, w, tcl);
    }

    virtual void fdwdx_colptrs(SUNMatrixWrapper &colptrs) override {        dwdx_colptrs_tmp0n6wklye(colptrs);
    }

    virtual void fdwdx_rowvals(SUNMatrixWrapper &rowvals) override {        dwdx_rowvals_tmp0n6wklye(rowvals);
    }


    virtual void fdwdw(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl) override {
        dwdw_tmp0n6wklye(dwdw, t, x, p, k, h, w, tcl);
    }

    virtual void fdwdw_colptrs(SUNMatrixWrapper &colptrs) override {        dwdw_colptrs_tmp0n6wklye(colptrs);
    }

    virtual void fdwdw_rowvals(SUNMatrixWrapper &rowvals) override {        dwdw_rowvals_tmp0n6wklye(rowvals);
    }


    virtual void fdxdotdw(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdw_tmp0n6wklye(dxdotdw, t, x, p, k, h, w);
    }

    virtual void fdxdotdw_colptrs(SUNMatrixWrapper &colptrs) override {        dxdotdw_colptrs_tmp0n6wklye(colptrs);
    }

    virtual void fdxdotdw_rowvals(SUNMatrixWrapper &rowvals) override {        dxdotdw_rowvals_tmp0n6wklye(rowvals);
    }


    virtual void fdxdotdp_explicit(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {}

    virtual void fdxdotdp_explicit_colptrs(SUNMatrixWrapper &colptrs) override {}

    virtual void fdxdotdp_explicit_rowvals(SUNMatrixWrapper &rowvals) override {}


    virtual void fdxdotdx_explicit(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {}

    virtual void fdxdotdx_explicit_colptrs(SUNMatrixWrapper &colptrs) override {}

    virtual void fdxdotdx_explicit_rowvals(SUNMatrixWrapper &rowvals) override {}


    virtual void fdydx(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        dydx_tmp0n6wklye(dydx, t, x, p, k, h, w, dwdx);
    }


    virtual void fdydp(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dtcldp) override {}


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

    virtual void froot(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        root_tmp0n6wklye(root, t, x, p, k, h);
    }


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

    virtual void fsigmay(realtype *sigmay, const realtype t, const realtype *p, const realtype *k) override {
        sigmay_tmp0n6wklye(sigmay, t, p, k);
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

    virtual void fstau(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip, const int ie) override {
        stau_tmp0n6wklye(stau, t, x, p, k, h, sx, ip, ie);
    }

    virtual void fsx0(realtype *sx0, const realtype t,const realtype *x, const realtype *p, const realtype *k, const int ip) override {}

    virtual void fsx0_fixedParameters(realtype *sx0_fixedParameters, const realtype t, const realtype *x0, const realtype *p, const realtype *k, const int ip, gsl::span<const int> reinitialization_state_idxs) override {
        sx0_fixedParameters_tmp0n6wklye(sx0_fixedParameters, t, x0, p, k, ip,  reinitialization_state_idxs);
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

    virtual void fw(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl) override {
        w_tmp0n6wklye(w, t, x, p, k, h, tcl);
    }


    virtual void fx0(realtype *x0, const realtype t, const realtype *p, const realtype *k) override {
        x0_tmp0n6wklye(x0, t, p, k);
    }


    virtual void fx0_fixedParameters(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs) override {
        x0_fixedParameters_tmp0n6wklye(x0_fixedParameters, t, p, k,  reinitialization_state_idxs);
    }


    virtual void fxdot(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        xdot_tmp0n6wklye(xdot, t, x, p, k, h, w);
    }


    virtual void fy(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        y_tmp0n6wklye(y, t, x, p, k, h, w);
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
        x_solver_tmp0n6wklye(x_solver, x_rdata);
    }


    virtual void ftotal_cl(realtype *total_cl, const realtype *x_rdata) override {}


    std::string getName() const override {
        return "tmp0n6wklye";
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
     * @brief Get names of model expressions
     * @return Expression names
     */
    virtual std::vector<std::string> getExpressionNames() const override {
        return std::vector<std::string>(expressionNames.begin(),
                                        expressionNames.end());
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
     * @brief Get IDs of model expressions
     * @return Expression IDs
     */
    virtual std::vector<std::string> getExpressionIds() const override {
        return std::vector<std::string>(expressionIds.begin(),
                                        expressionIds.end());
    }

    /**
     * @brief function indicating whether reinitialization of states depending on
     fixed parameters is permissible
     * @return flag indicating whether reinitialization of states depending on
     fixed parameters is permissible
     */
    virtual bool isFixedParameterStateReinitializationAllowed() const override {
        return true;
    }

    /**
     * @brief returns the AMICI version that was used to generate the model
     * @return AMICI version string
     */
    virtual std::string getAmiciVersion() const override {
        return "0.11.16";
    }

    /**
     * @brief returns the amici version that was used to generate the model
     * @return AMICI git commit hash
     */
    virtual std::string getAmiciCommit() const override {
        return "unknown";
    }

    virtual bool hasQuadraticLLH() const override {
        return true;
    }
};


} // namespace model_tmp0n6wklye

} // namespace amici

#endif /* _amici_TPL_MODELNAME_h */

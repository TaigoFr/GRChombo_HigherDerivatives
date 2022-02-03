#ifndef SCHWARZSCHILDISOTROPIC
#define SCHWARZSCHILDISOTROPIC

#include "BSSNVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"

#define MIN_CUT_OFF(VAR, MIN_VAL)                                              \
    VAR = simd_conditional(simd_compare_lt(VAR, MIN_VAL), MIN_VAL, VAR)
#define MAX_CUT_OFF(VAR, MAX_VAL)                                              \
    VAR = simd_conditional(simd_compare_gt(VAR, MAX_VAL), MAX_VAL, VAR)

class SchwarzschildIsotropic
{
  public:
    //! Stuct for the params of the Schw BH
    struct params_t
    {
        double mass;                            //!<< The mass of the Schw BH
        std::array<double, CH_SPACEDIM> center; //!< The center of the Schw BH
    };

  protected:
    double m_dx;
    params_t m_params;

  public:
    SchwarzschildIsotropic(params_t a_params, const double a_dx)
        : m_params(a_params), m_dx(a_dx)
    {
    }

    static double get_bh_radius(params_t a_params)
    {
        return a_params.mass / 2.;
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        BSSNVars::VarsWithGauge<data_t> vars;

        data_t r = coords.get_radius();
        data_t r_reg = MIN_CUT_OFF(
            r,
            get_bh_radius(m_params) /
                10.); // TODO: decide whether this regularisation is necessary

        data_t f = m_params.mass / (2. * r_reg);

        FOR(i) { vars.h[i][i] = 1.; }
        vars.chi = 1. / pow(1. + f, 4.);

        // vars.lapse = sqrt(vars.chi);
        vars.lapse = 1.;
        // vars.lapse = abs(1. - f) / (1. + f);

        current_cell.store_vars(vars);
    }
};
#endif /* SCHWARZSCHILDISOTROPIC*/

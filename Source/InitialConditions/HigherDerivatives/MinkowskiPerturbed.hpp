#ifndef MINKOWSKIPERTURBED_HPP_
#define MINKOWSKIPERTURBED_HPP_

#include "BSSNVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"

class MinkowskiPerturbed
{
  public:
    //! Stuct for the params of the Schw BH
    struct params_t
    {
        double amplitude;
        double r0;
        std::array<double, CH_SPACEDIM> center; //!< The center
    };

  protected:
    double m_dx;
    params_t m_params;

  public:
    MinkowskiPerturbed(params_t a_params, const double a_dx)
        : m_params(a_params), m_dx(a_dx)
    {
    }

    static double get_bh_radius(params_t a_params)
    {
        return a_params.r0; // one full loop of the sine perturbation
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        BSSNVars::VarsWithGauge<data_t> vars;

        data_t r = coords.get_radius();
        data_t r_n = r / m_params.r0;

        data_t pi2r = 2. * M_PI * r_n;
        data_t f = 1. + m_params.amplitude * sin(pi2r) / pi2r;

        vars.lapse = 1. / sqrt(f);
        FOR1(i) { vars.h[i][i] = 1.; }
        vars.chi = 1. / f;

        current_cell.store_vars(vars);
    }
};
#endif /* MINKOWSKIPERTURBED_HPP_ */

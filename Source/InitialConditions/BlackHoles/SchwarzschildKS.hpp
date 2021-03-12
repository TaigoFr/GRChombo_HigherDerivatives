#ifndef SCHWARZSCHILDKS
#define SCHWARZSCHILDKS

#include "BSSNVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"

#define MIN_CUT_OFF(VAR, MIN_VAL)                                              \
    VAR = simd_conditional(simd_compare_lt(VAR, MIN_VAL), MIN_VAL, VAR)
#define MAX_CUT_OFF(VAR, MAX_VAL)                                              \
    VAR = simd_conditional(simd_compare_gt(VAR, MAX_VAL), MAX_VAL, VAR)

class SchwarzschildKS
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
    SchwarzschildKS(params_t a_params, const double a_dx)
        : m_params(a_params), m_dx(a_dx)
    {
    }

    static double get_bh_radius(params_t a_params)
    {
        return a_params.mass * 2.;
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        BSSNVars::VarsWithGauge<data_t> vars;

        data_t r = coords.get_radius();
        data_t r_reg = MIN_CUT_OFF(
            r,
            get_bh_radius(m_params) /
                2.); // TODO: decide whether this regularisation is necessary

        Tensor<1, data_t> xyz;
        xyz[0] = coords.x;
        xyz[1] = coords.y;
        xyz[2] = coords.z;

        data_t f = 2. * m_params.mass / r_reg;

        Tensor<2, data_t> Kij;

        vars.chi = 1. / pow(1. + f, 1. / 3.);
        FOR1(i)
        {
            vars.shift[i] = f * xyz[i] / (r * (1. + f));
            FOR1(j)
            {
                vars.h[i][j] = f * xyz[i] * xyz[j] * vars.chi / (r_reg * r_reg);
                Kij[i][j] = -f * (4. + f) * xyz[i] * xyz[j] /
                            (2. * sqrt(1. + f) * r_reg * r_reg * r_reg);
            }

            vars.h[i][i] += vars.chi;
            Kij[i][i] += f / (sqrt(1. + f) * r_reg);

            vars.Gamma[i] = 2. * f * (2. + 3. * f) * xyz[i] *
                            pow(vars.chi, 5.) / (3. * r_reg * r_reg);
        }

        auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
        vars.K = vars.chi * TensorAlgebra::compute_trace(Kij, h_UU);
        FOR2(i, j)
        {
            vars.A[i][j] =
                vars.chi * Kij[i][j] - vars.K * vars.h[i][j] / GR_SPACEDIM;
        }

        // vars.lapse = sqrt(vars.chi);
        vars.lapse = 1.;
        // vars.lapse = 1. / sqrt(1. + f);

        current_cell.store_vars(vars);
    }
};
#endif /* SCHWARZSCHILDKS*/

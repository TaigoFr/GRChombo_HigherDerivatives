/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPUTECDOT
#define COMPUTECDOT

#include "C2EFT.hpp"
#include "CSystem.hpp"
#include "Cell.hpp"
#include "EBSystem.hpp"
#include "FourthOrderDerivatives.hpp"
#include "MatterCCZ4RHS.hpp"
#include "MovingPunctureGauge.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components

class ComputeCdot
{

    // Use the variable definitions in MatterCCZ4RHS
    template <class data_t>
    using Vars = typename MatterCCZ4RHS<C2EFT<EBSystem>>::template Vars<data_t>;

    template <class data_t>
    using Diff2Vars =
        typename MatterCCZ4RHS<C2EFT<EBSystem>>::template Diff2Vars<data_t>;

  public:
    ComputeCdot(double m_dx, int a_formulation,
                const CCZ4_params_t<> &a_ccz4_params,
                bool use_last_index_raised, int a_C_dot_comp)
        : m_formulation(a_formulation), m_ccz4_params(a_ccz4_params),
          m_deriv(m_dx), m_use_last_index_raised(use_last_index_raised),
          m_C_dot_comp(a_C_dot_comp)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto vars = current_cell.template load_vars<Vars>();

        /*
        double m_dx = m_deriv.m_dx;
        Coordinates<data_t> coords(current_cell, m_dx, {64., 0., 0});

        if (m_dx < 1)
            pout() << coords << std::endl;
        if (m_dx < 1 && abs(coords.x - 0 - 5) < m_dx &&
            abs(coords.y - 0) < m_dx)
        {
            if (m_use_last_index_raised)
            {
                pout() << coords << "; Eij_LU[0][0]=" << vars.Eij[0][0]
                       << "; dt_Eij_LU[0][0]=" << vars.Eaux[0][0]
                       << "; Bij_LU[0][0]=" << vars.Bij[0][0]
                       << "; dt_Bij_LU[0][0]=" << vars.Baux[0][0]
                       << "; Eij_LU[0][1]=" << vars.Eij[0][1]
                       << "; dt_Eij_LU[0][1]=" << vars.Eaux[0][1]
                       << "; Bij_LU[0][1]=" << vars.Bij[0][1]
                       << "; dt_Bij_LU[0][1]=" << vars.Baux[0][1]
                       << "; Eij_LU[1][1]=" << vars.Eij[1][1]
                       << "; dt_Eij_LU[1][1]=" << vars.Eaux[1][1]
                       << "; Bij_LU[1][1]=" << vars.Bij[1][1]
                       << "; dt_Bij_LU[1][1]=" << vars.Baux[1][1] << std::endl;
            }
            else
            {
                const auto d1 = m_deriv.template diff1<Vars>(current_cell);
                const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);
                const auto advec =
                    m_deriv.template advection<Vars>(current_cell, vars.shift);
                GeometricQuantities<data_t, Vars, Diff2Vars> gq(
                    vars, d1, d2, "ComputeCdot::compute");
                gq.set_formulation(m_formulation, m_ccz4_params);
                gq.set_advection_and_gauge(advec, EmptyGauge());
                const auto &h_UU = gq.get_metric_UU_spatial();

                Tensor<2, data_t> dt_Eij_LU = {0.};
                Tensor<2, data_t> dt_Bij_LU = {0.};

                Tensor<2, data_t> Eij_LU =
                    TensorAlgebra::compute_dot_product(vars.Eij, h_UU);
                Tensor<2, data_t> Bij_LU =
                    TensorAlgebra::compute_dot_product(vars.Bij, h_UU);

                Tensor<2, data_t> Eij_no_time_der = vars.Eij;
                Tensor<2, data_t> Bij_no_time_der = vars.Bij;

                Tensor<2, data_t> Eij_no_time_der_LU =
                    TensorAlgebra::compute_dot_product(Eij_no_time_der, h_UU);
                Tensor<2, data_t> Bij_no_time_der_LU =
                    TensorAlgebra::compute_dot_product(Bij_no_time_der, h_UU);

                Vars<data_t> rhs;
                gq.compute_rhs_equations_no_gauge(rhs);

                Tensor<2, data_t> dt_h;
                data_t chi2 = vars.chi * vars.chi;
                FOR(i, j)
                {
                    dt_h[i][j] =
                        rhs.h[i][j] / vars.chi - rhs.chi * vars.h[i][j] / chi2;
                }

                // (recall Eij and Bij are used to store the time derivatives LU
                // when 'm_compute_time_derivatives == true')
                FOR(i, j, k)
                {
                    dt_Eij_LU[i][j] += vars.Eaux[i][k] * h_UU[k][j];
                    dt_Bij_LU[i][j] += vars.Baux[i][k] * h_UU[k][j];
                    FOR(l)
                    {
                        dt_Eij_LU[i][j] +=
                            -Eij_no_time_der_LU[i][l] * h_UU[j][k] * dt_h[k][l];
                        dt_Bij_LU[i][j] +=
                            -Bij_no_time_der_LU[i][l] * h_UU[j][k] * dt_h[k][l];
                    }
                }

                pout() << coords << "; Eij_LU[0][0]=" << Eij_LU[0][0]
                       << "; dt_Eij_LU[0][0]=" << dt_Eij_LU[0][0]
                       << "; Bij_LU[0][0]=" << Bij_LU[0][0]
                       << "; dt_Bij_LU[0][0]=" << dt_Bij_LU[0][0]
                       << "; Eij_LU[0][1]=" << Eij_LU[0][1]
                       << "; dt_Eij_LU[0][1]=" << dt_Eij_LU[0][1]
                       << "; Bij_LU[0][1]=" << Bij_LU[0][1]
                       << "; dt_Bij_LU[0][1]=" << dt_Bij_LU[0][1]
                       << "; Eij_LU[1][1]=" << Eij_LU[1][1]
                       << "; dt_Eij_LU[1][1]=" << dt_Eij_LU[1][1]
                       << "; Bij_LU[1][1]=" << Bij_LU[1][1]
                       << "; dt_Bij_LU[1][1]=" << dt_Bij_LU[1][1] << std::endl;

                data_t C_dot2 = 0.;
                FOR(i, j)
                {
                    C_dot2 += 16. * (Eij_LU[i][j] * dt_Eij_LU[j][i] -
                                     Bij_LU[i][j] * dt_Bij_LU[j][i]);
                }
                pout() << "C_dot2 =" << C_dot2 << std::endl;
            }
        }
        */

        data_t C_dot = 0.;
        if (m_use_last_index_raised)
        {
            FOR(i, j)
            {
                C_dot += 16. * (vars.Eij[i][j] * vars.Eaux[j][i] -
                                vars.Bij[i][j] * vars.Baux[j][i]);
            }
        }
        else
        {
            const auto d1 = m_deriv.template diff1<Vars>(current_cell);
            const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);
            const auto advec =
                m_deriv.template advection<Vars>(current_cell, vars.shift);

            GeometricQuantities<data_t, Vars, Diff2Vars> gq(
                vars, d1, d2, "ComputeCdot::compute");
            gq.set_formulation(m_formulation, m_ccz4_params);
            gq.set_advection_and_gauge(advec, EmptyGauge());

            // raise last index
            const auto &h_UU = gq.get_metric_UU_spatial();
            Tensor<2, data_t> Eij_LU =
                TensorAlgebra::compute_dot_product(vars.Eij, h_UU);
            Tensor<2, data_t> Bij_LU =
                TensorAlgebra::compute_dot_product(vars.Bij, h_UU);
            // raise first index
            Tensor<2, data_t> Eij_UU =
                TensorAlgebra::compute_dot_product(h_UU, Eij_LU);
            Tensor<2, data_t> Bij_UU =
                TensorAlgebra::compute_dot_product(h_UU, Bij_LU);

            // now compute time derivative of spatial (not conformal!) metric
            Vars<data_t> rhs;
            gq.compute_rhs_equations_no_gauge(rhs);

            Tensor<2, data_t> dt_h;
            data_t chi2 = vars.chi * vars.chi;
            FOR(i, j)
            {
                dt_h[i][j] =
                    rhs.h[i][j] / vars.chi - rhs.chi * vars.h[i][j] / chi2;
            }

            FOR(i, j)
            {
                C_dot += 16. * (Eij_UU[i][j] * vars.Eaux[i][j] -
                                Bij_UU[i][j] * vars.Baux[i][j]);

                FOR(k)
                {
                    C_dot += -16. *
                             (Eij_UU[i][j] * Eij_LU[i][k] -
                              Bij_UU[i][j] * Bij_LU[i][k]) *
                             dt_h[j][k];
                }
            }
        }

        current_cell.store_vars(C_dot, m_C_dot_comp);

        /*
        if (m_dx < 1 && abs(coords.x - 0 - 5) < m_dx &&
            abs(coords.y - 0) < m_dx)
        {
            pout() << "C_dot =" << C_dot << std::endl;
            MayDay::Error("STOP HERE");
        }
        */
    }

  protected:
    int m_formulation;
    const CCZ4_params_t<> &m_ccz4_params;
    FourthOrderDerivatives m_deriv;

    bool m_use_last_index_raised;
    int m_C_dot_comp;
};

#endif /* COMPUTECDOT */

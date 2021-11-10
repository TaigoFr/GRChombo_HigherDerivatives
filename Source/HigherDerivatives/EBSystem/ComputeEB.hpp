/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPUTEEB
#define COMPUTEEB

#include "C2EFT.hpp"
#include "Cell.hpp"
#include "EBSystem.hpp"
#include "FourthOrderDerivatives.hpp"
#include "MatterCCZ4RHS.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components

class ComputeEB
{

    // Use the variable definitions in MatterCCZ4RHS
    template <class data_t>
    using Vars = typename MatterCCZ4RHS<C2EFT<EBSystem>>::template Vars<data_t>;

    template <class data_t>
    using Diff2Vars =
        typename MatterCCZ4RHS<C2EFT<EBSystem>>::template Diff2Vars<data_t>;

  public:
    ComputeEB(double m_dx, int a_formulation,
              const CCZ4_params_t<> &a_ccz4_params, const Interval &E_comps,
              const Interval &B_comps, bool use_last_index_raised,
              bool compute_time_derivatives = false)
        : m_formulation(a_formulation), m_ccz4_params(a_ccz4_params),
          m_deriv(m_dx), m_E_comps(E_comps), m_B_comps(B_comps),
          m_compute_time_derivatives(compute_time_derivatives),
          m_use_last_index_raised(use_last_index_raised)
    {
        CH_assert(m_E_comps.size() == (GR_SPACEDIM * (GR_SPACEDIM + 1) / 2));
        CH_assert(m_B_comps.size() == (GR_SPACEDIM * (GR_SPACEDIM + 1) / 2));
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);

        GeometricQuantities<data_t, Vars, Diff2Vars> gq(vars, d1, d2,
                                                        "ComputeEB::compute");

        gq.set_formulation(m_formulation, m_ccz4_params);

        // Eij and Bij are a bit everything here
        // they will be LL or LU depending on 'm_use_last_index_raised'
        // they will be the actual tensors OR THEIR TIME DERIVATIVE if
        // 'm_compute_time_derivatives == true'
        Tensor<2, data_t> Eij = {0.};
        Tensor<2, data_t> Bij = {0.};
        const auto &h_UU = gq.get_metric_UU_spatial();

        if (m_compute_time_derivatives)
        {
            // this is a pain
            // to compute the time derivative of LL via
            // GeometricQuantities::compute_dt_weyl_electric_part, so all
            // vars.E, d1.E, advec.E (same for B) that have one index raised
            // need to be lowered to have the right input (and raised
            // afterwards)
            const auto advec =
                m_deriv.template advection<Vars>(current_cell, vars.shift);

            Tensor<2, Tensor<1, data_t>> d1_Eij_LL, d1_Bij_LL;
            Tensor<2, data_t> Eij_LL, Bij_LL, advec_Eij_LL, advec_Bij_LL;

            if (m_use_last_index_raised)
            {
                const auto &metric_spatial = gq.get_metric_spatial();

                Tensor<2, data_t> advec_h;
                Tensor<2, Tensor<1, data_t>> d1_h;
                data_t chi2 = vars.chi * vars.chi;
                FOR(i, j)
                {
                    advec_h[i][j] = advec.h[i][j] / vars.chi -
                                    advec.chi * vars.h[i][j] / chi2;
                    FOR(k)
                    {
                        d1_h[i][j][k] = d1.h[i][j][k] / vars.chi -
                                        d1.chi[k] * vars.h[i][j] / chi2;
                    }
                }

                FOR(i, j)
                {
                    advec_Eij_LL[i][j] = 0.;
                    advec_Bij_LL[i][j] = 0.;

                    FOR(k)
                    {
                        advec_Eij_LL[i][j] +=
                            advec.Eij[i][k] * metric_spatial[k][j] +
                            vars.Eij[i][k] * advec_h[k][j];
                        advec_Bij_LL[i][j] +=
                            advec.Bij[i][k] * metric_spatial[k][j] +
                            vars.Bij[i][k] * advec_h[k][j];

                        d1_Eij_LL[i][j][k] = 0.;
                        d1_Bij_LL[i][j][k] = 0.;
                        FOR(l)
                        {
                            d1_Eij_LL[i][j][k] +=
                                d1.Eij[i][l][k] * metric_spatial[l][j] +
                                vars.Eij[i][l] * d1_h[l][j][k];
                            d1_Bij_LL[i][j][k] +=
                                d1.Bij[i][l][k] * metric_spatial[l][j] +
                                vars.Bij[i][l] * d1_h[l][j][k];
                        }
                    }
                }
                Eij_LL = TensorAlgebra::compute_dot_product(vars.Eij,
                                                            metric_spatial);
                Bij_LL = TensorAlgebra::compute_dot_product(vars.Bij,
                                                            metric_spatial);
            }
            else
            {
                d1_Eij_LL = d1.Eij;
                d1_Bij_LL = d1.Bij;
                Eij_LL = vars.Eij;
                Bij_LL = vars.Bij;
                advec_Eij_LL = advec.Eij;
                advec_Bij_LL = advec.Bij;
            }

            Tensor<2, data_t> dt_Eij = gq.compute_dt_weyl_electric_part(
                d1_Eij_LL, d1_Bij_LL, Eij_LL, Bij_LL, advec_Eij_LL,
                advec_Bij_LL);
            Tensor<2, data_t> dt_Bij = gq.compute_dt_weyl_magnetic_part(
                d1_Eij_LL, d1_Bij_LL, Eij_LL, Bij_LL, advec_Eij_LL,
                advec_Bij_LL);

            if (m_use_last_index_raised)
            {
                // now need to raise them back to get the time derivative LU
                // this needs the rhs of 'h', dt_h

                // gauge does not matter as we only compute
                // 'compute_rhs_equations_no_gauge' (we don't use the full RHS
                // equations)
                gq.set_advection_and_gauge(advec, EmptyGauge());

                Tensor<2, data_t> Eij_no_time_der = gq.get_weyl_electric_part();
                Tensor<2, data_t> Bij_no_time_der = gq.get_weyl_magnetic_part();

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
                    Eij[i][j] += dt_Eij[i][k] * h_UU[k][j];
                    Bij[i][j] += dt_Bij[i][k] * h_UU[k][j];
                    FOR(l)
                    {
                        Eij[i][j] +=
                            -Eij_no_time_der_LU[i][l] * h_UU[j][k] * dt_h[k][l];
                        Bij[i][j] +=
                            -Bij_no_time_der_LU[i][l] * h_UU[j][k] * dt_h[k][l];
                    }
                }
            }
            else
            {
                // replace tensors
                // (recall Eij and Bij are used to store the time derivatives
                // when 'm_compute_time_derivatives == true')
                Eij = dt_Eij;
                Bij = dt_Bij;
            }
        }
        else
        {
            Eij = gq.get_weyl_electric_part();
            Bij = gq.get_weyl_magnetic_part();

            if (m_use_last_index_raised)
            {
                // raise last index
                Eij = TensorAlgebra::compute_dot_product(Eij, h_UU);
                Bij = TensorAlgebra::compute_dot_product(Bij, h_UU);
            }
        }

        current_cell.store_vars(Eij[0][0], m_E_comps.begin() + 0);
        current_cell.store_vars(Eij[0][1], m_E_comps.begin() + 1);
        current_cell.store_vars(Eij[0][2], m_E_comps.begin() + 2);
        current_cell.store_vars(Eij[1][1], m_E_comps.begin() + 3);
        current_cell.store_vars(Eij[1][2], m_E_comps.begin() + 4);
        current_cell.store_vars(Eij[2][2], m_E_comps.begin() + 5);

        current_cell.store_vars(Bij[0][0], m_B_comps.begin() + 0);
        current_cell.store_vars(Bij[0][1], m_B_comps.begin() + 1);
        current_cell.store_vars(Bij[0][2], m_B_comps.begin() + 2);
        current_cell.store_vars(Bij[1][1], m_B_comps.begin() + 3);
        current_cell.store_vars(Bij[1][2], m_B_comps.begin() + 4);
        current_cell.store_vars(Bij[2][2], m_B_comps.begin() + 5);
    }

  protected:
    int m_formulation;
    const CCZ4_params_t<> &m_ccz4_params;
    FourthOrderDerivatives m_deriv;
    const Interval m_E_comps, m_B_comps;
    bool m_compute_time_derivatives;
    bool m_use_last_index_raised;
};

#endif /* COMPUTEEB */

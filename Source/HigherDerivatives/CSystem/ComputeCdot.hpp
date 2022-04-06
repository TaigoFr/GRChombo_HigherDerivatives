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
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);

        GeometricQuantities<data_t, Vars, Diff2Vars> gq(vars, d1, d2,
                                                        "ComputeCdot::compute");

        gq.set_formulation(m_formulation, m_ccz4_params);

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
                    C_dot += 16. *
                             (Eij_UU[i][j] * Eij_LU[i][k] -
                              Bij_UU[i][j] * Bij_LU[i][k]) *
                             dt_h[j][k];
                }
            }
        }

        current_cell.store_vars(C_dot, m_C_dot_comp);
    }

  protected:
    int m_formulation;
    const CCZ4_params_t<> &m_ccz4_params;
    FourthOrderDerivatives m_deriv;

    bool m_use_last_index_raised;
    int m_C_dot_comp;
};

#endif /* COMPUTECDOT */

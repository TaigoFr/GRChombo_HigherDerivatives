/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(C2EFT_HPP_)
#error "This file should only be included through C2EFT.hpp"
#endif

#ifndef C2EFT_IMPL_HPP_
#define C2EFT_IMPL_HPP_
#include "DimensionDefinitions.hpp"

// Calculate the stress energy tensor elements
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
emtensor_t<data_t> C2EFT::compute_emtensor(
    GeometricQuantities<data_t, vars_t, diff2_vars_t> &gq) const
{
    emtensor_t<data_t> out;

    Tensor<2, data_t, CH_SPACEDIM + 1> Tmn; // 4D

    compute_emtensor_4D(Tmn, gq);

    const auto &n_U = gq.get_normal_U_ST();
    const auto &proj_LU = gq.get_projector_LU_ST();

    Tensor<1, data_t, CH_SPACEDIM + 1> Tmn_dot_normal_U =
        TensorAlgebra::compute_dot_product(Tmn, n_U);

    Tensor<1, data_t, CH_SPACEDIM + 1> Si_4D =
        TensorAlgebra::compute_dot_product(Tmn_dot_normal_U, proj_LU);

    Tensor<2, data_t, CH_SPACEDIM + 1> Sij_4D =
        TensorAlgebra::compute_dot_product(
            proj_LU, TensorAlgebra::compute_dot_product(Tmn, proj_LU), 1, 0);

    out.rho = TensorAlgebra::compute_dot_product(Tmn_dot_normal_U, n_U);
    FOR(i)
    {
        out.Si[i] = -Si_4D[i + 1];
        FOR(j) { out.Sij[i][j] = Sij_4D[i + 1][j + 1]; }
    }
    out.S = TensorAlgebra::compute_trace(out.Sij, gq.get_metric_UU_spatial());

    return out;
}

template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
void C2EFT::compute_emtensor_4D(
    Tensor<2, data_t, CH_SPACEDIM + 1> &Tmn,
    GeometricQuantities<data_t, vars_t, diff2_vars_t> &gq) const
{
    const auto &vars = gq.get_vars();
    const auto &d1 = gq.get_d1_vars();
    const auto &d2 = gq.get_d2_vars();
    const auto &g = gq.get_metric_ST();
    const auto &riemann_LLLL = gq.get_riemann_LLLL_ST();
    const auto &riemann_LLLU = gq.get_riemann_LLLU_ST();
    const auto &riemann_LULU = gq.get_riemann_LULU_ST();
    const auto &chris_ST = gq.get_chris_ST();

    Tensor<1, data_t, CH_SPACEDIM + 1> d1_C_4D;
    Tensor<2, data_t, CH_SPACEDIM + 1> d2_C_4D;

    Vars<data_t> rhs;
    add_matter_rhs(rhs, gq);

    d1_C_4D[0] = vars.dCdt;
    d2_C_4D[0][0] = rhs.dCdt;
    FOR(i)
    {
        d1_C_4D[i + 1] = d1.C[i];

        d2_C_4D[0][i + 1] = d1.dCdt[i];
        d2_C_4D[i + 1][0] = d1.dCdt[i];

        FOR(j) { d2_C_4D[i + 1][j + 1] = d2.C[i][j]; }
    }

    const Tensor<2, data_t, CH_SPACEDIM + 1> covd2_C =
        TensorAlgebra::covariant_derivative(d2_C_4D, d1_C_4D, chris_ST);

    FOR_ST(a, b)
    {
        Tmn[a][b] = -0.5 * vars.C * vars.C * g[a][b];
        FOR_ST(c, d)
        {
            Tmn[a][b] += 8. * riemann_LULU[a][c][b][d] * covd2_C[c][d];

            FOR_ST(e)
            {
                Tmn[a][b] += -4. * vars.C * riemann_LULU[a][c][d][e] *
                             riemann_LLLU[b][c][e][d];
            }
        }
    }
}

// Adds in the RHS for the matter vars
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void C2EFT::add_matter_rhs(
    rhs_vars_t<data_t> &total_rhs,
    GeometricQuantities<data_t, vars_t, diff2_vars_t> &gq) const
{
    const auto &vars = gq.get_vars();
    const auto &advec = gq.get_advection();
    const auto &kretschmann = gq.get_kretschmann();

    CH_assert(m_params.sigma != 0.);
    // if (m_params.sigma != 0.)
    {
        const auto &d1 = gq.get_d1_vars();
        const auto &d2 = gq.get_d2_vars();
        const auto &g_UU = gq.get_metric_UU_ST();
        const auto &Gamma_ST = gq.get_Gamma_ST();

        total_rhs.C = vars.dCdt;
        total_rhs.dCdt = -m_params.tau / vars.lapse * (vars.dCdt - advec.C) +
                         kretschmann - vars.C -
                         m_params.sigma * Gamma_ST[0] * vars.dCdt;

        FOR(i)
        {
            total_rhs.dCdt +=
                m_params.sigma *
                (2. * g_UU[0][i + 1] * d1.dCdt[i] - Gamma_ST[i + 1] * d1.C[i]);
            FOR(j)
            {
                total_rhs.dCdt +=
                    m_params.sigma * g_UU[i + 1][j + 1] * d2.C[i][j];
            }
        }

        total_rhs.dCdt /= (-m_params.sigma * g_UU[0][0]);
    }
    // else
    // {
    //     total_rhs.C =
    //         (kretschmann - vars.C) * vars.lapse / m_params.tau + advec.C;
    //     total_rhs.dCdt = 0.;
    // }
}

#endif /* C2EFT_IMPL_HPP_ */

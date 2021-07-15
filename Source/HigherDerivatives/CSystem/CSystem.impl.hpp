/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(CSYSTEM_HPP_)
#error "This file should only be included through CSystem.hpp"
#endif

#ifndef CSYSTEM_IMPL_HPP_
#define CSYSTEM_IMPL_HPP_
#include "DimensionDefinitions.hpp"

// Calculate the stress energy tensor elements
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t, class gauge_t>
void CSystem::compute_C(
    data_t &C, Tensor<1, data_t, CH_SPACEDIM + 1> &d1_C,
    Tensor<2, data_t, CH_SPACEDIM + 1> &d2_C,
    GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq) const
{
    const auto &vars = gq.get_vars();
    const auto &d1 = gq.get_d1_vars();
    const auto &d2 = gq.get_d2_vars();

    C = vars.C;

    Vars<data_t> rhs;
    add_matter_rhs(rhs, gq);

    d1_C[0] = vars.dCdt;
    d2_C[0][0] = rhs.dCdt;
    FOR(i)
    {
        d1_C[i + 1] = d1.C[i];

        d2_C[0][i + 1] = d1.dCdt[i];
        d2_C[i + 1][0] = d1.dCdt[i];

        FOR(j) { d2_C[i + 1][j + 1] = d2.C[i][j]; }
    }
}

template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t, class gauge_t>
void CSystem::compute_Riemann(
    Tensor<4, data_t, CH_SPACEDIM + 1> &riemann_LLLU,
    Tensor<4, data_t, CH_SPACEDIM + 1> &riemann_LULU,
    GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq) const
{
    riemann_LLLU = gq.get_riemann_LLLU_ST();
    riemann_LULU = gq.get_riemann_LULU_ST();
}

// Adds in the RHS for the matter vars
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t, class gauge_t,
          template <typename> class rhs_vars_t>
void CSystem::add_matter_rhs(
    rhs_vars_t<data_t> &total_rhs,
    GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq) const
{
    const auto &vars = gq.get_vars();
    const auto &kretschmann = gq.get_kretschmann();

    CH_assert(m_params.sigma != 0.);

    total_rhs.C = vars.dCdt;

    if (m_params.use_only_time_derivatives)
    {
        data_t tau = m_params.tau;
        data_t sigma = m_params.sigma;
        data_t dCdt = vars.dCdt;

        if (m_params.rescale_tau_by_lapse)
            tau /= vars.lapse;
        if (m_params.rescale_sigma_by_lapse == 2)
            sigma /= (vars.lapse * vars.lapse);
        else if (m_params.rescale_sigma_by_lapse == 1)
            sigma /= vars.lapse;

        if (m_params.add_advection)
        {
            const auto &advec = gq.get_advection();
            dCdt -= advec.C;
        }
        total_rhs.dCdt = (-tau * dCdt + kretschmann - vars.C) / sigma;
    }
    else
    {
        const auto &d1 = gq.get_d1_vars();
        const auto &d2 = gq.get_d2_vars();
        const auto &advec = gq.get_advection();
        const auto &g_UU = gq.get_metric_UU_ST();
        const auto &Gamma_ST = gq.get_Gamma_ST();

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
}

template <class data_t, template <typename> class rhs_vars_t,
          template <typename> class vars_t,
          template <typename> class diff2_vars_t, class gauge_t>
void CSystem::add_diffusion_terms(
    rhs_vars_t<data_t> &rhs, //!< Reference to the variables into which the
                             //! output right hand side is written
    GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq,
    data_t diffCoeffSafe) const
{
    const auto &d2 = gq.get_d2_vars();

    data_t space_laplace_C = 0.;
    FOR(k) { space_laplace_C += d2.C[k][k]; }

    rhs.C += diffCoeffSafe * space_laplace_C;
}

#endif /* CSYSTEM_IMPL_HPP_ */

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
    data_t &C, Tensor<1, data_t, CH_SPACETIMEDIM> &d1_C,
    Tensor<2, data_t, CH_SPACETIMEDIM> &d2_C,
    GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq,
    const C2EFT<CSystem>::params_t &pm) const
{
    const auto &vars = gq.get_vars();
    const auto &d1 = gq.get_d1_vars();
    const auto &d2 = gq.get_d2_vars();

    C = vars.C;

    Vars<data_t> rhs;
    add_matter_rhs(rhs, gq, pm);

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
    Tensor<4, data_t, CH_SPACETIMEDIM> &riemann_LLLU,
    Tensor<4, data_t, CH_SPACETIMEDIM> &riemann_LULU,
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
    GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq,
    const C2EFT<CSystem>::params_t &pm) const
{
    const auto &vars = gq.get_vars();
    const auto &advec = gq.get_advection();
    const auto &kretschmann = gq.get_kretschmann();

    CH_assert(m_params.sigma != 0.);

    total_rhs.C = vars.dCdt;

    data_t tau = m_params.tau;
    data_t sigma = m_params.sigma;

    if (m_params.use_tau_radial_decay)
    {
        const Coordinates<data_t> &coords = gq.get_coordinates();
        tau = tau_from_radius(coords.get_radius());
    }
    
    if (m_params.use_tau_radial_decay)
    {
        const Coordinates<data_t> &coords = gq.get_coordinates();
        sigma = sigma_from_radius(coords.get_radius());
    }

    if (m_params.use_tau_chi_decay)
    {
       tau = (m_params.tau_asymptotic-m_params.tau)/(1.0 +exp(-(vars.chi/m_params.tau_decay_length -1.0)/m_params.tau_decay_width)) + m_params.tau;
    }    

   if (m_params.use_sigma_chi_decay)
    {
       sigma = (m_params.sigma_asymptotic-m_params.sigma)/(1.0 +exp(-(vars.chi/m_params.sigma_decay_length -1.0)/m_params.sigma_decay_width) ) + m_params.sigma;
    }


    data_t tau_rescaled = tau;
    data_t sigma_rescaled = sigma;
    if (m_params.rescale_tau_by_lapse)
        tau_rescaled /= vars.lapse;
    if (m_params.rescale_sigma_by_lapse == 2)
        sigma_rescaled /= (vars.lapse * vars.lapse);
    else if (m_params.rescale_sigma_by_lapse == 1)
        sigma_rescaled /= vars.lapse;

    if (m_params.version == 1)
    {
        data_t advfac = m_params.advection_type ? 1.0 : 0.0;

        data_t factor = 1.0;
        if (m_params.Box_transition)
        {
            factor = 0.0;
            data_t weak_field_var = 10.0;
            factor =
                C2EFT<CSystem>::weak_field_condition(weak_field_var, gq, pm);
        }
        const auto &d1 = gq.get_d1_vars();
        const auto &d2 = gq.get_d2_vars();
        const auto &g_UU = gq.get_metric_UU_ST();
        const auto &Gamma_ST = gq.get_Gamma_ST();

        total_rhs.dCdt =
            -tau / vars.lapse * (vars.dCdt - advfac * factor * advec.C) +
            kretschmann - vars.C - factor * sigma * Gamma_ST[0] * vars.dCdt;

        FOR(i)
        {
            total_rhs.dCdt +=
                factor * sigma *
                (2. * g_UU[0][i + 1] * d1.dCdt[i] - Gamma_ST[i + 1] * d1.C[i]);
            FOR(j)
            {
                total_rhs.dCdt +=
                    factor * sigma * g_UU[i + 1][j + 1] * d2.C[i][j];
            }
        }
        total_rhs.dCdt /= (-sigma * g_UU[0][0]);
    }
    else if (m_params.version == 2)
    {
        data_t dCdt = vars.dCdt;

        if (m_params.advection_type == 1 || m_params.advection_type == 2)
        {
            dCdt -= m_params.advection_coeff * advec.C;
        }
        total_rhs.dCdt =
            (-tau_rescaled * dCdt + kretschmann - vars.C) / sigma_rescaled;

        if (m_params.advection_type == 2)
        {
            const auto &d2 = gq.get_d2_vars();

            total_rhs.dCdt += 2. * m_params.advection_coeff * advec.dCdt;

            FOR(i, j)
            {
                total_rhs.dCdt += -m_params.advection_coeff *
                                  m_params.advection_coeff * vars.shift[i] *
                                  vars.shift[j] * d2.C[i][j];
            }
        }
    }
    else
    {
        MayDay::Error("Version not implemented");
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

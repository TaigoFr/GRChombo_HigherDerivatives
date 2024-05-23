/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(EBSYSTEM_HPP_)
#error "This file should only be included through EBSystem.hpp"
#endif

#ifndef EBSYSTEM_IMPL_HPP_
#define EBSYSTEM_IMPL_HPP_

#include "DimensionDefinitions.hpp"

// Calculate the stress energy tensor elements
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t, class gauge_t>
void EBSystem::compute_C(
    data_t &C, Tensor<1, data_t, CH_SPACETIMEDIM> &d1_C,
    Tensor<2, data_t, CH_SPACETIMEDIM> &d2_C,
    GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq,
    const C2EFT<EBSystem>::params_t &pm) const
{
    const auto &vars = gq.get_vars();
    const auto &d1 = gq.get_d1_vars();
    const auto &d2 = gq.get_d2_vars();
    const auto &Kij = gq.get_extrinsic_curvature();
    // needed if 'use_last_index_raised' = false
    const auto &metric_UU_spatial = gq.get_metric_UU_spatial();

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // 0th order first
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Tensor<2, data_t> Eij_LU, Bij_LU;

    if (m_params.use_last_index_raised)
    {
        // the evolved variables are already LU
        Eij_LU = vars.Eij;
        Bij_LU = vars.Bij;
    }
    else
    {
        // raise the last index
        Eij_LU =
            TensorAlgebra::compute_dot_product(vars.Eij, metric_UU_spatial);
        Bij_LU =
            TensorAlgebra::compute_dot_product(vars.Bij, metric_UU_spatial);
    }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Now 1st order
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // indices are only spatial, derivatives are 4D
    Tensor<2, Tensor<1, data_t, CH_SPACETIMEDIM>> d1_Eij, d1_Bij, d1_Eij_LU,
        d1_Bij_LU, d1_h, d1_h_dot_hUU;
    Tensor<2, Tensor<1, data_t, CH_SPACETIMEDIM>> d1_Kij; // needed for d2_h
    // 2nd mixed derivatives of shift needed for d2_h[i][j][0][0]
    Tensor<2, data_t> dtd1_shift; // mixed derivative

    vars_t<data_t> rhs;
    add_matter_rhs(rhs, gq, pm);
    gq.compute_rhs_equations(rhs);

    // first compute d1_h and d1_Kij
    // here 'h' is the spatial metric, not the conformal one
    auto &ccz4_params = gq.get_formulation_params();
    data_t chi2 = vars.chi * vars.chi;
    FOR(i, j)
    {
        d1_h[i][j][0] = rhs.h[i][j] / vars.chi - rhs.chi * vars.h[i][j] / chi2;
        d1_Kij[i][j][0] =
            -Kij[i][j] * rhs.chi / vars.chi +
            (rhs.A[i][j] +
             (vars.h[i][j] * rhs.K + vars.K * rhs.h[i][j]) / GR_SPACEDIM) /
                vars.chi;

        dtd1_shift[i][j] = ccz4_params.shift_Gamma_coeff * d1.B[i][j];

        FOR(k)
        {
            d1_h[i][j][k + 1] =
                d1.h[i][j][k] / vars.chi - d1.chi[k] * vars.h[i][j] / chi2;
            d1_Kij[i][j][k + 1] = -Kij[i][j] * d1.chi[k] / vars.chi +
                                  (d1.A[i][j][k] + (vars.h[i][j] * d1.K[k] +
                                                    vars.K * d1.h[i][j][k]) /
                                                       GR_SPACEDIM) /
                                      vars.chi;

            dtd1_shift[i][j] += ccz4_params.shift_advec_coeff *
                                (vars.shift[k] * d2.shift[i][k][j] +
                                 d1.shift[k][j] * d1.shift[i][k]);
        }
    }

    // compute d1_h_dot_hUU
    FOR(i, j)
    {
        FOR_ST(a)
        {
            d1_h_dot_hUU[i][j][a] = 0.;
            FOR(k)
            {
                d1_h_dot_hUU[i][j][a] +=
                    d1_h[i][k][a] * metric_UU_spatial[k][j];
            }
        }
    }

    // d1_Eij, d1_Bij refer to the evolved ones (could be raised or lowered
    // depending on the option)
    FOR(i, j)
    {
        d1_Eij[i][j][0] = rhs.Eij[i][j]; 
        d1_Bij[i][j][0] = rhs.Bij[i][j]; 

        FOR(k)
        {
            d1_Eij[i][j][k + 1] = d1.Eij[i][j][k];
            d1_Bij[i][j][k + 1] = d1.Bij[i][j][k];
        }
    }

    // now d1_Eij_LU and d1_Bij_LU
    if (m_params.use_last_index_raised)
    {
        // already LU
        d1_Eij_LU = d1_Eij;
        d1_Bij_LU = d1_Bij;
    }
    else
    {
        // raise last index, that needs 'd1_h_dot_hUU'
        FOR_ST(a)
        {
            FOR(i, j)
            {
                d1_Eij_LU[i][j][a] = 0.;
                d1_Bij_LU[i][j][a] = 0.;

                FOR(k)
                {
                    d1_Eij_LU[i][j][a] +=
                        d1_Eij[i][k][a] * metric_UU_spatial[k][j] -
                        Eij_LU[i][k] * d1_h_dot_hUU[k][j][a];
                    d1_Bij_LU[i][j][a] +=
                        d1_Bij[i][k][a] * metric_UU_spatial[k][j] -
                        Bij_LU[i][k] * d1_h_dot_hUU[k][j][a];
                }
            }
        }
    }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Finally 2nd order
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // indices are only spatial, derivatives are 4D
    // d2_Eij will be 'LU' if 'use_last_index_raised' = true
    Tensor<2, Tensor<2, data_t, CH_SPACETIMEDIM>> d2_Eij, d2_Bij, d2_h;
    Tensor<2, data_t, CH_SPACETIMEDIM> Eij_dot_d2_Eij, Bij_dot_d2_Bij;

    // compute d2_h
    FOR(i, j)
    {
        FOR(k)
        {
            // need to fill in reverse order in order to re-use them!

            FOR(l) // 2nd spatial ders
            {
                d2_h[i][j][k + 1][l + 1] =
                    d2.h[i][j][k][l] / vars.chi +
                    2. / (chi2 * vars.chi) * vars.h[i][j] * d1.chi[k] *
                        d1.chi[l] -
                    (d1.chi[k] * d1.h[i][j][l] + d1.chi[l] * d1.h[i][j][k] +
                     vars.h[i][j] * d2.chi[k][l]) /
                        chi2;
            }

            // mixed ders

            d2_h[i][j][0][k + 1] = -2. * vars.lapse * d1_Kij[i][j][k + 1] -
                                   2. * Kij[i][j] * d1.lapse[k];
            FOR(l)
            {
                d2_h[i][j][0][k + 1] +=
                    vars.shift[l] * d2_h[i][j][k + 1][l + 1] +
                    d1_h[i][j][l + 1] * d1.shift[l][k] +
                    d1_h[l][i][k + 1] * d1.shift[l][j] +
                    d1_h[l][j][k + 1] * d1.shift[l][i] +
                    vars.h[l][i] / vars.chi * d2.shift[l][j][k] +
                    vars.h[l][j] / vars.chi * d2.shift[l][i][k];
            }

            d2_h[i][j][k + 1][0] = d2_h[i][j][0][k + 1];
        }

        // 2nd time ders

        d2_h[i][j][0][0] =
            -2. * vars.lapse * d1_Kij[i][j][0] - 2. * Kij[i][j] * rhs.lapse;
        FOR(l)
        {
            d2_h[i][j][0][0] += vars.shift[l] * d2_h[i][j][0][l + 1] +
                                d1_h[i][j][l + 1] * rhs.shift[l] +
                                d1_h[l][i][0] * d1.shift[l][j] +
                                d1_h[l][j][0] * d1.shift[l][i] +
                                vars.h[l][i] / vars.chi * dtd1_shift[l][j] +
                                vars.h[l][j] / vars.chi * dtd1_shift[l][i];
        }
    }

    // compute d2_Eij, d2_Bij (raised or lowered depending on what we are
    // evolving)
    compute_d2_Eij_and_Bij(d2_Eij, d2_Bij, gq, rhs, d1_h);

    // now compute E_i^{~j} d_m d_n E_j^{~i}, B_i^{~j} d_m d_n B_j^{~i}
    if (m_params.use_last_index_raised)
    {
        FOR_ST(a, b)
        {
            Eij_dot_d2_Eij[a][b] = 0.;
            Bij_dot_d2_Bij[a][b] = 0.;
            FOR(i, j)
            {
                Eij_dot_d2_Eij[a][b] += Eij_LU[i][j] * d2_Eij[j][i][a][b];
                Bij_dot_d2_Bij[a][b] += Bij_LU[i][j] * d2_Bij[j][i][a][b];
            }
        }
    }
    else
    {
        // recast 'E_i^{~j} d_m d_n E_j^{~i}' to some expression that instead
        // uses d_m d_n E_{ij}
        // (see overleaf notes for this formula, it looks worse than it is)
        Tensor<2, data_t> Eij_UU =
            TensorAlgebra::compute_dot_product(metric_UU_spatial, Eij_LU, 0, 0);
        Tensor<2, data_t> Bij_UU =
            TensorAlgebra::compute_dot_product(metric_UU_spatial, Bij_LU, 0, 0);

        // computing some helping variables
        Tensor<2, Tensor<1, data_t, CH_SPACETIMEDIM>> d1_Eij_LL, d1_Bij_LL;
        Tensor<2, data_t> Eij_dot_Eij_UU, Bij_dot_Bij_UU;

        FOR(i, j)
        {
            d1_Eij_LL[i][j][0] = rhs.Eij[i][j];
            d1_Bij_LL[i][j][0] = rhs.Bij[i][j];

            Eij_dot_Eij_UU[i][j] = 0.;
            Bij_dot_Bij_UU[i][j] = 0.;

            FOR(k)
            {
                d1_Eij_LL[i][j][k + 1] = d1.Eij[i][j][k];
                d1_Bij_LL[i][j][k + 1] = d1.Bij[i][j][k];

                Eij_dot_Eij_UU[i][j] += Eij_UU[i][k] * Eij_LU[k][j];
                Bij_dot_Bij_UU[i][j] += Bij_UU[i][k] * Bij_LU[k][j];
            }
        }

        // now the real deal
        FOR_ST(a, b)
        {
            Eij_dot_d2_Eij[a][b] = 0.;
            Bij_dot_d2_Bij[a][b] = 0.;
            FOR(i, j)
            {
                Eij_dot_d2_Eij[a][b] += Eij_UU[i][j] * d2_Eij[i][j][a][b] -
                                        Eij_dot_Eij_UU[i][j] * d2_h[i][j][a][b];
                Bij_dot_d2_Bij[a][b] += Bij_UU[i][j] * d2_Bij[i][j][a][b] -
                                        Bij_dot_Bij_UU[i][j] * d2_h[i][j][a][b];

                FOR(k)
                {
                    Eij_dot_d2_Eij[a][b] +=
                        -Eij_UU[k][j] * d1_Eij_LU[j][i][b] * d1_h[k][i][a] -
                        Eij_UU[k][j] * d1_Eij_LU[j][i][a] * d1_h[k][i][b];
                    Bij_dot_d2_Bij[a][b] +=
                        -Bij_UU[k][j] * d1_Bij_LU[j][i][b] * d1_h[k][i][a] -
                        Bij_UU[k][j] * d1_Bij_LU[j][i][a] * d1_h[k][i][b];
                }
            }
        }
    }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Now make C, d1_C, d2_C
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // with all the above C, d1_C, d2_C are pretty straightforward
    C = 0.;
    FOR_ST(a)
    {
        d1_C[a] = 0.;
        FOR_ST(b)
        {
            d2_C[a][b] = 16. * (Eij_dot_d2_Eij[a][b] - Bij_dot_d2_Bij[a][b]);
        }
    }

    FOR(i, j)
    {
        C += 8. * (Eij_LU[i][j] * Eij_LU[j][i] - Bij_LU[i][j] * Bij_LU[j][i]);

        FOR_ST(a)
        {
            d1_C[a] += 16. * (Eij_LU[i][j] * d1_Eij_LU[j][i][a] -
                              Bij_LU[i][j] * d1_Bij_LU[j][i][a]);

            FOR_ST(b)
            {
                d2_C[a][b] += 16. * (d1_Eij_LU[i][j][b] * d1_Eij_LU[j][i][a] -
                                     d1_Bij_LU[i][j][b] * d1_Bij_LU[j][i][a]);
            }
        }
    }
}

template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t, class gauge_t>
void EBSystem::compute_Riemann(
    Tensor<4, data_t, CH_SPACETIMEDIM> &riemann_LLLU,
    Tensor<4, data_t, CH_SPACETIMEDIM> &riemann_LULU,
    GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq) const
{
    const auto &vars = gq.get_vars();
    const auto &g_UU = gq.get_metric_UU_ST();

    Tensor<2, data_t> Eij_LL;
    Tensor<2, data_t> Bij_LL;

    if (m_params.use_last_index_raised)
    {
        // evolved is LU, need to lower last index
        const auto &metric_spatial = gq.get_metric_spatial();
        Eij_LL = TensorAlgebra::compute_dot_product(vars.Eij, metric_spatial);
        Bij_LL = TensorAlgebra::compute_dot_product(vars.Bij, metric_spatial);
    }
    else
    {
        // evolved is already LL
        Eij_LL = vars.Eij;
        Bij_LL = vars.Bij;
    }

    // valid in vacuum!
    Tensor<4, data_t, CH_SPACETIMEDIM> riemann_LLLL =
        gq.compute_weyl_tensor_LLLL(Eij_LL, Bij_LL);

    FOR_ST(a, b, c, d, e)
    {
        riemann_LLLU[a][b][c][d] += riemann_LLLL[a][b][c][e] * g_UU[e][d];
    }

    FOR_ST(a, b, c, d, e)
    {
        riemann_LULU[a][b][c][d] += riemann_LLLU[a][e][c][d] * g_UU[e][b];
    }
}

// Adds in the RHS for the matter vars
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t, class gauge_t,
          template <typename> class rhs_vars_t>
void EBSystem::add_matter_rhs(
    rhs_vars_t<data_t> &total_rhs,
    GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq,
    const C2EFT<EBSystem>::params_t &pm) const
{
    const auto &vars = gq.get_vars();

    data_t tau = m_params.tau;
    data_t sigma = m_params.sigma;

    if (m_params.use_tau_radial_decay)
    {
        const Coordinates<data_t> &coords = gq.get_coordinates();
        tau = tau_from_radius(coords.get_radius());
    }

    if (m_params.use_sigma_radial_decay)
    {
        const Coordinates<data_t> &coords = gq.get_coordinates();
        sigma = sigma_from_radius(coords.get_radius());
    }

    if (m_params.use_tau_chi_decay)
    {
        tau = (m_params.tau_asymptotic - m_params.tau) /
                  (1.0 + exp(-(vars.chi / m_params.tau_decay_length - 1.0) /
                             m_params.tau_decay_width)) +
              m_params.tau;
    }

    if (m_params.use_sigma_chi_decay)
    {
        sigma = (m_params.sigma_asymptotic - m_params.sigma) /
                    (1.0 + exp(-(vars.chi / m_params.sigma_decay_length - 1.0) /
                               m_params.sigma_decay_width)) +
                m_params.sigma;
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
        // recall: Eaux is the physical E, Baux the physical B
        // (computed numerically from the metric on each cell)
        FOR(i, j)
        {
            total_rhs.Eij[i][j] =
                -1. / tau_rescaled * (vars.Eij[i][j] - vars.Eaux[i][j]);
            total_rhs.Bij[i][j] =
                -1. / tau_rescaled * (vars.Bij[i][j] - vars.Baux[i][j]);
            total_rhs.Eaux[i][j] = 0.;
            total_rhs.Baux[i][j] = 0.;
        }
    }
    else if (m_params.version == 2)
    {
        // recall: Eaux and Baux are the 1st time derivative of E and B

        const auto &d2 = gq.get_d2_vars();
        const auto &advec = gq.get_advection();

        // first compute the LL ones, but when evolving the raised E and B
        // these need to be raised
        Tensor<2, data_t> Eij = gq.get_weyl_electric_part();
        Tensor<2, data_t> Bij = gq.get_weyl_magnetic_part();

        if (m_params.use_last_index_raised)
        {
            // raising them here
            const auto &metric_UU_spatial = gq.get_metric_UU_spatial();
            Eij = TensorAlgebra::compute_dot_product(Eij, metric_UU_spatial);
            Bij = TensorAlgebra::compute_dot_product(Bij, metric_UU_spatial);
        }

        FOR(i, j)
        {
            data_t Eaux_with_advec = vars.Eaux[i][j];
            data_t Baux_with_advec = vars.Baux[i][j];

            if (m_params.advection_type == 1 || m_params.advection_type == 2)
            {
                Eaux_with_advec -= m_params.advection_coeff * advec.Eij[i][j];
                Baux_with_advec -= m_params.advection_coeff * advec.Bij[i][j];
            }

            total_rhs.Eij[i][j] = vars.Eaux[i][j];
            total_rhs.Bij[i][j] = vars.Baux[i][j];
            total_rhs.Eaux[i][j] =
                (-tau_rescaled * Eaux_with_advec + Eij[i][j] - vars.Eij[i][j]) /
                sigma_rescaled;
            total_rhs.Baux[i][j] =
                (-tau_rescaled * Baux_with_advec + Bij[i][j] - vars.Bij[i][j]) /
                sigma_rescaled;

            if (m_params.advection_type == 2)
            {
                total_rhs.Eaux[i][j] +=
                    2. * m_params.advection_coeff * advec.Eaux[i][j];
                total_rhs.Baux[i][j] +=
                    2. * m_params.advection_coeff * advec.Baux[i][j];

                FOR(k, l)
                {
                    total_rhs.Eaux[i][j] +=
                        -m_params.advection_coeff * m_params.advection_coeff *
                        vars.shift[k] * vars.shift[l] * d2.Eij[i][j][k][l];
                    total_rhs.Baux[i][j] +=
                        -m_params.advection_coeff * m_params.advection_coeff *
                        vars.shift[k] * vars.shift[l] * d2.Bij[i][j][k][l];
                }
            }
        }
    }
    else if (m_params.version == 3)
    {
        // recall: Eaux and Baux are the 1st time derivative of E and B
        const auto &d1 = gq.get_d1_vars();
        const auto &d2 = gq.get_d2_vars();
        const auto &advec = gq.get_advection();
        const auto &g_UU = gq.get_metric_UU_ST();
        const auto &Gamma_ST = gq.get_Gamma_ST();

        // first compute the LL ones, but when evolving the raised E and B
        // these need to be raised
        Tensor<2, data_t> Eij = gq.get_weyl_electric_part();
        Tensor<2, data_t> Bij = gq.get_weyl_magnetic_part();

        if (m_params.use_last_index_raised)
        {
            // raising them here
            const auto &metric_UU_spatial = gq.get_metric_UU_spatial();
            Eij = TensorAlgebra::compute_dot_product(Eij, metric_UU_spatial);
            Bij = TensorAlgebra::compute_dot_product(Bij, metric_UU_spatial);
        }

        data_t advfac = m_params.advection_type ? 1.0 : 0.0;
        FOR(i, j)
        {

            total_rhs.Eij[i][j] = vars.Eaux[i][j];
            total_rhs.Bij[i][j] = vars.Baux[i][j];

            total_rhs.Eaux[i][j] =
                -tau / vars.lapse *
                    (vars.Eaux[i][j] - advfac * advec.Eij[i][j]) +
                Eij[i][j] - vars.Eij[i][j] -
                sigma * Gamma_ST[0] * vars.Eaux[i][j];

            total_rhs.Baux[i][j] =
                -tau / vars.lapse *
                    (vars.Baux[i][j] - advfac * advec.Bij[i][j]) +
                Bij[i][j] - vars.Bij[i][j] -
                sigma * Gamma_ST[0] * vars.Baux[i][j];

            FOR(k)
            {
                total_rhs.Eaux[i][j] +=
                    sigma * (2. * g_UU[0][k + 1] * d1.Eaux[i][j][k] -
                             Gamma_ST[k + 1] * d1.Eij[i][j][k]);

                total_rhs.Baux[i][j] +=
                    sigma * (2. * g_UU[0][k + 1] * d1.Baux[i][j][k] -
                             Gamma_ST[k + 1] * d1.Bij[i][j][k]);
                FOR(l)
                {
                    total_rhs.Eaux[i][j] +=
                        sigma * g_UU[k + 1][l + 1] * d2.Eij[i][j][k][l];
                    total_rhs.Baux[i][j] +=
                        sigma * g_UU[k + 1][l + 1] * d2.Bij[i][j][k][l];
                }
            }
            total_rhs.Eaux[i][j] /= (-sigma * g_UU[0][0]);
            total_rhs.Baux[i][j] /= (-sigma * g_UU[0][0]);
        }
    }
    else if (m_params.version == 4)
    {
        // recall: Eaux and Baux are the 1st time derivative of E and B
        const auto &d1 = gq.get_d1_vars();
        const auto &d2 = gq.get_d2_vars();
        const auto &advec = gq.get_advection();
        const auto &g_UU = gq.get_metric_UU_ST();
        const auto &g_LL = gq.get_metric_ST();        
        const auto &CDn = gq.get_CD_n_UL_ST();  
        const auto &n = gq.get_normal_U_ST();
        const auto &Chris = gq.get_Chris_ULL_ST();        
        const auto &lieK = gq.get_lie_extrinsic_curvature();
        const auto &metric_spatial = gq.get_metric_spatial();         
        const auto &metric_UU_spatial = gq.get_metric_UU_spatial();
        const auto &levi_civita_spatial = gq.get_levi_civita_spatial();
        const auto &Kij = gq.get_extrinsic_curvature();
        const auto &Lie_acceleration = gq.get_LIE_acceleration_U_ST();        
        const auto &acceleration = gq.get_acceleration_U_ST();                        
        const auto &covD_K_tensor = gq.get_covd_extrinsic_curvature();
        const auto &chris_spatial = gq.get_chris_spatial();   ///this is ULL                             
        const auto &d1_chris_spatial_ULLL = gq.get_d1_chris_spatial_ULLL();
              
    	Tensor<1, data_t, CH_SPACETIMEDIM> acceleration_L;    	    	
        acceleration_L = TensorAlgebra::compute_dot_product(acceleration, g_LL);

        // first compute the LL ones, but when evolving the raised E and B
        // these need to be raised
        Tensor<2, data_t> Eij = gq.get_weyl_electric_part();
        Tensor<2, data_t> Bij = gq.get_weyl_magnetic_part();

        if (m_params.use_last_index_raised)
        {
            // raising them here
            const auto &metric_UU_spatial = gq.get_metric_UU_spatial();
            Eij = TensorAlgebra::compute_dot_product(Eij, metric_UU_spatial);
            Bij = TensorAlgebra::compute_dot_product(Bij, metric_UU_spatial);
        }
               
    	Tensor<2, data_t> Eij_LU, Bij_LU;    	    	
        Eij_LU = TensorAlgebra::compute_dot_product(vars.Eij, metric_UU_spatial);
        Bij_LU = TensorAlgebra::compute_dot_product(vars.Bij, metric_UU_spatial);
        
        Tensor<2, data_t> Eij_UU, Bij_UU ;
        Eij_UU = TensorAlgebra::compute_dot_product(metric_UU_spatial, Eij_LU, 0, 0);         
        Bij_UU = TensorAlgebra::compute_dot_product(metric_UU_spatial, Bij_LU, 0, 0);

        Tensor<3, data_t> D_Eij;        
        D_Eij = TensorAlgebra::covariant_derivative(d1.Eij, vars.Eij, chris_spatial);
        
        Tensor<3, data_t> D_Eij_LUL;  // this is as : D_Bij[i][j][k]  D_[k] Bij_[i]^[j]
        FOR(i, j, k)
        {
            FOR(l)
            {
        	D_Eij_LUL[i][j][k] += D_Eij[i][l][k]*metric_UU_spatial[l][j] ;
            }	 
        }        
        
        Tensor<3, data_t> D_Bij;        
        D_Bij = TensorAlgebra::covariant_derivative(d1.Bij, vars.Bij, chris_spatial);
        
        Tensor<3, data_t> D_Bij_LUL;  // this is as : D_Bij[i][j][k]  D_[k] Bij_[i]^[j]
        FOR(i, j, k)
        {
            FOR(l)
            {
        	D_Bij_LUL[i][j][k] += D_Bij[i][l][k]*metric_UU_spatial[l][j] ;
            }	 
        }          
        
        Tensor<4, data_t> DD_Eij; /// this is as : DD_Eij[i][j][l][m] D_m D_l Eij[i][j] 
        
        FOR(i, j , l , m)
        {
            DD_Eij[i][j][l][m] = d2.Eij[i][j][m][l] ;
            FOR(k)
            {
            DD_Eij[i][j][l][m] += - d1_chris_spatial_ULLL[k][i][l][m]*vars.Eij[k][j] - chris_spatial[k][i][l]*d1.Eij[k][j][m]
            		          - d1_chris_spatial_ULLL[k][l][j][m]*vars.Eij[i][k]  - chris_spatial[k][l][j]*d1.Eij[i][k][m] 
				  - chris_spatial[k][m][l]*D_Eij[i][j][k] - chris_spatial[k][m][i]*D_Eij[k][j][l] - chris_spatial[k][m][j]*D_Eij[i][k][l] ;            		                     
            }                         
        }
        
        Tensor<4, data_t> DD_Bij; /// this is as : DD_Bij[i][j][l][m] D_m D_l Bij[i][j] 
        
        FOR(i, j , l , m)
        {
            DD_Bij[i][j][l][m] = d2.Bij[i][j][m][l] ;
            FOR(k)
            {
            DD_Bij[i][j][l][m] += - d1_chris_spatial_ULLL[k][i][l][m]*vars.Bij[k][j] - chris_spatial[k][i][l]*d1.Bij[k][j][m]
            		          - d1_chris_spatial_ULLL[k][l][j][m]*vars.Bij[i][k]  - chris_spatial[k][l][j]*d1.Bij[i][k][m] 
				  - chris_spatial[k][m][l]*D_Bij[i][j][k] - chris_spatial[k][m][i]*D_Bij[k][j][l] - chris_spatial[k][m][j]*D_Bij[i][k][l] ;            		                     
            }
        }                 
        
    	Tensor<2, data_t> Kij_LU ;    	
        Kij_LU = TensorAlgebra::compute_dot_product(Kij, metric_UU_spatial);
        
        Tensor<2, data_t> Kij_UU ;
        Kij_UU = TensorAlgebra::compute_dot_product(metric_UU_spatial, Kij_LU, 0, 0);        
        
        Tensor<3, data_t> DKij_UUL ;
	FOR(i, j, k)
	{
	    FOR(l,m)
	    {
	        DKij_UUL[i][j][k] += covD_K_tensor[l][m][k]*metric_UU_spatial[i][l]*metric_UU_spatial[j][m] ;
	    }    
 	}
        
    	Tensor<2, data_t> Baux_LU, Eaux_LU ;    	    	
        Baux_LU = TensorAlgebra::compute_dot_product(vars.Baux, metric_UU_spatial);
        Eaux_LU = TensorAlgebra::compute_dot_product(vars.Eaux, metric_UU_spatial);                              
        
        FOR(i, j) // No contractions
        {

            total_rhs.Eij[i][j] = -vars.lapse*vars.Eaux[i][j];
            total_rhs.Bij[i][j] = -vars.lapse*vars.Baux[i][j];            
            FOR(k) // One spatial contraction
            {
            	total_rhs.Eij[i][j] += vars.shift[k]*d1.Eij[i][j][k] -vars.lapse*(vars.Eij[i][k]*CDn[k+1][j+1] + vars.Eij[k][j]*CDn[k+1][i+1]) ;
            	total_rhs.Bij[i][j] += vars.shift[k]*d1.Bij[i][j][k] -vars.lapse*(vars.Bij[i][k]*CDn[k+1][j+1] + vars.Bij[k][j]*CDn[k+1][i+1]) ;            	
            	
            	FOR_ST(a) // One ST contraction
            	{
             	    total_rhs.Eij[i][j] += vars.lapse*(n[a]*(Chris[k+1][i+1][a]*vars.Eij[k][j]  + Chris[k+1][j+1][a]*vars.Eij[i][k])) ;
             	    total_rhs.Bij[i][j] += vars.lapse*(n[a]*(Chris[k+1][i+1][a]*vars.Bij[k][j]  + Chris[k+1][j+1][a]*vars.Bij[i][k])) ;              	        
            	}
            }
            
            ////////////////////////WE START WITH THE SECOND ORDER EQUATION////////
            
            FOR(k) // One spatial contraction
            {
            	total_rhs.Eaux[i][j] += vars.shift[k]*d1.Eaux[i][j][k]/vars.lapse - (vars.Eaux[i][k]*CDn[k+1][j+1] + vars.Eaux[k][j]*CDn[k+1][i+1]) ;
            	total_rhs.Baux[i][j] += vars.shift[k]*d1.Baux[i][j][k]/vars.lapse - (vars.Baux[i][k]*CDn[k+1][j+1] + vars.Baux[k][j]*CDn[k+1][i+1]) ;
            	            	
            	FOR_ST(a) // One ST contraction
            	{
             	    total_rhs.Eaux[i][j] += (n[a]*(Chris[k+1][i+1][a]*vars.Eaux[k][j]  + Chris[k+1][j+1][a]*vars.Eaux[i][k])) ;
             	    total_rhs.Baux[i][j] += (n[a]*(Chris[k+1][i+1][a]*vars.Baux[k][j]  + Chris[k+1][j+1][a]*vars.Baux[i][k])) ;             	        
            	}
            	
	    }

	    if (m_params.epsilon > 0.0)
            {
	    		    
	    total_rhs.Eaux[i][j] += vars.K*vars.Eaux[i][j] ; // no contractions
	    total_rhs.Baux[i][j] += vars.K*vars.Baux[i][j] ; // no contractions
	    	    
	    FOR_ST(a) // One sT contraction
	    {
	        total_rhs.Eaux[i][j] -= -4.*acceleration[a]*acceleration_L[a]*vars.Eij[i][j];
	        total_rhs.Baux[i][j] -= -4.*acceleration[a]*acceleration_L[a]*vars.Bij[i][j];	        
	    }
	    	    
	    FOR(k) // One spatial contraction
	    {
	        total_rhs.Eaux[i][j] -= - Eij_LU[i][k]*lieK[j][k] - Eij_LU[j][k]*lieK[i][k] 
					+ 2.*Kij_LU[i][k]*vars.Eaux[j][k] + 2.*Kij_LU[j][k]*vars.Eaux[i][k] 
	        			+ 3.*acceleration[k+1]*acceleration_L[i+1]*vars.Eij[j][k] 
	        			+ 3.*acceleration[k+1]*acceleration_L[j+1]*vars.Eij[i][k]
	        			+ vars.K*Eij_LU[i][k]*Kij[j][k] + vars.K*Eij_LU[j][k]*Kij[i][k]	        			
	        			- 6.*(vars.Eij[k][i]*Eij_LU[j][k] - vars.Bij[k][i]*Bij_LU[j][k]);
	        			
	        total_rhs.Baux[i][j] -= - Bij_LU[i][k]*lieK[j][k] - Bij_LU[j][k]*lieK[i][k] 
					+ 2.*Kij_LU[i][k]*vars.Baux[j][k] + 2.*Kij_LU[j][k]*vars.Baux[i][k] 
	        			+ 3.*acceleration[k+1]*acceleration_L[i+1]*vars.Bij[j][k] 
	        			+ 3.*acceleration[k+1]*acceleration_L[j+1]*vars.Bij[i][k]
	        			+ vars.K*Bij_LU[i][k]*Kij[j][k] + vars.K*Bij_LU[j][k]*Kij[i][k]	        			
	        			- 6.*(vars.Eij[k][i]*Bij_LU[j][k] + vars.Eij[k][j]*Bij_LU[i][k]);	        			

	        			
	        FOR(l) // two spatial contractions
	        {
	            total_rhs.Eaux[i][j] -= -Bij_LU[i][l]*levi_civita_spatial[j][l][k]*Lie_acceleration[k+1] -Bij_LU[j][l]*levi_civita_spatial[i][l][k]*Lie_acceleration[k+1]
	                                    + DD_Eij[i][j][l][k]*metric_UU_spatial[l][k]
	            			    -2.*Baux_LU[i][l]*levi_civita_spatial[j][l][k]*acceleration[k+1] -2.*Baux_LU[j][l]*levi_civita_spatial[i][l][k]*acceleration[k+1]
	            			    -2.*acceleration[k+1]*acceleration[l+1]*vars.Eij[l][k]*metric_spatial[i][j]
	            			    +acceleration[k+1]*0.5*Bij_LU[i][l]*levi_civita_spatial[j][k][l]*vars.K
	            			    +acceleration[k+1]*0.5*Bij_LU[j][l]*levi_civita_spatial[i][k][l]*vars.K	            			    
					    +4.*Kij[l][k]*Kij_UU[l][k]*vars.Eij[i][j] 
					    -6.*Kij[l][k]*Eij_LU[i][l]*Kij_LU[j][k] -6.*Kij[l][k]*Eij_LU[j][l]*Kij_LU[i][k]
					    -2.*Eij_UU[l][k]*Kij[i][l]*Kij[j][k]
					    +2.*metric_spatial[i][j]*(vars.Eij[l][k]*Eij_UU[l][k] - vars.Bij[l][k]*Bij_UU[l][k]);
					    
	            total_rhs.Baux[i][j] -= -Eij_LU[i][l]*levi_civita_spatial[j][l][k]*Lie_acceleration[k+1] -Eij_LU[j][l]*levi_civita_spatial[i][l][k]*Lie_acceleration[k+1]
	                                    + DD_Bij[i][j][l][k]*metric_UU_spatial[l][k]
	            			    +2.*Eaux_LU[i][l]*levi_civita_spatial[j][l][k]*acceleration[k+1] +2.*Eaux_LU[j][l]*levi_civita_spatial[i][l][k]*acceleration[k+1]
	            			    -2.*acceleration[k+1]*acceleration[l+1]*vars.Bij[l][k]*metric_spatial[i][j]
	            			    -acceleration[k+1]*0.5*Eij_LU[i][l]*levi_civita_spatial[j][k][l]*vars.K
	            			    -acceleration[k+1]*0.5*Eij_LU[j][l]*levi_civita_spatial[i][k][l]*vars.K	            			    
					    +4.*Kij[l][k]*Kij_UU[l][k]*vars.Bij[i][j] 
					    -6.*Kij[l][k]*Bij_LU[i][l]*Kij_LU[j][k] -6.*Kij[l][k]*Bij_LU[j][l]*Kij_LU[i][k]
					    -2.*Bij_UU[l][k]*Kij[i][l]*Kij[j][k]
					    +4.*metric_spatial[i][j]*vars.Eij[l][k]*Bij_UU[l][k] ;					    
					    
		   FOR(m) // three spatial contractions
		   {
		       total_rhs.Eaux[i][j] -= Bij_LU[i][l]*levi_civita_spatial[j][l][k]*DKij_UUL[m][k][m]
		       			     + Bij_LU[j][l]*levi_civita_spatial[i][l][k]*DKij_UUL[m][k][m]
		       			     + 2.*D_Bij_LUL[i][k][l]*levi_civita_spatial[j][k][m]*Kij_UU[l][m]
		       			     + 2.*D_Bij_LUL[j][k][l]*levi_civita_spatial[i][k][m]*Kij_UU[l][m] 
					     - 2.*acceleration[k+1]*Bij_UU[l][m]*Kij[l][i]*levi_civita_spatial[j][k][m]
		       			     - 2.*acceleration[k+1]*Bij_UU[l][m]*Kij[l][j]*levi_civita_spatial[i][k][m]
		       			     - acceleration[k+1]*0.5*Bij_LU[i][l]*(levi_civita_spatial[j][l][m]*Kij_LU[k][m] + 2.*levi_civita_spatial[j][k][m]*Kij_LU[l][m])
		       			     - acceleration[k+1]*0.5*Bij_LU[j][l]*(levi_civita_spatial[i][l][m]*Kij_LU[k][m] + 2.*levi_civita_spatial[i][k][m]*Kij_LU[l][m])
		       			     + 2.*Eij_UU[l][k]*Kij_LU[l][m]*Kij[k][m]*metric_spatial[i][j];
		       			     
		       total_rhs.Baux[i][j] -= -Eij_LU[i][l]*levi_civita_spatial[j][l][k]*DKij_UUL[m][k][m]
		       			     + -Eij_LU[j][l]*levi_civita_spatial[i][l][k]*DKij_UUL[m][k][m]
		       			     - 2.*D_Eij_LUL[i][k][l]*levi_civita_spatial[j][k][m]*Kij_UU[l][m]
		       			     - 2.*D_Eij_LUL[j][k][l]*levi_civita_spatial[i][k][m]*Kij_UU[l][m] 
					     + 2.*acceleration[k+1]*Eij_UU[l][m]*Kij[l][i]*levi_civita_spatial[j][k][m]
		       			     + 2.*acceleration[k+1]*Eij_UU[l][m]*Kij[l][j]*levi_civita_spatial[i][k][m]
		       			     + acceleration[k+1]*0.5*Eij_LU[i][l]*(levi_civita_spatial[j][l][m]*Kij_LU[k][m] + 2.*levi_civita_spatial[j][k][m]*Kij_LU[l][m])
		       			     + acceleration[k+1]*0.5*Eij_LU[j][l]*(levi_civita_spatial[i][l][m]*Kij_LU[k][m] + 2.*levi_civita_spatial[i][k][m]*Kij_LU[l][m])
		       			     + 2.*Bij_UU[l][k]*Kij_LU[l][m]*Kij[k][m]*metric_spatial[i][j];		       			     
		   }			    			    
	        }				
	    }
	    }	    
	    total_rhs.Eaux[i][j] *= vars.lapse;
	    total_rhs.Baux[i][j] *= vars.lapse;	 
	    
	    
	    // NOW THE FIXING PART
	    
	   // total_rhs.Eaux[i][j] += -vars.lapse*tau_rescaled/sigma_rescaled*vars.Eaux[i][j] + vars.lapse/sigma_rescaled*(vars.Eij[i][j] - Eij[i][j] );
	   //total_rhs.Baux[i][j] += -vars.lapse*tau_rescaled/sigma_rescaled*vars.Baux[i][j] + vars.lapse/sigma_rescaled*(vars.Bij[i][j] - Bij[i][j] );
	    
	    
	    total_rhs.Eaux[i][j] += -tau_rescaled/sigma_rescaled*(-total_rhs.Eij[i][j]) + sigma_rescaled*(vars.Eij[i][j] - Eij[i][j] );
	    total_rhs.Baux[i][j] += -tau_rescaled/sigma_rescaled*(-total_rhs.Bij[i][j]) + sigma_rescaled*(vars.Bij[i][j] - Bij[i][j] );	       
        }
    }
    else
    {
        MayDay::Error("Version not implemented");
    }
}

// Adds in the RHS for the matter vars
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t, class gauge_t,
          template <typename> class rhs_vars_t>
void EBSystem::compute_d2_Eij_and_Bij(
    Tensor<2, Tensor<2, data_t, CH_SPACETIMEDIM>> &d2_Eij,
    Tensor<2, Tensor<2, data_t, CH_SPACETIMEDIM>> &d2_Bij,
    GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq,
    rhs_vars_t<data_t> &rhs,
    Tensor<2, Tensor<1, data_t, CH_SPACETIMEDIM>> &d1_h) const
{
    if (m_params.version == 1)
    {
        const auto &vars = gq.get_vars();
        const auto &d1 = gq.get_d1_vars();
        const auto &d2 = gq.get_d2_vars();
        const auto &advec = gq.get_advection();
        const auto &metric_UU_spatial = gq.get_metric_UU_spatial();
        const auto &Kij = gq.get_extrinsic_curvature();

        Tensor<2, Tensor<1, data_t>> d1_Eaux_LL, d1_Baux_LL;
        Tensor<2, data_t> Eaux_LL, Baux_LL, advec_Eaux_LL, advec_Baux_LL;

        if (m_params.use_last_index_raised)
        {
            // this is a pain
            // to compute the time derivative of LL via
            // GeometricQuantities::compute_dt_weyl_electric_part, so all
            // vars.E, d1.E, advec.E (same for B) that have one index raised
            // need to be lowered to have the right input (and raised
            // afterwards)
            const auto &metric_spatial = gq.get_metric_spatial();

            Tensor<2, data_t> advec_h;
            data_t chi2 = vars.chi * vars.chi;
            FOR(i, j)
            {
                advec_h[i][j] =
                    advec.h[i][j] / vars.chi - advec.chi * vars.h[i][j] / chi2;
            }

            FOR(i, j)
            {
                advec_Eaux_LL[i][j] = 0.;
                advec_Baux_LL[i][j] = 0.;

                FOR(k)
                {
                    advec_Eaux_LL[i][j] +=
                        advec.Eaux[i][k] * metric_spatial[k][j] +
                        vars.Eaux[i][k] * advec_h[k][j];
                    advec_Baux_LL[i][j] +=
                        advec.Baux[i][k] * metric_spatial[k][j] +
                        vars.Baux[i][k] * advec_h[k][j];

                    d1_Eaux_LL[i][j][k] = 0.;
                    d1_Baux_LL[i][j][k] = 0.;
                    FOR(l)
                    {
                        d1_Eaux_LL[i][j][k] +=
                            d1.Eaux[i][l][k] * metric_spatial[l][j] +
                            vars.Eaux[i][l] * d1_h[l][j][k + 1];
                        d1_Baux_LL[i][j][k] +=
                            d1.Baux[i][l][k] * metric_spatial[l][j] +
                            vars.Baux[i][l] * d1_h[l][j][k + 1];
                    }
                }
            }
            Eaux_LL =
                TensorAlgebra::compute_dot_product(vars.Eaux, metric_spatial);
            Baux_LL =
                TensorAlgebra::compute_dot_product(vars.Baux, metric_spatial);
        }
        else
        {
            // all good in this case, we have what we need
            d1_Eaux_LL = d1.Eaux;
            d1_Baux_LL = d1.Baux;
            Eaux_LL = vars.Eaux;
            Baux_LL = vars.Baux;
            advec_Eaux_LL = advec.Eaux;
            advec_Baux_LL = advec.Baux;
        }

        // compute time derivatives
        Tensor<2, data_t> dt_Eij = gq.compute_dt_weyl_electric_part(
            d1_Eaux_LL, d1_Baux_LL, Eaux_LL, Baux_LL, advec_Eaux_LL,
            advec_Baux_LL);
        Tensor<2, data_t> dt_Bij = gq.compute_dt_weyl_magnetic_part(
            d1_Eaux_LL, d1_Baux_LL, Eaux_LL, Baux_LL, advec_Eaux_LL,
            advec_Baux_LL);

        if (m_params.use_last_index_raised)
        {
            // now need to raise them back to get the time derivative LU
            // this needs the rhs of 'h', so d1_h[i][j][0]
            Tensor<2, data_t> dt_Eij_LU, dt_Bij_LU;

            Tensor<2, data_t> Eij_LU = TensorAlgebra::compute_dot_product(
                gq.get_weyl_electric_part(), metric_UU_spatial);
            Tensor<2, data_t> Bij_LU = TensorAlgebra::compute_dot_product(
                gq.get_weyl_magnetic_part(), metric_UU_spatial);

            FOR(i, j)
            {
                dt_Eij_LU[i][j] = 0.;
                dt_Bij_LU[i][j] = 0.;

                FOR(k)
                {
                    dt_Eij_LU[i][j] += dt_Eij[i][k] * metric_UU_spatial[k][j];
                    dt_Bij_LU[i][j] += dt_Bij[i][k] * metric_UU_spatial[k][j];

                    FOR(l)
                    {
                        dt_Eij_LU[i][j] += -Eij_LU[i][k] * d1_h[l][k][0] *
                                           metric_UU_spatial[j][l];
                        dt_Bij_LU[i][j] += -Bij_LU[i][k] * d1_h[l][k][0] *
                                           metric_UU_spatial[j][l];
                    }
                }
            }
            // replace old
            dt_Eij = dt_Eij_LU;
            dt_Bij = dt_Bij_LU;
        }

        data_t tau = m_params.tau;

        if (m_params.use_tau_radial_decay)
        {
            const Coordinates<data_t> &coords = gq.get_coordinates();
            tau = tau_from_radius(coords.get_radius());
        }

        if (m_params.use_tau_chi_decay)
        {
            tau = (m_params.tau_asymptotic - m_params.tau) /
                      (1.0 + exp(-(vars.chi / m_params.tau_decay_length - 1.0) /
                                 m_params.tau_decay_width)) +
                  m_params.tau;
        }

        if (m_params.rescale_tau_by_lapse)
            tau /= vars.lapse;

        FOR(i, j)
        {
            FOR(k)
            {
                FOR(l) // 2nd spatial ders
                {
                    d2_Eij[i][j][k + 1][l + 1] = d2.Eij[i][j][k][l];
                    d2_Bij[i][j][k + 1][l + 1] = d2.Bij[i][j][k][l];
                }

                // mixed ders

                d2_Eij[i][j][0][k + 1] =
                    -1. / tau * (d1.Eij[i][j][k] - d1.Eaux[i][j][k]);
                d2_Eij[i][j][k + 1][0] = d2_Eij[i][j][0][k + 1];

                d2_Bij[i][j][0][k + 1] =
                    -1. / tau * (d1.Bij[i][j][k] - d1.Baux[i][j][k]);
                d2_Bij[i][j][k + 1][0] = d2_Bij[i][j][0][k + 1];
            }

            // 2nd time ders

            d2_Eij[i][j][0][0] = -1. / tau * (rhs.Eij[i][j] - dt_Eij[i][j]);
            d2_Bij[i][j][0][0] = -1. / tau * (rhs.Bij[i][j] - dt_Bij[i][j]);
        }
    }
    else if (m_params.version == 2 || m_params.version == 3 )
    {
        const auto &d1 = gq.get_d1_vars();
        const auto &d2 = gq.get_d2_vars();

        FOR(i, j)
        {
            FOR(k)
            {
                FOR(l) // 2nd spatial ders
                {
                    d2_Eij[i][j][k + 1][l + 1] = d2.Eij[i][j][k][l];
                    d2_Bij[i][j][k + 1][l + 1] = d2.Bij[i][j][k][l];
                }

                // mixed ders

                d2_Eij[i][j][0][k + 1] = d1.Eaux[i][j][k];
                d2_Eij[i][j][k + 1][0] = d2_Eij[i][j][0][k + 1];

                d2_Bij[i][j][0][k + 1] = d1.Baux[i][j][k];
                d2_Bij[i][j][k + 1][0] = d2_Bij[i][j][0][k + 1];
            }

            // 2nd time ders

            d2_Eij[i][j][0][0] = rhs.Eaux[i][j];
            d2_Bij[i][j][0][0] = rhs.Baux[i][j];
        }
    }
    else if (m_params.version == 4 )
    {
        const auto &vars = gq.get_vars();    
        const auto &d1 = gq.get_d1_vars();
        const auto &d2 = gq.get_d2_vars();
        const auto &n = gq.get_normal_U_ST();
        const auto &d1n = gq.get_d1_n_UL_ST();        
        const auto &CDn = gq.get_CD_n_UL_ST();
        const auto &d1CDn = gq.get_d1CD_n_ULL_ST();
        const auto &Chris = gq.get_Chris_ULL_ST();
        const auto &d1Chris = gq.get_d1_Chris_ULLL_ST();        
        
        FOR(i, j)
        {
            FOR(k)
            {
                FOR(l) // 2nd spatial ders
                {
                    d2_Eij[i][j][k + 1][l + 1] = d2.Eij[i][j][k][l];
                    d2_Bij[i][j][k + 1][l + 1] = d2.Bij[i][j][k][l];
                }

                // mixed ders

                d2_Eij[i][j][0][k + 1] = -vars.lapse*d1.Eaux[i][j][k] - vars.Eaux[i][j]*d1.lapse[k] ; 
                d2_Bij[i][j][0][k + 1] = -vars.lapse*d1.Baux[i][j][k] - vars.Baux[i][j]*d1.lapse[k] ;
                                
                FOR(l)
                {
                    d2_Eij[i][j][0][k + 1] += d1.shift[l][k]*d1.Eij[i][j][l] 
                    			      + vars.shift[l]*d2.Eij[i][j][l][k]
                    			      + d1.lapse[k]*(-vars.Eij[i][l]*CDn[l+1][j+1] -vars.Eij[l][j]*CDn[l+1][i+1]) 
                    			      - vars.lapse*(d1.Eij[i][l][k]*CDn[l+1][j+1] + vars.Eij[i][l]*d1CDn[l+1][j+1][k+1] + d1.Eij[l][j][k]*CDn[l+1][i+1] + vars.Eij[l][j]*d1CDn[l+1][i+1][k+1]);
                    			      
                    d2_Bij[i][j][0][k + 1] += d1.shift[l][k]*d1.Bij[i][j][l] 
                    			      + vars.shift[l]*d2.Bij[i][j][l][k]
                    			      + d1.lapse[k]*(-vars.Bij[i][l]*CDn[l+1][j+1] -vars.Bij[l][j]*CDn[l+1][i+1]) 
                    			      - vars.lapse*(d1.Bij[i][l][k]*CDn[l+1][j+1] + vars.Bij[i][l]*d1CDn[l+1][j+1][k+1] + d1.Bij[l][j][k]*CDn[l+1][i+1] + vars.Bij[l][j]*d1CDn[l+1][i+1][k+1]);	      
                    			                      
                    FOR_ST(a)
                    {
                        d2_Eij[i][j][0][k + 1] += d1.lapse[k]*n[a]*(Chris[l+1][i+1][a]*vars.Eij[l][j] + Chris[l+1][j+1][a]*vars.Eij[i][l])
                    			          + vars.lapse*d1n[a][k+1]*(Chris[l+1][i+1][a]*vars.Eij[l][j] + Chris[l+1][j+1][a]*vars.Eij[i][l])
                    			          + vars.lapse*n[a]*(d1Chris[l+1][i+1][a][k+1]*vars.Eij[l][j] + Chris[l+1][i+1][a]*d1.Eij[l][j][k] 
                                                  + d1Chris[l+1][j+1][a][k+1]*vars.Eij[i][l] + Chris[l+1][j+1][a]*d1.Eij[i][l][k]) ;
                                              
                        d2_Bij[i][j][0][k + 1] += d1.lapse[k]*n[a]*(Chris[l+1][i+1][a]*vars.Bij[l][j] + Chris[l+1][j+1][a]*vars.Bij[i][l])
                    			          + vars.lapse*d1n[a][k+1]*(Chris[l+1][i+1][a]*vars.Bij[l][j] + Chris[l+1][j+1][a]*vars.Bij[i][l])
                    	   		          + vars.lapse*n[a]*(d1Chris[l+1][i+1][a][k+1]*vars.Bij[l][j] + Chris[l+1][i+1][a]*d1.Bij[l][j][k] 
                                                  + d1Chris[l+1][j+1][a][k+1]*vars.Bij[i][l] + Chris[l+1][j+1][a]*d1.Bij[i][l][k]) ;                                              
                    }                
                } 
                d2_Eij[i][j][k + 1][0] = d2_Eij[i][j][0][k + 1];
                d2_Bij[i][j][k + 1][0] = d2_Bij[i][j][0][k + 1];
            }

        }
        

            // 2nd time ders
        FOR(i, j)
        {            
                d2_Eij[i][j][0][0] = -vars.lapse*rhs.Eaux[i][j] - vars.Eaux[i][j]*rhs.lapse ;
                d2_Bij[i][j][0][0] = -vars.lapse*rhs.Baux[i][j] - vars.Baux[i][j]*rhs.lapse ;                
                
                FOR(l)
                {
                    d2_Eij[i][j][0][0] += rhs.shift[l]*d1.Eij[i][j][l] 
                    			  + vars.shift[l]*d2_Eij[i][j][l+1][0]
                    			  + rhs.lapse*(-vars.Eij[i][l]*CDn[l+1][j+1] -vars.Eij[l][j]*CDn[l+1][i+1]) 
                    			  - vars.lapse*(rhs.Eij[i][l]*CDn[l+1][j+1] + vars.Eij[i][l]*d1CDn[l+1][j+1][0] + rhs.Eij[l][j]*CDn[l+1][i+1] + vars.Eij[l][j]*d1CDn[l+1][i+1][0]);
                    			      
                    d2_Bij[i][j][0][0] += rhs.shift[l]*d1.Bij[i][j][l] 
                    			  + vars.shift[l]*d2_Bij[i][j][l+1][0]
                    			  + rhs.lapse*(-vars.Bij[i][l]*CDn[l+1][j+1] -vars.Bij[l][j]*CDn[l+1][i+1]) 
                    			  - vars.lapse*(rhs.Bij[i][l]*CDn[l+1][j+1] + vars.Bij[i][l]*d1CDn[l+1][j+1][0] + rhs.Bij[l][j]*CDn[l+1][i+1] + vars.Bij[l][j]*d1CDn[l+1][i+1][0]);
                
                    FOR_ST(a)
                    {
                        d2_Eij[i][j][0][0] += rhs.lapse*n[a]*(Chris[l+1][i+1][a]*vars.Eij[l][j] + Chris[l+1][j+1][a]*vars.Eij[i][l])
                    			      + vars.lapse*d1n[a][0]*(Chris[l+1][i+1][a]*vars.Eij[l][j] + Chris[l+1][j+1][a]*vars.Eij[i][l])
                    			      + vars.lapse*n[a]*(d1Chris[l+1][i+1][a][0]*vars.Eij[l][j] + Chris[l+1][i+1][a]*rhs.Eij[l][j] 
                                              + d1Chris[l+1][j+1][a][0]*vars.Eij[i][l] + Chris[l+1][j+1][a]*rhs.Eij[i][l]) ;
                                              
                        d2_Bij[i][j][0][0] += rhs.lapse*n[a]*(Chris[l+1][i+1][a]*vars.Bij[l][j] + Chris[l+1][j+1][a]*vars.Bij[i][l])
                    			      + vars.lapse*d1n[a][0]*(Chris[l+1][i+1][a]*vars.Bij[l][j] + Chris[l+1][j+1][a]*vars.Bij[i][l])
                    			      + vars.lapse*n[a]*(d1Chris[l+1][i+1][a][0]*vars.Bij[l][j] + Chris[l+1][i+1][a]*rhs.Bij[l][j] 
                                              + d1Chris[l+1][j+1][a][0]*vars.Bij[i][l] + Chris[l+1][j+1][a]*rhs.Bij[i][l]) ;
                    }
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
void EBSystem::add_diffusion_terms(
    rhs_vars_t<data_t> &rhs, //!< Reference to the variables into which the
                             //! output right hand side is written
    GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq,
    data_t diffCoeffSafe) const
{
    const auto &vars = gq.get_vars();
    const auto &d2 = gq.get_d2_vars();
    const auto &h_UU = gq.get_h_UU();

    Tensor<2, data_t> space_laplace_Eij = {0.};
    Tensor<2, data_t> space_laplace_Bij = {0.};

    FOR(i, j, k)
    {
        space_laplace_Eij[i][j] += d2.Eij[i][j][k][k];
        space_laplace_Bij[i][j] += d2.Bij[i][j][k][k];
    }

    data_t tr_space_laplace_Eij = 0.;
    data_t tr_space_laplace_Bij = 0.;
    if (m_params.use_last_index_raised)
    {
        // this is all very non-covariant, but for consistency I replaced the
        // metric by a delta function, since one of the indexes is already
        // raised
        FOR(i)
        {
            tr_space_laplace_Eij += space_laplace_Eij[i][i];
            tr_space_laplace_Bij += space_laplace_Bij[i][i];
        }

        FOR(i, j)
        {
            rhs.Eij[i][j] +=
                diffCoeffSafe * (space_laplace_Eij[i][j] -
                                 tr_space_laplace_Eij *
                                     TensorAlgebra::delta(i, j) / GR_SPACEDIM);
            rhs.Bij[i][j] +=
                diffCoeffSafe * (space_laplace_Bij[i][j] -
                                 tr_space_laplace_Bij *
                                     TensorAlgebra::delta(i, j) / GR_SPACEDIM);
        }
    }
    else
    {
        FOR(i, j)
        {
            tr_space_laplace_Eij += h_UU[i][j] * space_laplace_Eij[i][j];
            tr_space_laplace_Bij += h_UU[i][j] * space_laplace_Bij[i][j];
        }

        FOR(i, j)
        {
            rhs.Eij[i][j] += diffCoeffSafe * (space_laplace_Eij[i][j] -
                                              tr_space_laplace_Eij *
                                                  vars.h[i][j] / GR_SPACEDIM);
            rhs.Bij[i][j] += diffCoeffSafe * (space_laplace_Bij[i][j] -
                                              tr_space_laplace_Bij *
                                                  vars.h[i][j] / GR_SPACEDIM);
        }
    }
}

#endif /* EBSYSTEM_IMPL_HPP_ */

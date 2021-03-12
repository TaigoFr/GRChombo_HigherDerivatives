/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SYSTEMEB_HPP_)
#error "This file should only be included through SystemEB.hpp"
#endif

#ifndef SYSTEMEB_IMPL_HPP_
#define SYSTEMEB_IMPL_HPP_
#include "DimensionDefinitions.hpp"

// Calculate the stress energy tensor elements
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
void SystemEB::compute_C(
    data_t &C, Tensor<1, data_t, CH_SPACEDIM + 1> &d1_C,
    Tensor<2, data_t, CH_SPACEDIM + 1> &d2_C,
    GeometricQuantities<data_t, vars_t, diff2_vars_t> &gq) const
{
    const auto &vars = gq.get_vars();
    const auto &d1 = gq.get_d1_vars();
    const auto &d2 = gq.get_d2_vars();
    const auto &advec = gq.get_advection();
    const auto &metric_UU_spatial = gq.get_metric_UU_spatial();
    const auto &Kij = gq.get_extrinsic_curvature();

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // C is easy
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    C = 0.;
    FOR(i, j, k, l)
    {
        C +=
            8. * metric_UU_spatial[i][k] * metric_UU_spatial[j][l] *
            (vars.Eij[i][j] * vars.Eij[k][l] - vars.Bij[i][j] * vars.Bij[k][l]);
    }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // d1_C is ok
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // indices are only spatial, derivatives are 4D
    Tensor<2, Tensor<1, data_t, CH_SPACEDIM + 1>> d1_Eij;
    Tensor<2, Tensor<1, data_t, CH_SPACEDIM + 1>> d1_Bij;
    Tensor<2, Tensor<1, data_t, CH_SPACEDIM + 1>> d1_h;
    Tensor<2, Tensor<1, data_t, CH_SPACEDIM + 1>> d1_Kij; // needed for d2_h

    // 2nd mixed derivatives of shift needed for d2_h[i][j][0][0]
    Tensor<2, data_t> dtd1_shift; // mixed derivative
    auto ccz4_params = gq.get_formulation_params();

    vars_t<data_t> rhs;
    add_matter_rhs(rhs, gq);
    gq.compute_rhs_equations(rhs);

    data_t chi2 = vars.chi * vars.chi;

    FOR(i, j)
    {
        d1_Eij[i][j][0] = rhs.Eij[i][j];
        d1_Bij[i][j][0] = rhs.Bij[i][j];
        d1_h[i][j][0] = rhs.h[i][j] / vars.chi - rhs.chi * vars.h[i][j] / chi2;
        d1_Kij[i][j][0] =
            -Kij[i][j] * rhs.chi / vars.chi +
            (rhs.A[i][j] +
             (vars.h[i][j] * rhs.K + vars.K * rhs.h[i][j]) / GR_SPACEDIM) /
                vars.chi;

        dtd1_shift[i][j] = ccz4_params.shift_Gamma_coeff * d1.B[i][j];

        FOR(k)
        {
            d1_Eij[i][j][k + 1] = d1.Eij[i][j][k];
            d1_Bij[i][j][k + 1] = d1.Bij[i][j][k];
            d1_h[i][j][k + 1] =
                d1.h[i][j][k] / vars.chi - d1.chi[k] * vars.h[i][j] / chi2;
            d1_Kij[i][j][k + 1] = -Kij[i][j] * d1.chi[k] / vars.chi +
                                  (d1.A[i][j][k] + (vars.h[i][j] * d1.K[k] +
                                                    vars.K * d1.h[i][j][k]) /
                                                       GR_SPACEDIM) /
                                      vars.chi;

            dtd1_shift[i][j] += ccz4_params.shift_advec_coeff *
                                (vars.shift[k] * d2.shift[i][j][k] +
                                 d1.shift[i][k] * d1.shift[k][j]);
        }
    }

    Tensor<2, data_t> Eij_LU =
        TensorAlgebra::compute_dot_product(vars.Eij, metric_UU_spatial);
    Tensor<2, data_t> Eij_UU =
        TensorAlgebra::compute_dot_product(metric_UU_spatial, Eij_LU, 0, 0);
    Tensor<2, data_t> Eij_dot_Eij_UU =
        TensorAlgebra::compute_dot_product(Eij_UU, Eij_LU, 0, 0);

    Tensor<2, data_t> Bij_LU =
        TensorAlgebra::compute_dot_product(vars.Bij, metric_UU_spatial);
    Tensor<2, data_t> Bij_UU =
        TensorAlgebra::compute_dot_product(metric_UU_spatial, Bij_LU, 0, 0);
    Tensor<2, data_t> Bij_dot_Bij_UU =
        TensorAlgebra::compute_dot_product(Bij_UU, Bij_LU, 0, 0);

    FOR_ST(a)
    {
        d1_C[a] = 0.;
        FOR(i, j)
        {
            d1_C[a] += 16. * (Eij_UU[i][j] * d1_Eij[i][j][a] -
                              Eij_dot_Eij_UU[i][j] * d1_h[i][j][a]) -
                       16. * (Bij_UU[i][j] * d1_Bij[i][j][a] +
                              Bij_dot_Bij_UU[i][j] * d1_h[i][j][a]);
        }
    }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // d2_C is bad
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // indices are only spatial, derivatives are 4D
    Tensor<2, Tensor<2, data_t, CH_SPACEDIM + 1>> d2_Eij;
    Tensor<2, Tensor<2, data_t, CH_SPACEDIM + 1>> d2_Bij;
    Tensor<2, Tensor<2, data_t, CH_SPACEDIM + 1>> d2_h;

    Tensor<2, data_t> dt_Eij = gq.compute_dt_weyl_electric_part(
        d1.Eij, d1.Bij, vars.Eij, vars.Bij, advec.Eij, advec.Bij);
    Tensor<2, data_t> dt_Bij = gq.compute_dt_weyl_magnetic_part(
        d1.Eij, d1.Bij, vars.Eij, vars.Bij, advec.Eij, advec.Bij);

    FOR(i, j)
    {

        FOR(k)
        {
            // need to fill in reverse order in order to re-use them!

            FOR(l) // 2nd spatial ders
            {
                d2_Eij[i][j][k + 1][l + 1] = d2.Eij[i][j][k][l];
                d2_Bij[i][j][k + 1][l + 1] = d2.Bij[i][j][k][l];
                d2_h[i][j][k + 1][l + 1] =
                    d2.h[i][j][k][l] / vars.chi +
                    2. / (chi2 * vars.chi) * vars.h[i][j] * d1.chi[k] *
                        d1.chi[l] -
                    (d1.chi[k] * d1.h[i][j][l] + d1.chi[l] * d1.h[i][j][k] +
                     vars.h[i][j] * d2.chi[k][l]) /
                        chi2;
            }

            // mixed ders

            d2_Eij[i][j][0][k + 1] =
                -1. / m_params.tau * (d1.Eij[i][j][k] - d1.Ephys[i][j][k]);
            d2_Eij[i][j][k + 1][0] = d2_Eij[i][j][0][k + 1];

            d2_Bij[i][j][0][k + 1] =
                -1. / m_params.tau * (d1.Bij[i][j][k] - d1.Bphys[i][j][k]);
            d2_Bij[i][j][k + 1][0] = d2_Bij[i][j][0][k + 1];

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

        d2_Eij[i][j][0][0] =
            -1. / m_params.tau * (rhs.Eij[i][j] - dt_Eij[i][j]);
        d2_Bij[i][j][0][0] =
            -1. / m_params.tau * (rhs.Bij[i][j] - dt_Bij[i][j]);

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

    FOR_ST(a, b)
    {
        d2_C[a][b] = 0.;
        FOR(i, j)
        {
            d2_C[a][b] += 16. * (Eij_UU[i][j] * d2_Eij[i][j][a][b] -
                                 Eij_dot_Eij_UU[i][j] * d2_h[i][j][a][b]) -
                          16. * (Bij_UU[i][j] * d2_Bij[i][j][a][b] -
                                 Bij_dot_Bij_UU[i][j] * d2_h[i][j][a][b]);
            FOR(k, l)
            {
                d2_C[a][b] +=
                    16. *
                        (-2. * d1_Eij[i][j][a] * d1_h[k][l][b] *
                             metric_UU_spatial[i][k] * Eij_UU[j][l] -
                         2. * d1_Eij[i][j][b] * d1_h[k][l][a] *
                             metric_UU_spatial[i][k] * Eij_UU[j][l] +
                         d1_Eij[i][j][a] * d1_Eij[k][l][b] *
                             metric_UU_spatial[i][k] * metric_UU_spatial[j][l] +
                         d1_h[i][j][a] * d1_h[k][l][b] * Eij_UU[i][k] *
                             Eij_UU[j][l] +
                         d1_h[i][j][a] * d1_h[k][l][b] *
                             metric_UU_spatial[i][k] * Eij_dot_Eij_UU[j][l]) -
                    16. *
                        (-2. * d1_Bij[i][j][a] * d1_h[k][l][b] *
                             metric_UU_spatial[i][k] * Bij_UU[j][l] -
                         2. * d1_Bij[i][j][b] * d1_h[k][l][a] *
                             metric_UU_spatial[i][k] * Bij_UU[j][l] +
                         d1_Bij[i][j][a] * d1_Bij[k][l][b] *
                             metric_UU_spatial[i][k] * metric_UU_spatial[j][l] +
                         d1_h[i][j][a] * d1_h[k][l][b] * Bij_UU[i][k] *
                             Bij_UU[j][l] +
                         d1_h[i][j][a] * d1_h[k][l][b] *
                             metric_UU_spatial[i][k] * Bij_dot_Bij_UU[j][l]);
            }
        }
    }
}

template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
void SystemEB::compute_Riemann(
    Tensor<4, data_t, CH_SPACEDIM + 1> &riemann_LLLU,
    Tensor<4, data_t, CH_SPACEDIM + 1> &riemann_LULU,
    GeometricQuantities<data_t, vars_t, diff2_vars_t> &gq) const
{
    const auto &vars = gq.get_vars();
    const auto &g_UU = gq.get_metric_UU_ST();

    // valid in vacuum!
    Tensor<4, data_t, CH_SPACEDIM + 1> riemann_LLLL =
        gq.compute_weyl_tensor_LLLL(vars.Eij, vars.Bij);

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
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void SystemEB::add_matter_rhs(
    rhs_vars_t<data_t> &total_rhs,
    GeometricQuantities<data_t, vars_t, diff2_vars_t> &gq) const
{
    const auto &vars = gq.get_vars();
    const auto &d1 = gq.get_d1_vars();

    FOR(i, j)
    {
        total_rhs.Eij[i][j] =
            -1. / m_params.tau * (vars.Eij[i][j] - vars.Ephys[i][j]);
        total_rhs.Bij[i][j] =
            -1. / m_params.tau * (vars.Bij[i][j] - vars.Bphys[i][j]);
        total_rhs.Ephys[i][j] = 0.;
        total_rhs.Bphys[i][j] = 0.;
    }
}

#endif /* SYSTEMEB_IMPL_HPP_ */

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
    const auto &metric_UU_spatial = gq.get_metric_UU_spatial();

    C = 0.;
    FOR(i, j, k, l)
    {
        C +=
            8. * metric_UU_spatial[i][k] * metric_UU_spatial[j][l] *
            (vars.Eij[i][j] * vars.Eij[k][l] - vars.Bij[i][j] * vars.Bij[k][l]);
    }

    Tensor<2, Tensor<1, data_t, CH_SPACEDIM + 1>> d_Eij;
    Tensor<2, Tensor<1, data_t, CH_SPACEDIM + 1>> d_Bij;
    Tensor<2, Tensor<1, data_t, CH_SPACEDIM + 1>> d_h;

    Tensor<2, Tensor<2, data_t, CH_SPACEDIM + 1>> d2_Eij;
    Tensor<2, Tensor<2, data_t, CH_SPACEDIM + 1>> d2_Bij;
    Tensor<2, Tensor<2, data_t, CH_SPACEDIM + 1>> d2_h;

    vars_t<data_t> rhs;
    add_matter_rhs(rhs, gq);
    gq.compute_rhs_equations(rhs);

    FOR(i, j)
    {
        d_Eij[i][j][0] = rhs.Eij[i][j];
        d_Bij[i][j][0] = rhs.Bij[i][j];
        d_h[i][j][0] = vars.chi * rhs.h[i][j] + rhs.chi * vars.h[i][j];

        FOR(k)
        {
            d_Eij[i][j][k + 1] = d1.Eij[i][j][k];
            d_Bij[i][j][k + 1] = d1.Bij[i][j][k];
            d_h[i][j][k + 1] =
                vars.chi * d1.h[i][j][k] + d1.chi[k] * vars.h[i][j];
        }
    }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // MISSING DEFINING d2_Eij, d2_Bij, d2_h, d2_C
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Tensor<2, data_t> Eij_LU =
        TensorAlgebra::compute_dot_product(vars.Eij, metric_UU_spatial);
    Tensor<2, data_t> Eij_UU =
        TensorAlgebra::compute_dot_product(metric_UU_spatial, Eij_LU, 0, 0);

    Tensor<2, data_t> Bij_LU =
        TensorAlgebra::compute_dot_product(vars.Bij, metric_UU_spatial);
    Tensor<2, data_t> Bij_UU =
        TensorAlgebra::compute_dot_product(metric_UU_spatial, Bij_LU, 0, 0);

    FOR_ST(a)
    {
        d1_C[a] = 0.;
        FOR(i, j)
        {
            d1_C[a] += 16. * (Eij_UU[i][j] * d_Eij[i][j][a] -
                              Bij_UU[i][j] * d_Bij[i][j][a]);
            FOR(k)
            {
                d1_C[a] += 16. * (Eij_UU[i][j] * Eij_LU[i][k] * d_h[j][k][a] -
                                  Bij_UU[i][j] * Bij_LU[i][k] * d_h[j][k][a]);
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
    // const auto &advec = gq.get_advection();
    const auto &Eij = gq.get_weyl_electric_part();
    const auto &Bij = gq.get_weyl_magnetic_part();

    /*
        FOR(i, j)
        {
            total_rhs.Eij[i][j] =
                advec.Eij[i][j] -
                1. / m_params.tau * (vars.Eij[i][j] - Eij[i][j]) * vars.lapse;
            total_rhs.Bij[i][j] =
                advec.Bij[i][j] -
                1. / m_params.tau * (vars.Bij[i][j] - Bij[i][j]) * vars.lapse;

            FOR1(k)
            {
                total_rhs.Eij[i][j] += vars.Eij[i][k] * d1.shift[k][j] +
                                       vars.Eij[j][k] * d1.shift[k][i];
                total_rhs.Bij[i][j] += vars.Bij[i][k] * d1.shift[k][j] +
                                       vars.Bij[j][k] * d1.shift[k][i];
            }
        }
        */
    FOR(i, j)
    {
        total_rhs.Eij[i][j] = 1. / m_params.tau * (vars.Eij[i][j] - Eij[i][j]);
        total_rhs.Bij[i][j] = 1. / m_params.tau * (vars.Bij[i][j] - Bij[i][j]);
    }
}

#endif /* SYSTEMEB_IMPL_HPP_ */

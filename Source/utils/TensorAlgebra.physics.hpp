
/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(TENSORALGEBRA_HPP_)
#error "This file should only be included through TensorAlgebra.hpp"
#endif

#ifndef TENSORALGEBRA_PHYSICS_HPP_
#define TENSORALGEBRA_PHYSICS_HPP_

#include "VarsTools.hpp"

template <class data_t, int size = CH_SPACEDIM> struct chris_t
{
    Tensor<3, data_t, size> ULL;        //!< standard christoffel symbols
    Tensor<3, data_t, size> LLL;        //!< 3 lower indices
    Tensor<1, data_t, size> contracted; //!< contracted christoffel
};

template <class data_t, int size = CH_SPACEDIM> struct adm_metric_t
{
    data_t lapse;
    Tensor<1, data_t, size> shift;
    Tensor<1, data_t, size> shift_L;
    Tensor<2, data_t, size> metric_spatial;
    Tensor<2, data_t, size> metric_spatial_UU;
    Tensor<2, data_t, size> K_LL;
};

namespace TensorAlgebra
{

/// Computes the levi-civita symbol (3D, NB, symbol, not the Tensor)
inline Tensor<3, double, 3> epsilon()
{
    Tensor<3, double, 3> epsilon = {0.};
    epsilon[0][1][2] = 1.0;
    epsilon[1][2][0] = 1.0;
    epsilon[2][0][1] = 1.0;
    epsilon[0][2][1] = -1.0;
    epsilon[2][1][0] = -1.0;
    epsilon[1][0][2] = -1.0;

    return epsilon;
}

/// Computes the levi-civita symbol (4D, NB, symbol, not the Tensor)
inline Tensor<4, double, 4> epsilon4D()
{
    Tensor<4, double, 4> epsilon4D = {0.0};
    epsilon4D[0][1][2][3] = 1.0;
    epsilon4D[0][1][3][2] = -1.0;
    epsilon4D[0][3][1][2] = 1.0;
    epsilon4D[0][3][2][1] = -1.0;
    epsilon4D[0][2][1][3] = -1.0;
    epsilon4D[0][2][3][1] = 1.0;

    epsilon4D[1][0][2][3] = -1.0;
    epsilon4D[1][2][0][3] = 1.0;
    epsilon4D[1][2][3][0] = -1.0;
    epsilon4D[1][3][2][0] = 1.0;
    epsilon4D[1][3][0][2] = -1.0;
    epsilon4D[1][0][3][2] = 1.0;

    epsilon4D[2][0][1][3] = 1.0;
    epsilon4D[2][0][3][1] = -1.0;
    epsilon4D[2][3][0][1] = 1.0;
    epsilon4D[2][3][1][0] = -1.0;
    epsilon4D[2][1][3][0] = 1.0;
    epsilon4D[2][1][0][3] = -1.0;

    epsilon4D[3][0][1][2] = -1.0;
    epsilon4D[3][1][0][2] = 1.0;
    epsilon4D[3][1][2][0] = -1.0;
    epsilon4D[3][2][1][0] = 1.0;
    epsilon4D[3][2][0][1] = -1.0;
    epsilon4D[3][0][2][1] = 1.0;

    return epsilon4D;
}

/// Computes the conformal christoffel symbol
template <class data_t>
chris_t<data_t>
compute_christoffel(const Tensor<2, Tensor<1, data_t>> &d1_metric,
                    const Tensor<2, data_t> &h_UU)
{
    chris_t<data_t> out;

    FOR(i, j, k)
    {
        out.LLL[i][j][k] = 0.5 * (d1_metric[j][i][k] + d1_metric[k][i][j] -
                                  d1_metric[j][k][i]);
    }
    FOR(i, j, k)
    {
        out.ULL[i][j][k] = 0;
        FOR(l) { out.ULL[i][j][k] += h_UU[i][l] * out.LLL[l][j][k]; }
    }
    FOR(i)
    {
        out.contracted[i] = 0;
        FOR(j, k) { out.contracted[i] += h_UU[j][k] * out.ULL[i][j][k]; }
    }

    return out;
}

template <class data_t>
Tensor<3, data_t> compute_phys_chris(const Tensor<1, data_t> &d1_chi,
                                     const data_t &vars_chi,
                                     const Tensor<2, data_t> &vars_h,
                                     const Tensor<2, data_t> &h_UU,
                                     const Tensor<3, data_t> &chris_ULL)
{
    Tensor<3, data_t> chris_phys;
    FOR(i, j, k)
    {
        chris_phys[i][j][k] =
            chris_ULL[i][j][k] -
            0.5 / vars_chi *
                (delta(i, k) * d1_chi[j] + delta(i, j) * d1_chi[k]);
        FOR(m)
        {
            chris_phys[i][j][k] +=
                0.5 / vars_chi * vars_h[j][k] * h_UU[i][m] * d1_chi[m];
        }
    }
    return chris_phys;
}

// Note final indices are time for spacetime metric
template <class data_t>
adm_metric_t<data_t>
adm_vars_from_metric_ST(const Tensor<2, data_t, CH_SPACETIMEDIM> &g,
                        const Tensor<3, data_t, CH_SPACETIMEDIM> &dg)
{
    adm_metric_t<data_t> out;

    FOR(i)
    {
        // shift
        out.shift_L[i] = g[i + 1][0];
        // spatial metric
        FOR(j) { out.metric_spatial[i][j] = g[i + 1][j + 1]; }
    }

    out.metric_spatial_UU = compute_inverse_sym(out.metric_spatial);
    out.shift = raise_all(out.shift_L, out.metric_spatial_UU);

    // lapse
    data_t lapse_squared = -g[0][0];
    FOR(k) { lapse_squared += out.shift_L[k] * out.shift[k]; }
    out.lapse = sqrt(lapse_squared);

    // Kij
    Tensor<2, data_t> dt_metric_spatial;
    Tensor<1, Tensor<1, data_t>> d1_shift_L;
    Tensor<2, Tensor<1, data_t>> d1_metric_spatial;
    FOR2(i, j)
    {
        dt_metric_spatial[i][j] = dg[i + 1][j + 1][0];
        d1_shift_L[j][i] = dg[0][j + 1][i + 1];
        FOR1(k) { d1_metric_spatial[i][j][k] = dg[i + 1][j + 1][k + 1]; }
    }

    const auto chris_phys =
        compute_christoffel(d1_metric_spatial, out.metric_spatial_UU);

    Tensor<1, Tensor<1, data_t>> covd_shift_L = d1_shift_L;

    FOR(i, j, k)
    {
        covd_shift_L[i][j] += -chris_phys.ULL[k][i][j] * out.shift_L[k];
    }

    const data_t lapse_inverse = 1.0 / out.lapse;

    FOR(i, j)
    {
        out.K_LL[i][j] =
            -0.5 * lapse_inverse *
            (dt_metric_spatial[i][j] - covd_shift_L[i][j] - covd_shift_L[j][i]);
    }

    return out;
}

// Note final indices are time for spacetime metric
template <class data_t, template <typename> class vars_t>
void conformal_vars_from_adm_vars(vars_t<data_t> &vars,
                                  adm_metric_t<data_t> adm_vars)
{
    VarsTools::assign(vars, 0.); // Set only the non-zero components below

    data_t det_metric =
        TensorAlgebra::compute_determinant_sym(adm_vars.metric_spatial);
    vars.lapse = adm_vars.lapse;
    vars.shift = adm_vars.shift;
    vars.chi = pow(det_metric, -1. / GR_SPACEDIM);
    FOR(i, j) { vars.h[i][j] = vars.chi * adm_vars.metric_spatial[i][j]; }

    vars.K =
        TensorAlgebra::compute_trace(adm_vars.K_LL, adm_vars.metric_spatial_UU);
    FOR(i, j)
    {
        vars.A[i][j] =
            vars.chi * (adm_vars.K_LL[i][j] -
                        vars.K * adm_vars.metric_spatial[i][j] / GR_SPACEDIM);
    }

    // Theta = 0

    // Does NOT compute Gamma^i -> use GammaCalculator in another BoxTools in
    // initialData to compute them
}

template <class data_t, template <typename> class vars_t>
adm_metric_t<data_t> adm_vars_superposition(const vars_t<data_t> &vars1,
                                            const vars_t<data_t> &vars2)
{
    adm_metric_t<data_t> out;

    // first 3-metric
    FOR2(i, j)
    {
        out.metric_spatial[i][j] = vars1.metric_spatial[i][j] +
                                   vars2.metric_spatial[i][j] - delta(i, j);
    }
    out.metric_spatial_UU = compute_inverse_sym(out.metric_spatial);

    // gauge
    out.lapse =
        sqrt(vars1.lapse * vars1.lapse + vars2.lapse * vars2.lapse - 1.);
    FOR2(i, j)
    {
        out.shift[i] +=
            out.metric_spatial_UU[i][j] * (vars1.shift_L[j] + vars2.shift_L[j]);
    }
    out.shift_L = lower_all(out.shift, out.metric_spatial);

    // extrinsic curvature
    // raise an index before superposing
    Tensor<2, data_t> K_UL1, K_UL2, K_UL;
    FOR2(i, j)
    {
        K_UL1[i][j] = 0.;
        K_UL2[i][j] = 0.;
        K_UL[i][j] = 0.;
        out.K_LL[i][j] = 0.;
    }
    FOR3(i, j, k)
    {
        K_UL1[i][j] += vars1.metric_spatial_UU[i][k] * vars1.K_LL[k][j];
        K_UL2[i][j] += vars2.metric_spatial_UU[i][k] * vars2.K_LL[k][j];
    }
    FOR2(i, j) { K_UL[i][j] = K_UL1[i][j] + K_UL2[i][j]; }
    FOR3(i, j, k)
    {
        out.K_LL[i][j] += 0.5 * (out.metric_spatial[i][k] * K_UL[k][j] +
                                 out.metric_spatial[j][k] * K_UL[k][i]);
    }

    return out;
}

} // namespace TensorAlgebra

#endif /* TENSORALGEBRA_PHYSICS_HPP_ */

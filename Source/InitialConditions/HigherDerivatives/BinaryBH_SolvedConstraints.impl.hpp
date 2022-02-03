/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BINARYBH_SOLVEDCONSTRAINTS_HPP_)
#error                                                                         \
    "This file should only be included through BinaryBH_SolvedConstraints.hpp"
#endif

#ifndef BINARYBH_SOLVEDCONSTRAINTS_IMPL_HPP_
#define BINARYBH_SOLVEDCONSTRAINTS_IMPL_HPP_

#include "LorentzBoosts.hpp"
#include "VarsTools.hpp"

BinaryBH_SolvedConstraints::BinaryBH_SolvedConstraints(params_t a_bh1_params,
                                                       params_t a_bh2_params,
                                                       double a_dx)
    // : m_dx(a_dx), bh1(a_bh1_params, a_dx, "_1"), bh2(a_bh2_params, a_dx,
    // "_2")
    : m_dx(a_dx), m_bh1(a_bh1_params, a_dx), m_bh2(a_bh2_params, a_dx)
{
}

template <class data_t>
void BinaryBH_SolvedConstraints::compute(Cell<data_t> current_cell) const
{
    BSSNVars::VarsWithGauge<data_t> vars;
    VarsTools::assign(vars, 0.); // Set only the non-zero components below
    Coordinates<data_t> coords1_boosted(current_cell, m_dx,
                                        m_bh1.m_params.center);
    Coordinates<data_t> coords2_boosted(current_cell, m_dx,
                                        m_bh2.m_params.center);

    // compute the ADM variables from each star
    auto adm_vars1_boosted = m_bh1.compute_adm_vars(coords1_boosted);
    auto adm_vars2_boosted = m_bh2.compute_adm_vars(coords2_boosted);

    compute_vars_superposition(vars, adm_vars1_boosted, adm_vars2_boosted);

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

template <class data_t, template <typename> class vars_t>
void BinaryBH_SolvedConstraints::compute_vars_superposition(
    vars_t<data_t> &vars, const vars_t<data_t> &vars1,
    const vars_t<data_t> &vars2) const
{
    // first 3-metric
    Tensor<2, double> metric_spatial;
    FOR2(i, j)
    {
        metric_spatial[i][j] = adm_vars1.metric_spatial[i][j] +
                               adm_vars2.metric_spatial[i][j] - delta(i, j);
    }
    Tensor<2, double> metric_spatial_UU = compute_inverse_sym(metric_spatial);
    double det_metric = compute_determinant_sym(metric_spatial);
    vars.chi = pow(det_metric, -1. / GR_SPACEDIM);
    FOR2(i, j) { vars.h[i][j] = vars.chi * metric_spatial[i][j]; }

    // gauge
    vars.lapse = sqrt(adm_vars1.lapse * adm_vars1.lapse +
                      adm_vars2.lapse * adm_vars2.lapse - 1.);
    FOR2(i, j)
    {
        vars.shift[i] += metric_spatial_UU[i][j] *
                         (adm_vars1.shift_L[j] + adm_vars2.shift_L[j]);
    }

    // extrinsic curvature
    // raise an index before superposing
    Tensor<2, double> K_UL1, K_UL2, K_UL, K_LL;
    FOR2(i, j)
    {
        K_UL1[i][j] = 0.;
        K_UL2[i][j] = 0.;
        K_UL[i][j] = 0.;
        K_LL[i][j] = 0.;
    }
    FOR3(i, j, k)
    {
        K_UL1[i][j] += adm_vars1.metric_spatial_UU[i][k] * adm_vars1.K_LL[k][j];
        K_UL2[i][j] += adm_vars2.metric_spatial_UU[i][k] * adm_vars2.K_LL[k][j];
    }
    FOR2(i, j) { K_UL[i][j] = K_UL1[i][j] + K_UL2[i][j]; }
    FOR3(i, j, k)
    {
        K_LL[i][j] += 0.5 * (metric_spatial[i][k] * K_UL[k][j] +
                             metric_spatial[j][k] * K_UL[k][i]);
    }

    vars.K = compute_trace(K_LL, metric_spatial_UU);
    FOR2(i, j)
    {
        vars.A[i][j] = vars.chi * (K_LL[i][j] -
                                   vars.K * metric_spatial[i][j] / GR_SPACEDIM);
    }

    // TODO double check
}

#endif /* BINARYBH_SOLVEDCONSTRAINTS_IMPL_HPP_ */

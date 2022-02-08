/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOOSTEDSCHWARZSCHILD_SOLVEDCONSTRAINTS)
#error                                                                         \
    "This file should only be included through BoostedSchwarzschild_SolvedConstraints.hpp"
#endif

#ifndef BOOSTEDSCHWARZSCHILD_SOLVEDCONSTRAINTS_IMPL_HPP_
#define BOOSTEDSCHWARZSCHILD_SOLVEDCONSTRAINTS_IMPL_HPP_

#include "LorentzBoosts.hpp"

// Compute the value of the initial vars on the grid
template <class data_t>
void BoostedSchwarzschild_SolvedConstraints::compute(
    Cell<data_t> current_cell) const
{
    BSSNVars::VarsWithGauge<data_t> vars;
    VarsTools::assign(vars, 0.); // Set only the non-zero components below
    Coordinates<data_t> boosted_coords(current_cell, this->m_dx,
                                       this->m_params.center);

    // get boosted ADM variables
    auto adm_vars_boosted = compute_adm_boosted_vars(boosted_coords);

    // convert ADM variables to conformal variables
    compute_conformal_variables(vars, adm_vars_boosted);

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

// Compute the value of the initial vars on the grid
template <class data_t>
adm_metric_t<data_t>
BoostedSchwarzschild_SolvedConstraints::compute_adm_boosted_vars(
    const Coordinates<data_t> &boosted_coords) const
{
    // first calculate the rest frame coordinates
    Coordinates<data_t> unboosted_coords =
        LorentzBoosts::unboost_coords(boosted_coords, m_boost_velocity);

    // compute spacetime metric and derivatives
    Tensor<2, data_t, CH_SPACETIMEDIM> g_rest = {0.};
    Tensor<3, data_t, CH_SPACETIMEDIM> dg_rest = {0.};
    compute_rest_spacetime_metric(g_rest, dg_rest, unboosted_coords);

    // convert rest-frame tensors to boosted frame
    auto g_boosted = LorentzBoosts::boost_LL(g_rest, m_boost_velocity);
    auto dg_boosted = LorentzBoosts::boost_LLL(dg_rest, m_boost_velocity);

    return TensorAlgebra::adm_vars_from_metric_ST(g_boosted, dg_boosted);
}

template <class data_t>
void BoostedSchwarzschild_SolvedConstraints::compute_rest_spacetime_metric(
    Tensor<2, data_t, CH_SPACETIMEDIM> &g_rest,
    Tensor<3, data_t, CH_SPACETIMEDIM> &dg_rest,
    const Coordinates<data_t> &rest_coords) const
{
    // assumes shift = 0, lapse = sqrt(chi), h_ij=1
    // also that Arr=0 ==> K_ij=0 => d_t g_ij = 0
    // (the file doesn't assume they are 0, but it is always 0...)

    BSSNVars::VarsWithGauge<data_t> vars;
    VarsTools::assign(vars, 0.); // Set only the non-zero components below
    this->compute_vars(vars, rest_coords);

    data_t r = rest_coords.get_radius();
    data_t r_inv = 1. / r;
    data_t chi_rest_inv = 1. / vars.chi;

    // full unboosted spacetime metric (at boosted coordinates)
    FOR(i) { g_rest[i + 1][i + 1] = chi_rest_inv; }
    g_rest[0][0] = -vars.lapse * vars.lapse;

    // unit vector \hat{x}^i = x^i / r
    Tensor<1, data_t> xOr;
    xOr[0] = rest_coords.x * r_inv;
    xOr[1] = rest_coords.y * r_inv;
    xOr[2] = rest_coords.z * r_inv;

    // partial derivatives of full unboosted spacetime metric (at boosted
    // coordinates)
    // all time derivatives are 0
    data_t psi_rest = this->file_psi.interpolate(r, 0);
    data_t dpsi_rest_dr = this->file_psi.interpolate(r, 1); // 1st derivative
    data_t dchi_rest_dr = -4. * vars.chi * dpsi_rest_dr / psi_rest;
    data_t dlapse_rest_dr = 0.5 * vars.lapse * dchi_rest_dr / vars.chi;
    data_t chi_rest2_inv = chi_rest_inv * chi_rest_inv;

    FOR(i)
    {
        dg_rest[0][0][i + 1] = -2.0 * vars.lapse * dlapse_rest_dr * xOr[i];
        FOR(j)
        {
            dg_rest[j + 1][j + 1][i + 1] =
                -dchi_rest_dr * chi_rest2_inv * xOr[i];
        }
    }
}

template <class data_t, template <typename> class vars_t>
void BoostedSchwarzschild_SolvedConstraints::compute_conformal_variables(
    vars_t<data_t> &vars, adm_metric_t<data_t> adm_vars) const
{
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
}

#endif /* BOOSTEDSCHWARZSCHILD_SOLVEDCONSTRAINTS_IMPL_HPP_ */

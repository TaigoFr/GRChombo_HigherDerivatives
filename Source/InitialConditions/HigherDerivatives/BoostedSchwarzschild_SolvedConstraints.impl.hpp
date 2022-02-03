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
    auto adm_vars_boosted = compute_adm_vars(boosted_coords);

    // convert ADM variables to conformal variables
    compute_conformal_variables(vars, adm_vars_boosted);

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

// Compute the value of the initial vars on the grid
template <class data_t>
TensorAlgebra::adm_metric_t<data_t>
BoostedSchwarzschild_SolvedConstraints::compute_adm_boosted_vars(
    const Coordinates<data_t> &boosted_coords) const
{
    // first calculate the rest frame coordinates
    Coordinates<data_t> unboosted_coords =
        LorentzBoosts::unboost_coords(boosted_coords, m_boost_velocity);

    // compute spacetime metric and derivatives
    Tensor<2, data_t, CH_SPACETIMEDIM> &g;
    Tensor<3, data_t, CH_SPACETIMEDIM> &dg;
    compute_rest_spacetime_metric(g, dg, unboosted_coords);

    // convert rest-frame tensors to boosted frame
    auto g_boosted = boost_LL(g_rest, m_boost_velocity);
    auto dg_boosted = boost_LLL(dg_rest, m_boost_velocity);

    return adm_vars_from_metric_ST(g_boosted, dg_boosted);
}

template <class data_t>
void BoostedSchwarzschild_SolvedConstraints::compute_rest_spacetime_metric(
    Tensor<2, data_t, CH_SPACETIMEDIM> &g,
    Tensor<3, data_t, CH_SPACETIMEDIM> &dg,
    const Coordinates<data_t> &rest_coords) const
{
    // TODO

    // // compute the non-conformal ADM vars (metric, extrinsic K, lapse, shift)
    // this->compute_vars(vars, coords);
}

template <class data_t, template <typename> class vars_t>
void BoostedSchwarzschild_SolvedConstraints::compute_conformal_variables(
    vars_t<data_t> &vars, TensorAlgebra::adm_metric_t<data_t> adm_vars) const
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

    // TODO:
    // Compute Gamas?

    // TODO
    // check all variables are computed
}

#endif /* BOOSTEDSCHWARZSCHILD_SOLVEDCONSTRAINTS_IMPL_HPP_ */

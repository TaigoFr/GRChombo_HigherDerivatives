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
    Coordinates<data_t> coords1_boosted(current_cell, m_dx,
                                        m_bh1.m_params.center);
    Coordinates<data_t> coords2_boosted(current_cell, m_dx,
                                        m_bh2.m_params.center);

    // compute the ADM variables from each star
    auto adm_vars1_boosted = m_bh1.compute_adm_boosted_vars(coords1_boosted);
    auto adm_vars2_boosted = m_bh2.compute_adm_boosted_vars(coords2_boosted);

    auto adm_vars_superposed = TensorAlgebra::adm_vars_superposition(
        adm_vars1_boosted, adm_vars2_boosted);

    // convert ADM variables to conformal variables
    BSSNVars::VarsWithGauge<data_t> vars;
    TensorAlgebra::conformal_vars_from_adm_vars(vars, adm_vars_superposed);

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}
#endif /* BINARYBH_SOLVEDCONSTRAINTS_IMPL_HPP_ */

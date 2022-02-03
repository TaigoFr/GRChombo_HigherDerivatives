/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SCHWARZSCHILD_SOLVEDCONSTRAINTS)
#error                                                                         \
    "This file should only be included through Schwarzschild_SolvedConstraints.hpp"
#endif

#ifndef SCHWARZSCHILD_SOLVEDCONSTRAINTS_IMPL_HPP_
#define SCHWARZSCHILD_SOLVEDCONSTRAINTS_IMPL_HPP_

// Compute the value of the initial vars on the grid
template <class data_t>
void Schwarzschild_SolvedConstraints::compute(Cell<data_t> current_cell) const
{
    BSSNVars::VarsWithGauge<data_t> vars;
    VarsTools::assign(vars, 0.); // Set only the non-zero components below
    Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

    compute_vars(vars, coords);

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

template <class data_t, template <typename> class vars_t>
void Schwarzschild_SolvedConstraints::compute_vars(
    vars_t<data_t> &vars, const Coordinates<data_t> &coords) const
{
    // start with unit lapse and flat metric (must be relaxed for chi)
    vars.lapse = 1;

    // conformal metric is flat
    FOR(i) vars.h[i][i] = 1.;

    fill_from_files(vars.chi, vars.A, coords);
}

template <class data_t>
void Schwarzschild_SolvedConstraints::fill_from_files(
    data_t &chi, Tensor<2, data_t> &A, const Coordinates<data_t> &coords) const
{
    data_t r = coords.get_radius();
    static const double minimum_r = 1e-6;
    r = simd_max(r, minimum_r);

    double psi = file_psi.interpolate(r);
    if (psi >= 0.)
        chi = pow(psi, -4);
    else
        chi = 1. / pow(1. + m_params.mass / (2. * r), 4.);

    // assuming K = 0, then Arr = Krr

    data_t Arr = file_Krr.interpolate(r);

    Tensor<1, data_t> xyz = {coords.x / r, coords.y / r, coords.z / r};
    FOR(i)
    {
        FOR(j) { A[i][j] = Arr * 1.5 * xyz[i] * xyz[j]; }
        A[i][i] += -0.5 * Arr;
    }
}

#endif /* SCHWARZSCHILD_SOLVEDCONSTRAINTS_IMPL_HPP_ */

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

inline Schwarzschild_SolvedConstraints::Schwarzschild_SolvedConstraints(
    params_t a_params, double a_dx, const std::string &append)
    : m_params(a_params), m_dx(a_dx),
      file_psi("Npoints" + append, "rs" + append, "psi" + append, 1.),
      file_Krr("Npoints" + append, "rs" + append, "Krr" + append, 0.)
{
}

// Compute the value of the initial vars on the grid
template <class data_t>
void Schwarzschild_SolvedConstraints<matter_t>::compute(
    Cell<data_t> current_cell) const
{
    BSSNVars::VarsWithGauge<data_t> vars;
    VarsTools::assign(vars, 0.); // Set only the non-zero components below
    Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

    // start with unit lapse and flat metric (must be relaxed for chi)
    vars.lapse = 1;

    // conformal metric is flat
    FOR1(i) vars.h[i][i] = 1.;

    fill_from_data(vars.chi, vars.A, coords);

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

template <class data_t>
void Schwarzschild_SolvedConstraints<matter_t>::fill_from_data(
    data_t &chi, Tensor<2, data_t> &A, const Coordinates<data_t> &coords) const
{
    data_t r = coords.get_radius();
    static const double minimum_r = 1e-6;
    r = simd_max(r, minimum_r);

    double psi = file_psi.interpolate(r);
    chi = pow(psi, -4);

    // assuming K = 0, then Arr = Krr

    data_t Arr = file_Arr.interpolate(r);

    Tensor<1, data_t> xyz = {coords.x / r, coords.y / r, coords.z / r};
    FOR1(i)
    {
        FOR1(j) { A[i][j] = Arr * 1.5 * xyz[i] * xyz[j]; }
        A[i][i] += -0.5 * Arr;
    }
}

#endif /* SCHWARZSCHILD_SOLVEDCONSTRAINTS_IMPL_HPP_ */

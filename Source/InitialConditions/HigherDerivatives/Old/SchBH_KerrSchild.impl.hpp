/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SCHBH_KERRSCHILD_HPP_)
#error "This file should only be included through SchBH_KerrSchild.hpp"
#endif

#ifndef SCHBH_KERRSCHILD_IMPL_HPP_
#define SCHBH_KERRSCHILD_IMPL_HPP_

#include "BSSNVars.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "VarsTools.hpp"

template <class data_t> void SchBH::compute(Cell<data_t> current_cell) const
{
    BSSNVars::VarsWithGauge<data_t> vars;
    VarsTools::assign(vars,
                      0.); // Set only the non-zero components explicitly below
    Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

    data_t R = coords.get_radius();
    data_t M = m_params.mass;

    data_t x = coords.x;
    double y = coords.y;
    double z = coords.z;
    Tensor<1, data_t> xyz = {x, y, z};

    data_t factor = 1. + 2. * M / R;

    vars.lapse = 1. / sqrt(factor);

    FOR(i) { vars.shift[i] = 2. * M * xyz[i] / (factor * R * R); }

    vars.chi = pow(factor, -1. / GR_SPACEDIM);

    FOR(i, j)
    {
        vars.h[i][j] =
            ((i == j ? 1. : 0.) + 2. * M * xyz[i] * xyz[j] / (R * R * R));

        vars.A[i][j] = (i == j ? 2. * M / (sqrt(factor) * R * R) : 0.) -
                       2. * M * (M + 2. * R) * xyz[i] * xyz[j] /
                           (sqrt(factor) * pow(R, 5));
    }

    auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);

    vars.K = TensorAlgebra::compute_trace(vars.A, h_UU);
    TensorAlgebra::make_trace_free(vars.A, vars.h, h_UU);

    // Make conformal
    FOR(i, j)
    {
        vars.h[i][j] *= vars.chi;
        vars.A[i][j] *= vars.chi;
    }

    current_cell.store_vars(vars);
}

#endif /* SCHBH_KERRSCHILD_IMPL_HPP_ */

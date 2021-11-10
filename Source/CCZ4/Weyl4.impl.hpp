/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(WEYL4_HPP_)
#error "This file should only be included through Weyl4.hpp"
#endif

#ifndef WEYL4_IMPL_HPP_
#define WEYL4_IMPL_HPP_

template <class data_t> void Weyl4::compute(Cell<data_t> current_cell) const
{

    // copy data from chombo gridpoint into local variables
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);

    GeometricQuantities<data_t, Vars, Diff2Vars> gq(vars, d1, d2,
                                                    "Weyl4::compute");

    gq.set_formulation(m_formulation,
                       CCZ4::params_t() /*params don't matter here*/);

    // Compute the E and B fields
    EBFields_t<data_t> ebfields;
    ebfields.E = gq.get_weyl_electric_part();
    ebfields.B = gq.get_weyl_magnetic_part();

    // work out the Newman Penrose scalar
    const Coordinates<data_t> coords(current_cell, m_dx, m_center);
    NPScalar_t<data_t> out =
        compute_Weyl4(ebfields, vars, d1, d2, gq.get_h_UU(), coords);

    // Write the rhs into the output FArrayBox
    current_cell.store_vars(out.Real, c_Weyl4_Re);
    current_cell.store_vars(out.Im, c_Weyl4_Im);
}

// Calculation of the Weyl4 scalar
template <class data_t>
NPScalar_t<data_t> Weyl4::compute_Weyl4(const EBFields_t<data_t> &ebfields,
                                        const Vars<data_t> &vars,
                                        const Vars<Tensor<1, data_t>> &d1,
                                        const Diff2Vars<Tensor<2, data_t>> &d2,
                                        const Tensor<2, data_t> &h_UU,
                                        const Coordinates<data_t> &coords) const
{
    NPScalar_t<data_t> out;

    // Calculate the tetrads
    const Tetrad_t<data_t> tetrad = compute_null_tetrad(vars, h_UU, coords);

    // Projection of Electric and magnetic field components using tetrads
    out.Real = 0.0;
    out.Im = 0.0;
    FOR(i, j)
    {
        out.Real += 0.5 * (ebfields.E[i][j] * (tetrad.w[i] * tetrad.w[j] -
                                               tetrad.v[i] * tetrad.v[j]) -
                           2.0 * ebfields.B[i][j] * tetrad.w[i] * tetrad.v[j]);
        out.Im += 0.5 * (ebfields.B[i][j] * (-tetrad.w[i] * tetrad.w[j] +
                                             tetrad.v[i] * tetrad.v[j]) -
                         2.0 * ebfields.E[i][j] * tetrad.w[i] * tetrad.v[j]);
    }

    return out;
}

// Calculation of the null tetrad
// Defintions from gr-qc/0104063
// "The Lazarus project: A pragmatic approach to binary black hole evolutions",
// Baker et al.
template <class data_t>
Tetrad_t<data_t>
Weyl4::compute_null_tetrad(const Vars<data_t> &vars,
                           const Tensor<2, data_t> &h_UU,
                           const Coordinates<data_t> &coords) const
{
    Tetrad_t<data_t> out;

    // compute coords
    const data_t x = coords.x;
    const double y = coords.y;
    const double z = coords.z;

    // the alternating levi civita symbol
    const Tensor<3, double> epsilon = TensorAlgebra::epsilon();

    // calculate the tetrad
    out.u[0] = x;
    out.u[1] = y;
    out.u[2] = z;

    out.v[0] = -y;
    out.v[1] = x;
    out.v[2] = 0.0;

    out.w[0] = 0.0;
    out.w[1] = 0.0;
    out.w[2] = 0.0;

    // floor on chi
    const data_t chi = simd_max(vars.chi, 1e-4);

    FOR(i, j, k, m)
    {
        out.w[i] += 1. / sqrt(chi) * h_UU[i][j] * epsilon[j][k][m] * out.v[k] *
                    out.u[m];
    }

    // Gram Schmitt orthonormalisation
    // Choice of orthonormalisaion to avoid frame-dragging
    data_t omega_11 = 0.0;
    FOR(i, j) { omega_11 += out.v[i] * out.v[j] * vars.h[i][j] / chi; }
    FOR(i) { out.v[i] = out.v[i] / sqrt(omega_11); }

    data_t omega_12 = 0.0;
    FOR(i, j) { omega_12 += out.v[i] * out.u[j] * vars.h[i][j] / chi; }
    FOR(i) { out.u[i] += -omega_12 * out.v[i]; }

    data_t omega_22 = 0.0;
    FOR(i, j) { omega_22 += out.u[i] * out.u[j] * vars.h[i][j] / chi; }
    FOR(i) { out.u[i] = out.u[i] / sqrt(omega_22); }

    data_t omega_13 = 0.0;
    data_t omega_23 = 0.0;
    FOR(i, j)
    {
        omega_13 += out.v[i] * out.w[j] * vars.h[i][j] / chi;
        omega_23 += out.u[i] * out.w[j] * vars.h[i][j] / chi;
    }
    FOR(i) { out.w[i] += -(omega_13 * out.v[i] + omega_23 * out.u[i]); }

    data_t omega_33 = 0.0;
    FOR(i, j) { omega_33 += out.w[i] * out.w[j] * vars.h[i][j] / chi; }
    FOR(i) { out.w[i] = out.w[i] / sqrt(omega_33); }

    return out;
}

#endif /* WEYL4_HPP_ */


/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(TENSORALGEBRA_HPP_)
#error "This file should only be included through TensorAlgebra.hpp"
#endif

#ifndef TENSORALGEBRA_PHYSICS_HPP_
#define TENSORALGEBRA_PHYSICS_HPP_

template <class data_t, int size = CH_SPACEDIM> struct chris_t
{
    Tensor<3, data_t, size> ULL;        //!< standard christoffel symbols
    Tensor<3, data_t, size> LLL;        //!< 3 lower indices
    Tensor<1, data_t, size> contracted; //!< contracted christoffel
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

} // namespace TensorAlgebra

#endif /* TENSORALGEBRA_PHYSICS_HPP_ */

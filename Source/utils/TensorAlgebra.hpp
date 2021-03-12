/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef TENSORALGEBRA_HPP_
#define TENSORALGEBRA_HPP_

#include "AlwaysInline.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "simd.hpp"
#include <array>
#include <iostream>
#include <sstream>

#include "TensorAlgebra.aux.hpp"

template <class data_t, int size = CH_SPACEDIM> struct chris_t
{
    Tensor<3, data_t, size> ULL;        //!< standard christoffel symbols
    Tensor<3, data_t, size> LLL;        //!< 3 lower indices
    Tensor<1, data_t, size> contracted; //!< contracted christoffel
};

namespace TensorAlgebra
{
/// Computes determinant of a symmetric 3x3 matrix
template <class data_t>
ALWAYS_INLINE data_t compute_determinant_sym(const Tensor<2, data_t, 3> &matrix)
{
    data_t det = matrix[0][0] * matrix[1][1] * matrix[2][2] +
                 2 * matrix[0][1] * matrix[0][2] * matrix[1][2] -
                 matrix[0][0] * matrix[1][2] * matrix[1][2] -
                 matrix[1][1] * matrix[0][2] * matrix[0][2] -
                 matrix[2][2] * matrix[0][1] * matrix[0][1];

    return det;
}

/// Computes the determinant of a general 3x3 matrix.
/// Note: for a symmetric matrix use the simplified function
template <class data_t>
ALWAYS_INLINE data_t compute_determinant(const Tensor<2, data_t, 3> &matrix)
{
    data_t det =
        matrix[0][0] *
            (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
        matrix[0][1] *
            (matrix[2][2] * matrix[1][0] - matrix[1][2] * matrix[2][0]) +
        matrix[0][2] *
            (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
    return det;
}

/// Computes the inverse of a symmetric 3x3 matrix directly.
template <class data_t>
Tensor<2, data_t> compute_inverse_sym(const Tensor<2, data_t, 3> &matrix)
{
    data_t deth = compute_determinant_sym(matrix);
    data_t deth_inverse = 1. / deth;
    Tensor<2, data_t> h_UU;
    h_UU[0][0] = (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[1][2]) *
                 deth_inverse;
    h_UU[0][1] = (matrix[0][2] * matrix[1][2] - matrix[0][1] * matrix[2][2]) *
                 deth_inverse;
    h_UU[0][2] = (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) *
                 deth_inverse;
    h_UU[1][1] = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[0][2]) *
                 deth_inverse;
    h_UU[1][2] = (matrix[0][1] * matrix[0][2] - matrix[0][0] * matrix[1][2]) *
                 deth_inverse;
    h_UU[2][2] = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[0][1]) *
                 deth_inverse;
    h_UU[1][0] = h_UU[0][1];
    h_UU[2][0] = h_UU[0][2];
    h_UU[2][1] = h_UU[1][2];

    return h_UU;
}

/// Computes the inverse of a general 3x3 matrix.
/// Note: for a symmetric matrix use the simplified function
template <class data_t>
Tensor<2, data_t> compute_inverse(const Tensor<2, data_t, 3> &matrix)
{
    data_t deth = compute_determinant(matrix);
    data_t deth_inverse = 1. / deth;
    Tensor<2, data_t> h_UU;
    h_UU[0][0] = (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) *
                 deth_inverse;
    h_UU[1][1] = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) *
                 deth_inverse;
    h_UU[2][2] = (matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1]) *
                 deth_inverse;
    h_UU[1][0] = (matrix[2][0] * matrix[1][2] - matrix[1][0] * matrix[2][2]) *
                 deth_inverse;
    h_UU[0][1] = (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) *
                 deth_inverse;
    h_UU[2][0] = (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) *
                 deth_inverse;
    h_UU[0][2] = (matrix[0][1] * matrix[1][2] - matrix[1][1] * matrix[0][2]) *
                 deth_inverse;
    h_UU[2][1] = (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) *
                 deth_inverse;
    h_UU[1][2] = (matrix[1][0] * matrix[0][2] - matrix[0][0] * matrix[1][2]) *
                 deth_inverse;

    return h_UU;
}
/// Raises the index of a covector
template <class data_t, int size>
ALWAYS_INLINE Tensor<1, data_t, size>
raise_all(const Tensor<1, data_t, size> &tensor_L,
          const Tensor<2, data_t, size> &inverse_metric)
{
    Tensor<1, data_t, size> tensor_U = 0.;
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            tensor_U[i] += inverse_metric[i][j] * tensor_L[j];
    return tensor_U;
}

/// Raises the indices of a 2-Tensor
template <class data_t, int size>
ALWAYS_INLINE Tensor<2, data_t, size>
raise_all(const Tensor<2, data_t, size> &tensor_LL,
          const Tensor<2, data_t, size> &inverse_metric)
{
    Tensor<2, data_t, size> tensor_UU = 0.;

    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            for (int k = 0; k < size; ++k)
                for (int l = 0; l < size; ++l)
                    tensor_UU[i][j] += inverse_metric[i][k] *
                                       inverse_metric[j][l] * tensor_LL[k][l];
    return tensor_UU;
}

/// Lowers the indices of a vector
/// Note: same functionality as raise; included to improve readibility
template <class data_t, int size>
ALWAYS_INLINE Tensor<1, data_t, size>
lower_all(const Tensor<1, data_t, size> &tensor_U,
          const Tensor<2, data_t, size> &metric)
{ // The code for lowering is exactly the same as for raising
    return raise_all(tensor_U, metric);
}

/// Lowers the indices of a 2-Tensor
/// Note: same functionality as raise; included to improve readibility
template <class data_t, int size>
ALWAYS_INLINE Tensor<2, data_t, size>
lower_all(const Tensor<2, data_t, size> &tensor_UU,
          const Tensor<2, data_t, size> &metric)
{ // The code for lowering is exactly the same as for raising
    return raise_all(tensor_UU, metric);
}

/// Computes the (i,j) component of the Kronecker delta
constexpr int delta(int i, int j) { return (i == j); }

/// Computes the levi-civita symbol (3D, NB, symbol, not the Tensor)
inline Tensor<3, double> epsilon()
{
    Tensor<3, double> epsilon = {0.};
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

    FOR3(i, j, k)
    {
        out.LLL[i][j][k] = 0.5 * (d1_metric[j][i][k] + d1_metric[k][i][j] -
                                  d1_metric[j][k][i]);
    }
    FOR3(i, j, k)
    {
        out.ULL[i][j][k] = 0;
        FOR1(l) { out.ULL[i][j][k] += h_UU[i][l] * out.LLL[l][j][k]; }
    }
    FOR1(i)
    {
        out.contracted[i] = 0;
        FOR2(j, k) { out.contracted[i] += h_UU[j][k] * out.ULL[i][j][k]; }
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
    FOR3(i, j, k)
    {
        chris_phys[i][j][k] =
            chris_ULL[i][j][k] -
            0.5 / vars_chi *
                (delta(i, k) * d1_chi[j] + delta(i, j) * d1_chi[k]);
        FOR1(m)
        {
            chris_phys[i][j][k] +=
                0.5 / vars_chi * vars_h[j][k] * h_UU[i][m] * d1_chi[m];
        }
    }
    return chris_phys;
}

/// Makes a 2-Tensor symmetric
template <class data_t, int size>
ALWAYS_INLINE void make_symmetric(Tensor<2, data_t, size> &tensor_LL)
{
    for (int i = 0; i < size; ++i)
        for (int j = i + 1; j < size; ++j)
        {
            tensor_LL[i][j] = 0.5 * (tensor_LL[i][j] + tensor_LL[j][i]);
            tensor_LL[j][i] = tensor_LL[i][j];
        }
}

// print to any output stream 'os' (std::cout or pout() or a file)
template <class data_t, int size, int rank>
ALWAYS_INLINE void print(const Tensor<rank, data_t, size> &t,
                         std::ostream &os = std::cout)
{
    aux::print<data_t, size, rank>(
        (const typename Tensor<rank, data_t, size>::arr_t &)t, os);
}

// set a tensor to 0
template <class data_t, int size, int rank>
void set_to_zero(Tensor<rank, data_t, size> &tensor)
{
    aux::set_to_zero<data_t, size, rank>(
        (typename Tensor<rank, data_t, size>::arr_t &)tensor);
}

// copy any tensor 'src' to 'dest'
template <class data_t, int size, int rank>
void hard_copy(Tensor<rank, data_t, size> &dest,
               const Tensor<rank, data_t, size> &src)
{
    aux::hard_copy<data_t, size, rank>(
        (typename Tensor<rank, data_t, size>::arr_t &)dest,
        (const typename Tensor<rank, data_t, size>::arr_t &)src);
}

// copy any tensor 'src'*factor to 'dest'
template <bool reset, class data_t, int size, int rank>
void copy(Tensor<rank, data_t, size> &dest,
          const Tensor<rank, data_t, size> &src, data_t factor = 1.)
{
    aux::copy<data_t, size, rank, reset>(
        (typename Tensor<rank, data_t, size>::arr_t &)dest,
        (const typename Tensor<rank, data_t, size>::arr_t &)src, factor);
}

/*
// add any two tensors of same size and rank
template <class data_t, int size, int rank, long unsigned N>
ALWAYS_INLINE Tensor<rank, data_t, size>
add(const Tensor<rank, data_t, size> &tensor1,
    const Tensor<rank, data_t, size> &tensor2)
{
    Tensor<rank, data_t, size> result;
    aux::add<data_t, size, rank, false>(
        (const typename Tensor<rank, data_t, size>::arr_t &)tensor1,
        (const typename Tensor<rank, data_t, size>::arr_t &)tensor2);
    return result;
}

// add an arbitrary number of tensors of same size and rank
template <class data_t, int size, int rank, long unsigned N>
ALWAYS_INLINE Tensor<rank, data_t, size>
add(const std::array<const Tensor<rank, data_t, size> *, N> &tensors)
{
    Tensor<rank, data_t, size> result;
    aux::add<data_t, size, rank, N, false>(
        tensors, (typename Tensor<rank, data_t, size>::arr_t &)result);
    return result;
}
*/

// returns the same tensor with one rank less, by replacing index 'idx' for
// 'dir' e.g.: reduce_tensor(t[a][b][c], 1, 2) would give the tensor
// t2[a][b] = t[a][2][b]
template <class data_t, int size, int rank>
ALWAYS_INLINE Tensor<rank - 1, data_t, size>
reduce_tensor(const Tensor<rank, data_t, size> &tensor, int idx, int dir)
{
    Tensor<rank - 1, data_t, size> result;
    aux::reduce_tensor<data_t, size, rank, false>(
        (const typename Tensor<rank, data_t, size>::arr_t &)tensor,
        (typename Tensor<rank - 1, data_t, size>::arr_t &)result, idx, dir);
    return result;
}

// transpose indices 'idx1' and 'idx2' of tensor
template <class data_t, int size, int rank>
ALWAYS_INLINE Tensor<rank, data_t, size>
transpose(const Tensor<rank, data_t, size> &tensor, int idx1, int idx2)
{
    Tensor<rank, data_t, size> result;
    aux::transpose<data_t, size, rank>(
        (const typename Tensor<rank, data_t, size>::arr_t &)tensor,
        (typename Tensor<rank, data_t, size>::arr_t &)result, idx1, idx2);
    return result;
}

// compute external product between two arbitrary tensors
// e.g. t1[a][b] and t2[a][b][c] would result in a tensor
// t[a][b][c][d][e] = t1[a][b]*t2[c][d][e]
template <class data_t, int size, int rank1, int rank2>
ALWAYS_INLINE Tensor<rank1 + rank2, data_t, size>
compute_external_product(const Tensor<rank1, data_t, size> &t1,
                         const Tensor<rank2, data_t, size> &t2)
{
    Tensor<rank1 + rank2, data_t, size> result;
    aux::compute_external_product<data_t, size, rank1, rank2>(
        (const typename Tensor<rank1, data_t, size>::arr_t &)t1,
        (const typename Tensor<rank2, data_t, size>::arr_t &)t2,
        (typename Tensor<rank1 + rank2, data_t, size>::arr_t &)result);
    return result;
}

// dot product between t1 and t2, contracting index 'idx1' of t1 and 'idx2' of
// t2 (by default the last)
template <class data_t, int size, int rank1, int rank2>
ALWAYS_INLINE Tensor<rank1 + rank2 - 2, data_t, size>
compute_dot_product(const Tensor<rank1, data_t, size> &t1,
                    const Tensor<rank2, data_t, size> &t2, int idx1 = rank1 - 1,
                    int idx2 = rank2 - 1)
{
    Tensor<rank1 + rank2 - 2, data_t, size> result;
    aux::compute_dot_product<data_t, size, rank1, rank2>(
        (const typename Tensor<rank1, data_t, size>::arr_t &)t1,
        (const typename Tensor<rank2, data_t, size>::arr_t &)t2,
        (typename Tensor<rank1 + rank2 - 2, data_t, size>::arr_t &)result, idx1,
        idx2);
    return result;
}

// dot product between t1 and t2, contracting index 'idx1' of t1 and 'idx2' of
// t2 (by default the last)
template <class data_t, int size, int rank1, int rank2, int der2>
ALWAYS_INLINE Tensor<rank1 + rank2 + der2 - 2, data_t, size>
compute_dot_product(const Tensor<rank1, data_t, size> &t1,
                    const Tensor<rank2, Tensor<der2, data_t, size>, size> &t2,
                    int idx1 = rank1 - 1, int idx2 = rank2 - 1)
{
    Tensor<rank1 + rank2 + der2 - 2, data_t, size> result;
    aux::compute_dot_product<data_t, size, rank1, rank2 + der2>(
        (const typename Tensor<rank1, data_t, size>::arr_t &)t1,
        (const typename Tensor<rank2 + der2, data_t, size>::arr_t &)t2,
        (typename Tensor<rank1 + rank2 + der2 - 2, data_t, size>::arr_t &)
            result,
        idx1, idx2);
    return result;
}

// dot product between t1 and t2, contracting index 'idx1' of t1 and 'idx2' of
// t2 (by default the last) USING a 'contractor'
// this can be either the metric or the inverse metric, for UU or LL indices
template <class data_t, int size, int rank1, int rank2>
ALWAYS_INLINE Tensor<rank1 + rank2 - 2, data_t, size>
compute_dot_product(const Tensor<rank1, data_t, size> &t1,
                    const Tensor<rank2, data_t, size> &t2,
                    const Tensor<2, data_t, size> &contractor,
                    int idx1 = rank1 - 1, int idx2 = rank2 - 1)
{
    Tensor<rank1 + rank2 - 2, data_t, size> result;
    aux::compute_dot_product<data_t, size, rank1, rank2>(
        (const typename Tensor<rank1, data_t, size>::arr_t &)t1,
        (const typename Tensor<rank2, data_t, size>::arr_t &)t2,
        (typename Tensor<2, data_t, size>::arr_t &)contractor,
        (typename Tensor<rank1 + rank2 - 2, data_t, size>::arr_t &)result, idx1,
        idx2);
    return result;
}

// computes the trace of 'tensor' over indeces 'idx1' and 'idx2' (by default the
// last two)
template <class data_t, int size, int rank>
ALWAYS_INLINE Tensor<rank - 2, data_t, size>
compute_trace(const Tensor<rank, data_t, size> &tensor, int idx1 = rank - 2,
              int idx2 = rank - 1)
{
    Tensor<rank - 2, data_t, size> result;
    aux::compute_trace<data_t, size, rank>(
        (const typename Tensor<rank, data_t, size>::arr_t &)tensor,
        (typename Tensor<rank - 2, data_t, size>::arr_t &)result, idx1, idx2);
    return result;
}

// computes the trace of 'tensor' over indeces 'idx1' and 'idx2' (by default the
// last two)
// useful for derivative tensors 'Tensor<rank1, Tensor<rank2, data_t, size>'
template <class data_t, int size, int rank1, int rank2>
ALWAYS_INLINE Tensor<rank1 + rank2 - 2, data_t, size>
compute_trace(const Tensor<rank1, Tensor<rank2, data_t, size>, size> &tensor,
              int idx1 = rank1 + rank2 - 2, int idx2 = rank1 + rank2 - 1)
{
    Tensor<rank1 + rank2 - 2, data_t, size> result;
    aux::compute_trace<data_t, size, rank1 + rank2>(
        (const typename Tensor<rank1 + rank2, data_t, size>::arr_t &)tensor,
        (typename Tensor<rank1 + rank2 - 2, data_t, size>::arr_t &)result, idx1,
        idx2);
    return result;
}

// computes the trace of 'tensor' over indeces 'idx1' and 'idx2' (by default the
// last two)
template <class data_t, int size, int rank>
ALWAYS_INLINE Tensor<rank - 2, data_t, size>
compute_trace(const Tensor<rank, data_t, size> &tensor,
              const Tensor<2, data_t, size> &contractor, int idx1 = rank - 2,
              int idx2 = rank - 1)
{
    Tensor<rank - 2, data_t, size> result;
    aux::compute_trace<data_t, size, rank>(
        (const typename Tensor<rank, data_t, size>::arr_t &)tensor,
        (const typename Tensor<2, data_t, size>::arr_t &)contractor,
        (typename Tensor<rank - 2, data_t, size>::arr_t &)result, idx1, idx2);
    return result;
}

// computes the trace of 'tensor' over indeces 'idx1' and 'idx2' (by default the
// last two)
// useful for derivative tensors 'Tensor<rank1, Tensor<rank2, data_t, size>'
template <class data_t, int size, int rank1, int rank2>
ALWAYS_INLINE Tensor<rank1 + rank2 - 2, data_t, size>
compute_trace(const Tensor<rank1, Tensor<rank2, data_t, size>, size> &tensor,
              const Tensor<2, data_t, size> &contractor,
              int idx1 = rank1 + rank2 - 2, int idx2 = rank1 + rank2 - 1)
{
    Tensor<rank1 + rank2 - 2, data_t, size> result;
    aux::compute_trace<data_t, size, rank1 + rank2>(
        (const typename Tensor<rank1 + rank2, data_t, size>::arr_t &)tensor,
        (const typename Tensor<2, data_t, size>::arr_t &)contractor,
        (typename Tensor<rank1 + rank2 - 2, data_t, size>::arr_t &)result, idx1,
        idx2);
    return result;
}

/// Removes the trace of a 2-Tensor with lower indices given a metric and an
/// inverse metric.  Or a Tensor with upper indices given an inverse metric and
/// a metric.
template <class data_t, int size>
ALWAYS_INLINE void
    make_trace_free(Tensor<2, data_t, size> &tensor_LL,
                    const Tensor<2, data_t, size> &metric,
                    const Tensor<2, data_t, size> &inverse_metric)
{
    auto trace = compute_trace(tensor_LL, inverse_metric);
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            tensor_LL[i][j] -= metric[i][j] * trace / GR_SPACEDIM;
}

// moves index 'src' to position 'dest'
// e.g. move_index(t[a][b][c], 2, 0) would produce a tensor t2[a][b][c] =
// t[c][a][b]
template <class data_t, int size, int rank>
ALWAYS_INLINE Tensor<rank, data_t, size>
move_index(const Tensor<rank, data_t, size> &tensor, int src, int dest)
{
    Tensor<rank, data_t, size> result;
    aux::move_index<data_t, size, rank, false>(
        (const typename Tensor<rank, data_t, size>::arr_t &)tensor,
        (typename Tensor<rank, data_t, size>::arr_t &)result, src, dest);
    return result;
}

// compute covariant derivative of any 'tensor', given it's derivative
// 'd_tensor' and christoffel symbols 'chris_ULL' provide 'indices_up_or_down'
// if tensor contains raised indices (e.g. {1,0,0} for a tensor ULL)
// Derivative indices are last, tensor indices are first
template <class data_t, int size, int rank>
ALWAYS_INLINE Tensor<rank + 1, data_t, size>
covariant_derivative(const Tensor<rank + 1, data_t, size> &d_tensor,
                     const Tensor<rank, data_t, size> &tensor,
                     const Tensor<3, data_t, size> &chris_ULL,
                     const std::array<bool, rank> &indices_up_or_down = {false})
{
    Tensor<rank + 1, data_t, size> result;
    aux::covariant_derivative<data_t, size, rank>(
        (const typename Tensor<rank + 1, data_t, size>::arr_t &)d_tensor,
        (const typename Tensor<rank, data_t, size>::arr_t &)tensor,
        (const typename Tensor<3, data_t, size>::arr_t &)chris_ULL,
        (typename Tensor<rank + 1, data_t, size>::arr_t &)result,
        indices_up_or_down);
    return result;
}

// compute covariant derivative of any 'tensor', given it's derivative
// 'd_tensor' and christoffel symbols 'chris_ULL' provide 'indices_up_or_down'
// if tensor contains raised indices (e.g. {1,0,0} for a tensor ULL)
// useful for derivative tensors 'Tensor<rank1, Tensor<rank2, data_t, size>'
template <class data_t, int size, int rank>
ALWAYS_INLINE Tensor<rank + 1, data_t, size> covariant_derivative(
    const Tensor<rank, Tensor<1, data_t, size>, size> &d_tensor,
    const Tensor<rank, data_t, size> &tensor,
    const Tensor<3, data_t, size> &chris_ULL,
    const std::array<bool, rank> &indices_up_or_down = {false})
{
    Tensor<rank + 1, data_t, size> result;
    aux::covariant_derivative<data_t, size, rank>(
        (const typename Tensor<rank + 1, data_t, size>::arr_t &)d_tensor,
        (const typename Tensor<rank, data_t, size>::arr_t &)tensor,
        (const typename Tensor<3, data_t, size>::arr_t &)chris_ULL,
        (typename Tensor<rank + 1, data_t, size>::arr_t &)result,
        indices_up_or_down);
    return result;
}

// compute the maximum value (in absolute value) of a tensor 't' of any
// dimension or rank
template <class data_t, int size, int rank>
ALWAYS_INLINE data_t compute_tensor_max(const Tensor<rank, data_t, size> &t)
{
    return aux::compute_tensor_max<data_t, size, rank>(
        (const typename Tensor<rank, data_t, size>::arr_t &)t);
}

// assuming tensor has all indices down, creates a ST version
// (further generatilized for any dimension)
// assuming the tensor is spatial (such that time indices are simply 0)
template <class data_t, int size, int rank>
ALWAYS_INLINE Tensor<rank, data_t, size + 1>
make_spatial_tensor_up_ST(const Tensor<rank, data_t, size> &tensor)
{
    Tensor<rank, data_t, size + 1> tensor_ST;
    aux::make_spatial_tensor_up_ST<data_t, size, rank>(
        (const typename Tensor<rank, data_t, size>::arr_t &)tensor,
        (typename Tensor<rank, data_t, size + 1>::arr_t &)tensor_ST);
    return tensor_ST;
}

// assuming tensor has all indices down, creates a ST version
// (further generatilized for any dimension)
// assuming the tensor is spatial (such that time indices are simply
// contractions between the shift and the spatial ones)
template <class data_t, int size, int rank>
Tensor<rank, data_t, size + 1>
make_spatial_tensor_ST(const Tensor<rank, data_t, size> &tensor,
                       const Tensor<1, data_t, size + 1> &shift_ST)
{
    Tensor<rank, data_t, size + 1> tensor_ST;
    aux::make_spatial_tensor_ST<data_t, size, rank>(
        (const typename Tensor<rank, data_t, size>::arr_t &)tensor,
        (const typename Tensor<1, data_t, size + 1>::arr_t &)shift_ST,
        (typename Tensor<rank, data_t, size + 1>::arr_t &)tensor_ST);
    return tensor_ST;
}

template <class data_t, int size, int rank>
Tensor<rank, data_t, size + 1>
make_spatial_tensor_ST(const Tensor<rank, data_t, size> &tensor,
                       const Tensor<1, data_t, size> &shift)
{
    Tensor<1, data_t, size + 1> shift_ST;
    shift_ST[0] = 0.;
    for (int i = 0; i < size; ++i)
        shift_ST[i + 1] = shift[i];
    return make_spatial_tensor_ST(tensor, shift_ST);
}

// Lie derivative for scalar densities
template <class data_t>
ALWAYS_INLINE data_t lie_derivative(const data_t &advection_term,
                                    const data_t &tensor)
{
    return lie_derivative(advection_term, tensor, 0., 0.);
}

template <class data_t, int size>
ALWAYS_INLINE data_t lie_derivative(
    const data_t &advection_term, const data_t &tensor,
    const Tensor<1, Tensor<1, data_t, size>, size> &d_vector, double density)
{
    data_t div_vector;
    if (density == 0.)
        div_vector = 0.;
    else
        div_vector = compute_trace(d_vector);

    return lie_derivative(advection_term, tensor, div_vector, density);
}

template <class data_t>
ALWAYS_INLINE data_t lie_derivative(const data_t &advection_term,
                                    const data_t &tensor,
                                    const data_t &div_vector, double density)
{
    return advection_term + density * tensor * div_vector;
}

// Lie derivative for tensor densities
template <class data_t, int size, int rank>
ALWAYS_INLINE Tensor<rank, data_t, size>
lie_derivative(const Tensor<rank, data_t, size> &advection_term,
               const Tensor<rank, data_t, size> &tensor,
               const Tensor<1, Tensor<1, data_t, size>, size> &d_vector,
               const Tensor<1, data_t, size> &vector,
               const std::array<bool, rank> &indices_up_or_down = {false})
{
    return lie_derivative(advection_term, tensor, d_vector, vector, 0., 0.,
                          indices_up_or_down);
}

template <class data_t, int size, int rank>
ALWAYS_INLINE Tensor<rank, data_t, size>
lie_derivative(const Tensor<rank, data_t, size> &advection_term,
               const Tensor<rank, data_t, size> &tensor,
               const Tensor<1, Tensor<1, data_t, size>, size> &d_vector,
               const Tensor<1, data_t, size> &vector, double density,
               const std::array<bool, rank> &indices_up_or_down = {false})
{
    data_t div_vector;
    if (density == 0.)
        div_vector = 0.;
    else
        div_vector = compute_trace(d_vector);

    return lie_derivative(advection_term, tensor, d_vector, vector, div_vector,
                          density, indices_up_or_down);
}

template <class data_t, int size, int rank>
ALWAYS_INLINE Tensor<rank, data_t, size>
lie_derivative(const Tensor<rank, data_t, size> &advection_term,
               const Tensor<rank, data_t, size> &tensor,
               const Tensor<1, Tensor<1, data_t, size>, size> &d_vector,
               const Tensor<1, data_t, size> &vector, const data_t &div_vector,
               double density,
               const std::array<bool, rank> &indices_up_or_down = {false})
{
    Tensor<rank, data_t, size> result;
    aux::lie_derivative<data_t, size, rank>(
        (const typename Tensor<rank, data_t, size>::arr_t &)advection_term,
        (const typename Tensor<rank, data_t, size>::arr_t &)tensor,
        (const typename Tensor<2, data_t, size>::arr_t &)d_vector,
        (const typename Tensor<1, data_t, size>::arr_t &)vector, div_vector,
        (typename Tensor<rank, data_t, size>::arr_t &)result, density,
        indices_up_or_down);
    return result;
}

} // namespace TensorAlgebra

// template <class data_t, int rank, int size>
// ALWAYS_INLINE std::ostream &operator<<(std::ostream &os,
//                                        const Tensor<rank, data_t, size> &t)
// {
//     TensorAlgebra::print(t, os);
//     return os;
// }

#endif /* TENSORALGEBRA_HPP_ */

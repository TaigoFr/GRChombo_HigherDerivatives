
/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(TENSORALGEBRA_HPP_)
#error "This file should only be included through TensorAlgebra.hpp"
#endif

#ifndef TENSORALGEBRA_AUX_HPP_
#define TENSORALGEBRA_AUX_HPP_

#include "CH_assert.H"
#include <cstring> // memcpy

namespace TensorAlgebra
{

namespace aux
{

// TF: in the functions below, we use Tensor<>::arr_t explicitly instead of just
// tensor. This is mainly due to the fact that we don't want conversions from
// arr_t to Tensor when we access some element t[i] (it's what I actually tried
// first and the compiler actually complained)
// the extra benefit of it is that these functions will also work for weird
// mixed of Tensors that we get in derivatives e.g. A first derivative of a
// vector is a Tensor<1, Tensor<1, data_t, size>, size>. This is weird. But it's
// ::arr_t type is merely a data_t[size][size], just like the one of a Tensor<2,
// data_t, size>
// TF: I purposely didn't add operators to Tensors because that would make
// people use them too much, and additions/subtractions between tensors are
// costly and should be aggregated as much as possible inside loops

// for rank==0, return print
template <class data_t, int size, int rank>
ALWAYS_INLINE typename std::enable_if<(rank == 0), void>::type
print(const typename Tensor<rank, data_t, size>::arr_t &t,
      std::ostream &os = std::cout)
{
    os << t << "\t";
}

// for rank>0, print recursively through ranks and add a new line at the end of
// each rank termination
template <class data_t, int size, int rank>
ALWAYS_INLINE typename std::enable_if<(rank > 0), void>::type
print(const typename Tensor<rank, data_t, size>::arr_t &t,
      std::ostream &os = std::cout)
{
    for (int i = 0; i < size; ++i)
        print<data_t, size, rank - 1>(t[i], os);
    os << std::endl;
}

// auxiliary function for TensorAlgebra::set_to_zero for rank 0
template <class data_t, int size, int rank>
ALWAYS_INLINE typename std::enable_if<(rank == 0), void>::type
set_to_zero(typename Tensor<rank, data_t, size>::arr_t &tensor)
{
    tensor = (data_t)0.;
}

// auxiliary function for TensorAlgebra::set_to_zero for rank > 0
template <class data_t, int size, int rank>
ALWAYS_INLINE typename std::enable_if<(rank > 0), void>::type
set_to_zero(typename Tensor<rank, data_t, size>::arr_t &tensor)
{
    memset(tensor, 0, sizeof(tensor));
}

// auxiliary function for TensorAlgebra::hard_copy for rank 0
template <class data_t, int size, int rank>
ALWAYS_INLINE typename std::enable_if<(rank == 0), void>::type
hard_copy(typename Tensor<rank, data_t, size>::arr_t &dest,
          const typename Tensor<rank, data_t, size>::arr_t &src)
{
    dest = src;
}

// auxiliary function for TensorAlgebra::hard_copy for rank > 0
template <class data_t, int size, int rank>
ALWAYS_INLINE typename std::enable_if<(rank > 0), void>::type
hard_copy(typename Tensor<rank, data_t, size>::arr_t &dest,
          const typename Tensor<rank, data_t, size>::arr_t &src)
{
    memcpy(dest, src, sizeof(dest));
}

// auxiliary function for TensorAlgebra::copy for rank 0
template <class data_t, int size, int rank, bool reset>
ALWAYS_INLINE typename std::enable_if<(rank == 0 && reset == true), void>::type
copy(typename Tensor<rank, data_t, size>::arr_t &dest,
     const typename Tensor<rank, data_t, size>::arr_t &src, data_t factor = 1.)
{
    dest = src * factor;
}

// auxiliary function for TensorAlgebra::copy for rank 0
template <class data_t, int size, int rank, bool reset>
ALWAYS_INLINE typename std::enable_if<(rank == 0 && reset == false), void>::type
copy(typename Tensor<rank, data_t, size>::arr_t &dest,
     const typename Tensor<rank, data_t, size>::arr_t &src, data_t factor = 1.)
{
    dest += src * factor;
}

// auxiliary function for TensorAlgebra::copy for rank > 0
template <class data_t, int size, int rank, bool reset>
ALWAYS_INLINE typename std::enable_if<(rank > 0), void>::type
copy(typename Tensor<rank, data_t, size>::arr_t &dest,
     const typename Tensor<rank, data_t, size>::arr_t &src, data_t factor = 1.)
{
    for (int i = 0; i < size; ++i)
        copy<data_t, size, rank - 1, reset>(dest[i], src[i], factor);
}

/*
// auxiliary function for TensorAlgebra::add for rank 0
// reset == true means 'don't reset whatever is already there'
template <class data_t, int size, int rank, long unsigned N, bool reset>
typename std::enable_if<(rank == 0 && reset == true), void>::type
add(const typename Tensor<rank, data_t, size>::arr_t &tensor1,
    const typename Tensor<rank, data_t, size>::arr_t &tensor2,
    typename Tensor<rank, data_t, size>::arr_t &result)
{
    result += tensor1 + tensor2;
}

// auxiliary function for TensorAlgebra::add for rank 0
// reset == false means 'reset whatever is already there'
template <class data_t, int size, int rank, long unsigned N, bool reset>
typename std::enable_if<(rank == 0 && reset == false), void>::type
add(const typename Tensor<rank, data_t, size>::arr_t &tensor1,
    const typename Tensor<rank, data_t, size>::arr_t &tensor2,
    typename Tensor<rank, data_t, size>::arr_t &result)
{
    result = tensor1 + tensor2;
}

// auxiliary function for TensorAlgebra::add for rank > 0
template <class data_t, int size, int rank, long unsigned N, bool reset>
typename std::enable_if<(rank > 0), void>::type
add(const typename Tensor<rank, data_t, size>::arr_t &tensor1,
    const typename Tensor<rank, data_t, size>::arr_t &tensor2,
    typename Tensor<rank, data_t, size>::arr_t &result)
{
    for (int i = 0; i < size; ++i)
        add<data_t, size, rank, N, reset>(tensor1[i], tensor2[i], result[i]);
}

// auxiliary function for TensorAlgebra::add for many tensors for rank 0
template <class data_t, int size, int rank, long unsigned N, bool reset>
ALWAYS_INLINE typename std::enable_if<(rank == 0 && reset == true), void>::type
add(const std::array<const typename Tensor<rank, data_t, size>::arr_t *, N>
        &tensors,
    typename Tensor<rank, data_t, size>::arr_t &result)
{
    for (int n = 0; n < N; ++n)
        result += *tensors[n];
}

// auxiliary function for TensorAlgebra::add for many tensors for rank 0
template <class data_t, int size, int rank, long unsigned N, bool reset>
ALWAYS_INLINE typename std::enable_if<(rank == 0 && reset == false), void>::type
add(const std::array<const typename Tensor<rank, data_t, size>::arr_t *, N>
        &tensors,
    typename Tensor<rank, data_t, size>::arr_t &result)
{
    result = 0.;
    for (int n = 0; n < N; ++n)
        result += *tensors[n];
}

// auxiliary function for TensorAlgebra::add for many tensors for rank > 0
template <class data_t, int size, int rank, long unsigned N, bool reset>
typename std::enable_if<(rank > 0), void>::type
add(const std::array<const typename Tensor<rank, data_t, size>::arr_t *, N>
        &tensors,
    typename Tensor<rank, data_t, size>::arr_t &result)
{
    for (int i = 0; i < size; ++i)
    {
        std::array<const typename Tensor<rank - 1, data_t, size>::arr_t *, N>
            sub_tensors;
        for (int n = 0; n < N; ++n)
            sub_tensors[n] = &((*(tensors[n]))[i]);
        add<data_t, size, rank - 1, N, reset>(sub_tensors, result[i]);
    }
}

// auxiliary function for TensorAlgebra::add for many tensors for rank > 0
template <class data_t, int size, int rank, long unsigned N, bool reset>
typename std::enable_if<(rank > 0), void>::type
add(const std::array<const Tensor<rank, data_t, size> *, N> &tensors,
    typename Tensor<rank, data_t, size>::arr_t &result)
{
    for (int i = 0; i < size; ++i)
    {
        std::array<const typename Tensor<rank - 1, data_t, size>::arr_t *, N>
            sub_tensors;
        for (int n = 0; n < N; ++n)
            sub_tensors[n] = &((*(tensors[n]))[i]);
        add<data_t, size, rank - 1, N, reset>(sub_tensors, result[i]);
    }
}
*/

// auxiliary function for TensorAlgebra::reduce_tensor for rank 1
// that for 'reset==true', equals and replaces the old
// factor is useful when computing a trace and contracting with a metric for
// reset==false, which brings extra factors to the calculation
template <class data_t, int size, int rank, bool reset>
ALWAYS_INLINE typename std::enable_if<(rank == 1), void>::type
reduce_tensor(const typename Tensor<rank, data_t, size>::arr_t &tensor,
              typename Tensor<rank - 1, data_t, size>::arr_t &result, int idx,
              int dir, data_t factor = 1.)
{
    copy<data_t, size, rank - 1, reset>(result, tensor[dir], factor);
}

// auxiliary function for TensorAlgebra::reduce_tensor for rank > 0
// that for 'reset==true', equals and replaces the old
// factor is useful when computing a trace and contracting with a metric for
// reset==false, which brings extra factors to the calculation
template <class data_t, int size, int rank, bool reset>
typename std::enable_if<(rank > 1), void>::type
reduce_tensor(const typename Tensor<rank, data_t, size>::arr_t &tensor,
              typename Tensor<rank - 1, data_t, size>::arr_t &result, int idx,
              int dir, data_t factor = 1.)
{
    CH_assert(idx < rank && dir < size);

    if (idx > 0)
        for (int i = 0; i < size; ++i)
            reduce_tensor<data_t, size, rank - 1, reset>(tensor[i], result[i],
                                                         idx - 1, dir, factor);
    else
        copy<data_t, size, rank - 1, reset>(result, tensor[dir]);
}

// auxiliary function for TensorAlgebra::transpose for rank 1
// swaps index 0 with index 'idx', knowing that 'tensor' and 'result' already
// had their original first swapping index to 'idx1_fixed' and 'idx2_fixed'
// template <class data_t, int size, int rank>
// ALWAYS_INLINE typename std::enable_if<(rank == 1), void>::type
// transpose(const typename Tensor<rank, data_t, size>::arr_t &tensor,
//           typename Tensor<rank, data_t, size>::arr_t &result, int idx,
//           int idx1_fixed, int idx2_fixed)
// {
//     result[idx1_fixed] = tensor[idx2_fixed];
// }

// auxiliary function for TensorAlgebra::transpose for rank > 0
template <class data_t, int size, int rank>
typename std::enable_if<(rank == 0), void>::type
transpose(const typename Tensor<rank, data_t, size>::arr_t &tensor,
          typename Tensor<rank, data_t, size>::arr_t &result, int idx,
          int idx1_fixed, int idx2_fixed)
{
    CH_assert(false); // should never get here, but function needs to exist
}

// auxiliary function for TensorAlgebra::transpose for rank > 0
// swaps index 0 with index 'idx', knowing that 'tensor' and 'result' already
// had their original first swapping index to 'idx1_fixed' and 'idx2_fixed'
template <class data_t, int size, int rank>
typename std::enable_if<(rank > 0), void>::type
transpose(const typename Tensor<rank, data_t, size>::arr_t &tensor,
          typename Tensor<rank, data_t, size>::arr_t &result, int idx,
          int idx1_fixed, int idx2_fixed)
{
    if (idx == 0)
        hard_copy<data_t, size, rank - 1>(result[idx1_fixed],
                                          tensor[idx2_fixed]);
    else
        for (int i = 0; i < size; ++i)
            transpose<data_t, size, rank - 1>(tensor[i], result[i], idx - 1,
                                              idx1_fixed, idx2_fixed);
}

// auxiliary function for TensorAlgebra::transpose for rank 1
// needed for templates not to complain that it is missing
template <class data_t, int size, int rank>
ALWAYS_INLINE typename std::enable_if<(rank == 1), void>::type
transpose(const typename Tensor<rank, data_t, size>::arr_t &tensor,
          typename Tensor<rank, data_t, size>::arr_t &result, int idx1,
          int idx2)
{
    CH_assert(idx1 == 0 && idx2 == 0);
    hard_copy<data_t, size, rank>(result, tensor);
}

// auxiliary function for TensorAlgebra::transpose for rank >1
// swaps index 'idx1' with index 'idx2'
template <class data_t, int size, int rank>
typename std::enable_if<(rank >= 2), void>::type
transpose(const typename Tensor<rank, data_t, size>::arr_t &tensor,
          typename Tensor<rank, data_t, size>::arr_t &result, int idx1,
          int idx2)
{
    CH_assert(idx1 < rank && idx1 >= 0);
    CH_assert(idx2 < rank && idx2 >= 0);
    if (idx1 > idx2) // repeat with swapped indices
        transpose<data_t, size, rank>(tensor, result, idx2, idx1);
    else if (idx1 == idx2)
        hard_copy<data_t, size, rank>(result, tensor);
    else if (idx1 > 0)
        for (int i = 0; i < size; ++i)
            transpose<data_t, size, rank - 1>(tensor[i], result[i], idx1 - 1,
                                              idx2 - 1);
    else // idx1 == 0
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                transpose<data_t, size, rank - 1>(tensor[i], result[j],
                                                  idx2 - 1, i, j);
}

// auxiliary function for TensorAlgebra::compute_external_product for rank 0
// this is the final destination of the templating, that for 'reset==true',
// equals and replaces the old
// factor is useful when contracting with a metric for reset==false, which
// brings extra factors to the calculation
template <class data_t, int size, int rank1, int rank2, bool reset>
ALWAYS_INLINE typename std::enable_if<(rank1 == 0 && rank2 == 0), void>::type
compute_external_product(
    const typename Tensor<rank1, data_t, size>::arr_t &t1,
    const typename Tensor<rank2, data_t, size>::arr_t &t2,
    typename Tensor<rank1 + rank2, data_t, size>::arr_t &result,
    data_t factor = 1.)
{
    copy<data_t, size, 0, reset>(result, t1 * t2, factor);
}

// auxiliary function for TensorAlgebra::compute_external_product for rank1==0
// and rank2>0
template <class data_t, int size, int rank1, int rank2, bool reset>
ALWAYS_INLINE typename std::enable_if<(rank1 == 0 && rank2 > 0), void>::type
compute_external_product(
    const typename Tensor<rank1, data_t, size>::arr_t &t1,
    const typename Tensor<rank2, data_t, size>::arr_t &t2,
    typename Tensor<rank1 + rank2, data_t, size>::arr_t &result,
    data_t factor = 1.)
{
    for (int i = 0; i < size; ++i)
        compute_external_product<data_t, size, rank1, rank2 - 1, reset>(
            t1, t2[i], result[i], factor);
}

// auxiliary function for TensorAlgebra::compute_external_product for rank1>0
// and rank2>=0
// factor is useful when contracting with a metric for reset==false, which
// brings extra factors to the calculation
template <class data_t, int size, int rank1, int rank2, bool reset>
ALWAYS_INLINE typename std::enable_if<(rank1 > 0 && rank2 >= 0), void>::type
compute_external_product(
    const typename Tensor<rank1, data_t, size>::arr_t &t1,
    const typename Tensor<rank2, data_t, size>::arr_t &t2,
    typename Tensor<rank1 + rank2, data_t, size>::arr_t &result,
    data_t factor = 1.)
{
    for (int i = 0; i < size; ++i)
        compute_external_product<data_t, size, rank1 - 1, rank2, reset>(
            t1[i], t2, result[i], factor);
}

// auxiliary function for TensorAlgebra::compute_external_product for rank1==0
// and rank2==1
// computes external product between t1 and t2, fixing the index of t2 to
// 'fix_dir2' (in this case, fix_idx2 must ==0, but needed for the template
// structure)
template <class data_t, int size, int rank1, int rank2, bool reset>
ALWAYS_INLINE typename std::enable_if<
    (rank1 == 0 && rank2 == 1 /*&& rank1 + rank2 - 1 == 0*/), void>::type
compute_external_product(
    const typename Tensor<rank1, data_t, size>::arr_t &t1,
    const typename Tensor<rank2, data_t, size>::arr_t &t2,
    typename Tensor<rank1 + rank2 - 1, data_t, size>::arr_t &result,
    int fix_idx2, int fix_dir2, data_t factor = 1.)
{
    // CH_assert(fix_idx2 == 0 && fix_dir2 >= 0);
    // CH_assert(fix_dir2 < size);
    compute_external_product<data_t, size, rank1, rank2 - 1, reset>(
        t1, t2[fix_dir2], result, factor);
}

// auxiliary function for TensorAlgebra::compute_external_product for rank1==0
// and rank2>1
// computes external product between t1 and t2, fixing the index 'fix_idx2' of
// t2 to 'fix_dir2'
template <class data_t, int size, int rank1, int rank2, bool reset>
typename std::enable_if<(rank1 == 0 && rank2 > 1 /*&& rank1 + rank2 - 1 > 0*/),
                        void>::type
compute_external_product(
    const typename Tensor<rank1, data_t, size>::arr_t &t1,
    const typename Tensor<rank2, data_t, size>::arr_t &t2,
    typename Tensor<rank1 + rank2 - 1, data_t, size>::arr_t &result,
    int fix_idx2, int fix_dir2, data_t factor = 1.)
{
    // CH_assert(fix_idx2 >= 0 && fix_dir2 >= 0);
    // CH_assert(fix_idx2 < rank2 && fix_dir2 < size);
    if (fix_idx2 == 0)
        compute_external_product<data_t, size, rank1, rank2 - 1, reset>(
            t1, t2[fix_dir2], result, factor);
    else
        for (int i = 0; i < size; ++i)
            compute_external_product<data_t, size, rank1, rank2 - 1, reset>(
                t1, t2[i], result[i], fix_idx2 - 1, fix_dir2, factor);
}

// auxiliary function for TensorAlgebra::compute_external_product for rank1>0
// and rank2>0
// computes external product between t1 and t2, fixing the index 'fix_idx2' of
// t2 to 'fix_dir2'
// factor is useful when contracting with a metric for reset==false, which
// brings extra factors to the calculation
template <class data_t, int size, int rank1, int rank2, bool reset>
ALWAYS_INLINE
    typename std::enable_if<(rank1 > 0 && rank2 > 0 && rank1 + rank2 - 1 > 0),
                            void>::type
    compute_external_product(
        const typename Tensor<rank1, data_t, size>::arr_t &t1,
        const typename Tensor<rank2, data_t, size>::arr_t &t2,
        typename Tensor<rank1 + rank2 - 1, data_t, size>::arr_t &result,
        int fix_idx2, int fix_dir2, data_t factor = 1.)
{
    CH_assert(fix_idx2 >= 0 && fix_dir2 >= 0);
    CH_assert(fix_idx2 < rank2 && fix_dir2 < size);
    for (int i = 0; i < size; ++i)
        compute_external_product<data_t, size, rank1 - 1, rank2, reset>(
            t1[i], t2, result[i], fix_idx2, fix_dir2, factor);
}

// auxiliary function for TensorAlgebra::compute_external_product for rank1==0
// and rank2>0
// needed for templates not to complain that it is missing
template <class data_t, int size, int rank1, int rank2, bool reset>
ALWAYS_INLINE
    typename std::enable_if<(rank1 == 0 && rank2 > 0 && rank1 + rank2 - 2 >= 0),
                            void>::type
    compute_external_product(
        const typename Tensor<rank1, data_t, size>::arr_t &t1,
        const typename Tensor<rank2, data_t, size>::arr_t &t2,
        typename Tensor<rank1 + rank2 - 2, data_t, size>::arr_t &result,
        int fix_idx1, int fix_dir1, int fix_idx2, int fix_dir2, data_t factor)
{
    CH_assert(false); // should never get here, but function needs to exist
}

// auxiliary function for TensorAlgebra::compute_external_product for rank1==1
// ad rank2==1
// computes external product between t1 and t2, fixing the index of t1 to
// 'fix_dir1' and the index of t2 to 'fix_dir2' (in this case, fix_idx1 and
// fix_idx2 must ==0, but needed for the template structure)
template <class data_t, int size, int rank1, int rank2, bool reset>
ALWAYS_INLINE typename std::enable_if<
    (rank1 == 1 && rank2 == 1 /*&& rank1+rank2-2 == 0*/), void>::type
compute_external_product(
    const typename Tensor<rank1, data_t, size>::arr_t &t1,
    const typename Tensor<rank2, data_t, size>::arr_t &t2,
    typename Tensor<rank1 + rank2 - 2, data_t, size>::arr_t &result,
    int fix_idx1, int fix_dir1, int fix_idx2, int fix_dir2, data_t factor = 1.)
{
    // CH_assert(fix_idx1 == 0 && fix_idx2 == 0 && fix_dir1 >= 0 && fix_dir2 >=
    // 0);
    // CH_assert(fix_dir1 < size && fix_dir2 < size);
    compute_external_product<data_t, size, rank1 - 1, rank2, reset>(
        t1[fix_dir1], t2, result, fix_idx2, fix_dir2, factor);
}

// auxiliary function for TensorAlgebra::compute_external_product for rank1>0
// ad rank2>0
// computes external product between t1 and t2, fixing the index 'fix_idx1' of
// t1 to 'fix_dir1' and the index 'fix_idx2' of t2 to 'fix_dir2'
// factor is useful when contracting with a metric for reset==false, which
// brings extra factors to the calculation
template <class data_t, int size, int rank1, int rank2, bool reset>
typename std::enable_if<(rank1 > 0 && rank2 > 0 && rank1 + rank2 - 2 > 0),
                        void>::type
compute_external_product(
    const typename Tensor<rank1, data_t, size>::arr_t &t1,
    const typename Tensor<rank2, data_t, size>::arr_t &t2,
    typename Tensor<rank1 + rank2 - 2, data_t, size>::arr_t &result,
    int fix_idx1, int fix_dir1, int fix_idx2, int fix_dir2, data_t factor = 1.)
{
    CH_assert(fix_idx1 >= 0 && fix_idx2 >= 0 && fix_dir1 >= 0 && fix_dir2 >= 0);
    CH_assert(fix_idx1 < rank1 && fix_idx2 < rank2 && fix_dir1 < size &&
              fix_dir2 < size);
    if (fix_idx1 == 0)
        compute_external_product<data_t, size, rank1 - 1, rank2, reset>(
            t1[fix_dir1], t2, result, fix_idx2, fix_dir2, factor);
    else
        for (int i = 0; i < size; ++i)
            compute_external_product<data_t, size, rank1 - 1, rank2, reset>(
                t1[i], t2, result[i], fix_idx1 - 1, fix_dir1, fix_idx2,
                fix_dir2, factor);
}

// dot product between t1 and t2, contracting index 'idx1' of t1 and 'idx2' of
// t2 (by default the last)
// set factor to -1 to subtract (useful for covariant derivative)
template <class data_t, int size, int rank1, int rank2>
ALWAYS_INLINE void compute_dot_product(
    const typename Tensor<rank1, data_t, size>::arr_t &t1,
    const typename Tensor<rank2, data_t, size>::arr_t &t2,
    typename Tensor<rank1 + rank2 - 2, data_t, size>::arr_t &result,
    int idx1 = rank1 - 1, int idx2 = rank2 - 1, double factor = 1.)
{
    // TF: I tried many ways of doing this (3 or 4) and this one seemed the
    // fastest across multiple tests to me
    CH_assert(idx1 < rank1 && idx2 < rank2);
    set_to_zero<data_t, size, rank1 + rank2 - 2>(result);
    for (int i = 0; i < size; ++i)
        compute_external_product<data_t, size, rank1, rank2, false>(
            t1, t2, result, idx1, i, idx2, i, factor);
}

// dot product between t1 and t2, contracting index 'idx1' of t1 and 'idx2' of
// t2 (by default the last) USING a 'contractor'
// this can be either the metric or the inverse metric, for UU or LL indices
template <class data_t, int size, int rank1, int rank2>
ALWAYS_INLINE void compute_dot_product(
    const typename Tensor<rank1, data_t, size>::arr_t &t1,
    const typename Tensor<rank2, data_t, size>::arr_t &t2,
    const typename Tensor<2, data_t, size>::arr_t &contractor,
    typename Tensor<rank1 + rank2 - 2, data_t, size>::arr_t &result,
    int idx1 = rank1 - 1, int idx2 = rank2 - 1)
{
    // TF: I tried many ways of doing this (3 or 4) and this one seemed the
    // fastest across multiple tests to me
    CH_assert(idx1 < rank1 && idx2 < rank2);
    set_to_zero<data_t, size, rank1 + rank2 - 2>(result);
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            compute_external_product<data_t, size, rank1, rank2, false>(
                t1, t2, result, idx1, i, idx2, j, contractor[i][j]);
}

// auxiliary function for TensorAlgebra::compute_trace for rank==1
template <class data_t, int size, int rank>
typename std::enable_if<(rank == 2), void>::type
compute_trace(const typename Tensor<rank, data_t, size>::arr_t &tensor,
              typename Tensor<rank - 2, data_t, size>::arr_t &result,
              int idx1 = rank - 2, int idx2 = rank - 1)
{
    CH_assert(idx1 == 0 && idx2 == 1);
    set_to_zero<data_t, size, rank - 2>(result);
    for (int i = 0; i < size; ++i)
        reduce_tensor<data_t, size, rank - 1, false>(tensor[i], result,
                                                     idx2 - 1, i);
}

// auxiliary function for TensorAlgebra::compute_trace for rank >=2
// computes the trace of 't' over indeces 'idx1' and 'idx2' (by default the last
// two)
template <class data_t, int size, int rank>
typename std::enable_if<(rank > 2), void>::type
compute_trace(const typename Tensor<rank, data_t, size>::arr_t &tensor,
              typename Tensor<rank - 2, data_t, size>::arr_t &result,
              int idx1 = rank - 2, int idx2 = rank - 1)
{
    CH_assert(idx1 >= 0 && idx2 >= 0 && idx1 != idx2);
    CH_assert(idx1 < rank && idx2 < rank);
    if (idx1 > idx2)
        compute_trace<data_t, size, rank>(tensor, result, idx2, idx1);
    else if (idx1 == 0)
    {
        set_to_zero<data_t, size, rank - 2>(result);
        for (int i = 0; i < size; ++i)
            reduce_tensor<data_t, size, rank - 1, false>(tensor[i], result,
                                                         idx2 - 1, i);
    }
    else
        for (int i = 0; i < size; ++i)
            compute_trace<data_t, size, rank - 1>(tensor[i], result[i],
                                                  idx1 - 1, idx2 - 1);
}

// auxiliary function for TensorAlgebra::compute_trace for rank==2
// computes the trace of 't' over indeces 'idx1' and 'idx2' (by default the last
// two)
// USING a 'contractor' this can be either the metric or the inverse metric, for
// LL or UU indices
template <class data_t, int size, int rank>
typename std::enable_if<(rank == 2), void>::type
compute_trace(const typename Tensor<rank, data_t, size>::arr_t &tensor,
              const typename Tensor<2, data_t, size>::arr_t &contractor,
              typename Tensor<rank - 2, data_t, size>::arr_t &result,
              int idx1 = rank - 2, int idx2 = rank - 1)
{
    CH_assert(idx1 == 0 && idx2 == 1);
    set_to_zero<data_t, size, rank - 2>(result);
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            reduce_tensor<data_t, size, rank - 1, false>(
                tensor[i], result, idx2 - 1, j, contractor[i][j]);
}

// auxiliary function for TensorAlgebra::compute_trace
// computes the trace of 't' over indeces 'idx1' and 'idx2' (by default the last
// two)
// USING a 'contractor' this can be either the metric or the inverse metric, for
// LL or UU indices
template <class data_t, int size, int rank>
typename std::enable_if<(rank > 2), void>::type
compute_trace(const typename Tensor<rank, data_t, size>::arr_t &tensor,
              const typename Tensor<2, data_t, size>::arr_t &contractor,
              typename Tensor<rank - 2, data_t, size>::arr_t &result,
              int idx1 = rank - 2, int idx2 = rank - 1)
{
    CH_assert(idx1 >= 0 && idx2 >= 0 && idx1 != idx2);
    CH_assert(idx1 < rank && idx2 < rank);
    if (idx1 > idx2)
        compute_trace<data_t, size, rank>(tensor, contractor, result, idx2,
                                          idx1);
    else if (idx1 == 0)
    {
        set_to_zero<data_t, size, rank - 2>(result);
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                reduce_tensor<data_t, size, rank - 1, false>(
                    tensor[i], result, idx2 - 1, j, contractor[i][j]);
    }
    else
        for (int i = 0; i < size; ++i)
            compute_trace<data_t, size, rank - 1>(
                tensor[i], contractor, result[i], idx1 - 1, idx2 - 1);
}

// auxiliary function for TensorAlgebra::move_index for rank 1
// moves index 0 with index 'idx', knowing that 'result' already
// had their original first moving index 'dest' fixed to 'dest_fixed'
template <class data_t, int size, int rank, bool reset>
ALWAYS_INLINE typename std::enable_if<(rank == 1), void>::type
move_index(const typename Tensor<rank, data_t, size>::arr_t &tensor,
           typename Tensor<rank - 1, data_t, size>::arr_t &result, int dest,
           int dest_fixed)
{
    copy<data_t, size, rank - 1, reset>(result, tensor[dest_fixed]);
}

// auxiliary function for TensorAlgebra::move_index for rank >1
// moves index 0 with index 'idx', knowing that 'result' already
// had their original first moving index 'dest' fixed to 'dest_fixed'
template <class data_t, int size, int rank, bool reset>
typename std::enable_if<(rank > 1), void>::type
move_index(const typename Tensor<rank, data_t, size>::arr_t &tensor,
           typename Tensor<rank - 1, data_t, size>::arr_t &result, int dest,
           int dest_fixed)
{
    if (dest == 0)
        copy<data_t, size, rank - 1, reset>(result, tensor[dest_fixed]);
    else
        for (int i = 0; i < size; ++i)
            move_index<data_t, size, rank - 1, reset>(tensor[i], result[i],
                                                      dest - 1, dest_fixed);
}

// auxiliary function for TensorAlgebra::move_index for rank 1
// needed for templates not to complain that it is missing
template <class data_t, int size, int rank, bool reset>
ALWAYS_INLINE typename std::enable_if<(rank == 1), void>::type
move_index(const typename Tensor<rank, data_t, size>::arr_t &tensor,
           typename Tensor<rank, data_t, size>::arr_t &result, int src,
           int dest)
{
    CH_assert(src == 0 && dest == 0);
    copy<data_t, size, rank, reset>(result, tensor);
}

// auxiliary function for TensorAlgebra::move_index for rank > 1
// moves index 'src' to position 'dest'
// e.g. move_index(t[a][b][c], 2, 0) would produce a tensor t2[a][b][c] =
// t[c][a][b]
template <class data_t, int size, int rank, bool reset>
typename std::enable_if<(rank >= 2), void>::type
move_index(const typename Tensor<rank, data_t, size>::arr_t &tensor,
           typename Tensor<rank, data_t, size>::arr_t &result, int src,
           int dest)
{
    CH_assert(src < rank && src >= 0);
    CH_assert(dest < rank && dest >= 0);
    if (src > dest) // repeat with swapped indices
        move_index<data_t, size, rank, reset>(tensor, result, dest, src);
    else if (src == dest)
        copy<data_t, size, rank, reset>(result, tensor);
    else if (src > 0)
        for (int i = 0; i < size; ++i)
            move_index<data_t, size, rank - 1, reset>(tensor[i], result[i],
                                                      src - 1, dest - 1);
    else
    {
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                move_index<data_t, size, rank - 1, reset>(
                    tensor[j], result[i][j], dest - 1, i);
    }
}

template <class data_t, int size, int rank>
void covariant_derivative(
    const typename Tensor<rank + 1, data_t, size>::arr_t &d_tensor,
    const typename Tensor<rank, data_t, size>::arr_t &tensor,
    const typename Tensor<3, data_t, size>::arr_t &chris_ULL,
    typename Tensor<rank + 1, data_t, size>::arr_t &result,
    const std::array<bool, rank> &indices_up_or_down = {false})
{
    hard_copy<data_t, size, rank + 1>(result, d_tensor);

    for (int i = 0; i < rank; ++i)
    {
        typename Tensor<rank + 1, data_t, size>::arr_t tmp;
        if (indices_up_or_down[i])
        {
            compute_dot_product<data_t, size, rank, 3>(tensor, chris_ULL, tmp,
                                                       i, 2, 1.);
            move_index<data_t, size, rank + 1, false>(tmp, result, i, 0);
        }
        else
        {
            compute_dot_product<data_t, size, rank, 3>(tensor, chris_ULL, tmp,
                                                       i, 0, -1.);
            move_index<data_t, size, rank + 1, false>(tmp, result, i, rank - 1);
        }
    }
}

// for rank==0, return its value
template <class data_t, int size, int rank>
ALWAYS_INLINE typename std::enable_if<(rank == 0), data_t>::type
compute_tensor_max(const typename Tensor<rank, data_t, size>::arr_t &t)
{
    return t;
}

// for rank>=1, compute the maximum value of the tensor as the maximum of its
// parts
template <class data_t, int size, int rank>
typename std::enable_if<(rank > 0), data_t>::type
compute_tensor_max(const typename Tensor<rank, data_t, size>::arr_t &t)
{
    // looking for maximum in absolute value
    data_t max = 0.;
    for (int i = 0; i < size; ++i)
    {
        data_t tmp = compute_tensor_max<data_t, size, rank - 1>(t[i]);
        tmp = simd_conditional(simd_compare_gt(tmp, 0.), tmp, -tmp);
        max = simd_max(max, tmp);
    }
    return max;
}

// assuming tensor has all indices down, creates a ST version
// (further generatilized for any dimension)
// assuming the tensor is spatial (such that time indices are simply 0)
template <class data_t, int size, int rank>
ALWAYS_INLINE typename std::enable_if<(rank == 0), void>::type
make_spatial_tensor_up_ST(
    const typename Tensor<rank, data_t, size>::arr_t &tensor,
    typename Tensor<rank, data_t, size + 1>::arr_t &result)
{
    result = tensor;
}

template <class data_t, int size, int rank>
ALWAYS_INLINE typename std::enable_if<(rank > 0), void>::type
make_spatial_tensor_up_ST(
    const typename Tensor<rank, data_t, size>::arr_t &tensor,
    typename Tensor<rank, data_t, size + 1>::arr_t &result)
{
    set_to_zero<data_t, size + 1, rank - 1>(result[0]);
    for (int i = 0; i < size; ++i)
        make_spatial_tensor_up_ST<data_t, size, rank - 1>(tensor[i],
                                                          result[i + 1]);
}

// assuming tensor has all indices down, creates a ST version
// (further generatilized for any dimension)
// assuming the tensor is spatial (such that time indices are simply
// contractions between the shift and the spatial ones)
template <class data_t, int size, int rank>
ALWAYS_INLINE typename std::enable_if<(rank == 1), void>::type
make_spatial_tensor_ST(
    const typename Tensor<rank, data_t, size>::arr_t &tensor,
    const typename Tensor<1, data_t, size + 1>::arr_t &shift_ST,
    typename Tensor<rank, data_t, size + 1>::arr_t &result)
{
    for (int i = 0; i < size; ++i)
        result[i + 1] = tensor[i];
    set_to_zero<data_t, size + 1, rank - 1>(result[0]);
    compute_dot_product<data_t, size + 1, rank, 1>(result, shift_ST, result[0],
                                                   0, 0);
}

template <class data_t, int size, int rank>
ALWAYS_INLINE typename std::enable_if<(rank >= 2), void>::type
make_spatial_tensor_ST(
    const typename Tensor<rank, data_t, size>::arr_t &tensor,
    const typename Tensor<1, data_t, size + 1>::arr_t &shift_ST,
    typename Tensor<rank, data_t, size + 1>::arr_t &result)
{
    for (int i = 0; i < size; ++i)
        make_spatial_tensor_ST<data_t, size, rank - 1>(tensor[i], shift_ST,
                                                       result[i + 1]);
    set_to_zero<data_t, size + 1, rank - 1>(result[0]);
    compute_dot_product<data_t, size + 1, rank, 1>(result, shift_ST, result[0],
                                                   0, 0);
}

template <class data_t, int size, int rank>
ALWAYS_INLINE typename std::enable_if<(rank > 0), void>::type
lie_derivative(const typename Tensor<rank, data_t, size>::arr_t &advection_term,
               const typename Tensor<rank, data_t, size>::arr_t &tensor,
               const typename Tensor<2, data_t, size>::arr_t &d_vector,
               const typename Tensor<1, data_t, size>::arr_t &vector,
               const data_t &div_vector,
               typename Tensor<rank, data_t, size>::arr_t &result,
               double density,
               const std::array<bool, rank> &indices_up_or_down = {false})
{
    hard_copy<data_t, size, rank>(result, advection_term);

    for (int i = 0; i < rank; ++i)
    {
        typename Tensor<rank, data_t, size>::arr_t tmp;
        if (indices_up_or_down[i])
        {
            compute_dot_product<data_t, size, 2, rank>(d_vector, tensor, tmp, 1,
                                                       i, -1.);
            move_index<data_t, size, rank, false>(tmp, result, 0, i);
        }
        else
        {
            compute_dot_product<data_t, size, 2, rank>(d_vector, tensor, tmp, 0,
                                                       i, 1.);
            move_index<data_t, size, rank, false>(tmp, result, 0, i);
        }
    }

    if (density != 0.)
        copy<data_t, size, rank, false>(result, tensor, density * div_vector);
}

} // namespace aux

} // namespace TensorAlgebra

#endif /* TENSORALGEBRA_AUX_HPP_ */

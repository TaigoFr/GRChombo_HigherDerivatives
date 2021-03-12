/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef TENSOR_HPP_
#define TENSOR_HPP_

#include "AlwaysInline.hpp"
#include "DimensionDefinitions.hpp"

/// This class implements a Tensor with given rank, element data type, and
/// dimension.  By default the dimension is equal to DEFAULT_TENSOR_DIM.
template <int rank, class data_t, int size = DEFAULT_TENSOR_DIM> class Tensor
{
  public:
    template <int, class, int> friend class Tensor;
    typedef typename Tensor<rank - 1, data_t, size>::arr_t arr_t[size];

  private:
    arr_t arr;

  public:
    ALWAYS_INLINE
    Tensor() {}

    //    ALWAYS_INLINE
    //    Tensor(std::initializer_list<data_t> list) :
    //        arr (list)
    //    {}

    template <typename... T> ALWAYS_INLINE Tensor(T... data) : arr{data...} {}

    // to convert to 'arr_t&' implicitly
    operator arr_t &() { return arr; }
    operator const arr_t &() const { return arr; }

    // to convert to 'arr_t&' explicitly (e.g. t() gives t.arr)
    arr_t &operator()() { return arr; }
    const arr_t &operator()() const { return arr; }
};

template <class data_t, int size> class Tensor<0, data_t, size>
{
  public:
    template <int, class, int> friend class Tensor;
    typedef data_t arr_t;

  private:
    arr_t arr;

  public:
    ALWAYS_INLINE
    Tensor() {}

    ALWAYS_INLINE
    Tensor(data_t val) : arr(val) {}

    // to convert to 'arr_t&' implicitly
    operator arr_t &() { return arr; }
    operator const arr_t &() const { return arr; }

    // to convert to 'arr_t&' explicitly (e.g. t() gives t.arr)
    arr_t &operator()() { return arr; }
    const arr_t &operator()() const { return arr; }
};

#endif /* TENSOR_HPP_ */

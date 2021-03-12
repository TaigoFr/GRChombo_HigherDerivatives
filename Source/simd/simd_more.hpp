
#ifndef SIMD_MORE_HPP_
#define SIMD_MORE_HPP_

#if !defined(SIMD_HPP_)
#error "This file should only be included through simd.hpp"
#endif

template <typename t>
ALWAYS_INLINE bool simd_compare_lt_any(const t &a, const t &b)
{
    return a < b;
}

// true if any of the simd is less than 'b'
template <typename t>
ALWAYS_INLINE bool simd_compare_lt_any(const simd<t> &a, const t &b)
{
    t in_arr[simd_traits<t>::simd_len];
    simd<t>::store(in_arr, a);

    bool out_bool = false;
    for (int i = 0; i < simd_traits<t>::simd_len; ++i)
        out_bool = out_bool || (in_arr[i] < b);

    return out_bool;
}

template <typename t>
ALWAYS_INLINE bool simd_compare_gt_any(const t &a, const t &b)
{
    return a > b;
}

// true if any of the simd is greater than 'b'
template <typename t>
ALWAYS_INLINE bool simd_compare_gt_any(const simd<t> &a, const t &b)
{
    t in_arr[simd_traits<t>::simd_len];
    simd<t>::store(in_arr, a);

    bool out_bool = false;
    for (int i = 0; i < simd_traits<t>::simd_len; ++i)
        out_bool = out_bool || (in_arr[i] > b);

    return out_bool;
}

template <typename t> ALWAYS_INLINE bool simd_isnan_any(const t &a)
{
    return std::isnan(a);
}

// true if any of the simd is nan
template <typename t> ALWAYS_INLINE bool simd_isnan_any(const simd<t> &a)
{
    t in_arr[simd_traits<t>::simd_len];
    simd<t>::store(in_arr, a);

    bool out_bool = false;
    for (int i = 0; i < simd_traits<t>::simd_len; ++i)
        out_bool = out_bool || std::isnan(in_arr[i]);

    return out_bool;
}

#endif /* SIMD_MORE_HPP_ */
